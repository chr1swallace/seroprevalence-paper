source("common.R")
library(ggplot2)
library(cowplot)

(load(file_rdata_v5)) #m
m[, type := make.names(type)]
#type: historical controls = controls, covid = cases.

#only keep historical and COVID
dat <- m[!(type %in% c("No.sample.control", "Patient.4")),]
dat[, c("RBD.log", "SPIKE.log") := list(log(RBD), log(SPIKE))]
train.ind <- dat$type %in% c("COVID", "Historical.controls")

rbd <- dat[type == "Historical.controls"]$RBD
spk <- dat[type == "Historical.controls"]$SPIKE

pred.rbd.3sd <- ifelse(dat[!train.ind]$RBD < mean(rbd) + 3 * sd(rbd), "Negative.pred", "COVID.pred")
pred.spk.3sd <- ifelse(dat[!train.ind]$SPIKE < mean(spk) + 3 * sd(spk), "Negative.pred", "COVID.pred")
pred.rbd.6sd <- ifelse(dat[!train.ind]$RBD < mean(rbd) + 6 * sd(rbd), "Negative.pred", "COVID.pred")
pred.spk.6sd <- ifelse(dat[!train.ind]$SPIKE < mean(spk) + 6 * sd(spk), "Negative.pred", "COVID.pred")


table(pred.rbd.3sd, dat[!train.ind]$type)
table(pred.spk.3sd, dat[!train.ind]$type)
table(pred.rbd.6sd, dat[!train.ind]$type)
table(pred.spk.6sd, dat[!train.ind]$type)


###ML ANALYSIS
run.analysis <- function(dat) {
  train.ind <- dat$type %in% c("COVID", "Historical.controls")

  #----------------------------------------------------------------------------------------
  ##LOGISTIC REGRESSION
  message("startig logreg...")
  mod.log.rbd <- glm(data = dat[train.ind,], factor(type) ~ log(RBD), family = binomial, maxit = 100)
  prob.log.rbd <- 1 - predict(mod.log.rbd, dat[!train.ind], type = "response") #prob of covid
  pred.log.rbd <- ifelse(prob.log.rbd > 0.5, "COVID.pred", "Negative.pred")

  mod.log.spk <- glm(data = dat[train.ind,], factor(type) ~ log(SPIKE), family = binomial, maxit = 100)
  prob.log.spk <- 1 - predict(mod.log.spk, dat[!train.ind], type = "response") #prob of covid
  pred.log.spk <- ifelse(prob.log.spk > 0.5, "COVID.pred", "Negative.pred")

  message("logreg done!")
  #----------------------------------------------------------------------------------------
  ##SVM: linear kernel
  message("starting linear svm...")
  ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 5, classProbs = TRUE, allowParallel = TRUE)
  grd <- data.frame(C = c(0.001, 0.01, 0.5, 1, 2, 5, 10))

  mod.svm.rbd <- train(data = dat[train.ind,], factor(type) ~ log(RBD), method = "svmLinear", trControl = ctrl, tuneGrid = grd)
  prob.svm.rbd <- predict(mod.svm.rbd, dat[!train.ind,], type = "prob")$COVID
  pred.svm.rbd <- ifelse(prob.svm.rbd > 0.5, "COVID.pred", "Negative.pred")

  mod.svm.spk <- train(data = dat[train.ind,], factor(type) ~ log(SPIKE), method = "svmLinear", trControl = ctrl, tuneGrid = grd)
  prob.svm.spk <- predict(mod.svm.spk, dat[!train.ind,], type = "prob")$COVID
  pred.svm.spk <- ifelse(prob.svm.spk > 0.5, "COVID.pred", "Negative.pred")

  message("linear svm done!")
  #----------------------------------------------------------------------------------------
  #LDA
  message("starting lda...")
  mod.lda.rbd <- train(data = dat[train.ind,], factor(type) ~ log(RBD), method = "lda")
  prob.lda.rbd <- predict(mod.lda.rbd, dat[!train.ind,], type = "prob")$COVID
  pred.lda.rbd <- ifelse(prob.lda.rbd > 0.5, "COVID.pred", "Negative.pred")

  mod.lda.spk <- train(data = dat[train.ind,], factor(type) ~ log(SPIKE), method = "lda")
  prob.lda.spk <- predict(mod.lda.spk, dat[!train.ind,], type = "prob")$COVID
  pred.lda.spk <- ifelse(prob.lda.spk > 0.5, "COVID.pred", "Negative.pred")

  message("lda done!")
  #----------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------
  #SAVE
  res <- data.table(dat[!train.ind, 1 : 8], prob.SVM.rbd = prob.svm.rbd, status.SVM.rbd = pred.svm.rbd,
                                            prob.LDA.rbd = prob.lda.rbd, status.LDA.rbd = pred.lda.rbd,
                                            prob.LOG.rbd = prob.log.rbd, status.LOG.rbd = pred.log.rbd,
                                            prob.SVM.spk = prob.svm.spk, status.SVM.spk = pred.svm.spk,
                                            prob.LDA.spk = prob.lda.spk, status.LDA.spk = pred.lda.spk,
                                            prob.LOG.spk = prob.log.spk, status.LOG.spk = pred.log.spk)

  message("all done!")
  list(res = res, mods.spk = list(mod.lda = mod.lda.spk, mod.svm = mod.svm.spk, mod.log = mod.log.spk), mods.rbd = list(mod.lda = mod.lda.rbd, mod.svm = mod.svm.rbd, mod.log = mod.log.rbd))
}


res <- run.analysis(dat)
res.x <- run.analysis(dat[!(dat$Sample.ID %in% c("1", "2")), ])

mod.rbd <- res$mods.rbd
mod.spk <- res$mod.spk
mod.rbd.x <- res.x$mods.rbd
mod.spk.x <- res.x$mod.spk

save(mod.rbd, mod.spk, mod.rbd.x, mod.spk.x, file = "~/rds/rds-cew54-wallace-share/Data/elisa/full_univar_mods.Rdata")

res <- res$res
res.x <- res.x$res

res[ , c("prob.ENS.rbd", "prob.ENS.spk") := list((prob.LDA.rbd + prob.SVM.rbd)/2, (prob.LDA.spk + prob.SVM.spk)/2)]
res[, c("status.ENS.rbd", "status.ENS.spk") := list(ifelse(prob.ENS.rbd > 0.5, "COVID.pred", "Negative.pred"), ifelse(prob.ENS.spk > 0.5, "COVID.pred", "Negative.pred"))]

res[, c("status.SD3.rbd", "status.SD3.spk", "status.6SD.rbd", "status.6SD.spk") := list(pred.rbd.3sd, pred.spk.3sd, pred.rbd.6sd, pred.spk.6sd)]

res.x[ , c("prob.ENS.rbd", "prob.ENS.spk") := list((prob.LDA.rbd + prob.SVM.rbd)/2, (prob.LDA.spk + prob.SVM.spk)/2)]
res.x[, c("status.ENS.rbd", "status.ENS.spk") := list(ifelse(prob.ENS.rbd > 0.5, "COVID.pred", "Negative.pred"), ifelse(prob.ENS.spk > 0.5, "COVID.pred", "Negative.pred"))]


save(res, res.x, file = "~/rds/rds-cew54-wallace-share/Data/elisa/full_univar_preds.Rdata")




















##
