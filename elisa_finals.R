source("common.R")
library(ggplot2)
library(cowplot)
library(kernlab)
library(caret)

# d <- "/rds/user/ng414/hpc-work/elisa/results"

(load(file_rdata)) #m
m[, type := make.names(type)]

dat <- m[!(type %in% c("No.sample.control", "Patient.4")),]

run.analysis <- function(dat) {
  train.ind <- dat$type %in% c("COVID", "Historical.controls")

  #----------------------------------------------------------------------------------------
  ##LOGISTIC REGRESSION
  message("startig logreg...")
  mod.log <- glm(data = dat[train.ind,], factor(type) ~ log(SPIKE) + log(RBD), family = binomial, maxit = 100)
  prob.log <- 1 - predict(mod.log, dat[!train.ind], type = "response") #prob of covid
  pred.log <- ifelse(prob.log > 0.5, "COVID.pred", "Negative.pred")

  message("logreg done!")
  #----------------------------------------------------------------------------------------
  ##SVM: linear kernel
  message("starting linear svm...")
  ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 5, classProbs = TRUE, allowParallel = TRUE)
  grd <- data.frame(C = c(0.001, 0.01, 0.5, 1, 2, 5, 10))

  mod.svm <- train(data = dat[train.ind,], factor(type) ~ log(SPIKE) + log(RBD), method = "svmLinear", trControl = ctrl, tuneGrid = grd)
  prob.svm <- predict(mod.svm, dat[!train.ind,], type = "prob")$COVID
  pred.svm <- ifelse(prob.svm > 0.5, "COVID.pred", "Negative.pred")

  message("linear svm done!")
  #----------------------------------------------------------------------------------------
  ##SVM: quadratic kernel
  message("starting quadratic svm...")
  ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 5, classProbs = TRUE, allowParallel = TRUE)
  grd <- expand.grid(degree = 2, scale = c(0.1, 1, 2), C = c(0.001, 0.01, 0.5, 1, 2, 5, 10))

  mod.svm2 <- train(data = dat[train.ind,], factor(type) ~ log(SPIKE) + log(RBD), method = "svmPoly", trControl = ctrl, tuneGrid = grd)
  prob.svm2 <- predict(mod.svm2, dat[!train.ind,], type = "prob")$COVID
  pred.svm2 <- ifelse(prob.svm2 > 0.5, "COVID.pred", "Negative.pred")

  message("quadratic svm done!")
  #----------------------------------------------------------------------------------------
  #LDA
  message("starting lda...")
  mod.lda <- train(data = dat[train.ind,], factor(type) ~ log(SPIKE) + log(RBD), method = "lda")
  prob.lda <- predict(mod.lda, dat[!train.ind,], type = "prob")$COVID
  pred.lda <- ifelse(prob.lda > 0.5, "COVID.pred", "Negative.pred")

  message("lda done!")
  #----------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------
  #SAVE
  res <- data.table(dat[!train.ind, 1 : 8], prob.SVM = prob.svm, status.SVM = pred.svm, prob.LDA = prob.lda, status.LDA = pred.lda, prob.LOG = prob.log, status.LOG = pred.log, prob.SVM2 = prob.svm2, status.SVM2 = pred.svm2)

  message("all done!")
  list(res = res, mod.lda = mod.lda, mod.svm = mod.svm, mod.log = mod.log, mod.svm2 = mod.svm2)
}

res <- run.analysis(dat)
lda.mod <- res$mod.lda
svm.mod <- res$mod.svm
log.mod <- res$mod.log
svm2.mod <- res$mod.svm2

res <- res$res[, prob.ENS := (prob.SVM + prob.LDA)/2]
res[, status.ENS := ifelse(prob.ENS > 0.5, "COVID.pred", "Negative.pred")]

# res <- res$res[, prob.ENS.LOG := (prob.SVM + prob.LOG)/2]
# res[, status.ENS.LOG := ifelse(prob.ENS.LOG > 0.5, "COVID.pred", "Negative.pred")]

saveRDS(res, file = "~/rds/rds-cew54-wallace-share/Data/elisa/predictions_on_test_set.rds")
save(lda.mod, svm.mod, log.mod, svm2.mod, file = "~/rds/rds-cew54-wallace-share/Data/elisa/full_models.RData")
rm(lda.mod, svm.mod, log.mod, svm2.mod)


res.x <- run.analysis(dat[!(dat$Sample.ID %in% c("1", "2")), ])
lda.mod <- res.x$mod.lda
svm.mod <- res.x$mod.svm
log.mod <- res.x$mod.log
svm2.mod <- res.x$mod.svm2

res.x <- res.x$res[, prob.ENS := (prob.SVM + prob.LDA)/2]
res.x[, status.ENS := ifelse(prob.ENS > 0.5, "COVID.pred", "Negative.pred")]

# res.x <- res.x$res[, prob.ENS.LOG := (prob.SVM + prob.LOG)/2]
# res.x[, status.ENS.LOG := ifelse(prob.ENS.LOG > 0.5, "COVID.pred", "Negative.pred")]

saveRDS(res.x, file = "~/rds/rds-cew54-wallace-share/Data/elisa/predictions_on_test_set_no_12.rds")
save(lda.mod, svm.mod, log.mod, svm2.mod, file = "~/rds/rds-cew54-wallace-share/Data/elisa/no_12_models.RData")

#END OF MODEL FITTING


## #PLOT
## res <- melt(res, measure.vars = c("status.SVM", "status.LDA", "status.LOG"))
## res[, method := sub("status\\.", "", variable)]
## setnames(res, "value", "status")

## ddat <- dat[train.ind, 1 : 8]
## ddat[, status := type]
## ddat <- rbind(ddat, ddat, ddat)
## ddat[, method := rep(c("LDA", "SVM", "LOG"), len = nrow(ddat))]

## X <- rbind(res[, c("SPIKE", "RBD", "method", "status")], ddat[, c("SPIKE", "RBD", "method", "status")])

## ggplot(X, aes(x = log(SPIKE), y = log(RBD), colour = status)) + geom_point(alpha = 0.6, size = 2) + facet_grid(.~method) + theme_cowplot(font_size = 25) + background_grid()



































##
