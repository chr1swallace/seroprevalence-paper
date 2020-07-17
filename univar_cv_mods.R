source("common.R")
library(ggplot2)
library(cowplot)
library(kernlab)
library(caret)
library(parallel)

file_rdata <- file_rdata_v5
(load(file_rdata)) #m
m[, type:=make.names(type)]
#type: historical controls = controls, covid = cases.

#only keep historical and COVID
dat <- m[!is.na(fold1)]
table(m$type)
dat.x <- dat[!(Sample.ID %in% c("1", "2"))]
table(m$type)

#----------------------------------------------------------------------------------------------------------
#LOGISTIC REGRESSION
#logistic: 0 = COVID, 1 = historical
#prot = {RBD, SPIKE}
run.logreg <- function(dat, k, prot) {
  fold <- paste0("fold", k)
  out <- data.table(dat[, c("Sample.ID", "type")], fold = dat[, get(fold)], method = "LOG", fold.n = k, protein = prot, prob.covid = vector("numeric", nrow(dat)), pred = vector("character", nrow(dat)))

  mod.output <- lapply(1 : 10, FUN = function(i) {
    message(sprintf("fold%s: %s", k, i))
    if(prot == "RBD") {
      mod <- glm(data = dat[get(fold) != i,], factor(type) ~ log(RBD), family = binomial, maxit = 100)
    } else {
        mod <- glm(data = dat[get(fold) != i,], factor(type) ~ log(SPIKE), family = binomial, maxit = 100)
  }
    prob <- 1 - predict(mod, dat[get(fold) == i], type = "response") #prob of covid
    list(mod = mod, prob = prob)
  })

  for(i in 1 : 10) {
    message(i)
    out[fold == i, prob.covid := mod.output[[i]]$prob]
    out[fold == i, pred := ifelse(prob.covid > 0.5, "COVID", "Historical.controls")]
  }

  list(mod.output = mod.output, pred = out)
}

#----------------------------------------------------------------------------------------------------------
#ctrl = options for parameter tuning, leave as NULL if don't want to tune
#grd = parameter grid for tuning, leave as NULL if don't want to tune
#log_ = log SPIKE and RDS?
#k = which fold set?
#prot = {RBD, SPIKE}
run.mod <- function(dat, k, method, grd = NULL, ctrl = NULL, prot) {
  fold <- paste0("fold", k)
  out <- data.table(dat[, c("Sample.ID", "type")], fold = dat[, get(fold)], method = method, fold.n = k, protein = prot, prob.covid = vector("numeric", nrow(dat)), pred = vector("character", nrow(dat)))

  mod.output <- lapply(1 : 10, FUN = function(i) {
    message(sprintf("%s: fold%s %s", method, k, i))
    dat.train <- dat[get(fold) != i,]
    dat.test <- dat[get(fold) == i,]

    if(all(is.null(grd), is.null(ctrl))) {
      if(prot == "RBD") {
        mod <- train(data = dat.train, factor(type) ~ log(RBD), method = method) } else {
                mod <- train(data = dat.train, factor(type) ~ log(SPIKE), method = method)
        }
    } else {
      if(prot == "RBD") {
        mod <- train(data = dat.train, factor(type) ~ log(RBD), method = method, trControl = ctrl, tuneGrid = grd) } else {
                mod <- train(data = dat.train, factor(type) ~ log(SPIKE), method = method, trControl = ctrl, tuneGrid = grd)
        }
    }

    pred <- predict(mod, dat.test, type = "prob")
    list(mod = mod, prob = pred$COVID)
  })

  for(i in 1 : 10) {
    message(i)
    out[fold == i, prob.covid := mod.output[[i]]$prob]
    out[fold == i, pred := ifelse(prob.covid > 0.5, "COVID", "Historical.controls")]
  }

  list(mod.output = mod.output, pred = out)
}

#----------------------------------------------------------------------------------------------------------
#3/6 SD

run.sd <- function(dat, k, prot, d) {
  fold <- paste0("fold", k)
  out <- data.table(dat[,c("Sample.ID", "type")], fold = dat[,get(fold)], method = paste0("SD", k), fold.n = k, protein = prot, prob.covid = vector("numeric", nrow(dat)), pred = vector("character", nrow(dat)))

  mod.output <- lapply(1 : 10, FUN = function(i) {
    message(sprintf("fold%s: %s", k, i))
    if(prot == "RBD") {
      rbd <- dat[get(fold) != i & type == "Historical.controls",]$RBD
      thr <- mean(rbd) + d * sd(rbd)
      pred <- ifelse(dat[get(fold) == i]$RBD < thr, "Negative.pred", "COVID.pred")
        } else {
          spk <- dat[get(fold) != i & type == "Historical.controls",]$SPIKE
          thr <- mean(spk) + d * sd(spk)
          pred <- ifelse(dat[get(fold) == i]$SPIKE < thr, "Negative.pred", "COVID.pred")
        }
      pred
  })

  for(i in 1 : 10) {
    message(i)
    out[fold == i, pred := mod.output[[i]]]
  }

  out
}


#----------------------------------------------------------------------------------------------------------
##LOGREG

logreg.rbd <- lapply(1 : NFOLDS, FUN = function(i) run.logreg(dat, i, "RBD"))
logreg.spk <- lapply(1 : NFOLDS, FUN = function(i) run.logreg(dat, i, "SPIKE"))
logreg.rbd.x <- lapply(1 : NFOLDS, FUN = function(i) run.logreg(dat.x, i, "RBD"))
logreg.spk.x <- lapply(1 : NFOLDS, FUN = function(i) run.logreg(dat.x, i, "SPIKE"))

#SD
sd3.rbd <- lapply(1 : NFOLDS, FUN = function(i) run.sd(dat, i, "RBD", 3))
sd3.spk <- lapply(1 : NFOLDS, FUN = function(i) run.sd(dat, i, "SPIKE", 3))
sd6.rbd <- lapply(1 : NFOLDS, FUN = function(i) run.sd(dat, i, "RBD", 6))
sd6.spk <- lapply(1 : NFOLDS, FUN = function(i) run.sd(dat, i, "SPIKE", 6))

##SVM: LINEAR KERNEL
ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 5, classProbs = TRUE, allowParallel = TRUE)
grd <- data.frame(C = c(0.001, 0.01, 0.5, 1, 2, 5, 10))

svm.rbd <- mclapply(1 : NFOLDS, FUN = function(i) run.mod(dat, i, method = "svmLinear", ctrl = ctrl, grd = grd, prot = "RBD"), mc.cores = NFOLDS)
svm.spk <- mclapply(1 : NFOLDS, FUN = function(i) run.mod(dat, i, method = "svmLinear", ctrl = ctrl, grd = grd, prot = "SPIKE"), mc.cores = NFOLDS)
svm.rbd.x <- mclapply(1 : NFOLDS, FUN = function(i) run.mod(dat.x, i, method = "svmLinear", ctrl = ctrl, grd = grd, prot = "RBD"), mc.cores = NFOLDS)
svm.spk.x <- mclapply(1 : NFOLDS, FUN = function(i) run.mod(dat.x, i, method = "svmLinear", ctrl = ctrl, grd = grd, prot = "SPIKE"), mc.cores = NFOLDS)


##LDA
lda.rbd <- mclapply(1 : NFOLDS, FUN = function(i) run.mod(dat, i, method = "lda", prot = "RBD"), mc.cores = NFOLDS)
lda.spk <- mclapply(1 : NFOLDS, FUN = function(i) run.mod(dat, i, method = "lda", prot = "SPIKE"), mc.cores = NFOLDS)
lda.rbd.x <- mclapply(1 : NFOLDS, FUN = function(i) run.mod(dat.x, i, method = "lda", prot = "RBD"), mc.cores = NFOLDS)
lda.spk.x <- mclapply(1 : NFOLDS, FUN = function(i) run.mod(dat.x, i, method = "lda", prot = "SPIKE"), mc.cores = NFOLDS)



#----------------------------------------------------------------------------------------------------------
#SAVING AND ASSEMBLING

save(logreg.rbd, logreg.spk, logreg.rbd.x, logreg.spk.x, file = file.path(d, "logreg_univar_output.RData"))
save(svm.rbd, svm.spk, svm.rbd.x, svm.spk.x, file = file.path(d, "svm_lin_univar_output.RData"))
save(lda.rbd, lda.spk, lda.rbd.x, lda.spk.x, file = file.path(d, "lda_univar_output.RData"))

##assembly

foo.rbd <- function(i) {
  data.table(pr.logreg = logreg.rbd[[i]]$pred$prob.covid, pred.log = logreg.rbd[[i]]$pred$pred,
                          pr.svm = svm.rbd[[i]]$pred$prob.covid, pred.svm = svm.rbd[[i]]$pred$pred,
                          pr.lda = lda.rbd[[i]]$pred$prob.covid, pred.lda = lda.rbd[[i]]$pred$pred,
                          pred.3sd = sd3.rbd[[i]]$pred, pred.6sd = sd6.rbd[[i]]$pred)
}

foo.spk <- function(i) {
  data.table(pr.logreg = logreg.spk[[i]]$pred$prob.covid, pred.log = logreg.spk[[i]]$pred$pred,
                          pr.svm = svm.spk[[i]]$pred$prob.covid, pred.svm = svm.spk[[i]]$pred$pred,
                          pr.lda = lda.spk[[i]]$pred$prob.covid, pred.lda = lda.spk[[i]]$pred$pred,
                          pred.3sd = sd3.spk[[i]]$pred, pred.6sd = sd6.spk[[i]]$pred)
}

foo.rbd.x <- function(i) {
  data.table(pr.logreg = logreg.rbd.x[[i]]$pred$prob.covid, pred.log = logreg.rbd.x[[i]]$pred$pred,
                          pr.svm = svm.rbd.x[[i]]$pred$prob.covid, pred.svm = svm.rbd.x[[i]]$pred$pred,
                          pr.lda = lda.rbd.x[[i]]$pred$prob.covid, pred.lda = lda.rbd.x[[i]]$pred$pred)
                        }

foo.spk.x <- function(i) {
  data.table(pr.logreg = logreg.spk.x[[i]]$pred$prob.covid, pred.log = logreg.spk.x[[i]]$pred$pred,
                          pr.svm = svm.spk.x[[i]]$pred$prob.covid, pred.svm = svm.spk.x[[i]]$pred$pred,
                          pr.lda = lda.spk.x[[i]]$pred$prob.covid, pred.lda = lda.spk.x[[i]]$pred$pred)
}


res.rbd <- lapply(1 : NFOLDS, FUN = function(i) data.table(dat, foo.rbd(i)))
res.spk <- lapply(1 : NFOLDS, FUN = function(i) data.table(dat, foo.spk(i)))

res.rbd.x <- lapply(1 : NFOLDS, FUN = function(i) data.table(dat.x, foo.rbd.x(i)))
res.spk.x <- lapply(1 : NFOLDS, FUN = function(i) data.table(dat.x, foo.spk.x(i)))

res.rbd <- lapply(res.rbd, FUN = function(z) {
  z[, pr.ens := 1/2 * (pr.svm + pr.lda)]
  z[, pred.ens := ifelse(pr.ens > 0.5, "COVID", "Historical.controls")]

  z[, pr.ens.log := 1/2 * (pr.svm + pr.logreg)]
  z[, pred.ens.log := ifelse(pr.ens.log > 0.5, "COVID", "Historical.controls")]
})
res.rbd.x <- lapply(res.rbd.x, FUN = function(z) {
  z[, pr.ens := 1/2 * (pr.svm + pr.lda)]
  z[, pred.ens := ifelse(pr.ens > 0.5, "COVID", "Historical.controls")]

  z[, pr.ens.log := 1/2 * (pr.svm + pr.logreg)]
  z[, pred.ens.log := ifelse(pr.ens.log > 0.5, "COVID", "Historical.controls")]
})
res.spk <- lapply(res.spk, FUN = function(z) {
  z[, pr.ens := 1/2 * (pr.svm + pr.lda)]
  z[, pred.ens := ifelse(pr.ens > 0.5, "COVID", "Historical.controls")]

  z[, pr.ens.log := 1/2 * (pr.svm + pr.logreg)]
  z[, pred.ens.log := ifelse(pr.ens.log > 0.5, "COVID", "Historical.controls")]
})
res.spk.x <- lapply(res.spk.x, FUN = function(z) {
  z[, pr.ens := 1/2 * (pr.svm + pr.lda)]
  z[, pred.ens := ifelse(pr.ens > 0.5, "COVID", "Historical.controls")]

  z[, pr.ens.log := 1/2 * (pr.svm + pr.logreg)]
  z[, pred.ens.log := ifelse(pr.ens.log > 0.5, "COVID", "Historical.controls")]
})



save(res.rbd, res.spk, res.rbd.x, res.spk.x, file = "~/rds/rds-cew54-wallace-share/Data/elisa/univar_results.RData")




























##
