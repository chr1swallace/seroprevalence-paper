source("common.R")
library(ggplot2)
library(cowplot)
library(kernlab)
library(caret)
library(parallel)

# d <- "/rds/user/ng414/hpc-work/elisa/results"

(load(file_rdata)) #m
m[, type:=make.names(type)]
#type: historical controls = controls, covid = cases.
f <- 9

#only keep historical and COVID
dat <- m[!is.na(fold1)]
table(m$type)
dat.x <- dat[!(Sample.ID %in% c("1", "2"))]
table(m$type)

g1 <- ggplot(data = dat, aes(x = SPIKE, y = RBD, colour = type)) + geom_point(alpha = 0.5) + geom_point(data = m[Sample.ID %in% c("1", "2", "10", "Sp2019 158", "Sp2019 3")], aes(x = SPIKE, y = RBD), shape = 8, size = 12)
g2 <- ggplot(data = dat, aes(x = log(SPIKE), y = log(RBD), colour = type)) + geom_point(alpha = 0.5) + geom_point(data = m[Sample.ID %in% c("1", "2", "10", "Sp2019 158", "Sp2019 3")], aes(x = log(SPIKE), y = log(RBD)), shape = 8, size = 12)
plot_grid(g1, g2)

# dat[, c("SPIKE", "RBD") := list(log(SPIKE), log(RBD))]

#----------------------------------------------------------------------------------------------------------
#LOGISTIC REGRESSION
#logistic: 0 = COVID, 1 = historical
run.logreg <- function(dat, k, log_ = TRUE) {
  fold <- paste0("fold", k)
  out <- data.table(dat[,c("Sample.ID", "type")], fold = dat[, get(fold)], method = "LOG", fold.n = k, prob.covid = vector("numeric", nrow(dat)), pred = vector("character", nrow(dat)))

  mod.output <- lapply(1 : 10, FUN = function(i) {
    message(sprintf("fold%s: %s", k, i))
    if(log_) {
    mod <- glm(data = dat[get(fold) != i,], factor(type) ~ log(SPIKE) + log(RBD), family = binomial, maxit = 100)
  } else {
        mod <- glm(data = dat[get(fold) != i,], factor(type) ~ SPIKE + RBD, family = binomial, maxit = 100)
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

logreg <- lapply(1 : f, FUN = function(i) run.logreg(dat, i))
for(i in 1 : f) {
  message("i\n")
  print(table(logreg[[i]]$pred$type, logreg[[i]]$pred$pred))
  print(logreg[[i]]$pred[type != pred])
}

#----
logreg.x <- lapply(1 : f, FUN = function(i) run.logreg(dat.x, i))
for(i in 1 : f) {
  message("i\n")
  print(table(logreg.x[[i]]$pred$type, logreg.x[[i]]$pred$pred))
  print(logreg.x[[i]]$pred[type != pred])
}
#----------------------------------------------------------------------------------------------------------
#ctrl = options for parameter tuning, leave as NULL if don't want to tune
#grd = parameter grid for tuning, leave as NULL if don't want to tune
#log_ = log SPIKE and RDS?
#k = which fold set?
run.mod <- function(dat, k, method, grd = NULL, ctrl = NULL, log_ = TRUE) {
  fold <- paste0("fold", k)
  out <- data.table(dat[,c("Sample.ID", "type")], fold = dat[,get(fold)], method = method, fold.n = k, prob.covid = vector("numeric", nrow(dat)), pred = vector("character", nrow(dat)))

  mod.output <- lapply(1 : 10, FUN = function(i) {
    message(sprintf("%s: fold%s %s", method, k, i))
    dat.train <- dat[get(fold) != i,]
    dat.test <- dat[get(fold) == i,]

    if(all(is.null(grd), is.null(ctrl))) {
      if(log_) {
        mod <- train(data = dat.train, factor(type) ~ log(SPIKE) + log(RBD), method = method) } else {
                mod <- train(data = dat.train, factor(type) ~ SPIKE + RBD, method = method)
        }
    } else {

      if(log_) {
        mod <- train(data = dat.train, factor(type) ~ log(SPIKE) + log(RBD), method = method, trControl = ctrl, tuneGrid = grd) } else {
                mod <- train(data = dat.train, factor(type) ~ SPIKE + RBD, method = method, trControl = ctrl, tuneGrid = grd)
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
#SVM
##LINEAR KERNEL
ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 5, classProbs = TRUE, allowParallel = TRUE)
grd <- data.frame(C = c(0.001, 0.01, 0.5, 1, 2, 5, 10))

svm.lin <- mclapply(1 : f, FUN = function(i) run.mod(dat, i, method = "svmLinear", ctrl = ctrl, grd = grd), mc.cores = f)
for(i in 1 : f) {
  message("i\n")
  print(table(svm.lin[[i]]$pred$type, svm.lin[[i]]$pred$pred))
  print(svm.lin[[i]]$pred[type != pred])
}

#----
svm.lin.x <- mclapply(1 : f, FUN = function(i) run.mod(dat.x, i, method = "svmLinear", ctrl = ctrl, grd = grd), mc.cores = f)
for(i in 1 : f) {
  message("i\n")
  print(table(svm.lin.x[[i]]$pred$type, svm.lin.x[[i]]$pred$pred))
  print(svm.lin.x[[i]]$pred[type != pred])
}

##POLYNIMIAL (QUADRATIC) KERNEL
ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 5, classProbs = TRUE, allowParallel = TRUE)
grd <- expand.grid(degree = 2, scale = c(0.1, 1, 2), C = c(0.001, 0.01, 0.5, 1, 2, 5, 10))

svm.quad <- mclapply(1 : f, FUN = function(i) run.mod(dat, i, method = "svmPoly", ctrl = ctrl, grd = grd), mc.cores = f)
for(i in 1 : f) {
  message("i\n")
  print(table(svm.quad[[i]]$pred$type, svm.quad[[i]]$pred$pred))
  print(svm.quad[[i]]$pred[type != pred])
}

#----
svm.quad.x <- mclapply(1 : f, FUN = function(i) run.mod(dat.x, i, method = "svmPoly", ctrl = ctrl, grd = grd), mc.cores = f)
for(i in 1 : f) {
  message("i\n")
  print(table(svm.quad.x[[i]]$pred$type, svm.quad.x[[i]]$pred$pred))
  print(svm.quad.x[[i]]$pred[type != pred])
}

#----------------------------------------------------------------------------------------------------------
#LINEAR DISCRIMINANT ANALYSIS

lda <- mclapply(1 : f, FUN = function(i) run.mod(dat, i, method = "lda"), mc.cores = f)
for(i in 1 : f) {
  message("i\n")
  print(table(lda[[i]]$pred$type, lda[[i]]$pred$pred))
  print(lda[[i]]$pred[type != pred])
}

#----
lda.x <- mclapply(1 : f, FUN = function(i) run.mod(dat.x, i, method = "lda"), mc.cores = f)
for(i in 1 : f) {
  message("i\n")
  print(table(lda.x[[i]]$pred$type, lda.x[[i]]$pred$pred))
  print(lda.x[[i]]$pred[type != pred])
}


#----------------------------------------------------------------------------------------------------------
#SAVING AND ASSEMBLING

save(logreg, logreg.x, file = file.path(d, "logreg_output.RData"))
save(svm.lin, svm.lin.x, file = file.path(d, "svm_lin_output.RData"))
save(svm.quad, svm.quad.x, file = file.path(d, "svm_quad_output.RData"))
save(lda, lda.x, file = file.path(d, "lda_output.RData"))

##assembly

foo <- function(i) {
  data.table(pr.logreg = logreg[[i]]$pred$prob.covid, pred.logreg = logreg[[i]]$pred$pred,
                          pr.svmlin = svm.lin[[i]]$pred$prob.covid, pred.svmlin = svm.lin[[i]]$pred$pred,
                          pr.svm2 = svm.quad[[i]]$pred$prob.covid, pred.svm2 = svm.quad[[i]]$pred$pred,
                          pr.lda = lda[[i]]$pred$prob.covid, pred.lda = lda[[i]]$pred$pred)
}

foo.x <- function(i) {
  data.table(pr.logreg = logreg.x[[i]]$pred$prob.covid, pred.logreg = logreg.x[[i]]$pred$pred,
                          pr.svmlin = svm.lin.x[[i]]$pred$prob.covid, pred.svmlin = svm.lin.x[[i]]$pred$pred,
                          pr.svm2 = svm.quad.x[[i]]$pred$prob.covid, pred.svm2 = svm.quad.x[[i]]$pred$pred,
                          pr.lda = lda.x[[i]]$pred$prob.covid, pred.lda = lda.x[[i]]$pred$pred)
}

res <- lapply(1 : f, FUN = function(i) data.table(dat, foo(i)))
res.x <- lapply(1 : f, FUN = function(i) data.table(dat.x, foo.x(i)))

#----------------------------------------------------------------------------------------------------------
#ENSEMBLE: SVM + LDA

res <- lapply(res, FUN = function(z) {
  z[, pr.ens := 1/2 * (pr.svmlin + pr.lda)]
  z[, pred.ens := ifelse(pr.ens > 0.5, "COVID", "Historical.controls")]

  z[, pr.ens.log := 1/2 * (pr.svmlin + pr.logreg)]
  z[, pred.ens.log := ifelse(pr.ens.log > 0.5, "COVID", "Historical.controls")]
})
res.x <- lapply(res.x, FUN = function(z) {
  z[, pr.ens := 1/2 * (pr.svmlin + pr.lda)]
  z[, pred.ens := ifelse(pr.ens > 0.5, "COVID", "Historical.controls")]

  z[, pr.ens.log := 1/2 * (pr.svmlin + pr.logreg)]
  z[, pred.ens.log := ifelse(pr.ens.log > 0.5, "COVID", "Historical.controls")]
})

save(res, res.x, file = "~/rds/rds-cew54-wallace-share/Data/elisa/ml_results.RData")

 # q('no')
#----------------------------------------------------------------------------------------------------------

#
# ###analysis
# ind <- which(dat$Sample.ID %in% c("1", "2", "10", "Sp2019 158"))[-1] #double Sp2019 158, only need the second one
#
# ###fixing folds
# (load(file_rdata)) #m
# m <- m[,!grepl("fold", colnames(m)), with = FALSE]
#
# X <- m[!duplicated(Sample.ID), c('Sample.ID', 'type', 'group', 'rowid'), with = FALSE]
# set.seed(42)
# X[type %in% c("Historical controls","COVID"),fold1:=1:.N %% 10 + 1,by="type"]
# set.seed(42)
# X[type %in% c("Historical controls","COVID"),fold2:=1:.N %% 10 + 1,by="type"]
# set.seed(42)
# X[type %in% c("Historical controls","COVID"),fold3:=1:.N %% 10 + 1,by="type"]
# m <- data.table(m, X[match(m$Sample.ID, Sample.ID), c("fold1", "fold2", "fold3")])
#
#
# #check
# mm <- m[type %in% c("Historical controls","COVID")]
# tapply(mz$Sample.ID, mz$fold1, FUN = function(x) sum(table(x) == 2)) %>% sum
# # [1] 295
# table(table(mm$Sample.ID))
# #   1   2
# # 131 295
#
# mz <- mz[order('Sample.ID', 'type', 'group', 'rowid')]
# m <- m[order('Sample.ID', 'type', 'group', 'rowid')]
# all.equal(mz[, 1 : 8], m)
#

































##
