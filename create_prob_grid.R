source("common.R")
library(ggplot2)
library(cowplot)
library(kernlab)
library(caret)
library(viridis)
library(ggnewscale)

# d <- "/rds/user/ng414/hpc-work/elisa/results"

(load(file_rdata)) #m
m[, type:=make.names(type)]

dat <- m[!(type %in% c("No.sample.control", "Patient.4")),]
train.ind <- dat$type %in% c("COVID", "Historical.controls")
dat[!train.ind, week := substr(gsub("^[P|p] ", "", Sample.ID), 1, 4)]

rbd <- dat[type == "Historical.controls"]$RBD
spk <- dat[type == "Historical.controls"]$SPIKE

thr.rbd3 <- mean(rbd) + 3 * sd(rbd)
thr.rbd6 <- mean(rbd) + 6 * sd(rbd)
thr.spk3 <- mean(spk) + 3 * sd(spk)
thr.spk6 <- mean(spk) + 6 * sd(spk)

spike.seq <- seq(min(log(dat$SPIKE)), max(log(dat$SPIKE)), len = 20)
rbd.seq <- seq(min(log(dat$RBD)), max(log(dat$RBD)), len = 20)

grd <- expand.grid(SPIKE = exp(spike.seq), RBD = exp(rbd.seq))


(load("~/rds/rds-cew54-wallace-share/Data/elisa/no_12_models.RData")) #lda.mod, svm.mod

prob.grd <- lapply(1 : nrow(grd), FUN = function(i) {
  message(i)
  dt=data.table(LOG = 1- predict(log.mod, grd[i, ], type = "response"), # control = 1, covid = 0
                LDA = predict(lda.mod, grd[i, ], type = "prob")$COVID,
                SVM = predict(svm.mod, grd[i, ], type = "prob")$COVID,
                SVM2=predict(svm2.mod, grd[i, ], type = "prob")$COVID)
  dt[,ENS := (LDA + SVM)/2]
  dt[,ENS2 := (LDA + SVM2)/2]
  dt
}) %>% do.call(rbind, .) %>% data.table(., grd)

prob.grd <- melt(prob.grd, c("SPIKE", "RBD"))
setnames(prob.grd, c("variable", "value"), c("method", "covid.prob"))
saveRDS(prob.grd,file.path(d, "probability_grid_no12.rds"))
