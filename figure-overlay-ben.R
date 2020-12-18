source("common.R")
library(Hmisc) # binconf
(load(file=file_rdata))

################################################################################

## from https://www.tandfonline.com/doi/abs/10.1080/00031305.2018.1473796?af=R&journalCode=utas20
## functions for upper and lower limits of MI-Wilson interval
#Upper limit function for the MI-Wilson interval for multiple imputation
upper.limit= function(z, n, r.m, mean_q) {
  ul = ((((2*mean_q) + ((z^2)/n) + (((z^2) * r.m)/n))/
           (2*(1 + ((z^2)/n) + (((z^2)*r.m)/n)))) +
          sqrt(((((2*mean_q) + ((z^2)/n) + (((z^2)*r.m)/n))^2)/
                  (4*(1 + ((z^2)/n) + (((z^2)*r.m)/n))^2)) -
                 ((mean_q^2)/(1 + ((z^2)/n) + (((z^2)*r.m)/n)))))
  return(ul)
}

#Lower limit function for the MI-Wilson interval for multiple imputation
lower.limit= function(z, n, r.m, mean_q) {
  ll = ((((2*mean_q) + ((z^2)/n) + (((z^2) * r.m)/n))/
           (2*(1 + ((z^2)/n) + (((z^2)*r.m)/n)))) -
          sqrt(((((2*mean_q) + ((z^2)/n) + (((z^2)*r.m)/n))^2)/
                  (4*(1 + ((z^2)/n) + (((z^2)*r.m)/n))^2)) -
                 ((mean_q^2)/(1 + ((z^2)/n) + (((z^2)*r.m)/n)))))
  return(ll)
}

x <- file.path(d,"predictions_on_test_set_no_12.rds")  %>% readRDS()
with(x,table(group,type))
x[,type:=sub("\\."," ",type)]
x[,Sample.ID:=gsub("wk|WK","Wk",Sample.ID)] # typo
x[grep("Wk",Sample.ID), Week:=sub(".*(Wk[0-9]+).*","\\1",Sample.ID) ]
x[,prob.ENS:=(prob.SVM+prob.LDA)/2]

################################################################################

## estimate proportion + CI

fitone <- function(x,probvar="prob.log",byvars=c("type","Week"),M=1000) {
  set.seed(42) # reprocibility
  summ <- lapply(1:M, function(j) {
    x[,y:=runif(.N) <= x[[probvar]] ] # simulate case status conditional on fit
    x[,.(yhat=sum(y),phat=mean(y),varp=mean(y)*mean(1-y)/.N,N=.N,j=j),by=byvars]
  })  %>% rbindlist()
  ## wilson's interval with observed data for comparison
  wilsons <- with(summ, binconf(yhat,N))
  summ <- cbind(summ,wilsons)
  summ2 <- summ[,.(N=mean(N),
                   qm=mean(phat),um=mean(varp),bm=var(phat), # MI
                   Est=mean(PointEst),Lower=mean(Lower),Upper=mean(Upper)), # observed data
                by=byvars]
  summ2[,Tm:=(1+1/M)*bm + um]
  if(any(summ2$bm==0)) { # all realisations the same - doesn't happen, just in case
    summ2[bm==0, c("v","rm"):=list(1000000000,0)]
  }
  summ2[bm!=0,
        c("v","rm"):=list((M-1)*(1+(um/((1+ (1/M))*bm)))^2, (1 + (1/M))*(bm/um))]
  summ2[,ll:=lower.limit(z = qt(0.975, df = v), n = N, r.m = rm,
                             mean_q = qm),by=byvars]
  summ2[,ul:=upper.limit(z = qt(0.975, df = v), n = N, r.m = rm,
                             mean_q = qm),by=byvars]
  options(digits=2)
  summ2[,c(byvars,c("N","Est","Lower","Upper","qm","ll","ul")),with=FALSE]
}

m <- "ENS"
## per group and week
probs <- fitone(x,paste0("prob.",m),byvars=c("type","Week"))
probs[order(type,Week),.(type,Week,N,estimate=qm,lower.ci=ll,upper.ci=ul)]

## per week
probs2 <- fitone(x,paste0("prob.",m),byvars="Week")[,type:="Combined"]
probs2[order(Week),.(Week,N,estimate=qm,lower.ci=ll,upper.ci=ul)]

## combine
probs <- rbind(probs,probs2)[,method:=m]
probs[,week_num:=as.numeric(sub("Wk","",Week))]
probs[,type:=factor(type,levels=c("Blood donors","Pregnant volunteers","Combined"))]

## Ben's
b <- fread(file.path(d,"prev_intervals_10000.csv"))
setnames(b,"timepoints","week_num")
b <- melt(b,"week_num")
b[,group:=sub("low|high","",variable)]
b[group!="median",group:=paste0(group,"%")]

pred <- merge(probs[type=="Combined"],b,by="week_num",all=TRUE)

col <- "#5D7A9B"
col <- "dodgerblue"
theme_set(theme_cowplot(font_size=10))
p <- ggplot(pred,aes(x=week_num)) +
  geom_pointrange(aes(y=qm,ymin=ll,ymax=ul)) +
  geom_path(aes(y=value,lty=group,group=variable),data=b,col=col) +
  labs(x="Week",y="Machine learning prediction + 95% conf. int.") +
  scale_linetype_manual("Bayesian posterior",values=c(median="solid","70%"="dashed","95%"="dotted")) +
 background_grid(major="y") +
  theme(legend.title=element_text(colour=col,face="bold"),
        legend.text=element_text(colour=col),
        legend.position=c(0.1,0.9))
p

ggsave("figure-overlay.png",p,height=6,width=6)
