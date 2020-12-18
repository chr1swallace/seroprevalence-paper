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

## add 3sd double positives
spike.3sd <- with(m[type=="Historical controls"], mean(SPIKE) + 3*sd(SPIKE))
rbd.3sd <- with(m[type=="Historical controls"], mean(RBD) + 3*sd(RBD))
x[,prob.3sd:=as.numeric(SPIKE > spike.3sd & RBD > rbd.3sd)]
m[,prob.3sd:=as.numeric(SPIKE > spike.3sd & RBD > rbd.3sd)]

mx <- m[type %in% c("Historical controls","COVID")]
library(viridis)

p3 <- ggplot(mapping=aes(y=log2(SPIKE),x=log2(RBD))) +
  geom_point(aes(pch=type),data=mx,col="grey") +
  geom_point(aes(col=prob.SVM),data=x,pch=17) +
  theme(legend.position="bottom") +
  ## scale_colour_viridis() +
  scale_colour_gradient2(midpoint=0.5,mid="orange") + 
  scale_shape_manual(values=c("Historical controls"=21,"COVID"=20))
p2 <- ggplot(mapping=aes(y=log2(SPIKE),x=log2(RBD))) +
  ## geom_point(aes(pch=type),data=mx,col="grey") +
  geom_point(aes(col=prob.SVM),data=x,pch=17) +
  theme(legend.position="bottom") +
  ## scale_colour_viridis() +
  scale_colour_gradient2(midpoint=0.5,mid="orange") + 
  scale_shape_manual(values=c("Historical controls"=21,"COVID"=20))
p1 <- ggplot(mapping=aes(y=log2(SPIKE),x=log2(RBD))) +
  geom_point(aes(pch=type),data=mx,col="grey") +
  theme(legend.position="bottom") +
  ## scale_colour_viridis() +
  scale_colour_gradient2(midpoint=0.5,mid="orange") + 
  ## geom_point(aes(col=prob.SVM),data=x,pch=4) +
  scale_shape_manual(values=c("Historical controls"=21,"COVID"=20))
svm <- plot_grid(p1,p2,p3,nrow=1)
svm

p3 <- ggplot(mapping=aes(y=log2(SPIKE),x=log2(RBD))) +
  geom_point(aes(pch=type),data=mx,col="grey") +
  geom_point(aes(col=prob.ENS),data=x,pch=17) +
  theme(legend.position="bottom") +
  ## scale_colour_viridis() +
  scale_colour_gradient2(midpoint=0.5,mid="orange") + 
  scale_shape_manual(values=c("Historical controls"=21,"COVID"=20))
p2 <- ggplot(mapping=aes(y=log2(SPIKE),x=log2(RBD))) +
  ## geom_point(aes(pch=type),data=mx,col="grey") +
  geom_point(aes(col=prob.ENS),data=x,pch=17) +
  theme(legend.position="bottom") +
  ## scale_colour_viridis() +
  scale_colour_gradient2(midpoint=0.5,mid="orange") + 
  scale_shape_manual(values=c("Historical controls"=21,"COVID"=20))
p1 <- ggplot(mapping=aes(y=log2(SPIKE),x=log2(RBD))) +
  geom_point(aes(pch=type),data=mx,col="grey") +
  theme(legend.position="bottom") +
  ## scale_colour_viridis() +
  scale_colour_gradient2(midpoint=0.5,mid="orange",
    limits=c(0,1),breaks=c(0,1),labels=c(0,1),
                         guide=guide_colorbar(labels=NULL,draw.ulim=TRUE,draw.llim=TRUE)) +
  ## geom_point(aes(col=prob.ENS),data=x,pch=4) +
  scale_shape_manual(values=c("Historical controls"=21,"COVID"=20))
ens <- plot_grid(p1,p2,p3,nrow=1,align="hv")
ens

p3 <- ggplot(mapping=aes(y=log2(SPIKE),x=log2(RBD))) +
  geom_point(aes(pch=type),data=mx,col="grey") +
  geom_point(aes(col=prob.LOG),data=x,pch=17) +
  theme(legend.position="bottom") +
  ## scale_colour_viridis() +
  scale_colour_gradient2(midpoint=0.5,mid="orange") + 
  scale_shape_manual(values=c("Historical controls"=21,"COVID"=20))
p2 <- ggplot(mapping=aes(y=log2(SPIKE),x=log2(RBD))) +
  ## geom_point(aes(pch=type),data=mx,col="grey") +
  geom_point(aes(col=prob.LOG),data=x,pch=17) +
  theme(legend.position="bottom") +
  ## scale_colour_viridis() +
  scale_colour_gradient2(midpoint=0.5,mid="orange") + 
  scale_shape_manual(values=c("Historical controls"=21,"COVID"=20))
p1 <- ggplot(mapping=aes(y=log2(SPIKE),x=log2(RBD))) +
  geom_point(aes(pch=type),data=mx,col="grey") +
  theme(legend.position="bottom") +
  ## scale_colour_viridis() +
  scale_colour_gradient2(midpoint=0.5,mid="orange") + 
  ## geom_point(aes(col=prob.LOG),data=x,pch=4) +
  scale_shape_manual(values=c("Historical controls"=21,"COVID"=20))
logistic <- plot_grid(p1,p2,p3,nrow=1)


p3 <- ggplot(mapping=aes(y=log2(SPIKE),x=log2(RBD))) +
  geom_point(aes(pch=type),data=mx,col="grey") +
  geom_point(aes(col=prob.3sd),data=x,pch=17) +
  theme(legend.position="bottom") +
  ## scale_colour_viridis() +
  scale_colour_gradient2(midpoint=0.5,mid="orange") + 
  scale_shape_manual(values=c("Historical controls"=21,"COVID"=20)) +
geom_vline(xintercept=log2(rbd.3sd)) +
  geom_hline(yintercept=log2(spike.3sd))
p2 <- ggplot(mapping=aes(y=log2(SPIKE),x=log2(RBD))) +
  ## geom_point(aes(pch=type),data=mx,col="grey") +
  geom_point(aes(col=prob.3sd),data=x,pch=17) +
  theme(legend.position="bottom") +
  ## scale_colour_viridis() +
  scale_colour_gradient2(midpoint=0.5,mid="orange") + 
  scale_shape_manual(values=c("Historical controls"=21,"COVID"=20)) +
  geom_vline(xintercept=log2(rbd.3sd)) +
  geom_hline(yintercept=log2(spike.3sd))
p1 <- ggplot(mapping=aes(y=log2(SPIKE),x=log2(RBD))) +
  geom_point(aes(pch=type),data=mx,col="grey") +
  theme(legend.position="bottom") +
  ## scale_colour_viridis() +
  scale_colour_gradient2(midpoint=0.5,mid="orange") + 
  ## geom_point(aes(col=prob.LOG),data=x,pch=4) +
  scale_shape_manual(values=c("Historical controls"=21,"COVID"=20))
sd3 <- plot_grid(p1,p2,p3,nrow=1)

plot_grid(ens,svm,logistic,sd3,ncol=1)

ggplot(x, aes(x=prob.LOG)) + geom_histogram(col="grey")
ggplot(x, aes(x=prob.ENS)) + geom_histogram(col="grey")
ggplot(x, aes(x=SPIKE,y=RBD,col=prob.LOG)) + geom_point()
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

## per group and week
m <- "ENS"
probs <- fitone(x,paste0("prob.",m),byvars=c("type","Week"))
probs[order(type,Week),.(type,Week,N,estimate=qm,lower.ci=ll,upper.ci=ul)]

## per week
probs2 <- fitone(x,paste0("prob.",m),byvars="Week")[,type:="Combined"]
probs2[order(Week),.(Week,N,estimate=qm,lower.ci=ll,upper.ci=ul)]

## per group
probs3 <- fitone(x,paste0("prob.",m),byvars="type")[,Week:="Combined"]
probs3[order(type),.(type,N,estimate=qm,lower.ci=ll,upper.ci=ul)]

## combine
probs <- rbind(probs,probs2,probs3)[,method:=m]
probs[,week_num:=as.numeric(sub("Wk","",Week))]
probs[,type:=factor(type,levels=c("Blood donors","Pregnant volunteers","Combined"))]

pred <- probs
pred[type=="Combined" & Week=="Wk14",.(Week,method,qm)]
pred[type=="Combined" & Week=="Wk25",.(Week,method,qm)]

tmp <- dcast(pred[type!="Combined" & method=="ENS" & Week!="Wk14"],
      Week ~ type, value.var="Est")
tmp
setnames(tmp,make.names(names(tmp)))
tmp[,.(bd=mean(Blood.donors),pv=mean(Pregnant.volunteers))]


ggplot(pred[type=="Combined"],aes(x=week_num,y=qm,ymin=ll,ymax=ul#,col=type
                                  )) +
  geom_smooth(se=FALSE,
              ## method="lm",formula=y~poly(x,3),
                                        span=0.8,
    alpha=0.2,
              linetype="dotted") +
  geom_pointrange() +
  ## facet_wrap(~type) +
  facet_grid(type ~ method) +
  background_grid(major="y")

options(digits=6)
probs[order(type,Week),.(type,Week,N,estimate=qm,lower.ci=ll,upper.ci=ul)]
fwrite(probs[order(type,Week),.(type,Week,N,estimate=qm,lower.ci=ll,upper.ci=ul)],
       file=file.path(d,"estimates.csv"))
fwrite(probs[order(type,Week),.(type,Week,N,estimate=qm,lower.ci=ll,upper.ci=ul)],
       file=file.path("estimates.csv"))


## compare with sweden deaths from ecdc
library(utils)
#read the Dataset sheet into “R”. The dataset will be called "data".
data <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM")
data%<>%as.data.table

pl=ggplot(pred[type=="Combined"],aes(x=week_num,y=qm,ymin=ll,ymax=ul#,col=type
                                  )) +
  geom_smooth(se=FALSE,
              ## method="lm",formula=y~poly(x,3),
                                        span=0.8,
    alpha=0.2,
              linetype="dotted") +
  geom_pointrange() +
  ## facet_wrap(~type) +
  facet_grid(type ~ method) +
  background_grid(major="y")

ggplot(deaths,aes(x=week_num,y=cumsum(per100k))) +
  geom_smooth(method="lm",formula=y~poly(x,3),) +
  geom_point()

deaths=data[countriesAndTerritories=="Sweden"]
deaths[,week_num:=as.numeric(sub("2020-","",year_week))]
deaths=deaths[order(week_num)]
deaths[,perd:=100*deaths_weekly/popData2019]
deaths[,perc:=cases_weekly/popData2019]
comb=merge(unique(pred[type=="Combined",.(week_num,qm,ll,ul)]),
           deaths[,.(week_num,perd,perc)],
           by="week_num", all=TRUE)
comb=comb[order(week_num)]

comb[,perd.norm:=median(qm,na.rm=TRUE) * perd / median(cumsum(perd),na.rm=TRUE)]

ggplot(comb,aes(x=week_num)) +
  geom_pointrange(aes(y=qm,ymin=ll,ymax=ul),col="darkblue") +
  geom_smooth(aes(y=qm),col="darkblue",se=FALSE,linetype="dotted") +#,linetype=="dashed") +
  ## geom_point(aes(y=cumsum(10*perc)),col="purple") +
  ## geom_path(aes(y=cumsum(10*perc)),col="purple",linetype="dashed") +
  geom_point(aes(x=week_num,y=cumsum(perd.norm)),col="black") +
  geom_path(aes(x=week_num,y=cumsum(perd.norm)),col="black",linetype="dashed")
