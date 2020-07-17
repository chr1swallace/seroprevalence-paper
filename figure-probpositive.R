source("common.R")
library(Hmisc) # binconf
(load(file=file_rdata_v5)) # data

################################################################################

#predictions
x <- file.path(d,"predictions_on_test_set_no_12.rds")  %>% readRDS()
with(x,table(group,type))
x[,type:=sub("\\."," ",type)]
x[,Sample.ID:=gsub("wk|WK","Wk",Sample.ID)] # typo
x[grep("Wk",Sample.ID), Week:=sub(".*(Wk[0-9]+).*","\\1",Sample.ID) ]
x[,prob.ENS:=(prob.SVM+prob.LDA)/2]

spike.3sd <- with(m[type=="Historical controls"], mean(SPIKE) + 3*sd(SPIKE))
rbd.3sd <- with(m[type=="Historical controls"], mean(RBD) + 3*sd(RBD))
x[,prob.3sd:=as.numeric(SPIKE > spike.3sd & RBD > rbd.3sd)]
m[,prob.3sd:=as.numeric(SPIKE > spike.3sd & RBD > rbd.3sd)]

spike.6sd <- with(m[type=="Historical controls"], mean(SPIKE) + 6*sd(SPIKE))
rbd.6sd <- with(m[type=="Historical controls"], mean(RBD) + 6*sd(RBD))
x[,prob.6sd:=as.numeric(SPIKE > spike.6sd & RBD > rbd.6sd)]
m[,prob.6sd:=as.numeric(SPIKE > spike.6sd & RBD > rbd.6sd)]

mx <- m[type %in% c("Historical controls","COVID") & !(Sample.ID %in% 1:2)]

x[,Prob:=prob.ENS]
mx[,Label:=sub("Historical controls","Hist. ctrl",type)]

library(cowplot)
theme_set(theme_cowplot(font_size=8))

library(viridis)
p <- ggplot(mapping=aes(y=SPIKE,x=RBD)) +
  theme(legend.position="bottom") +
  scale_x_log10() +
  scale_y_log10() +
  background_grid() +
  geom_vline(xintercept=rbd.3sd,col="black",linetype="solid") +
  geom_hline(yintercept=spike.3sd,col="black",linetype="solid") +
  geom_vline(xintercept=rbd.6sd,col="black",linetype="dashed") +
  geom_hline(yintercept=spike.6sd,col="black",linetype="dashed") +
  ## scale_colour_viridis() +
  scale_colour_gradient2(#"Prob",
    low=scales::muted("blue"),high=scales::muted("red"),
    midpoint=0.5,mid="orange",
    limits=c(0,1),breaks=c(0,1),labels=c(0,1),
    guide=guide_colorbar(labels=NULL,draw.ulim=TRUE,draw.llim=TRUE)) +
  scale_fill_gradient2(#"Prob",
    low=scales::muted("blue"),high=scales::muted("red"),
    midpoint=0.5,mid="orange",
    limits=c(0,1),breaks=c(0,1),labels=c(0,1),
    guide=guide_colorbar(labels=NULL,draw.ulim=TRUE,draw.llim=TRUE)) +
  scale_shape_manual(values=c("Hist. ctrl"=21,
                              "COVID"=19)) +
  draw_text(x=1,y=spike.3sd,text="3SD",vjust=1.2,size=8) +
  draw_text(x=1,y=spike.6sd,text="6SD",vjust=1.2,size=8)

library(ggrepel)
## p3 <- p +
##   geom_point(aes(pch=Label),data=mx,col="grey",size=3) +
##   geom_point(aes(fill=Prob),data=x,pch=23,size=3,alpha=0.7) +
##   ggtitle("Overlay")
x[Prob>0.4 & Prob < 0.6]
x[Prob>0.9 & SPIKE <0.31]
x[Prob<0.1 & RBD > 0.1 & SPIKE <0.07]

if("lab" %in% names(x))
  x[,lab:=NULL]
x[Sample.ID %in% c("P Wk25 25","Wk20 33","Wk24 45") |
 SPIKE==min(x$SPIKE) | RBD==max(x$RBD),lab:=paste0("Pr=",signif(Prob,3))]
p2 <- p +
  geom_point(aes(fill=Prob),data=x,pch=23,size=3,alpha=0.7) +
  geom_label_repel(aes(label=lab),data=x,size=2,min.segment.length=0) + #nudge_x=-0.05,nudge_y=-0.1) +
  ggtitle("Test data")
p1 <- p +
  geom_point(aes(pch=Label),data=mx,col="grey",size=3)+
  ggtitle("Training data")
ens <- plot_grid(p1,p2#,p3
                ,nrow=1,align="h",axis="h")
ens

head(x)
w <- 1.1
ggsave("~/ens-probabilities.png",ens,height=4*w,width=8*w)

xx <- x[,.(group,type,Sample.ID,SPIKE,RBD,prob.ENS,prob.3sd,prob.6sd)]
fwrite(xx, file="~/ens-probabilities.csv")
