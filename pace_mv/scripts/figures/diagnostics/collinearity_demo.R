## collinearity_demo.R — visualise E^tech vs E^bio collinearity for BC (Stromal->Tumour).
## Panel A: density kernels K_tech[5um] vs K_bio[30um]; B: expression kernels (EPCAM);
## C: cor vs h_bio bandwidth sweep (density + EPCAM), with current/peak/separability marked.
suppressPackageStartupMessages({library(dbscan); library(ggplot2); library(patchwork)})
setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv")
b<-readRDS("data/breast_cancer/Y_df_for_mcsd.rds"); d<-b$df; Y<-b$Y; d$celltype<-as.character(d$celltype)
HBIO<-c(5,8,12,16,20,30,40,50,75,100,150); HT<-5; EPS<-250
Kt<-Et<-numeric(0); Kb<-Eb<-lapply(HBIO,function(x)numeric(0))
for(s in unique(d$slide)){
  si<-which(d$slide==s); if(length(si)<50)next
  ds<-d[si,]; fr<-dbscan::frNN(as.matrix(ds[,c("x","y")]),eps=EPS); isT<-ds$celltype=="Tumour"; yE<-Y[si,"EPCAM"]
  for(i in which(ds$celltype=="Stromal")){
    id<-fr$id[[i]]; if(!length(id))next; di<-fr$dist[[i]]; m<-isT[id]; if(!any(m))next
    dd<-di[m]; ye<-yE[id[m]]; wt<-exp(-dd/HT); Kt<-c(Kt,sum(wt)); Et<-c(Et,sum(wt*ye))
    for(k in seq_along(HBIO)){wb<-exp(-dd^2/HBIO[k]^2); Kb[[k]]<-c(Kb[[k]],sum(wb)); Eb[[k]]<-c(Eb[[k]],sum(wb*ye))}
  }
}
i30<-which(HBIO==30)
cordK<-sapply(Kb,function(x)cor(Kt,x)); cordE<-sapply(Eb,function(x)cor(Et,x))
TB<-"#4E79A7"; OR<-"#F28E2B"; RD<-"#E15759"; GY<-"grey50"
thm<-theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=11))

dfA<-data.frame(x=Kt,y=Kb[[i30]])
pA<-ggplot(dfA,aes(x,y))+geom_hex(bins=50)+scale_fill_gradient(low="grey88",high=TB,guide="none")+
  geom_smooth(method="lm",se=FALSE,color=RD,linewidth=.6,linetype="dashed")+
  labs(title="A  Density kernels",x="K_tech  (exp(-d/5))",y="K_bio  (exp(-d^2/30^2))")+
  annotate("text",x=min(Kt),y=max(dfA$y),hjust=0,vjust=1,label=sprintf("r = %.2f",cor(Kt,Kb[[i30]])),fontface="bold")+thm

dfB<-data.frame(x=Et,y=Eb[[i30]])
pB<-ggplot(dfB,aes(x,y))+geom_hex(bins=50)+scale_fill_gradient(low="grey88",high=OR,guide="none")+
  geom_smooth(method="lm",se=FALSE,color=RD,linewidth=.6,linetype="dashed")+
  labs(title="B  Expression kernels (EPCAM, tumour->stroma)",x="E_tech  (exp(-d/5))",y="E_bio  (exp(-d^2/30^2))")+
  annotate("text",x=min(Et),y=max(dfB$y),hjust=0,vjust=1,label=sprintf("r = %.2f",cor(Et,Eb[[i30]])),fontface="bold")+thm

dfC<-rbind(data.frame(h=HBIO,cor=cordK,kind="density (K)"),data.frame(h=HBIO,cor=cordE,kind="EPCAM expr (E)"))
pC<-ggplot(dfC,aes(h,cor,color=kind))+
  geom_hline(yintercept=0.5,linetype="dotted",color=GY)+
  geom_vline(xintercept=30,linetype="dashed",color=RD,linewidth=.5)+
  geom_vline(xintercept=12,linetype="dotted",color=GY,linewidth=.4)+
  geom_line(linewidth=.8)+geom_point(size=1.6)+
  scale_color_manual(values=c("density (K)"=TB,"EPCAM expr (E)"=OR),name=NULL)+
  scale_x_continuous(breaks=c(5,12,30,50,75,100,150))+ylim(0,1)+
  annotate("text",x=30,y=0.04,label="current\nh_bio=30",color=RD,size=3,hjust=-0.05)+
  annotate("text",x=12,y=0.97,label="peak\n~12um",color=GY,size=3,hjust=-0.1)+
  labs(title="C  Collinearity vs bio-kernel bandwidth (h_tech=5 fixed)",x="h_bio (um)",y="cor(E_tech, E_bio)")+
  thm+theme(legend.position=c(.98,.98),legend.justification=c(1,1),legend.background=element_rect(fill=alpha("white",.7)))

fig<-(pA|pB)/pC + plot_annotation(title=sprintf("E_tech vs E_bio collinearity  (BC, Stromal focal w/ Tumour neighbours, n=%d cells)",length(Kt)),
       theme=theme(plot.title=element_text(face="bold",size=12)))
ggsave("plots/diagnostics/collinearity_demo.pdf",fig,width=9,height=7.5,device=cairo_pdf)
cat(sprintf("saved plots/diagnostics/collinearity_demo.pdf  (n=%d; current h=30: K r=%.2f, E r=%.2f)\n",length(Kt),cordK[i30],cordE[i30]))
