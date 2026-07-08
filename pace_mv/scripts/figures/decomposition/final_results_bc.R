## final_results_bc.R — FINAL BC 313 decomposition results (spec^2-weighted ANOVA SS).
suppressPackageStartupMessages({library(Matrix); library(dbscan); library(dplyr)})
say <- function(...) cat(sprintf(...),"\n")
mv<-readRDS("data/breast_cancer/sweeps/mvpql_percell_hc.rds"); fit<-mv$fit; U<-fit$U; seU<-fit$se_U
TYPES<-sort(unique(as.character(mv$cell_meta$celltype))); ct<-as.character(mv$cell_meta$celltype); nC<-mv$cell_meta$nCount
gf<-mv$decomposition$gene_focal_4block; gs<-mv$gene_set
ymv<-readRDS("data/breast_cancer/Y_df_for_mcsd.rds"); df<-ymv$df
xy<-as.matrix(df[,c("x","y")]); n<-nrow(df); H<-30
fr<-dbscan::frNN(xy,eps=90); iv<-rep.int(seq_len(n),lengths(fr$id)); jv<-unlist(fr$id,use.names=FALSE); dv<-unlist(fr$dist,use.names=FALSE)
jc<-match(ct[jv],TYPES); ok<-!is.na(jc); A<-as.matrix(sparseMatrix(i=iv[ok],j=jc[ok],x=exp(-(dv[ok])^2/H^2),dims=c(n,length(TYPES)))); colnames(A)<-TYPES
rm(fr,iv,jv,dv,jc,ok); invisible(gc())
share<-matrix(0,length(TYPES),length(TYPES),dimnames=list(neighbour=TYPES,focal=TYPES)); drivers<-list()
for(F in TYPES){cF<-which(ct==F); spec<-gf$spec[gf$focal==F][match(gs,gf$gene[gf$focal==F])]; spec[is.na(spec)]<-0
  geneV<-numeric(length(gs)); Vpair<-setNames(numeric(length(TYPES)-1),setdiff(TYPES,F))
  for(k in setdiff(TYPES,F)){row<-paste0(F,"::",k); if(!(row%in%rownames(U)))next
    a<-A[cF,k];ex<-a>1e-8;vc<-if(sum(ex)>2)var(a[ex])else 0; contr<-(U[row,]^2+seU[row,]^2)*vc
    geneV<-geneV+contr; Vpair[k]<-sum(spec^2*contr)}
  if(sum(Vpair)>0)share[names(Vpair),F]<-Vpair/sum(Vpair)
  drivers[[F]]<-gs[order(spec^2*geneV,decreasing=TRUE)][1:8]
}
## emit final tables
sh<-as.data.frame(round(100*share,1)); write.csv(sh,"streaming/FINAL_bc313_neighbour_shares.csv")
dv2<-do.call(rbind,lapply(TYPES,function(F)data.frame(focal=F,rank=1:8,gene=drivers[[F]]))); write.csv(dv2,"streaming/FINAL_bc313_driver_genes.csv",row.names=FALSE)

say("================ FINAL BC 313 DECOMPOSITION (spec^2-weighted ANOVA SS) ================")
say("\n--- PAIRWISE: top neighbour shares per focal (columns sum to 100%%) ---")
for(F in TYPES){v<-sort(share[,F],decreasing=TRUE);v<-v[v>0.03][1:min(3,sum(v>0.03))]
  say("  %-15s <- %s", F, paste(sprintf("%s %.0f%%",names(v),100*v),collapse=", "))}
say("\n--- DRIVER PROGRAMS: top genes per focal (spec^2-weighted, contamination removed) ---")
for(F in c("Stromal","Macrophage","Endothelial","Myoepithelial","Tumour","T_Cell"))
  say("  %-15s : %s", F, paste(drivers[[F]],collapse=", "))
say("\nsaved: streaming/FINAL_bc313_neighbour_shares.csv  +  FINAL_bc313_driver_genes.csv")
say("figure: streaming/figures/decomp_anova_ss_weighted_bc.pdf")
