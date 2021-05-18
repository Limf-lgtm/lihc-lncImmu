########## inflammatory ##########

setwd("E:/jt/0_other/LIHC_immlnc")

library(fgsea)

GO_all <- gmtPathways("5_inflam/c5.bp.v7.0.symbols.gmt")
GO_all_names <- names(GO_all)
GO_inflamm <- grep("INFLAMM",GO_all_names,value = T)
GO_inflamm_pathway <- GO_all[GO_inflamm]

######## cor phyper #######
N=19406
func_phper <- function(v1,v2){
  m <- length(v1)
  k <- length(v2)
  n <- N-m
  x <- length(intersect(v1,v2))
  p <- ifelse(x>=3,1-phyper(x-1,m,n,k),10)
}

exp_mRNA_cancer_no0 <- as.matrix(read.table("0_data/TCGA/0_deal/FPKM_mRNA_cancer_no0.txt",header=T,row.names = 1,sep = "\t",check.names = F))
exp_lncRNA_cancer_no0 <- as.matrix(read.table("0_data/TCGA/0_deal/FPKM_lncRNA_cancer_no0.txt",header=T,row.names = 1,sep = "\t",check.names = F))
mRNA_sort <- rownames(exp_mRNA_cancer_no0)
lncRNA_sort <- rownames(exp_lncRNA_cancer_no0)
rm(exp_mRNA_cancer_no0);rm(exp_lncRNA_cancer_no0)

con1 <- file("1_GSEA/cor_res/p_all_cancer.txt","r")
con2 <- file("1_GSEA/cor_res/r_all_cancer.txt","r")
pairs_phyper_all <- c()
i <- 0
time1= Sys.time()
while(1){
  oneline_p = readLines(con1, n = 1)
  oneline_r = readLines(con2, n = 1)
  if((length(oneline_p)+length(oneline_r)) == 0){
    break
  }
  i <- i+1
  lnc_i <- lncRNA_sort[i]
  p_i <- as.numeric(unlist(strsplit(oneline_p,"\t")))
  r_i <- as.numeric(unlist(strsplit(oneline_r,"\t")))
  ind_sig <- which((p_i<0.05)&(abs(r_i)>0.3))
  if(length(ind_sig)>=3){
    mRNA_corSig <- mRNA_sort[ind_sig]
    p_lnci <- sapply(GO_inflamm_pathway,func_phper,mRNA_corSig)
    pairs_phyper_all <- rbind(pairs_phyper_all,c(lnc_i,p_lnci))
  }
}
close(con1);close(con2)
time2= Sys.time()
time2-time1

lncs_sig_num <- apply(pairs_phyper_all, 1, function(v){
  sum_sig <- sum(as.numeric(v[2:length(v)])<0.05)
  return(sum_sig)
})
lncs_sigInflam_1 <- pairs_phyper_all[which(lncs_sig_num>=1),1]
intersect(lncs_sigInflam_1,"WWTR1-AS1")

setwd("E:/jt/0_other/LIHC_immlnc")
save.image("5_inflam/all.RData")
load("all.RData")
