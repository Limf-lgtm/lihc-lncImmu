
setwd("E:/jt/0_other/LIHC_immlnc")

###############  deal exp #################
id2symbol <- function(v,m_list,n_col,new_col){
  m_list <- as.matrix(m_list)
  v_inter <- intersect(v,m_list[,n_col])
  m_list_inter <- m_list[which(m_list[,n_col]%in%v_inter),]
  rownames(m_list_inter) <- m_list_inter[,n_col]
  ind_v <- which(v%in%v_inter)
  v_name <- m_list_inter[v[ind_v],new_col]
  return(list(ind_v,v_name))
}

mRNA_idList <- as.matrix(read.table("0_data/idList/mRNA_unique.txt",header = F,sep = "\t"))
lncRNA_idList <- as.matrix(read.table("0_data/idList/lncRNA_unique.txt",header = F,sep = "\t"))
mRNA_idList[,2] <- substr(mRNA_idList[,2],1,nchar("ENSG00000242268"))
lncRNA_idList[,2] <- substr(lncRNA_idList[,2],1,nchar("ENSG00000242268"))

exp_all <- read.table("0_data/TCGA/HTSeq - FPKM.merge.txt",header = T,sep = "\t",row.names = 1,check.names = F)
genes_all <- rownames(exp_all) ; genes_all_id <- substr(genes_all,1,nchar("ENSG00000242268"))

res_mRNA <- id2symbol(genes_all_id,mRNA_idList,2,1)
exp_mRNA <- exp_all[res_mRNA[[1]],] ; rownames(exp_mRNA) <- res_mRNA[[2]]
res_lncRNA <- id2symbol(genes_all_id,lncRNA_idList,2,1)
exp_lncRNA <- exp_all[res_lncRNA[[1]],] ; rownames(exp_lncRNA) <- res_lncRNA[[2]]
ind_0row_mRNA <- apply(exp_mRNA, 1, function(v){all(v == 0)})
ind_0row_lncRNA <- apply(exp_lncRNA, 1, function(v){all(v == 0)})
exp_mRNA_no0 <- exp_mRNA[!ind_0row_mRNA,] # 19421
exp_lncRNA_no0 <- exp_lncRNA[!ind_0row_lncRNA,] # 13532

samples_cancer <- colnames(exp_all)[(as.numeric(substr(colnames(exp_all),14,15))<=9)]
samples_normal <- colnames(exp_all)[!((as.numeric(substr(colnames(exp_all),14,15))<=9))]
exp_mRNA_cancer <- exp_mRNA_no0[,samples_cancer] ; exp_mRNA_normal <- exp_mRNA_no0[,samples_normal]
# 19421 374 ; 19421 50
exp_lncRNA_cancer <- exp_lncRNA_no0[,samples_cancer] ; exp_lncRNA_normal <- exp_lncRNA_no0[,samples_normal]
# 13532 374 ; 13532 50

write.table(exp_mRNA_cancer,"0_data/TCGA/0_deal/FPKM_mRNA_cancer.txt",quote = F,row.names = T,col.names=T,sep="\t")
write.table(exp_mRNA_normal,"0_data/TCGA/0_deal/FPKM_mRNA_normal.txt",quote = F,row.names = T,col.names=T,sep="\t")
write.table(exp_lncRNA_cancer,"0_data/TCGA/0_deal/FPKM_lncRNA_cancer.txt",quote = F,row.names = T,col.names=T,sep="\t")
write.table(exp_lncRNA_normal,"0_data/TCGA/0_deal/FPKM_lncRNA_normal.txt",quote = F,row.names = T,col.names=T,sep="\t")

############## read data ######################
exp_mRNA_cancer <- as.matrix(read.table("0_data/TCGA/0_deal/FPKM_mRNA_cancer.txt",header=T,row.names = 1,sep = "\t",check.names = F))
exp_mRNA_normal <- as.matrix(read.table("0_data/TCGA/0_deal/FPKM_mRNA_normal.txt",header=T,row.names = 1,sep = "\t",check.names = F))
exp_lncRNA_cancer <- as.matrix(read.table("0_data/TCGA/0_deal/FPKM_lncRNA_cancer.txt",header=T,row.names = 1,sep = "\t",check.names = F))
exp_lncRNA_normal <- as.matrix(read.table("0_data/TCGA/0_deal/FPKM_lncRNA_normal.txt",header=T,row.names = 1,sep = "\t",check.names = F))

############ cor ##########################

ind_0row_mRNA_cancer <- apply(exp_mRNA_cancer, 1, function(v){all(v == 0)}) # 15
ind_0row_mRNA_normal <- apply(exp_mRNA_normal, 1, function(v){all(v == 0)}) # 516
ind_0row_lncRNA_cancer <- apply(exp_lncRNA_cancer, 1, function(v){all(v == 0)})#25
ind_0row_lncRNA_normal <- apply(exp_lncRNA_normal, 1, function(v){all(v == 0)})#1762

# cor for cancer samples
exp_mRNA_cancer_no0 <- exp_mRNA_cancer[!ind_0row_mRNA_cancer,] # 19406 374
exp_lncRNA_cancer_no0 <- exp_lncRNA_cancer[!ind_0row_lncRNA_cancer,] # 13507 374
write.table(exp_mRNA_cancer_no0,"0_data/TCGA/0_deal/FPKM_mRNA_cancer_no0.txt",quote = F,row.names = T,col.names=T,sep="\t")
write.table(exp_lncRNA_cancer_no0,"0_data/TCGA/0_deal/FPKM_lncRNA_cancer_no0.txt",quote = F,row.names = T,col.names=T,sep="\t")


func_cor <- function(v,m){
  res_i <- apply(m, 1, function(m_i){
    test_i=cor.test(m_i,v)
    test_i_p <- test_i$p.value
    test_i_r <- test_i$estimate
    test_i_rp <- -log10(test_i$p.value)*sign(test_i$estimate)
    return(c(test_i_p,test_i_r,test_i_rp))
  })
  write.table(matrix(res_i[1,],nrow = 1),"1_GSEA/cor_res/p_all_cancer.txt",row.names = F,col.names = F,quote = F,append = T,sep = "\t")
  write.table(matrix(res_i[2,],nrow = 1),"1_GSEA/cor_res/r_all_cancer.txt",row.names = F,col.names = F,quote = F,append = T,sep = "\t")
  write.table(matrix(res_i[3,],nrow = 1),"1_GSEA/cor_res/rp_all_cancer.txt",row.names = F,col.names = F,quote = F,append = T,sep = "\t")
  return(NULL)
}
test_time1 <- Sys.time()
test_res <- apply(log2(exp_lncRNA_cancer_no0+0.01),1,func_cor,log2(exp_mRNA_cancer_no0+0.01))
test_time2 <- Sys.time()
test_time2-test_time1



# cor for normal samples
exp_mRNA_normal_no0 <- exp_mRNA_normal[!ind_0row_mRNA_normal,] # 18905 50
exp_lncRNA_normal_no0 <- exp_lncRNA_normal[!ind_0row_lncRNA_normal,] # 11770 50

func_cor <- function(v,m){
  res_i <- apply(m, 1, function(m_i){
    test_i=cor.test(m_i,v)
    test_i_p <- test_i$p.value
    test_i_r <- test_i$estimate
    test_i_rp <- -log10(test_i$p.value)*sign(test_i$estimate)
    return(c(test_i_p,test_i_r,test_i_rp))
  })
  write.table(matrix(res_i[1,],nrow = 1),"1_GSEA/cor_res/p_all_normal.txt",row.names = F,col.names = F,quote = F,append = T,sep = "\t")
  write.table(matrix(res_i[2,],nrow = 1),"1_GSEA/cor_res/r_all_normal.txt",row.names = F,col.names = F,quote = F,append = T,sep = "\t")
  write.table(matrix(res_i[3,],nrow = 1),"1_GSEA/cor_res/rp_all_normal.txt",row.names = F,col.names = F,quote = F,append = T,sep = "\t")
  return(NULL)
}
test_time1 <- Sys.time()
test_res <- apply(log2(exp_lncRNA_normal_no0+0.01),1,func_cor,log2(exp_mRNA_normal_no0+0.01))
test_time2 <- Sys.time()
test_time2-test_time1



################# GSEA ######################
func.unlist.paste <- function(v){
  rr <- unlist(v)
  return(paste(rr,collapse = ","))
}

library(fgsea)
pathways <- gmtPathways("0_data/immune/immu_pathway/Immune.gmt")

##### GSEA for cancer
con <- file("1_GSEA/cor_res/rp_all_cancer.txt","r")
lineCnt = 0
test_time1 <- Sys.time()
while(1){
  oneline = readLines(con, n = 1)
  if(length(oneline) == 0){
    break
  }
  lineCnt = lineCnt+1
  rp_i <- as.numeric(unlist(strsplit(oneline,"\t")))
  names(rp_i) <- rownames(exp_mRNA_cancer_no0)
  if(sum(is.infinite(rp_i))>0){
    rp_i[which(is.infinite(rp_i))] <- 100
  }
  fgseaRes<-fgsea(pathways,stats=rp_i,minSize=1, maxSize=5000, nperm=1000)
  P_ES_Edge<-fgseaRes[,c(1,2,3,4,8)]
  leadingedge<-apply(P_ES_Edge[,5],1,func.unlist.paste)
  P_ES_Edge_end<-cbind(as.matrix(P_ES_Edge[,1:4]),leadingedge)
  lncRES <- 1-2*as.numeric(P_ES_Edge_end[,2])
  GSEA_res_i <- cbind(rownames(exp_lncRNA_cancer_no0)[lineCnt],cbind(P_ES_Edge_end[,1:4],lncRES,P_ES_Edge_end[,5]))
  write.table(GSEA_res_i,"1_GSEA/GSEA_res/cancer_all.txt",row.names = F,col.names = F,quote = F,sep = "\t",append = T)
  GSEA_res_fdr005_lncRES0995 <- GSEA_res_i[which((as.numeric(GSEA_res_i[,4])<0.05)&(as.numeric(GSEA_res_i[,6])>0.995)),]
  write.table(matrix(GSEA_res_fdr005_lncRES0995,ncol = ncol(GSEA_res_i)),"1_GSEA/GSEA_res/cancer_fdr005_lncRES0995.txt",row.names = F,col.names = F,quote = F,sep = "\t",append = T)
}
close(con)
test_time2 <- Sys.time()
test_time2-test_time1


#### GSEA for normal
con <- file("1_GSEA/cor_res/rp_all_normal.txt","r")
lineCnt = 0
test_time1 <- Sys.time()
while(1){
  oneline = readLines(con, n = 1)
  if(length(oneline) == 0){
    break
  }
  lineCnt = lineCnt+1
  rp_i <- as.numeric(unlist(strsplit(oneline,"\t")))
  names(rp_i) <- rownames(exp_mRNA_normal_no0)
  if(sum(is.infinite(rp_i))>0){
    rp_i[which(is.infinite(rp_i))] <- 100
  }
  fgseaRes<-fgsea(pathways,stats=rp_i,minSize=1, maxSize=5000, nperm=1000)
  P_ES_Edge<-fgseaRes[,c(1,2,3,4,8)]
  leadingedge<-apply(P_ES_Edge[,5],1,func.unlist.paste)
  P_ES_Edge_end<-cbind(as.matrix(P_ES_Edge[,1:4]),leadingedge)
  lncRES <- 1-2*as.numeric(P_ES_Edge_end[,2])
  GSEA_res_i <- cbind(rownames(exp_lncRNA_normal_no0)[lineCnt],cbind(P_ES_Edge_end[,1:4],lncRES,P_ES_Edge_end[,5]))
  write.table(GSEA_res_i,"1_GSEA/GSEA_res/normal_all.txt",row.names = F,col.names = F,quote = F,sep = "\t",append = T)
  GSEA_res_fdr005_lncRES0995 <- GSEA_res_i[which((as.numeric(GSEA_res_i[,4])<0.05)&(as.numeric(GSEA_res_i[,6])>0.995)),]
  write.table(matrix(GSEA_res_fdr005_lncRES0995,ncol = ncol(GSEA_res_i)),"1_GSEA/GSEA_res/normal_fdr005_lncRES0995.txt",row.names = F,col.names = F,quote = F,sep = "\t",append = T)
}
close(con)
test_time2 <- Sys.time()
test_time2-test_time1






setwd("E:/jt/0_other/LIHC_immlnc")
save.image("1_GSEA/all.RData")
load("1_GSEA/all.RData")



##pathway
pathway_check <- as.matrix(read.table("E:/jt/0_other/LIHC_immlnc/0_data/immune/immu_pathway/Gene-from_ImmPort.txt",header = T,quote = "",fill = T,sep = "\t"))
genes_all_pathway <- unique(pathway_check[,1]) # 1811

## sig cor number
con_r_cancer <- file("1_GSEA/cor_res/r_all_cancer.txt","r")
con_p_cancer <- file("1_GSEA/cor_res/p_all_cancer.txt","r")
lineCnt = 0
sigCor_num = 0
test_time1 <- Sys.time()
while(1){
  oneline_r = readLines(con_r_cancer, n = 1)
  oneline_p = readLines(con_p_cancer, n = 1)
  if(length(oneline_r) == 0){
    break
  }
  lineCnt = lineCnt+1
  r_i <- as.numeric(unlist(strsplit(oneline_r,"\t")))
  p_i <- as.numeric(unlist(strsplit(oneline_p,"\t")))
  sig_i <- which((abs(r_i)>0.3)&(p_i<0.05))
  sigCor_num <- sum(sigCor_num,length(sig_i))
}
close(con_r_cancer)
close(con_p_cancer)
test_time2 <- Sys.time()
test_time2-test_time1


con_r_normal <- file("1_GSEA/cor_res/r_all_normal.txt","r")
con_p_normal <- file("1_GSEA/cor_res/p_all_normal.txt","r")
lineCnt_normal = 0
sigCor_num_normal = 0
test_time1 <- Sys.time()
while(1){
  oneline_r = readLines(con_r_normal, n = 1)
  oneline_p = readLines(con_p_normal, n = 1)
  if(length(oneline_r) == 0){
    break
  }
  lineCnt_normal = lineCnt_normal+1
  r_i <- as.numeric(unlist(strsplit(oneline_r,"\t")))
  p_i <- as.numeric(unlist(strsplit(oneline_p,"\t")))
  sig_i <- which((abs(r_i)>0.3)&(p_i<0.05))
  sigCor_num_normal <- sum(sigCor_num_normal,length(sig_i))
}
close(con_r_normal)
close(con_p_normal)
test_time2 <- Sys.time()
test_time2-test_time1
