##### GSEA & diffexp ######
setwd("E:/jt/0_other/LIHC_immlnc")

##### GSEA lncRNA #####
GSEA_sigRes_cancer <- (read.table("1_GSEA/GSEA_res/cancer_fdr005_lncRES0995.txt",header = F,sep = "\t"))
GSEA_sigRes_normal <- (read.table("1_GSEA/GSEA_res/normal_fdr005_lncRES0995.txt",header = F,sep = "\t"))
pairs_cancer <- apply(GSEA_sigRes_cancer,1,function(v){paste(v[1:2],collapse = ",")})
pairs_normal <- apply(GSEA_sigRes_normal,1,function(v){paste(v[1:2],collapse = ",")})
pairs_inter <- intersect(pairs_cancer,pairs_normal) # 3134
length(pairs_inter)/nrow(GSEA_sigRes_cancer) # 0.2269
length(pairs_inter)/nrow(GSEA_sigRes_normal) # 0.2204
lncs_inter <- intersect(as.character(GSEA_sigRes_cancer[,1]),
                        as.character(GSEA_sigRes_normal[,1])) # 1748
length(lncs_inter)/length(unique(as.character(GSEA_sigRes_cancer[,1]))) # /3732 = 0.4684
length(lncs_inter)/length(unique(as.character(GSEA_sigRes_normal[,1]))) # /5322 = 0.3284
lncs_cancer_only <- setdiff(as.character(GSEA_sigRes_cancer[,1]),
                            as.character(GSEA_sigRes_normal[,1])) # 1984 -use

pairs_cancer_only <- setdiff(pairs_cancer,pairs_normal) # 10676
lncs_inPairs_cancer_only <- unique(unlist(strsplit(pairs_cancer_only,","))[rep(c(TRUE,FALSE),length(pairs_cancer_only))])
length(lncs_inPairs_cancer_only) # 3478

################## enrich pathway analysis
pathway_names <- as.matrix(read.table("0_data/immune/immune_names.txt",header = F,sep = "\t"))
pathways <- pathway_names[,2]
names(pathways) <- pathway_names[,1]

pathways_enrich_cancer <- as.data.frame(table(GSEA_sigRes_cancer[,2]))
pathways_enrich_cancer$names <- pathways[as.character(pathways_enrich_cancer[,1])]
pathways_enrich_normal <- as.data.frame(table(GSEA_sigRes_normal[,2]))
pathways_enrich_normal$names <- pathways[as.character(pathways_enrich_normal[,1])]

pathwayNum_plot <- data.frame(pathways=c(as.character(pathways_enrich_cancer$names),as.character(pathways_enrich_normal$names)),
                              num_lncs=c(as.numeric(pathways_enrich_cancer$Freq),-(as.numeric(pathways_enrich_normal$Freq))),
                              class=c(rep("cancer",nrow(pathways_enrich_cancer)),rep("normal",nrow(pathways_enrich_normal))))
ggplot(pathwayNum_plot,aes(x=pathways,y=num_lncs,fill=as.factor(class)))+geom_bar(stat ="identity")+scale_fill_manual(values=alpha(c("#6495ED","#FFA500"), 0.5))+
  coord_flip()+theme( legend.key = element_blank())+ theme_set(theme_bw())


############### cancerImmLnc analysis
lncs_cancer_only
###### pairs
length(lncs_cancer_only)/length(unique(unique(GSEA_sigRes_cancer[,1]))) # 0.5316 53.16%
lncs_cancer_only_Pairs <- GSEA_sigRes_cancer[which(GSEA_sigRes_cancer[,1] %in% lncs_cancer_only),]
nrow(lncs_cancer_only_Pairs)/nrow(GSEA_sigRes_cancer) # 0.4084721  40.85% 5641/13810
###### pathways
pathwayNum_perLncAll <- nrow(GSEA_sigRes_cancer)/length(unique(GSEA_sigRes_cancer[,1])) # 3.70
pathwayNum_perLncOnly <- nrow(lncs_cancer_only_Pairs)/length(unique(lncs_cancer_only_Pairs[,1])) # 2.84

num_pathways_cancer <- table(GSEA_sigRes_cancer[,2])
pathways_ratio_cancer <- num_pathways_cancer/sum(num_pathways_cancer)
num_pathways_cancerOnly <- table(lncs_cancer_only_Pairs[,2])
pathways_ratio_cancerOnly <- num_pathways_cancerOnly/sum(num_pathways_cancerOnly)
pathways_ratio2_cancer <- pathways_ratio_cancer/(pathways_ratio_cancer+pathways_ratio_cancerOnly)
pathways_ratio_df <- data.frame(number=c(num_pathways_cancer,num_pathways_cancerOnly),
                                freq=c(pathways_ratio_cancer,pathways_ratio_cancerOnly),
                                freq_2=c(pathways_ratio2_cancer,1-pathways_ratio2_cancer),
                                pathwyas=c(pathways[names(num_pathways_cancer)],pathways[names(num_pathways_cancerOnly)]),
                                type=c(rep("lncRNA_cancer_all",length(num_pathways_cancer)),rep("lncRNA_cancer_only",length(num_pathways_cancerOnly))))
ggplot(pathways_ratio_df, aes(x=pathwyas, y=freq, fill=type)) +geom_bar(stat ="identity",position = "dodge")+
  scale_fill_manual(values=alpha(c("#6495ED","#2a5caa"), 0.5))+coord_flip()+theme( legend.key = element_blank())+ theme_set(theme_bw())
ggplot(pathways_ratio_df, aes(x=pathwyas, y=freq_2, fill=type)) +geom_bar(stat ="identity",position = "stack")+
  scale_fill_manual(values=alpha(c("#6495ED","#2a5caa"), 0.5))+coord_flip()+theme( legend.key = element_blank())+ theme_set(theme_bw())



##### diffexp lncRNA #####
diffLnc_all <- (read.table("2_diffexp/glmLRT.csv",header = T,sep = ","))
diffLnc_fdr005 <- diffLnc_all[which(diffLnc_all$FDR<0.05),] # 5244 - use
diffLnc_fdr005_logFC1 <- diffLnc_fdr005[which(abs(diffLnc_fdr005$logFC)>1),] # 3420
sum(diffLnc_fdr005[,2]<0) # 4218(upReg)
sum(diffLnc_fdr005[,2]>0) # 1026(downReg)

ratio_diff_allLnc <- nrow(diffLnc_fdr005)/nrow(diffLnc_all) # 0.393 39.3%  5244/13342
ratio_diff_allimmuneLnc_cancer <- length(intersect(diffLnc_fdr005[,1],unique(GSEA_sigRes_cancer[,1])))/length(unique(GSEA_sigRes_cancer[,1]))
#1193/3732 0.32 32%
ratio_diff_allimmuneLnc_noraml <- length(intersect(diffLnc_fdr005[,1],unique(GSEA_sigRes_normal[,1])))/length(unique(GSEA_sigRes_normal[,1]))
#2359/5322 0.4432 44.32%%
ratio_diff_onlyimmuneLnc_cancer <- length(intersect(diffLnc_fdr005[,1],lncs_cancer_only))/length(lncs_cancer_only)
#498/1984 0.251

##### diffexp cancer imm-lncRNA #####
lncs_cancer_imm_diffFDR <- intersect(lncs_cancer_only,diffLnc_fdr005[,1]) # 498 - use
lncs_cancer_imm_diffFDRFC <- intersect(lncs_cancer_only,diffLnc_fdr005_logFC1[,1]) # 312
lncs_cancer_imm_diffFDR_diffres <- diffLnc_fdr005[which(diffLnc_fdr005[,1]%in%lncs_cancer_imm_diffFDR),]

###### pairs
lncs_immDiff_Pairs <- GSEA_sigRes_cancer[which(GSEA_sigRes_cancer[,1] %in% lncs_cancer_imm_diffFDR),]
dim(lncs_immDiff_Pairs) # 1482 7
###### pathways
lncNum_perpathway_immDiff <- nrow(lncs_immDiff_Pairs)/length(unique(lncs_immDiff_Pairs[,2])) # 92.63

num_pathways_immDiffPairs <- table(lncs_immDiff_Pairs[,2])
names(num_pathways_immDiffPairs) <- pathways[names(num_pathways_immDiffPairs)]
ggplot(as.data.frame(num_pathways_immDiffPairs), aes(x=as.character(Var1), y=Freq)) +geom_bar(stat ="identity",fill="#deab8a")+
  coord_flip()+theme( legend.key = element_blank())+ theme_set(theme_bw())

pathways_ratio_immDiff <- num_pathways_immDiffPairs/sum(num_pathways_immDiffPairs)
pathways_ratio_cancerOnly
pathways_ratio3_immDiff <- pathways_ratio_immDiff/(pathways_ratio_immDiff+pathways_ratio_cancerOnly)
pathways_ratio_df_new <- data.frame(
                                freq=c(pathways_ratio3_immDiff,1-pathways_ratio3_immDiff),
                                pathwyas=c(pathways[names(pathways_ratio3_immDiff)],pathways[names(pathways_ratio_cancerOnly)]),
                                type=c(rep("lncRNA1_immDiff",length(pathways_ratio3_immDiff)),rep("lncRNA2_cancer_only",length(pathways_ratio_cancerOnly))))
ggplot(pathways_ratio_df_new, aes(x=pathwyas, y=freq, fill=type)) +geom_bar(stat ="identity",position = "stack")+
  scale_fill_manual(values=alpha(c("#deab8a","#2a5caa"), 0.5),breaks=c("lncRNA1_immDiff","lncRNA2_cancer_only"),labels=c("lncRNA_immDiff","lncRNA_cancer_only"))+
  coord_flip()+theme( legend.key = element_blank())+ theme_set(theme_bw())

#### GSEA plot
library(fgsea)
cat(rownames(exp_mRNA_cancer_no0),file = "2_diffexp/res_analysis/rp_immDiffLncs.txt",sep = "\t",append = T)
cat("\n",file = "2_diffexp/res_analysis/rp_immDiffLncs.txt",sep = "\t",append = T)
all_lncs <- rownames(exp_lncRNA_cancer_no0)
con <- file("1_GSEA/cor_res/rp_all_cancer.txt","r")
lineCnt = 0
test_time1 <- Sys.time()
while(1){
  oneline = readLines(con, n = 1)
  if(length(oneline) == 0){
    break
  }
  lineCnt = lineCnt+1
  lnc_i <- all_lncs[lineCnt]
  if(lnc_i %in% lncs_cancer_imm_diffFDR){
    rp_i <- as.numeric(unlist(strsplit(oneline,"\t")))
    if(sum(is.infinite(rp_i))>0){
      rp_i[which(is.infinite(rp_i))] <- 100
    }
    cat(c(lnc_i,rp_i),file = "2_diffexp/res_analysis/rp_immDiffLncs.txt",sep = "\t",append = T)
    cat("\n",file = "2_diffexp/res_analysis/rp_immDiffLncs.txt",sep = "\t",append = T)
    
  }else{
    next
  }
}
close(con)
test_time2 <- Sys.time()
test_time2-test_time1

pathways_list <- gmtPathways("0_data/immune/immu_pathway/Immune.gmt")
rp_immDiffLncs <- as.matrix(read.table("2_diffexp/res_analysis/rp_immDiffLncs.txt",header = T,row.names = 1,sep = "\t"))
fgseaRes<-fgsea(pathways,stats=rp_i,minSize=1, maxSize=5000, nperm=1000)
lncs_i <- "LINC02550"
exampleRanks <- rp_immDiffLncs[lncs_i,]
pathways_i <- pathways["Immune6"]
GSEA_res_i <- lncs_immDiff_Pairs[which((lncs_immDiff_Pairs[,1]==lncs_i)&(lncs_immDiff_Pairs[,2]==names(pathways_i))),]
plotEnrichment(pathways_list[[names(pathways_i)]],
               exampleRanks) + labs(title=paste0(lncs_i," enriched in ",pathways_i,", FDR=",GSEA_res_i[4],", ES=",GSEA_res_i[5]))



####### immCell ###########

con <- file("3_immCell/immCell_geneSet.txt","r")
immCells <- c()
geneSets <- c()
while(1){
  oneline = readLines(con, n = 1)
  if(length(oneline) == 0){
    break
  }
  geneSet_i <- as.character(unlist(strsplit(oneline,"\t")))
  geneSet_i <- geneSet_i[!(geneSet_i=="")]
  immCells <- c(immCells,geneSet_i[2:length(geneSet_i)])
  geneSets <- c(geneSets,rep(geneSet_i[1],(length(geneSet_i)-1)))
}
close(con)
immCell_geneSet <- cbind(geneSets,immCells)
colnames(immCell_geneSet) <- c("immCell","geneSet")
write.table(immCell_geneSet,"3_immCell/immCell_geneSet_df.txt",row.names = F,quote = F,sep = "\t")

##### immCell cor ######
library(psych)
lncs_cancer_imm_diffFDR # 498 character
lncs_cancer_imm_diffFDRFC # 312
exp_mRNA_cancer_no0 <- as.matrix(read.table("0_data/TCGA/0_deal/FPKM_mRNA_cancer_no0.txt",header=T,row.names = 1,sep = "\t",check.names = F))
exp_lncRNA_cancer_no0 <- as.matrix(read.table("0_data/TCGA/0_deal/FPKM_lncRNA_cancer_no0.txt",header=T,row.names = 1,sep = "\t",check.names = F))

immDiff_lncExp_cancer <- exp_lncRNA_cancer_no0[lncs_cancer_imm_diffFDR,]

func_mtx2df <- function(m){
  row_name <- rep(rownames(m),ncol(m))
  col_name <- rep(colnames(m),each=nrow(m))
  df <- data.frame(row_name,col_name,as.vector(m))
  return(df)
}

sigPairs <- c()
immCells_uniq <- unique(immCell_geneSet[,1])
for (i in 1:length(immCells_uniq)) {
  immCell_i <- immCells_uniq[i]
  immCell_genes_i <- immCell_geneSet[which(immCell_geneSet[,1]==immCell_i),2]
  imm_mRNAexp_cancer <- exp_mRNA_cancer_no0[intersect(immCell_genes_i,rownames(exp_mRNA_cancer_no0)),]
  res_cor_cancer <- corr.test(t(log2(immDiff_lncExp_cancer+0.01)),t(log2(imm_mRNAexp_cancer+0.01)))
  ind_absR_0.03 <- (abs(res_cor_cancer$r)>0.3)
  ind_p_0.05 <- res_cor_cancer$p<0.05
  ind_sigPairs <- ind_absR_0.03&ind_p_0.05
  Pairs_all <- func_mtx2df(ind_sigPairs)
  sigPairs_i <- Pairs_all[which(Pairs_all[,3]),1:2]
  sigPairs <- rbind(sigPairs,cbind(immCell_i,sigPairs_i))
}

###### all lnc cor #######
exp_lncRNA_cancer_no0
exp_mRNA_cancer_no0
con1 <- file("1_GSEA/cor_res/p_all_cancer.txt","r")
con2 <- file("1_GSEA/cor_res/r_all_cancer.txt","r")
pairs_sig_all <- c()
i <- 0
time1= Sys.time()
while(1){
  oneline_p = readLines(con1, n = 1)
  oneline_r = readLines(con2, n = 1)
  if((length(oneline_p)+length(oneline_r)) == 0){
    break
  }
  i <- i+1
  if(rownames(exp_lncRNA_cancer_no0)[i] %in% lncs_cancer_imm_diffFDR){
    p_i <- as.numeric(unlist(strsplit(oneline_p,"\t")))
    r_i <- as.numeric(unlist(strsplit(oneline_r,"\t")))
    ind_sig <- which((p_i<0.05)&(abs(r_i)>0.3))
    if(length(ind_sig)>0){
      pairs_sig_all <- rbind(pairs_sig_all,
                             cbind(rownames(exp_lncRNA_cancer_no0)[i],rownames(exp_mRNA_cancer_no0)[ind_sig]))
      
    }
  }
}
close(con1);close(con2)
time2= Sys.time()
time2-time1
#write.table(immCell_geneSet,"3_immCell/immCell_geneSet_df.txt",row.names = F,quote = F,sep = "\t")

####### immlnc_cor_immcell v.s. immlnc_cor_all###########
immCells_uniq
immCells_geneNumber <- table(immCell_geneSet[,1])
immlnc_cor_inter <- intersect(unique(pairs_sig_all[,1]),unique(sigPairs[,2]))
N=19406
immLnc_phyperRes_all <- c()
for (i in 1:length(immlnc_cor_inter)) {
  lnc_i <- immlnc_cor_inter[i]
  m <- sum(pairs_sig_all[,1]==lnc_i)
  for (j in 1:length(immCells_uniq)) {
    immCells_j <- immCells_uniq[j]
    k <- immCells_geneNumber[immCells_j]
    x <- sum((sigPairs[,2]==lnc_i)&(sigPairs[,1]==immCells_j))
    p <- ifelse(x>=3,1-phyper(x-1,m,N-m,k),NA)
    immLnc_phyperRes_ij <- c(lnc_i,immCells_j,x,m,k,p)
    immLnc_phyperRes_all <- rbind(immLnc_phyperRes_all,immLnc_phyperRes_ij)
  }
}

rownames(immLnc_phyperRes_all) <- NULL
immLnc_phyperRes_all_df <- data.frame(immlncs=immLnc_phyperRes_all[,1],
                                      immCells=immLnc_phyperRes_all[,2],
                                      x=as.numeric(immLnc_phyperRes_all[,3]),
                                      m=as.numeric(immLnc_phyperRes_all[,4]),
                                      k=as.numeric(immLnc_phyperRes_all[,5]),
                                      p=as.numeric(immLnc_phyperRes_all[,6]))
write.table(immLnc_phyperRes_all_df,"3_immCell/immCells/immlnc_phyper_info_all.txt",row.names = F,quote = F,sep = "\t")
immLnc_phyper_p <- matrix(immLnc_phyperRes_all_df$p,nrow = 315,byrow = T)
rownames(immLnc_phyper_p) <- immlnc_cor_inter
colnames(immLnc_phyper_p) <- immCells_uniq

immlnc_sigPhyper_num <- apply(immLnc_phyper_p, 1, function(v){
  ind_na <- which(is.na(v))
  if(length(ind_na)>0){v[ind_na]<-10}
  ind_p005 <- which(v<0.05)
  return(c(length(ind_p005),length(ind_na)))
})
immlnc_sigPhyper_num <- t(immlnc_sigPhyper_num)

##### sigRes
lncs_sigImm <- rownames(immlnc_sigPhyper_num[which(immlnc_sigPhyper_num[,1]>=10),])
exp_immRelatedLnc <- exp_lncRNA_cancer_no0[lncs_sigImm,]
write.table(exp_immRelatedLnc,"4_cancerSubtype/kidney/exp_immRelatedLnc_merge.txt",quote = F,sep = "\t")

lncs_sigImmDiff_Pairs <- GSEA_sigRes_cancer[which(GSEA_sigRes_cancer[,1] %in% lncs_sigImm),]
dim(lncs_sigImmDiff_Pairs) # 120 7


##### heatmap
library(pheatmap)
log10P <- -log10(immLnc_phyper_p[lncs_sigImm,])
cell_sigSum <- colSums(log10P,na.rm = T)
ind_cell_sort <- order(cell_sigSum,decreasing = T)
lncs_sigSum <- rowSums(log10P,na.rm = T)
ind_lncs_sort <- order(lncs_sigSum,decreasing = T)
pheatmap(log10P[ind_lncs_sort,ind_cell_sort], cluster_row = F, cluster_col = F,show_rownames = T,na_col = "#BEBEBE",
         col = colorRampPalette(c("white","#FFA500"))(40),display_numbers = TRUE, number_format="%.2f")


####### subtype ########

library(ConsensusClusterPlus)
setwd("E:/jt/0_other/LIHC_immlnc/4_subtype")

results = ConsensusClusterPlus(exp_immRelatedLnc,maxK=6,reps=50,
                               pItem=0.8,pFeature=1,title = "pam_pearson",
                               clusterAlg="pam",distance="pearson",
                               seed=1262118388.71279,plot="png",writeTable = T)
##### heatmap exp  
sample_class <- as.matrix(read.table("pam_pearson/pam.k=2.consensusClass.csv",sep=","))
sample_class[1:3,]
samples_1 <- sample_class[which(sample_class[,2]=="1"),1]
samples_2 <- sample_class[which(sample_class[,2]=="2"),1]
exp_lncs_class <- exp_lncRNA_cancer_no0[lncs_sigImm,c(samples_1,samples_2)]
annotation_col = data.frame(sampleGroup = factor(rep(c("group1", "group2"), table(sample_class[,2]))))
rownames(annotation_col) <- c(samples_1,samples_2)
heatmap_res <- pheatmap(log10(exp_lncs_class+0.1), cluster_row = T, cluster_col = F,scale = "column",
         show_rownames = T,show_colnames = F,na_col = "#BEBEBE",
         annotation_col = annotation_col,annotation_colors=list(sampleGroup=c(group1="red",group2="blue")),
         col = colorRampPalette(c("#6495ED","white","#FFA500"))(100))
lncs_sort <- rownames(exp_lncs_class)[heatmap_res$tree_row$order]
lncs_sort_use <- c(lncs_sort[1:3],lncs_sort[16:18],lncs_sort[4:15])
pheatmap(log10(exp_lncs_class[lncs_sort_use,]+0.1), cluster_row = F, cluster_col = F,scale = "column",
         show_rownames = T,show_colnames = F,na_col = "#BEBEBE",cellheight = 12,gaps_col = length(samples_1),
         annotation_col = annotation_col,annotation_colors=list(sampleGroup=c(group1="#f8aba6",group2="#afb4db")),
         col = colorRampPalette(c("#6495ED","white","#FFA500"))(100))

write.table(lncs_sort,"3_immCell/immCells/res_keyLncs.txt",row.names = F,col.names = F,quote = F,sep = "\t")



setwd("E:/jt/0_other/LIHC_immlnc")
save.image("3_immCell/all.RData")
load("3_immCell/all.RData")


#### choose k ####
Kvec = 2:6
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)}#end for i# The optimal K
optK = Kvec[which.min(PAC)]


#### test diffexp ####
exp_lncRNA_cancer <- as.matrix(read.table("0_data/TCGA/0_deal/FPKM_lncRNA_cancer.txt",header=T,row.names = 1,sep = "\t",check.names = F))
exp_lncRNA_normal <- as.matrix(read.table("0_data/TCGA/0_deal/FPKM_lncRNA_normal.txt",header=T,row.names = 1,sep = "\t",check.names = F))

sum(exp_lncRNA_cancer["WWTR1-AS1",],na.rm = T)/ncol(exp_lncRNA_cancer)
sum(exp_lncRNA_normal["WWTR1-AS1",],na.rm = T)/ncol(exp_lncRNA_normal)
