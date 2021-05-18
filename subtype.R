##### subtype #####
library(ggplot2)

chisq_p <- function(df){
  class1_num <- table(as.character(df[,2]))
  class2_num <- table(as.character(df[,3]))
  class1_num_noless <- class1_num[which(class1_num>=5)]
  if(length(class1_num_noless)>1){
    m <- c()
    for (i in 1:length(class1_num_noless)) {
      a_i <- length(intersect(df[which(df[,2]==names(class1_num_noless)[i]),1],
                              df[which(df[,3]==names(class2_num)[1]),1]))
      b_i <- length(intersect(df[which(df[,2]==names(class1_num_noless)[i]),1],
                              df[which(df[,3]==names(class2_num)[2]),1]))
      m <- rbind(m,c(a_i,b_i))
    }
    rownames(m) <- names(class1_num_noless)
    colnames(m) <- names(class2_num)
    chiaqRes <- chisq.test(m)
    return(list(m=m,p=chiaqRes$p.value))
  }else{
    return(list(m="noRes",p=10))
  }
}

##### subtype TCGA clinical analysis #####

sample_class <- read.table("4_subtype/1_immCells/pam_pearson/pam.k=2.consensusClass.csv",sep=",")
sample_info <- read.table("E:/jt/0_other/LIHC_immlnc/0_data/TCGA/Clinical BCR XML.merge.txt",sep="\t",header = T,fill=T,quote = "")

sample_class[,1] <- substr(sample_class[,1],1,nchar("TCGA-MI-A75G"))
sample_class_rep <- sample_class[which(duplicated(sample_class[,1])),]
sample_class[which(sample_class[,1]==sample_class_rep[3,1]),]
sample_class_noRep <- sample_class[-which(duplicated(sample_class[,1])),]
rownames(sample_class_noRep) <- sample_class_noRep[,1]

classes <- colnames(sample_info)
rownames(sample_info) <- sample_info[,1]

for (i in 2:length(classes)) {
  samples_i <- sample_info[,c(1,i)]
  samples_i_noNA <- samples_i[(!is.na(samples_i[,2]))&(samples_i[,2]!=""),]
  samples_i_inter <- intersect(sample_class_noRep[,1],samples_i_noNA[,1])
  df <- cbind(samples_i_noNA[samples_i_inter,],sample_class_noRep[samples_i_inter,2])
  colnames(df) <- c("samples","classes","type")
  
  type_i <- "numeric" # character
  
  if(type_i=="character"){
    table(as.character(df$classes))
    {df$classes <- substr(as.character(df$classes),1,2)}
    {
      df$classes[which(as.character(df$classes)=="Stage IIIA")]="Stage III"
      df$classes[which(as.character(df$classes)=="Stage IIIB")]="Stage III"
      df$classes[which(as.character(df$classes)=="Stage IIIC")]="Stage III"
      df$classes[which(as.character(df$classes)=="Stage IVA")]="Stage IV"
      df$classes[which(as.character(df$classes)=="Stage IVB")]="Stage IV"
    }
    {df <- df[which(df$classes!="Not Available"),]}
    table(as.character(df$classes))
    res <- chisq_p(df)
    p <- res$p
    m <- res$m
    df_plot <- data.frame(classes=rep(rownames(m),2),
                          type=rep(colnames(m),each=nrow(m)),
                          num=as.numeric(m))
    ggplot(df_plot, aes(x = classes,y = num,fill = type))+
      geom_bar(stat ="identity",width = 0.6,position = "dodge")+
      scale_fill_manual(values = c("red","blue"))+
      labs(x = "",y = "", title = classes[i])+
      geom_text(aes(label = df_plot$num),position=position_dodge(width = 0.5),size = 3,vjust = -0.25) ###########设置柱子上的标签文字，文字的position_dodge(width=0.5)设置，保证分隔宽度�?
      #guides(fill = guide_legend(reverse = F))##############图例顺序反转
    
    ggplot(df_plot, aes(x =type ,y = num,fill = classes ))+
      geom_bar(stat ="identity",width = 0.6,position = "dodge")+
      labs(x = "",y = "", title = classes[i])+
      geom_text(aes(label = df_plot$num),position=position_dodge(width = 0.5),size = 3,vjust = -0.25) ###########设置柱子上的标签文字，文字的position_dodge(width=0.5)设置，保证分隔宽度�?
    #guides(fill = guide_legend(reverse = F))##############图例顺序反转
    
  }else{
    p <- wilcox.test(classes ~ type,df)$p.value
    ggplot(df)+
      geom_boxplot(aes(x=as.factor(type),y=classes,fill=as.factor(type)),width=0.6,outlier.size = 0,outlier.color = "white")+ #position = position_dodge(0.8)
      scale_fill_manual(values = c("red", "blue"),breaks=c("1","2"),labels=c("Group 1","Group 2"))+
      xlab("")+
      ylab("")+
      scale_y_continuous(limits = c(0,50))#,breaks=seq(0,110,5)
  }
  
}

######## immune signature ##############33
sample_class_noRep
signatures_all <- read.table("3_immCell/immSignature/LIHC_signature.txt",header = T,row.names = 1,fill = T,quote = "")
signatures <- colnames(signatures_all)
samples_inter <- length(intersect(rownames(signatures_all),rownames(sample_class_noRep))) # 368
for (i in 1:length(signatures)) {
  {
    samples_i <- signatures_all[,i]
    names(samples_i) <- rownames(signatures_all)
    samples_i_noNA <- samples_i[(!is.na(samples_i))&(samples_i!="")]
    samples_i_inter <- intersect(sample_class_noRep[,1],names(samples_i_noNA))
    df <- cbind(names(samples_i_noNA[samples_i_inter]),samples_i_noNA[samples_i_inter],sample_class_noRep[samples_i_inter,2])
    colnames(df) <- c("samples","signatures","type")
  }
  
  type_i <- "numeric" # character
  
  if(type_i=="character"){
    table(as.character(df$signatures))
    {df$signatures <- substr(as.character(df$signatures),1,2)}
    {df <- df[which(df$signatures!="Not Available"),]}
    table(as.character(df$signatures))
    res <- chisq_p(df)
    p <- res$p
    m <- res$m
    df_plot <- data.frame(signatures=rep(rownames(m),2),
                          type=rep(colnames(m),each=nrow(m)),
                          num=as.numeric(m))
    ggplot(df_plot, aes(x = signatures,y = num,fill = type))+
      geom_bar(stat ="identity",width = 0.6,position = "dodge")+
      scale_fill_manual(values = c("red","blue"))+
      labs(x = "",y = "", title = signatures[i])+
      geom_text(aes(label = df_plot$num),position=position_dodge(width = 0.5),size = 3,vjust = -0.25) ###########设置柱子上的标签文字，文字的position_dodge(width=0.5)设置，保证分隔宽度�?
    
    ggplot(df_plot, aes(x =type ,y = num,fill = signatures ))+
      geom_bar(stat ="identity",width = 0.6,position = "dodge")+
      labs(x = "",y = "", title = signatures[i])+
      geom_text(aes(label = df_plot$num),position=position_dodge(width = 0.5),size = 3,vjust = -0.25) ###########设置柱子上的标签文字，文字的position_dodge(width=0.5)设置，保证分隔宽度�?
    
  }else{
    df <- data.frame(samples=df[,1],signatures=as.numeric(df[,2]),type=as.factor(df[,3]))
    table(df$type)
    p <- wilcox.test(signatures ~ type,df)$p.value
    ggplot(df,aes(x=as.factor(type),y=signatures,fill=as.factor(type)))+
      geom_boxplot(notch = T,width=0.6,outlier.size = 0,outlier.color = "black")+
      geom_point(aes(color=as.factor(type)),position="jitter",pch=16,cex=1)+
      labs(x = "",y = "", title = paste0(signatures[i],"   p=",p),subtitle=paste0(table(df$type)["1"]," group1 ,",table(df$type)["2"]," group2" ))+
      scale_fill_manual(values = alpha(c("red", "blue"),0.5),breaks=c("1","2"),labels=c("Group 1","Group 2"))+
      scale_colour_manual(values = c("red", "blue"),breaks=c("1","2"),labels=c("Group 1","Group 2"))+
      theme(plot.title=element_text(hjust=0.5),plot.subtitle =element_text(hjust=0.5) )
  }
}





setwd("E:/jt/0_other/LIHC_immlnc")
save.image("4_subtype/all.RData")
load("4_subtype/all.RData")
