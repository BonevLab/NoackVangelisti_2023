source('scripts/config.R')
source('scripts/aux_functions.R')
source('scripts/plot_functions.R')

source(paste0(main_dir,'results/hg38/hic/config.R'))

library(MPRAnalyze)
library(ggplot2)
library(ggpubr)
library(pals)
library(cowplot)
library(patchwork)
library(circlize)
library(Hmisc)
library(Matrix)
library(matrixStats)
require(LSD)
library(ggrepel)
library(ggpointdensity)
library(grid)
library(Rcpp)
library(ComplexHeatmap)
library(monaLisa)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(stringr)
library(plyr)
library(dplyr)
library(reshape2)
multicoreParam <- MulticoreParam(workers = 30)
theme_set(theme_cowplot())
theme_border <- theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line=element_blank(),panel.border = element_rect(color = "black", size = 0.5))
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}

FigureS7A<-function(MPRA_data,MPRA_design,cutoff,groups,newnames,out_f,colours,height= 5, width= 4){
  file <- read.table(MPRA_data)
  colnames(file) <- c("seq", "number")
  file <- file[file$number > cutoff,]
  file$label <- str_split_fixed(file$seq, "_", 2)[,1]
  out<-as.data.frame(str_split_fixed(as.character(file$seq), "_", 3))
  out$anno<-"WT"
  out$anno[grepl("mut",as.character(out$V1))]<-"Mut"
  out$anno[grepl("scr",as.character(out$V1))]<-"Scr"
  out$anno[grepl("NEUROG2",as.character(out$V1))]<-"Mut"
  file$label<-out$anno
  stats <- file %>% 
    group_by(label) %>% 
    dplyr::summarize(mean = mean(number), median = median(number))
  design <- read.table(MPRA_design, sep="\t")
  stats <- file %>% 
    group_by(label) %>% 
    dplyr::summarize(mean = mean(number), median = median(number))
  out<-as.data.frame(str_split_fixed(as.character(design$V1), "_", 3))
  out$anno<-"WT"
  out$anno[grepl("mut",as.character(out$V1))]<-"Mut"
  out$anno[grepl("scr",as.character(out$V1))]<-"Scr"
  out$anno[grepl("NEUROG2",as.character(out$V1))]<-"Mut"
  design$V2<-out$anno
  counter_design <- table(design$V2)
  counter_file <- table(file$label)
  counters <- data.frame(counter_file,counter_design[names(counter_design) %in% names(counter_file)],names = names(counter_file))
  counters <- counters[,c(2,4)]
  colnames(counters) <- c("obs", "exp")
  counters$fraction <- round(counters$obs / counters$exp, 3)
  stats <- cbind.data.frame(counters, stats)
  stats_subset <- stats
  stats_subset <- stats[stats$label %in% groups,]
  stats_subset$label <- factor(stats_subset$label, levels = groups, ordered=T)
  p <- ggplot(data = stats_subset, aes(x = label, y = fraction, fill = label)) + geom_bar(stat="identity",colour="black")+
    ylab("Fraction of recovered CREs") +  
    scale_y_continuous(breaks = seq(0,1,0.1),expand = c(0, 0))+scale_fill_manual(values=colours) + 
    theme(axis.title.x=element_blank(),legend.title = element_blank(), legend.position = 'none')+
    geom_hline(yintercept=1.0, linetype="dashed",color = "blue", size=1)
  pdf(out_f, height= height, width= width)
  print(p)
  dev.off()
}

FigureS7B<-function(MPRA_data,MPRA_design,cutoff,groups,newnames,out_f,colours,height= 5, width= 4){
  file <- read.table(MPRA_data)
  design <- read.table(MPRA_design, sep="\t")
  colnames(file) <- c("seq", "number")
  file <- file[file$number > cutoff,]
  out<-as.data.frame(str_split_fixed(as.character(file$seq), "_", 3))
  out$anno<-"WT"
  out$anno[grepl("mut",as.character(out$V1))]<-"Mut"
  out$anno[grepl("scr",as.character(out$V1))]<-"Scr"
  out$anno[grepl("NEUROG2",as.character(out$V1))]<-"Mut"
  file$label<-out$anno
  file_subset <- file[file$label %in% groups,]
  file_subset$label <- factor(file_subset$label, levels =groups, ordered=T)
  p <- ggplot(file_subset, aes(x = label, y = log2(number), fill=label)) +
    geom_violin() +
    geom_boxplot(width=0.1, outlier.shape = NA) + scale_fill_manual(values=colours) +
    theme(axis.title.x=element_blank(),legend.title = element_blank(), legend.position = 'none')+ylab("log2 (Number of Barcodes per CRE)")+
    scale_y_continuous(breaks = seq(0,10,2)) 
  pdf(out_f,height= height, width= width)
  plot(p)
  dev.off()
}

FigureS7D <- function(file_f,conditions,reps,which_genes,cols,ylim=NULL){
  mat <- read.csv(file_f,header=T,row.names = 'Gene')
  row.names(mat) <- gsub("TBR2","EOMES",row.names(mat))
  s_names <- paste0(rep(conditions,each=2),reps)
  df <- mat[which_genes,]
  df$Gene <- factor(which_genes,levels=which_genes)
  df <- melt(df,varnames = c('Gene'),value.name = 'FC')
  df$Condition <- factor(rep(conditions,each=2*length(which_genes)),levels=conditions)
  df_summ <- data_summary(df, varname="FC",groupnames=c('Condition','Gene'))
  p <- ggplot(df_summ, aes(x=Gene, y=FC, fill=Condition)) + scale_fill_manual(values = cols) +
    geom_bar(stat="identity", color="black",position=position_dodge()) + xlab('') + ylab('Relative Fold Change') +
    geom_errorbar(aes(ymin=FC-sd, ymax=FC+sd), width=.2,position=position_dodge(.9))
  p <- p+ geom_point(data = df,aes(x=Gene, y=FC, fill=Condition),position = position_jitterdodge(dodge.width =.9,jitter.width = 0.25,seed = 42),size=2) 
  if(!is.null(ylim)){
    p <- p+ylim(ylim)
  } 
  return(p)
}

FigureS7E<-function(MPRA_data_f,MPRA_design,cutoff,out_f,cols,celltypes,organ, width=4,height=5){ 
  names<-c("Scr", "WT","Mut")
  lables<- read.table(MPRA_design, sep="\t")
  out<-as.data.frame(str_split_fixed(as.character(lables$V1), "_", 3))
  out$anno<-"WT"
  out$anno[grepl("mut",as.character(out$V1))|grepl("NEUROG2",as.character(out$V1))]<-"Mut"
  out$anno[grepl("scr",as.character(out$V1))]<-"Scr"
  lables$V2<-out$anno
  to_plot<-data.frame()
  for (c in 1:length(celltypes)) {
    celltype<-celltypes[c]
    dnaCounts <- read.delim(paste0(MPRA_data_f,celltype,'_',organ,'/MPRA_3DRAM_',celltype,'/dna_counts.tsv'),header=T,row.names = 'seq_id')
    rep1<-dnaCounts[1:(1+1)!=(1+1)]
    rep1$barcodes<-rowSums(rep1>=1, na.rm = T)
    rep1<-subset(rep1,rep1$barcodes >= cutoff)
    out<-as.data.frame(str_split_fixed(as.character(row.names(rep1)), "_", 3))
    out$anno<-"WT"
    out$anno[grepl("mut",as.character(out$V1))|grepl("NEUROG2",as.character(out$V1))]<-"Mut"
    out$anno[grepl("scr",as.character(out$V1))]<-"Scr"
    rep1$X1<-out$anno
    rep2<-dnaCounts[1:(1+1)==(1+1)]
    rep2$barcodes<-rowSums(rep2>=1, na.rm = T)
    rep2<-subset(rep2,rep2$barcodes >= cutoff)
    out<-as.data.frame(str_split_fixed(as.character(row.names(rep2)), "_", 3))
    out$anno<-"WT"
    out$anno[grepl("mut",as.character(out$V1))|grepl("NEUROG2",as.character(out$V1))]<-"Mut"
    out$anno[grepl("scr",as.character(out$V1))]<-"Scr"
    rep2$X1<-out$anno
    name_Vec<-unique(grep(lables$V2, pattern="",value = TRUE))
    rep_vec<-c("rep1","rep2")
    for (k in 1:length(rep_vec)) {
      df<-data.frame()
      filtered_barcodes<-get(rep_vec[k])
      for (i in 1:length(name_Vec)) {
        term<-name_Vec[i]
        number<-as.data.frame(nrow(subset(filtered_barcodes,filtered_barcodes$X1 == term)) / nrow(subset(lables,lables$V2 == term)))*100
        colnames(number)<-"Percentage"
        number$type<-term
        df<-rbind(df,number)
      }
      df$replicate<-rep_vec[k]
      df$celltype<-celltype
      to_plot<-rbind(to_plot,df)
    }
  }
  to_plot[grepl("Scr",to_plot$type),'class']<-"Scr"
  to_plot[grepl("WT",to_plot$type),'class']<-"WT"
  to_plot[grepl("Mut",to_plot$type),'class']<-"Mut"
  to_plot<-na.omit(to_plot)
  to_plot<-to_plot %>% dplyr::group_by(class,celltype) %>% mutate(Mean=mean(Percentage))
  to_plot<-to_plot %>% dplyr::group_by(class,celltype,replicate) %>% mutate(Mean_replicate=mean(Percentage))
  to_plot$class <- factor(to_plot$class,levels = c("Scr","WT","Mut"))
  to_plot$celltype <- factor(to_plot$celltype,levels = c("Pax6", "Tbr2", "PN"))
  levels(to_plot$celltype) <- c('RGC','IPC','N')
  
  to_plot$Mean<-to_plot$Mean/100
  to_plot$Percentage<-to_plot$Percentage/100
  to_plot$Mean_replicate<-to_plot$Mean_replicate/100
  to_plot$type<-NULL
  to_plot$Percentage<-NULL
  to_plot<-unique(to_plot)
  p<-ggplot(to_plot, aes(y=Mean,x=class,fill=celltype))+geom_bar(stat="identity",position="dodge", color="black")+
    xlab(label = 'CRS Type')+ylab(label = paste0('Percentage recovery')) +
    ylab("Fraction of recovered CREs") +  
    scale_y_continuous(breaks = seq(0,1,0.1))+
    theme(axis.title.x=element_blank(),legend.title = element_blank())+
    geom_hline(yintercept=1.0, linetype="dashed",color = "blue", size=1)
  p<-p+ geom_point(data = to_plot,aes(x=class, y=Mean_replicate, fill=celltype),position = position_jitterdodge(dodge.width =.9,jitter.width = 0.25,seed = 42),size=2)+
    theme(legend.position='top',legend.justification = "center")+scale_fill_manual(values = cols) 
  pdf(file=out_f, width=width,height=height)
  print(p)
  dev.off()  
}

FigureS7E_2<-function(MPRA_data_f,out_f,cols,names,celltypes,organ){
  to_plot<-data.frame()
  for (c in 1:length(celltypes)) {
    celltype<-celltypes[c]
    dnaCounts <- read.delim(paste0(MPRA_data_f,celltype,'_',organ,'/MPRA_3DRAM_',celltype,'/dna_counts.tsv'),header=T,row.names = 'seq_id')
    rep1<-dnaCounts[1:(1+1)!=(1+1)]
    rep1$barcodes<-rowSums(rep1>=1, na.rm = T)
    rep1$replicate<-'Rep1'
    rep2<-dnaCounts[1:(1+1)==(1+1)]
    rep2$barcodes<-rowSums(rep2>=1, na.rm = T)
    rep2$replicate<-'Rep2'
    df<-rbind(rep1[,c('barcodes','replicate')],rep2[,c('barcodes','replicate')])
    df$celltype<-names[c]
    to_plot<-rbind(to_plot,df)
  }
  to_plot$celltype <- factor(to_plot$celltype,levels = names)
  p <- ggplot(to_plot, aes(x = celltype, y = log2(barcodes), fill=celltype))+
    geom_violin() +geom_boxplot(width=0.1, outlier.shape = NA) +
    theme_classic() +ylab("log2 (Number of Barcodes per CRE)") +
    theme(axis.title.x=element_blank(),legend.position = "none") +scale_y_continuous(breaks = seq(0,10,2))+scale_fill_manual(values = cols) 
  pdf(file=out_f, height = 5, width=3)
  print(p)
  dev.off() 
}

FigureS7F<-function(MPRA_data_f,remove=T,out_f,cols,organ="hg",width=5,height=5){
  celltypes<-c('Pax6','Tbr2','PN')
  merged_df<-data.frame
  for (i in 1:length(celltypes)) {
    celltype<-celltypes[i]
    dnaCounts <- read.delim(paste0(MPRA_data_f,celltype,'_',organ,'/MPRA_3DRAM_',celltype,'/dna_counts.tsv'),header=T,row.names = 'seq_id')
    rnaCounts <- read.delim(paste0(MPRA_data_f,celltype,'_',organ,'/MPRA_3DRAM_',celltype,'/rna_counts.tsv'),header=T,row.names = 'seq_id')
    #split in replicates + normalize by total number of reads
    rep1_DNA<-dnaCounts[1:(1+1)!=(1+1)]
    rep1_DNA[is.na(rep1_DNA)]<-0
    rep1_RNA<-rnaCounts[1:(1+1)!=(1+1)]
    rep1_RNA[is.na(rep1_RNA)]<-0
    if (remove==T) {
      rep1_RNA[rep1_DNA==0]<-0
      rep1_DNA[rep1_RNA==0]<-0
    }
    rep1_DNA<-rep1_DNA[]/sum(rep1_DNA)
    rep1_DNA$DNA_read_sum_rep1<-rowSums(rep1_DNA,na.rm=T)
    rep1_DNA$DNA_barcodes_rep1<-rowSums(rep1_DNA>0, na.rm = T)-1
    rep1_RNA<-rep1_RNA[]/sum(rep1_RNA)
    rep1_RNA$RNA_read_sum_rep1<-rowSums(rep1_RNA,na.rm=T)
    rep1_RNA$RNA_barcodes_rep1<-rowSums(rep1_RNA>0, na.rm = T)-1
    ###for replicate 2
    rep2_DNA<-dnaCounts[1:(1+1)==(1+1)]
    rep2_DNA[is.na(rep2_DNA)]<-0
    rep2_RNA<-rnaCounts[1:(1+1)==(1+1)]
    rep2_RNA[is.na(rep2_RNA)]<-0
    if (remove==T) {
      rep2_RNA[rep2_DNA==0]<-0
      rep2_DNA[rep2_RNA==0]<-0
    }
    rep2_DNA<-rep2_DNA[]/sum(rep2_DNA)
    rep2_DNA$DNA_read_sum_rep2<-rowSums(rep2_DNA,na.rm=T)
    rep2_DNA$DNA_barcodes_rep2<-rowSums(rep2_DNA>0, na.rm = T)-1
    rep2_RNA<-rep2_RNA[]/sum(rep2_RNA)
    rep2_RNA$RNA_read_sum_rep2<-rowSums(rep2_RNA,na.rm=T)
    rep2_RNA$RNA_barcodes_rep2<-rowSums(rep2_RNA>0, na.rm = T)-1
    df<-as.data.frame(cbind(rep1_DNA[,c('DNA_read_sum_rep1','DNA_barcodes_rep1')],rep1_RNA[,c('RNA_read_sum_rep1','RNA_barcodes_rep1')],rep2_DNA[,c('DNA_read_sum_rep2','DNA_barcodes_rep2')],rep2_RNA[,c('RNA_read_sum_rep2','RNA_barcodes_rep2')]))
    row.names(df)<-row.names(dnaCounts)
    #####calculate the sum ratio 
    df$ratio_sums_rep1<-(df$RNA_read_sum_rep1+1)/(df$DNA_read_sum_rep1+1)
    df$ratio_sums_rep2<-(df$RNA_read_sum_rep2+1)/(df$DNA_read_sum_rep2+1)
    df$mean_ration_sums<-rowMeans(df[,c("ratio_sums_rep1", "ratio_sums_rep2")], na.rm=TRUE)
    df$name<-row.names(df)
    colnames(df)<-c(paste0(celltype,'_DNA_read_sum_rep1'),paste0(celltype,'_DNA_barcodes_rep1'),
                    paste0(celltype,'_RNA_read_sum_rep1'),paste0(celltype,'_RNA_barcodes_rep1'),
                    paste0(celltype,'_DNA_read_sum_rep2'),paste0(celltype,'_DNA_barcodes_rep2'),
                    paste0(celltype,'_RNA_read_sum_rep2'),paste0(celltype,'_RNA_barcodes_rep2'),
                    paste0(celltype,'_ratio_sums_rep1'),paste0(celltype,'_ratio_sums_rep2'),paste0(celltype,'_mean_ration_sums'),paste0(celltype,'_name')
    )
    ifelse(i==1,merged_df<-df,merged_df<-merge(merged_df,df,by.x=paste0(celltypes[1],'_name'), by.y=paste0(celltype,'_name'), all.x=T)) ###keep only CRS which are in all celltypes !!! otherwise all=T 
  }
  out<-as.data.frame(str_split_fixed(as.character(merged_df$Pax6_name), "_", 3))
  out$anno<-"WT"
  out$anno[grepl("mut",as.character(out$V1))|grepl("NEUROG2",as.character(out$V1))]<-"Mut"
  out$anno[grepl("scr",as.character(out$V1))]<-"Scr"
  merged_df$CDS_type<-out$anno
  #####Generate PCA plot
  df<- merged_df[,c('Pax6_ratio_sums_rep1','Pax6_ratio_sums_rep2','Tbr2_ratio_sums_rep1','Tbr2_ratio_sums_rep2','PN_ratio_sums_rep1','PN_ratio_sums_rep2')]
  row.names(df)<-row.names(merged_df)
  df<-na.omit(df)
  df<-as.data.frame(as.matrix(t(df)))
  df_pca <- prcomp(df)
  df_pca_out<-as.data.frame(df_pca$x)
  df_pca_out$Celltype<-rep(names,each=2)
  df_pca_out$replicate<-c('rep1','rep2','rep1','rep2','rep1','rep2')
  percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
  percentage <- paste(colnames(df_pca_out), "(", paste( as.character(percentage), "%", ")", sep="") )
  df_pca_out$Celltype<- factor(df_pca_out$Celltype, levels = names)
  p<-ggplot(df_pca_out,aes(x=PC1,y=PC2,color=Celltype))+geom_point(size = 5)+xlab(percentage[1]) + ylab(percentage[2])+ 
    theme(legend.position = c(0.75,0.15),legend.title = element_blank(),legend.text = element_text( size = 15))+scale_color_manual(values = cols)
  pdf(file=out_f, width=width,height=height)
  print(p)
  dev.off()
}

FigureS7G<-function(v_enhancer,MPRA_data,group,out_f,ylim,p.ylim =6,width=5,height=5){
  v_enhancer<-read.table(v_enhancer, sep='\t', header=F)
  gr_vista_enhancer<-makeGRangesFromDataFrame(v_enhancer, keep.extra.columns=F,ignore.strand=T,seqnames.field="V1",start.field="V2",end.field="V3")
  df=vroom::vroom(paste0(MPRA_data))
  mut_cor<-subset(df, df$enh_typ == 'scr') 
  mut_cor$type<-'Scr'
  df<-subset(df,df$enh_typ=="WT")
  peaks <- gsub('WT_','',df$names)
  featureSplit <- as.data.frame(stringr::str_split(paste0(peaks), pattern =':' , n = 2, simplify = TRUE))
  featureSplit$V3<-as.numeric(paste0(featureSplit$V2))+135
  featureSplit$V2<-as.numeric(paste0(featureSplit$V2))-135
  peaks <- GRanges(featureSplit$V1,IRanges(as.integer(featureSplit$V2),as.integer(featureSplit$V3)),names=df$names,Pax6_mad.score=df$Pax6_mad.score,
                   Tbr2_mad.score=df$Tbr2_mad.score,PN_mad.score=df$PN_mad.score)
  pos_cor<-as.data.frame(peaks[peaks %over% gr_vista_enhancer])
  pos_cor$type<-'Vista'
  for_ploting<-rbind(pos_cor[,c("Pax6_mad.score","Tbr2_mad.score",'PN_mad.score',"type")],mut_cor[,c("Pax6_mad.score","Tbr2_mad.score",'PN_mad.score',"type")])
  colnames(for_ploting)<-c('RGC','IPC','N','type')
  res_df<-melt(for_ploting)
  res_df <- res_df[complete.cases(res_df),]
  p <- ggplot(res_df,aes(x=type,y=value,fill=type)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend =T,width=0.8) + scale_fill_manual(name='',values=c("#F1F1F1","#6B6B6B"))
  p <- p + coord_cartesian(ylim=ylim) + ylab('MPRA signal') + xlab('Cell type')+ theme(legend.position = "none")+ facet_wrap(~variable) + xlab("")
  p <- p + stat_compare_means(comparisons = list(c('Scr','Vista')),label = "p.format",method = 'wilcox',paired = F,label.y =p.ylim,tip.length = c(0.003,0.001),inherit.aes = T)
  pdf(file=out_f, width=width,height=height)
  print(p)
  dev.off()
}


### Plot Figures


FigureS7A(MPRA_data='/home/fnoack/projects/MPRA/MPRAflow/outs/3DRAM_BC_asso/merged_BC_CRS_asso/filtered_combined.tsv',
          MPRA_design='/home/fnoack/projects/MPRA/dataMPRA_3DRAM/BC-CRS/labels.tsv',
          cutoff=1,colours=c("white",'grey75','grey35'),
          groups=c("Scr", "WT","Mut"),
          out_f="figures/FigureS7A.pdf")

FigureS7B(MPRA_data='/home/fnoack/projects/MPRA/MPRAflow/outs/3DRAM_BC_asso/merged_BC_CRS_asso/filtered_combined.tsv',
          MPRA_design='/home/fnoack/projects/MPRA/dataMPRA_3DRAM/BC-CRS/labels.tsv',
          cutoff=1,colours=c("white",'grey75','grey35'),
          groups=c("Scr", "WT","Mut"),
          out_f="figures/FigureS7B.pdf")

p1 <- FigureS7D(file_f="/home/fnoack/projects/MPRA/dataMPRA_3DRAM/qPCR_results_mm/qPCR_results_human.csv",conditions=c('RGC','IPC','N'),
          reps=1:2,cols=c(cell_colors,'grey50'),which_genes=c('SOX2')) #'Hes1','Pax6','Sox2','Btg2','Eomes','Sox5','Tubb3'
p2 <- FigureS7D(file_f="/home/fnoack/projects/MPRA/dataMPRA_3DRAM/qPCR_results_mm/qPCR_results_human.csv",conditions=c('RGC','IPC','N'),
                reps=1:2,cols=c(cell_colors,'grey50'),which_genes=c('EOMES')) #'Hes1','Pax6','Sox2','Btg2','Eomes','Sox5','Tubb3'
p3 <- FigureS7D(file_f="/home/fnoack/projects/MPRA/dataMPRA_3DRAM/qPCR_results_mm/qPCR_results_human.csv",conditions=c('RGC','IPC','N'),
                reps=1:2,cols=c(cell_colors,'grey50'),which_genes=c('SOX5')) #'Hes1','Pax6','Sox2','Btg2','Eomes','Sox5','Tubb3'
p <- p1+plot_spacer() + (p2+ylab(""))+plot_spacer() +(p3+ylab(""))
pdf('figures/FigureS7D.pdf',height=5.5,width=6)
print(p+plot_layout(guides="collect",widths = c(10,-5,10,-5,10)) & theme(legend.position = "top",legend.title=element_blank(),legend.justification = "center"))
dev.off()

FigureS7E(MPRA_data_f='/home/fnoack/projects/MPRA/MPRAflow/outs/MPRA_3DRAM/',
          MPRA_design='/home/fnoack/projects/MPRA/dataMPRA_3DRAM/BC-CRS/labels.tsv',
          cutoff=1,cols=c(cell_colors,'grey50'),
          celltypes<-c('Pax6','Tbr2','PN'),organ='hg',
          out_f="figures/FigureS7E.pdf",height=5.5,width=5)

FigureS7F(MPRA_data_f='/home/fnoack/projects/MPRA/MPRAflow/outs/MPRA_3DRAM/',
          out_f="figures/FigureS7F.pdf",remove=T,
          cols=c(cell_colors,'grey50'))

FigureS7G(v_enhancer="/home/fnoack/temp/vista_enhancer/ensemble_vista_enhancer/all_forebrain_enhancer_hg38.bed",
          MPRA_data="/home/fnoack/projects/MPRA/dataMPRA_3DRAM/mad_scores_hORG_10BC_threshold.txt",
          out_f="figures/FigureS7G.pdf",
          ylim = c(-3,7),p.ylim =6,height=5,width=5)
