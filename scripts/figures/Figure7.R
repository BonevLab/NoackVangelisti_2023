source('scripts/config.R')
source('scripts/aux_functions.R')
source('scripts/plot_functions.R')

source(paste0(main_dir,'results/hg38/hic/config.R'))
source(paste0(main_f,'scripts/main_functions.R'))
source(paste0(main_f,'scripts/aux_functions.R'))
source(paste0(main_f,'scripts/plot_functions.R'))
source(paste0(main_f,'scripts/temp_functions.R'))

library(doParallel)
registerDoParallel(cores=12)
library(JASPAR2020)
library(Matrix)
library(TFBSTools)
library(motifmatchr)
library(Hmisc)
library(BSgenome.Hsapiens.UCSC.hg38)
library(monaLisa)
library(MPRAnalyze)
library(ggplot2)
library(ggpubr)
library(pals)
library(cowplot)
library(patchwork)
library(circlize)
library(dplyr)
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
library(Biostrings)

theme_set(theme_cowplot())

# Import data
mpra_annot <- read.table('/home/hpc/bonev/projects/mpra/mpra_3dram/data/mpra_loop_anno.tsv',header=T)
res <- read.table('/home/hpc/bonev/projects/mpra/mpra_3dram/results/horg_immunoMPRA_res3.tsv',header=T)


Figure7B<-function(MPRA_data_f,out_f,cols,names,celltypes,organ='hg', height = 5, width=3){
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
  pdf(file=out_f, height = height, width=width)
  print(p)
  dev.off() 
}


Figure7C <- function(res,metric='mad.score',out_f,height,width,ylim=c(-3,30)){
 res_sub <- res[res$NSC_pval.mad<=0.1|res$IPC_pval.mad<=0.1,]
 # res_sub <- res[((res$enh_class=='NSC')&(res$NSC_pval.mad<=0.05))|((res$enh_class=='IPC')&(res$IPC_pval.mad<=0.05)),]
  res_sub_WT <- as.data.frame(res_sub[res_sub$enh_typ%in%'WT',c(paste0(c('NSC','IPC','N'),'_',metric),'enh_typ','enh_class')])
  res_sub_scr <- as.data.frame(res[res$enh_typ=='scr',c(paste0(c('NSC','IPC','N'),'_',metric),'enh_typ','enh_class')])
  res_sub_scr$enh_class <- 'Scr'
  res_sub <- rbind(res_sub_WT,res_sub_scr)
  colnames(res_sub) <- c('RGC','IPC','N','mpra_type','enh_class')
  res_sub$enh_class <- factor(res_sub$enh_class,levels=c('Scr','NSC','IPC'))
  levels(res_sub$enh_class)[2] <- 'RGC'
  res_df <- melt(res_sub,id.vars = c('enh_class'),measure.vars = c('RGC','IPC','N'))
  p <- ggplot(res_df,aes(x=enh_class,y=value,fill=enh_class)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend =F,width=0.8) + scale_fill_manual(name='',values=c("#F1F1F1",cell_colors,"#6B6B6B"))
  p <- p + coord_cartesian(ylim=ylim) + ylab('MPRA signal') + facet_wrap(~variable) + xlab('')
  pdf(paste0('figures/',out_f,'.pdf'),height=height,width=width)
  print(p+xlab(''))
  dev.off()
}

Figure7D <- function(res,k=5,cols,out_f,width,height,seed){
  df <- as.data.frame(res[res$enh_typ=='WT',c(paste0(c('NSC','IPC','N'),'_','mad.score'),'enh_typ','enh_class','names')])
  colnames(df) <- gsub('_mad.score','',colnames(df))
  df$sig <- ifelse(rowMins(as.matrix(res[res$enh_typ=='WT',c('NSC_pval.mad','IPC_pval.mad','N_pval.mad')]),na.rm=T)<=0.1,'yes','no')
  df_sig <- subset(df,sig=='yes')
  mat <- t(scale(t(df_sig[,1:3])))
  mat_k <- tglkmeans::TGL_kmeans(mat,k = k,reorder_func = hclust,parallel = F,seed=seed,max_iter = 10000,id_column = F)
  df_sig$cluster <- mat_k$cluster
  la <- rowAnnotation(foo = anno_block(labels = 1:k,labels_rot = 0,gp = gpar(fill = rev(colorpalette('brbg',k)),col = "black", border = "black", fontsize = 12)))
  
  hm <- Heatmap(df_sig[,1:3],name='MPRA signal',row_split = df_sig$cluster,cluster_columns = F,show_row_names = F,cluster_row_slices = F,cluster_rows = F,left_annotation = la,
                col = cols,column_names_rot = 0,column_names_centered = T,border=T,row_title = NULL,use_raster = F)
  pdf(paste0('figures/',out_f,'.pdf'),width=width,height=height)
  print(hm)
  dev.off()
  df$sig[df$sig=='no'] <- 'fuck'
  df$sig[df$sig=='fuck'&rowMins(as.matrix(res[res$enh_typ=='WT',c('NSC_pval.mad','IPC_pval.mad','N_pval.mad')]),na.rm=T)>0.1] <- 'no'
  df_Nosig <- df[df$sig=='no',]
  df_Nosig$cluster <- k+1
  df_test <- rbind(df_sig,df_Nosig)
  write.table(df_test,'/home/hpc/bonev/projects/mpra/mpra_3dram/results/horg_immunoMPRA_sigDF_clustered2.tsv',col.names=T,row.names=F,sep='\t',quote=F)
}

Figure7D2 <- function(df,pwms,cre_size,genome,mcparams,FDR.cutoff,enr.cutoff,out_f,height,width){
  df <- res[res$enh_typ=='WT'&(res$NSC_pval.mad<=.1|res$IPC_pval.mad<=.1),]
  peaks <- gsub('WT_','',df$names)
  peaks <- as.data.frame(stringr::str_split(paste0(peaks), pattern =':' , n = 2, simplify = TRUE))
  peaks$V2 <- as.numeric(as.character(peaks$V2))
  peaks[,3] <- peaks[,2]+1
  colnames(peaks) <- c('chrom','start','end')
  peaks <- makeGRangesFromDataFrame(peaks,keep.extra.columns = T)
  names(peaks) <- 1:length(peaks)
  peaks <- resize(peaks,cre_size,fix='center')          #Change coords to actual MPRA coords
  peakseqs <- getSeq(genome, peaks)
  #bins <- monaLisa::bin(x = factor(df_test$cluster), binmode = "equalN")
  se <- calcBinnedMotifEnr(seqs = peakseqs, bins = factor(df$cluster),
                           min.score=10,motifs = pwms,BPPARAM = mcparam)
  sel1 <- apply(assay(se, "negLog10Padj"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > FDR.cutoff
  sel2 <- apply(assay(se, "log2enr"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > enr.cutoff
  seSel <- se[sel1&sel2, ]
  #SimMatSel <- motifSimilarity(rowData(seSel)$motif.pfm,BPPARAM = mcparams)
  hcl <- hclust(dist(assay(seSel,"log2enr")), method = "single")
  pdf(paste0('figures/',out_f,'.pdf'),height=height,width=width)
  plotMotifHeatmaps(x = seSel, which.plots = c("log2enr"), width = 1.8,cluster=hcl,
                    show_dendrogram = TRUE, show_seqlogo = TRUE,
                    width.seqlogo = 1.2)
  dev.off()
}

Figure7D3 <- function(df,pwms,cre_size,genome,mcparams,FDR.cutoff,enr.cutoff,out_f,height,width){
  df <- res[res$enh_typ=='WT',]
  peaks <- gsub('WT_','',df$names)
  peaks <- as.data.frame(stringr::str_split(paste0(peaks), pattern =':' , n = 2, simplify = TRUE))
  peaks$V2 <- as.numeric(as.character(peaks$V2))
  peaks[,3] <- peaks[,2]+1
  colnames(peaks) <- c('chrom','start','end')
  peaks <- makeGRangesFromDataFrame(peaks)
  names(peaks) <- 1:length(peaks)
  peaks <- resize(peaks,cre_size,fix='center')          #Change coords to actual MPRA coords
  peakseqs <- getSeq(genome, peaks)
  suppressWarnings({
    hits <- findMotifHits(query = pwms, subject = peakseqs, min.score = 10.0,
                          BPPARAM = multicoreParam)
  })
  TFBSmatrix <- unclass(table(factor(seqnames(hits), levels = seqlevels(hits)),
                              factor(hits$pwmname, levels = name(pwms))))
  zero_TF <- colSums(TFBSmatrix) == 0
  TFBSmatrix <- TFBSmatrix[, !zero_TF]
  set.seed(79)
  se <- randLassoStabSel(x = TFBSmatrix, y = df$IPC_mad.score-df$NSC_mad.score,
                         cutoff = 0.8)
  pdf(paste0('figures/',out_f,'.pdf'),height=height,width=width)
  plotSelectionProb(se, directional = TRUE,selProbMin = 0.4,showSelProbMin = F,selProbMinPlot = 0.3,col = rep('black',3))
  dev.off()
  
  
}

Figure7E <- function(df,pwms,cre_size,genome,mcparams,motif1,motif2,y_lim=NULL,vals=c('NSC','IPC'),out_f,height,width,plot.sig=T){
  peaks <- gsub('WT_','',df$names)
  peaks <- as.data.frame(stringr::str_split(paste0(peaks), pattern =':' , n = 2, simplify = TRUE))
  peaks$V2 <- as.numeric(as.character(peaks$V2))
  peaks[,3] <- peaks[,2]+1
  colnames(peaks) <- c('chrom','start','end')
  peaks <- makeGRangesFromDataFrame(peaks)
  names(peaks) <- 1:length(peaks)
  peaks <- resize(peaks,cre_size,fix='center')          #Change coords to actual MPRA coords
  peakseqs <- getSeq(genome, peaks)
  suppressWarnings({
    hits <- findMotifHits(query = pwms, subject = peakseqs, min.score = 10.0,
                          BPPARAM = mcparams)
  })
  TFBSmatrix <- unclass(table(factor(seqnames(hits), levels = seqlevels(hits)),
                              factor(hits$pwmname, levels = name(pwms))))
  mat <- df[,c(paste0(vals,'_mad.score'),'enh_typ')]
  colnames(mat) <- c(vals,'enh_type')
  mat$motif1 <- TFBSmatrix[,motif1]
  mat$motif2 <- TFBSmatrix[,motif2]
  mat$motif1_motif2 <- mat$motif1+mat$motif2
  mat$motif1_motif2[mat$motif1==0|mat$motif2==0] <- 0
  mat$motif <- factor(case_when(mat$motif1_motif2>=2 ~ 'Both',
                         mat$motif1>=1 ~ motif1,
                         mat$motif2>=1 ~ motif2,
                         TRUE ~ 'None'),levels=c('None',motif1,motif2,'Both'))
  mat$NSC_sig <- df$NSC_mad.score<=0.1
  mat$IPC_sig <- df$IPC_mad.score<=0.1
  
  mat_m <- melt(mat,measure.vars = vals)  
  # mat_m$motif1 <- factor(mat_m$motif1,levels=c('scr',0:max(as.numeric(mat$motif1),na.rm=T)))
  # mat_m$motif2 <- factor(mat_m$motif2,levels=c('scr',0:max(as.numeric(mat$motif2),na.rm=T)))
  # mat_m$motif1_motif2 <- factor(mat_m$motif1_motif2,levels=c('scr',0:max(as.numeric(mat$motif1_motif2),na.rm=T)))
  
  if(is.null(y_lim)){
    y_lim <- c(quantile(mat_m$value,0.01)*1.5,quantile(mat_m$value,0.95)*1.5)
  }
  p <- ggplot(mat_m,aes(x=motif,y=value,fill=motif)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8,outlier.shape = NA)
  p <- p + scale_fill_npg() + xlab('') + ylab('MPRA signal') + theme(legend.position = "none") + coord_cartesian(ylim=y_lim) + facet_wrap(~variable)
  if(plot.sig){  
  p <- p + stat_compare_means(comparisons = list(c('None',motif1),c('None',motif2),c('None','Both')),label = "p.format",method='wilcox',paired = F,label.y=label.y,tip.length = rep(0.005,8))
  }
  # 
  # p1 <- ggplot(mat_m,aes(x=motif1,y=value,fill=motif1)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8,outlier.shape = NA)
  # p1 <- p1 + scale_fill_npg() + xlab('') + ylab('MPRA signal') + theme(legend.position = "none") + coord_cartesian(ylim=y_lim) + facet_wrap(~variable)
  # p1 <- p1 + stat_compare_means(comparisons = list(c('scr','0'),c('0','1'),c('1',2)),label = "p.format",method='wilcox',paired = F,label.y=label.y,tip.length = 0.005)
  # p2 <- ggplot(mat_m,aes(x=motif2,y=value,fill=motif2)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8,outlier.shape = NA)
  # p2 <- p2 + scale_fill_npg() + xlab('') + ylab('MPRA signal') + theme(legend.position = "none") + coord_cartesian(ylim=y_lim) + facet_wrap(~variable)
  # p2 <- p2 + stat_compare_means(comparisons = list(c('scr','0'),c('0','1'),c('1',2)),label = "p.format",method='wilcox',paired = F,label.y=label.y,tip.length = 0.005)
  # p3 <- ggplot(mat_m,aes(x=motif1_motif2,y=value,fill=motif1_motif2)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8,outlier.shape = NA)
  # p3 <- p3 + scale_fill_npg() + xlab('') + ylab('MPRA signal') + theme(legend.position = "none") + coord_cartesian(ylim=y_lim) + facet_wrap(~variable)
  # p3 <- p3 + stat_compare_means(comparisons = list(c('scr','0'),c('0','2'),c('2',3)),label = "p.format",method='wilcox',paired = F,label.y=label.y,tip.length = 0.005)
  # 
  pdf(paste0('figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
  return(mat)
}


Figure7E <- function(res,pwms,cre_size,genome,mcparams,FDR.cutoff,features,anno.size=4,out_f,height,width,cols){
  df <- res[(res$NSC_pval.mad<=0.1|res$IPC_pval.mad<=0.1)&res$enh_typ=='WT',]
  peaks <- gsub('WT_','',df$names)
  peaks <- as.data.frame(stringr::str_split(paste0(peaks), pattern =':' , n = 2, simplify = TRUE))
  peaks$V2 <- as.numeric(as.character(peaks$V2))
  peaks[,3] <- peaks[,2]+1
  colnames(peaks) <- c('chrom','start','end')
  peaks <- makeGRangesFromDataFrame(peaks)
  names(peaks) <- 1:length(peaks)
  peaks <- resize(peaks,cre_size,fix='center')          #Change coords to actual MPRA coords
  peakseqs <- getSeq(genome, peaks)
  motif_ix <- matchMotifs(pwms, peakseqs,p.cutoff=5e-05) 
  TFBSmatrix <- motifMatches(motif_ix)
  colnames(TFBSmatrix) <- as.vector(name(pwms))
  row.names(TFBSmatrix) <- df$names
  tf_mat <- adply(TFBSmatrix[,colSums(TFBSmatrix)>=10],2,function(x){
    mat <- df[which(x>=1),grep('mad.score',colnames(df))]
    if(nrow(mat)>=10){
    return(c(colMedians(as.matrix(mat[,1:3]),na.rm=T),wilcox.test(mat[,1],mat[,2],paired=T)$p.value,nrow(mat)))
    }
  })
  colnames(tf_mat) <- c('TF','NSC MPRA signal','IPC MPRA signal','N MPRA signal','p.value','N')
  
  permut_df <- as.data.frame(matrix(NA,nrow=1000,ncol = 2))
  for (i in 1:1000){
    idx <- sample(1:nrow(df),50)
    permut_df[i,] <- colMedians(as.matrix(df[idx,c('NSC_mad.score','IPC_mad.score')]),na.rm=T)
  }
  nsc_score_q <- quantile((permut_df[,2]-permut_df[,1]),0.05)
  ipc_score_q <- quantile((permut_df[,2]-permut_df[,1]),0.95)
  
 tf_mat$perm <- 'ns' 
 tf_mat$perm[tf_mat$`IPC MPRA signal`-tf_mat$`NSC MPRA signal`>=ipc_score_q] <- 'ipc'
 tf_mat$perm[tf_mat$`IPC MPRA signal`-tf_mat$`NSC MPRA signal`<=nsc_score_q] <- 'nsc'
 
 
 p <-  ggplot(tf_mat,aes(x=`NSC MPRA signal`,y=`IPC MPRA signal`)) + geom_point(data = tf_mat[tf_mat$p.value>FDR.cutoff,],col='grey80') + geom_point(data = tf_mat[tf_mat$p.value<=FDR.cutoff&(tf_mat$`NSC MPRA signal`-tf_mat$`IPC MPRA signal`>0),],col=cols[1]) + geom_point(data = tf_mat[tf_mat$p.value<=FDR.cutoff&(tf_mat$`NSC MPRA signal`-tf_mat$`IPC MPRA signal`<0),],col=cols[2]) 
 p <- p  + geom_abline(slope = 1,intercept = 0,linetype = 2,col='black') + coord_cartesian(xlim=c(0,7),ylim=c(0,7))
 p <- p + geom_text_repel(
   data = tf_mat, size = anno.size,box.padding = 0.5,segment.alpha = 0.5, min.segment.length = 0,max.iter = 20000,
   aes(color=NULL,label=ifelse(TF%in%features, as.character(TF), "")),force=10)
 
 pdf(paste0('figures/',out_f,'.pdf'),height=height,width=width)
 print(p)
 dev.off()  
 return(TFBSmatrix)
}

Figure7F <- function(res,TFBSmatrix,motif,y_lim=NULL,cols,label.y=c(20,22,18),out_f,height,width){
  #df <- res[(res$NSC_pval.mad<=0.1|res$IPC_pval.mad<=0.1)&res$enh_typ=='WT',]
  df <- res[res$enh_typ=='WT',]
  df_scr <- res[res$class=='scr',]
  mat <- df[TFBSmatrix[,motif]>=1,grep('mad.score',colnames(df))]
 # scr_mat <- df_scr[,grep('mad.score',colnames(df_scr))]
  colnames(mat) <- c('NSC','IPC','N')
#  colnames(scr_mat) <- c('NSC','IPC','N')
  mat_m <- melt(mat)
  if(is.null(y_lim)){
    y_lim <- c(quantile(mat_m[,2],0.01)*1.5,quantile(mat_m[,2],0.95))
  }
  p1 <- ggplot(mat_m,aes(x=variable,y=value,fill=variable)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8,outlier.shape = NA) 
  p1 <- p1 + scale_fill_manual(values=cols) + xlab('') + ylab('MPRA signal') + theme(legend.position = "none") + coord_cartesian(ylim=y_lim)
  p1 <- p1 + stat_compare_means(comparisons = list(c('NSC','IPC'),c('NSC','N'),c('IPC','N')),label = "p.format",method='wilcox',paired = T,label.y=label.y,tip.length = 0.005)
  pdf(paste0('figures/',out_f,'.pdf'),height=height,width=width)
  print(p1)
  dev.off()
}

Figure7G <- function(res,y_lim=NULL,cols=c('grey80','red','blue'),out_f,height,width){
  mat <- res[res$enh_typ=='scr'|res$enh_typ=='UCON'|res$enh_typ=='MER130',c(paste0(c('NSC','IPC','N'),'_mad.score'),'enh_typ','names')]
  mat_sig <- res[res$enh_typ=='UCON'|res$enh_typ=='MER130',c(paste0(c('NSC','IPC','N'),'_pval.mad'),'enh_typ','names')]
  colnames(mat) <- c('NSC','IPC','N','Class','Name')
  mat$Class <- factor(mat$Class,levels=c('scr','UCON','UCON31','MER130'))
  mat$Class[mat$Class=='UCON'] <- 'UCON31'
  mat$Class <- droplevels(mat$Class)
  mat_m <- melt(mat,id.vars = c('Class','Name'))
  if(is.null(y_lim)){
    y_lim <- c(quantile(mat_m$value,0.01)*2,quantile(mat_m$value,0.95))
  }
  p1 <- ggplot(mat_m,aes(x=variable,y=value,fill=Class)) + geom_boxplot(outlier.size=1,show.legend = T,width=0.8,outlier.shape = NA) 
  p1 <- p1 + scale_fill_manual(name='',values=cols) + xlab('') + ylab('MPRA signal') + theme(legend.position = c(.75,.9)) + coord_cartesian(ylim=y_lim)
  sig_df <- ddply(mat_sig,.(enh_typ),function(x){
    test <- round(colSums(x[,1:3]<=0.1)/nrow(x)*100,2)
    return(test)
  })
  colnames(sig_df) <- c('Class','NSC_sig_%','IPC_sig_%','N_sig_%')
  #p1 <- p1 + stat_compare_means(comparisons = list(c('scr','NSC'),c('scr','IPC'),c('NSC','IPC')),,label = "p.format",method='wilcox',paired = T,label.y=label.y,tip.length = 0.005)
  pdf(paste0('figures/',out_f,'.pdf'),height=height,width=width)
  print(p1)
  dev.off()
  return(sig_df)
}

Figure7H <- function(df,mut_motif,y_lim,cols=c("#3B4992FF","#EE0000FF"),label.y=25,out_f,height,width){
  
  df_mut <- res[res$enh_typ==mut_motif,]
  df$coord <- stringr::str_split(df$names, pattern ='WT_' , n = 2, simplify = TRUE)[,2]
  df_mut$coord <- stringr::str_split(df_mut$names, pattern ='WT_' , n = 2, simplify = TRUE)[,2]
  df_mut <- df_mut[df_mut$coord%in%df$coord,]
  df <- df[df$coord%in%df_mut$coord,]
  df_mut <- df_mut[match(df$coord,df_mut$coord),]
  mat <- rbind(df,df_mut)[,c('NSC_mad.score','IPC_mad.score','N_mad.score','class','coord')]
  colnames(mat) <- c('NSC','IPC','N','Class','Coord')
  mat$Class <- factor(mat$Class,levels=c('WT','mut'))
  levels(mat$Class) <- c('WT','MUT')
  mat_m <- melt(mat,id.vars = c('Class','Coord'))
  if(is.null(y_lim)){
    y_lim <- c(quantile(mat_m[,4],0.01)*1.5,quantile(mat_m[,4],0.95)*1.5)
  }
  p1 <- ggplot(mat_m,aes(x=Class,y=value,fill=Class)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8,outlier.shape = NA)
  p1 <- p1 + scale_fill_manual(values=cols) + xlab('') + ylab('MPRA signal') + theme(legend.position = "none") + coord_cartesian(ylim=y_lim) + facet_wrap(~variable)
  p1 <- p1 + stat_compare_means(comparisons = list(c('WT','MUT')),label = "p.format",method='wilcox',paired = T,label.y=label.y,tip.length = 0.005)
  pdf(paste0('figures/',out_f,'.pdf'),height=height,width=width)
  print(p1)
  dev.off()
  message(nrow(df))
}


Figure7M <- function(facs_res_f,cols=c('#aaa9ad','#aaa9ad')) {
  df <- read.table(facs_res_f,header=T,sep='\t')
  df <- melt(df,varnames = c('Species'),measure.vars = 'Positive_Cells')
  df_summ <- data_summary(data=df, varname="value",groupnames=c("Species"))
  p <- ggplot(df_summ, aes(x=Species, y=value, fill=Species)) + scale_fill_grey(start = 0.35,end = 1) +
    geom_bar(stat="identity", color="black",position=position_dodge()) + xlab('') + ylab('') +
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.9)) + ggtitle('')
  p <- p + theme(legend.position='none',plot.title = element_text(hjust = 0.5)) + geom_point(data = df,position=position_jitter(w = 0.2, h = 0),size=1) + annotate("text", x = c(1,2), y = c(17,3), 
                                                                                                                                                                    label = round(df_summ$value,2))
  return(p)
}

### Plot Figures

Figure7B(MPRA_data_f='/home/fnoack/projects/MPRA/MPRAflow/outs/MPRA_3DRAM/',
          cols=c(cell_colors,'grey50'),
          celltypes<-c('Pax6','Tbr2','PN'),names=c('RGC','IPC','N'),
          out_f="figures/Figure7B.pdf",height=3.5,width=2.5)

Figure7C(res=res,metric = 'mad.score',out_f='Figure7C',ylim = c(-3,23),height=5,width=6)
Figure7D(res=res,k=5,cols=colorRamp2(c(seq(0,20,length.out=15)), rev(colorpalette('reds',15))),out_f='Figure7D',width=5,height=6,seed=123)
df <- Figure7E(df=res[res$enh_typ=='WT'&(res$NSC_pval.mad<=.1|res$IPC_pval.mad<=.1),],
               mcparams = MulticoreParam(10),motif1='NEUROG2(var.2)',motif2 = 'EOMES',vals=c('IPC'),y_lim = c(-2,16),plot.sig = F,
               pwms=readRDS(paste0(main_dir,'data/hg38/pwms_FPKM1.RDS')),cre_size=266,genome=BSgenome.Hsapiens.UCSC.hg38,
               out_f='Figure7E_1',height=4,width=3.5)
df <- Figure7E(df=res[res$enh_typ=='WT'&(res$NSC_pval.mad<=.1|res$IPC_pval.mad<=.1),],
         mcparams = MulticoreParam(10),motif1='SOX2',motif2 = 'LHX2',vals=c('NSC'),y_lim = NULL,plot.sig = F,
         pwms=readRDS(paste0(main_dir,'data/hg38/pwms_FPKM1.RDS')),cre_size=266,genome=BSgenome.Hsapiens.UCSC.hg38,
         out_f='Figure7E_2_TEST',height=4,width=3.5)

df <- Figure7E(df=res[res$enh_typ=='WT'&(res$NSC_pval.mad<=.1|res$IPC_pval.mad<=.1),],
               mcparams = MulticoreParam(10),motif1='NEUROG2(var.2)',motif2 = 'EOMES',vals=c('IPC'),y_lim = c(-2,16),plot.sig = F,
               pwms=readRDS(paste0(main_dir,'data/hg38/pwms_FPKM1.RDS')),cre_size=266,genome=BSgenome.Hsapiens.UCSC.hg38,
               out_f='Figure7E_bla',height=4,width=3.5)

df <- Figure7E(df=res[res$enh_typ=='WT'&(res$NSC_pval.mad<=.1|res$IPC_pval.mad<=.1|res$N_pval.mad<=.1),],
               mcparams = MulticoreParam(10),motif1='NEUROD1',motif2 = 'POU3F3',vals=c('IPC'),y_lim = c(-2,16),plot.sig = F,
               pwms=readRDS(paste0(main_dir,'data/hg38/pwms_FPKM1.RDS')),cre_size=266,genome=BSgenome.Hsapiens.UCSC.hg38,
               out_f='Figure7E_POU3F3',height=4,width=3.5)

Figure7F(res=res,TFBSmatrix=TFBSmatrix,motif='NEUROG2(var.2)',y_lim=c(-3,25),cols=c(cell_colors,'grey80'),label.y=c(20,22,18),out_f='Figure7F_Neurog2',height=4,width=3.7)
Figure7F(res=res,TFBSmatrix=TFBSmatrix,motif='HES1',y_lim=c(-3,25),cols=c(cell_colors,'grey80'),label.y=c(20,22,18),out_f='Figure7F_HES1',height=4,width=3.7)
Figure7G(res=res,cols=c('grey80','red','blue'),out_f='Figure7G_UCON',height=6,width=5,y_lim = c(-2.7,5))
Figure7H(df=res[(res$IPC_pval.mad<=0.1)&res$enh_typ=='WT',],
         mut_motif='NEUROG2',y_lim=c(-2,25),cols=c("#3B4992FF","#EE0000FF"),label.y=24,out_f='Figure7H_Neurog2',height=5,width=4.5)
Figure7H(df=res[(res$IPC_pval.mad<=0.1)&res$enh_typ=='WT',],
         mut_motif='EOMESmut',y_lim=NULL,cols=c("#3B4992FF","#EE0000FF"),label.y=23,out_f='Figure7H_EOMES',height=5,width=4.5)
Figure7H(df=res[(res$NSC_pval.mad<=0.1)&res$enh_typ=='WT',],
         mut_motif='SOX2mut',y_lim=NULL,cols=c("#3B4992FF","#EE0000FF"),label.y=12,out_f='Figure7H_SOX2',height=5,width=4.5)
Figure7H(df=res[(res$NSC_pval.mad<=0.1)&res$enh_typ=='WT',],
         mut_motif='LHX2mut',y_lim=NULL,cols=c("#3B4992FF","#EE0000FF"),label.y=12,out_f='Figure7H_LHX2',height=5,width=4.5)

Figure7H(df=res[(res$IPC_pval.mad<=0.1)&res$enh_typ%in%c('MER130'),],
         mut_motif='NEUROG2',y_lim=c(-2,30),cols=c("#3B4992FF","#EE0000FF"),label.y=29,out_f='Figure7H_MER130_NEUROG2',height=5,width=4.5)
Figure7H(df=res[(res$IPC_pval.mad<=0.1)&res$enh_typ%in%c('UCON'),],
         mut_motif='NEUROG2',y_lim=c(-2,30),cols=c("#3B4992FF","#EE0000FF"),label.y=29,out_f='Figure7H_UCON_NEUROG2',height=5,width=4.5)

source(paste0(main_dir,'results/hg38/hic/config.R'))

#####
df <- res[(res$IPC_pval.mad<=0.1)&res$enh_typ=='WT',]
peaks <- gsub('WT_','',df$names)
peaks <- as.data.frame(stringr::str_split(paste0(peaks), pattern =':' , n = 2, simplify = TRUE))
peaks$V2 <- as.numeric(as.character(peaks$V2))
peaks[,3] <- peaks[,2]+1
colnames(peaks) <- c('chrom','start','end')
mpra_coords <- intervals.normalize(peaks,1000)      #For Visualisation in big regions
mpra_tracks <- c("mpra.horgD45_WT_NSC","mpra.horgD45_WT_IPC","mpra.horgD45_WT_N")

anno <- gintervals('chr13',66781409,66781410)
anno <- intervals.normalize(anno,1000)   #Makes all intervals 1kb long
misha_meth_tracks <- c("methylation.3DRAM_Pax6_CpG_merged_10x","methylation.3DRAM_Tbr2_CpG_merged_10x","methylation.3DRAM_Pax6_GpC_merged_10x","methylation.3DRAM_Tbr2_GpC_merged_10x") 

plotMisha(main_f=main_f,targetGene='PCDH9',outDir='figures/',out_f='Figure7_PCDH9_MPRA',upstream=6e5,downstream=1e5,
          chipNames=c('','','',''),window_scale=1.8,chipRes=20,pointCEX=2,conditions=score_tracks,binSize=5e3,radius=2e4,annIntervals=anno,
          chipTracksToExtract=c(misha_meth_tracks[grep('GpC',misha_meth_tracks)],'rnaseq.Pax6_HMGU1_Day45_RNA_merged','rnaseq.Tbr2_HMGU1_Day45_RNA_merged'),chipColors=rep(cell_colors,2), 
          mpra_coords=mpra_coords,mpra_tracks=mpra_tracks,mpra_colors=colorpalette('matlablike2',11),mpraNames=c('','',''),mpra_ylim=c(0,160),
          chipYlim=matrix(c(0,100,0,100,0,13.5,0,13.5),ncol = 2,byrow = T),img_factor=1.2,
          plotOrder=list(scores=TRUE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,anno=FALSE,MPRA=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.2, VP=1.5, loops=2.2, rna=0.6, chip=0.5,MPRA=0.3,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.2))

df <- res[grep('chr13:66781409',res$names),]
df <- df[,c('NSC_mad.score','IPC_mad.score','N_mad.score')]
colnames(df) <- c('NSC','IPC','PN')
df$enh_type <- factor(c('WT','NEUROG2mut','EOMESmut'),levels=c('WT','NEUROG2mut','EOMESmut'))
mat <- melt(df,id.vars = 'enh_type')
p <- ggplot(mat,aes(x=enh_type,y=value,fill=enh_type)) + geom_bar(stat="identity",col='black') + facet_wrap(~variable,ncol=1)
p <- p + scale_fill_manual(values=c("#3B4992FF","#008B45FF","#EE0000FF")) + xlab('') + ylab('MPRA signal') + theme(legend.position = "none") 
pdf('figures/Figure7I.pdf',height=6,width=3)
p
dev.off()


p <- Figure7M(paste0(main_dir,'results/hg38/FBXO32_FACS.tsv'))
pdf('figures/Figure7M.pdf',height=4,width=2,useDingbats = F)
print(p)
dev.off()

# library(Signac)
# NSC_motifs <- enrichedMotifs(peaks='../../results/hg38/mpra/cluster1.bed',bg_peaks = '../../results/hg38/mpra/ALL_WT.bed',genome = 'hg38',genome_BS = BSgenome.Hsapiens.UCSC.hg38,pwm=readRDS(paste0(main_dir,'data/hg38/pwms_FPKM1.RDS')),
#                              features=c('Neurog2(var.2)','Neurod1','Neurod2','Eomes','Ctcf','Tcf12(var.2)','Tcf4','Nrf1','Fos::jun','Lhx2','Sox2','Tgif2','Nfia','Rfx4','Zeb1','Zbtb18','Tead2','Mybl2','Emx2','Tcf7l2','Irf1','Stat1::stat2','Pou3f2','Rest'),
#                              cols=c("blue",'grey80','red'),logFC=0.2,logP=1,point.size=4,anno.size=6)
#         FDR.cutoff=0.05,cols = cell_colors,anno.size = 4,
#         features=c('NEUROG2(var.2)','NEUROD1','BHLHE40','PAX6','POU3F2','TCF7L2','DMRT3','NFIA','NFIX','FOSB::JUN','ZNF135','SOX2','TEAD2','RBPJ','RXRA::VDR','CTCF','NEUROD2','NFIA','TGIF2','LHX2'),

