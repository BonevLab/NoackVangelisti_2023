source('scripts/config.R')
source('scripts/aux_functions.R')
source('scripts/plot_functions.R')

source(paste0(main_dir,'results/mm10/hic/config.R'))
source(paste0(main_f,'scripts/main_functions.R'))
source(paste0(main_f,'scripts/aux_functions.R'))
source(paste0(main_f,'scripts/plot_functions.R'))
source(paste0(main_f,'scripts/temp_functions.R'))

########################
###Functions###
########################

#FigS2A

FigureS2A <-function(DRAM_RNA,BB_RNA,cormethod,main='',filter_gtf='protein_coding',gtf_f,fig_f,height,width) {
  RAM_fpkm<-read.table(file =DRAM_RNA, header =T , sep='\t')
  bonev_fpkm<-read.table(file =BB_RNA, header =T , sep='\t')
  df <-merge(x=RAM_fpkm, y=bonev_fpkm, by='gene_name')
  colnames(df)<-c('gene','rep1','rep2','rep3','3DRAM-seq','RNA-seq mESC','RNA-seq NPC','RNA-seq CN','RNA-seq ncxNPC','RNA-seq ncxCN')
  df[,2:ncol(df)] <- log10(df[,2:ncol(df)]+1)
  if(!is.null(filter_gtf)){
    gtf <- rtracklayer::import(gtf_f)
    filtered_genes <- unique(gtf$gene_name[gtf$gene_type==filter_gtf])
    df <- df[df$gene%in%filtered_genes,]
  }
  jpeg(fig_f,height=height,width=width,res = 300)   
  print(heatpairs(as.matrix(df[,5:8]),method=cormethod,main = main))
  dev.off()
} 

FigureS2E <- function(){
  chromHHM<-read.table('/home/fnoack/temp/RAM_rebottle/ChromHMM/mESC_mm10_chromHMM.bed', header=F)
  states<-unique(chromHHM$V4)
  tracks <- c("methylation.ES_3DRAM_GpC_10x")
  df<-data.frame()
  for (i in 2:length(states)) {
    peaks<-subset(chromHHM,chromHHM$V4 == states[i])
    peaks<-peaks[,1:3]
    colnames(peaks)<-colnames(gintervals.all())
    res <- misha_extract(tracks=tracks,regions=peaks,iterator=peaks) 
    res$intervalID<-states[i]
    res <- res[,4:5]
    df<-rbind(df,res)
  }
  colnames(df)<-c('GpC_meth','ChromHMM')
  df$ChromHMM<-as.character(df$ChromHMM)
  df[grepl("E11",df$ChromHMM),2]<-'Weak Enhancer'
  df[grepl("E10",df$ChromHMM),2]<-'Transcription Elongation'
  df[grepl("E1",df$ChromHMM),2]<-'CTCF'
  df[grepl("E2",df$ChromHMM),2]<-'Intergenic'
  df[grepl("E3",df$ChromHMM),2]<-'Heterochromatin'
  df[grepl("E4",df$ChromHMM),2]<-'Enhancer'
  df[grepl("E5",df$ChromHMM),2]<-'Repressed Chromatin'
  df[grepl("E6",df$ChromHMM),2]<-'Bivalent Promoter'
  df[grepl("E7",df$ChromHMM),2]<-'Active Promoter'
  df[grepl("E8",df$ChromHMM),2]<-'Strong Enhancer'
  df[grepl("E9",df$ChromHMM),2]<-'Transcriptional Transition'
  p<-ggplot(df, aes(x=ChromHMM, y=GpC_meth)) + xlab("") + 
    geom_boxplot(outlier.size=1,outlier.shape=NA,show.legend = T,width=0.8)+theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+ylab('GpC Methylation (%)')
  return(p)
}

########################
###Call Functions###
########################

#FigS2A -> Correlation of RNA-seq data with 3DRAM seq data
FigureS2A(DRAM_RNA='/home/hpc/bonev/projects/ram/data/mm10/RAM_RNA/FPKM_means.txt',
       BB_RNA='/home/hpc/bonev/projects/ram/data/mm10/RAM_RNA/Bonev_fpkm_means.tsv',cormethod='spearman',
       main='',filter_gtf='protein_coding',gtf_f='/home/hpc/bonev/annotations/mm10/mm10.gencode.gtf',
       fig_f='figures/FigureS2A.jpg',height=2400,width=2400)


#FigS2B Pearson correlation between datasets -> use deeptools to calculate the values  
tracks <- c("methylation.ES_3DRAM_GpC_10x","atac.ES_GSE113592","atac.ES_DNAseI_ENCSR000CMW")
res <- misha_extract(tracks=tracks,regions=gintervals.all(),iterator=10000,blacklist = mm10_bl)
colnames(res)[4:6] <- c('3DRAM-seq','ATAC','DNAseI')
cor_res <- cor(res[,4:6],use = 'complete.obs',method='spearman')
   
jpeg('figures/FigureS2B.jpg',height=1600,width=1600,res = 300)           #pdf file becomes too big
heatpairs(as.matrix(res[,4:6]),method='spearman',main = '',xlim = c(0,50),ylim=c(0,50))
dev.off()

#pheatmap(cor_res,color = colorRampPalette(heatmap_colors)(101), breaks = seq(-1,1,length=100),display_numbers = TRUE, number_color = "black",fontsize_number = 25,
#         border_color = "black",cluster_rows = F,cluster_cols = F,angle_col = 0,fontsize = 14,filename ='figures/FigureS2B.pdf',width=7,height=6)

#FigS2C VennDiagram overlaps 
require(ChIPpeakAnno)
gNOME_peaks <- toGRanges(paste0(main_dir,'data/mm10/beds/gNOMEpeaks/gNOME_peaks.bed'), format="BED", header=FALSE) 
ATAC_peaks <- toGRanges(paste0(main_dir,'data/mm10/ATAC_seq/ATAC_peaks.bed'), format="BED", header=FALSE) 
DNAse_peaks <- toGRanges(paste0(main_dir,'data/mm10/DNase_ENCSR000CMW/DNA_E14_peaks_encode.bed'), format="BED", header=FALSE)

ol <- findOverlapsOfPeaks(ATAC_peaks,DNAse_peaks,gNOME_peaks)
pdf(file='figures/FigureS2C.pdf',height=6,width=6)
print(makeVennDiagram(ol, NameOfPeaks=c('ATAC','DHS','3DRAM'),
                     # fill=c("#E69F00", "#56B4E9",'grey57'), # circle fill color
                     fill=c('blue','green','red'),
                      col=c('black'), #circle border color
                      cat.col=c("black"),scale=T,cat.cex=2,cex=2)
)
dev.off()

#FigS2D
tracks <- c("methylation.ES_WGBS_10x","methylation.ES_3DRAM_CpG_10x","methylation.ES_MethylHiC_CpG_merged_10x","methylation.ES_Methyl3C_CpG_merged_10x")
res <- misha_extract(tracks=tracks,regions=gintervals.all(),iterator=10000,blacklist = mm10_bl)
colnames(res)[4:7] <- c('WGBS','3DRAM-seq','Methyl-HiC','Methyl-3C')
cor_res <- cor(res[,4:7],use = 'complete.obs',method='pearson')

#pheatmap(cor_res,color = colorRampPalette(heatmap_colors)(101), breaks = seq(-1,1,length=100),display_numbers = TRUE, number_color = "black",fontsize_number = 25,
#         border_color = "black",cluster_rows = F,cluster_cols = F,angle_col = 0,fontsize = 14,filename ='figures/FigureS2D.pdf',width=7,height=6)

jpeg('figures/FigureS2D.jpg',height=1600,width=1600,res = 300)           #pdf file becomes too big
heatpairs(as.matrix(res[,4:7]),method='pearson',main = '')
dev.off()

#Fig2E Heatmap signals across GpC peaks
tracks <- c(paste0(main_dir,'data/mm10/RAM_bigwigs/3DRAM_ES_250k_merged_CpG_cov10x.bw'),paste0(main_dir,'data/mm10/Methyl3C_Ecker/3C_mESC_ecker_CpG_cov10x.bw'),paste0(main_dir,'data/mm10/MethylHiC_Ren/MethylHiC_merged_CpG_cov10x.bw'),paste0(main_dir,'data/mm10/WGBS/ES_WGBS_CpG_cov10x.bw'))
peaks <- paste0(main_dir,'data/mm10/beds/ES_CTCF_top30k_motifs.bed')
res <- getPlotSetArray(tracks=tracks,features=peaks,refgenome='mm10',type = 'mf',add_heatmap=T,xmin=1000,xmax=1000,bin = 50,ignore_strand = T)
jpeg('figures/FigureS2E.jpg',width=3200,height=1600,res=300)
plotHeatmap(res,labels = c('3DRAM-seq','Methyl-3C','Methyl-HiC','WGBS'),sortrows = 'increasing',clspace = rev(colorpalette('blues',12)),
            clstmethod=FALSE,raster = T,autoscale = F,zmin = 0,zmax = 100,indi=F)
dev.off()


#FigS2E 
p <- FigureS2E()
pdf('figures/FigureS2E.pdf',width=5,height=5,useDingbats = F)
print(p+coord_cartesian(ylim=c(0,65)))
dev.off()



#FigS2F CisDecay of different datasets generated with the HiC scripts 

file_name <- 'cisDecay_combined'
run_cisDecay(tracks=all_tracks,cellsToPlot=c("ES_3DRAM_250k","ES_Methylhic_Ecker","ES_Methylhic_Ren","ES_unsorted2017"),
                    path=paste0(main_f,'analysis/cis_decay/'),log_base=10,cells=c("ES_3DRAM_250k","ES_Methylhic_Ecker","ES_Methylhic_Ren","ES_unsorted2017"),
                       file_f='cisDecay_combined',minDist=3e3,maxDist=1e8,
                       colors=c("red1",'green1','royalblue1','black'),width=8,height=8,
                       labels=c('3DRAM-seq','Methyl-HiC','Methyl-3C','Hi-C'),alpha=0.3,out_f='figures/FigureS2F.pdf')  

#FigureS2G

hic_files <- c('/home/fnoack/projects/3D_RAMseq_250k/merged_replicates/merged_HiC/3DRAM_ES_250k.hic','/home/fnoack/projects/3D_RAMseq_250k/3C_mESC_eckert/3C_Ecker_frag.hic','/home/fnoack/projects/3D_RAMseq_250k/3C_mESC_ren/3C_mESC_ren_frag.hic',)
res_df <- correlate_hic(hic_files,1e4)            #in aux_functions.R
dimnames(res_df) <- list(paste0('rep',1:3),paste0('rep',1:3))
pdf('figures/FigureS1G.pdf',height=4,width=4)
corrplot(res_df, method="color", col=colorRampPalette(heatmap_colors)(100), order="hclust",number.digits = 3, 
         addCoef.col = "black",cl.lim=c(-1,1),cl.length=3,addgrid.col = 'black',
         tl.col="black", tl.srt=0,number.font=1,bg='transparent'
)
dev.off()

  
## Figure 3A
intervals_to_test <- c("ES_Nanog.narrowPeak","ES_Sox2.bed","ES_Nrf1.bed","ES_Klf4.bed","ES_Rest.bed","ES_CTCF_top30k_motifs.bed","ES_activeTss.bed","ES_inactiveTss.bed","mES_ATAC_top50k_distal.bed","ES_CTCF_top30k_shuffled.bed")
res <- matrix(rep(NA,length(intervals_to_test)*3),ncol=3)
row.names(res) <- c('Nanog','Nrf1','Klf4','Rest','Sox2','Ctcf','activeTSS','inactiveTSS','CREs','random')
colnames(res) <- c('HiCvsCpG','HiCvsGpC','CpGvsGpC')
if (!file.exists(paste0(main_dir,'results/mm10/TF_cors.tsv'))){
  for (i in seq_along(intervals_to_test)){
    res[i,] <- extract_misha_pairs(bed1=paste0(main_dir,"data/mm10/beds/",intervals_to_test[i]),
                                   bed2=paste0(main_dir,"data/mm10/beds/",intervals_to_test[i]),
                                   interval_window=200,min_dist=5e4,max_dist=2e6,expand=c(-5000,5000),
                                   domains_f="hic.ES_Bonev2017.ins_250_domains_expanded",score_tracks=score_tracks[1],
                                   hic_names='HiC',other_tracks=c("methylation.ES_3DRAM_CpG_10x","methylation.ES_3DRAM_GpC_10x"),
                                   other_names=c('CpG','GpC'),intra_only = F)
  }
  write.table(res,paste0(main_dir,'results/mm10/TF_cors.tsv'),col.names=T,row.names=T,quote=F,sep='\t')
} else {
  res <- read.table(paste0(main_dir,'results/mm10/TF_cors.tsv'),header=T)
}

res$Name <- factor(res$Name,levels=c('Nanog','Nrf1','Klf4','Rest','Sox2','Ctcf','activeTSS','inactiveTSS','CREs','random'))

#Figure S2G
pdf('figures/FigureS2G.pdf',width=6,height=6)
p <- ggplot(res,aes(x=CpGvsGpC,y=Name)) + geom_bar(stat = "identity",fill='#a62c2b',color='black') 
p <- p + scale_y_discrete(limits=rev(levels(res$Name))) + xlab('Correlation') + ylab('') + xlim(c(-0.75,0.25)) + ggtitle('CpG Methylation vs GpC accessibility')
p <- p + geom_vline(xintercept = 0,linetype=2) 
print(p)
dev.off()

#Figure S2H
pdf('figures/FigureS2H.pdf',width=6,height=6)
p <- ggplot(res,aes(x=HiCvsCpG,y=Name)) + geom_bar(stat = "identity",fill='#a62c2b',color='black') 
p <- p + scale_y_discrete(limits=rev(levels(res$Name))) + xlab('Correlation') + ylab('') + xlim(c(-0.25,0.25)) + ggtitle('HiC score vs CpG Methylation')
p <- p + geom_vline(xintercept = 0,linetype=2) 
print(p)
dev.off()

#Figure S2I
pdf('figures/FigureS2I.pdf',width=6,height=6)
p <- ggplot(res,aes(x=HiCvsGpC,y=Name)) + geom_bar(stat = "identity",fill='#a62c2b',color='black') 
p <- p + scale_y_discrete(limits=rev(levels(res$Name))) + xlab('Correlation') + ylab('') + xlim(c(-0.25,0.25)) + ggtitle('HiC score vs GpC Accessibility')
p <- p + geom_vline(xintercept = 0,linetype=2) 
print(p)
dev.off()

###To fix - incorporate into figure S2

Figure3B <- function(res1,res2,zlim=c(-1,1),extra_cols,intervals1,intervals2,plot_res,ylim1,ylim2,fig_f,height,width){
  pdf(fig_f,height=height,width=width)
  
  layout(mat=matrix(c(1:12),byrow=F,ncol=4),widths = c(4,0.5,0.5,0.5),heights=c(0.5,0.5,4),respect = T)
  par(mar=c(1.5,0.5,1.5,0.5),mgp=c(1,0.5,0)) 
  plotMext(INPUTS=res1[1]$data, xlim = NULL,
           ylim = ylim1, main = '', xlab = "", ylab = '',
           plotScale = "linear", type = "full", error.estimates = T,
           legend = F, legend_ext = F,bty='n',
           xaxt='n',yaxt='n',yaxs='i',xaxs='i',
           colvec = extra_cols[1],par_mar=c(0.2,0.5,0.2,0.5),
           ln.v = FALSE, ln.h = NULL)
  axis(side = 4,labels=ylim1,at=ylim1,pos=4e4,tick = T,las=2,padj = c(0,1))
  plotMext(INPUTS=res1[2]$data, xlim = NULL,
           ylim = ylim2, main = '', xlab = "", ylab = '',
           plotScale = "linear", type = "full", error.estimates = T,
           legend = F, legend_ext = F,bty='n',
           xaxt='n',yaxt='n',yaxs='i',xaxs='i',
           colvec = extra_cols[2],par_mar=c(0.2,0.5,0.2,0.5),
           ln.v = FALSE, ln.h = NULL)
  axis(side = 4,labels=ylim2,at=ylim2,pos=4e4,tick = T,las=2,padj = c(0,1))
  params <- plot_aggregateHiC(cells=cells[1],pool=T,range_f=40000,filter_f=0,res_f=1000,plot_res=plot_res,grid_mode='loops',zlim=zlim,which_plot=c(2),
                              intervals1=intervals1,intervals2=intervals2,
                              interval1_name = '',interval2_name = '',
                              add_plot=T,plot_mean=T)
  par(mar=c(0,0,0,0))
  plot(0,type='n',axes=FALSE,ann=FALSE)
  plot(0,type='n',axes=FALSE,ann=FALSE)
  plotMext(INPUTS=res2[2]$data, xlim = NULL,
           ylim = ylim2, main = '', xlab = "", ylab = '',
           plotScale = "linear", type = "full", error.estimates = T,
           legend = F, legend_ext = F,bty='n',
           xaxt='n',yaxt='n',yaxs='i',xaxs='i',
           colvec = extra_cols[2],par_mar=c(1.5,0.2,1.5,0.2),
           ln.v = FALSE, ln.h = NULL,invert=T)
  #axis(side = 3,labels=ylim2,at=ylim2,pos=4e4,tick = T,las=1)
  plot(0,type='n',axes=FALSE,ann=FALSE)
  plot(0,type='n',axes=FALSE,ann=FALSE)
  plotMext(INPUTS=res2[1]$data, xlim = NULL,
           ylim = ylim1, main = '', xlab = "", ylab = '',
           plotScale = "linear", type = "full", error.estimates = T,
           legend = F, legend_ext = F,bty='n',
           xaxt='n',yaxt='n',yaxs='i',xaxs='i',
           colvec = extra_cols[1],par_mar=c(1.5,0.2,1.5,0.2),
           ln.v = FALSE, ln.h = NULL,invert=T)
  #axis(side = 3,labels=ylim1,at=ylim1,pos=4e4,tick = T,las=1)
  plot(0,type='n',axes=FALSE,ann=FALSE)
  plot(0,type='n',axes=FALSE,ann=FALSE)
  par(mar=c(1.5,0.5,1.5,2))
  image.scale(as.matrix(1),zlim=params$zlim, col=params$cols,axis.pos=4,adj=0.5,cex.axis = 1.2)
  dev.off()
}

###Plot Figures

tracks <- list.files(paste0(main_dir,'data/mm10/RAM_bigwigs'),pattern = 'merged',full.names = T)
res1 <- getPlotSetArray(tracks=tracks,features=paste0(main_dir,'data/mm10/microc/res1000/ES_loops_L.bed'),refgenome='mm10',type = 'mf',add_heatmap=F,xmin=40000,xmax=40000,bin = 50,ignore_strand = T)  #Extract CpG and CpG from loop coordinates
res2 <- getPlotSetArray(tracks=tracks,features=paste0(main_dir,'data/mm10/microc/res1000/ES_loops_R.bed'),refgenome='mm10',type = 'mf',add_heatmap=F,xmin=40000,xmax=40000,bin = 100,ignore_strand = T)

Figure3B(res1=res1,res2=res2,extra_cols=c(cpg_color,gpc_color),
         intervals1='ES_loops_L.bed',intervals2='ES_loops_R.bed',
         plot_res=1000,
         ylim1=c(45,80),ylim2=c(10,20),zlim=c(-0.75,0.75),fig_f='figures/Figure3B.pdf',height=6,width=7)

res1 <- getPlotSetArray(tracks=tracks,features=paste0(main_dir,'data/mm10/microc/res1000/ES_loops_L_noCTCF_either.bed'),refgenome='mm10',type = 'mf',add_heatmap=F,xmin=40000,xmax=40000,bin = 100,ignore_strand = T)  #Extract CpG and CpG from loop coordinates
res2 <- getPlotSetArray(tracks=tracks,features=paste0(main_dir,'data/mm10/microc/res1000/ES_loops_R_noCTCF_either.bed'),refgenome='mm10',type = 'mf',add_heatmap=F,xmin=40000,xmax=40000,bin = 100,ignore_strand = T)
Figure3B(res1=res1,res2=res2,extra_cols=c(cpg_color,gpc_color),
         intervals1='ES_loops_L_noCTCF_either.bed',intervals2='ES_loops_R_noCTCF_either.bed',
         plot_res=2000,
         ylim1=c(35,80),ylim2=c(10,25),zlim=c(-0.75,0.75),fig_f='figures/Figure3B_1.pdf',height=6,width=7)

#Fig3C


res1 <- getPlotSetArray(tracks=tracks,features=paste0(main_dir,'data/mm10/microc/res1000/ES_loops_L.bed'),refgenome='mm10',type = 'mf',add_heatmap=T,xmin=40000,xmax=40000,bin = 2000,ignore_strand = T)  #Extract CpG and CpG from loop coordinates
res2 <- getPlotSetArray(tracks=tracks,features=paste0(main_dir,'data/mm10/microc/res1000/ES_loops_R.bed'),refgenome='mm10',type = 'mf',add_heatmap=T,xmin=40000,xmax=40000,bin = 2000,ignore_strand = T)  #Extract CpG and CpG from loop coordinates
pdf('figures/Figure3C.pdf',width=6,height=8)
par(mfrow=c(1,2))
plotHeatmap(res1[2],labels = c(''),sortrows = FALSE,clspace = rev(colorpalette('reds',12)),
            clstmethod=FALSE,raster = T,embed=T)
plotHeatmap(res2[2],labels = c(''),sortrows = FALSE,clspace = rev(colorpalette('reds',12)),
            clstmethod=FALSE,raster = T,embed=T)
dev.off()

###Fig3D


peaks <- c(paste0(main_dir,'data/mm10/beds/gNOMEpeaks/sig_NOMEPeaks_bedform.bed'),paste0(main_dir,'data/mm10/ATAC_seq/ATAC_peaks.bed'),paste0(main_dir,'data/mm10/DNase_ENCSR000CMW/DNA_E14_peaks_encode.bed'),paste0(main_dir,'data/mm10/ATAC_seq/ATAC_peaks_shuffle.bed'))
res <- getPlotSetArray(tracks="/home/hpc/bonev/projects/ram/data/mm10/mnase_GSM1400766/mES_MNAse_SRR1302810.bw",features=peaks,refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 5,ignore_strand = F)
pdf('figures/FigureS2J.pdf',width=5,height=5)
plotAverage(plotset=res, labels = c('gNOME','ATAC','DHS','shuffled'), xlim = NULL,
            ylim = NULL, main = '', xlab = "", ylab = 'MNase-seq (GSE58101)',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = T, legend_ext = F, legend_pos = "bottomright",
            legend_ext_pos = "bottomright", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12)
axis(side = 1,labels=c('-1kb','0','+1kb'),at=c(-1000,0,1000),pos=10,tick = T)
dev.off()



source('/home/faye/Misc_requests/hic/config.R')
script_dir <- '/home/hpc/bonev/projects/ram/results/mm10/hic/'
source(paste0(script_dir,'scripts/main_functions.R'))
source(paste0(script_dir,'scripts/aux_functions.R'))
source(paste0(script_dir,'scripts/plot_functions.R'))
path <- paste0(main_f,'analysis/compartments/')
binSize=2.5e5
#rank_matrix(eigen_tracks=gtrack.ls('eigen',250),cells=cells[1:3],chrs=c(1:19,'X'),path=path,binSize=binSize,min_dist=1e7,ranks=100)
file_f='/home/faye/Misc_requests/hic/data/extractedBins/test'
plot_rankMatrix(file_f=file_f,zlim=c(-1.2,1.2),cells=cells,col=wide_red_blue_pal(1000),plot_chip=FALSE,out_f = 'figures/FigureS2K.pdf',)





#Fig1B Merged plotted across CTCF
peaks <- paste0(main_dir,'data/mm10/beds/ES_CTCF_top30k_motifs.bed')
tracks <- c(paste0(main_dir,'data/mm10/RAM_bigwigs/3DRAM_ES_250k_merged_GpC_cov10x.bw'),paste0(main_dir,'data/mm10/WGBS/WGBS_mESC_merged.NOMe.GpC_cov10x.bw'))
res <- getPlotSetArray(tracks=tracks,features=peaks,refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 5,ignore_strand = F)
pdf('figures/Rev_Figure1B.pdf',width=5,height=5)
plotAverage(plotset=res, labels = c('3DRAM-seq GpC','WGBS GpC'), xlim = NULL,
            ylim = c(-5,60), main = '', xlab = "", ylab = '% Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = T, legend_ext = F, legend_pos = "topright",
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = c(cpg_color,gpc_color))
axis(side = 1,labels=c('-1kb','CTCF->','+1kb'),at=c(-1000,0,1000),pos =-5,tick = T)
dev.off()


#Fig1C Merged plotted across Klf4
peaks <- paste0(main_dir,'data/mm10/beds/ES_Nrf1_motifs.bed')
res <- getPlotSetArray(tracks=tracks,features=peaks,refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 2,ignore_strand = F)
pdf('figures/Rev_Figure1C.pdf',width=5,height=5)
plotAverage(plotset=res, labels = c('3DRAM-seq GpC','WGBS GpC'), xlim = NULL,
            ylim = c(-5,50), main = '', xlab = "", ylab = '% Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = T, legend_ext = F, legend_pos = "topright",
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = c(cpg_color,gpc_color))
axis(side = 1,labels=c('-1kb','Nrf1','+1kb'),at=c(-1000,0,1000),pos =-5,tick = T)
dev.off()

res <- getPlotSetArray(tracks=tracks[2],features=peaks,refgenome='mm10',type = 'mf',add_heatmap=F,xmin=50,xmax=50,bin = 2,ignore_strand = F)
pdf('figures/Figure1C_i.pdf',width=3,height=3)
plotAverage(plotset=res, labels = , xlim = NULL,
            ylim = c(25,50), main = '', xlab = "", ylab = '',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = F, legend_ext = F, legend_pos = "topright",
            legend_ext_pos = "topleft", cex.axis = 10, cex.lab = 10,xaxt='n',
            cex.main = 10, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = gpc_color)
axis(side = 1,labels=c('-100bp','Nrf1','+100bp'),at=c(-100,0,100),pos =25,tick = T,cex.axis=0.8)
dev.off()


