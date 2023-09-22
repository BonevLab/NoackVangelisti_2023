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
library(Signac)
library(JASPAR2020)
library(Matrix)
library(TFBSTools)
library(motifmatchr)
library(Hmisc)
library(GenomicRanges)
library(ChIPpeakAnno)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GenomicRanges)
library(SummarizedExperiment)
library(seqplots)
require(plyr)
require(dplyr)
library(ChIPseeker)
library("methylKit", lib.loc="/home/hpc/bonev/software/R/library")
library(ggpointdensity)
library(gridExtra)
library(pheatmap)
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(motifmatchr)
library(Hmisc)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(ggpubr)
theme_border <- theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line=element_blank(),panel.border = element_rect(color = "black", size = 0.5),legend.margin = margin(0,0.5,0,0.5,unit='inch'))

tracks <- list.files(paste0(main_dir,'data/hg38/RAM_bigwigs/'),pattern = 'merged',full.names = T)

########################
###Functions###
########################


########################
###Call Functions###
########################
region <- toGRanges(paste0(main_dir,'data/hg38/beds/gNOME_peaks/mergedgNOME_sig_peaks.bed'))
files_f <- as.list(list.files(paste0(main_dir,'data/hg38/RAM_covfiles/replicates/'),pattern = 'GpC',full.names = T))
DMR_call(con = 'GpC',files = files_f,genome = 'hg38',threshold=10,difference=0,min_Cov=10,qval=0.05,
         high_cov_perc=99.9,statistical_method='over_disp_chi',region=region,output_path=paste0(main_dir,'data/hg38/DMR_analysis/'))
# To generate shuffled regions 
# for i in /home/hpc/bonev/projects/ram/data/hg38/DMR_analysis/*.bed; do bedtools shuffle -i $i -g /home/hpc/bonev/annotations/hg38/hg38.chrom.sizes -noOverlapping -excl /home/hpc/bonev/annotations/hg38/encode/hg38.blacklist.bed >${i/.bed/_shuffled.bed};done  

#Figure 5A -> DMR GpC and CpG levels 
p <- plot_DMR(con='GpC',files=files_f,genome='hg38',threshold=10,difference=0,min_Cov=10,qval=0.05,labels=c('RGC GpC Accessibility','IPC GpC Accessibility'),
         high_cov_perc=99.9,statistical_method='over_disp_chi',region=region,output_path=paste0(main_dir,'data/hg38/DMR_analysis/'),colours=cell_colors) #adjust context rest as above ! 
pdf('figures/Figure5A.pdf',width = 6,height=6,useDingbats = F)
print(p + xlab("RGC GpC Accessibility (%)") + ylab("IPC GpC Accessibility (%)") + theme(legend.position = "none"))
dev.off()

#Fig5B <- average methylation levels across DAR sites 
peaks <-paste0(main_dir,'data/hg38/DMR_analysis/GpC_hypo_0.05_0_over_disp_chi.bed')
res <- getPlotSetArray(tracks=tracks[grep('CpG',tracks)],features=peaks,refgenome='hg38',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 20,ignore_strand = F)
pdf('figures/Figure5B_1.pdf',width=4.5,height=4.55)
plotAverage(plotset=res, labels = c('RGC','IPC'), xlim = NULL,
            ylim = c(15,95), main = '', xlab = "", ylab = "",
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = T, legend_ext = F, legend_pos = "bottomright",
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = cell_colors)
title(ylab = "CpG Methylation (%)", line = 2.5)
axis(side = 1,labels=c('-1kb','RGC DAR','+1kb'),at=c(-1000,0,1000),pos =15,tick = T,padj=-0.75)
dev.off()

peaks <- paste0(main_dir,'data/hg38/DMR_analysis/GpC_hyper_0.05_0_over_disp_chi.bed')
res <- getPlotSetArray(tracks=tracks[grep('CpG',tracks)],features=peaks,refgenome='hg38',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 20,ignore_strand = F)
pdf('figures/Figure5B_2.pdf',width=4.5,height=4.5)
plotAverage(plotset=res, labels = c('RGC','IPC'), xlim = NULL,
            ylim = c(15,95), main = '', xlab = "", ylab = '',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = T, legend_ext = F, legend_pos = "bottomright",
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = cell_colors)
title(ylab = "CpG Methylation (%)", line = 2.5)
axis(side = 1,labels=c('-1kb','IPC DAR','+1kb'),at=c(-1000,0,1000),pos =15,tick = T,padj=-0.75)
dev.off()


#Figure 5C -> GOterm analysis


nsc_mat <- extract_misha_ep_pairs(enhancer_bed = paste0(main_dir,'data/hg38/DMR_analysis/NSC_distal_DAR.bed'),
                                  tss_bed = paste0(main_dir,'data/hg38/beds/hg38_proteinCoding_TSS.bed'),
                                  interval_window = 200,min_dist = 5e3,max_dist = 2e6,expand = c(-1e4,1e4),
                                  domains_f = "hic.3DRAM_D45_Pax6.ins_250_domains_expanded",
                                  score_tracks = score_tracks,hic_names = c('NSCscore','IPCscore'),
                                  other_tracks =gtrack.ls('methylation','GpC','10x'),other_names = c('NSC GpC accessbility','IPC GpC Accessibility'))

nsc_p1 <- GOterm_enrichment(distal_mat=nsc_mat,promoter_bed=paste0(main_dir,'data/hg38/DMR_analysis/NSC_promoter_DAR.bed'),
                            all_TSS=paste0(main_dir,'data/hg38/beds/hg38_proteinCoding_TSS.bed'),
                            closest_TSS = T)
nsc_p2 <- GOterm_enrichment(distal_mat=nsc_mat,promoter_bed=paste0(main_dir,'data/hg38/DMR_analysis/NSC_promoter_DAR.bed'),
                        all_TSS=paste0(main_dir,'data/hg38/beds/hg38_proteinCoding_TSS.bed'),which_HiC = 'NSCscore',
                        closest_TSS = F)
pdf(paste0('figures/Figure5C_1.pdf'),width=9,height=4,useDingbats = F)
print(dotplot(nsc_p2$p, showCategory=10,orderBy = "GeneRatio", x='GeneRatio')+theme_cowplot())
dev.off()




ipc_mat <- extract_misha_ep_pairs(enhancer_bed = paste0(main_dir,'data/hg38/DMR_analysis/IPC_distal_DAR.bed'),
                                  tss_bed = paste0(main_dir,'data/hg38/beds/hg38_proteinCoding_TSS.bed'),
                                  interval_window = 200,min_dist = 5e3,max_dist = 2e6,expand = c(-1e4,1e4),
                                  domains_f = "hic.3DRAM_D45_Tbr2.ins_250_domains_expanded",
                                  score_tracks = score_tracks,hic_names = c('NSCscore','IPCscore'),
                                  other_tracks =gtrack.ls('methylation','GpC','10x'),other_names = c('NSC GpC accessbility','IPC GpC Accessibility'))

ipc_p1 <- GOterm_enrichment(distal_mat=ipc_mat,promoter_bed=paste0(main_dir,'data/hg38/DMR_analysis/IPC_promoter_DAR.bed'),
                        all_TSS=paste0(main_dir,'data/hg38/beds/hg38_proteinCoding_TSS.bed'),
                        closest_TSS = T)
ipc_p2 <- GOterm_enrichment(distal_mat=ipc_mat,promoter_bed=paste0(main_dir,'data/hg38/DMR_analysis/IPC_promoter_DAR.bed'),
                        all_TSS=paste0(main_dir,'data/hg38/beds/hg38_proteinCoding_TSS.bed'),which_HiC = 'IPCscore',
                        closest_TSS = F)
pdf(paste0('figures/Figure5C_2.pdf'),width=9,height=4,useDingbats = F)
print(dotplot(ipc_p2$p, showCategory=10,orderBy = "GeneRatio", x='GeneRatio'))
dev.off()

# Comparison to gene expression 

res = read.table('/home/hpc/bonev/projects/rna/ram/IPCvsNSC_res_proteinCoding_all.tsv')
nsc_mat_res <- nsc_p2$mat
nsc_mat_res$log2FoldChange <- res$log2FoldChange[match(nsc_mat_res$geneName,row.names(res))]
ipc_mat_res <- ipc_p2$mat
ipc_mat_res$log2FoldChange <- res$log2FoldChange[match(ipc_mat_res$geneName,row.names(res))]

nsc_mat_res <- nsc_mat_res[!duplicated(nsc_mat_res$geneName),'log2FoldChange',drop=F]
ipc_mat_res <- ipc_mat_res[!duplicated(ipc_mat_res$geneName),'log2FoldChange',drop=F]
fc_res <- data.frame(Condition=factor(c(rep('RGC',length(nsc_mat_res$log2FoldChange)),rep('IPC',length(ipc_mat_res$log2FoldChange))),levels=c('RGC','IPC')),Value=c(nsc_mat_res$log2FoldChange,ipc_mat_res$log2FoldChange))
p1 <- ggplot(fc_res,aes(x=Condition,y=Value,fill=Condition)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8,outlier.shape = NA) + xlab('')
p1 <- p1 + scale_fill_manual(values=cell_colors) + ylab('log2 Fold Change') + theme(legend.position = "none") + stat_summary(fun=mean, geom="point", shape=8, size=3, color="darkred", fill="darkred") + coord_cartesian(ylim=c(-3.5,4.5))
p1 <- p1 + stat_compare_means(comparisons = list(c('RGC','IPC')),label = "p.format",method='wilcox',label.y = c(4.5,4.5),tip.length = c(0.01,0.03))
pdf('figures/Figure5CD.pdf',height=4.5,width=4)
p1 + geom_hline(yintercept = 0,linetype = 2,col='black') + ggtitle('Expression of linked genes')
dev.off()
#Figure 5D-> motif enrichment analysis######
background_peaks <- paste0(main_dir,'data/hg38/beds/gNOME_peaks/mergedgNOME_sig_peaks.bed')
pwm <- readRDS(paste0(main_dir,'/data/hg38/pwms_FPKM1.RDS')) 

peaks <- paste0(main_dir,'data/hg38/DMR_analysis/GpC_hyper_0.05_0_over_disp_chi.bed')
IPC_motifs <- enrichedMotifs(peaks=peaks,bg_peaks = background_peaks,genome = 'hg38',genome_BS = BSgenome.Hsapiens.UCSC.hg38,
                             features= c('Neurog2(var.2)','Neurod1','Neurod2','Eomes','Ctcf','Nrf1','Fos::jun','Lhx2','Sox2','Nfia','Rfx4','Zbtb18','Tead2','Irf1','Pou3f3','Rest','Hes1','Nr2f1'),
                         cols=c("blue",'grey80','red'),logFC=0.20,logP=2,point.size=4,anno.size=6)
pdf(paste0('figures/Figure5D_1.pdf'),width=5,height=5,useDingbats = F)
print(IPC_motifs$p)
dev.off()


peaks <- paste0(main_dir,'data/hg38/DMR_analysis/GpC_hypo_0.05_0_over_disp_chi.bed')
NSC_motifs <- enrichedMotifs(peaks=peaks,bg_peaks = background_peaks,genome_f = 'hg38',genome_BS = BSgenome.Hsapiens.UCSC.hg38,
                             features=c('Neurog2(var.2)','Neurod1','Neurod2','Eomes','Ctcf','Nrf2','Fos::jun','Lhx2','Sox2','Nfia','Rfx4','Zbtb18','Tead2','Irf1','Pou3f3','Rest','Hes1','Nr2f1'),
                             cols=c("blue",'grey80','red'),logFC=0.2,logP=2,point.size=4,anno.size=6)
pdf(paste0('figures/Figure5D_2.pdf'),width=5,height=5,useDingbats = F)
print(NSC_motifs$p)
dev.off()

#Figure 5E-> AggregatedHiC of Ngn2  with/without motif 
#source('/home/fnoack/projects/3DRAM/scripts/figures/Flo_aux_functions.R')
#chip_f <- "/home/hpc/bonev/projects/ram/data/hg38/DMR_analysis/GpC_hyper_0.05_0_over_disp_chi.bed" #Full path to peaks
#pwm_f <- readRDS('/home/hpc/bonev/projects/ram/data/mm10/combined_pwm.RDS') ### motif matrix
#motif.names <- as.data.frame(capitalize(tolower(as.vector(name(pwm_f))))) ### Explore this to see how exacly to spell the motif name you need for the function
#motif.names$ID<-capitalize(tolower(as.vector(ID(pwm_f))))
#motif_name<-'Neurog2(var.2)' 
#genome='hg38'
#cutoff=0.0005
#bed_dir='/home/hpc/bonev/projects/ram/data/hg38/DMR_analysis/Peaks_TF_Centered/GpC_hyper_' #dir and name where to safe the motif centered beds 
#DMR_with_motifs(chip_f,pwm_f,motif.names,motif_name,genome,cutoff,bed_dir) #generates bed files 
#for (cell in cells){
#  submit_aggregateHiC(cells=cell,tracks=all_tracks,range_f=40000,res_f=1000,filter_f=0,
#                      intervals1='/home/hpc/bonev/projects/ram/data/hg38/DMR_analysis/Peaks_TF_Centered/GpC_hyper_Neurog2var.2_0.0005_Neurog2var.2.bed', 
#                      intervals2='/home/hpc/bonev/projects/ram/data/hg38/DMR_analysis/Peaks_TF_Centered/GpC_hyper_Neurog2var.2_0.0005_Neurog2var.2.bed',
#                      grid_mode='1D',mem=80) 
#}
# for (cell in cells){
#   pdf(paste0('/home/fnoack/projects/3DRAM/Figures/Figure5_AggregatedHiC_GpC_hyper_Ngn2_genomePeaks_',cell,'.pdf'),width=4,height=4)
#   layout(matrix(c(1:2),nrow=1,byrow=F),widths = c(4,2),heights=c(4),respect = T)
#   params <- plot_aggregateHiC(cells=cell,pool=T,range_f=40000,filter_f=0,res_f=1000,plot_res=4000,grid_mode='1D',zlim=c(-0.4,0.4),which_plot=c(2),
#                               intervals1='GpC_hyper_Neurog2var.2_0.0005_Ngn2.bed',intervals2='GpC_hyper_Neurog2var.2_0.0005_Ngn2.bed',
#                               interval1_name = 'NEUROG2',interval2_name = 'NEUROG2',
#                               add_plot=T,plot_mean=T)   #plot_res=1000 does not work because of ->  Error in (floor(ncol(o)/2)):(floor(ncol(o)/2) + 2) : argument of length 0 -> do with 2000 and can be after the pool file is generated with 1000 
#   par(mar=c(1.5,2,1.5,3),mgp=c(1,0.5,0)) 
#   image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
#   dev.off()
# }

pdf('figures/Figure5E_1.pdf',width=4,height=8.5)
layout(matrix(c(1:3,3),nrow=2,byrow=F),widths = c(4,1.5),heights=c(4,4),respect = T)
for (cell in cells){
  params <- plot_aggregateHiC(cells=cell,pool=T,range_f=40000,filter_f=0,res_f=1000,plot_res=4000,grid_mode='1D',zlim=c(-0.4,0.4),which_plot=c(2),
                              intervals1='GpC_hyper_Neurog2var.2_0.0005_Ngn2.bed',intervals2='GpC_hyper_Neurog2var.2_0.0005_Ngn2.bed',
                              interval1_name = '',interval2_name = '',
                              add_plot=T,plot_mean=T) 
}
par(mar=c(1.5,1.5,1.5,3),mgp=c(1,0.5,0)) 
image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
dev.off()

pdf('figures/Figure5E_2.pdf',width=4,height=8.5)
layout(matrix(c(1:3,3),nrow=2,byrow=F),widths = c(4,1.5),heights=c(4,4),respect = T)
for (cell in cells){
  params <- plot_aggregateHiC(cells=cell,pool=T,range_f=40000,filter_f=0,res_f=1000,plot_res=4000,grid_mode='1D',zlim=c(-0.4,0.4),which_plot=c(2),
                              intervals1='GpC_hypo_Lhx2_0.0005_Lhx2.bed',intervals2='GpC_hypo_Lhx2_0.0005_Lhx2.bed',
                              interval1_name = '',interval2_name = '',
                              add_plot=T,plot_mean=T) 
}
par(mar=c(1.5,1.5,1.5,3),mgp=c(1,0.5,0)) 
image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
dev.off()

#Figure 5X plot example regions
misha_meth_tracks <- c("methylation.3DRAM_Pax6_CpG_merged_10x","methylation.3DRAM_Tbr2_CpG_merged_10x","methylation.3DRAM_Pax6_GpC_merged_10x","methylation.3DRAM_Tbr2_GpC_merged_10x") 
tss <- gintervals.load(tss_f)

anno <- read.table(paste0(main_dir,'data/hg38/DMR_analysis/Peaks_TF_Centered/GpC_hypo_Lhx2_0.0005.bed'))
anno <- gintervals(anno[,1],anno[,2],anno[,3])
anno <- intervals.normalize(anno,1000)   #Makes all intervals 1kb long

plotMisha(main_f=main_f,targetGene='GAS1',outDir='figures/',out_f='Figure5F_1',upstream=6e5,downstream=5e4,chipYlim=matrix(c(0,100,0,100,0,4,0,4),nrow = 4,ncol = 2,byrow = T),
          chipNames=c('','','',''),window_scale=1.8,chipRes=20,pointCEX=1,conditions=score_tracks,binSize=5e3,radius=2e4,
          chipTracksToExtract=c(misha_meth_tracks[grep('GpC',misha_meth_tracks)],'rnaseq.Pax6_HMGU1_Day45_RNA_merged','rnaseq.Tbr2_HMGU1_Day45_RNA_merged'),chipColors=rep(cell_colors,2), 
          annIntervals=anno,
          plotOrder=list(scores=TRUE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,anno=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.2, VP=1.5, loops=2.2, rna=0.6, chip=0.6,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.2))

anno <- read.table(paste0(main_dir,'data/hg38/DMR_analysis/Peaks_TF_Centered/GpC_hyper_Neurog2var.2_0.0005_Neurog2var.2.bed'))
anno <- gintervals(anno[,1],anno[,2],anno[,3])
anno <- intervals.normalize(anno,1000)   #Makes all intervals 1kb long

plotMisha(main_f=main_f,targetGene='NFIA',outDir='figures/',out_f='Figure5F_2',upstream=1.3e6,downstream=2e5,chipYlim=matrix(c(0,100,0,100,0,10.8,0,10.8),nrow = 4,ncol = 2,byrow = T),
          chipNames=c('','','',''),window_scale=1.8,chipRes=20,pointCEX=1,conditions=score_tracks,binSize=5e3,radius=2e4,
          chipTracksToExtract=c(misha_meth_tracks[grep('GpC',misha_meth_tracks)],'rnaseq.Pax6_HMGU1_Day45_RNA_merged','rnaseq.Tbr2_HMGU1_Day45_RNA_merged'),chipColors=rep(cell_colors,2), 
          annIntervals=anno,
          plotOrder=list(scores=TRUE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,anno=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.2, VP=1.5, loops=2.2, rna=0.6, chip=0.6,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.2))


smf_files <- list.files(paste0(main_dir,'results/hg38/pair_smf/'),pattern = 'fst',full.names = T)
anchor_files <- list.files(paste0(main_dir,'data/hg38/beds/gNOME_peaks/motif_centered/'),pattern = 'bed',full.names = T)
bed_files <- list.files(paste0(main_dir,'data/mm10/beds/'),pattern = 'bed',full.names = T)

pair_smf_f <- read.fst(smf_files[grep('Neurog2_GpCPax6',smf_files)])
binSize=100   #window +- interval center
res <- filter_smf(df=pair_smf_f,binSize=binSize,filter_dist=c(1000,4e6),ncalls=1,out_res_coord = TRUE) 

plot_res <- suppressWarnings(plot_paired_smf(res,TF_name='NEUROG2',met_lims=c(0.2,0.5),
                                             cluster_names=c('C1','C2','C3','C4'),
                                             cols = colorRamp2(c(0,100), c("grey","red")),
                                             mat_name='%GpC',
                                             ggplot_cols=hue_pal()(4),
                                             window=1,
                                             window2=10,
                                             binSize=binSize,
                                             pair_smf_f=pair_smf_f,
                                             meth_colors=c('blue','red'),
                                             point.size=16,
                                             plot_binSize=200))

pdf('figures/Figure5J_1.pdf',height=8,width=4)
draw(plot_res$hm,heatmap_legend_side = "bottom",annotation_legend_side = "bottom")
dev.off()

pdf('figures/Figure5K_1.pdf',height=3,width=2.5,pointsize = 30,useDingbats = F)
plot_res$p2$p_la|plot_res$p2$p_ra
dev.off()


pair_smf_f <- read.fst(smf_files[grep('Neurog2_GpCTbr2',smf_files)])
binSize=100   #window +- interval center
res <- filter_smf(df=pair_smf_f,binSize=binSize,filter_dist=c(1000,4e6),ncalls=1,out_res_coord = TRUE) 

plot_res <- suppressWarnings(plot_paired_smf(res,TF_name='NEUROG2',met_lims=c(0.2,0.5),
                                             cluster_names=c('C1','C2','C3','C4'),
                                             cols = colorRamp2(c(0,100), c("grey","red")),
                                             mat_name='%GpC',
                                             ggplot_cols=hue_pal()(4),
                                             window=1,
                                             window2=10,
                                             binSize=binSize,
                                             pair_smf_f=pair_smf_f,
                                             meth_colors=c('blue','red'),
                                             point.size=16,
                                             plot_binSize=200))

pdf('figures/Figure5J_2.pdf',height=8,width=4)
draw(plot_res$hm,heatmap_legend_side = "bottom",annotation_legend_side = "bottom")
dev.off()
pdf('figures/Figure5K_2.pdf',height=3,width=2.5,pointsize = 30,useDingbats = F)
plot_res$p2$p_la|plot_res$p2$p_ra
dev.off()


pair_smf_f <- read.fst(smf_files[grep('Pax6_CpG_GpC_Neurog2v2_500R1.',smf_files)])

res <- filter_smf(df=pair_smf_f, binSize=100, filter_dist=c(0,300), ncalls=1, 
                  out_res_coord = TRUE) 
plot_res <- suppressWarnings(plot_paired_smf(res,TF_name='NEUROG2',met_lims=c(0.2,0.5),
                                             cluster_names=c('C1','C2','C3','C4'),
                                             cols = colorRamp2(c(0,100), c("grey","red")),
                                             mat_name='% methylation',
                                             ggplot_cols=hue_pal()(4),
                                             window=1,
                                             window2=10,
                                             binSize=binSize,
                                             pair_smf_f=pair_smf_f,
                                             meth_colors=c('blue','red'),
                                             point.size=16,
                                             plot_binSize=120))
pdf('figures/Figure5L_1.pdf',height=8,width=4)
draw(plot_res$hm,heatmap_legend_side = "bottom",annotation_legend_side = "bottom")
dev.off()
pdf('figures/Figure5M_1.pdf',height=3,width=2.5,pointsize = 30,useDingbats = F)
plot_res$p2$p_la|plot_res$p2$p_ra
dev.off()

pair_smf_f <- read.fst(smf_files[grep('Tbr2_CpG_GpC_Neurog2v2_500R1.',smf_files)])

res <- filter_smf(df=pair_smf_f, binSize=100, filter_dist=c(0,300), ncalls=1, 
                  out_res_coord = TRUE) 
plot_res <- suppressWarnings(plot_paired_smf(res,TF_name='NEUROG2',met_lims=c(0.2,0.5),
                                             cluster_names=c('C1','C2','C3','C4'),
                                             cols = colorRamp2(c(0,100), c("grey","red")),
                                             mat_name='% methylation',
                                             ggplot_cols=hue_pal()(4),
                                             window=1,
                                             window2=10,
                                             binSize=binSize,
                                             pair_smf_f=pair_smf_f,
                                             meth_colors=c('blue','red'),
                                             point.size=16,
                                             plot_binSize=120))
pdf('figures/Figure5L_2.pdf',height=8,width=4)
draw(plot_res$hm,heatmap_legend_side = "bottom",annotation_legend_side = "bottom")
dev.off()
pdf('figures/Figure5M_2.pdf',height=3,width=2.5,pointsize = 30,useDingbats = F)
plot_res$p2$p_la|plot_res$p2$p_ra
dev.off()

misha_meth_tracks <- c("methylation.3DRAM_Pax6_CpG_merged_10x","methylation.3DRAM_Tbr2_CpG_merged_10x","methylation.3DRAM_Pax6_GpC_merged_10x","methylation.3DRAM_Tbr2_GpC_merged_10x") 
tss <- gintervals.load(tss_f)

source(paste0(main_dir,'results/hg38/hic/config.R'))
anno <- read.table(paste0(main_dir,'data/hg38/DMR_analysis/Peaks_TF_Centered/GpC_hypo_Lhx2_0.0005.bed'))
anno <- gintervals(anno[,1],anno[,2],anno[,3])
anno <- intervals.normalize(anno,1000)   #Makes all intervals 1kb long

plotMisha(main_f=main_f,targetGene='FAM107A',outDir='figures/',out_f='FAM107A_human',upstream=2e5,downstream=5e4,chipYlim=matrix(c(0,100,0,100,0,3.2,0,3.2),nrow = 4,ncol = 2,byrow = T),
          chipNames=c('','','',''),window_scale=1.8,chipRes=20,pointCEX=2,conditions=score_tracks,binSize=5e3,radius=2e4,
          chipTracksToExtract=c(misha_meth_tracks[grep('GpC',misha_meth_tracks)],'rnaseq.Pax6_HMGU1_Day45_RNA_merged','rnaseq.Tbr2_HMGU1_Day45_RNA_merged'),chipColors=rep(cell_colors,2), 
          annIntervals=anno,
          plotOrder=list(scores=TRUE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,anno=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.2, VP=1.5, loops=2.2, rna=0.6, chip=0.6,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.2))


source(paste0(main_dir,'results/mm10/hic/config.R'))
plotMisha(main_f=main_f,targetGene='Fam107a',outDir='figures/',out_f='FAM107A_mouse',upstream=2e5,downstream=5e4,
          chipNames=c('','',''),window_scale=1.8,chipRes=20,pointCEX=2,conditions=score_tracks[1:2],binSize=5e3,radius=2e4,
          chipTracksToExtract=c('scATAC.E14_NSC','scATAC.E14_IPC'),chipColors=rep('black',2), 
          plotOrder=list(scores=TRUE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,anno=FALSE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.2, VP=1.5, loops=2.2, rna=0.6, chip=0.6,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.2))


source(paste0(main_dir,'results/hg38/hic/config.R'))
anno <- gintervals('chr3',58581032,58581042)
anno <- intervals.normalize(anno,10)   #Makes all intervals 1kb long

plotMisha(main_f=main_f,targetGene='chr3,58580782,58581442',outDir='figures/',out_f='FAM107A_enh_human',upstream=2e3,downstream=2e3,
          chipNames=c('NSC ATAC','IPC ATAC'),window_scale=1.8,chipRes=10,pointCEX=2,conditions=score_tracks,binSize=5e3,radius=2e4,
          chipTracksToExtract=c("atac.HMGU1_D45_Pax6_ATAC_merged","atac.HMGU1_D45_Tbr2_ATAC_merged"),chipColors=rep(cell_colors,2),
          methTracksToExtract=c(misha_meth_tracks),methColors=rep(cell_colors,2),methNames=c('NSC CpG','IPC CpG','NSC GpC','IPC GpC'),
          annIntervals=anno,
          plotOrder=list(scores=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,anno=TRUE,meth=TRUE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.2, VP=1.5, loops=2.2, rna=0.6, chip=0.6,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.2))


source(paste0(main_dir,'results/hg38/hic/config.R'))
anno <- read.table(paste0(main_dir,'data/hg38/DMR_analysis/Peaks_TF_Centered/GpC_hypo_Lhx2_0.0005.bed'))
anno <- gintervals(anno[,1],anno[,2],anno[,3])
anno <- intervals.normalize(anno,1000)   #Makes all intervals 1kb long

plotMisha(main_f=main_f,targetGene='PCDH9',outDir='figures/',out_f='PCDH9_test',upstream=6e5,downstream=1e5,
          chipNames=c('','','',''),window_scale=1.8,chipRes=20,pointCEX=1,conditions=score_tracks,binSize=5e3,radius=2e4,
          chipTracksToExtract=c(misha_meth_tracks[grep('GpC',misha_meth_tracks)],'rnaseq.Pax6_HMGU1_Day45_RNA_merged','rnaseq.Tbr2_HMGU1_Day45_RNA_merged'),chipColors=rep(cell_colors,2), 
          #chipYlim=matrix(c(0,100,0,100,0,1.2,0,1.2),ncol = 2,byrow = T),
          plotOrder=list(scores=TRUE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,anno=FALSE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.2, VP=1.5, loops=2.2, rna=0.6, chip=0.6,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.2))

source(paste0(main_dir,'results/mm10/hic/config.R'))
plotMisha(main_f=main_f,targetGene='Fbxo32',outDir='figures/',out_f='Fbxo32_mouse',upstream=1e5,downstream=2.3e5,
          chipNames=c('','',''),window_scale=1.8,chipRes=20,pointCEX=2,conditions=score_tracks[1:2],binSize=5e3,radius=2e4,
          chipTracksToExtract=c('scATAC.E14_NSC','scATAC.E14_IPC'),chipColors=rep('black',2), 
          plotOrder=list(scores=TRUE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,anno=FALSE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.2, VP=1.5, loops=2.2, rna=0.6, chip=0.6,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.2))







### Test CpG/GpC per motif
library(doParallel)
registerDoParallel(cores=8)

misha_meth_tracks <- c("methylation.3DRAM_Pax6_CpG_merged_10x","methylation.3DRAM_Tbr2_CpG_merged_10x","methylation.3DRAM_Pax6_GpC_merged_10x","methylation.3DRAM_Tbr2_GpC_merged_10x") 
bed_fs <- list.files('/home/hpc/bonev/projects/ram/results/hg38/beds/motif_ALLgNOME_all/',pattern = 'bed',full.names = T)

res_list <- foreach (bed_f=bed_fs) %do% {
  motif_name <- gsub("\\.bed.*","",basename(bed_f))
  bed <- read.table(bed_f)
  colnames(bed) <- c('chrom','start','end')
  df <- misha_extract(misha_meth_tracks,bed,iterator = bed,mode = 'avg')
  colnames(df)[4:8] <- c('NSC_CpG','IPC_CpG','NSC_GpC','IPC_GpC','motif')
  df[,8] <- motif_name
  return(df)
}
names(res_list) <- gsub("\\.bed.*","",basename(bed_fs))
saveRDS(res_list,'/home/hpc/bonev/projects/ram/results/hg38/beds/motif_ALLgNOME_all/all_CpG_GpC.RDS')

# res <- sapply(res_list,function(x){
#   return(colMeans(x[,4:7],na.rm=T))
# })
# 
# res <- sapply(res_list,function(x){
#   return(as.numeric(cor.test(x$IPC_GpC-x$NSC_GpC,x$IPC_CpG-x$NSC_CpG)$estimate))
# })
mat <- res_list[['NEUROG2-var.2-']]
mat$dCpG <- mat$IPC_CpG-mat$NSC_CpG
mat$dGpC <- mat$IPC_GpC-mat$NSC_GpC

p <- ggplot(mat, aes( x = dGpC, y = dCpG)) + geom_pointdensity(shape=19,size=0.6,alpha=1) + scale_color_gradientn(colours = c('black','red','orange'),name='')
p <- p + geom_vline(xintercept = 0,linetype = 2,col='darkblue')+ geom_hline(yintercept = 0,linetype = 2,col='darkblue') + xlim(-70,70) + ylim(-100,100) + theme(legend.position = c(0.7,0.85),legend.direction = 'horizontal')
p <- p + xlab('Accessibility change (IPC-NSC)') + ylab('Methylation change (IPC-NSC)')
pdf('figures/FigureS5J.pdf',height=5,width=5)
p + geom_text(x=-50, y=85,size=6, label=paste0('r=',round(as.numeric(cor.test(x=mat$dGpC,y=mat$dCpG,method = 'pearson',use='complete')$estimate),2)))
dev.off()

### Calculate change in Hi-C as function of expression

nsc_mat <- extract_misha_ep_pairs(enhancer_bed = paste0(main_dir,'data/hg38/DMR_analysis/NSC_distal_DAR.bed'),
                                  tss_bed = paste0(main_dir,'data/hg38/beds/hg38_proteinCoding_TSS.bed'),
                                  interval_window = 200,min_dist = 5e3,max_dist = 2e6,expand = c(-1e4,1e4),
                                  domains_f = "hic.3DRAM_D45_Pax6.ins_250_domains_expanded",
                                  score_tracks = score_tracks,hic_names = c('NSCscore','IPCscore'),
                                  other_tracks =gtrack.ls('methylation','GpC','10x'),other_names = c('NSC GpC accessbility','IPC GpC Accessibility'))

ipc_mat <- extract_misha_ep_pairs(enhancer_bed = paste0(main_dir,'data/hg38/DMR_analysis/IPC_distal_DAR.bed'),
                                  tss_bed = paste0(main_dir,'data/hg38/beds/hg38_proteinCoding_TSS.bed'),
                                  interval_window = 200,min_dist = 5e3,max_dist = 2e6,expand = c(-1e4,1e4),
                                  domains_f = "hic.3DRAM_D45_Tbr2.ins_250_domains_expanded",
                                  score_tracks = score_tracks,hic_names = c('NSCscore','IPCscore'),
                                  other_tracks =gtrack.ls('methylation','GpC','10x'),other_names = c('NSC GpC accessbility','IPC GpC Accessibility'))


res = read.table('/home/hpc/bonev/projects/rna/ram/IPCvsNSC_res_proteinCoding_all.tsv')
ipc_mat$log2FoldChange <- res$log2FoldChange[match(ipc_mat$geneName,row.names(res))]
ipc_mat$RNA_padjust <- res$padj[match(ipc_mat$geneName,row.names(res))]

nsc_mat$log2FoldChange <- res$log2FoldChange[match(nsc_mat$geneName,row.names(res))]
nsc_mat$RNA_padjust <- res$padj[match(nsc_mat$geneName,row.names(res))]

test <- ipc_mat[ipc_mat$RNA_padjust<=0.05&ipc_mat$domain=='intraTAD',]
test <- test[!is.na(test$log2FoldChange),]
test_up <- test[test$log2FoldChange>=0.25,]
test_down <- test[test$log2FoldChange<=(-0.25),]
summary(test_down$NSCscore)
summary(test_down$IPCscore)

summary(test_up$NSCscore)
summary(test_up$IPCscore)

features_IPC <- c('chr3:27994033-EOMES','chr7:31603221-NEUROD6','chr7:31736471-PPP1R17','chr9:22326727-ELAVL2','chr1:60235197-NFIA','chr3:78498190-ROBO2')
features_NSC <- c('chr18:78721808-SALL3','chr10:16986658-VIM','chr9:86624737-GAS1','chr6:19644247-ID4','chr8:123720319-FBXO32','chr8:123658629-FBXO32','chr3:58581398-FAM107A','chr3:194296611-HES1')

xscore='NSCscore'
yscore='IPCscore'
df <- test_up
df_b <- melt(test_down[,c(xscore,yscore)])
df_b$variable <- gsub('score','',df_b$variable)
df_b$variable <- factor(df_b$variable,levels=c("NSC",'IPC'))
levels(df_b$variable) <- c('RGC','IPC')
p1 <- ggplot(df_b,aes(x=variable,y=value,fill=variable)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8,outlier.shape = NA) 
p1 <- p1 + scale_fill_manual(values=cell_colors) + xlab('') + ylab('Hi-C Score') + theme(legend.position = "none")
p1 <- p1 + stat_compare_means(comparisons = list(c('RGC','IPC')),label = "p.format",method='wilcox',paired = F)
pdf('figures/Figure5N_NSC_box.pdf',height=5,width=4)
p1 
dev.off()


p <- ggplot( df, aes_( x = as.name(xscore), y = as.name(yscore) )) + geom_pointdensity(shape=19,size=1,alpha=1) + scale_color_gradientn(colours = c('black','red','orange'),name='')
p <- p + geom_abline(slope = 1,intercept = 0,linetype = 2,col='darkblue') + xlim(0,100) + ylim(0,100) + theme(legend.position='none')
df <- df[order(df[,xscore]-df[,yscore],df[,xscore],decreasing=T),]
df$labels <- paste0(df$distalIDs,'-',df$geneName)
if (is.null(features)){
  df1 <- as.data.frame(df %>% group_by(geneName) %>% dplyr::slice(2))
  df1 <- df1[order(df1[,yscore]-df1[,xscore],df1[,yscore],decreasing=T),]
  df1 <- df1[df1[,yscore]>0&df1[,xscore]>0,]
  features1 <- c(head(df1$labels,5),features1)
} else {
  features1 <- features_NSC
}
p <- p + theme(legend.direction = "horizontal",legend.justification=c(0,0), legend.position=c(0.85, 0.05),legend.background =element_blank(),legend.key = element_blank(),legend.text = element_blank())+ guides(color = guide_colorbar(barwidth = 3, barheight = 1))
p <- p + geom_text_repel(
  data = df, size = 3,box.padding = 0.5,segment.alpha = 0.5, min.segment.length = 0,max.iter = 20000,
  aes(color=NULL,label=ifelse(labels%in%features1, as.character(labels), "")),force=10)
pdf('figures/Figure5N_NSC_sc.pdf',height=6,width=6)
p + xlab('RGC') + ylab('IPC')
dev.off()
