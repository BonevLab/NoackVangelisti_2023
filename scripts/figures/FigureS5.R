source('scripts/config.R')
source('scripts/aux_functions.R')
source('scripts/plot_functions.R')

source(paste0(main_dir,'results/hg38/hic/config.R'))
source(paste0(main_f,'scripts/main_functions.R'))
source(paste0(main_f,'scripts/aux_functions.R'))
source(paste0(main_f,'scripts/plot_functions.R'))
source(paste0(main_f,'scripts/temp_functions.R'))



tracks <- list.files(paste0(main_dir,'data/hg38/RAM_bigwigs/'),pattern = 'merged',full.names = T)

peaks <- paste0(main_dir,'data/hg38/DMR_analysis/GpC_hypo_0.05_0_over_disp_chi.bed')
res <- getPlotSetArray(tracks=tracks[grep('GpC',tracks)],features=peaks,refgenome='hg38',type = 'mf',add_heatmap=T,xmin=500,xmax=500,bin = 20,ignore_strand = T)
pdf('figures/FigureS5A.pdf',width=6,height=8)
plotHeatmap(res,labels = c('',''),sortrows = 'decreasing',clspace = rev(colorpalette('reds',12)),
            clstmethod=FALSE,raster = T)
dev.off()

peaks <- paste0(main_dir,'data/hg38/DMR_analysis/GpC_hyper_0.05_0_over_disp_chi.bed')
res <- getPlotSetArray(tracks=tracks[grep('GpC',tracks)],features=peaks,refgenome='hg38',type = 'mf',add_heatmap=T,xmin=500,xmax=500,bin = 20,ignore_strand = T)
pdf('figures/FigureS5B.pdf',width=6,height=8)
plotHeatmap(res,labels = c('',''),sortrows = 'decreasing',clspace = rev(colorpalette('reds',12)),
            clstmethod=FALSE,raster = T)
dev.off()


region <- toGRanges(paste0(main_dir,'data/hg38/beds/gNOME_peaks/mergedgNOME_sig_peaks.bed'))
files_f <- as.list(list.files(paste0(main_dir,'data/hg38/RAM_covfiles/replicates/'),pattern = 'CpG',full.names = T))
DMR_call(con = 'CpG',files = files_f,genome = 'hg38',threshold=10,difference=0,min_Cov=10,qval=0.05,
         high_cov_perc=99.9,statistical_method='over_disp_chi',region=region,output_path=paste0(main_dir,'data/hg38/DMR_analysis/'))


p <- plot_DMR(con='CpG',files=files_f,genome='hg38',threshold=10,difference=0,min_Cov=10,qval=0.05,labels=c('RGC CpG Methylation (%)','IPC CpG Methylation (%)'),
              high_cov_perc=99.9,statistical_method='over_disp_chi',region=region,output_path=paste0(main_dir,'data/hg38/DMR_analysis/'),colours=cell_colors) #adjust context rest as above ! 
pdf('figures/FigureS5C.pdf',width = 5,height=5,useDingbats = F)
print(p + xlab("RGC CpG Methylation (%)") + ylab("IPC CpG Methylation (%)") + theme(legend.position = "none"))
dev.off()

peaks <- paste0(main_dir,'data/hg38/DMR_analysis/CpG_hypo_0.05_0_over_disp_chi.bed')
background_peaks <- paste0(main_dir,'data/hg38/beds/gNOME_peaks/mergedgNOME_sig_peaks.bed')

NSC_motifs <- enrichedMotifs(peaks=peaks,bg_peaks = background_peaks,genome = 'hg38',genome_BS = BSgenome.Hsapiens.UCSC.hg38,
                             features=c('Neurog2(var.2)','Neurod1','Neurod2','Eomes','Nfia','Nfib','Fosb::jun','Lhx2','Sox2','Sox10','Tead2'),
                             cols=c("blue",'grey80','red'),logFC=0.25,logP=2,point.size=4,anno.size=6)
pdf(paste0('figures/FigureS5D.pdf'),width=5,height=5,useDingbats = F)
print(NSC_motifs$p)
dev.off()

peaks <- paste0(main_dir,'data/hg38/DMR_analysis/Peaks_TF_Centered/GpC_hyper_Neurog2var.2_0.0005_Neurog2var.2.bed')
res <- getPlotSetArray(tracks=tracks[grep('CpG',tracks)],features=peaks,refgenome='hg38',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 20,ignore_strand = F)
pdf('figures/FigureS5E_1.pdf',width=6,height=6)
plotAverage(plotset=res, labels = c('NSC','IPC'), xlim = NULL,
            ylim = c(40,90), main = '', xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = T, legend_ext = F, legend_pos = "bottomright",
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = cell_colors)
axis(side = 1,labels=c('-1kb','NEUROG2','+1kb'),at=c(-1000,0,1000),pos =40,tick = T)
dev.off()

peaks <- paste0(main_dir,'data/hg38/DMR_analysis/Peaks_TF_Centered/GpC_hyper_Neurog2var.2_0.0005_noNeurog2var.2.bed')
res <- getPlotSetArray(tracks=tracks[grep('CpG',tracks)],features=peaks,refgenome='hg38',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 20,ignore_strand = F)
pdf('figures/FigureS5E_2.pdf',width=6,height=6)
plotAverage(plotset=res, labels = c('NSC','IPC'), xlim = NULL,
            ylim = c(40,90), main = '', xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = T, legend_ext = F, legend_pos = "bottomright",
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = cell_colors)
axis(side = 1,labels=c('-1kb','no NEUROG2','+1kb'),at=c(-1000,0,1000),pos =40,tick = T)
dev.off()

pdf('figures/FigureS5F.pdf',width=7.5,height=3)
layout(matrix(c(1:3),nrow=1,byrow=F),widths = c(4,4,1.5),heights=c(4),respect = T)
for (cell in cells){
  params <- plot_aggregateHiC(cells=cell,pool=T,range_f=40000,filter_f=0,res_f=1000,plot_res=4000,grid_mode='1D',zlim=c(-0.4,0.4),which_plot=c(2),
                              intervals1='GpC_hyper_Neurog2var.2_0.0005_noNgn2.bed',intervals2='GpC_hyper_Neurog2var.2_0.0005_noNgn2.bed',
                              interval1_name = '',interval2_name = '',
                              add_plot=T,plot_mean=F)    
  par(mar=c(1.5,2,1.5,4),mgp=c(1,0.5,0)) 
}
image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
dev.off()

peaks <- paste0(main_dir,'data/hg38/DMR_analysis/Peaks_TF_Centered/GpC_hypo_Lhx2_0.0005.bed')
res <- getPlotSetArray(tracks=tracks[grep('CpG',tracks)],features=peaks,refgenome='hg38',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 20,ignore_strand = F)
pdf('figures/FigureS5G_1.pdf',width=6,height=6)
plotAverage(plotset=res, labels = c('NSC','IPC'), xlim = NULL,
            ylim = c(0,100), main = '', xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = T, legend_ext = F, legend_pos = "bottomright",
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = cell_colors)
axis(side = 1,labels=c('-1kb','LHX2','+1kb'),at=c(-1000,0,1000),pos =0,tick = T)
dev.off()

peaks <- paste0(main_dir,'data/hg38/DMR_analysis/Peaks_TF_Centered/GpC_hypo_Lhx2_0.0005_noLhx2.bed')
res <- getPlotSetArray(tracks=tracks[grep('CpG',tracks)],features=peaks,refgenome='hg38',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 20,ignore_strand = F)
pdf('figures/FigureS5G_2.pdf',width=6,height=6)
plotAverage(plotset=res, labels = c('NSC','IPC'), xlim = NULL,
            ylim = c(0,100), main = '', xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = T, legend_ext = F, legend_pos = "bottomright",
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = cell_colors)
axis(side = 1,labels=c('-1kb','no LHX2','+1kb'),at=c(-1000,0,1000),pos =00,tick = T)
dev.off()

pdf('figures/FigureS5H.pdf',width=7.5,height=3)
layout(matrix(c(1:3),nrow=1,byrow=F),widths = c(4,4,1.5),heights=c(4),respect = T)
for (cell in cells){
  params <- plot_aggregateHiC(cells=cell,pool=T,range_f=40000,filter_f=0,res_f=1000,plot_res=4000,grid_mode='1D',zlim=c(-0.4,0.4),which_plot=c(2),
                              intervals1='GpC_hypo_Lhx2_0.0005_noLhx2.bed',intervals2='GpC_hypo_Lhx2_0.0005_noLhx2.bed',
                              interval1_name = '',interval2_name = '',
                              add_plot=T,plot_mean=F)    
  par(mar=c(1.5,2,1.5,4),mgp=c(1,0.5,0)) 
}
image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
dev.off()

### Test CpG/GpC per motif
misha_meth_tracks <- c("methylation.3DRAM_Pax6_CpG_merged_10x","methylation.3DRAM_Tbr2_CpG_merged_10x","methylation.3DRAM_Pax6_GpC_merged_10x","methylation.3DRAM_Tbr2_GpC_merged_10x") 

res_list <- list()

bed_fs <- list.files('/home/hpc/bonev/projects/ram/results/hg38/beds/motif_ALLgNOME_all/',pattern = 'bed',full.names = T)
for (bed_f in bed_fs){
  motif_name <- gsub("\\.bed.*","",basename(bed_f))
  bed <- read.table(bed_f)
  colnames(bed) <- c('chrom','start','end')
  df <- misha_extract(misha_meth_tracks,bed,iterator = bed,mode = 'avg')
  
  colMeans(df[,4:7],na.rm=T)
  cor(df[,4],df[,6],use = 'complete.obs')
  cor(df[,5],df[,7],use = 'complete.obs')
}


chip_f <- paste0(main_dir,"data/hg38/beds/gNOME_peaks/mergedgNOME_sig_peaks.bed") #Full path to peaks
pwm_f <- readRDS(paste0(main_dir,'data/mm10/combined_pwm.RDS')) ### motif matrix
bed_dir=paste0(main_dir,'data/hg38/beds/gNOME_peaks/motif_centered/all_gNOME_')
tracks <- list.files(paste0(main_dir,'data/hg38/RAM_bigwigs'),pattern = c('merged'),full.names = T)

centeredTF_plot(pwm=pwm_f,chip_f=chip_f,motif.name='Neurog2(var.2)',genome='hg38',tracks=tracks[grep('GpC',tracks)],   #Figure 4F
                min=500, max=500, bin=10,lab=c('NSC','IPC'),con='GpC',ylim=c(5,35),height=6,width=6,
                colours=cell_colors,save_plot=T,bed_dir=bed_dir,out_file='figures/FigureS5I_1.pdf',label='Neurog2',ylab='%GpC Accessibility')
centeredTF_plot(pwm=pwm_f,chip_f=chip_f,motif.name='Neurog2(var.2)',genome='hg38',tracks=tracks[grep('GpC',tracks)],   #Figure 4F
                min=200, max=200, bin=4,lab=c('NSC','IPC'),con='GpC',ylim=NULL,height=6,width=6,
                colours=cell_colors,save_plot=F,bed_dir=bed_dir,out_file='figures/FigureS5I.pdf',label='Neurog2',ylab='%GpC Accessibility')

centeredTF_plot(pwm=pwm_f,chip_f=chip_f,motif.name='Neurog2(var.2)',genome='hg38',tracks=tracks[grep('CpG',tracks)],   #Figure 4F
                min=500, max=500, bin=20,lab=c('NSC','IPC'),con='CpG',ylim=c(42,82),height=6,width=6,
                colours=cell_colors,save_plot=T,bed_dir=bed_dir,out_file='figures/FigureS5I_2.pdf',label='Neurog2',ylab='%CpG Methylation')






