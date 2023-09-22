source('scripts/config.R')
source('scripts/aux_functions.R')
source('scripts/plot_functions.R')

source(paste0(main_dir,'results/hg38/hic/config.R'))
source(paste0(main_f,'scripts/main_functions.R'))
source(paste0(main_f,'scripts/aux_functions.R'))
source(paste0(main_f,'scripts/plot_functions.R'))
source(paste0(main_f,'scripts/temp_functions.R'))

library(motifmatchr)
library(TFBSTools)

library(doParallel)
library(corrplot)
registerDoParallel(cores=8)

misha_meth_rep_tracks <- c(paste0("methylation.3DRAM_Pax6_CpG_5x_rep",1:2),paste0("methylation.3DRAM_Tbr2_CpG_5x_rep",1:2),paste0("methylation.3DRAM_Pax6_GpC_5x_rep",1:2),paste0("methylation.3DRAM_Tbr2_GpC_5x_rep",1:2))
########################
###Functions###
########################

FigureS4B <- function(control_files,samples,cols=c('#aaa9ad','#aaa9ad')) {
  df <- read.table(control_files[grep(samples[1],control_files)],header=T,sep='\t')
  df <- melt(df,varnames = c('sample','Context'),measure.vars = 'methylation')
  df_summ <- data_summary(data=df, varname="value",groupnames=c("Context"))
  p <- ggplot(df_summ, aes(x=Context, y=value, fill=Context)) + scale_fill_manual(values = cols) +
    geom_bar(stat="identity", color="black",position=position_dodge()) + xlab('') + ylab('Methylation (%)') +
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.9)) + ggtitle(samples[1])
  p1 <- p + theme(legend.position='none',plot.title = element_text(hjust = 0.5)) + geom_point(data = df,position=position_jitter(w = 0.2, h = 0),size=1) + annotate("text", x = c(1,2), y = c(7,105), 
                                                                                                                                                                    label = round(df_summ$value,2))
  
  df <- read.table(control_files[grep(samples[2],control_files)],header=T,sep='\t')
  df <- melt(df,varnames = c('sample','Context'),measure.vars = 'methylation')
  df_summ <- data_summary(data=df, varname="value",groupnames=c("Context"))
  p <- ggplot(df_summ, aes(x=Context, y=value, fill=Context)) + scale_fill_manual(values = cols) +
    geom_bar(stat="identity", color="black",position=position_dodge()) + xlab('') + ylab('Methylation (%)') +
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.9)) + ylab('') + ggtitle(samples[2])
  p2 <- p + theme(legend.position='none',plot.title = element_text(hjust = 0.5)) + geom_point(data = df,position=position_jitter(w = 0.2, h = 0),size=1) + annotate("text", x = c(1,2), y = c(105,105), label = round(df_summ$value,2))
  return(p1|p2)
}


########################
###Call Functions###
########################

control_files <- list.files(paste0(main_dir,'data/hg38/RAM_controls/'),full.names = T)
p <- FigureS4B(control_files = control_files,samples=c('lambda','puc19'))
pdf('figures/FigureS4B.pdf',width=4,height=4,useDingbats = F)
print(p)
dev.off()

#FigS4E -> coverage in 100bp bins 
kb_data<-read.table(file=paste0(main_dir,'data/hg38/misc/100bp_coverage.txt'), header=T, sep='\t')
kb_data$read_number<-kb_data$read_number*1000
kb_data$Celltype <- factor(kb_data$Celltype, levels = c("PAX6+", "EOMES+"))
kb_data$Celltype <- revalue(kb_data$Celltype, c("PAX6+"="RGC", "EOMES+"="IPC"))

p <- ggplot(data = kb_data, aes(x = as.numeric(as.character(read_number)), y = mean, color= Celltype)) + geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.2)+ 
  labs(x='Number of reads', y='Coverage of 100bp bins (%)', color="Celltype") + geom_line(size=1)+ geom_point(size=1, colour='black')+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position=c(0.05,0.85)) + scale_color_manual(values=cell_colors)

pdf(file='figures/FigureS4C.pdf',width=4.5,height=4,useDingbats = F)
print(p)
dev.off()

#FigS4E+F Correlation 5mC and CpG
res <- misha_extract(tracks=misha_meth_rep_tracks,regions=gintervals.all(),iterator=1e4)
colnames(res)[4:11] <- c(paste0('NSC rep',1:2),paste0('IPC rep',1:2),paste0('NSC rep',1:2),paste0('IPC rep',1:2))

jpeg('figures/FigureS4D.jpg',height=1600,width=1600,res = 300)            #pdf file becomes too big
heatpairs(as.matrix(res[,4:7]),method='pearson',main = 'CpG methylation')
dev.off()

jpeg('figures/FigureS4E.jpg',height=1600,width=1600,res = 300)           #pdf file becomes too big
heatpairs(as.matrix(res[,8:11]),method='pearson',main = 'GpC accessibility')
dev.off()


#FigS4F Correlation HiC ->  Correlation matrix calculated with HiCRep SCC 10kb bins 
hic_files <- c('/home/fnoack/projects/3DRAM_horg/Pax6_rep1/3DRAM_Pax6_rep1.hic','/home/fnoack/projects/3DRAM_horg/Pax6_rep2/3DRAM_Pax6_rep2.hic','/home/fnoack/projects/3DRAM_horg/Tbr2_rep1/3DRAM_Tbr2_rep1.hic', '/home/fnoack/projects/3DRAM_horg/Tbr2_rep2/3DRAM_Tbr2_rep2.hic')
res_df <- correlate_hic(hic_files,bin=1e4,chrs = 1:21)            #adapted for aux funktion add that error get script otherwise didnt worked... .errorhandling = c("remove")
#dimnames(res_df) <- list(paste0('Pax6 rep',1:2),paste0('Tbr2 rep',1:2))
colnames(res_df) <-c('PAX6+ rep1', 'PAX6+ rep2','EOMES+ rep1','EOMES+ rep2')
rownames(res_df) <-c('PAX6+ rep1', 'PAX6+ rep2','EOMES+ rep1','EOMES+ rep2')
pdf('Figures/FigureS4F.pdf',height=5,width=5)
corrplot(res_df, method="color", col=colorRampPalette(heatmap_colors)(100), order="hclust",number.digits = 3, 
         addCoef.col = "black",cl.lim=c(-1,1),cl.length=3,addgrid.col = 'black',
         tl.col="black", tl.srt=45,number.font=1,bg='transparent'
)
dev.off()


#FigS4G plot GpC accessbility across Song et. al. 2020 ATAC peaks from sorted RG and IPC

peaks <- list.files(paste0(main_dir,'data/hg38/beds/song2020_IPC_RG/'),pattern = 'ATAC-seq.bed',full.names = T) #Full path to peaks
tracks <-list.files('/home/hpc/bonev/projects/ram/data/hg38/RAM_bigwigs/',full.names = T)



#FigS4G GpC and CpG methylation at CTCF peaks centered around CTCT motif
chip_f <- paste0(main_dir,"data/hg38/beds/fetal_CTCF_peaks_GSE116825/hglft_fetal_CTCF_merged.bed") #Full path to peaks
pwm_f <- readRDS(paste0(main_dir,'data/mm10/combined_pwm.RDS')) ### motif matrix
bed_dir=paste0(main_dir,'data/hg38/beds/fetal_CTCF_peaks_GSE116825/fetal_CTCF')
tracks <- list.files(paste0(main_dir,'data/hg38/RAM_bigwigs/'),pattern = c('cov10x.bw'),full.names = T)

centeredTF_plot(pwm=pwm_f,chip_f=chip_f,motif.name='Ctcf',genome='hg38',tracks=tracks[grep('CpG',tracks)],   #Figure 4F
                min=1000, max=1000, bin=10,lab=c('RGC','IPC'),con='CpG',ylim=c(0,100),height=5,width=5,
                colours=cell_colors,save_plot=T,bed_dir=bed_dir,out_file='figures/FigureS4G_1.pdf',label='CTCF',ylab='CpG Methylation (%')

centeredTF_plot(pwm=pwm_f,chip_f=chip_f,motif.name='Ctcf',genome='hg38',tracks=tracks[grep('GpC',tracks)],   #Figure 4F
                min=1000, max=1000, bin=10,lab=c('RGC','IPC'),con='GpC',ylim=c(0,60),height=5,width=5,
                colours=cell_colors,save_plot=T,bed_dir=bed_dir,out_file='figures/FigureS4G_2.pdf',label='CTCF',ylab='GpC Methylation (%)')




#FigureS4G_1
res <- getPlotSetArray(tracks=tracks[grep('GpC',tracks)],features=peaks[grep('RG\\.',peaks)], refgenome='hg38', type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 10,ignore_strand = T)
pdf(file='figures/FigureS4G_1.pdf',height=6,width=6)
plotAverage(plotset=res, labels = c('NSC','IPC'), xlim = NULL,
            ylim =c(7,23), main = NULL, xlab = "", ylab = paste0('%GpC Methylation'),
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = T, legend_ext = F, legend_pos = 'topright',
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = cell_colors)
axis(side = 1,labels=c('-1kb','NSC/RG peaks','+1kb'),at=c(-1000,0,1000),pos =7,tick = T)
dev.off()

res <- getPlotSetArray(tracks=tracks[grep('GpC',tracks)],features=peaks[grep('IPC\\.',peaks)], refgenome='hg38', type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 10,ignore_strand = T)
pdf(file='figures/FigureS4G_2.pdf',height=6,width=6)
plotAverage(plotset=res, labels = c('NSC','IPC'), xlim = NULL,
            ylim =c(7,23), main = NULL, xlab = "", ylab = paste0('%GpC Methylation'),
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = T, legend_ext = F, legend_pos = 'topright',
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = cell_colors)
axis(side = 1,labels=c('-1kb','IPC peaks','+1kb'),at=c(-1000,0,1000),pos =7,tick = T)
dev.off()

#FigS4H AggregatedHiC across CTCF chip-seq peaks 
# for (cell in cells){
#   submit_aggregateHiC(cells=cell,tracks=all_tracks,range_f=40000,res_f=1000,filter_f=0,
#                       intervals1='/home/hpc/bonev/projects/ram/data/hg38/beds/fetal_CTCF_peaks_GSE116825/all_CTCF_motif_centered_f.bed',
#                       intervals2='/home/hpc/bonev/projects/ram/data/hg38/beds/fetal_CTCF_peaks_GSE116825/all_CTCF_motif_centered_r.bed',
#                       grid_mode='1D',mem=80) 
# }

pdf('figures/FigureS4I_test.pdf',width=7.5,height=3)
layout(matrix(c(1:3),nrow=1,byrow=F),widths = c(4,4,1.5),heights=c(4),respect = T)
for (cell in cells){
  params <- plot_aggregateHiC(cells=cell,pool=T,range_f=40000,filter_f=0,res_f=1000,plot_res=2000,grid_mode='1D',zlim=c(-1,1),which_plot=c(2),
                              intervals1='all_gNOME_CTCF_F.bed',intervals2='all_gNOME_CTCF_R.bed',
                              interval1_name = '',interval2_name = '',
                              add_plot=T,plot_mean=T)    
  par(mar=c(1.5,2,1.5,4),mgp=c(1,0.5,0)) 
}
image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
dev.off()

source(paste0(main_dir,'results/mm10/hic/config.R'))

plotMisha(main_f=main_f,targetGene='Idh1',outDir='figures/',out_f='FigureS4J',upstream=2.2e5,downstream=7e4,chipYlim=matrix(c(0,7,0,7),nrow = 4,ncol = 2,byrow = T),
          chipNames=c('',''),window_scale=1.8,chipRes=10,pointCEX=5,conditions=c("hic.E14_NSC.score_k100","hic.E14_IPC.score_k100"),binSize=1e3,radius=1e4,
          chipTracksToExtract=c("scATAC.E14_NSC","scATAC.E14_IPC"),chipColors=c(cell_colors), 
          plotOrder=list(scores=TRUE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,anno=FALSE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.2, VP=1.5, loops=2.2, rna=0.6, chip=0.6,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.2))

run_cisDecay(tracks=all_tracks,cellsToPlot=cells,
             path=paste0(main_f,'analysis/cis_decay/'),log_base=10,cells=cells,
             file_f='cisDecay_Pax6_Tbr2',minDist=3e3,maxDist=1e8,width=6,height=6.5,labels=c('RGC','IPC'),
             colors=cell_colors,alpha=0.3,out_f='figures/FigureS4K.pdf')
