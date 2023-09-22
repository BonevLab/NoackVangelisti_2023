source('scripts/config.R')
source('scripts/aux_functions.R')
source('scripts/plot_functions.R')

source(paste0(main_dir,'results/mm10/hic/config.R'))
source(paste0(main_f,'scripts/main_functions.R'))
source(paste0(main_f,'scripts/aux_functions.R'))
source(paste0(main_f,'scripts/plot_functions.R'))
source(paste0(main_f,'scripts/temp_functions.R'))

###Input ####
rep_tracks <- list.files(paste0(main_dir,'data/mm10/RAM_bigwigs'),pattern = 'cov5x',full.names = T)
misha_meth_tracks <- c("methylation.ES_3DRAM_CpG_10x","methylation.ES_3DRAM_GpC_10x")
misha_meth_rep_tracks <- c(paste0("methylation.ES_3DRAM_CpG_5x_rep",1:3),paste0(paste0("methylation.ES_3DRAM_GpC_5x_rep",1:3)))
misha_rna_tracks <- "rnaseq_RPM.ES_3DRAM_merged"



########################
###Functions###
########################

Figure2A <- function(peaks1,peaks2,peaks3,peak_names,window=100,extend=50,cols=c("#E69F00", "green4", "#56B4E9",'grey57'),fig_f,key.size=1,width=4,height=4){
  peaks <- c(peaks1,peaks2,peaks3)
  peak_list <- lapply(peaks,function(x){
    df <- read.table(x)
    df <- df[df[,1] %in% gintervals.all()$chrom,]
    df <- gintervals(df[,1],df[,2],df[,3])
    df <- intervals.expand(df,extend)
    df <- gintervals.canonic(df)
    df$names <- paste0(basename(x),1:nrow(df))
    return(df)
  })
  names(peak_list) <- peak_names
  overlaps <- data.frame(Overlap=c('none','DHS','ATAC','both'),Values=NA,Percentage=NA)
  overlaps_1vs2 <- gintervals.neighbors(peak_list[[1]],peak_list[[2]])
  colnames(overlaps_1vs2)[8:9] <- c('names1','dist1')
  overlaps_1vs2vs3 <- gintervals.neighbors(overlaps_1vs2,peak_list[[3]])
  colnames(overlaps_1vs2vs3)[13:14] <- c('names2','dist2')
  overlaps[overlaps$Overlap=='both','Values'] <- sum(overlaps_1vs2vs3$dist1<=window&overlaps_1vs2vs3$dist2<=window)
  overlaps[overlaps$Overlap=='ATAC','Values'] <- sum(overlaps_1vs2vs3$dist1<=window&overlaps_1vs2vs3$dist2>window)
  overlaps[overlaps$Overlap=='DHS','Values'] <- sum(overlaps_1vs2vs3$dist1>window&overlaps_1vs2vs3$dist2<=window)
  overlaps[overlaps$Overlap=='none','Values'] <- sum(overlaps_1vs2vs3$dist1>window&overlaps_1vs2vs3$dist2>window)
  overlaps$Percentage <- round(overlaps$Values/sum(overlaps$Values)*100,2)
  overlaps$Overlap <- factor(overlaps$Overlap,levels=c('both','ATAC','DHS','none'))
  overlaps <-overlaps %>% arrange(desc(Overlap)) %>%
    mutate(lab.ypos = cumsum(Percentage) - 0.5*Percentage)
  p<- ggplot(overlaps, aes(x="", y=Percentage, fill=Overlap))+geom_bar(width = 1, stat = "identity", color = "black")+ scale_fill_manual(name='% overlap',values=cols)
  p <- p + coord_polar("y", start=0)+geom_text(aes(y = lab.ypos, label = Percentage), color = "white", size=6) + theme_void()
  p<- p + theme(legend.position=c(0.5,0.95),legend.direction = 'horizontal', legend.box = "horizontal") + guides(colour=guide_legend(override.aes = list(size = key.size),nrow = 1)) 
  pdf(fig_f,height=height,width=width)
  print(p)
  dev.off()
  }

########################
###Call Functions###
########################

#Fig2A 
Figure2A(peaks1=paste0(main_dir,'data/mm10/beds/gNOMEpeaks/gNOME_peaks.bed'),
         peaks2=paste0(main_dir,'data/mm10/ATAC_seq/ATAC_peaks.bed'),
         peaks3=paste0(main_dir,'data/mm10/DNase_ENCSR000CMW/DNA_E14_peaks_encode.bed'),
         peak_names=c('gNOME','ATAC','DHS'),window=100,extend=50,
         cols=c("#E69F00", "green4", "#56B4E9",'grey57'),fig_f='figures/Figure2A.pdf',key.size=1,width=4,height=4)

#Fig2B 
peaks <- c(paste0(main_dir,'data/mm10/beds/gNOMEpeaks/sig_NOMEPeaks_bedform.bed'),paste0(main_dir,'data/mm10/ATAC_seq/ATAC_peaks.bed'),paste0(main_dir,'data/mm10/DNase_ENCSR000CMW/DNA_E14_peaks_encode.bed'),paste0(main_dir,'data/mm10/ATAC_seq/ATAC_peaks_shuffle.bed'))
tracks <- list.files(paste0(main_dir,'data/mm10/RAM_bigwigs'),pattern = 'merged',full.names = T)
res <- getPlotSetArray(tracks=tracks[2],features=peaks,refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 5,ignore_strand = F)
pdf('figures/Figure2B.pdf',width=5,height=4.5)
plotAverage(plotset=res, labels = c('3DRAM','ATAC','DHS','Shuffled'), xlim = NULL,
            ylim = c(10,45), main = '', xlab = "", ylab = 'GpC accessibility (%)',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = F, legend_ext = F, legend_pos = "topright",
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 12, ln.v = FALSE, ln.h = NULL, pointsize = 12)
axis(side = 1,labels=c('-1kb','0','+1kb'),at=c(-1000,0,1000),pos=10,tick = T)
dev.off()

#Fig2C Heatmap signals across GpC peaks
tracks <- c(paste0(main_dir,'data/mm10/RAM_bigwigs/3DRAM_ES_250k_merged_GpC_cov10x.bw'),paste0(main_dir,'data/mm10/ATAC_seq/mES_GSE113592.bw'),paste0(main_dir,'data/mm10/DNase_ENCSR000CMW/DNAse_E14_merged.bw'),"/home/hpc/bonev/projects/ram/data/mm10/mnase_GSM1400766/mES_MNAse_SRR1302810.bw")
peaks <- paste0(main_dir,'data/mm10/beds/gNOMEpeaks/gNOME_peaks.bed')
res <- getPlotSetArray(tracks=tracks,features=peaks,refgenome='mm10',type = 'mf',add_heatmap=T,xmin=500,xmax=500,bin = 20,ignore_strand = T)
pdf('figures/Figure2C_alt.pdf',width=8,height=8)
plotHeatmap(res,labels = c('','','',''),sortrows = 'decreasing',clspace = rev(colorpalette('reds',12)),
            clstmethod=FALSE,raster = T)
dev.off()

#Fig2D -> coverage in 100bp bins 
kb_data<-read.table(file=paste0(main_dir,'data/mm10/misc/100bb_coverage_Ren_Ecker_WGBS.txt'), header=T, sep='\t',nrows = 63)
kb_data$read_number<-kb_data$read_number*1000
pdf(file='figures/Figure2D.pdf',width=4.5,height=4,useDingbats = F)
ggplot(data = kb_data, aes(x = as.numeric(as.character(read_number)), y = mean, color= reorder(Method,holder))) + geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.2)+ 
  labs(x='Number of reads', y='Coverage of 100bp bins (in %)', color="Method") + geom_line(size=1.5)+ geom_point(size=1, colour='black')+
  scale_color_npg(name='') + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position=c(0.05,0.9))
dev.off()

#Fig2E -> HiC stats of Ren Ecker and 3DRAM-seq -> put them together  
hic_stat<-as.data.frame(read.table(paste0(main_dir,'data/mm10/misc/HiC_stat_3DRAM_Ecker_Ren.txt'), fill=T, header=T,sep = "\t"))
levels(hic_stat$Method) <- c('3DRAM-seq','Methyl-3C','Methyl-HiC')
pdf(file='figures/Figure2E.pdf',height=4,width=4.5,useDingbats = F)
ggplot(hic_stat, aes( y=Mean, x=reorder(Feature, +Holder), fill=Method))+ 
        geom_bar(position="dodge", stat="identity",col='black')+ 
        geom_dotplot(aes(x= reorder(Feature, +Holder),y=Percentage),binaxis='y', stackdir='center' ,dotsize=1,position="dodge",binwidth = 1)+ 
        geom_errorbar(aes(ymin=Mean-STD, ymax=Mean+STD), width=.2,position=position_dodge(.9))+
        ylim(0,40) +
        labs(y= 'Percentage of Total Reads',x= "", fill= '')+
        scale_fill_grey(start = 0.25,end = 1)+ theme(legend.position=c(0.65,0.95)) + scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))
dev.off()



#Fig2F
pdf('figures/Figure2F.pdf',width=7,height=6)
layout(matrix(c(1,2,5,3,4,5),nrow=2,byrow=T),widths = c(4,4,1,4,4),heights=c(4,4,8,4,4),respect = T)
for (cell in cells){
  params <- plot_aggregateHiC(cells=cell,pool=T,range_f=40000,filter_f=0,res_f=1000,plot_res=1000,grid_mode='loops',zlim=c(-0.75,0.75),which_plot=c(2),
                            intervals1='ES_Bonev2017_loops_5000_pval0_05_L.bed',intervals2='ES_Bonev2017_loops_5000_pval0_05_R.bed',
                            interval1_name = '',interval2_name = '',
                            add_plot=T,plot_mean=F)
}
par(mar=c(1.5,1,1.5,3),mgp=c(1,0.5,0)) 
image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
dev.off()

#
plotMisha(main_f=main_f,targetGene='Zfp42',out_f='Figure2G',upstream=1e5,downstream=8e5,
          chipNames=c(''),window_scale=1.8,chipRes=20,pointCEX=2.5,conditions=score_tracks[1:4],chipColors='black',methNames=c('','','',''),
          methTracksToExtract=c("methylation.ES_3DRAM_CpG_10x","methylation.ES_Methyl3C_CpG_merged_10x","methylation.ES_MethylHiC_CpG_merged_10x","methylation.ES_WGBS_10x"),
          methColors=pal_npg()(4),
          chipTracksToExtract=c("chipseq_RPM.ES_CTCF"),
          plotOrder=list(scores=TRUE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,meth=TRUE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.5,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15))
