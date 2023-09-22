source('scripts/config.R')
source('scripts/aux_functions.R')
source('scripts/plot_functions.R')

source(paste0(main_dir,'results/mm10/hic/config.R'))

theme_set(theme_cowplot())
theme_border <- theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line=element_blank(),panel.border = element_rect(color = "black", size = 0.5),legend.margin = margin(0,0.5,0,0.5,unit='inch'))

library(doParallel)
library(corrplot)
registerDoParallel(cores=8)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

####Input 
tracks <- list.files(paste0(main_dir,'data/mm10/RAM_bigwigs'),pattern = 'merged',full.names = T)
rep_tracks <- list.files(paste0(main_dir,'data/mm10/RAM_bigwigs'),pattern = 'cov5x',full.names = T)
misha_meth_tracks <- c("methylation.ES_3DRAM_CpG_10x","methylation.ES_3DRAM_GpC_10x")
misha_meth_rep_tracks <- c(paste0("methylation.ES_3DRAM_CpG_5x_rep",1:3),paste0(paste0("methylation.ES_3DRAM_GpC_5x_rep",1:3)))
misha_rna_tracks <- "rnaseq_RPM.ES_3DRAM_merged"

########################
###Functions###
########################
#FigS1A
FigureS1A<-function(amplicons,reads_f,cov_threshold,center_dist=50) {
  reads_list <- lapply(reads_f,read.table)
  names(reads_list) <- gsub('\\..*','',basename(reads_f))
  df <- do.call(rbind,reads_list)
  colnames(df)[1:6] <- c('chrom','start','end','MethPer','Meth','UnMeth')
  df$start <- df$start-1   #account for -1 pos in coverage files
  df$condition <- factor(gsub('\\..*','',row.names(df)),levels=gsub('\\..*','',basename(reads_f)))
  levels(df$condition) <- c('10min_live','10min','1h','2h','3h','4h')
  colnames(amplicons) <- c('chrom','start','end','amplicon')
  if(center_dist!=0){
    amplicons <- intervals.centers(amplicons)
  }
  df <- gintervals.neighbors(df,amplicons,maxdist = center_dist)
  df <- df[rowSums(df[,c('Meth','UnMeth')])>=cov_threshold,c(1:7,11)]
  df_m <- melt(df[,c('condition','amplicon','MethPer')])
  p <- ggplot(df, aes(x=condition, y=MethPer)) + geom_violin(draw_quantiles = c(0.25, 0.75),linetype = "dashed",fill='#aaa9ad') +
    geom_violin(fill="transparent",draw_quantiles = 0.5) +
    stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red") + ylim(c(0,100)) + xlab('') + ylab('GpC Methylation (%)')
  return(p)
}

#FigS1B -> Per gene coverage was calculated by RSeQC 
FigureS1B<- function (rseqc_f,cols,anno.size=12) {
  df <- read.table(rseqc_f,header = T,row.names = 1)
  row.names(df) <- paste0('rep',1:3)
  colnames(df) <- 1:ncol(df)
  df <- t(df/rowSums(df))
  df <- melt(df)
  p <- ggplot(df,aes(x=Var1,y=value,color=Var2)) + geom_line() + scale_y_continuous(labels = scales::rescale) + ylab('Coverage') + xlab("Gene body percentile (5'−>3')") + scale_color_manual(name='',values=as.character(cols))  
  p <- p + theme(legend.text = element_text(size=anno.size)) + theme(legend.position=c(0.85,0.95), legend.box = "horizontal",legend.spacing.x = unit(0.1, 'inch')) 
  return(p)
}

#FigS1C
FigureS1C <- function(control_files,samples,cols=c('#aaa9ad','#aaa9ad')) {
  df <- read.table(control_files[grep(samples[1],control_files)],header=T,sep='\t')
  df <- melt(df,varnames = c('sample','Context'),measure.vars = 'methylation')
  df_summ <- data_summary(df, varname="value",groupnames=c("Context"))
  p <- ggplot(df_summ, aes(x=Context, y=value, fill=Context)) + scale_fill_manual(values = cols) +
    geom_bar(stat="identity", color="black",position=position_dodge()) + xlab('') + ylab('Methylation (%)') +
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.9)) + ggtitle(samples[1])
  p1 <- p + theme(legend.position='none',plot.title = element_text(hjust = 0.5)) + geom_point(data = df,position=position_jitter(w = 0.1, h = 0),size=2) + annotate("text", x = c(1,2), y = c(7,105), label = round(df_summ$value,2))
  
  df <- read.table(control_files[grep(samples[2],control_files)],header=T,sep='\t')
  df <- melt(df,varnames = c('sample','Context'),measure.vars = 'methylation')
  df_summ <- data_summary(df, varname="value",groupnames=c("Context"))
  p <- ggplot(df_summ, aes(x=Context, y=value, fill=Context)) + scale_fill_manual(values = cols) +
    geom_bar(stat="identity", color="black",position=position_dodge()) + xlab('') + ylab('Methylation (%)') +
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.9)) + ylab('') + ggtitle(samples[2])
  p2 <- p + theme(legend.position='none',plot.title = element_text(hjust = 0.5)) + geom_point(data = df,position=position_jitter(w = 0.1, h = 0),size=2) + annotate("text", x = c(1,2), y = c(105,105), label = round(df_summ$value,2))
  return(p1|p2)
}

#FigureS1D
FigureS1D<- function (RNA_data,cormethod,main='',fig_f,height,width) {
  df <- read.table(file =RNA_data, header =T , sep='\t')
  df <- as.matrix(log10(df[,1:3]+1))
  colnames(df) <- c('rep1','rep2','rep3')
  jpeg(fig_f,height=height,width=width,res = 300)   
  print(heatpairs(as.matrix(df),method=cormethod,main = main))
  dev.off()
}

#FigureS1G

########################
###Call Functions###
########################

#FigS1A -> BSamp-seq timecourse

amplicons <- read.table(paste0(main_dir,'data/mm10/NOME_BSamp_timecourse/analysis/amplicons.txt'), fill=T, header=F,sep = "\t")
reads_f <- list.files(paste0(main_dir,'data/mm10/NOME_BSamp_timecourse/analysis/'), pattern = paste0('GpC','.cov'), all.files = T,full.names = T)
reads_f <- reads_f[c(2,1,3:6)]           #Reordering for factor names
p <- FigureS1A(amplicons,reads_f,cov_threshold=50,center_dist = 50) + geom_vline(xintercept = 1.5,lty=2)
p <- p + stat_compare_means(comparisons = list(c(2,3),c(3,4),c(4,5),c(5,6)),label = "p.format",method='t.test',label.y=c(65,80,88,95))
pdf('figures/FigureS1A.pdf',width=5,height=4,useDingbats = F)
print(p)
dev.off()

#FigS1B -> Per gene coverage was calculated by RSeQC 
p <- FigureS1B(rseqc_f=paste0(main_dir,'data/mm10/RAM_RNA/RNA_quality_3DRAM_RNA.pdf.geneBodyCoverage.txt'),cols=rep_colors,anno.size=12)
ggsave('figures/FigureS1B.pdf',p,width=6,height=4)

#FigS1C puc and lambda control plots 
control_files <- list.files(paste0(main_dir,'data/mm10/RAM_controls/'),full.names = T)
p <- FigureS1C(control_files = control_files,samples=c('lambda','puc19'))
pdf('figures/FigureS1C.pdf',width=4,height=4,useDingbats = F)
print(p)
dev.off()

#FigS1D Correlation RNA_Seq
FigureS1D(RNA_data=paste0(main_dir,'data/mm10/RAM_RNA/FPKM_means.txt'),cormethod='spearman',fig_f = 'figures/FigureS1D.jpg',width=1600,height=1600,main='Gene expression')

#FigS1E+F Correlation 5mC and CpG
res <- misha_extract(tracks=misha_meth_rep_tracks,regions=gintervals.all(),iterator=1e4)
colnames(res)[4:9] <- c(paste0('rep',1:3),paste0('rep',1:3))

jpeg('figures/FigureS1E.jpg',height=1600,width=1600,res = 300)            #pdf file becomes too big
heatpairs(as.matrix(res[,7:9]),method='pearson',main = 'GpC accessibility')
dev.off()

jpeg('figures/FigureS1F.jpg',height=1600,width=1600,res = 300)           #pdf file becomes too big
heatpairs(as.matrix(res[,4:6]),method='pearson',main = 'CpG methylation')
dev.off()

#FigS1G Correlation HiC ->  Correlation matrix calculated with HiCRep SCC 10kb bins 
hic_files <- c('/home/fnoack/projects/3D_RAMseq_250k/rep1/3DRAM_ES_250k_rep1.hic','/home/fnoack/projects/3D_RAMseq_250k/rep2/3DRAM_ES_250k_rep2.hic','/home/fnoack/projects/3D_RAMseq_250k/rep3/3DRAM_ES_250k_rep3.hic')
res_df <- correlate_hic(hic_files,1e4)            #in aux_functions.R
dimnames(res_df) <- list(paste0('rep',1:3),paste0('rep',1:3))
pdf('figures/FigureS1G.pdf',height=4,width=4)
corrplot(res_df, method="color", col=colorRampPalette(heatmap_colors)(100), order="hclust",number.digits = 3, 
         addCoef.col = "black",cl.lim=c(-1,1),cl.length=3,addgrid.col = 'black',
         tl.col="black", tl.srt=0,number.font=1,bg='transparent'
)
dev.off()

#FigS1H Plot methylation and accessbility levels across TSS of repressed and 25% high expressed genes filtered for FPKM>1

generate_rna_quantile(tss_f=paste0(main_dir,'data/mm10/beds/mm10_TSS.bed'),                     #in aux_functions.R
                      fpkm_f=paste0(main_dir,'data/mm10/RAM_RNA/FPKM_means.txt'))

peaks <- c(paste0(main_dir,'data/mm10/beds/ES_TSS_inactive.bed'),paste0(main_dir,'data/mm10/beds/ES_TSS_Q4.bed'))
res <- getPlotSetArray(tracks=tracks[1],features=peaks,refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 20,ignore_strand = F)
pdf('figures/FigureS1H_1.pdf',width=5,height=5)
plotAverage(plotset=res, labels = c('repressed','highly expressed'), xlim = NULL,
            ylim = c(0,80), main = '', xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = F, legend_ext = F, legend_pos = "topright",
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = c('darkblue','darkred'))
axis(side = 1,labels=c('-2kb','TSS','+2kb'),at=c(-2000,0,2000),pos =0,tick = T)
dev.off()


res <- getPlotSetArray(tracks=tracks[2],features=peaks,refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 20,ignore_strand = F)
pdf('figures/FigureS1H_2.pdf',width=5,height=5)
plotAverage(plotset=res, labels = c('repressed','highly expressed'), xlim = NULL,
            ylim = c(10,70), main = '', xlab = "", ylab = '% GpC Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = T, legend_ext = F, legend_pos = "topright",
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = c('darkblue','darkred'))
axis(side = 1,labels=c('-2kb','TSS','+2kb'),at=c(-2000,0,2000),pos =10,tick = T)
dev.off()

peaks <- paste0(main_dir,'data/mm10/beds/ES_CTCF_top30k_motifs.bed')
res <- getPlotSetArray(tracks=rep_tracks[grep('CpG',rep_tracks)],features=peaks,refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 5,ignore_strand = F)
pdf('figures/FigureS1I_1.pdf',width=5,height=5)
plotAverage(plotset=res, labels = paste0('rep',1:3), xlim = NULL,
            ylim = c(0,80), main = '', xlab = "", ylab = '% CpG Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = T, legend_ext = F, legend_pos = "bottomright",
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12)
axis(side = 1,labels=c('-1kb','CTCF->','+1kb'),at=c(-1000,0,1000),pos =0,tick = T)
dev.off()

res <- getPlotSetArray(tracks=rep_tracks[grep('GpC',rep_tracks)],features=peaks,refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 5,ignore_strand = F)
pdf('figures/FigureS1I_2.pdf',width=5,height=5)
plotAverage(plotset=res, labels = paste0('rep',1:3), xlim = NULL,
            ylim = c(5,60), main = '', xlab = "", ylab = '% GpC Methylation',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = T, legend_ext = F, legend_pos = "topright",
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12)
axis(side = 1,labels=c('-1kb','CTCF->','+1kb'),at=c(-1000,0,1000),pos =5,tick = T)
dev.off()


peaks <- paste0(main_dir,'data/mm10/beds/ES_Nrf1_motifs.bed')
res <- getPlotSetArray(tracks=paste0(main_dir,'data/mm10/ATAC_seq/mES_GSE113592.bw'),features=peaks,refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 2,ignore_strand = F)
pdf('figures/FigureS1J.pdf',width=5,height=5)
plotAverage(plotset=res, labels = c('ATAC accessibility'), xlim = NULL,
            ylim = c(2,40), main = 'ATAC Accessibility', xlab = "", ylab = 'RPM',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = F, legend_ext = F, legend_pos = "topright",
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = 'black')
axis(side = 1,labels=c('-1kb','Nrf1','+1kb'),at=c(-1000,0,1000),pos =2,tick = T)
dev.off()

res <- getPlotSetArray(tracks=paste0(main_dir,'data/mm10/ATAC_seq/mES_GSE113592.bw'),features=peaks,refgenome='mm10',type = 'mf',add_heatmap=F,xmin=50,xmax=50,bin = 2,ignore_strand = F)
pdf('figures/FigureS1J_i.pdf',width=3,height=3)
plotAverage(plotset=res, labels = c('ATAC accessibility'), xlim = NULL,
            ylim = c(31,35), main = '', xlab = "", ylab = '',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = F, legend_ext = F, legend_pos = "topright",
            legend_ext_pos = "topleft", cex.axis = 10, cex.lab = 10,xaxt='n',
            cex.main = 10, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = 'black')
axis(side = 1,labels=c('-50bp','Nrf1','+50bp'),at=c(-50,0,50),pos =31,tick = T,cex.axis=0.8)
dev.off()

