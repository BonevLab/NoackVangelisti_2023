source('scripts/config.R')
source('scripts/aux_functions.R')
source('scripts/plot_functions.R')

source(paste0(main_dir,'results/mm10/hic/config.R'))
source(paste0(main_f,'scripts/main_functions.R'))
source(paste0(main_f,'scripts/aux_functions.R'))
source(paste0(main_f,'scripts/plot_functions.R'))
source(paste0(main_f,'scripts/temp_functions.R'))

set.seed(123)
theme_set(theme_cowplot())
theme_border <- theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line=element_blank(),panel.border = element_rect(color = "black", size = 0.5),legend.margin = margin(0,0.5,0,0.5,unit='inch'))

tracks <- list.files(paste0(main_dir,'data/mm10/RAM_bigwigs'),pattern = 'merged',full.names = T)
misha_meth_tracks <- c("methylation.ES_3DRAM_CpG_10x","methylation.ES_3DRAM_GpC_10x")
misha_rna_tracks <- "rnaseq_RPM.ES_3DRAM_merged"
########################
###Functions###
########################

Figure1D <- function(res,zlim=c(-1,1),cols,extra_cols,path,file_f,ylim1,ylim2,xlim,fig_f,height,width){
  pdf(fig_f,height=height,width=width)
  layout(mat=matrix(c(1,2,3,5,4,6),byrow=T,ncol=2),widths = c(4,1,4,1,4,1),heights=c(1.5,0.5,0.5),respect = F)
  plot_averageTAD(cells=cells[1],fig_name = '',path=path,file_f=file_f,stats_f=F,zlim=zlim,
                  flip=T,z_colors=cols,add_plot = T,par_mar=c(1,2,1,0))
  par(mar=c(1,2,1,5))
  image.scale(as.matrix(1),zlim=zlim, col=cols,axis.pos=4,adj=0.5,cex.axis = 1.2)
  par(mar=c(1,2,1,0))
  plotMext(INPUTS=res[1]$data, xlim = xlim,
           ylim = ylim1, main = '', xlab = "", ylab = '',
           plotScale = "linear", type = "full", error.estimates = T,
           legend = F, legend_ext = F,bty='n',
           xaxt='n',yaxt='n',yaxs='i',xaxs='i',
           colvec = extra_cols[1],par_mar=c(0.5,2,0.5,0),
           ln.v = FALSE, ln.h = NULL)
  axis(side = 2,labels=ylim1,at=ylim1,pos=0,tick = T,las=2)
  plotMext(INPUTS=res[2]$data, xlim = xlim,
           ylim = ylim2, main = '', xlab = "", ylab = '',
           plotScale = "linear", type = "full", error.estimates = T,
           legend = F, legend_ext = F,bty='n',
           xaxt='n',yaxt='n',yaxs='i',xaxs='i',
           colvec = extra_cols[2],par_mar=c(0.5,2,0.5,0),
           ln.v = FALSE, ln.h = NULL)
  axis(side = 2,labels=ylim2,at=ylim2,pos=0,tick = T,las=2)
  dev.off()
}

########################

#Fig1B 
peaks <- paste0(main_dir,'data/mm10/beds/ES_CTCF_top30k_motifs.bed')
res <- getPlotSetArray(tracks=tracks,features=peaks,refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 5,ignore_strand = F)
pdf('figures/Figure1B.pdf',width=5,height=5)
plotAverage(plotset=res, labels = c('CpG methylation','GpC accessibility'), xlim = NULL,
            ylim = c(0,100), main = '', xlab = "", ylab = '',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = T, legend_ext = F, legend_pos = "topright",
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = c(cpg_color,gpc_color))
axis(side = 1,labels=c('-1kb','CTCF->','+1kb'),at=c(-1000,0,1000),pos =0,tick = T)
title(ylab="Methylation (%)",line = 2.5)
dev.off()


#Fig1C 
peaks <- paste0(main_dir,'data/mm10/beds/ES_Nrf1_motifs.bed')
res <- getPlotSetArray(tracks=tracks,features=peaks,refgenome='mm10',type = 'mf',add_heatmap=F,xmin=1000,xmax=1000,bin = 2,ignore_strand = F)
pdf('figures/Figure1C.pdf',width=5,height=5)
plotAverage(plotset=res, labels = c('CpG methylation','GpC accessibility'), xlim = NULL,
            ylim = c(0,100), main = '', xlab = "", ylab = '',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = T, legend_ext = F, legend_pos = "topright",
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = c(cpg_color,gpc_color))
axis(side = 1,labels=c('-1kb','Nrf1','+1kb'),at=c(-1000,0,1000),pos =0,tick = T)
title(ylab="Methylation (%)",line = 2.5)
dev.off()

res <- getPlotSetArray(tracks=tracks[2],features=peaks,refgenome='mm10',type = 'mf',add_heatmap=F,xmin=50,xmax=50,bin = 2,ignore_strand = F)
pdf('figures/Figure1C_i.pdf',width=4,height=4)
plotAverage(plotset=res, labels = , xlim = NULL,
            ylim = c(27,47), main = '', xlab = "", ylab = '',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',xaxs='i',
            legend = F, legend_ext = F, legend_pos = "topright",
            legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 14,xaxt='n',yaxt='n',
            cex.main = 14, cex.legend = 14, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = gpc_color)
axis(side = 1,labels=c('-50bp','Nrf1','+50bp'),at=c(-50,0,50),pos =27,tick = T,cex.axis=1.5)
axis(side = 2,labels=c(27,35,47),at=c(27,35,47),pos =-50,tick = T,cex.axis=1.5,las = 2)
dev.off()

#Fig1D 
zlim=c(-1,1)
domains <- paste0(main_dir,'results/mm10/hic/data/3DRAM_TADs_lengthExpanded.bed')
res <- getPlotSetArray(tracks=tracks,features=domains,refgenome='mm10',type = 'af',add_heatmap=F,xmin=1000,xmax=1000,bin = 100,ignore_strand = T,xanchored=10000)  #Extract CpG and CpG from TAD coordinates
averageTAD_colors <- c(blue_white_pal(length(seq(zlim[1],0.01,by=0.01))),'white',white_red_pal(length(seq(0.01,zlim[2],by=0.01))))

Figure1D(res=res,zlim=zlim,cols=averageTAD_colors,extra_cols=c(cpg_color,gpc_color),
         path=paste0(main_f,'analysis/averageTAD/'),file_f='averageTAD_3DRAM',
         ylim1=c(70,80),ylim2=c(10,13),xlim=c(0,10000),fig_f='figures/Figure1D.pdf',height=3,width=6)

#Fig1E 

pdf('figures/Figure1E.pdf',width=4,height=4)
layout(matrix(c(1:2),nrow=1,byrow=F),widths = c(4,2),heights=c(4),respect = T)
params <- plot_aggregateHiC(cells=cells[1],pool=T,range_f=40000,filter_f=0,res_f=1000,plot_res=1000,grid_mode='1D',zlim=c(-0.5,0.5),which_plot=c(2),
                            intervals1='ES_CTCF_for_top30k.bed',intervals2='ES_CTCF_rev_top30k.bed',
                            interval1_name = 'CTCF->',interval2_name = '<-CTCF',
                            add_plot=T,plot_mean=T)
par(mar=c(1.5,2,1.5,3),mgp=c(1,0.5,0)) 
image.scale.aggregate(params$input,zlim=params$zlim, col=params$cols,axis.pos=4,label='') 
dev.off()

#Fig1F

plotMisha(main_f=main_f,targetGene='Sox2',out_f='Figure1F',upstream=5e4,downstream=5e5,chipYlim=matrix(c(0,100,0,100,rep(NA,8)),nrow = 6,ncol = 2,byrow = T),
          chipNames=c('CpG methylation','GpC accessibility','ATAC','H3K27ac','RNA'),window_scale=1.8,chipRes=50,pointCEX=1,conditions=score_tracks[1],
          chipTracksToExtract=c(misha_meth_tracks,"atac.ES_GSE113592",'chipseq_RPM.ES_H3K27ac',misha_rna_tracks),chipColors=c(cpg_color,gpc_color,rep('black',4)),
          plotOrder=list(scores=TRUE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.5,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15))

plotMisha(main_f=main_f,targetGene='Sox2',out_f='Figure1G',upstream=1e4,downstream=1.2e5,
          chipNames=c('ATAC','H3K27ac','RNA'),window_scale=1.8,chipRes=10,pointCEX=5,conditions=score_tracks[1],methNames=c('CpG','GpC'),
          chipTracksToExtract=c("atac.ES_GSE113592",'chipseq_RPM.ES_H3K27ac',misha_rna_tracks),chipColors=c(rep('black',4)),methTracksToExtract=misha_meth_tracks,methColors=c(cpg_color,gpc_color),
          plotOrder=list(scores=TRUE,anno=FALSE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE,meth=TRUE, chip=TRUE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.5,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15))


