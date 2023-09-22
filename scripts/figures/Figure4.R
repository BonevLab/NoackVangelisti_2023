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
library(ggrastr)
library(org.Hs.eg.db)
library(clusterProfiler)
########################
###Functions###
#######################

Figure4C <- function(res,features,point.size=1,anno.size=4,sample_colors,out_f,height,width,ylim=NULL,xlim=NULL){
  res <- res[!is.na(row.names(res)),]
  res$gene_symbol <- row.names(res)
  p <- ggplot(res,aes(x=log2FoldChange,y=-log10(padj))) + geom_point_rast(color='grey',size = point.size) + geom_point(fill=sample_colors[2],alpha=0.8, pch = I(21),size = point.size,data = res[res$padj<0.05&res$log2FoldChange<0,]) + geom_point(fill=sample_colors[1],alpha=0.8, pch = I(21),size = point.size,data = res[res$padj<0.05&res$log2FoldChange>0,])  
  p <- p + xlab(expression(Log[2]~Fold~Change)) + ylab(expression(-Log[10]~(P)))
  p <- p + ggrepel::geom_text_repel(
    data = res, size = anno.size,seed = 42,
    box.padding =0.8, min.segment.length = 0,max.iter = 10000,
    aes(x=log2FoldChange,y=-log10(padj),color=NULL,label=ifelse(gene_symbol%in%features, as.character(gene_symbol), "")),force=10)
  if(!is.null(ylim)){
    p <- p+ ylim(ylim)
  }
  if(!is.null(xlim)){
    p <- p+ xlim(xlim)
  }
  pdf(paste0('figures/',out_f,'.pdf'),height=height,width=width)
  print(p)
  dev.off()
}  


Figure4D <- function(res,fc,direction,orgDB=org.Mm.eg.db,plot_what='GO',qvalue.cutoff=0.05,sample_colors,out_f,height,width){
  require(clusterProfiler)
  require(enrichplot)
  res <- res[!is.na(row.names(res)),]
  res$gene_symbol <- row.names(res)
  df <- res[res$padj<=0.05,]
  if(direction=='up'){
    df <- df[df$log2FoldChange>=fc,]
  } else if(direction=='down'){
    df <- df[df$log2FoldChange<=fc,]
  } else {
    df <- df[abs(df$log2FoldChange)>=fc,]
  }
  geneList <- as.numeric(df$log2FoldChange)
  names(geneList) <- as.character(df$gene_symbol)
  geneList <- geneList[order(geneList,decreasing = T)]
  names(geneList) <- bitr(names(geneList), fromType = "SYMBOL",
                          toType = c("ENTREZID"),
                          OrgDb = orgDB)$ENTREZID
  if(plot_what=='GO'){
    ego <- enrichGO(gene         = unique(names(geneList)),
                    OrgDb         = orgDB,
                    keyType       = 'ENTREZID',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = qvalue.cutoff)
    #p <- simplify(ego)
  } else if(plot_what=='GSEA'){
    ego <- gseGO(geneList     = geneList,
                 OrgDb        = orgDB,
                 ont          = "BP",
                 keyType       = 'ENTREZID',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
    #p <- simplify(ego)
  } else if(plot_what=='KEGG'){
    ego <- enrichKEGG(gene        = names(geneList),
                      organism     = 'mmu',
                      pvalueCutoff = 0.05)
  }
  return(list(ego=ego,geneList=geneList))
}



#Fig4D AggregatedTAd with GpC and CpG tracks !
Figure4H <- function(res,zlim=c(-1,1),z_cols,extra_cols,path,file_f,ylim1,ylim2,xlim,fig_f,height,width,cell){
  pdf(fig_f,height=height,width=width)
  layout(mat=matrix(c(1,2,3,5,4,6),byrow=T,ncol=2),widths = c(4,1,4,1,4,1),heights=c(1.5,0.5,0.5),respect = F)
  plot_averageTAD(cells=cell,fig_name = '',path=path,file_f=file_f,stats_f=F,zlim=zlim,
                  flip=T,z_colors=z_cols,add_plot = T,par_mar=c(1,2,1,0))
  par(mar=c(1,2.2,1,5))
  image.scale(as.matrix(1),zlim=zlim, col=z_cols,axis.pos=4,adj=0.5,cex.axis = 1.2)
  par(mar=c(1,2.2,1,0))
  plotMext(INPUTS=res[1]$data, xlim = xlim,
           ylim = ylim1, main = '', xlab = "", ylab = '',
           plotScale = "linear", type = "full", error.estimates = T,
           legend = F, legend_ext = F,bty='n',
           xaxt='n',yaxt='n',yaxs='i',xaxs='i',
           colvec = extra_cols[1],par_mar=c(0.5,2.2,0.5,0),
           ln.v = FALSE, ln.h = NULL)
  axis(side = 2,labels=ylim1,at=ylim1,pos=0,tick = T,las=2)
  plotMext(INPUTS=res[2]$data, xlim = xlim,
           ylim = ylim2, main = '', xlab = "", ylab = '',
           plotScale = "linear", type = "full", error.estimates = T,
           legend = F, legend_ext = F,bty='n',
           xaxt='n',yaxt='n',yaxs='i',xaxs='i',
           colvec = extra_cols[2],par_mar=c(0.5,2.2,0.5,0),
           ln.v = FALSE, ln.h = NULL)
  axis(side = 2,labels=ylim2,at=ylim2,pos=0,tick = T,las=2)
  dev.off()
}


########################
###Call Functions###
########################


#Fig4C
Figure4C(res = read.table('/home/hpc/bonev/projects/rna/ram/IPCvsNSC_res_proteinCoding_all.tsv'),xlim=c(-12,12),
         features <- c('SOX2','VIM','HES1','HOPX','GAS1','PAX6','NEUROG2','EOMES','NEUROD6','SSTR2','PCDH9','REST','LHX2','HES5','FOS','NEUROD2','NFIA','FBXO32'),
         point.size=2,sample_colors = rev(cell_colors),out_f = 'Figure4C',height=5,width=5)

#Fig4D
ego <- Figure4D(res = read.table('/home/hpc/bonev/projects/rna/ram/IPCvsNSC_res_proteinCoding_all.tsv'),fc = 0.5,plot_what='GO',
                direction = 'all',orgDB = org.Hs.eg.db,qvalue.cutoff = 0.05,sample_colors = rev(cell_colors))
p <- dotplot(simplify(ego$ego), showCategory=10,font.size=12)
p$theme <- theme_cowplot()
pdf('figures/Figure4D.pdf',height=6,width=8.75,useDingbats = F)
print(p+ theme(legend.position = c(0.8,0.35))) 
dev.off()


#Fig4 E-F -> Average plot of CpG and GpC methylation across motifs in all gnomePEAKS
chip_f <- paste0(main_dir,"data/hg38/beds/gNOME_peaks/mergedgNOME_sig_peaks.bed") #Full path to peaks
pwm_f <- readRDS(paste0(main_dir,'data/mm10/combined_pwm.RDS')) ### motif matrix
bed_dir=paste0(main_dir,'data/hg38/beds/gNOME_peaks/motif_centered/all_gNOME_')
tracks <- list.files(paste0(main_dir,'data/hg38/RAM_bigwigs'),pattern = c('merged'),full.names = T)

#Figure 4E
centeredTF_plot(pwm=pwm_f,chip_f=chip_f,motif.name='Ctcf',genome='hg38',tracks=tracks[grep('CpG',tracks)],   
                min=1000, max=1000, bin=10,lab=c('RGC','IPC'),con='CpG',ylim=c(0,100),height=5.5,width=5.5,
                colours=cell_colors,save_plot=T,bed_dir=bed_dir,out_file='figures/Figure4E.pdf',label='CTCF',ylab='CpG Methylation (%)')

#Figure 4F
centeredTF_plot(pwm=pwm_f,chip_f=chip_f,motif.name='Ctcf',genome='hg38',tracks=tracks[grep('GpC',tracks)],   
                min=1000, max=1000, bin=10,lab=c('RGC','IPC'),con='GpC',ylim=c(0,45),height=5.5,width=5.5,
                colours=cell_colors,save_plot=T,bed_dir=bed_dir,out_file='figures/Figure4F.pdf',label='CTCF',ylab='GpC Accessibility (%)')

#Fig4G
#generateBins <- extractBinned(tracks=all_tracks,cells=cells,chrs=chrs,path=path,binSize=binSize,file_f=file_f) #<- has to be done once 
# run over console since otherwisw it gives a connection error :(
plotBinned_CpG_GpC(extra_tracks=meth_tracks,cells=cells,init_cell=cells[1],chrs='chr3',binSize=2e5,
                   path=paste0(main_f,'analysis/compartments/'),file_f=paste0(main_f,'analysis/compartments/chr3_200kb'),
                   plot_what='obs',balance=T,fig_name='figures/Figure4G.pdf',width=4,height=4)

#Fig4H

celltypes=c('Pax6','Tbr2')
ylims=list(c(5.5,6.5),c(6.5,7.5))
zlim=c(-1.3,1)
z_cols <- c(blue_white_pal(length(seq(zlim[1],-0.01,by=0.01))),'white',white_red_pal(length(seq(0.01,zlim[2],by=0.01))))

for (i in 1:2){
  tracks <- list.files(paste0(main_dir,'data/hg38/RAM_bigwigs'),pattern = c(paste0(celltypes[i]),'merged'),full.names = T)[2:3]
  domains <- paste0(main_dir,'results/hg38/hic/data/3DRAM_',celltypes[i],'_TADs_lengthExpanded.bed')
  res <- getPlotSetArray(tracks=tracks,features=domains,refgenome='hg38',type = 'af',add_heatmap=F,xmin=1000,xmax=1000,bin = 100,ignore_strand = T,xanchored=10000)  #Extract CpG and CpG from TAD coordinates
  Figure4H(res=res,zlim=zlim,z_cols=z_cols,extra_cols=c(cpg_color,gpc_color),
           path=paste0(main_f,'analysis/averageTAD/'),file_f=paste0('averageTAD_3DRAM_all'),
           ylim1=c(83,90),ylim2=ylims[[i]],xlim=c(0,10000),fig_f=paste0('figures/Figure4H_',i,'.pdf'),height=3,width=6,cell=cells[i])
}


#Fig4I

misha_gpc_tracks <- c("methylation.3DRAM_Pax6_GpC_merged_10x","methylation.3DRAM_Tbr2_GpC_merged_10x") 
misha_cpc_tracks <- c("methylation.3DRAM_Pax6_CpG_merged_10x","methylation.3DRAM_Tbr2_CpG_merged_10x") 
misha_rna_tracks <- c("rnaseq.Pax6_HMGU1_Day45_RNA_merged","rnaseq.Tbr2_HMGU1_Day45_RNA_merged") 


anno <- read.table('/home/hpc/bonev/projects/ram/data/hg38/beds/gNOME_peaks/motif_centered/all_gNOME_Ctcf.bed')
anno <- gintervals(anno[,1],anno[,2],anno[,3])
anno <- intervals.normalize(anno,1000)   #Makes all intervals 1kb long
plotMisha(main_f=main_f,targetGene='SOX2',outDir='figures/',out_f='Figure4J_SOX2',upstream=1e5,downstream=1e6,
chipYlim=matrix(c(0,100,0,100,0,100,0,100,0,12,0,12),nrow = 12,ncol = 2,byrow = T),annIntervals=anno,
chipNames=c('','','','','',''),window_scale=1.8,chipRes=10,pointCEX=.5,conditions=score_tracks,binSize=1e3,
chipTracksToExtract=c(misha_gpc_tracks,misha_cpg_tracks,misha_rna_tracks),chipColors=rep(cell_colors,3),
plotOrder=list(scores=TRUE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,anno=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
plotRatios=list(unitHeight=120, scores=2.2, VP=1.5, loops=2.2, rna=0.6, chip=0.6,meth=0.5, domains=0.15, genes=1, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.2))
