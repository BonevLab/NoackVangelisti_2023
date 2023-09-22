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
library(clusterProfiler)
library(org.Hs.eg.db)
library(JASPAR2020)
library(Matrix)
library(TFBSTools)
library(motifmatchr)
library(Hmisc)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
library(monaLisa)

TE_beds <- paste0(main_dir,'data/hg38/repeats/bed_files/hg38_repeats_cleaned.bed')
peaks <- read.table(TE_beds)
colnames(peaks) <- c('chrom','start','end','Class','TE_class','strand')
hic_names=c('NSCscore','IPCscore')
other_tracks <- c("methylation.3DRAM_Pax6_CpG_merged_10x","methylation.3DRAM_Tbr2_CpG_merged_10x","methylation.3DRAM_Pax6_GpC_merged_10x","methylation.3DRAM_Tbr2_GpC_merged_10x","hic.3DRAM_D45_Pax6.ins_250","hic.3DRAM_D45_Tbr2.ins_250")
other_names <- c("NSC_CpG","IPC_CpG","NSC_GpC","IPC_GpC","NSC_ins","IPC_ins")

pwms <- getMatrixSet(JASPAR2020,
                     opts = list(matrixtype = "PWM",
                                 species = 9606))

monalisa_peaks <- function(peaks_f,pwms,p.cutoff){
  peaks <- read.table(peaks_f)
  colnames(peaks) <- c('chrom','start','end','family','subfamily','strand')
  peaks <- makeGRangesFromDataFrame(peaks)
  peakSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, peaks)
  
  
  se <- calcBinnedMotifEnrR(seqs = peakSeq,pwmL = pwms,background = 'genome',genome = BSgenome.Hsapiens.UCSC.hg38,genome.oversample = 10)
  sel <- apply(assay(se, "negLog10Padj"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > p.cutoff
  seSel <- se[sel, ]
  SimMatSel <- motifSimilarity(rowData(seSel)$motif.pfm)
  hcl <- hclust(as.dist(1 - SimMatSel), method = "average")
  
  # peaks <- read.table(peaks_f)[,1:3]
  # colnames(peaks) <- c('chrom','start','end')
  # #peaks <- intervals.normalize(peaks,200)
  # df <- misha_extract(tracks=c("methylation.3DRAM_Pax6_GpC_merged_10x","methylation.3DRAM_Tbr2_GpC_merged_10x"),regions = peaks,iterator = peaks,mode = 'avg')
  # ###Predict TFs responsible
  # peaks <- makeGRangesFromDataFrame(df)
  # peakSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, peaks)
  # hits <- findMotifHits(query = pwms, subject = peakSeq, min.score = 10.0,
  #                       BPPARAM = BiocParallel::MulticoreParam(2))
  # TFBSmatrix <- unclass(table(factor(seqnames(hits), levels = seqlevels(hits)),
  #                             factor(hits$pwmname, levels = name(pwms))))
  # zero_TF <- colSums(TFBSmatrix) == 0
  # TFBSmatrix <- TFBSmatrix[, !zero_TF]
  # fMono <- oligonucleotideFrequency(peakSeq, width = 1L, as.prob = TRUE)
  # fDi <- oligonucleotideFrequency(peakSeq, width = 2L, as.prob = TRUE)
  # fracGC <- fMono[,"G"] + fMono[,"C"]
  # oeCpG <- (fDi[,"CG"] + 0.01) / (fMono[,"G"] * fMono[,"C"] + 0.01)
  # TFBSmatrix <- cbind(fracGC, oeCpG, TFBSmatrix)
  # set.seed(42)
  # se <- randLassoStabSel(x = TFBSmatrix, y = df$v_methylation.3DRAM_Tbr2_GpC_merged_10x-df$v_methylation.3DRAM_Pax6_GpC_merged_10x,weakness = 1,cutoff = 0.6, mc.cores = 1)
  return(list(seSel=seSel,hcl=hcl))
}

if(!file.exists(paste0(main_dir,'results/hg38/repeats/repeats_analysis.tsv'))){
  res_f <- list()
  res_f <- foreach(te_class=unique(peaks$TE_class),
                   .export=c('peaks','hic_names','other_tracks','other_names'),
                   .errorhandling='remove') %dopar% {
                     res <- repeat_analysis(repeat_bed=peaks[peaks$TE_class==te_class,],interval_window=0,min_dist=5e4,max_dist=10e6,expand=c(-1e4,1e4),domains_f="hic.3DRAM_D45_Pax6.ins_250_domains_expanded",
                                            score_tracks=score_tracks,hic_names=hic_names,
                                            other_tracks=other_tracks,other_names=other_names)
                     intra_hic_df <- res$hicscores[res$hicscores$domain=='intraTAD',grep('score',colnames(res$hicscores))]
                     intra_hic_df <- intra_hic_df[complete.cases(intra_hic_df),]
                     inter_hic_df <- res$hicscores[res$hicscores$domain=='interTAD',grep('score',colnames(res$hicscores))]
                     inter_hic_df <- inter_hic_df[complete.cases(inter_hic_df),]
                     
                     other_df <- res$lin_scores[,-c(1:3,ncol(res$lin_scores))]
                     intra_res_hic <- c(colMeans(intra_hic_df,na.rm=T),hicsign=wilcox.test(intra_hic_df[,1],intra_hic_df[,2],paired = T)$p.value)
                     inter_res_hic <- c(colMeans(inter_hic_df,na.rm=T),hicsign=wilcox.test(inter_hic_df[,1],inter_hic_df[,2],paired = T)$p.value)
                     res_cpg <- c(colMeans(other_df[,grep('CpG',colnames(other_df))],na.rm=T),cpg_sign=wilcox.test(other_df[,1],other_df[,2],paired = T)$p.value)
                     res_gpc <- c(colMeans(other_df[,grep('GpC',colnames(other_df))],na.rm=T),gpc_sign=wilcox.test(other_df[,grep('GpC',colnames(other_df))[1]],other_df[,grep('GpC',colnames(other_df))[2]],paired = T)$p.value)
                     res_ins <- c(colMeans(other_df[,grep('ins',colnames(other_df))],na.rm=T),ins_sign=wilcox.test(other_df[,grep('ins',colnames(other_df))[1]],other_df[,grep('ins',colnames(other_df))[2]],paired = T)$p.value)
                     res <- c(as.character(te_class),intra_res_hic,nrow(intra_hic_df),inter_res_hic,nrow(inter_hic_df),res_cpg,sum(complete.cases(other_df[,grep('CpG',colnames(other_df))])),res_gpc,sum(complete.cases(other_df[,grep('GpC',colnames(other_df))])),res_ins,sum(complete.cases(other_df[,grep('ins',colnames(other_df))])))
                     names(res)[c(1,5,9,13,17,21)] <- c('TE','N_coveredIntraHiC_pairs','N_coveredInterHiC_pairs','N_coveredCpG','N_coveredGpC','N_coveredINS')
                     return(res)
                   }
  res <- as.data.frame(do.call("rbind",res_f))
  row.names(res) <- res$TE
  res <- res[,-1]
  res_path <- paste0(main_dir,'results/hg38/repeats/repeats_analysis.tsv')
  write.table(res,res_path,quote =F,col.names=T,row.names=T,sep='\t')
} else {
  res <- read.table(paste0(main_dir,'results/hg38/repeats/repeats_analysis.tsv'),header=T)
  res_ctcf <- read.table(paste0(main_dir,'results/hg38/repeats/repeats_analysis_ctcf.tsv'),header=T)
  res_shuffle <- read.table(paste0(main_dir,'results/hg38/repeats/repeats_analysis_shuffle.tsv'),header=T)
}

res <- rbind(res,res_shuffle)
res$TE_family <- peaks$Class[match(row.names(res),peaks$TE_class)]

#Figures

df <- plot_repeat_scatter(res,metric=c('GpC'),min_n = 50,diff.cutoff = 5,anno.size = 4,cols = c('grey',cell_colors),
                          features=c('MER130','UCON31'))
pdf('figures/Figure6A.pdf',height=5,width=5,useDingbats = F)
print(df$p + geom_smooth(method = "lm", se = FALSE) + xlab('RGC GpC Accessibility (%)') + ylab('IPC GpC Accessibility (%)'))
dev.off()


df <- plot_repeat_scatter(res,metric=c('CpG'),min_n = 30,diff.cutoff = 5,anno.size = 4,cols = c('grey',cell_colors),
                          features=c('MER130','UCON31'))
pdf('figures/Figure6B.pdf',height=5,width=5,useDingbats = F)
print(df$p + geom_smooth(method = "lm", se = FALSE) + xlab('RGC CpG Methylation (%)') + ylab('IPC CpG Methylation (%)'))
dev.off()

##Figure 6C###
peaks <- read.table(paste0(main_dir,'data/hg38/repeats/bed_files/MER130.bed'))[,1:3]
colnames(peaks) <- colnames(gintervals.all())
tracks=c("methylation.3DRAM_Pax6_GpC_merged_10x","methylation.3DRAM_Tbr2_GpC_merged_10x","methylation.3DRAM_Pax6_CpG_merged_10x","methylation.3DRAM_Tbr2_CpG_merged_10x")
df <- misha_extract(tracks = tracks,regions=peaks,iterator = peaks,mode = 'avg')
df <- df[,4:8]
colnames(df)[1:4] <- c('RGC_GpC','IPC_GpC','RGC_CpG','IPC_CpG')
df_melt <- melt(df,id.vars = 'intervalID')
df_melt$condition <- factor(gsub("_.*","",df_melt$variable),levels=c('RGC','IPC'))
df_melt$mark <- factor(gsub(".*_","",gsub('_distal','',df_melt$variable)),levels=c('GpC','CpG'))
p1 <- ggplot(df_melt,aes(x=mark,y=value,fill=condition)) + geom_boxplot(outlier.size=1,outlier.shape=NA,show.legend = T,width=0.8) 
p1 <- p1 + scale_fill_manual(values=cell_colors,name='') + xlab('') + ylab('Methylation (%)') + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.35, 0.9))
p1 <- p1 + stat_compare_means(label = "p.format",method='wilcox',paired = F)
pdf(paste0('figures/Figure6C.pdf'),width=5,height=5)
print(p1)
dev.off()
##Figure 6D###
peaks <- read.table(paste0(main_dir,'data/hg38/repeats/bed_files/UCON31.bed'))[,1:3]
colnames(peaks) <- colnames(gintervals.all())
tracks=c("methylation.3DRAM_Pax6_GpC_merged_10x","methylation.3DRAM_Tbr2_GpC_merged_10x","methylation.3DRAM_Pax6_CpG_merged_10x","methylation.3DRAM_Tbr2_CpG_merged_10x")
df <- misha_extract(tracks = tracks,regions=peaks,iterator = peaks,mode = 'avg')
df <- df[,4:8]
colnames(df)[1:4] <- c('RGC_GpC','IPC_GpC','RGC_CpG','IPC_CpG')
df_melt <- melt(df,id.vars = 'intervalID')
df_melt$condition <- factor(gsub("_.*","",df_melt$variable),levels=c('RGC','IPC'))
df_melt$mark <- factor(gsub(".*_","",gsub('_distal','',df_melt$variable)),levels=c('GpC','CpG'))
p1 <- ggplot(df_melt,aes(x=mark,y=value,fill=condition)) + geom_boxplot(outlier.size=1,outlier.shape=NA,show.legend = T,width=0.8) 
p1 <- p1 + scale_fill_manual(values=cell_colors,name='') + xlab('') + ylab('Methylation (%)') + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.35, 0.9))
p1 <- p1 + stat_compare_means(label = "p.format",method='wilcox',paired = F)
pdf(paste0('figures/Figure6D.pdf'),width=5,height=5,useDingbats = F)
print(p1)
dev.off()
##Figure 6E###

ucon31_mat <- extract_misha_ep_pairs(enhancer_bed = paste0(main_dir,'data/hg38/repeats/bed_files/UCON31.bed'),
                                  tss_bed = paste0(main_dir,'data/hg38/beds/hg38_proteinCoding_TSS.bed'),
                                  interval_window = 200,min_dist = 5e3,max_dist = 2e6,expand = c(-1e4,1e4),
                                  domains_f = "hic.3DRAM_D45_Tbr2.ins_250_domains_expanded",
                                  score_tracks = score_tracks,hic_names = c('NSCscore','IPCscore'),
                                  other_tracks =c(gtrack.ls('methylation','GpC','10x'),gtrack.ls('methylation','CpG','10x')),other_names = c('NSC_GpC','IPC_GpC','NSC_CpG','IPC_CpG'))

ucon31_p1 <- GOterm_enrichment(distal_mat=ucon31_mat,promoter_bed=paste0(main_dir,'data/hg38/repeats/bed_files/UCON31_prom.bed'),
                        all_TSS=paste0(main_dir,'data/hg38/beds/hg38_proteinCoding_TSS.bed'),
                        closest_TSS = F,orgDB = org.Hs.eg.db)
pdf(paste0('figures/Figure6C_UCON31.pdf'),width=5,height=4,useDingbats = F)
print(dotplot(ucon31_p1$p, showCategory=10,orderBy = "GeneRatio", x='GeneRatio') + theme_cowplot() +theme(legend.position = c(0.85,0.15)))
dev.off()

mer130_mat <- extract_misha_ep_pairs(enhancer_bed = paste0(main_dir,'data/hg38/repeats/bed_files/MER130.bed'),
                                     tss_bed = paste0(main_dir,'data/hg38/beds/hg38_proteinCoding_TSS.bed'),
                                     interval_window = 200,min_dist = 5e3,max_dist = 2e6,expand = c(-1e4,1e4),
                                     domains_f = "hic.3DRAM_D45_Tbr2.ins_250_domains_expanded",
                                     score_tracks = score_tracks,hic_names = c('NSCscore','IPCscore'),
                                     other_tracks =c(gtrack.ls('methylation','GpC','10x'),gtrack.ls('methylation','CpG','10x')),other_names = c('NSC_GpC','IPC_GpC','NSC_CpG','IPC_CpG'))

mer130_p1 <- GOterm_enrichment(distal_mat=mer130_mat,promoter_bed=paste0(main_dir,'data/hg38/repeats/bed_files/MER130_prom.bed'),
                               all_TSS=paste0(main_dir,'data/hg38/beds/hg38_proteinCoding_TSS.bed'),
                               closest_TSS = F,orgDB = org.Hs.eg.db)
p <- dotplot(mer130_p1$p, showCategory=10,orderBy = "GeneRatio", x='GeneRatio') + theme_cowplot()
p <- p + theme(legend.position = c( 0.85,0.25),axis.text.y = element_text(size=18)) + scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 25))
pdf(paste0('figures/Figure6C_MER130.pdf'),width=8,height=8)
print(p)
dev.off()


res = read.table('/home/hpc/bonev/projects/rna/ram/IPCvsNSC_res_proteinCoding_all.tsv')
fpkm_res <- read.table('/home/hpc/bonev/projects/rna/ram/fpkm_annotated_means.txt')
ucon31_mat_res <- ucon31_p1$mat
ucon31_mat_res$log2FoldChange <- res$log2FoldChange[match(ucon31_mat_res$geneName,row.names(res))]
ucon31_mat_res$NSC_FPKM <- fpkm_res$NSC[match(ucon31_mat_res$geneName,row.names(fpkm_res))]
ucon31_mat_res$IPC_FPKM <- fpkm_res$IPC[match(ucon31_mat_res$geneName,row.names(fpkm_res))]

mer130_mat_res <- mer130_p1$mat
mer130_mat_res$log2FoldChange <- res$log2FoldChange[match(mer130_mat_res$geneName,row.names(res))]
mer130_mat_res$NSC_FPKM <- fpkm_res$NSC[match(mer130_mat_res$geneName,row.names(fpkm_res))]
mer130_mat_res$IPC_FPKM <- fpkm_res$IPC[match(mer130_mat_res$geneName,row.names(fpkm_res))]

fc_res <- data.frame(Condition=factor(c(rep('UCON31',length(ucon31_mat_res$log2FoldChange)),rep('MER130',length(mer130_mat_res$log2FoldChange))),levels=c('UCON31','MER130')),Value=c(ucon31_mat_res$log2FoldChange,mer130_mat_res$log2FoldChange))
p1 <- ggplot(fc_res,aes(x=Condition,y=Value,fill=Condition)) + geom_boxplot(outlier.size=1,show.legend = F,width=0.8,outlier.shape = NA) + xlab('') + coord_cartesian(ylim=c(-4.5,4.5))
p1 <- p1 + scale_fill_manual(values=c('gold','cyan')) + ylab('log2 Fold Change') + theme(legend.position = "none")  
pdf('figures/Figure6J.pdf',height=5,width=3.5)
p1 + geom_hline(yintercept = 0,linetype = 2,col='black') 
dev.off()

cols=colorRampPalette(c('blue','white','red'))(101)
features <- c('NEUROD1','PRDM8','BCL11B','ROBO1')
p <- ggplot(mer130_mat_res,aes(x=IPC_GpC_distal-NSC_GpC_distal,y=IPC_CpG_distal-NSC_CpG_distal,fill=log2FoldChange)) + geom_point(pch = I(21),size = 2)  
p <- p + scale_fill_gradientn(colours = cols,name='',breaks=c(-2,-1,0,1,2),labels=c('-2','-1','0','+1','+2'),limits=c(-5,5))
p <- p + xlab('dGpC Accessibility (IPC-NSC)') + ylab('dCpG Methylation (IPC-NSC)')
p <- p + coord_cartesian(xlim=c(-70,70),ylim=c(-70,70)) + geom_vline(xintercept = 0,lty=2) + geom_hline(yintercept = 0,lty=2)
p <- p + geom_text_repel(
  data = mer130_mat_res, size = 2,box.padding = 0.8,segment.alpha = 0.8, min.segment.length = 0,max.iter = 10000,
  aes(col='black',label=ifelse(geneName%in%features, as.character(geneName), "")),force=10)



res <- monalisa_peaks(peaks_f=paste0(main_dir,'data/hg38/repeats/bed_files/UCON31.bed'),pwms = pwms,p.cutoff = 2)
pdf('figures/Figure6D.pdf',width=6,height=6)
p <- plotMotifHeatmaps(x = res$seSel, which.plots = c("log2enr","negLog10Padj"), width = 1,
                  cluster = res$hcl, maxEnr = 3, maxSig = 10,
                  show_dendrogram = TRUE, show_seqlogo = TRUE,
                  width.seqlogo = 1.2)
dev.off()

tracks <- list.files(paste0(main_dir,'data/hg38/RAM_bigwigs/'),pattern = 'NOMe',full.names = T)
peaks <- paste0(main_dir,'data/hg38/repeats/bed_files/UCON31_bed6.bed')
res <- getPlotSetArray(tracks=tracks,features=peaks,refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 40,ignore_strand = F)
pdf('figures/Figure6E.pdf',width=5,height=5)
plotAverage(plotset=res[c(2,4)], labels = c('NSC','IPC'), xlim = NULL,
            ylim = c(3,28), main = '', xlab = "", ylab = '% GpC Accessibility',
            plotScale = "linear", type = "full", error.estimates = T,yaxs='i',
            legend = T, legend_ext = F, legend_pos = "topright",
            legend_ext_pos = "topleft", cex.axis = 12, cex.lab = 12,xaxt='n',
            cex.main = 12, cex.legend = 10, ln.v = FALSE, ln.h = NULL, pointsize = 12, colvec = cell_colors)
axis(side = 1,labels=c('-2kb','UCON31','+2kb'),at=c(-2000,0,2000),pos =3,tick = T)
dev.off()

anno <- read.table(paste0(main_dir,'data/hg38/repeats/bed_files/UCON31.bed'))
anno <- gintervals(anno[,1],anno[,2],anno[,3])
anno <- intervals.normalize(anno,1000)   #Makes all intervals 1kb long

plotMisha(main_f=main_f,targetGene='ZBTB18',outDir='figures/',out_f='ZBTB18_UCON31',upstream=5e4,downstream=2.5e5,chipYlim=matrix(c(0,100,0,100),nrow = 10,ncol = 2,byrow = T),
          chipNames=c('NSC GpC','IPC GpC'),window_scale=1.8,chipRes=20,pointCEX=3,conditions=score_tracks,
          chipTracksToExtract=c("methylation.3DRAM_Pax6_GpC_merged_10x","methylation.3DRAM_Tbr2_GpC_merged_10x"),chipColors=c(cell_colors), #,rep('black',4)
          annIntervals=anno,
          plotOrder=list(scores=TRUE, domains=FALSE, loops=FALSE, VP=FALSE, arcs=FALSE, rna=FALSE, chip=TRUE,anno=TRUE,meth=FALSE,axis=FALSE,scATAC=FALSE, genes=TRUE,ideogram=TRUE),
          plotRatios=list(unitHeight=120, scores=2.5, VP=1.5, loops=2.2, rna=0.6, chip=0.5,meth=0.5, domains=0.15, genes=0.8, arcs=0.5, axis=0.4,ideogram=0.2,scATAC=4,anno=0.2))
dev.off()
