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
mcparams <- BiocParallel::MulticoreParam(10L)


library(JASPAR2020)
library(Matrix)
library(TFBSTools)
library(motifmatchr)
library(Hmisc)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
library(monaLisa)
library(clusterProfiler)
library(org.Hs.eg.db)

TE_beds <- paste0(main_dir,'data/hg38/repeats/bed_files/hg38_repeats_cleaned.bed')
peaks <- read.table(TE_beds)
colnames(peaks) <- c('chrom','start','end','Class','TE_class','strand')
hic_names=c('NSCscore','IPCscore')
other_tracks <- c("methylation.3DRAM_Pax6_CpG_merged_10x","methylation.3DRAM_Tbr2_CpG_merged_10x","methylation.3DRAM_Pax6_GpC_merged_10x","methylation.3DRAM_Tbr2_GpC_merged_10x","hic.3DRAM_D45_Pax6.ins_250","hic.3DRAM_D45_Tbr2.ins_250")
other_names <- c("NSC_CpG","IPC_CpG","NSC_GpC","IPC_GpC","NSC_ins","IPC_ins")

monalisa_peaks <- function(peaks_f,pwms,p.cutoff,log2enr.cutoff=0,mcparams=BiocParallel::MulticoreParam(1L)){
  peaks <- read.table(peaks_f)
  colnames(peaks) <- c('chrom','start','end','family','subfamily','strand')
  peaks <- makeGRangesFromDataFrame(peaks)
  peakSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, peaks)
  
  
  se <- calcBinnedMotifEnrR(seqs = peakSeq,pwmL = pwms,background = 'genome',genome = BSgenome.Hsapiens.UCSC.hg38,genome.oversample = 10,BPPARAM = mcparams)
  sel1 <- apply(assay(se, "negLog10Padj"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > p.cutoff
  sel2 <- apply(assay(se, "log2enr"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > log2enr.cutoff
  seSel <- se[sel1&sel2, ]
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

res <- rbind(rbind(res,res_ctcf),res_shuffle)
res$TE_family <- peaks$Class[match(row.names(res),peaks$TE_class)]

#Figures

df <- plot_repeat_scatter(res,metric=c('intraHiC'),min_n = 100,diff.cutoff = 5,anno.size = 3,cols = c('grey',cell_colors),
                          features=c('CTCF_convergent','random','LTR24C','LTR55','Charlie9','Charlie4','MER21-int','PABL_B-int','LTR66','Charlie17','ERVL47-int','MER68-int','L1M3f','X9b_DNA','HUERS-P1-int'))
pdf('figures/FigureS6A.pdf',height=6,width=6,useDingbats = F)
print(df$p+xlab('RGC'))
dev.off()

df <- plot_repeat_scatter(res,metric=c('interHiC'),min_n = 100,diff.cutoff = 5,anno.size = 3,cols = c('grey',cell_colors),
                          features=c('CTCF_convergent','random','HERVE-int','LTR44','MER70-int','LTR34','HUERS-P3-int','LTR19A'))
pdf('figures/FigureS6B.pdf',height=6,width=6,useDingbats = F)
print(df$p+xlab('RGC'))
dev.off()


df <- plot_repeat_scatter(res,metric=c('ins'),min_n = 50,diff.cutoff = 5,anno.size = 3,cols = c('grey',cell_colors),
                          features=c('CTCF_convergent','random','HERVE-int','HERVH-int','LTR10G','LTR43-int','LTR108a_Mam','LTR13'))
pdf('figures/FigureS6C.pdf',height=6,width=6,useDingbats = F)
print(df$p+ geom_smooth(method = "lm", se = TRUE))
dev.off()

 df <- plot_repeat_scatter(res,metric=c('CpG'),min_n = 30,diff.cutoff = 5,anno.size = 4,cols = c('grey',cell_colors),
                           features=c('MER130','UCON31','random'))
 pdf('figures/FigureS6D1.pdf',height=5,width=5,useDingbats = F)
 print(df$p + geom_smooth(method = "lm", se = TRUE)+xlab('RGC')+coord_cartesian(xlim=c(63,100),ylim = c(63,100)))
 dev.off()



pwms <- getMatrixSet(JASPAR2020,
                     opts = list(matrixtype = "PWM",
                                 species = 9606))
  

res <- monalisa_peaks(peaks_f=paste0(main_dir,'data/hg38/repeats/bed_files/LTR24C.bed'),pwms = pwms,p.cutoff = 8,log2enr.cutoff = 1)
pdf('figures/FigureS6F.pdf',width=8,height=8)
plotMotifHeatmaps(x = res$seSel, which.plots = c("log2enr","negLog10Padj"), width = 1.8,
                  cluster = res$hcl,
                  show_dendrogram = TRUE, show_seqlogo = TRUE,
                  width.seqlogo = 1.2)
dev.off()

res <- monalisa_peaks(peaks_f=paste0(main_dir,'data/hg38/repeats/bed_files/HERVE.bed'),pwms = pwms,p.cutoff = 8,log2enr.cutoff = 1,mcparams = mcparams)
pdf('figures/FigureS6G.pdf',width=8,height=8)
plotMotifHeatmaps(x = res$seSel, which.plots = c("log2enr","negLog10Padj"), width = 1.8,
                  cluster = res$hcl,
                  show_dendrogram = TRUE, show_seqlogo = TRUE,
                  width.seqlogo = 1.2)
dev.off()

res <- monalisa_peaks(peaks_f=paste0(main_dir,'data/hg38/repeats/bed_files/LTR43.bed'),pwms = pwms,p.cutoff = 8,log2enr.cutoff = 1)
pdf('figures/FigureS6H.pdf',width=8,height=8)
plotMotifHeatmaps(x = res$seSel, which.plots = c("log2enr","negLog10Padj"), width = 1.8,
                  cluster = res$hcl,
                  show_dendrogram = TRUE, show_seqlogo = TRUE,
                  width.seqlogo = 1.2)
dev.off()

mat <- extract_misha_ep_pairs(enhancer_bed = paste0(main_dir,'data/hg38/repeats/bed_files/MER130.bed'),
                               tss_bed = paste0(main_dir,'data/hg38/beds/hg38_proteinCoding_TSS.bed'),
                               interval_window = 200,min_dist = 5e3,max_dist = 2e6,expand = c(-1e4,1e4),
                               domains_f = "hic.3DRAM_D45_Tbr2.ins_250_domains_expanded",
                               score_tracks = score_tracks,hic_names = c('NSCscore','IPCscore'),
                               c(gtrack.ls('methylation','GpC','10x'),gtrack.ls('methylation','CpG','10x')),other_names = c('NSC_GpC','IPC_GpC','NSC_CpG','IPC_CpG'))


peaks <- read.table(paste0(main_dir,'data/hg38/repeats/bed_files/MER130.bed'))[,1:3]
colnames(peaks) <- colnames(gintervals.all())
tracks=c("methylation.3DRAM_Pax6_GpC_merged_10x","methylation.3DRAM_Tbr2_GpC_merged_10x","methylation.3DRAM_Pax6_CpG_merged_10x","methylation.3DRAM_Tbr2_CpG_merged_10x")
df <- misha_extract(tracks = tracks,regions=peaks,iterator = peaks,mode = 'avg')

df <- df[,4:8]
colnames(df)[1:4] <- c('NSC_GpC','IPC_GpC','NSC_CpG','IPC_CpG')
df_melt <- melt(df,id.vars = 'intervalID')
df_melt$condition <- factor(gsub("_.*","",df_melt$variable),levels=c('NSC','IPC'))
df_melt$mark <- factor(gsub(".*_","",gsub('_distal','',df_melt$variable)),levels=c('GpC','CpG'))

p1 <- ggplot(df_melt,aes(x=mark,y=value,fill=condition)) + geom_boxplot(outlier.size=1,outlier.shape=NA,show.legend = T,width=0.8) 
p1 <- p1 + scale_fill_manual(values=cell_colors,name='') + xlab('') + ylab('% Methylation') + theme(legend.direction = "horizontal",legend.justification=c(1,1), legend.position=c(0.35, 0.9))
p1 <- p1 + stat_compare_means(label = "p.format",method='wilcox',paired = F)
wilcox.test(df$NSC_GpC,df$IPC_GpC,paired = T)$p.value

pdf(paste0('figures/FigureS6D.pdf'),width=5,height=5)
print(p1)
dev.off()


p2 <- GOterm_enrichment(distal_mat=mat,promoter_bed=paste0(main_dir,'data/hg38/repeats/bed_files/MER130_prom.bed'),
                        all_TSS=paste0(main_dir,'data/hg38/beds/hg38_proteinCoding_TSS.bed'),
                        closest_TSS = F)
pdf(paste0('figures/Figure6E.pdf'),width=6,height=4)
print(dotplot(p2$p, showCategory=10,orderBy = "GeneRatio", x='GeneRatio'))
dev.off()
