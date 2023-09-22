source('scripts/config.R')
source('scripts/aux_functions.R')
source('scripts/plot_functions.R')

source(paste0(main_dir,'results/mm10/hic/config.R'))
source(paste0(main_f,'scripts/main_functions.R'))
source(paste0(main_f,'scripts/aux_functions.R'))
source(paste0(main_f,'scripts/plot_functions.R'))
source(paste0(main_f,'scripts/temp_functions.R'))

###Files
smf_files <- list.files(paste0(main_dir,'results/mm10/pair_smf/'),pattern = 'fst',full.names = T)
anchor_files <- list.files(paste0(main_dir,'data/mm10/beds/pairs/'),pattern = 'bed',full.names = T)
bed_files <- list.files(paste0(main_dir,'data/mm10/beds/'),pattern = 'bed',full.names = T)
###

### Read data and filter interv ########
pair_smf_f <- read.fst(smf_files[grep('ES_CTCF_top30k_motifs_GpC',smf_files)])
lanchor <- read.table(bed_files[grep('ES_CTCF_for_top30k_motifs',bed_files)])[,1:3]
ranchor <- read.table(bed_files[grep('ES_CTCF_rev_top30k_motifs',bed_files)])[,1:3]

binSize=100   #window +- interval center
res1 <- filter_smf(df=pair_smf_f,binSize=binSize,filter_dist=c(250,4e6),ncalls=1,                             #only reads with at least this many calls
                  filter_interv_l = intervals.centers(lanchor),ldist = 0,
                  filter_interv_r = intervals.centers(ranchor),rdist = 0,out_res_coord = TRUE) 
res2 <- filter_smf(df=pair_smf_f,binSize=binSize,filter_dist=c(250,4e6),ncalls=1,                             #only reads with at least this many calls
                   filter_interv_l = intervals.centers(ranchor),ldist = 0,
                   filter_interv_r = intervals.centers(lanchor),rdist = 0,out_res_coord = TRUE) 
res <- list(res_ratios=rbind(res1$res_ratios,res2$res_ratios),res_coords=rbind(res1$res_coords,res2$res_coords),
            lcoord_mat=rbindlist(list(res1$lcoord_mat,res2$lcoord_mat),fill=T),rcoord_mat=rbindlist(list(res1$rcoord_mat,res2$rcoord_mat),fill=T),
            lcall_mat=rbindlist(list(res1$lcall_mat,res2$lcall_mat),fill=T),rcall_mat=rbindlist(list(res1$rcall_mat,res2$rcall_mat),fill=T),
            linterv_mat=rbindlist(list(res1$linterv_mat,res2$linterv_mat),fill=T),rinterv_mat=rbindlist(list(res1$rinterv_mat,res2$rinterv_mat),fill=T),idx=unique(c(res1$idx,res2$idx)))

plot_res <- suppressWarnings(plot_paired_smf(res,TF_name='CTCF',met_lims=c(0.2,0.5),
                cluster_names=c('C1','C2','C3','C4'),
                cols = colorRamp2(c(0,100), c("grey","red")),
                mat_name='%GpC accessibility',
                ggplot_cols=hue_pal()(4),
                window=1,
                window2=10,
                binSize=binSize,
                pair_smf_f=pair_smf_f,
                meth_colors=c('blue','red'),
                point.size=16,
                plot_binSize=200))

saveRDS(list(res=res,plot_res=plot_res),paste0(main_dir,'results/mm10/pair_smf/res/CTCF.RDS'))

odds_ratio_mat <- matrix(rev(table(plot_res$cluster_id)),ncol=2,nrow=2,byrow = T)
fisher.test(odds_ratio_mat)


pdf('figures/Figure3B.pdf',height=8,width=4)
draw(plot_res$hm,heatmap_legend_side = "bottom",annotation_legend_side = "bottom")
dev.off()

pdf('figures/Figure3C.pdf',height=3,width=2.5,pointsize = 30,useDingbats = F)
plot_res$p2$p_la|plot_res$p2$p_ra 
dev.off()

peaks <- res$res_coords
peaks <- peaks[!duplicated(peaks$index),]

write.table(cbind(peaks[,1:3],'.',as.numeric(plot_res$cluster_id),'+'),paste0(main_dir,'results/mm10/pair_smf/res/Figure3E_LA.bed'),quote=F,col.names=F,row.names=F,sep='\t')
write.table(cbind(peaks[,4:6],'.',as.numeric(plot_res$cluster_id),'+'),paste0(main_dir,'results/mm10/pair_smf/res/Figure3E_RA.bed'),quote=F,col.names=F,row.names=F,sep='\t')

tracks <- c('/home/hpc/bonev/data/ChIPseq/norm_wigs/ES_CTCF.bw','/home/hpc/bonev/projects/ram/data/mm10/RAM_bigwigs/3DRAM_ES_250k_merged_GpC_cov10x.bw')
seqplots_res1 <- getPlotSetArray(tracks=tracks,features=paste0(main_dir,'results/mm10/pair_smf/res/Figure3E_LA.bed'),refgenome='mm10',type = 'mf',add_heatmap=T,xmin=500,xmax=500,bin = 5,ignore_strand = F)  #Extract CpG and CpG from loop coordinates
seqplots_res2 <- getPlotSetArray(tracks=tracks,features=paste0(main_dir,'results/mm10/pair_smf/res/Figure3E_RA.bed'),refgenome='mm10',type = 'mf',add_heatmap=T,xmin=500,xmax=500,bin = 10,ignore_strand = F)  #Extract CpG and CpG from loop coordinates


pdf('figures/Figure3D.pdf',width=5,height=8)
par(mfrow=c(1,2))
plotHeatmap(seqplots_res1[1],labels = c('CTCF'),sortrows = FALSE,clspace = rev(colorpalette('reds',12)),
            clstmethod='bed_scores',raster = T,embed=T,FO = order(plot_res$cluster_id),CL=plot_res$cluster_id)
plotHeatmap(seqplots_res1[2],labels = c('GpC accessibility'),sortrows = FALSE,clspace = rev(colorpalette('reds',12)),
            clstmethod='bed_scores',raster = T,embed=T)
dev.off()

### Calculate enrichment at clusters

lin_scores <- extract_linear_score(res,plot_res,pair_smf_f,tracks=c("chipseq_RPM.ES_CTCF",'chipseq_RPM.ES_Smc1'),expand=c(-100,100),add_shuffle = TRUE,shift_r=c(seq(-2e6,-5e4,by=10000),seq(5e4,2e6,by=10000)))

score_p <- plot_linear_scores(scores=lin_scores,cluster_cols=c(hue_pal()(4),'grey80'),cluster_names=c('low','mid1','mid2','high'),chip_name='ES_CTCF')
p1 <- score_p + ylim(c(0,25)) + ylab('Average ChIP signal (RPM)') + ggtitle('Ctcf')+ theme(plot.title = element_text(hjust = 0.5))
score_p <- plot_linear_scores(scores=lin_scores,cluster_cols=c(hue_pal()(4),'grey80'),cluster_names=c('low','mid1','mid2','high'),chip_name='ES_Smc1')
p2 <- score_p + ylim(c(0,6)) + ylab('Average ChIP signal (RPM)') + ggtitle('Smc1')+ theme(plot.title = element_text(hjust = 0.5))

pdf('figures/Figure3E.pdf',height=8.5,width=5)
print(p1/p2)
dev.off()

### Enhancer-Promoter
pair_smf_f <- read.fst(smf_files[grep('ES_ATAC_top50_Prom_GpC',smf_files)])
lanchor <- read.table(bed_files[grep('mES_ATAC_top50k_distal',bed_files)])[,1:3]
ranchor <- read.table(bed_files[grep('ES_activeTss',bed_files)])[,1:3]

binSize=100   #window +- interval center
res <- filter_smf(df=pair_smf_f,binSize=binSize,filter_dist=c(5000,4e6),ncalls=1,                             #only reads with at least this many calls
                   filter_interv_l = intervals.centers(lanchor),ldist = 0,
                   filter_interv_r = intervals.centers(ranchor),rdist = 0,out_res_coord = TRUE) 

plot_res <- suppressWarnings(plot_paired_smf(res,TF_name='E-P',met_lims=c(0.2,0.5),
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

odds_ratio_mat <- matrix(rev(table(plot_res$cluster_id)),ncol=2,nrow=2,byrow = T)
fisher.test(odds_ratio_mat)

pdf('figures/Figure3F.pdf',height=8,width=4)
draw(plot_res$hm,heatmap_legend_side = "bottom",annotation_legend_side = "bottom")
dev.off()

pdf('figures/Figure3G.pdf',height=3,width=2.5,pointsize = 30,useDingbats = F)
plot_res$p2$p_la|plot_res$p2$p_ra
dev.off()

peaks <- res$res_coords
peaks <- peaks[!duplicated(peaks$index),]

write.table(cbind(peaks[,1:3],'.',as.numeric(plot_res$cluster_id),'+'),paste0(main_dir,'results/mm10/pair_smf/res/Figure3H_LA.bed'),quote=F,col.names=F,row.names=F,sep='\t')
write.table(cbind(peaks[,4:6],'.',as.numeric(plot_res$cluster_id),'+'),paste0(main_dir,'results/mm10/pair_smf/res/Figure3H_RA.bed'),quote=F,col.names=F,row.names=F,sep='\t')

tracks <- c('/home/hpc/bonev/projects/atac/public/mES_GSE113592.bw','/home/hpc/bonev/projects/ram/data/mm10/RAM_bigwigs/3DRAM_ES_250k_merged_GpC_cov10x.bw')
seqplots_res1 <- getPlotSetArray(tracks=tracks,features=paste0(main_dir,'results/mm10/pair_smf/res/Figure3H_LA.bed'),refgenome='mm10',type = 'mf',add_heatmap=T,xmin=5000,xmax=5000,bin = 50,ignore_strand = F)  #Extract CpG and CpG from loop coordinates
seqplots_res2 <- getPlotSetArray(tracks=tracks,features=paste0(main_dir,'results/mm10/pair_smf/res/Figure3H_RA.bed'),refgenome='mm10',type = 'mf',add_heatmap=T,xmin=5000,xmax=5000,bin = 50,ignore_strand = F)  #Extract CpG and CpG from loop coordinates


pdf('figures/Figure3H_1.pdf',width=5,height=8)
par(mfrow=c(1,2))
plotHeatmap(seqplots_res1[1],labels = c('ATAC'),sortrows = FALSE,clspace = rev(colorpalette('reds',12)),
            clstmethod='bed_scores',raster = T,embed=T,FO = order(plot_res$cluster_id),CL=plot_res$cluster_id)
plotHeatmap(seqplots_res1[2],labels = c('GpC accessibility'),sortrows = FALSE,clspace = rev(colorpalette('reds',12)),
            clstmethod='bed_scores',raster = T,embed=T)
dev.off()


pdf('figures/Figure3H_2.pdf',width=5,height=8)
par(mfrow=c(1,2))
plotHeatmap(seqplots_res2[1],labels = c('ATAC'),sortrows = FALSE,clspace = rev(colorpalette('reds',12)),
            clstmethod='bed_scores',raster = T,embed=T,FO = order(plot_res$cluster_id),CL=plot_res$cluster_id)
plotHeatmap(seqplots_res2[2],labels = c('GpC accessibility'),sortrows = FALSE,clspace = rev(colorpalette('reds',12)),
            clstmethod='bed_scores',raster = T,embed=T)
dev.off()




lin_scores <- extract_linear_score(res,plot_res,pair_smf_f,tracks=c("chipseq_RPM.ES_H3K27ac","atac.ES_GSE113592"),expand=c(-250,250),add_shuffle = TRUE,shift_r=c(seq(-2e6,-5e4,by=10000),seq(5e4,2e6,by=10000)))

score_p <- plot_linear_scores(scores=lin_scores,cluster_cols=c(hue_pal()(4),'grey80'),cluster_names=c('low','mid1','mid2','high'),chip_name='ES_GSE113592')
p1 <- score_p + ylim(c(0,125)) + ylab('Average accessibility (RPM)') + ggtitle('ATAC')+ theme(plot.title = element_text(hjust = 0.5))
score_p <- plot_linear_scores(scores=lin_scores,cluster_cols=c(hue_pal()(4),'grey80'),cluster_names=c('low','mid1','mid2','high'),chip_name='ES_H3K27ac')
p2 <- score_p + ylim(c(0,2.5)) + ylab('Average ChIP signal (RPM)') + ggtitle('H3K27ac')+ theme(plot.title = element_text(hjust = 0.5))

pdf('figures/Figure3I.pdf',height=8.5,width=5)
print(p1/p2) 
dev.off()
