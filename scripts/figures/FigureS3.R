source('scripts/config.R')
source('scripts/aux_functions.R')
source('scripts/plot_functions.R')

source(paste0(main_dir,'results/mm10/hic/config.R'))
source(paste0(main_f,'scripts/main_functions.R'))
source(paste0(main_f,'scripts/aux_functions.R'))
source(paste0(main_f,'scripts/plot_functions.R'))
source(paste0(main_f,'scripts/temp_functions.R'))


###Files

smf_files <- list.files(paste0(main_dir,'results/mm10/pair_smf/controls/'),pattern = 'fst',full.names = T)
bed_files <- list.files(paste0(main_dir,'data/mm10/beds/'),pattern = 'bed',full.names = T)

### CTCF shuffled ####
pair_smf_f <- read.fst(smf_files[grep('GpC_CTCF_ShuffleR2_R',smf_files)])
lanchor <- read.table(paste0(main_dir,'results/mm10/pair_smf/res/Figure3E_LA.bed'))[,1:3]

res <- filter_smf(df=pair_smf_f, binSize=100, filter_dist=c(1000, 4e6), ncalls=1,
                  filter_interv_l = intervals.centers(lanchor),ldist = 0,
                  shufflebS_lim=5e3, shuffleSet="R2", 
                  out_res_coord = TRUE)
plot_res <- suppressWarnings(plot_smf_shuffle(res,TF_name='CTCF',met_lims=c(0.2,0.5),k=5,
                                              cluster_names=c('C1','C2','C3','C4','C5'),
                                              cols = colorRamp2(c(0,100), c("grey","red")),
                                              mat_name='%GpC accessibility',
                                              ggplot_cols=hue_pal()(5)))

saveRDS(list(res=res,plot_res=plot_res),paste0(main_dir,'results/mm10/pair_smf/res/CTCF_shuffleR2.RDS'))

odds_ratio_mat <- matrix(rev(table(plot_res$cluster_id)),ncol=2,nrow=2,byrow = T)
fisher.test(odds_ratio_mat)

pdf('figures/FigureS3A.pdf',height=8,width=4)
draw(plot_res$hm,heatmap_legend_side = "bottom",annotation_legend_side = "bottom")
dev.off()

###CTCF convergent overlapping with microC loops

pair_smf_f <- read.fst(smf_files[grep('ES_CTCF_top30k_motifs_GpC',smf_files)])
lanchor <- read.table(bed_files[grep('ES_CTCF_for_top30k_motifs',bed_files)])[,1:3]
ranchor <- read.table(bed_files[grep('ES_CTCF_rev_top30k_motifs',bed_files)])[,1:3]
loops <- read.table(paste0(main_dir,'data/mm10/microc/res2500/ES_loops.tsv'))
loops <- intervals.2d.centers(rbind(gintervals.2d(loops[,1],loops[,2],loops[,3],loops[,4],loops[,5],loops[,6]),
                                    gintervals.2d(loops[,4],loops[,5],loops[,6],loops[,1],loops[,2],loops[,3])))
binSize=100   #window +- interval center
res1 <- filter_smf(df=pair_smf_f,binSize=binSize,filter_dist=c(1000,4e6),ncalls=1,                             #only reads with at least this many calls
                   filter_interv_l = intervals.centers(lanchor),ldist = 0,
                   filter_interv_r = intervals.centers(ranchor),rdist = 0,filter_interv_2d = loops,dist=25000,
                   out_res_coord = TRUE) 
res2 <- filter_smf(df=pair_smf_f,binSize=binSize,filter_dist=c(1000,4e6),ncalls=1,                             #only reads with at least this many calls
                   filter_interv_l = intervals.centers(ranchor),ldist = 0,
                   filter_interv_r = intervals.centers(lanchor),rdist = 0,filter_interv_2d = loops,dist=25000,
                   out_res_coord = TRUE) 
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

pdf('figures/FigureS3B.pdf',height=8,width=4)
draw(plot_res$hm,heatmap_legend_side = "bottom",annotation_legend_side = "bottom")
dev.off()

###Files
smf_files <- list.files(paste0(main_dir,'results/mm10/pair_smf/'),pattern = 'fst',full.names = T)
anchor_files <- list.files(paste0(main_dir,'data/mm10/beds/pairs/'),pattern = 'bed',full.names = T)
bed_files <- list.files(paste0(main_dir,'data/mm10/beds/'),pattern = 'bed',full.names = T)
###

pair_smf_f <- read.fst(smf_files[grep('ES_CTCF_top30k_motifs_GpC',smf_files)])
lanchor <- read.table(bed_files[grep('ES_CTCF_for_top30k_motifs',bed_files)])[,1:3]
ranchor <- read.table(bed_files[grep('ES_CTCF_rev_top30k_motifs',bed_files)])[,1:3]

### CTCF for-for
binSize=100
res <- filter_smf(df=pair_smf_f,binSize=binSize,filter_dist=c(1000,4e6),ncalls=1,               
                  filter_interv_l = intervals.centers(lanchor),ldist = 0,
                  filter_interv_r = intervals.centers(lanchor),rdist = 0,
                  out_res_coord = TRUE) 

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

pdf('figures/FigureS3C1.pdf',height=8,width=4)
draw(plot_res$hm,heatmap_legend_side = "bottom",annotation_legend_side = "bottom")
dev.off()


### CTCF rev-rev

res <- filter_smf(df=pair_smf_f,binSize=binSize,filter_dist=c(1000,4e6),ncalls=1,                             #only reads with at least this many calls
                  filter_interv_l = intervals.centers(ranchor),ldist = 0,
                  filter_interv_r = intervals.centers(ranchor),rdist = 0,
                  out_res_coord = TRUE) 

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

pdf('figures/FigureS3C2.pdf',height=8,width=4)
draw(plot_res$hm,heatmap_legend_side = "bottom",annotation_legend_side = "bottom")
dev.off()

### CTCF GpC vs CpG #####
smf_files <- list.files(paste0(main_dir,'results/mm10/pair_smf/controls/'),pattern = 'fst',full.names = T)
pair_smf_f <- read.fst(smf_files[grep('CpG_GpC_CTCF_500R1.',smf_files)])
tad_borders <- gintervals.load("hic.ES_Bonev2017.ins_250_borders")

res <- filter_smf(df=pair_smf_f, binSize=100, filter_dist=c(0,300), ncalls=1, 
                            out_res_coord = TRUE) 

plot_res <- suppressWarnings(plot_paired_smf(res,TF_name='CTCF',met_lims=c(0.2,0.5),
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


pdf('figures/FigureS3D.pdf',height=8,width=4)
draw(plot_res$hm,heatmap_legend_side = "bottom",annotation_legend_side = "bottom")
dev.off()

pdf('figures/FigureS3E.pdf',height=3,width=2.5,pointsize = 30,useDingbats = F)
plot_res$p2$p_la|plot_res$p2$p_ra 
#(plot_res1$p2$p_la|plot_res1$p2$p_ra)/(plot_res2$p2$p_la|plot_res2$p2$p_ra)
dev.off()

lin_scores <- extract_linear_score(res,plot_res,pair_smf_f,tracks=c("chipseq_RPM.ES_CTCF",'chipseq_RPM.ES_Smc1'),expand=c(-100,100),add_shuffle = TRUE,shift_r=c(seq(-2e6,-5e4,by=10000),seq(5e4,2e6,by=10000)))
df1 <- lin_scores[,grep('\\.y',colnames(lin_scores),invert=T)]
colnames(df1) <- gsub('\\.x','',colnames(df1))

ctcf_coords <- res$res_coords[,1:3]
ctcf_coords <- read.table('../../data/mm10/beds/ES_CTCF_top30k_motifs.bed')[,1:3]
colnames(ctcf_coords) <- c('chrom','start','end')
ctcfs_at_TADbs <- gintervals.neighbors(ctcf_coords[!duplicated(ctcf_coords),],tad_borders,maxdist = 10000)[,1:3]
df2 <- misha_extract(tracks = c("chipseq_RPM.ES_CTCF",'chipseq_RPM.ES_Smc1'),regions = intervals.normalize(ctcfs_at_TADbs,200),iterator = intervals.normalize(ctcfs_at_TADbs,200),mode = 'avg')
df2$index <- 1:nrow(df2)
df2$cluster <- 'TAD Boundary'
df2 <- df2[,match(colnames(df1),colnames(df2))]
df <- rbind(df1,df2)
df$cluster <- factor(df$cluster,levels=c('C1','C2','C3','C4','TAD Boundary','control'))

score_p <- plot_linear_scores(scores=df,cluster_cols=c(hue_pal()(5),'grey80'),names = c('Ctcf'),chip_name='ES_CTCF')
p1 <- score_p + ylab('Average ChIP signal (RPM)') + scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))  
score_p <- plot_linear_scores(scores=df,cluster_cols=c(hue_pal()(5),'grey80'),names = c('Smc1'),chip_name='ES_Smc1')
p2 <- score_p + ylim(c(0,6)) + ylab('Average ChIP signal (RPM)') + scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))

pdf('figures/FigureS3F.pdf',height=8.5,width=4.5)
print(p1/p2) 
dev.off()



res <- filter_smf(df=pair_smf_f, binSize=100, filter_dist=c(0,300), ncalls=1,filter_interv_l = tad_borders,ldist = 5000, 
                  out_res_coord = TRUE) 

plot_res <- suppressWarnings(plot_paired_smf(res,TF_name='CTCF',met_lims=c(0.2,0.5),
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


pdf('figures/FigureS3G.pdf',height=8,width=4)
draw(plot_res$hm,heatmap_legend_side = "bottom",annotation_legend_side = "bottom")
dev.off()

