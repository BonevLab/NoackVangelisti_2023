main_f<- '/home/hpc/bonev/projects/ram/results/mm10/hic/'

genome <- 'mm10'
chrom_sizes_f <- '/home/hpc/bonev/annotations/mm10/encode/mm10.chrom.sizes'
map3c_f <- "/home/hpc/bonev/software/hicpipe/" 
cells <- c('ES_3DRAM_250k','ES_Methylhic_Ecker' ,'ES_Methylhic_Ren','ES_Bonev2017')
#cells <- c('E14_NSC','E14_IPC','E14_PN')
#cells <- c('S2_3DRAM_2259','S1_3DRAM_53646','S3_3DRAM_8082','S5_3DRAM_8570','S4_3DRAM_S559')
bowtieInd <- "/home/hpc/bonev/annotations/mm10/encode"
blacklist <-  read.table('/home/hpc/bonev/annotations/mm10/encode/mm10.blacklist.bed.gz')

sunGrid=F
trackdb <- paste0("/home/hpc/bonev/trackdb/",genome,"/")

wide_red_blue_pal = colorRampPalette(c("purple", "navy", "blue", "#87FFFF", "white", "#FF413D", "black", "orange", "yellow"))
reds = colorRampPalette(c("white","orange","red","darkRed"))
cworld_path <- '/home/hpc/bonev/software/cworld/'

require(misha)
gsetroot(trackdb)
gdb.reload()

require(reshape2)
require(plyr)
require(dplyr)
require(zoo)
options(scipen=1e9)
options(gmax.data.size=1.5e8)
options(gmultitasking=F)
require(plotrix)
require(matrixStats)
require(hashmap)
require(pheatmap)
require(tidyr)
source(paste0(main_f,'scripts/aux_functions.R'))
require(ggplot2)
require(magicaxis)
require(gridBase)
require(grid)
require(Gviz)

chrom_sizes <- gintervals.all()
all_tracks <- as.vector(unlist(sapply(paste0('hic.',cells),function(x){gtrack.ls(x)},simplify=T)))
all_tracks <- all_tracks[grep('shuffle',all_tracks,invert=T)]
score_tracks <- all_tracks[grep('score',all_tracks)]
all_tracks <- all_tracks[grep('score',all_tracks,invert=T)]
all_tracks <- all_tracks[grep('ins',all_tracks,invert=T)]
all_tracks <- all_tracks[grep('others',all_tracks,invert=T)]
all_tracks <- all_tracks[grep('HiChIP',all_tracks,invert=T)]
all_tracks <- all_tracks[grep('100kb',all_tracks,invert=T)]
all_tracks <- all_tracks[grep('INS',all_tracks,invert=T)]

eigen_tracks=gtrack.ls('eigen')
eigen_tracks <- eigen_tracks[grep('rep',eigen_tracks,invert=T)]
chip_tracks <- gtrack.ls('chipseq_RPM.')
repli_tracks <- gtrack.ls('repli.')
ins_tracks <- gtrack.ls('hic.','ins_250')
ins_rep_tracks <- gtrack.ls('hic.','ins_250')[grep('rep',gtrack.ls('hic.','ins_250'))]
ins_tracks <- ins_tracks[grep('rep',ins_tracks,invert=T)]
meth_tracks <- gtrack.ls('methylation.','10x')
atac_tracks <- gtrack.ls('scATAC')

tss_f <- "glp.intervals.ucscCanTSS"
tss_full <- "glp.intervals.genCodeUniqueTSS"
genes_f <- "glp.intervals.ucscCanGenes"
blacklist <- gintervals(chroms=blacklist$V1,starts=blacklist$V2,ends=blacklist$V3)

dir.create(paste0(main_f,'logs/'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'qsub/'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'data/'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'temp/'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'analysis/'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'analysis/figures/'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'analysis/figures/gene_plots/'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'analysis/aggregateHiC'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'analysis/aggregateDiagonal'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'analysis/averageTAD'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'analysis/cis_decay'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'analysis/compartments'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'analysis/insulation'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'data/aggregateHiC'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'data/aggregateDiagonal'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'data/cis_decay'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'data/diffHiC'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'data/extractedBins'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'data/fastq'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'logs/aggregateHiC'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'logs/averageTAD'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'logs/aggregateDiagonal'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'logs/cis_decay'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'logs/map3c'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'logs/diffHiC'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'logs/prepareTracks'), showWarnings = FALSE,recursive=T)
dir.create(paste0(main_f,'data/extractedGenes/'), showWarnings = FALSE,recursive=T)
