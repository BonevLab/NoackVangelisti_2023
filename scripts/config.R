### Libraries 

library(plyr)
library(dplyr)
library(tidyr)
library(matrixStats)
library(fst)
library(misha)
library(data.table)
library(circlize)
library(ComplexHeatmap)
library(ggplot2)
library(cowplot)
library(patchwork)
library(scales)
library(GenomicRanges)
library(ggpubr)
library(LSD)
library(pals)
library(Hmisc)
library(doParallel)
library(LSD)
library(ggrepel)
library(seqplots)
library(SummarizedExperiment)
library(factoextra)
library(ggrastr)
library(Sushi)
library(ggsci)

theme_set(theme_cowplot())
options(gmultitasking=FALSE)
options(scipen=1e9)


### Global Directories ###
main_dir <- '/home/hpc/bonev/projects/ram/'
mm10_trackdb <- '/home/hpc/bonev/trackdb/mm10/' 
mm10_bl <- paste0(main_dir,'data/mm10/beds/mm10_blacklist.bed')
hg38_trackdb <- '/home/hpc/bonev/trackdb/hg38/'

mm10_hic_f <- ''
hg38_hic_f <- ''

### Colors
archr_colors <- ArchR::ArchRPalettes

cpg_color <- 'red'
gpc_color <- 'blue'
gpc_hm_colors <- rev(colorpalette('reds',10))
cpg_hm_colors <- rev(colorpalette('blues',10))
rep_colors <- glasbey(3)
heatmap_colors <- colorpalette('matlablike',10)

blue_white_pal = colorRampPalette(c("purple", "navy", "blue", "#87FFFF", "white"))
white_red_pal = colorRampPalette(c("white","#FF413D", "black", "orange", "yellow"))

### Genomic databases ###

#txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
