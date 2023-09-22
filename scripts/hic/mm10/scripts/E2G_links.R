require(misha)
require(tidyr)
require(plyr)
require(dplyr)

gsetroot('/home/hpc/bonev/trackdb/mm10')

construct.grid = function(interv1,interv2,min_dist,max_dist){
  return(ddply(interv1, .(chrom), function(i1) {
    i2 = interv2[as.character(interv2$chrom) == as.character(i1$chrom[1]),]
    if (nrow(i2) ==0) {
      return(c())
    }
    g = expand.grid(1:nrow(i1), 1:nrow(i2))
    g = g[abs(i2$start[g$Var2]-i1$start[g$Var1]) > min_dist & abs(i2$start[g$Var2]-i1$start[g$Var1]) < max_dist,]
    grid = cbind(i1[g$Var1,c("chrom", "start", "end")], i2[g$Var2,c("chrom", "start", "end")])
    colnames(grid) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
    return(grid)
  })[,-1] )
}

min_dist <- 5000
max_dist <- 2e6
expand=c(-2500,2500)
domains <- gintervals.load("hic.ES.ins_250_domains_expanded")
tracks <- c('hic.ES.score_k200','hic.NPC.score_k200')

options(gmax.data.size=1e8)
options(gmultitasking=TRUE)
for (i in 1:length(tracks)){
  gvtrack.create(paste0('v_',tracks[i]),tracks[i],'max')
  gvtrack.iterator.2d(paste0('v_',tracks[i]), sshift1=min(expand), eshift1=max(expand), sshift2=min(expand), eshift2=max(expand))
}


tss <- gintervals.load("glp.intervals.ucscCanTSS")
enh <- read.table('/home/hpc/bonev/data/hic/data/enhancers/ES_all_enhancers.bed')

enh <- gintervals(enh[,1],enh[,2],enh[,3])

grid1 <- construct.grid(tss,enh,min_dist,max_dist)
grid2 <- construct.grid(enh,tss,min_dist,max_dist)
grid <- unique(rbind(grid1,grid2))



dom = gintervals.2d(chroms1=domains$chrom, starts1=domains$start, ends1=domains$end,chroms2=domains$chrom, starts2=domains$start, ends2=domains$end)
intra_grid <- gintervals.intersect(grid,dom)
grid_indx <- as.vector(unite(grid[,c(1,2,5)],col='peaks_ID',sep='_'))[,1]		
intra_indx = as.vector(unite(intra_grid[,c(1,2,5)],col='peaks_ID',sep='_'))[,1]		
inter_grid <- grid[!(grid_indx%in%intra_indx),-7]

intra_scores<- gextract(paste0('v_',tracks),intervals = intra_grid,iterator = intra_grid,band = -c(max_dist+max(expand),min_dist-min(expand)))
inter_scores <- gextract(paste0('v_',tracks),intervals = inter_grid,iterator = inter_grid,band = -c(max_dist+max(expand),min_dist-min(expand)))


intra_scores$domain <- 'intraTAD'
inter_scores$domain <- 'interTAD'
grid <- rbind(intra_scores,inter_scores)
grid <- grid[!is.na(grid$peakID),]
l_anchor <- grid[,1:3]
r_anchor <- grid[,4:6]
colnames(l_anchor) <- c('chrom','start','end')
colnames(r_anchor) <- c('chrom','start','end')
l_dist <- gintervals.neighbors(l_anchor,tss)
r_dist <- gintervals.neighbors(r_anchor,tss)



linksSig <- linksSig[match(grid$peakID,linksSig$IDs)]
linksSig$ESscore <- grid$v_hic.ES.score_k200
linksSig$NPCscore <- grid$v_hic.ncx_Hes5.score_k200
linksSig$CNscore <- grid$v_hic.ncx_Dcx.score_k200
linksSig$domain <- grid$domain
linksSig$deltaHiC <- rowMaxs(as.matrix(mcols(linksSig[,c('ESscore','NPCscore','CNscore')])))-rowMins(as.matrix(mcols(linksSig[,c('ESscore','NPCscore','CNscore')])))
linksSig$deltaNPC_CN_HiC <- rowMaxs(as.matrix(mcols(linksSig[,c('NPCscore','CNscore')])))-rowMins(as.matrix(mcols(linksSig[,c('NPCscore','CNscore')])))
linksSig <- linksSig[order(linksSig$FDR,linksSig$gene_name,100-abs(linksSig$deltaNPC_CN_HiC),decreasing=F),]
