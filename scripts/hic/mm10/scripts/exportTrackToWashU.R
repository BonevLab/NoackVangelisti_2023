genome <- 'mm10'
chrom_sizes_f <- '/home/hpc/bonev/annotations/mm10/mm10.chrom.sizes'
outdir_f <- '/home/hpc/bonev/projects/hic/test/data/extractedBins/washU/'
sample <- 'ES'
trackdb <- paste0("/home/hpc/bonev/trackdb/",genome,"/")

require(misha)
require(plyr)
require(dplyr)
gsetroot(trackdb)
gdb.reload()

chrs <- c(18)
trackmode <- 'avg'
binSize <- 5000
maxSize <- 4e6

tracks <- gtrack.ls("hic.",sample,".score_k200")
v_tracks <- paste0('v_',tracks)


all_intervs <- gintervals.2d.all()
for (i in 1:length(tracks)){
  gvtrack.create(v_tracks[i],tracks[i],func=trackmode)
}

for (chr in chrs){
  iter2d <- giterator.intervals(intervals=gintervals.2d(chr), iterator=c(binSize,binSize), band=-c(maxSize,binSize))
  grid <- gextract(v_tracks,iter2d,iterator = iter2d,colnames='score',band = -c(maxSize,binSize))
  grid <- grid[!is.na(grid$score),-8]
  grid$score <- round(grid$score,digits = 0)
  grid1 <- grid[,1:3]
  grid1$anchor <- do.call(sprintf, c(grid[c('chrom2','start2','end2','score')], '%s:%s-%s,%s'))
  grid2 <- grid[,4:6]
  colnames(grid1) <- c('chrom','start','end','inter')
  grid2$anchor <- do.call(sprintf, c(grid[c('chrom1','start1','end1','score')], '%s:%s-%s,%s'))
  colnames(grid2) <- c('chrom','start','end','inter')
  grid_f <- rbind(grid1,grid2)
  grid_f <- grid_f[order(grid_f$start),]
  grid_f$ID <- 1:nrow(grid_f)
  grid_f$dir <- '.'
  name_f <- paste0(outdir_f,sample,'_',chr,'.txt')
  write.table(grid_f[,1:4],name_f,quote=F,col.names=F,row.names=F,sep='\t')
}


