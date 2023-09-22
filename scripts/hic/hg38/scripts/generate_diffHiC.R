#!/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(diffHic, quietly = TRUE, verbose = FALSE)
library(GenomicRanges, quietly = TRUE, verbose = FALSE)
require(edgeR)

cell=args[1]
track=args[2]
message('Working on track: ',track)
n_track=as.numeric(args[3])
chr=args[4]
path=args[5]
mode=args[6]
filter_v=args[7]
main_f=args[8]

source(paste0(main_f,'config.R'))
options(gmax.data.size=2e08)

chrom_sizes <- read.table(paste0(trackdb,'/chrom_sizes.txt'))
chrom_sizes$V1 <- paste0('chr',chrom_sizes$V1)
chromInf <- Seqinfo(seqnames=as.character(chrom_sizes$V1), seqlengths=chrom_sizes$V2,genome = genome)

if (!file.exists(paste0(path,'/frags.G'))){
	frags <- read.table(paste0(trackdb,'/seq/redb/GATC.frags'))
	colnames(frags) <- c('idx','chrom','start','end')
	frags$chrom <- paste0('chr',frags$chrom)
	frags.G <- makeGRangesFromDataFrame(frags, seqnames.field="chrom", start.field="start", end.field="end", ignore.strand=T, seqinfo=chromInf)
	save(frags.G,file=paste0(path,'/frags.G'))
} else {
	frags.G <- get(load(paste0(path,'/frags.G')))
	frags.param<-pairParam(frags.G)
}


path <- paste0(path,cell,'/')
dir.create(path, showWarnings = FALSE,recursive=T)
dir.create(paste0(path,cell,'_',n_track,'/'), showWarnings = FALSE,recursive=T)


file_f <- paste0(path,cell,'_',n_track,'/',as.character(chr),'.h5')
df <- gextract(track,gintervals.2d(chr),band=(-c(gintervals(chr)$end,1000)))
df1 <- makeGRangesFromDataFrame(df, seqnames.field="chrom1", start.field="start1", end.field="end1", ignore.strand=T, seqinfo=chromInf)
df2 <- makeGRangesFromDataFrame(df, seqnames.field="chrom2", start.field="start2", end.field="end2", ignore.strand=T, seqinfo=chromInf)
res1 <- data.frame(anchor.id=nearest(df1,frags.param$fragments),target.id=nearest(df2,frags.param$fragments),count=df[,7])
res2 <- res1[res1$count==2,]
res3 <- res1[res1$count==3,]
res4 <- res1[res1$count==4,]
res5 <- res1[res1$count==5,]
res <- rbind(res1,res2,res3,res3,res4,res4,res4,res5,res5,res5,res5)
res <- res[,-3]
res <- res[order(res$anchor.id,res$target.id),]
colnames(res) <- c('anchor1.id','anchor2.id')
savePairs(res,file=file_f,frags.param)
if (mode=='exp'){
	file_f <- paste0(path,cell,'_',n_track,'/',as.character(chr),'_shuffle.h5')
	df <- gextract(paste0(track,'_shuffle'),gintervals.2d(chr),band=(-c(gintervals(chr)$end,1000)))
	df1 <- makeGRangesFromDataFrame(df, seqnames.field="chrom1", start.field="start1", end.field="end1", ignore.strand=T, seqinfo=chromInf)
	df2 <- makeGRangesFromDataFrame(df, seqnames.field="chrom2", start.field="start2", end.field="end2", ignore.strand=T, seqinfo=chromInf)
	res1 <- data.frame(anchor.id=nearest(df1,frags.param$fragments),target.id=nearest(df2,frags.param$fragments),count=df[,7])
	res2 <- res1[res1$count==2,]
	res3 <- res1[res1$count==3,]
	res4 <- res1[res1$count==4,]
	res5 <- res1[res1$count==5,]
	res <- rbind(res1,res2,res3,res3,res4,res4,res4,res5,res5,res5,res5)
	res <- res[,-3]
	res <- res[order(res$anchor.id,res$target.id),]
	colnames(res) <- c('anchor1.id','anchor2.id')
	savePairs(res,file=file_f,frags.param)
} else if (mode=='norm'){
	score_track <- gtrack.ls(cell,'score')
	df <- gscreen(paste0(score_track,'>=',filter_v),gintervals.2d(chr),iterator=track,band=(-c(gintervals(chr)$end,1000)))
	file_f <- paste0(path,cell,'_',n_track,'/',as.character(chr),'_score',filter_v,'.h5')
	df1 <- makeGRangesFromDataFrame(df, seqnames.field="chrom1", start.field="start1", end.field="end1", ignore.strand=T, seqinfo=chromInf)
	df2 <- makeGRangesFromDataFrame(df, seqnames.field="chrom2", start.field="start2", end.field="end2", ignore.strand=T, seqinfo=chromInf)
	res <- data.frame(anchor.id=nearest(df1,frags.param$fragments),target.id=nearest(df2,frags.param$fragments))
	res <- res[order(res$anchor.id,res$target.id),]
	colnames(res) <- c('anchor1.id','anchor2.id')
	savePairs(res,file=file_f,frags.param)
}