#!/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(diffHic, quietly = TRUE, verbose = FALSE)
library(GenomicRanges, quietly = TRUE, verbose = FALSE)
require(edgeR)
require(csaw)
require(statmod)

cells=args[1]
cells=as.vector(unlist(strsplit(cells,',')))
chr=as.character(args[2])
binSize=as.numeric(args[3])
path=args[4]
read_thresh=as.numeric(args[5])
shuffle_thresh=as.numeric(args[6])
shuffle_keep=as.logical(args[7])
name=args[8]
main_f=args[9]

# cells <- c('D0','D2','D4','D6','D10')
# chr='chr19'
# path=paste0(main_f,'analysis/diffHiC/',ifelse(shuff_norm,'shaman/','loess/'))
path <- paste0(path,name,'/',binSize/1000,'kb/')
dir.create(path, showWarnings = FALSE,recursive=T)

expand.grid.unique <- function(x, y, include.equals=FALSE)
{
    x <- unique(x)
    y <- unique(y)
    g <- function(i)
    {
        z <- setdiff(y, x[seq_len(i-include.equals)])
        if(length(z)) cbind(x[i], z, deparse.level=0)
    }
    do.call(rbind, lapply(seq_along(x), g))
}



frags.G <- get(load(paste0(main_f,'data/diffHiC/frags.G')))
frags.param <- pairParam(frags.G)
param <- reform(frags.param, restrict=chr)
rep_h5 <- list.files(path=paste0(main_f,'data/diffHiC/'),pattern=paste0(chr,'.h5'),full.names=T,recursive=T,include.dirs=T)
rep_h5 <- as.vector(unlist(sapply(cells,function(x){return(rep_h5[grep(paste0('/',x),rep_h5,fixed=T)])})))
shuffle_rep_h5 <- paste0(gsub('.h5','',rep_h5),'_shuffle.h5')
cells_idx <- as.vector(unlist(sapply(cells,function(x){return(rep(x,length(grep(paste0('/',x),rep_h5,fixed=T))))})))
reps_idx <- as.vector(unlist(sapply(cells,function(x){return(rep(1:length(grep(paste0('/',x),rep_h5,fixed=T))))})))
design_matrix <- data.frame(cells=cells_idx,reps=reps_idx,files=rep_h5)
print(design_matrix)
design <- model.matrix(~0+factor(cells_idx,levels=cells))
colnames(design) <- cells

print(paste("reading files",rep_h5))
#data <- squareCounts(c(rep_h5,shuffle_rep_h5), param, width=binSize, filter=read_thresh)
#shuffle_data <- data[,(length(rep_h5)+1):ncol(data)]
#data <- data[,1:length(rep_h5)]

#smaller.data <- squareCounts(c(rep_h5,shuffle_rep_h5), param, width=binSize/2, filter=read_thresh)
#smaller.shuffle_data <- smaller.data[,(length(rep_h5)+1):ncol(smaller.data)]
#smaller.data <- smaller.data[,1:length(rep_h5)]

#data_trended <- squareCounts(rep_h5, param, width=5e5, filter=read_thresh)


#### Testng with preset regions #######

p2glinks <- readRDS('/home/hpc/bonev/projects/SC/E14_analysis/results/P2G-Links.RDS')
p2glinks <- c(p2glinks$posCor,p2glinks$negCor,p2glinks$noCor)
DA_anchor <- GRanges(p2glinks)
mcols(DA_anchor) <- NULL
DA_anchor <- resize(x = DA_anchor,width = binSize+1,fix = 'center')
GA_anchor <- GRanges(seqnames = p2glinks$gene_chr,IRanges(start=p2glinks$gene_start,end = p2glinks$gene_start+1,names=p2glinks$gene_name))
GA_anchor <- resize(x = GA_anchor,width = binSize+1,fix = 'center')
print(paste("reading files",rep_h5))
data <- connectCounts(c(rep_h5,shuffle_rep_h5), param, regions=DA_anchor, second.regions=GA_anchor)
shuffle_data <- data[,(length(rep_h5)+1):ncol(data)]
data <- data[,1:length(rep_h5)]
DA_anchor <- resize(x = DA_anchor,width = binSize/2+1,fix = 'center')
GA_anchor <- resize(x = GA_anchor,width = binSize/2+1,fix = 'center')
smaller.data <- connectCounts(c(rep_h5,shuffle_rep_h5), param, regions=DA_anchor, second.regions=GA_anchor)
smaller.shuffle_data <- smaller.data[,(length(rep_h5)+1):ncol(smaller.data)]
smaller.data <- smaller.data[,1:length(rep_h5)]
data_trended <- connectCounts(c(rep_h5,shuffle_rep_h5), param, regions=DA_anchor, second.regions=GA_anchor, filter=read_thresh)


######################
### Filter data ######
######################
grid <- expand.grid.unique(cells,cells)

shuffle_idx <- 2*assay(data)/assay(shuffle_data)
is.na(shuffle_idx) <- sapply(shuffle_idx, is.infinite)
shuffle.keep <- abs(log2(rowMaxs(shuffle_idx,na.rm=T)))>log2(shuffle_thresh)
message('Shuffle kept:',paste0(table(shuffle.keep),collapse='-'))

smaller.shuffle_idx <- 2*assay(smaller.data)/assay(smaller.shuffle_data)
is.na(smaller.shuffle_idx) <- sapply(smaller.shuffle_idx, is.infinite)
smaller.shuffle.keep <- rowMaxs(smaller.shuffle_idx,na.rm=T)>shuffle_thresh

# margin.data <- marginCounts(rep_h5, param, width=binSize)
# adjc <- cpm(asDGEList(margin.data), log=TRUE, prior.count=1)
# colnames(adjc) <- paste0(design_matrix$cells,'_',design_matrix$reps)
# pdf(paste0(path,chr,'_marginals.pdf'),height=4,width=nrow(grid)*4)
# par(mfrow=c(1,nrow(grid)))
# for (i in 1:nrow(grid)){
# 	smoothScatter(rowSums(adjc[,grep(grid[i,1],colnames(adjc))],na.rm=T),rowSums(adjc[,grep(grid[i,2],colnames(adjc))],na.rm=T),xlab=grid[i,1], ylab=grid[i,2], main=paste0(grid[i,1],'vs',grid[i,2]))
# }
# dev.off()

ave.ab <- aveLogCPM(asDGEList(data))
avg_thresh <- aveLogCPM(read_thresh, lib.size=mean(data$totals))
count.keep <- ave.ab >= avg_thresh

pdf(paste0(path,chr,'_abundance.pdf'),height=5,width=10)
par(mfrow=c(1,2))
hist(ave.ab,xlab="Average abundance",main=paste0('all bins. Kept: ',summary(count.keep)[3],'; discarded: ',summary(count.keep)[2]))
abline(v=avg_thresh,col='red',lty=5)
hist(ave.ab,xlab="Average abundance",xlim=c(min(ave.ab,na.rm=T),quantile(ave.ab,0.95)),main='0.95 percentile')
abline(v=avg_thresh,col='red',lty=5)
dev.off()


rm(ave.ab)
message('Average read kept:',paste0(summary(count.keep),collapse=' '))

smaller.ave.ab <- aveLogCPM(asDGEList(smaller.data))
smaller.avg_thresh <- aveLogCPM(read_thresh, lib.size=mean(data$totals))
smaller.count.keep <- smaller.ave.ab >= smaller.avg_thresh
rm(smaller.ave.ab)

trended_idx <- matrix(FALSE,nrow=nrow(assay(data)),ncol=nrow(design_matrix))
pdf(paste0(path,chr,'_cisDecay.pdf'),height=5,width=5*length(cells))
par(mfrow=c(1,length(cells)))
for (cell in cells){
	idx <- as.numeric(row.names(design_matrix)[grep(cell,design_matrix$cell)])
	trended <- filterTrended(data[,idx],reference=data_trended[,idx])
	smoothScatter(trended$log.distance, trended$abundances, xlab="Log-Distance", ylab="Normalized abundance",main=cell)
	lines(trended$log.distance[order(trended$log.distance)], trended$threshold[order(trended$log.distance)], col="red", lwd=2)	
	trended_idx[,idx] <- trended$abundances > trended$threshold
}
dev.off()
trended.keep <- rowSums(trended_idx)>1
message('Trended kept:',paste0(summary(trended.keep),collapse='-'))


trended_idx <- matrix(FALSE,nrow=nrow(assay(smaller.data)),ncol=nrow(design_matrix))
for (cell in cells){
	idx <- as.numeric(row.names(design_matrix)[grep(cell,design_matrix$cell)])
	trended <- filterTrended(smaller.data[,idx],reference=data_trended[,idx])
	trended_idx[,idx] <- trended$abundances > trended$threshold
}
smaller.trended.keep <- rowSums(trended_idx)>1

if (shuffle_keep){
	total.keep <- count.keep&shuffle.keep&trended.keep
	smaller.total.keep <- smaller.count.keep&smaller.shuffle.keep&smaller.trended.keep
} else {
	total.keep <- count.keep&trended.keep
	smaller.total.keep <- smaller.count.keep&smaller.trended.keep
}
total.keep[is.na(total.keep)] <- FALSE
smaller.total.keep[is.na(smaller.total.keep)] <- FALSE
message('Total kept:',paste0(summary(total.keep),collapse='-'))

data_orig <- data
smaller.data_orig <- smaller.data

data <- data[total.keep,]
shuffle_idx <- shuffle_data[total.keep,]

smaller.data <- smaller.data[smaller.total.keep,]
smaller.shuffle_idx <- smaller.shuffle_data[smaller.total.keep,]


######################
### Normalize data ###
######################

data <- normOffsets(data, type="loess", se.out=TRUE)
# neardiag <- filterDiag(data, by.dist=binSize*1.5)
# nb.off <- matrix(0, nrow=nrow(data), ncol=ncol(data))
# nb.off[neardiag] <- normOffsets(data[neardiag,], type="loess")
# nb.off[!neardiag] <- normOffsets(data[!neardiag,], type="loess")

smaller.data <- normOffsets(smaller.data, type="loess", se.out=TRUE)
# smaller.neardiag <- filterDiag(smaller.data, by.dist=binSize*0.75)
# smaller.nb.off <- matrix(0, nrow=nrow(smaller.data), ncol=ncol(smaller.data))
# smaller.nb.off[smaller.neardiag] <- normOffsets(smaller.data[smaller.neardiag,], type="loess")
# smaller.nb.off[!smaller.neardiag] <- normOffsets(smaller.data[!smaller.neardiag,], type="loess")

temp <- log2(assay(shuffle_idx)/rowMeans(assay(shuffle_idx),na.rm=T))
is.na(temp) <- sapply(temp, is.infinite)
temp[is.na(temp)] <- 0

smaller.temp <- log2(assay(smaller.shuffle_idx)/rowMeans(assay(smaller.shuffle_idx),na.rm=T))
is.na(smaller.temp) <- sapply(smaller.temp, is.infinite)
smaller.temp[is.na(smaller.temp)] <- 0

ab <- aveLogCPM(asDGEList(data))
o <- order(ab)
pdf(paste0(path,chr,'_adjCounts',cells[1],'.pdf'),height=4,width=12)
par(mfrow=c(1,3))
adj.counts <- cpm(asDGEList(data), log=TRUE)
mval <- adj.counts[,3]-adj.counts[,1]
smoothScatter(ab, mval, xlab="A", ylab="M",  main=paste0(cells[1],"_Rep2vsRep1_orig"))
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")
adj.counts <- log2(assay(data) + 0.5) - assay(data, "offset")/log(2)
mval <- adj.counts[,3]-adj.counts[,1]
smoothScatter(ab, mval, xlab="A", ylab="M", main=paste0(cells[1],"_Rep2vsRep1_loess"))
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")
adj.counts <- log2(assay(data) + 0.5) - temp/log(2)
mval <- adj.counts[,2]-adj.counts[,1]
smoothScatter(ab, mval, xlab="A", ylab="M", main=paste0(cells[1],"_Rep2vsRep1_shaman"))
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")
dev.off()

# adj.counts <- log2(assay(data) + 0.5) - temp/log(2)
# colnames(adj.counts) <- paste0(design_matrix$cells,'_',design_matrix$reps)
# pdf(paste0(path,name,'_adjCounts_allShaman.pdf'),height=4,width=nrow(grid)*4)
# par(mfrow=c(1,nrow(grid)))
# for (i in 1:nrow(grid)){
# 	if ((!grepl('ncx',grid[i,1]))&(!grepl('ncx',grid[i,2]))) {
# 	idx <- c(as.numeric(row.names(design_matrix)[intersect(grep(grid[i,1],design_matrix$cells),grep('ncx',design_matrix$cells,invert=T))]),as.numeric(row.names(design_matrix)[intersect(grep(grid[i,2],design_matrix$cells),grep('ncx',design_matrix$cells,invert=T))]))
# 	 } else if ((!grepl('ncx',grid[i,1]))&(grepl('ncx',grid[i,2]))) {
# 	idx <- c(as.numeric(row.names(design_matrix)[intersect(grep(grid[i,1],design_matrix$cells),grep('ncx',design_matrix$cells,invert=T))]),as.numeric(row.names(design_matrix)[grep(grid[i,2],design_matrix$cells)]))
# 	 } else if ((grepl('ncx',grid[i,1]))&(!grepl('ncx',grid[i,2]))) {
# 	idx <- c(as.numeric(row.names(design_matrix)[grep(grid[i,2],design_matrix$cells)]),as.numeric(row.names(design_matrix)[intersect(grep(grid[i,2],design_matrix$cells),grep('ncx',design_matrix$cells,invert=T))]))
# 	 } else {
# 	idx <- c(as.numeric(row.names(design_matrix)[grep(grid[i,1],design_matrix$cells)]),as.numeric(row.names(design_matrix)[grep(grid[i,2],design_matrix$cells)]))
# 	 }
# 	ab <- aveLogCPM(asDGEList(data[,idx]))
# 	o <- order(ab)
# 	mval <- rowSums(adj.counts[,grep(grid[i,1],colnames(adj.counts))],na.rm=T)-rowSums(adj.counts[,grep(grid[i,2],colnames(adj.counts))],na.rm=T)
# 	smoothScatter(ab, mval, xlab="A", ylab="M", main=paste0(grid[i,1],'vs',grid[i,2]))
# 	fit <- loessFit(x=ab, y=mval)
# 	lines(ab[o], fit$fitted[o], col="red")
# }
# dev.off()

# adj.counts <- log2(assay(data) + 0.5) - nb.off/log(2)
# colnames(adj.counts) <- paste0(design_matrix$cells,'_',design_matrix$reps)
# pdf(paste0(path,name,'_adjCounts_allLoess.pdf'),height=4,width=nrow(grid)*4)
# par(mfrow=c(1,nrow(grid)))
# for (i in 1:nrow(grid)){
# 	if ((!grepl('ncx',grid[i,1]))&(!grepl('ncx',grid[i,2]))) {
# 	idx <- c(as.numeric(row.names(design_matrix)[intersect(grep(grid[i,1],design_matrix$cells),grep('ncx',design_matrix$cells,invert=T))]),as.numeric(row.names(design_matrix)[intersect(grep(grid[i,2],design_matrix$cells),grep('ncx',design_matrix$cells,invert=T))]))
# 	 } else if ((!grepl('ncx',grid[i,1]))&(grepl('ncx',grid[i,2]))) {
# 	idx <- c(as.numeric(row.names(design_matrix)[intersect(grep(grid[i,1],design_matrix$cells),grep('ncx',design_matrix$cells,invert=T))]),as.numeric(row.names(design_matrix)[grep(grid[i,2],design_matrix$cells)]))
# 	 } else if ((grepl('ncx',grid[i,1]))&(!grepl('ncx',grid[i,2]))) {
# 	idx <- c(as.numeric(row.names(design_matrix)[grep(grid[i,2],design_matrix$cells)]),as.numeric(row.names(design_matrix)[intersect(grep(grid[i,2],design_matrix$cells),grep('ncx',design_matrix$cells,invert=T))]))
# 	 } else {
# 	idx <- c(as.numeric(row.names(design_matrix)[grep(grid[i,1],design_matrix$cells)]),as.numeric(row.names(design_matrix)[grep(grid[i,2],design_matrix$cells)]))
# 	 }
# 	ab <- aveLogCPM(asDGEList(data[,idx]))
# 	o <- order(ab)
# 	mval <- rowSums(adj.counts[,grep(grid[i,1],colnames(adj.counts))],na.rm=T)-rowSums(adj.counts[,grep(grid[i,2],colnames(adj.counts))],na.rm=T)
# 	smoothScatter(ab, mval, xlab="A", ylab="M", main=paste0(grid[i,1],'vs',grid[i,2]))
# 	fit <- loessFit(x=ab, y=mval)
# 	lines(ab[o], fit$fitted[o], col="red")
# }
# dev.off()



# corrected <- correctedContact(data_orig, average=FALSE)
# anchor1.bias <- corrected$bias[anchors(data_orig, type="first", id=TRUE),]
# anchor2.bias <- corrected$bias[anchors(data_orig, type="second", id=TRUE),]
# iter.off <- log(anchor1.bias * anchor2.bias)[total.keep,]

######################
### Perform Testing ##
######################

for (shuff_norm in c(FALSE,TRUE)){
	path=paste0(main_f,'analysis/diffHiC/',name,'/',binSize/1000,'kb/',ifelse(shuff_norm,'/shaman/','/loess/'))
	dir.create(path, showWarnings = FALSE,recursive=T)

	if (shuff_norm){
		nb.off <- temp
		smaller.nb.off <- smaller.temp
	}
	y <- asDGEList(data)
#	y <- scaleOffset(y, nb.off)
	y <- estimateDisp(y, design)

	smaller.y <- asDGEList(smaller.data)
#	smaller.y <- scaleOffset(smaller.y, smaller.nb.off)
	smaller.y <- estimateDisp(smaller.y, design)

	pdf(paste0(path,chr,'_dispNB.pdf'),height=4,width=4)
	plotBCV(y)
	dev.off()

	fit <- glmQLFit(y, design, robust=TRUE)
	pdf(paste0(path,chr,'_dispQL.pdf'),height=4,width=4)
	plotQLDisp(fit)
	dev.off()

	smaller.fit <- glmQLFit(smaller.y, design, robust=TRUE)

	for (i in 1:nrow(grid)){
		name2 <- paste0(grid[i,2],'-',grid[i,1])
		message('Working on ',name2)
		result <- glmQLFTest(fit, contrast=makeContrasts(constrasts=name2, levels=design))
		adj.p <- p.adjust(result$table$PValue, method="BH")
		useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
		inter.frame <- data.frame(DA_anchor=names(regions(data))[anchorIds(data,type='first')],GA_anchor=names(regions(data))[anchorIds(data,type='second')],stringsAsFactors = F)
			#	inter.frame <- as.data.frame(interactions(data))[,useful.cols]
		results.r <- data.frame(inter.frame, result$table, FDR=adj.p,stringsAsFactors = F)
		results.r <- results.r[order(results.r$PValue),]
		results.r <- unique(results.r)
		results.r$FDR <- p.adjust(results.r$PValue, method="BH")
	#	results.r <- results.r[results.r$FDR<0.05,]
		path_t <- paste0(path,name2,'/direct/')
		dir.create(path_t, showWarnings = FALSE,recursive=T) 
		write.table(results.r, file=paste0(path_t,chr,'.txt'), sep="\t", quote=FALSE, row.names=FALSE)
    next
		clustered <- clusterPairs(data, tol=1, upper=binSize*10)
		tabcluster <- combineTests(clustered$indices[[1]], result$table)
		inter.frame <- as.data.frame(clustered$interactions)[,useful.cols]
		results.i <- data.frame(inter.frame, tabcluster)
		results.i <- results.i[order(results.i$PValue),]
		results.i <- results.i[results.i$FDR<0.05,]
		path_t <- paste0(path,name2,'/clusters/')
		dir.create(path_t, showWarnings = FALSE,recursive=T) 
		write.table(results.i, file=paste0(path_t,chr,'_independent.txt'), sep="\t", quote=FALSE, row.names=FALSE)
	
		clustered.sig <- diClusters(data, result$table, target=0.05, cluster.args=list(tol=1))
		tabcom <- combineTests(clustered.sig$indices[[1]], result$table)
		tabbest <- getBestTest(clustered.sig$indices[[1]], result$table)
		check <- (nrow(tabcom)==length(tabbest$logFC))&(nrow(tabcom)==length(clustered.sig$FDR))&(length(tabbest$logFC)==length(clustered.sig$FDR))
		if (check){
			tabstats <- data.frame(tabcom[,1:4], logFC=tabbest$logFC, FDR=clustered.sig$FDR)
			results.d <- data.frame(as.data.frame(clustered.sig$interactions)[,useful.cols], tabstats)	
			results.d <- results.d[order(results.d$PValue),]
			write.table(results.d, file=paste0(path_t,chr,'_DI.txt'), sep="\t", quote=FALSE, row.names=FALSE)
		}
		smaller.result <- glmQLFTest(smaller.fit, contrast=makeContrasts(constrasts=name2, levels=design))
		merged.data <- list(data, smaller.data)
		merged.results <- list(result$table, smaller.result$table)
		clustered.mult <- diClusters(merged.data, merged.results,target=0.05, cluster.args=list(tol=1))	
	
	#	boxed <- boxPairs(larger=data, smaller=smaller.data)
		cons <- consolidatePairs(clustered.mult$indices, merged.results)
	#	cons2 <- consolidatePairs(boxed$indices, merged.results)
		results.b <- data.frame(as.data.frame(clustered.mult$interactions)[,useful.cols], cons)
		o.b <- order(results.b$PValue)	
		results.b <- results.b[o.b,]
		path_t <- paste0(path,name2,'/nested/')
		dir.create(path_t, showWarnings = FALSE,recursive=T) 
		write.table(results.b, file=paste0(path_t,chr,'_boxed.txt'), sep="\t", quote=FALSE, row.names=FALSE)	
	
		inside <- getBestTest(clustered.mult$indices[[2]], smaller.result$table)
		best.interactions <- interactions(smaller.data)[inside$best,]
		inter.frame <- as.data.frame(best.interactions)[,useful.cols[-c(1,4)]]
		nested <- data.frame(inter.frame, inside[,c("logFC", "F")])
		expanded <- rep(NA, nrow(results.b)) # As some parents do not have nested elements.
		expanded[as.integer(rownames(inside))] <- seq_len(nrow(inside))
		results.b <- data.frame(results.b, best=nested[expanded,])	
		write.table(results.b[o.b,], file=paste0(path_t,chr,'_boxedBest.txt'), sep="\t", quote=FALSE, row.names=FALSE)
	}	
}	
	