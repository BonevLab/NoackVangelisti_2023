#!/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
main_f <- as.character(args[7])
source(paste0(main_f,'config.R'))

construct.grid = function(interv1,interv2,min_dist,max_dist){
  return(ddply(interv1, .(chrom), function(i1) {
	i2 = interv2[as.character(interv2$chrom) == as.character(i1$chrom[1]),]
	if (nrow(i2) ==0) {
		return(c())
	}
	g = expand.grid(1:nrow(i1), 1:nrow(i2))
  	g = g[i2$start[g$Var2]-i1$start[g$Var1] > min_dist & i2$start[g$Var2]-i1$start[g$Var1] < max_dist,]
	grid = cbind(i1[g$Var1,c("chrom", "start", "end")], i2[g$Var2,c("chrom", "start", "end")])
	colnames(grid) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
	return(grid)
   })[,-1] )
 }

construct.grid2 = function(interv1,interv2,min_dist,max_dist){
  return(ddply(interv1, .(chrom), function(i1) {
	i2 = interv2[as.character(interv2$chrom) == as.character(i1$chrom[1]),]
	if (nrow(i2) ==0) {
		return(c())
	}
	g = expand.grid(1:nrow(i1), 1:nrow(i2))
  	g = g[i2$start[g$Var2]-i1$start[g$Var1] > min_dist &
        	i2$start[g$Var2]-i1$start[g$Var1] < max_dist,]
	  grid = cbind(i1[g$Var1,c("chrom", "start", "end")], i2[g$Var2,c("chrom", "start", "end")])
	  colnames(grid) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
	 return(grid)
   })[,-1] )
}

construct.grid.trans = function(interv1,interv2){
  return(ddply(interv1, .(chrom), function(i1) {
	i2 = interv2[as.character(interv2$chrom) != as.character(i1$chrom[1]),]
	if (nrow(i2) ==0) {
		return(c())
	}
	g = expand.grid(1:nrow(i1), 1:nrow(i2))
	grid = cbind(i1[g$Var1,c("chrom", "start", "end")], i2[g$Var2,c("chrom", "start", "end")])
	colnames(grid) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
	return(grid)
   })[,-1] )
 }


file1 <- args[1]
file2 <- args[2]
name <- args[3]
grid_mode <- args[4]
interval_window <- as.numeric(args[5])
hic_f <- args[6]
domains_f <- args[8]
min_dist <- as.numeric(args[9])
max_dist <- as.numeric(args[10])
res <- as.numeric(args[11]) 
window_f <- as.numeric(args[12]) 
domain_mode <- args[13]
save_folder <- args[14]

dir.create(paste0(main_f,'data/APA/',save_folder), showWarnings = FALSE,recursive=T)


juicebox_jar='/home/hpc/bonev/juicer/scripts/common/juicer_tools.jar'

source_f <- gintervals.all()


if (grepl('Peak',file1)|grepl('bed',file1)) {
	peaks <- read.table(file1)
	interv1 <- gintervals(chroms=peaks$V1,starts=peaks$V2,ends=peaks$V3)
	interv1 <- interv1[interv1$chrom %in% gintervals.all()$chrom,]
	intervals1 <- interv1[(as.character(interv1$chrom)!='chrM'|as.character(interv1$chrom)!='chrY'),]
} else {
	intervals1 = get(load(file=file1))
}
intervals1 <- intervals1[complete.cases(intervals1),]
if (interval_window!=0) {intervals1 <- intervals.expand(intervals1, interval_window/2)} else {intervals1 <- intervals.expand(intervals1, 1000)}
if (grid_mode!="merged") {intervals1 <- gintervals.canonic(intervals1)}
if (interval_window!=0) {intervals1 <- intervals.normalize(intervals1,interval_window)}                 ## else { intervals1 <- intervals1[(intervals1[,3]-intervals1[,2])>=max(dist_f),] }
interv1 <- intervals1
message('working with ',nrow(interv1),' intervals')
interv1$start <- round_any(interv1$start,res)
interv1$end <- round_any(interv1$end,res)

if (grepl('Peak',file2)|grepl('bed',file2)) {
	peaks <- read.table(file2)
	interv2 <- gintervals(chroms=peaks$V1,starts=peaks$V2,ends=peaks$V3)
	interv2 <- interv2[interv2$chrom %in% gintervals.all()$chrom,]
	intervals2 <- interv2[(as.character(interv2$chrom)!='chrM'|as.character(interv2$chrom)!='chrY'),]
} else { intervals2 = get(load(file=file2))}
intervals2 <- intervals2[complete.cases(intervals2),]
if (interval_window!=0) { intervals2 <- intervals.expand(intervals2, interval_window/2) } else {intervals2 <- intervals.expand(intervals2, 1000)}           ##  else {intervals2 <- intervals.expand(intervals2, mean(intervals2[,3]-intervals2[,2]))}
if (grid_mode!="merged") {intervals2 <- gintervals.canonic(intervals2)}
if (interval_window!=0) {intervals2 <- intervals.normalize(intervals2,interval_window)} 
interv2 <- intervals2
message('working with ',nrow(interv2),' intervals')
interv2$start <- round_any(interv2$start,res)
interv2$end <- round_any(interv2$end,res)


if (!gintervals.exists(domains_f)) {
	domains <- gintervals.ls(domains_f)
	domains <- domains[grep(cell,domains)]
	domains <- gintervals.load(domains)[,1:3]
} else {	 
	domains <- gintervals.load(domains_f)[,1:3]
}

if (grid_mode=="cis"){
message("Grid mode = cis")
	grid <- construct.grid(interv1,interv2,min_dist,max_dist)
	if (domain_mode=='intra'){
		dom = gintervals.2d(chroms1=domains$chrom, starts1=domains$start, ends1=domains$end,chroms2=domains$chrom, starts2=domains$start, ends2=domains$end)
		grid <- gintervals.intersect(grid,dom)
	} else if (domain_mode=='inter'){
		inter_2d = construct.grid2(domains,domains,min_dist,max_dist)			
		grid <- gintervals.intersect(grid,inter_2d)
	}
}

if (grid_mode=="trans"){
message("Grid mode = trans")
grid <- construct.grid.trans(interv1,interv2)
}

colnames(grid) <- c("chr1", "x1", "x2", "chr2", "y1", "y2")

grid$color <- '0,255,0'

file_f <- paste0(main_f,"/data/APA/",name)
write.table(grid,file_f,quote=F,col.names=T,row.names=F,sep='\t')


command_f <- paste0('java -Djava.awt.headless=true -jar ',juicebox_jar,' apa -n ',round(min_dist/res),' -x ',round(max_dist/res),' -w ',round(window_f/res),' -r ',res,' -c ',paste(gsub('chr','',gintervals.all()$chrom),collapse=','),' -q 3 -k KR ',ifelse(grid_mode=='trans','-e ',''),'-u ',hic_f,' ',file_f,' ',paste0(main_f,'data/APA/',save_folder,'/')) 
print(command_f)
system(command_f)