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


dist_f = 10^seq(3.9,8,by=0.1)
dist = 10^seq(3.9,8,by=0.1)[-1]
min_dist <- min(dist_f)
max_dist <- max(dist_f)


file1 <- args[1]
file2 <- args[2]
name <- args[3]
grid_mode <- args[4]
interval_window <- as.numeric(args[5])
cells <- args[6]
cells <- c(unlist(strsplit(cells,'\\,')))
domains_f <- args[8]

message('Working on cells: ',paste0(cells,collapse=','))

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
#if (grid_mode!="merged") {intervals1 <- gintervals.canonic(intervals1)}
intervals1 <- gintervals.canonic(intervals1)
if (interval_window!=0) {intervals1 <- intervals.normalize(intervals1,interval_window)}                 ## else { intervals1 <- intervals1[(intervals1[,3]-intervals1[,2])>=max(dist_f),] }
interv1 <- intervals1
message('working with ',nrow(interv1),' intervals')
if (grepl('Peak',file2)|grepl('bed',file2)) {
	peaks <- read.table(file2)
	interv2 <- gintervals(chroms=peaks$V1,starts=peaks$V2,ends=peaks$V3)
	interv2 <- interv2[interv2$chrom %in% gintervals.all()$chrom,]
	intervals2 <- interv2[(as.character(interv2$chrom)!='chrM'|as.character(interv2$chrom)!='chrY'),]
} else { intervals2 = get(load(file=file2))}
intervals2 <- intervals2[complete.cases(intervals2),]
if (interval_window!=0) { intervals2 <- intervals.expand(intervals2, interval_window/2) } else {intervals2 <- intervals.expand(intervals2, 1000)}           ##  else {intervals2 <- intervals.expand(intervals2, mean(intervals2[,3]-intervals2[,2]))}
#if (grid_mode!="merged") {intervals2 <- gintervals.canonic(intervals2)}
intervals2 <- gintervals.canonic(intervals2)
if (interval_window!=0) {intervals2 <- intervals.normalize(intervals2,interval_window)} 
interv2 <- intervals2
message('working with ',nrow(interv2),' intervals')


	if (grid_mode=="1D"){
		grid <- construct.grid(interv1,interv2,min_dist,max_dist)
	}

	 if (grid_mode=="2D"){
	 message("Grid mode = 2D")
		grid1 <- construct.grid(interv1,interv2,min_dist,max_dist)
		grid2 <- construct.grid(interv2,interv1,min_dist,max_dist)
		grid <- unique(rbind(grid1,grid2))
	 }

if (grid_mode=="merged") {
	 message("Grid mode = merged")

		grid <- cbind(as.data.frame(interv1[,1:3]),as.data.frame(interv2[,1:3]))
		message(nrow(grid))
# 		for (i in 1:nrow(grid)){
# 		 if (grid[i,2] > grid[i,5]) {
# 			 grid[i,] <- grid[i,c(4,5,6,1,2,3)]
# 			}
# 		} 
		colnames(grid) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
		grid <- grid[((grid$start2 - grid$start1) >= min_dist) & ((grid$start2 - grid$start1) <= max_dist),]
}



source_f <- interv1
message('Working on grid with size: ',nrow(grid)) 

chroms <- gintervals.all()
interv1_marg <- adply(interv1,1,function(x){
	return(chroms[chroms$chrom%in%x$chrom,])
})
interv1_marg <- gintervals.2d(interv1[,1],interv1[,2],interv1[,3],interv1_marg[,1],interv1_marg[,2],interv1_marg[,3])
interv2_marg <- adply(interv2,1,function(x){
	return(chroms[chroms$chrom%in%x$chrom,])
})
interv2_marg <- gintervals.2d(interv2[,1],interv2[,2],interv2[,3],interv2_marg[,1],interv2_marg[,2],interv2_marg[,3])



stats_df <- as.data.frame(matrix(nrow=length(all_tracks),ncol=5))
colnames(stats_df) <- c("obs_intra","obs_inter","exp_intra","exp_inter",'coverage')

cis_decay <- list()
score_list <- list()
for (cell in cells){
	tracks = all_tracks[grep(cell,all_tracks)]
	shuffled_tracks = paste0(tracks,'_shuffle')
	if (length(tracks)!=length(shuffled_tracks)){
	message('Check your tracks input - something went wrong')
	stop()
	}
	if(grepl('chromosome',domains_f)) {
		domains <- gintervals.all()
	} else if (!gintervals.exists(domains_f)) {
		domains <- gintervals.ls(domains_f)
		domains <- domains[grep(cell,domains)]
		domains <- gintervals.load(domains)
	} else {	 
		domains <- gintervals.load(domains_f)
	}
	temp_list <- list()


	
	######### Score Calculation #####	
	if (grid_mode!="merged") {
	  min_dist <- 5e4
	  max_dist <- 2e6
		score_track <- gtrack.ls(paste0('hic.',cell),'.score_k')
		gvtrack.create("v_score", score_track, "max")
		score_grid <- construct.grid(interv1,interv2,min_dist,max_dist)
		dom = gintervals.2d(chroms1=domains$chrom, starts1=domains$start, ends1=domains$end,chroms2=domains$chrom, starts2=domains$start, ends2=domains$end)
		intra_grid <- gintervals.intersect(score_grid,dom)
		inter_2d = construct.grid2(domains,domains,min_dist,max_dist)			
		inter_grid <- gintervals.intersect(score_grid,inter_2d)
		intra_scores <- gextract("v_score", intra_grid, iterator=intra_grid)
		inter_scores <- gextract("v_score", inter_grid, iterator=inter_grid)
		gvtrack.rm("v_score")
		score_list[[cell]]$intra <- intra_scores
		score_list[[cell]]$inter <- inter_scores
	}

	save(score_list,file=paste0(main_f,"/data/cis_decay/",name,'_score'))
	######### Score Calculation #####
		
	for (i in 1:length(tracks)){
		message('Working on track: ',tracks[i])
		obs_cis_decay <- as.data.frame(gcis_decay(tracks[i], dist_f, source_f, domains, grid))
		obs_cis_decay$distance <- row.names(obs_cis_decay)
		exp_cis_decay <- as.data.frame(gcis_decay(shuffled_tracks[i], dist_f, source_f, domains, grid))/2
		exp_cis_decay$distance <- row.names(exp_cis_decay)
		merged <- merge(obs_cis_decay,exp_cis_decay,by="distance",sort=FALSE)
		colnames(merged) <- c("distance","obs_intra","obs_inter","exp_intra","exp_inter")
		temp_list[[i]] <- merged
		names(temp_list)[i] <- unlist(strsplit(tracks[i],'_'))[2]

		 gvtrack.create("v_sum", tracks[i], "weighted.sum")
		 cov1 <- sum(gextract("v_sum", interv1_marg, iterator=interv1_marg, band=c(-max(gintervals.all()$end), -min_dist))$v_sum,na.rm=T)
		 cov2 <- sum(gextract("v_sum", interv2_marg, iterator=interv2_marg, band=c(-max(gintervals.all()$end), -min_dist))$v_sum,na.rm=T)
		 gvtrack.rm("v_sum")

		stats_df[match(tracks[i],all_tracks),] <- c(t(colSums(merged[,-1])),(cov1+cov2))
		rownames(stats_df)[match(tracks[i],all_tracks)] <- tracks[i]
	}
	cis_decay[[cell]] <- temp_list
	cis_decay[['stats']] <- stats_df
	save(cis_decay,file=paste0(main_f,"/data/cis_decay/",name))
}
# 
# cis_decay[['stats']] <- stats_df
# 
# save(cis_decay,file=paste0(main_f,"/data/cis_decay/",name))
