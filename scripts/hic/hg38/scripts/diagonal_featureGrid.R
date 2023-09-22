diagonal_contactMatrix <- function(track, shuffled_track, v_scoreTracks, filter_value, interval, expand, expand_res, min_dist, max_dist,regularization=30,interval_name=NULL,consider_strand=FALSE,track_db=trackdb,tss=tss_f){
	expand1 <- expand
	expand2 <- expand
	o = matrix(0, nrow=length(expand1)-1, ncol=length(expand2)-1)
	e = matrix(0, nrow=length(expand1)-1, ncol=length(expand2)-1)
	filtered_int_name <- paste0('intervals.',interval_name,"_",as.integer(expand_res),"_",as.integer(min_dist),"_",as.integer(max_dist),"_",as.integer(max(expand2)),"_",filter_value)
	grid_f <- paste0(track_db,'/tracks/intervals/',filtered_int_name)
	message(filtered_int_name)
	band = c(-max_dist, -min_dist)

	if (!file.exists(grid_f)){
		if (grepl('Peak',interval_f)|grepl('bed',interval_f)) {
			peaks <- read.table(interval_f)
			interval <- gintervals(chroms=peaks$V1,starts=peaks$V2,ends=peaks$V3)
			interval <- interval[interval$chrom %in% gintervals.all()$chrom,]
			interval <- interval[(as.character(interval$chrom)!='chrM'|as.character(interval$chrom)!='chrY'),]
		} else {
			interval <- get(load(file=interval_f))
		}
		interval <- intervals.expand(interval, expand_res/2)
		interval <- gintervals.canonic(interval)
		interval <- unique(intervals.centers(interval))
		grid <- gintervals.2d(interval[,1],interval[,2],interval[,3],interval[,1],interval[,2],interval[,3])
		if (filter_value!=0) {
			filter_vtrack_contraint = paste(v_scoreTracks, ">", filter_value, collapse =" | ")
			if (nrow(grid)>100000){
				chunk <- 100000
				n <- nrow(grid)
				r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
				d <- split(grid,r)
				d1 <- lapply(d, function(x) gscreen(filter_vtrack_contraint, intervals=x, iterator=x, band=band))
				grid <- do.call("rbind", d1)
				row.names(grid) <- seq(1,nrow(grid))
				rm(d)
				rm(d1)
				rm(r)	
			} else {
				grid = gscreen(filter_vtrack_contraint, intervals=grid, iterator=grid, band=band)
			}
		 }
		 message(paste("found ", nrow(grid), " regions after score filter"))
		 save(grid,file=grid_f)
	} else {
		grid <- get(load(grid_f))
		message(paste("found ", nrow(grid), " regions after score filter"))
	}
	if (consider_strand){
  		tss <- gintervals.load(tss)
  		temp_tss <- grid[,1:3]
		colnames(temp_tss) <- c('chrom','start','end') 
		temp_tss <- gintervals.neighbors(temp_tss,tss) 
		grid$type=1
		grid$strand <- temp_tss$strand
		colnames(grid) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2","type","strand")
	}
	interv_set_out = gsub("+", "", paste(track,paste(rev(-(band)), collapse="_"), sep="_"), fixed=T)
	if (!gintervals.exists(interv_set_out)) {
		giterator.intervals(track, band=band, intervals.set.out=interv_set_out)
	}

	a <- gintervals.neighbors(interv_set_out, grid, mindist1=min(expand1), maxdist1=max(expand1), mindist2=min(expand2), maxdist2=max(expand2))
	colnames(a)[7:(ncol(a)-2)] = paste0("grid.", colnames(grid))
	gc(verbose=F)
	if (!consider_strand){
		o = table(cut(a$start1 - a$grid.start1, breaks=expand1, include.lowest=TRUE), cut(a$start2 - a$grid.start2, breaks=expand2, include.lowest=TRUE))
 	} else {
 		rotate <- function(x) t(apply(x, 2, rev))
		a1 <- a[a$grid.type==1&a$grid.strand==1,]
		o = table(cut(a1$start1 - a1$grid.start1, breaks=expand1, include.lowest=TRUE), cut(a1$start2 - a1$grid.start2, breaks=expand2, include.lowest=TRUE))
		message(' N of obs contacts type1 strand1:',nrow(a1))
		for (type in c(1,2)){
			for (strand in c(1,(-1))){
				a1 <- a[a$grid.type==type&a$grid.strand==strand,]
				message('type:',type,', strand:',strand,' N of obs contacts:',nrow(a1))
				o1 = table(cut(a1$start1 - a1$grid.start1, breaks=expand1, include.lowest=TRUE), cut(a1$start2 - a1$grid.start2, breaks=expand2, include.lowest=TRUE))
				if (type==1&strand==(-1)) { o1 <- t(rotate(rotate(o1)))} 
				else if (type==1&strand==1) { next}		
				colnames(o1) <- colnames(o)
				rownames(o1) <- rownames(o)
				o <- o + o1
			}
		}
	}
 

  	gc(verbose=F)
	interv_set_out = gsub("+", "", paste(shuffled_track,paste(rev(-(band)), collapse="_"), sep="_"), fixed=T)
	if (!gintervals.exists(interv_set_out)) {
		giterator.intervals(shuffled_track, band=band, intervals.set.out=interv_set_out)
	}

	a <- gintervals.neighbors(interv_set_out, grid, mindist1=min(expand1), maxdist1=max(expand1), mindist2=min(expand2), maxdist2=max(expand2))
  	colnames(a)[7:(ncol(a)-2)] = paste0("grid.", colnames(grid))
  	gc(verbose=T)
  	if (!consider_strand){
  	e = table(cut(a$start1 - a$grid.start1, breaks=expand1, include.lowest=TRUE), cut(a$start2 - a$grid.start2, breaks=expand2, include.lowest=TRUE))
	} else {
	 ### Accounting for strand - exp
 		a1 <- a[a$grid.type==1&a$grid.strand==1,]
		e = table(cut(a1$start1 - a1$grid.start1, breaks=expand1, include.lowest=TRUE), cut(a1$start2 - a1$grid.start2, breaks=expand2, include.lowest=TRUE))
		message(' N of exp contacts type1 strand1:',nrow(a1))
		for (type in c(1,2)){
			for (strand in c(1,-1)){
				a1 <- a[a$grid.type==type&a$grid.strand==strand,]
				e1 = table(cut(a1$start1 - a1$grid.start1, breaks=expand1, include.lowest=TRUE), cut(a1$start2 - a1$grid.start2, breaks=expand2, include.lowest=TRUE))
				if (type==1&strand==(-1)) { e1 <- t(rotate(rotate(e1)))} 
				else if (type==1&strand==1) { next}		
				colnames(e1) <- colnames(e)
				rownames(e1) <- rownames(e)
				e <- e + e1
				message('type:',type,', strand:',strand,' N of exp contacts:',nrow(a1))
			}
		}
	}	
 	
	o[o < regularization] = NA
	e[is.na(o) | e < regularization] = NA
	total_obs = sum(o, na.rm=TRUE)
	total_exp = sum(e, na.rm=TRUE)
	obs = o / total_obs
	exp = e / total_exp
	return(list(obs=obs, exp=exp, total_obs=total_obs, total_exp=total_exp, grid_size=nrow(grid)));
}


args = commandArgs(trailingOnly=TRUE)


main_f=args[1]
track = args[2]
shuffled_track = paste0(track, "_shuffle")
score_track = args[3] 
filter_value = as.numeric(args[4])
interval_f = args[5]
data_file = args[6]
min_dist <- as.numeric(args[7])
expand_range = as.numeric(args[8])
expand_res = as.numeric(args[9])
max_dist <- expand_range
grid_mode= args[10]
interval_name = args[11]
source(paste0(main_f,'config.R'))
options(gmax.data.size=2e8)
message(paste("expand_range=", expand_range))
expand = seq(0-expand_range, expand_range, by=expand_res)

scoreTracks <- paste0("hic.",cells,".score_",score_track)                       #
if (filter_value!=0){
	 v_scoreTracks <- paste("v_", scoreTracks, sep="")
	   for (i in 1:length(v_scoreTracks)) {
		gvtrack.create(v_scoreTracks[i], scoreTracks[i], 'max')
		gvtrack.iterator.2d(v_scoreTracks[i], sshift1=min(expand), eshift1=max(expand), sshift2=min(expand), eshift2=max(expand))
	   }
} else { v_scoreTracks <- '' }

obs = list()
exp = list()
grid_size = c()
total_obs = c()
total_exp = c()
print(interval_f)
grid = diagonal_contactMatrix(track=track, shuffled_track=shuffled_track, v_scoreTracks=v_scoreTracks, filter_value=filter_value,interval=interval_f, expand=expand, expand_res=expand_res, min_dist=min_dist, max_dist=max_dist,regularization=0,interval_name=interval_name,consider_strand=as.logical(grid_mode))
obs[[length(obs)+1]] = grid$obs
exp[[length(exp)+1]] = grid$exp
grid_size = c(grid_size, grid$grid_size)
total_obs = c(total_obs, grid$total_obs)
total_exp = c(total_exp, grid$total_exp)
message(paste("grid size=", grid$grid_size, ", max obs=", max(log10(grid$obs),na.rm=TRUE), ", max obs/exp=", max(log2(grid$obs/grid$exp),na.rm=TRUE)))

   
grid = list(obs = obs, exp=exp, grid_size=grid_size, total_obs=total_obs, total_exp = total_exp)
save(grid, file=data_file)
