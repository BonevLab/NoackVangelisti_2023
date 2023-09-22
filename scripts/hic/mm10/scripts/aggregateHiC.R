args = commandArgs(trailingOnly=TRUE)
if (length(args) < 10) {
stop("Usage: generate_feature_grid <trackdb> <track> <source track> <score filter> <feature1 data> <feature2 data>
	<output_data> <range> <resolution> <output_image> <domains>")
}

construct.grid = function(interv1,interv2,min_dist,max_dist){
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
	  grid$start1 = grid$start1 + floor((grid$end1-grid$start1)/2)
	  grid$start2 = grid$start2 + floor((grid$end2-grid$start2)/2)
	  grid$end1 = grid$start1+1
	  grid$end2 = grid$start2+1

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
 
norm_contact_grid_by_shuffled_filter_domains <- function(track, shuffled_track, v_scoreTracks, filter_value, interv1, interv2, expand1, expand2, min_dist, max_dist,regularization=30,dom,domains,grid_mode='1D',interval_name=NULL) {
	message(1)
	o = matrix(0, nrow=length(expand1)-1, ncol=length(expand2)-1)
	e = matrix(0, nrow=length(expand1)-1, ncol=length(expand2)-1)
	
	filtered_int_name <- paste0('intervals.',interval_name,"_",as.integer(min_dist),"_",as.integer(max_dist),"_",as.integer(max(expand2)),"_",domains,"_",grid_mode,"_",filter_value)
	grid_f <- paste0('/work/project/Cavalli-mammals/boyan/HiC/trackdb/hg19/tracks/intervals/',filtered_int_name)
	message(filtered_int_name)
	message(paste0("Grid mode = ",grid_mode))
	band = c(-(max_dist-min(expand1)+max(expand2)), -max(min_dist-max(expand1)+min(expand2), 10000))
	
	if (max_dist<=2000000) { 
			band=c(-2040000,-10000) 
		} else {
			band=c(-10040000,-1960000) 
		} 
	if (max_dist==50e6 & min_dist==10e6  ) { band=c(-50100000,-9600000) }
	
	if (max(expand1)!=20000){
		min_dist <- min_dist+max(expand1)
		max_dist <- max_dist-max(expand1)
	}
	if (filter_value==0){
		  message(paste0("Working on band: ",-band[1]," to ",-band[2]))
		  message("building grid...")
  
		if (grid_mode=="1D"){
			grid <- construct.grid(interv1,interv2,min_dist,max_dist)
		} else if (grid_mode=="2D"){
			message("Grid mode = 2D")
			grid1 <- construct.grid(interv1,interv2,min_dist,max_dist)
			grid1$type <- 1
			temp_tss <- grid1[,4:6]
			colnames(temp_tss) <- c('chrom','start','end') 
			temp_tss <- gintervals.neighbors(temp_tss,tss)
			grid1$strand <- temp_tss$strand
			grid2 <- construct.grid(interv2,interv1,min_dist,max_dist)
			grid2$type <- 2
			temp_tss <- grid2[,1:3]
			colnames(temp_tss) <- c('chrom','start','end') 
			temp_tss <- gintervals.neighbors(temp_tss,tss)
			grid2$strand <- temp_tss$strand
			grid <- unique(rbind(grid1,grid2))
			temp <- grid
		} else if (grid_mode=="2Dup"){
			message("Grid mode = 2Dup")
			grid2 <- construct.grid(interv2,interv1,min_dist,max_dist)
			grid2$type <- 2
			temp_tss <- grid2[,1:3]
			colnames(temp_tss) <- c('chrom','start','end') 
			temp_tss <- gintervals.neighbors(temp_tss,tss)
			grid2$strand <- temp_tss$strand
			grid <- unique(grid2)
			temp <- grid
		} else if (grid_mode=="2Ddown"){
			message("Grid mode = 2Ddown")
			grid1 <- construct.grid(interv1,interv2,min_dist,max_dist)
			grid1$type <- 1
			temp_tss <- grid1[,4:6]
			colnames(temp_tss) <- c('chrom','start','end') 
			temp_tss <- gintervals.neighbors(temp_tss,tss)
			grid1$strand <- temp_tss$strand
			grid <- unique(grid1)
			temp <- grid
		} else if (grid_mode=="merged") {
			message("Grid mode = merged")
			grid <- as.data.frame(cbind(interv1[,1:3],interv2[,1:3]))
			message(nrow(grid))
			#  grid <- grid[,-1]
			grid[,7]=1
			for (i in 1:nrow(grid)){
			 if (grid[i,2] > grid[i,5]) {
				 grid[i,] <- grid[i,c(4,5,6,1,2,3)]
				 grid[i,7] <- 2
				}
			} 
			colnames(grid) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2","type")		
			### Assumes that interv1 is Tss #####
			temp_tss <- interv1
			#########
			colnames(temp_tss) <- c('chrom','start','end') 
			temp_tss <- gintervals.neighbors(temp_tss,tss)
 			grid$strand <- temp_tss$strand
			grid <- grid[((grid$start2 - grid$start1) >= min_dist) & ((grid$start2 - grid$start1) <= max_dist),]
			message(nrow(grid))
			message(paste0('inverted ',nrow(grid[grid$type==2,]),' number of rows'))
			temp <- grid
		} else if (grid_mode=="loops") {
		  message("Grid mode = loops")
		  grid <- as.data.frame(cbind(interv1[,1:3],interv2[,1:3]))
		  message(nrow(grid))
		  #  grid <- grid[,-1]
		  grid[,7]=1
		  for (i in 1:nrow(grid)){
		    if (grid[i,2] > grid[i,5]) {
		      grid[i,] <- grid[i,c(4,5,6,1,2,3)]
		      grid[i,7] <- 2
		    }
		  } 
		  colnames(grid) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2","type")
		  grid <- grid[((grid$start2 - grid$start1) >= min_dist) & ((grid$start2 - grid$start1) <= max_dist),]
		  message(nrow(grid))
		  message(paste0('inverted ',nrow(grid[grid$type==2,]),' number of rows'))
		  grid <- unique(grid)
		}

		message(paste("starting with ", nrow(grid),"intervals"))
		if (min_dist<2e6){
			if (domains=='intra'){
				dom = gintervals.2d(chroms1=dom$chrom, starts1=dom$start, ends1=dom$end,chroms2=dom$chrom, starts2=dom$start, ends2=dom$end)
				grid <- gintervals.intersect(grid,dom)
			} else if (domains=='inter'){
				inter_2d = construct.grid2(dom,dom,min_dist,max_dist)			
				grid <- gintervals.intersect(grid,inter_2d)
	#			grid <- anti_join(grid,gintervals.intersect(grid,dom))
			}
			message(paste("found ", nrow(grid), " overlapping intervals"))
		} else {
			inter_2d = construct.grid2(dom,dom,min_dist,max_dist)			
			grid <- gintervals.intersect(grid,inter_2d)			
			print(head(grid))
		}
  	
		if (filter_value>0){
			message("screening for filter")
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
		#	save(grid,file=grid_f)
		}
		if (grid_mode!='1D'&grid_mode!='loops'){
			grid <- grid[,1:6]
			grid$id <- paste0(grid$chrom1,grid$start1,grid$end1,grid$start2,grid$end2)
			temp$id <- paste0(temp$chrom1,temp$start1,temp$end1,temp$start2,temp$end2)
			temp <- merge(grid,temp,by='id',all=FALSE)
			grid <- temp[,c(2:7,14,15)]
			colnames(grid) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2","type","strand")
		#	save(grid,file=grid_f)
		} 

	} else {
		#grid <- get(load(grid_f))
	}
  ######### End of modification
  grid <- grid[grid$chrom1!='chrM',]
	levels(grid$chrom1) <- levels(gintervals.all()$chrom)
	levels(grid$chrom2) <- levels(gintervals.all()$chrom)
	message(paste("working on",domains," intervals"))
	message(paste("found ", nrow(grid), " regions after score filter"))
	if(!use_score){
	  interv_set_out = gsub("+", "", paste(track,paste(rev(-(band)), collapse="_"), sep="_"), fixed=T)
	  if (!gintervals.exists(interv_set_out)) {
	    giterator.intervals(track, band=band, intervals.set.out=interv_set_out)
	  }
	  
	  a = gintervals.neighbors(interv_set_out, grid, 
	                           mindist1=min(expand1), maxdist1=max(expand1), 
	                           mindist2=min(expand2), maxdist2=max(expand2))
	  colnames(a)[7:(ncol(a)-2)] = paste0("grid.", colnames(grid)) 
	  message(' N of obs contacts in intervals:',nrow(a))
	  if (grid_mode=="1D"|grid_mode=="loops"){
	    o = table(cut(a$start1 - a$grid.start1, breaks=expand1, include.lowest=TRUE), 
	              cut(a$start2 - a$grid.start2, breaks=expand2, include.lowest=TRUE))
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
	        if (type==2&strand==(-1)) { o1 <- t(rotate(rotate(o1)))} 
	        else if (type==2&strand==(1)) { o1 <- apply(o1,1,rev)} 
	        else if (type==1&strand==(-1)) { o1 <- t(apply(o1,1,rev))}
	        else if (type==1&strand==1) { next}		
	        colnames(o1) <- colnames(o)
	        rownames(o1) <- rownames(o)
	        o <- o + o1
	      }
	    }
	  }
	  
	  interv_set_out = gsub("+", "", paste(shuffled_track,paste(rev(-(band)), collapse="_"), sep="_"), fixed=T)
	  if (!gintervals.exists(interv_set_out)) {
	    giterator.intervals(shuffled_track, band=band, intervals.set.out=interv_set_out)
	  }
	  a = gintervals.neighbors(interv_set_out, grid, mindist1=min(expand1), maxdist1=max(expand1), 
	                           mindist2=min(expand2), maxdist2=max(expand2))
	  colnames(a)[7:(ncol(a)-2)] = paste0("grid.", colnames(grid))
	  message(' N of exp contacts in intervals:',nrow(a))
	  if (grid_mode=="1D"|grid_mode=="loops"){
	    e = table(cut(a$start1 - a$grid.start1, breaks=expand1, include.lowest=TRUE), 
	              cut(a$start2 - a$grid.start2, breaks=expand2, include.lowest=TRUE))
	  } else {
	    a1 <- a[a$grid.type==1&a$grid.strand==1,]
	    e = table(cut(a1$start1 - a1$grid.start1, breaks=expand1, include.lowest=TRUE), cut(a1$start2 - a1$grid.start2, breaks=expand2, include.lowest=TRUE))
	    message(' N of exp contacts type1 strand1:',nrow(a1))
	    for (type in c(1,2)){
	      for (strand in c(1,-1)){
	        a1 <- a[a$grid.type==type&a$grid.strand==strand,]
	        e1 = table(cut(a1$start1 - a1$grid.start1, breaks=expand1, include.lowest=TRUE), cut(a1$start2 - a1$grid.start2, breaks=expand2, include.lowest=TRUE))
	        if (type==2&strand==(-1)) { e1 <- t(rotate(rotate(e1)))} 
	        else if (type==2&strand==(1)) { e1 <- apply(e1,1,rev)} 
	        else if (type==1&strand==(-1)) { e1 <- t(apply(e1,1,rev))}
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
	  return(list(obs=obs, exp=exp, total_obs=total_obs, total_exp=total_exp, grid_size=nrow(grid)))
	} else {
	  library(doParallel)
	  registerDoParallel(8)
	  res <- foreach(i=1:nrow(grid))%dopar%{
	    interv <- grid[i,]
	    o = matrix(0, nrow=length(expand1), ncol=length(expand2))
	    colnames(o) <- expand1
	    row.names(o) <- expand2
	    for (s in 1:nrow(o)){
	      for (k in 1:ncol(o)){
	        inter <- grid[i,]
	        inter[,c(2,3)] <- inter[,c(2,3)] + as.numeric(row.names(o)[s])
	        inter[,c(5,6)] <- inter[,c(5,6)] + as.numeric(colnames(o)[k])
	        inter[,7]<- s
	        inter[,8] <- k
	        interv <- rbind(interv,inter)
	      }
	    }
	    interv <- interv[-1,]
	    df <- gextract(track,intervals = interv,iterator = interv)
	    for (s in 1:nrow(o)){
	      for (k in 1:ncol(o)){
	        o[s,k] <- df[which(interv[,7]==s&interv[,8]==k),7]
	      }
	    }
	    type <- grid[i,'type']
	    strand <- grid[i,'strand']
	    if (type==2&strand==(-1)) { o <- t(rotate(rotate(o)))} 
	    else if (type==2&strand==(1)) { o <- apply(o,1,rev)} 
	    else if (type==1&strand==(-1)) { o <- t(apply(o,1,rev))}
	    return(o)
	  }
	  arr <- array( unlist(res) , c(length(expand),length(expand),length(res)) )
	  df <- apply( arr , 1:2 , mean,na.rm=T )
	  return(list(obs=df, exp=NA, total_obs=0, total_exp=0, grid_size=nrow(grid)))
	}
}

main_f=args[1]
track = args[2]
shuffled_track = paste0(track, "_shuffle")
#score_track = paste0(track, "_score")
score_track = args[3] 
filter_v = as.numeric(args[4])
i1 = args[5]
i2 = args[6]
data_file = args[7]
expand_range = as.numeric(args[8])
expand_res = as.numeric(args[9])
domains_f = args[10]
grid_mode = args[11]
interval_name = args[12]
use_score=as.logical(args[13])
long_range=FALSE


source(paste0(main_f,'config.R'))
message(paste("expand_range=", expand_range))
tss <- gintervals.load(tss_f)
expand = seq(0-expand_range, expand_range, by=expand_res)
if (grepl('narrowPeak',i1)|grepl('\\.bed',i1)) {
	peaks <- read.table(i1)
	interv1 <- gintervals(chroms=peaks$V1,starts=peaks$V2,ends=peaks$V3)
	interv1 <- interv1[interv1$chrom %in% gintervals.all()$chrom,]
	interv1 <- interv1[(as.character(interv1$chrom)!='chrM'|as.character(interv1$chrom)!='chrY'),]
} else {
	interv1 = get(load(file=i1))
}
if (grid_mode!='merged'&grid_mode!='loops'){
	interv1 <- intervals.expand(interv1, expand_res/2)
	interv1 <- gintervals.canonic(interv1)
	interv1 <- unique(intervals.centers(interv1))
} else {
	interv1 <- intervals.centers(interv1)
}

if (grepl('narrowPeak',i2)|grepl('\\.bed',i2)) {
	peaks <- read.table(i2)
	interv2 <- gintervals(chroms=peaks$V1,starts=peaks$V2,ends=peaks$V3)
	interv2 <- interv2[interv2$chrom %in% gintervals.all()$chrom,]
	interv2 <- interv2[(as.character(interv2$chrom)!='chrM'|as.character(interv2$chrom)!='chrY'),]
} else {
	interv2 = get(load(file=i2))
}
if (grid_mode!='merged'&grid_mode!='loops'){
	interv2 <- intervals.expand(interv2, expand_res/2)
	interv2 <- gintervals.canonic(interv2)
	interv2 <- unique(intervals.centers(interv2))
} else {
	interv2 <- intervals.centers(interv2)
}


if (!grepl('dm',genome)) {bands_matrix <- as.data.frame(cbind(c(5e4,5e4,2e5,2e6),c(2e5,2e6,2e6,10e6)))} else {bands_matrix <- as.data.frame(cbind(c(5e4,5e4,2e6),c(2e5,2e6,5e6)))}
if (long_range){bands_matrix <- rbind(bands_matrix,c(1e7,5e7))}
if(grepl('exonQ',interval_name)){bands_matrix <- as.data.frame(cbind(c(2e5,2e6,10e6),c(2e6,10e6,50e6)))}
if(grepl('Exon',interval_name)){bands_matrix <- as.data.frame(cbind(c(2e6,10e6),c(10e6,50e6)))}
if(expand_range==100000){bands_matrix <- as.data.frame(cbind(c(10e4,2e6),c(2e6,10e6)))}
if(expand_range==80000){bands_matrix <- as.data.frame(cbind(c(8e4,2e6),c(2e6,10e6)))}
if(expand_range==60000){bands_matrix <- as.data.frame(cbind(c(6e4,2e6),c(2e6,10e6)))}

#bands_matrix <- as.data.frame(cbind(c(10e6),c(50e6)))          # Temp testing

scoreTracks <- paste0("hic.",cells,".score_",score_track)                       #
if (filter_v!=0){
	 v_scoreTracks <- paste("v_", scoreTracks, sep="")
	   for (i in 1:length(v_scoreTracks)) {
		gvtrack.create(v_scoreTracks[i], scoreTracks[i], 'max')
		gvtrack.iterator.2d(v_scoreTracks[i], sshift1=min(expand), eshift1=max(expand), sshift2=min(expand), eshift2=max(expand))
	   }
} else { v_scoreTracks <- '' }

if(use_score){
  gvtrack.create(paste0('v_',track), track, 'avg')
  gvtrack.iterator.2d(paste0('v_',track), sshift1=-expand_res, eshift1=expand_res, sshift2=-expand_res, eshift2=expand_res)
  track=paste0('v_',track)
}

if (domains_f=='NULL'){
  domains <- gintervals.all()
} else {
  domains <- gintervals.load(domains_f) 
}

message(paste0("relying on ",nrow(domains)," domains"))
#dom_2d = gintervals.2d(chroms1=domains$chrom, starts1=domains$start, ends1=domains$end,chroms2=domains$chrom, starts2=domains$start, ends2=domains$end)
dom_2d <- unique(domains)

print(bands_matrix)
obs = list()
exp = list()
grid_size = c()
total_obs = c()
total_exp = c()
message('interv1:',nrow(interv1))
message('interv2:',nrow(interv2))
for (s in 1:nrow(bands_matrix)) {
	if (bands_matrix[s,1]<2e6){
		if (grid_mode!='merged'&grid_mode!='loops'){
			grid = norm_contact_grid_by_shuffled_filter_domains(track, shuffled_track, v_scoreTracks, filter_v, unique(interv1[interv1$chrom != "chrM",1:3]), unique(interv2[interv2$chrom != "chrM",1:3]),expand, expand, min_dist=bands_matrix[s,1], max_dist=bands_matrix[s,2], regularization=0,dom=dom_2d,domains="intra",grid_mode=grid_mode,interval_name=interval_name)
		} else {
			grid = norm_contact_grid_by_shuffled_filter_domains(track, shuffled_track, v_scoreTracks, filter_v, interv1, interv2, expand, expand,  min_dist=bands_matrix[s,1], max_dist=bands_matrix[s,2], regularization=0,dom=dom_2d,domains="intra",grid_mode=grid_mode,interval_name=interval_name)
		}
			obs[[length(obs)+1]] = grid$obs
			exp[[length(exp)+1]] = grid$exp
			grid_size = c(grid_size, grid$grid_size)
			total_obs = c(total_obs, grid$total_obs)
			total_exp = c(total_exp, grid$total_exp)
			message(paste("grid size=", grid$grid_size, ", max obs=", max(log10(grid$obs),na.rm=TRUE), ", max obs/exp=", max(log2(grid$obs/grid$exp),na.rm=TRUE)))
	}
	if (bands_matrix[s,2]>2e5){
		if ((grid_mode!='merged'&grid_mode!='loops')&(domains_f!='NULL')){
			grid = norm_contact_grid_by_shuffled_filter_domains(track, shuffled_track, v_scoreTracks, filter_v, unique(interv1[interv1$chrom != "chrM",1:3]), unique(interv2[interv2$chrom != "chrM",1:3]),expand, expand, min_dist=bands_matrix[s,1], max_dist=bands_matrix[s,2], regularization=0,dom=dom_2d,domains="inter",grid_mode=grid_mode,interval_name=interval_name)
		} else {
		  next
		 # grid = norm_contact_grid_by_shuffled_filter_domains(track, shuffled_track, v_scoreTracks, filter_v, interv1, interv2,expand, expand, min_dist=bands_matrix[s,1], max_dist=bands_matrix[s,2], regularization=0,dom=dom_2d,domains="inter",grid_mode=grid_mode,interval_name=interval_name)
		}
		obs[[length(obs)+1]] = grid$obs
		exp[[length(exp)+1]] = grid$exp
		grid_size = c(grid_size, grid$grid_size)
		total_obs = c(total_obs, grid$total_obs)
		total_exp = c(total_exp, grid$total_exp)
		message(paste("grid size=", grid$grid_size, ", max obs=", max(log10(grid$obs),na.rm=TRUE), ", max obs/exp=", max(log2(grid$obs/grid$exp),na.rm=TRUE)))
	}	
}
grid = list(obs = obs, exp=exp, grid_size=grid_size, total_obs=total_obs, total_exp = total_exp)
save(grid, file=data_file)




