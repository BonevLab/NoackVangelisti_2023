#!/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

main_f <- as.character(args[5])
source(paste0(main_f,'config.R'))

ins_window <- as.numeric(args[1])
k=as.numeric(args[2])
cell=args[3]
near_cis <- as.numeric(args[4])
eigenBin <- as.numeric(args[6])
expand_cis <- as.numeric(args[7])

work_dir <- paste0(main_f,'temp/')

library(shaman)

options("shaman.debug"=FALSE)
if (sunGrid) {options("shaman.sge_support"=1)} else {options("shaman.mc_support"=1)}
options("shaman.sge_flags"="-l mem=60G -l h_vmem=60G")

intervals.expand <- function(inv,expansion=100){
  inv[,2]<-inv[,2]-expansion
  inv[,3]<-inv[,3]+expansion
  return(gintervals.force_range(inv))
}

calculate_eigen <-function(trackdb=trackdb,tracks=all_tracks,cells=cells,chrs=gsub('chr','',gintervals.all()[,1]),path_f=paste0(main_f,'analysis/compartments/'),bin=1e5,gen=genome,cworld=cworld_path,chrom_f=chrom_sizes_f){
  for (t in 1:length(tracks)){
    gvtrack.create(paste0('v_',tracks[t]), tracks[t], "area")
  }
  for (cell in cells){
    for (chr in chrs){
      chr_int <- gintervals.2d(chr)
      ## Drosophila chromosomes are strange and eigenvector doesn't work well with the centromere regions. We fix it by avoiding 4MB from the centromere for chromosomes 2 and 3
      if(chr=='2R'|chr=='3R') {
        chr_int$start1 <- round_any(chr_int$start1+2e6,bin)
        chr_int$start2 <- round_any(chr_int$start2+2e6,bin)
      } else if (chr=='2L'|chr=='3L'){
        chr_int$end1 <- round_any(chr_int$end1-2e6,bin)
        chr_int$end2 <- round_any(chr_int$end2-2e6,bin)
      }
      ########################
      chr_int <- gintervals.force_range(chr_int)		
      bins_2d <- giterator.intervals(intervals=chr_int, iterator=c(bin,bin))
      message(cell,'_chr',chr)
      if(!file.exists(paste0(path_f,cell,'/',bin/1000,'kb/chr',chr,'_cis.txt'))|(file.info(paste0(path_f,cell,'/',bin/1000,'kb/chr',chr,'_cis.txt'))$size==0)){
        dir.create(paste0(path_f,cell,'/',bin/1000,'kb'),showWarnings = FALSE, recursive=T)
        df <- gextract(paste0('v_',tracks[grep(cell,tracks)],collapse='+'),intervals=bins_2d,iterator=bins_2d)
        df <- df[,-ncol(df)]
        df$bin1 <- paste0(df[,2]/bin,'|',gen,'|chr',chr,':',df[,2],'-',df[,3])
        df$bin2 <- paste0(df[,5]/bin,'|',gen,'|chr',chr,':',df[,5],'-',df[,6])
        df <- df[,7:ncol(df)]
        colnames(df)[1] <- 'count'
        df_cast <- dcast(df,factor(bin1,levels=unique(bin1))~factor(bin2,levels=unique(bin2)),value.var='count',)
        colnames(df_cast)[1] <- 'bin'
        write.table(df_cast,paste0(path_f,cell,'/',bin/1000,'kb/chr',chr,'_cis.txt'),col.names=T,row.names=F,quote=F,sep='\t')
        command_f <- paste0('cd ',paste0(path_f,cell,'/',bin/1000,'kb/ ; '),'perl -I ',cworld,'/lib/ ',cworld,'/scripts/perl/matrix2compartment.pl -i ',paste0(path_f,cell,'/',bin/1000,'kb/chr',chr,'_cis.txt'),' --et --minDist 2000000 -o ',path_f,cell,'/',bin/1000,'kb/chr',chr,'_eigen')
        system(command_f,wait = TRUE)
      } else if (!file.exists(paste0(path_f,cell,'/',bin/1000,'kb/chr',chr,'_eigen.zScore.eigen1.bedGraph'))){
        command_f <- paste0('cd ',paste0(path_f,cell,'/',bin/1000,'kb/ ; '),'perl -I ',cworld,'/lib/ ',cworld,'/scripts/perl/matrix2compartment.pl -i ',paste0(path_f,cell,'/',bin/1000,'kb/chr',chr,'_cis.txt'),' --et --minDist 2000000 -o ',path_f,cell,'/',bin/1000,'kb/chr',chr,'_eigen')
        system(command_f,wait = TRUE)
      } else {
        message('Eigen track exist. Skipping...')
      }
    }
    command_f <- paste0('cd ',paste0(path_f,cell,'/',bin/1000,'kb/ ;'),'cat *eigen1.bedGraph > combined.bedGraph')
    system(command_f,wait = TRUE)
    command_f <- paste0('cd ',paste0(path_f,cell,'/',bin/1000,'kb/ ; grep -v bedGraph combined.bedGraph > combined_pol.bedGraph; mv combined_pol.bedGraph combined.bedGraph; sort -k1,1 -k2,2n combined.bedGraph > combined_pol.bedGraph; mv combined_pol.bedGraph combined.bedGraph; '),'bedGraphToBigWig combined.bedGraph ',chrom_f,' ',paste0(path_f,cell,'.bw'))
    system(command_f,wait = TRUE)
    gtrack.import(track=paste0('eigen.',cell,'_',bin/1000,'kb'),description=paste0('eigen vector based on ',bin/1000,'kb'),file=paste0(path_f,cell,'/',bin/1000,'kb/combined.bedGraph'),binsize=bin)
  }
}

gtrack.2d.gen_insu_prof = function(track_nm, scale, res, min_diag_d=1000)
{
  #extract - using an iter on  	
  iter_1d = giterator.intervals(intervals=ALLGENOME[[1]], iterator=res)
  iter_2d = gintervals.2d(chroms1 = iter_1d$chrom, 
                          starts1=iter_1d$start, 
                          ends1=iter_1d$end, 
                          chroms2 = iter_1d$chrom, 
                          starts2=iter_1d$start, 
                          ends2=iter_1d$end)
  if(length(gvtrack.ls("obs_big")) == 1) {
    gvtrack.rm("obs_big")
  }
  if(length(gvtrack.ls("obs_ins")) == 1) {
    gvtrack.rm("obs_ins")
  }
  gvtrack.create("obs_big", track_nm, "weighted.sum")
  gvtrack.create("obs_ins", track_nm, "weighted.sum")
  gvtrack.iterator.2d("obs_big", 
                      eshift1=scale, sshift1=-scale, 
                      eshift2=scale, sshift2=-scale)
  gvtrack.iterator.2d("obs_ins", 
                      eshift1=0, sshift1=-scale, 
                      eshift2=scale, sshift2=0)
  message("will iter on ", dim(iter_2d)[1])
  ins = gextract("obs_big", "obs_ins", 
                 gintervals.2d.all(),
                 iterator=iter_2d, band=c(-scale*2,0))
  ins_diag = gextract("obs_big", "obs_ins", 
                      gintervals.2d.all(),
                      iterator=iter_2d, band=c(-min_diag_d,0))
  ###fixing NA bug
  ins[is.na(ins)] = 0
  ins_diag[is.na(ins_diag)] = 0
  ##
  ins$obs_big = ins$obs_big - ins_diag$obs_big
  ins$obs_ins = ins$obs_ins - ins_diag$obs_ins
  message("will retrun ins with ", dim(ins)[1], " rows")
  ##### Fix for edge of repeats ####
  ins$obs_big[ins$obs_ins <= 10] = NA
  ins$obs_ins[ins$obs_ins <= 10] = NA	     
  #### End of fix #####
  
  return(ins)
}

gtrack.2d.gen_joint_insu_prof = function(track_nms, scale, res, min_diag_d = 1000)
{
  ins = gtrack.2d.gen_insu_prof(track_nms[1], scale, res, min_diag_d)
  for (track in track_nms[-1]) {
    t_ins =  gtrack.2d.gen_insu_prof(track, scale, res, min_diag_d)
    ins$obs_big = rowSums(cbind(ins$obs_big, t_ins$obs_big), na.rm=T)
    ins$obs_ins = rowSums(cbind(ins$obs_ins, t_ins$obs_ins), na.rm=T)
  }
  ##### Fix for edge of repeats ####
  ins$obs_big[ins$obs_ins <= 10] = NA
  ins$obs_ins[ins$obs_ins <= 10] = NA	     
  #### End of fix #####
  return(ins)
}


gtrack.2d.gen_insu_track = function(track_nm, scale, res, min_diag_d=1000, new_track, description="")
{
  k_reg = 10
  if (length(track_nm) > 1) {
    prof = gtrack.2d.gen_joint_insu_prof(track_nm, scale, res, min_diag_d)
  } else {
    print(paste0("single source track: ",track_nm))
    prof = gtrack.2d.gen_insu_prof(track_nm, scale, res, min_diag_d)
    prof$obs_big[prof$obs_ins <= 10] = NA
    prof$obs_ins[prof$obs_ins <= 10] = NA	 
  }
  message("names ", paste(names(prof), collapse=","))
  names(prof)[1] = "chrom"
  names(prof)[2] = "start"
  names(prof)[3] = "end"
  
  ins <- as.data.frame(cbind(prof[,c(1,2,3)], log2(prof$obs_ins/(prof$obs_big+k_reg))))
  ins_na <- ins[is.na(ins[,4]),]
  inter_na <- gintervals.canonic(ins_na)
  inter_na <- intervals.expand(inter_na,expansion=scale)
  ins <- gintervals.neighbors(ins,inter_na)
  ins[ins$dist==0,4] <- NaN
  gtrack.create_sparse(track=new_track, ins[,1:3], ins[,4], description=description)
}


all_tracks <- all_tracks[grep(cell,all_tracks)] 

#####################
## Shuffle tracks ###
#####################

shuffled_tracks <- paste0(all_tracks,'_shuffle')
for (track in all_tracks[!gtrack.exists(shuffled_tracks)]){
  message('Shuffling: working on Track: ',track)
  shaman_shuffle_hic_track(track_db=trackdb, obs_track_nm=track,work_dir=work_dir,max_jobs=nrow(gintervals.all()))
}

#####################
##### Insulation ####
#####################

ins_rep_tracks <- paste0(all_tracks,'_ins_',ins_window)
for (track in all_tracks[!gtrack.exists(ins_rep_tracks)]){
  gtrack.2d.gen_insu_track(track, ins_window*1000, 1000, new_track=paste0(track,'_ins_',ins_window), description=paste0(ins_window,"k scale at 1kbp resolution"), min_diag_d=10000)
}

ins_track <- paste0('hic.',cell,'.ins_',ins_window)
if (!gtrack.exists(ins_track)){
  message('Insulation: working on Track: ',cell)
  gtrack.2d.gen_insu_track(all_tracks, ins_window*1000, 1000, new_track=ins_track, description=paste0(ins_window,"k scale at 1kbp resolution"), min_diag_d=10000)
}

#####################
##### Eigen ####
#####################

# gdir.create('eigen', showWarnings=F)
# eigen_rep_tracks <- paste0(gsub(paste0('hic.',cell,'.'),'eigen.',all_tracks),'_',eigenBin,'kb')
# for (track in all_tracks[!gtrack.exists(eigen_rep_tracks)]){
# 	calculate_eigen(tracks=track,cells=cell,bin=eigenBin*1000)
# }

eigen_track <- paste0('eigen.',cell,'_',eigenBin,'kb')
if (!gtrack.exists(eigen_track)){
  calculate_eigen(tracks=all_tracks,cells=cell,bin=eigenBin*1000)
}

#####################
##### HiC Scores ####
#####################	
gdb.reload()
score_track <- paste0('hic.',cell,'.score_k',k)
if (!gtrack.exists(score_track)){
  message('Calculating scores: working on Track: ',cell)
  shaman_score_hic_track(track_db=trackdb, work_dir=work_dir, score_track_nm=score_track, obs_track_nms=all_tracks,near_cis=near_cis, expand=expand_cis, k=k, max_jobs=16)
}






