library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(tidyfast)
library(doParallel)
library(fst)

intervals.centers <- function(inv){
  inv[,2:3]<-floor((inv[,2]+inv[,3])/2)
  inv[,3]<-inv[,3]+1
  return(inv)
}

#returns the nearest occurence of x in vec   credit(https://stackoverflow.com/questions/10160400/r-find-nearest-index)
nearest.vec <- function(x, vec)
{
  smallCandidate <- findInterval(x, vec, all.inside=TRUE)
  largeCandidate <- smallCandidate + 1
  #nudge is TRUE if large candidate is nearer, FALSE otherwise
  nudge <- 2 * x > vec[smallCandidate] + vec[largeCandidate]
  return(smallCandidate + nudge)
}

##################
args = commandArgs(trailingOnly=TRUE)

binSize <- as.numeric(args[1])
interv1_f <- args[2]
path <- args[3]
file_id <- args[4]
out_f <- args[5]
n_cores <- as.numeric(args[6])

registerDoParallel(cores=n_cores)

chr_sizes <- read.table('/home/hpc/bonev/annotations/mm10/mm10.chrom.sizes')
files_f <- dir(path, pattern=file_id,full.names = T)
##################


interv1 <- read.delim2(interv1_f, header = FALSE)[,1:3]
colnames(interv1) = c('chrom', 'start', 'end')
interv1 <- intervals.centers(interv1)

colnames(chr_sizes) <- c('chrom','end')
chr_sizes$start <- 1
chr_sizes <- chr_sizes[,c(1,3,2)]
chrs <- unique(as.character(chr_sizes$chrom))
chrs <- chrs[(chrs!='chrM')&(chrs!='chrY')]

res <- list()

Sys.time()
res <- foreach(chr=chrs, .export=c('chr_sizes','interv1','interv2'),
               .packages=c('dplyr','tidyr','tibble','stringr','tidyfast')) %dopar% {
                 temp_df1 <- list()
                 temp_df2 <- list()
                 
                 filter1_df <- interv1[interv1$chrom==chr,]     #Subset filter intervals1 by chr
                 filter1 <- unique(c(chr_sizes[chr_sizes$chrom==chr,'start'],sort(filter1_df$start),chr_sizes[chr_sizes$chrom==chr,'end']))          #Create sorted vector with interval starts, beginning and ending with chromosome coordinates. This is needed for the binary search later.  
                 
                 for (read_id in 1:2) {                  #Go through the CpG selected reads and retain only reads matching filter
                   reads <- read.fst(paste0(path,'/CpG/',file_id,'_CpG_',chr,'_',read_id,'.fst'))
                   temp_range <- str_split(reads$Range,pattern = ',',simplify = TRUE)
                   mode(temp_range) = "numeric"
                   
                   a_idx <- t(apply(temp_range,1,function(x){
                     #interv_bin <- filter1[findInterval(x, filter1)]             #execute binary search to identify the closest filter interval for each methylation call coordinate
                     interv_bin <-filter1[nearest.vec(x,filter1)]                #execute binary search to identify the closest filter interval for each methylation call coordinate
                     nearest_dist <- abs(x-interv_bin)                           #calculate the absolute distance between each methylation call and the nearest filter interval
                     interv_bin[nearest_dist>binSize] <- NA                      #transform all calls outside of the user defined window to NAs
                     if (min(nearest_dist,na.rm=T)<=binSize&!(is.infinite(min(nearest_dist,na.rm=T)))){      #Return a row (read) if at least one call was within the specified window. Check for non-ifinite rows otherwise min causes bugs when only NAs.  
                       return(interv_bin)
                     } else {
                       rep(NA,length(interv_bin))
                     }
                   }))
                   
                   
                   reg_a <- which(rowSums(a_idx,na.rm=T)>=1)               #find rows which have TRUE in them
                   interv_id <- a_idx[reg_a,]
                   count <- data.frame(Mcalls=rowSums(!is.na(a_idx[reg_a,]), na.rm=T))
                   temp_call <- str_split(reads$Call[reg_a],pattern = ',',simplify = TRUE)     #split actual methylation calls in individual columns for next step
                   mode(temp_call) = "numeric"
                   df <- cbind(as.matrix(reads[reg_a,c('V1','rid', 'A')]),count,round(rowMeans(temp_call,na.rm=T),2),temp_range[reg_a,],temp_call,interv_id)    #Create a summary matrix with read names, ids, .. as well as the filter interval ids(starts) and call coordinates per row 
                   colnames(df) <- c('read_name','read_id','chr','n_calls','perc_meth',paste0('call_coord',1:ncol(temp_range)),paste0('call',1:ncol(temp_call)),paste0('interv_coord',1:ncol(interv_id)))
                   temp_df1[[read_id]] <- df[order(df$read_name),]   
                 }
                 ### Switch to GpC and select matching reads
                 for (read_id in c(1,2)) {
                   reads <- read.fst(paste0(path,'/GpC/',file_id,'_GpC_',chr,'_',read_id,'.fst'))
                   df <- temp_df1[[read_id]]
                   
                   reg_a <- match(df$read_name,reads$V1)             #Match read IDs
                   reg_a <- reg_a[!is.na(reg_a)]
                   reads <- reads[reg_a,]
                   
                   temp_range <- str_split(reads$Range,pattern = ',',simplify = TRUE)
                   mode(temp_range) = "numeric"
                   
                   interv_id <- temp_range
                   count <- data.frame(Mcalls=rowSums(!is.na(temp_range), na.rm=T))
                   temp_call <- str_split(reads$Call,pattern = ',',simplify = TRUE)     #split actual methylation calls in individual columns for next step
                   mode(temp_call) = "numeric"
                   df <- cbind(as.matrix(reads[,c('V1','rid', 'A')]),count,round(rowMeans(temp_call,na.rm=T),2),temp_range,temp_call,interv_id)    #Create a summary matrix with read names, ids, .. as well as the filter interval ids(starts) and call coordinates per row 
                   colnames(df) <- c('read_name','read_id','chr','n_calls','perc_meth',paste0('call_coord',1:ncol(temp_range)),paste0('call',1:ncol(temp_call)),paste0('interv_coord',1:ncol(interv_id)))
                   temp_df2[[read_id]] <- df[order(df$read_name),]   
                 }
                 
                 gc()
                 #Now it is time to match read ids and ensure that only pairs overlapping both filter anchors simultaneously are retained
                 
                 df1 <- merge(temp_df1[[1]],temp_df2[[1]],by='read_name')    #First match read ids for when read1 overlaps left filter anchor and read2 overlaps right filter anchors
                 df2 <- merge(temp_df1[[2]],temp_df2[[2]],by='read_name')    #Analogous for read2 overlapping right anchor
                 df <- rbindlist(list(df1,df2), fill=TRUE)
                 return(df)
               }
Sys.time()
res <- rbindlist(res, fill=TRUE)
res <- Filter(function(x)!all(is.na(x)), res)     #filter all na columns

write.fst(res,out_f)


