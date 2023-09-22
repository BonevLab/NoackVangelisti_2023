library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(misha)


#####################
#Variables to define 
#####################
# main_f <- '/home/fnoack/projects/3D_RAMseq_250k/rep1/'
# bedgraph_file <- '3DRAM_ES_250k_rep1_R1.fastq.gz_trimmed.fastq_CpG_cov5x.bedGraph'
# rf <- c('CpG_CTOB_3DRAM_ES_250k_rep1_R1.fastq.gz_trimmed.fastq_all_merged_3split_dedupped.txt', 
#                 'CpG_CTOT_3DRAM_ES_250k_rep1_R1.fastq.gz_trimmed.fastq_all_merged_3split_dedupped.txt',
#                 'CpG_OB_3DRAM_ES_250k_rep1_R1.fastq.gz_trimmed.fastq_all_merged_3split_dedupped.txt',
#                 'CpG_OT_3DRAM_ES_250k_rep1_R1.fastq.gz_trimmed.fastq_all_merged_3split_dedupped.txt')



#gene <- c('Pax6')
#tss_f <- tss_f[tss_f$geneName %in% gene,]
# 
# #####################
# #1 - Filtering
# #####################
# 
# bG_file <- fread(paste0(main_f, bedgraph_file))
# bG_file$id1 <- paste(bG_file$V1, bG_file$V2, sep="_")
# bG_file$id2 <- paste(bG_file$V1, bG_file$V3, sep="_")
# 
# c = 1 
# for (i in rf) {
#   #myenv <- new.env()
#   x <- fread(paste0(main_f, i))
#   assign(paste0('read_file', c, sep =""), x)
#   c = c + 1
# }
# l <- ls(pattern='read_file')
# g = 1
# 
# for (i in l) {
#   z <- get(i)
#   colnames(z) = c('V1', 'V2', 'V3', 'V4', 'V5')
#   z$id <- paste(z$V3, z$V4, sep="_")
#   assign(paste0('read_file', g, sep =""), z)
#   g = g+1
#   
# }
# 
# R1 <- rbind(read_file3, read_file4)
# R2 <- rbind(read_file1, read_file2)
# 
# filtered_R1 <- R1 %>% filter(id %in% bG_file$id1 | id %in% bG_file$id2)
# filtered_R2 <- R2 %>% filter(id %in% bG_file$id1 | id %in% bG_file$id2)
# read_file1 <- readRDS('/home/faye/scripts/rf1.RDS')
# read_file2 <- readRDS('/home/faye/scripts/rf2.RDS')
# read_file3 <- readRDS('/home/faye/scripts/rf3.RDS')
# read_file4 <- readRDS('/home/faye/scripts/rf4.RDS')
# 
# 
# ####################
# #1a - checking read similarity
# ####################
# readN1 <- read_file1$V1
# read_name1 <- unique(sapply(strsplit(readN1, '\\_'), "[", 1))
# readN2 <- read_file2$V1
# read_name2 <- unique(sapply(strsplit(readN2, '\\_'), "[", 1))
# readN3 <- read_file3$V1
# read_name3 <- unique(sapply(strsplit(readN3, '\\_'), "[", 1))
# readN4 <- read_file4$V1
# read_name4 <- unique(sapply(strsplit(readN4, '\\_'), "[", 1))
# Rtwo <- unique(c(read_name1, read_name2))
# Rone <- unique(c(read_name3, read_name4))
# Rtwo1 <- c(read_name1, read_name2)
# Rone1 <- c(read_name3, read_name4)
# saveRDS(Rtwo, '/home/faye/scripts/Rtwo.RDS')
# saveRDS(Rone, '/home/faye/scripts/Rone.RDS')
# count(Rtwo[Rtwo %in% Rone])
# length(intersect(Rtwo1, Rone1))
# 
# 
# ####################
# #2 - Per read info
# ####################
# filtered_CpG <- readRDS('/home/faye/scripts/3DRAM/filtered_CpG.RDS')
files <- readRDS('/home/faye/scripts/filteredR1.RDS')

# #filtered_R2 <- readRDS('/home/faye/scripts/filteredR2.RDS')
# files <- bind_rows(filtered_CpG)
# saveRDS(files, '/home/faye/scripts/3DRAM/Rfiltered_CpG.RDS')

files <- readRDS('/home/faye/scripts/3DRAM/Rfiltered_CpG.RDS')
write.csv2(reads,'/home/faye/scripts/3DRAM/reads.csv')
reads <- readRDS('/home/faye/scripts/3DRAM/reads_CpG.RDS')
files$V1 <- NULL
colnames(files) <- c('V1', 'V2', 'V3', 'V4', 'V5', 'id')

#reads <- foreach(i=1:nrow(files)) %dopar% {
#for (i in files) {
  read <- files$V1
  read_name <- sapply(strsplit(read, '\\_'), "[", 1)
  id <-  sapply(strsplit(read, '\\_'), "[", 2)
  id_n <- sapply(strsplit(id, ':'), "[", 1)
  files$V1 = paste0(read_name, '_', id_n)
  files <- mutate_at(files, "V2",
                    list(~ ifelse(. == '+', 1, 0))) 
  
  new <- files %>% group_by(V1) %>%
      summarise(Chr=paste0(V3, collapse = ','),
      #Positions = paste0(V2, V4, V3, collapse = ','),
      Range=paste0(V4, collapse = ','),
      Call =paste0(V2, collapse = ','))
  
  new_df <- new %>% separate(col = Chr, into = LETTERS[1:20], sep = ",")
  new_df <- new_df[,c('V1', 'A', 'Range', 'Call')]
  #new_df <- new_df %>% separate(col = Positions, into = LETTERS[2:25], sep = ",")
  #df_final <- new_df[,colSums(is.na(new_df))<nrow(new_df)]
  
  r <- new_df$V1
  X1 <- sapply(strsplit(r, "\\_"), "[", 2)
  V1 <- sapply(strsplit(r, "\\_"), "[", 1)
  new_df$V1 <- V1
  new_df <- add_column(new_df, X1, .after = 1)
  new_df <- as.data.frame(new_df)



readRDS(new_df, '/home/faye/scripts/3DRAM/reads_CpG.RDS')

