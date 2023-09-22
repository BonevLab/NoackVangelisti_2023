#!/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

main_f <- as.character(args[1])
track <- as.character(args[2])
config_f <- as.character(args[3])
descr <- as.character(args[4])
import_adj <- as.logical(args[5])

source(paste0(main_f,'config.R'))
setwd(paste0(map3c_f,'/map3c/'))
Sys.setenv(PIPELINE_HOME = map3c_f)
source("TG3C/imp3c.r")

gtrack.2d.create_from_3Cseq(track, config_f, descr, groot=trackdb, verbose=TRUE, combine_adj=FALSE, skip.if.exist=FALSE, overwrite.if.exist=FALSE, vars_dir="track_vars", interim_step=F,import_adj=import_adj)
