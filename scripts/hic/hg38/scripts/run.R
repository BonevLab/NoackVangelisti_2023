#!/bin/env Rscript
source('/home/hpc/bonev/projects/hic/sc/config.R')

source(paste0(main_f,'scripts/main_functions.R'))
source(paste0(main_f,'scripts/aux_functions.R'))
source(paste0(main_f,'scripts/plot_functions.R'))
source(paste0(main_f,'scripts/temp_functions.R'))

library(doParallel)
registerDoParallel(8)

### Inititalize genome parameters ### 
chrom_sizes <- gintervals.all()
#### 


########################
##### Prepare Tracks ###
########################

submit_prepareTracks(ins_window=250,k=100,cell=cells[2],eigenBin=250,mem=40,near_cis=5e6,expand_cis=5e5)

########################



### Calculate matrix bins and output square matrix

path <- paste0(main_f,'analysis/compartments/')
dir.create(paste0(main_f,'analysis/compartments/'), showWarnings = FALSE,recursive=T)
binSize=2e5
chrs=c('chr3')
file_f <- paste0(path,'chr3_',binSize/1000,'kb')

test <- extractBinned(tracks=all_tracks,cells=cells,chrs=chrs,path=path,binSize=binSize,file_f=file_f)
plotBinned(extra_tracks=meth_tracks,cells=cells,init_cell="E14_NSC",chrs='chr3',binSize=binSize,path=path,file_f=file_f,plot_what='obs',balance=T)

#to get the compartment strenth and on top the chips
rank_matrix(eigen_tracks=gtrack.ls('eigen','250'),cells=cells,chrs=c(1:19,'X'),path=path,binSize=2.5e5,min_dist=1e7,ranks=100)

file_f <- paste0(path,'rankMatrix_',binSize/1000,'kb_ranks',100)
plot_rankMatrix(file_f=file_f,zlim=c(-1.2,1.2),cells=cells,col=wide_red_blue_pal(1000),plot_chip=FALSE)

########################
##### Compartments #####
########################

res <- analyzeCompartments(cells=cells,hic_tracks = all_tracks,insTracks = ins_tracks,ins_dom_thresh=0.1,insulation=250,min_domainSize = 5e4)

########################

res_comp <- comp_boundaries(type='eigen',cells=cells, write_borders=T)
res_comp <- comp_boundaries(type='amos',cells=cells, write_borders=T)
# run in R and print res_comp to get the transition in terms of % for amos and in terms of number for eigen
# the fines generated van be used withintersect to get the CTCF sites overthe transition an then use them in contact enrichment figure


res_cor <- comp_correlation(cells=cells,e_tracks=eigen_tracks,tracks=atac_tracks[-c(4)],binSize=1e5,type='chip',cutoff=NULL,method='spearman')
res_cor <- comp_correlation(cells=cells,e_tracks=eigen_tracks,binSize=1e5,marks=c('OIS_0hrs_Replicate1'),type='repli',cutoff=0.75)

test <- corrSamples(tracks=all_tracks,binSizes=1e5,minDist=1e4,maxDist=max(gintervals.all()$end),cutoff=5,cells_vector=rep(cells,each=3),path=paste0(main_f,'analysis/compartments/'))

comp_expression(cells=cells,K=6,beanplot_f=T,cutoff=NULL,cell_colors=rainbow(length(cells)),path_f=paste0(main_f,'analysis/compartments/'))

########################
##### Comp cisDecay ####
########################

path <- paste0(main_f,'analysis/cis_decay/')
dir.create(path, showWarnings = FALSE,recursive=T)
file_name <- 'IUE_cisDecay'
run_cisDecay(tracks=all_tracks,cells=cells[4:5],path=path,log_base=10,file_f=file_name,maxDist=1e8,colors=c('darkgreen','darkred'),alpha=0.3,disp_colors=c('lightgreen','pink'))


path <- paste0(main_f,'analysis/cis_decay/')
dir.create(path, showWarnings = FALSE,recursive=T)
file_name <- 'comp_domains'
cells <- c('D2','D4')
comp_cisDecay(comp_tracks=all_tracks,cells=cells,path=path,ins=250,file_f=file_name)
plot_comp_cisDecay(plot_data=paste0(path,file_name),plot_cell=c('D2','D4'),plot_what='BvsB',min_dist=10e6,max_dist=1.5e8,fig.height=6,fig.width=6,smooth=T,legend_pos="topright")    #ylim=c(value1,value2)
plot_comp_cisDecay(plot_data=paste0(path,file_name),plot_cell=c('D2','D4'),plot_what='total_AvsB',min_dist=10e6,max_dist=1.5e8,fig.height=6,fig.width=6,smooth=T,legend_pos="topleft",returnData=T)    #ylim=c(value1,value2)

#####################
##### Insulation ####
#####################

path <- paste0(main_f,'analysis/insulation/')
dir.create(path, showWarnings = FALSE,recursive=T)

res <- analyzeInsulation(tracks=ins_tracks,min_nreps = 2,rep_tracks=ins_rep_tracks,path=paste0(main_f,'analysis/insulation/'),cells=cells,ins_dom_thresh=0.1,border_max_distance=1e4,sig_thresh=0.8,no_sig_thresh=0.2,use.normFACS=TRUE,use.sigThresh=TRUE,write.insTracks=T,write.borders=T,anyDiff=T)
analyzeDiffBorders(input=res$differential,expression_f='',path=paste0(main_f,'analysis/insulation/'),cells=cells,scale=F,clustering='hierarchical',method='ward.D2',K=4,write.borders=T)
write.table(res$constant[,1:3],paste0(path,'conserved_borders.bed'),quote=F,row.names=F,col.names=F,sep='\t')


#####################
##### AverageTAD ####
#####################

path <- paste0(main_f,'analysis/averageTAD/')
dir.create(path, showWarnings = FALSE,recursive=T)

averageTAD(tracks=all_tracks,cells=cells,score_tracks = score_tracks,domain_size=c(5e4,4e6),file_f='averageTAD',path=paste0(main_f,'analysis/averageTAD/'),bins=100,stats_f=F)

averageTAD(tracks=all_tracks,cells=cells,score_tracks = score_tracks,domain_size=c(5e4,4e6),file_f='averageTAD',path=paste0(main_f,'analysis/averageTAD/'),bins=100,stats_f=T)

zlim=c(-1.5,1)

blue_white_pal = colorRampPalette(c("purple", "navy", "blue", "#87FFFF", "white"))
white_red_pal = colorRampPalette(c("white","#FF413D", "black", "orange", "yellow"))
averageTAD_colors <- c(blue_white_pal(1000*length(seq(zlim[1],0,by=0.1))/length(seq(zlim[1],zlim[2],by=0.1))),white_red_pal(1000*length(seq(0,zlim[2],by=0.1))/length(seq(zlim[1],zlim[2],by=0.1))))

plot_averageTAD(tracks=all_tracks,add_matrix=res1,mat_lim=c(60,90),cells=cells[1],path=path,file_f='averageTAD',stats_f=F,plot_what='',zlim=zlim,flip=T,z_colors=averageTAD_colors,height=500,width=800)


res1 <- averageTrack(tracks="methylation.E14_NSC_10x",regions="hic.E14_NSC.ins_250_domains_expanded",bins=100,anchor_middle=F)
res2 <- averageTrack(tracks="methylation.E14_IPC_10x",regions="hic.E14_IPC.ins_250_domains_expanded",bins=100,anchor_middle=F)
res3 <- averageTrack(tracks="methylation.E14_PN_10x",regions="hic.E14_PN.ins_250_domains_expanded",bins=100,anchor_middle=F)


#######################################
##### aggregateHiC_quantification #####
#######################################
intervals_f <- list.files('/home/hpc/bonev/projects/SC/E14_analysis/results/HiC/motif_beds_top5000/',pattern = 'bed',full.names = T)
shuffled_int <- intervals_f[grep('shuffle',intervals_f)]
intervals_f <- intervals_f[grep('shuffle',intervals_f,invert=T)]
for (i in 951:1000){
submit_cisDecayIntervals(cells=cells,window_i=10000,
                           intervals1='/home/hpc/bonev/data/hic/data/peaks/final/NPC_Pax6',
                           intervals2='/home/hpc/bonev/data/hic/data/peaks/final/NPC_Pax6',
                           grid_mode='1D',domains = "hic.ncx_Hes5.ins_250_domains_expanded")
}


path <- paste0(main_f,'analysis/cis_decay/')
files <- list.files(paste0(main_f,'data/cis_decay/'),full.names=T)
files <- files[grep('score',files,invert=T)]


### intra composite
dist = 10^seq(4,8,by=0.1)
min_dist=1e4
max_dist=1e7
ylim=c(0,2)
name <- 'Ngn2'	
min_ind <- findInterval(min_dist,dist)
max_ind <- findInterval(max_dist,dist)
dist <- dist[min_ind:max_ind]
pdf(paste0(path,name,'_intra.pdf'),width=8,height=6)
plot(log10(dist),y=rep(NA,length(dist)),ylim=ylim,type='n',main=name,ylab="log2(obs/exp)",xaxt='n',xlab='')
magaxis(1,majorn=5,unlog=TRUE,grid=T)
plot_cis_decay(plot_data=files[grep('NPC_Neurog2.bed_NPC_Neurog2.bed.1D.10000',files)],plot_cell=cells,plot_what='intra',min_dist=min_dist,max_dist=max_dist,composite=T,smooth=T,returnData=T)
plot_cis_decay(plot_data=files[grep('D2_CTCF',files)],plot_cell=c('D2'),plot_what='intra',min_dist=min_dist,max_dist=max_dist,composite=T,colors=c('green'),disp_colors=c("lightgreen"),smooth=T,returnData=F)
plot_cis_decay(plot_data=files[grep('D10_CTCF',files)],plot_cell=c('D10'),plot_what='intra',min_dist=min_dist,max_dist=max_dist,composite=T,colors=c('purple'),disp_colors=c("darkred"),smooth=T,returnData=F)
legend(x="topleft",legend=c('D0','D2','D10'),cex=1.5,lty=rep(1,2),col=c('blue','green','purple'))
grid(nx=0,ny=NULL)
dev.off()


#to plot the cis decay on the chip marks###
name <- 'D0_CTCF'	
plot_cis_decay(plot_data=files[grep('D0_CTCF_For',files)],plot_cell=c('D0'),plot_what='intra',composite=F,smooth=T,returnData=F,min_dist=1e4,max_dist=1e7,fig.height=6,fig.width=10,legend_pos='topleft')


dist = 10^seq(4,8,by=0.1)
min_dist=1e7
max_dist=1e8
ylim=c(-0.5,0.5)
name <- 'CTCF'	
min_ind <- findInterval(min_dist,dist)
max_ind <- findInterval(max_dist,dist)
dist <- dist[min_ind:max_ind]
pdf(paste0(path,name,'inter.pdf'),width=8,height=6)
plot(log10(dist),y=rep(NA,length(dist)),ylim=ylim,type='n',main=name,ylab="log2(obs/exp)",xaxt='n',xlab='')
magaxis(1,majorn=5,unlog=TRUE,grid=T)
plot_cis_decay(plot_data=files[grep('D0_CTCF',files)],plot_cell=c('D0'),plot_what='inter',min_dist=min_dist,max_dist=max_dist,composite=T,colors=c('blue'),disp_colors=c("lightblue"),smooth=T,returnData=F)
plot_cis_decay(plot_data=files[grep('D2_CTCF',files)],plot_cell=c('D2'),plot_what='inter',min_dist=min_dist,max_dist=max_dist,composite=T,colors=c('red'),disp_colors=c("pink"),smooth=T,returnData=F)
legend(x="topleft",legend=c('D0','D2'),cex=1.5,lty=rep(1,2),col=c('blue','red'))
grid(nx=0,ny=NULL)
grid(nx=0,ny=NULL)
dev.off()

plot_cis_decay(plot_data=files[grep('D0_CTCF',files)],plot_cell=c('D0'),plot_what='inter',composite=F,smooth=T,returnData=F,min_dist=1e6,max_dist=1e8,fig.height=6,fig.width=10,legend_pos='topleft')

name <- 'all'
plot_cis_decay(plot_data=files[grep('D0_H3K27me3',files)],plot_cell=c('D0'),plot_what='inter',composite=T,smooth=T,returnData=T,min_dist=1e6,max_dist=1e8,fig.height=6,fig.width=10,legend_pos='topleft')
plot_cis_decay(plot_data=files[grep('D2_H3K27me3',files)],plot_cell=c('D2'),plot_what='inter',composite=T,smooth=T,returnData=T,min_dist=1e6,max_dist=1e8,fig.height=6,fig.width=10,legend_pos='topleft')
plot_cis_decay(plot_data=files[grep('D6_H3K27me3',files)],plot_cell=c('D6'),plot_what='inter',composite=T,smooth=T,returnData=F,min_dist=1e6,max_dist=1e8,fig.height=6,fig.width=10,legend_pos='topleft')
plot_cis_decay(plot_data=files[grep('D10_H3K27me3',files)],plot_cell=c('D10'),plot_what='inter',composite=T,smooth=T,returnData=F,min_dist=1e6,max_dist=1e8,fig.height=6,fig.width=10,legend_pos='topleft')

grid(nx=0,ny=NULL)
dev.off()






########################
##### aggregateHiC #####
########################

for (cell in cells){
  submit_aggregateHiC(cells=cell,tracks=all_tracks,range_f=40000,res_f=1000,filter_f=0,
                      intervals1='/home/hpc/bonev/projects/ram/data/hg38/beds/IPC_ATAC_NGN2_5000.bed',
                      intervals2='/home/hpc/bonev/projects/ram/data/hg38/beds/IPC_ATAC_NGN2_5000.bed',
                      domains=paste0("hic.HMGU1_D45_Tbr2.ins_250_domains_expanded"),
                      grid_mode='1D',mem=80)
  }

for (cell in cells){
  for (cluster in c('NSC','IPC','PN')){
    for (corre in c('neg','no')){
      submit_aggregateHiC(cells=cell,tracks=all_tracks,range_f=60000,res_f=1000,filter_f=0,
                          intervals1=paste0("/home/hpc/bonev/projects/SC/E14_analysis/results/beds/",cluster,"_",corre,"Cor_GA"),
                          intervals2=paste0("/home/hpc/bonev/projects/SC/E14_analysis/results/beds/",cluster,"_",corre,"Cor_DA"),
                          domains=paste0("hic.E14_",cluster,".ins_250_domains_expanded"),
                          grid_mode='merged',mem=80)
    #  submit_aggregateHiC(cells=cell,tracks=score_tracks,range_f=100000,res_f=10000,filter_f=0,
    #                      intervals1=paste0("/home/hpc/bonev/projects/SC/E14_analysis/results/beds/",cluster,"_",corre,"Cor_GA"),
    #                     intervals2=paste0("/home/hpc/bonev/projects/SC/E14_analysis/results/beds/",cluster,"_",corre,"Cor_DA"),
    #                      domains=paste0("hic.E14_",cluster,".ins_250_domains_expanded"),use_score=TRUE,
    #                      grid_mode='merged',mem=80)
    }
  }
}


for (cell in cells){
  for (cluster in c('NSC','IPC','PN')){
    for (corre in c('pos')){
  plot_aggregateHiC(cells=cell,pool=T,intervals1=paste0(cluster,"_",corre,"Cor_GA"),intervals2=paste0(cluster,"_",corre,"Cor_DA"),range_f=80000,filter_f=0,res_f=1000,plot_res=6000,grid_mode='merged',zlim=c(-0.75,0.75),which_plot=c(1),plot_scale = F,interval1_name = 'TSS',interval2_name = 'Distal')
    }
  }
}



for (cell in cells){
  plot_aggregateHiC(cells=cell,pool=T,intervals1='Pou3f2.bed',intervals2='Pou3f2.bed',range_f=40000,filter_f=0,res_f=1000,plot_res=2000,grid_mode='1D',zlim=c(-0.5,0.5),which_plot=c(2,4),plot_scale = F,interval1_name = '',interval2_name = '',add_plot = F)
}



cells <- c('D0','D2')
for (cell in cells){
  submit_aggregateDiagonal(cells=cell,tracks=all_tracks,intervals_f=paste0('/work/project/Cavalli-mammals/satish/HiC/data/peaks/processed/',cell,'_CTCF.bed'),consider_strand=F,mem=50)
}


cells <- c('D0','D2')
for (cell in cells){
  plot_aggregateDiagonal(cells=cell,pool=T,intervals1=paste0(cell,'_CTCF.bed'),res=500,consider_strand='FALSE',zlim=c(-1,1))
}



#################################
##### process_peaks and rna #####
#################################

cells=c('D2')
for (cell in cells){
  process_peaks(peaks=paste0(cell,'_CTCF.bed'),chip_max_peaks=NULL,peaks_mode='sharp',binSize=1000,cut_off=0.98)
  # 	process_peaks(peaks=paste0(cell,'_H3K4me3.bed'),chip_max_peaks=NULL,peaks_mode='sharp',binSize=1000,cut_off=0.98)
  # 	process_peaks(peaks=paste0(cell,'_H3K27Ac.bed'),chip_max_peaks=20000,peaks_mode='sharp',binSize=1000,cut_off=0.98)
  
}	
cells <- c('D6','D10','prolif','resen')
rank_Tss(cells=cells,tss=gintervals.load(tss_f),fpkm_thresh=1,ranks=4)
rank_Exons(cells=cells,fpkm_thresh=1,fpkm_ranks=4,exon_ranks=6)

## Match expression 

path <- '/work/project/Cavalli-mammals/satish/HiC/data/rna/'
file_f <- '/work/project/Cavalli-mammals/satish/HiC/data/peaks/senescence_genes.txt'
which_column <- 1
cellToMatch <- 'D0'
int_name <- paste0('SASP_matched',cellToMatch)

tss <- gintervals.load('glp.intervals.ucscCanTSS')
allGenes <- gintervals.load('glp.intervals.ucscCanTSS')
if (grepl('txt',file_f)|grepl('bed',file_f)) {
  peaks <- read.table(file_f,header=F)
  genes <- tss[tss$geneName %in% peaks[,which_column],]
} else {
  peaks <- get(load(file_f))
  genes <- tss[tss$geneName %in% peaks[,which_column],]
}
expression <- read.table("/work/project/Cavalli-mammals/satish/HiC/data/fpkm_genes.txt",header=T)
expression$FPKM <- rowMeans(expression[,grep(cellToMatch,colnames(expression))],na.rm=T)
df <- sampleMatchedExpression(genes,allGenes,expression,0.2)
temp <- df[,1:5]
save(temp, file=paste0(path,int_name))

##### 

###################
##### diffHiC #####
###################


for (cell in cells[3]){
  for (chr in gintervals.all()$chrom){
    generate_diffHiC(cell=cell,chr=chr,mode='exp',mem=30)
  }
}

binSizes=c(1e4,2.5e4)
for (chr in chrom_sizes[,1]){
  for (binSize in binSizes){
    analyze_diffHiC(cells=cells,binSize=binSize,chr=chr,name='Test',mem=40,read_thresh=1)
  }
}


#######################
##### Methylation #####
#######################

#for mC to get quantiles ##
path <- paste0(main_f,'analysis/methylation/100kb/')
dir.create(path, showWarnings = FALSE,recursive=T)
# methylation_quantiles(tracks=meth_tracks,cells=c('D0','D2','D6','prolif','resen'),gpath=path,bin=1e4,inters=gintervals.all(), quantile_cutoff=c(60,80),file_f='all_60_40.pdf',fig.height=5,fig.width=6)
methylation_quantiles(tracks=meth_tracks,cells=c('D0','D2','D6','prolif','resen'),gpath=path,bin=1e4,inters=gintervals.all(), quantile_cutoff=c(80,100),file_f='all_80_100.pdf',fig.height=7,fig.width=6)
# methylation_quantiles(tracks=meth_tracks,cells=c('D0','D2','D6','prolif','resen'),gpath=path,bin=1e4,inters=gintervals.all(), quantile_cutoff=c(95,100),file_f='test.pdf',fig.height=5,fig.width=6)


#for repliseq to get quantiles ##
path <- paste0(main_f,'analysis/repliseq/100kb/')
dir.create(path, showWarnings = FALSE,recursive=T)
# methylation_quantiles(tracks=repli_tracks,cells=c('0hrs','12hrs','24hrs','36hrs'),gpath=path,bin=1e5,inters=gintervals.all(), quantile_cutoff=c(0,20),file_f='new.pdf',fig.height=5,fig.width=6,cex.legend=0.5)
# methylation_quantiles(tracks=repli_tracks,cells=c('0hrs','12hrs','24hrs','36hrs'),gpath=path,bin=1e5,inters=gintervals.all(), quantile_cutoff=c(20,40),file_f='new1.pdf',fig.height=8,fig.width=7,cex.legend=0.5)
# methylation_quantiles(tracks=repli_tracks,cells=c('0hrs','12hrs','24hrs','36hrs'),gpath=path,bin=1e5,inters=gintervals.all(), quantile_cutoff=c(40,60),file_f='test1.pdf',fig.height=6,fig.width=7,cex.legend=0.5)
# methylation_quantiles(tracks=repli_tracks,cells=c('0hrs','12hrs','24hrs','36hrs'),gpath=path,bin=1e5,inters=gintervals.all(), quantile_cutoff=c(60,80),file_f='test.pdf',fig.height=5,fig.width=6,cex.legend=0.5)
methylation_quantiles(tracks=repli_tracks,cells=c('0hrs','12hrs','24hrs','36hrs'),gpath=path,bin=1e5,inters=gintervals.all(), quantile_cutoff=c(80,100),file_f='test.pdf',fig.height=5,fig.width=6,cex.legend=0.5)



######## Testing ###############

submit_aggregateHiC2(cells=cells[4],tracks=all_tracks,domains="hic.D6.A_domains",range_f=40000,res_f=1000,filter_f=0,intervals1='/work/project/Cavalli-mammals/satish/HiC/data/rna/D6_activeTss',intervals2='/work/project/Cavalli-mammals/satish/HiC/data/rna/D6_activeTss',grid_mode='1D',mem=40,long_range=T)
submit_aggregateHiC2(cells=cells[4],tracks=all_tracks,domains="hic.D6.B_domains",range_f=40000,res_f=1000,filter_f=0,intervals1='/work/project/Cavalli-mammals/satish/HiC/data/rna/D6_activeTss',intervals2='/work/project/Cavalli-mammals/satish/HiC/data/rna/D6_activeTss',grid_mode='1D',mem=40,long_range=T)





submit_APA(hic_file='/work/project/Cavalli-mammals/boyan/HiC/juicer/projects/CN_mapq30.hic',domains='hic.CN.ins_250_domains_expanded',window_i=55000,intervals1=paste0('/work/project/Cavalli-mammals/boyan/HiC/data/geneBodies/CN_activeTss_exonQ2_fpkmQ4'),intervals2=paste0('/work/project/Cavalli-mammals/boyan/HiC/data/geneBodies/CN_activeTss_exonQ2_fpkmQ4'),grid_mode='trans',domain_mode='inter',min_dist=5e5,max_dist=2e6,window_f=5000000,res=50000,save_folder='CN_activeTSS_exonQ2_fpkmQ4_trans_50kb',mem=20)





### Aggregate diffHiC Hi-C results ####

contrast <- 'E14_IPC-E14_NSC'
path <- paste0("/home/hpc/bonev/projects/hic/sc/analysis/diffHiC/Test/25kb//loess/",contrast,"/direct/")
files <- list.files(path,recursive=T,full.names=T)
df <- read.table(files[1],header=T)
for (i in 2:length(files)){
  df <- rbind(df,read.table(files[i],header=T))
}
write.table(df,'/home/hpc/bonev/projects/hic/sc/analysis/diffHiC/Test/25kb/aggregate_results.tsv',quote=F,sep='\t',row.names=F,col.names=T)






