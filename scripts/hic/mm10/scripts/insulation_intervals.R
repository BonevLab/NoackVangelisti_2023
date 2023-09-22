library(misha)
library(ggplot2)
library(cowplot)
library(ggrastr)
theme_set(theme_cowplot())

gsetroot('/home/hpc/bonev/trackdb/mm10/')
options(gmultitasking = F)
options(gmax.data.size = 5e7)

ins_tracks <- c("hic.ES.ins_250","hic.NPC.ins_250")
for (track in c(ins_tracks)){                                ### Create virtual intervals based on average signal so that we can work with variable bin sizes
  gvtrack.create(paste0('v_',track),track,'avg')
}

plot_pair <- function(borders,binSize=1e4,intervals1=NULL,intervals2=NULL,ins_tracks,cols=c('red','blue'),point.size=1){
  iter <- giterator.intervals(intervals=gintervals.all(), iterator=binSize)        #####   Creates intervals from the chosen bin size 
  df <- gextract(c(paste0('v_',ins_tracks,'*(-1)')),intervals = iter ,iterator=iter)        ### Main misha function: extracts values from the given (virtual tracks), given set of intervals and iterators
  df <- df[complete.cases(df),]
  df$labels <- paste0(df$chrom,':',df$start,'-',df$end)
  if(!is.null(intervals1)){
  df$intervals1 <- gintervals.neighbors(df,intervals1)$dist
  } else  {
    df$intervals1 <- 1e8
  }
  if(!is.null(intervals2)){
    df$intervals2 <- gintervals.neighbors(df,intervals2)$dist
  } else  {
    df$intervals_2 <- 1e8
  }
  colnames(df)[4:5] <- c('ES_ins','NPC_ins')
  p <- ggplot(df, aes(x = ES_ins,y = NPC_ins)) + geom_point_rast(size=point.size,alpha=1,col='grey') + geom_point_rast(size=point.size,alpha=1,col=cols[1],data=df[df$intervals1<=1,]) + geom_point_rast(size=point.size,alpha=1,col=cols[2],data=df[df$intervals2<=1,])
  p <- p + geom_abline(slope = 1,intercept = 0,linetype = 2,col='black') + xlab('ES ins') + ylab('NPC ins') + xlim(c(0,5)) + ylim(c(0,5))
  return(p)
}

borders <- gintervals.load('intervals.ES_NPC_unityBorders')
es_borders <- gintervals.load("hic.ES.ins_250_borders")
npc_borders <- gintervals.load("hic.NPC.ins_250_borders")
tss <- gintervals.load("glp.intervals.ucscCanTSS")
ctcf_30k <- gintervals.load('intervals.ES_CTCF_30k')

pdf('example_plot.pdf',height=8,width=8)
p <- plot_pair(borders=borders,binSize=1e4,intervals1=tss,intervals2=NULL,ins_tracks=ins_tracks,cols=c('red','blue'),point.size=1)    ### In this function replace either intervals1 or intervals2 (or both) with any intervals from the selection above. If you want no, or only one set leave the entry as =NULL
print(p)
dev.off()
