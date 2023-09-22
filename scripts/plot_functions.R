
palette.breaks = function(n, colors, breaks){
  colspec = colorRampPalette(c(colors[1],colors[1]))(breaks[1])
  
  for(i in 2:(length(colors)) ){
    colspec = c(colspec, colorRampPalette(c(colors[i-1], colors[i]))(abs(breaks[i]-breaks[i-1])))
  }
  colspec = c( colspec,
               colorRampPalette(c(colors[length(colors)],colors[length(colors)]))(n-breaks[length(colors)])
  )
  colspec
}



plot_paired_smf <- function(res,TF_name,met_lims,cluster_names,cols,mat_name,ggplot_cols,window,binSize,pair_smf_f,meth_colors,window2=10,point.size=4,plot_binSize){
  mat <- as.data.frame(res$res_ratios)
  row.names(mat) <- 1:nrow(mat)
  test <- t(scale(t(mat)))
  test[!complete.cases(test),1] <- mat[!complete.cases(test),1]
  test[!complete.cases(test),2] <- mat[!complete.cases(test),2]
  test[mat[,1]>=met_lims[2]&mat[,2]>=met_lims[2],] <- c(1,1)
  test[mat[,1]<=met_lims[1]&mat[,2]<=met_lims[1],] <- c(-1,-1) 
  mat$cluster = kmeans(mat,centers = 4,iter.max = 10000,nstart = 100)$cluster
  #mat$cluster = cluster::pam(mat,4,cluster.only=T,metric = 'euclidean')
  cluster_means <- ddply(mat,.(cluster),function(x){
    return(c(mean(x[,1],na.rm=T),mean(x[,2],na.rm=T)))
  })
  cluster_means <- cluster_means[order(round(rowMeans(as.matrix(cluster_means[,2:3]),na.rm=T),2),round(cluster_means[,2],1),cluster_means[,3]),]
  cluster_means$names <- cluster_names
  mat$cluster <- factor(cluster_names[match(mat$cluster,cluster_means$cluster)],levels=cluster_names)
  left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = ggplot_cols),
                                                   labels = paste0(round(table(mat$cluster)/length(mat$cluster)*100),' %'),
                                                   labels_gp = gpar(col = "white", fontsize = 10)))
  
  hm <- Heatmap(mat[,1:2]*100, name = mat_name,row_split =mat$cluster,cluster_columns = F,cluster_rows=F,left_annotation = left_annotation,
                cluster_row_slices = FALSE,show_row_names = F,show_column_names = F,show_row_dend = F,col=cols,use_raster=F,row_title = NULL,heatmap_legend_param = list(direction = "horizontal",legend_width = unit(3, "inch")))
  
  test <- data.frame(cluster=as.factor(mat$cluster),TF=TF_name)
  test <- melt(test)
  total <- ddply(test, .(cluster), function(x) nrow(x))[, 2]
  total <- round(total/sum(total)*100,0)
  p <- ggplot(test, aes(x=TF,y=1,fill=cluster)) +  geom_bar(position="fill", stat="identity") + scale_fill_manual(name='cluster',values = as.character(ggplot_cols))
  p <- p + annotate(geom="text", x=1, y=total[4]/200, label=paste0(total[4],'%'),color='black')
  p <- p + annotate(geom="text", x=1, y=total[4]/100+total[3]/200, label=paste0(total[3],'%'),color='black')
  p <- p + annotate(geom="text", x=1, y=total[4]/100+total[3]/100+total[2]/200, label=paste0(total[2],'%'),color='black')
  p <- p + annotate(geom="text", x=1, y=total[4]/100+total[3]/100+total[2]/100+total[1]/200, label=paste0(total[1],'%'),color='black')
  p1 <- p + ylab('Percentage of Total') + ggtitle(paste0(nrow(mat)," total pairs"))+ scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('') + theme(legend.position = 'none',plot.title = element_text(hjust = 0.5))
  
  p2 <- smfCoveragePlot(
    res,
    pair_smf_f,
    binSize,
    mat,
    window,
    meth_colors = meth_colors,
    window2 = window2,
    point.size=point.size,
    plot_binSize=plot_binSize)
  return(list(hm=hm,p1=p1,p2=p2,cluster_id=mat$cluster,mat=mat))
}

plot_smf_shuffle <- function(res, TF_name, met_lims,k=4, cluster_names, cols, mat_name, ggplot_cols){
  mat <- as.data.frame(res$res_ratios)
  row.names(mat) <- 1:nrow(mat)
  test <- t(scale(t(mat)))
  test[!complete.cases(test),1] <- mat[!complete.cases(test),1]
  test[!complete.cases(test),2] <- mat[!complete.cases(test),2]
  test[mat[,1]>=met_lims[2]&mat[,2]>=met_lims[2],] <- c(1,1)
  test[mat[,1]<=met_lims[1]&mat[,2]<=met_lims[1],] <- c(-1,-1) 
  mat$cluster = kmeans(mat,centers = k,iter.max = 1000,nstart = 10)$cluster
  #mat$cluster = cluster::pam(mat,4,cluster.only=T,metric = 'euclidean')
  cluster_means <- ddply(mat,.(cluster),function(x){
    return(c(mean(x[,1],na.rm=T),mean(x[,2],na.rm=T)))
  })
  cluster_means <- cluster_means[order(round(rowMeans(as.matrix(cluster_means[,2:3]),na.rm=T),2),round(cluster_means[,2],1),cluster_means[,3]),]
  cluster_means <- cluster_means[c(1,3,2,4,5),]
  cluster_means$names <- cluster_names
  mat$cluster <- factor(cluster_names[match(mat$cluster,cluster_means$cluster)],levels=cluster_names)
  left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = ggplot_cols),
                                                   labels = paste0(round(table(mat$cluster)/length(mat$cluster)*100),' %'), 
                                                   labels_gp = gpar(col = "white", fontsize = 10)))
  
  hm <- Heatmap(mat[,1:2]*100, name = mat_name,row_split =mat$cluster,cluster_columns = F,cluster_rows=F,left_annotation = left_annotation,
                cluster_row_slices = FALSE,show_row_names = F,show_column_names = F,show_row_dend = F,col=cols,use_raster=F,row_title = NULL,heatmap_legend_param = list(direction = "horizontal",legend_width = unit(3, "inch")))
  
  test <- data.frame(cluster=as.factor(mat$cluster),TF=TF_name)
  test <- melt(test)
  total <- ddply(test, .(cluster), function(x) nrow(x))[, 2]
  total <- round(total/sum(total)*100,0)
  p <- ggplot(test, aes(x=TF,y=1,fill=cluster)) +  geom_bar(position="fill", stat="identity") + scale_fill_manual(name='cluster',values = as.character(ggplot_cols))
  p <- p + annotate(geom="text", x=1, y=total[4]/200, label=paste0(total[4],'%'),color='black')
  p <- p + annotate(geom="text", x=1, y=total[4]/100+total[3]/200, label=paste0(total[3],'%'),color='black')
  p <- p + annotate(geom="text", x=1, y=total[4]/100+total[3]/100+total[2]/200, label=paste0(total[2],'%'),color='black')
  p <- p + annotate(geom="text", x=1, y=total[4]/100+total[3]/100+total[2]/100+total[1]/200, label=paste0(total[1],'%'),color='black')
  p1 <- p + ylab('Percentage of Total') + ggtitle(paste0(nrow(mat)," total pairs"))+ scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab('') + theme(legend.position = 'none',plot.title = element_text(hjust = 0.5))
  return(list(hm=hm, p1=p1, cluster_id=mat$cluster))
  
}

plot_paired_smf_line <- function(pair_smf_f,res,plot_res,lmat,rmat,cluster_id,ylim=NULL,rollmean_k=5){
  df_id <- data.frame(read_id=pair_smf_f$read_name[res$idx],cluster=plot_res$cluster_id)
  lplot_means <- lmat
  lplot_means_clust <- lplot_means[row.names(lplot_means)%in%df_id$read_id[df_id$cluster==cluster_id],]
  lplot <- data.frame(x=as.numeric(names(colMeans(lplot_means_clust,na.rm=T))),y=colMeans(lplot_means_clust,na.rm=T))
  if(rollmean_k!=0){
    lplot$y <- zoo::rollmean(lplot$y,k = rollmean_k, fill = NA)
  }
  lp <- ggplot(lplot,aes(x=x,y=y,group=1)) + geom_line() + xlab('') + ylab('%Meth')
  
  rplot_means <- rmat
  rplot_means_clust <- rplot_means[row.names(rplot_means)%in%df_id$read_id[df_id$cluster==cluster_id],]
  rplot <- data.frame(x=as.numeric(names(colMeans(rplot_means_clust,na.rm=T))),y=colMeans(rplot_means_clust,na.rm=T))
  if(rollmean_k!=0){
    rplot$y <- zoo::rollmean(rplot$y,k = rollmean_k, fill = NA)
  }
  rp <- ggplot(rplot,aes(x=x,y=y,group=1)) + geom_line() + xlab('') + ylab('%Meth')
  if(!is.null(ylim)){
    return(lp+ ylim(ylim)|rp+ ylim(ylim) )
  } else {
    return(lp|rp)
  }
}


plot_hic_scores <- function(scores,cluster_cols,cluster_names,add_stats=F){
  df_m <- melt(scores[,c(grep('score',colnames(scores)),ncol(scores))])
  p1 <- ggplot(df_m,aes(x=cluster,y=value,fill=cluster)) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = F,width=0.8) 
  p1 <- p1 + scale_fill_manual(values=cluster_cols) + xlab('') + ylab('Hi-C Score') + theme(legend.position = "none")
  if(add_stats){
    p1 <- p1 + stat_compare_means(comparisons = list(cluster_names[1:2],cluster_names[2:3],cluster_names[3:4]),label = "p.format",method='wilcox')
  }
  return(p1)
}

plot_linear_scores <- function(scores,cluster_cols,cluster_names,names=c('read1','read2'),chip_name,add_stats=F){
  df_m <- melt(scores[,c(grep(chip_name,colnames(scores)),ncol(scores))])
  df_m$variable <- names[as.numeric(df_m$variable)]
  p1 <- ggplot(df_m,aes(x=cluster,y=value,fill=cluster))+facet_grid(~variable) + geom_boxplot(outlier.size=1,outlier.shape = NA,show.legend = F,width=0.8) 
  p1 <- p1 + scale_fill_manual(values=cluster_cols) + xlab('') + ylab(paste0('Average ',chip_name,' ChIP signal')) + theme(legend.position = "none") 
  if(add_stats){
    p1 <- p1 + stat_compare_means(comparisons = list(cluster_names[1:2],cluster_names[2:3],cluster_names[3:4]),label = "p.signif",method='wilcox')
  }
  return(p1)
}


plotMext <- function(
  INPUTS, xlim=NULL, ylim=NULL, main=NULL, xlab='', ylab='', 
  plotScale='linear', type='full', error.estimates=TRUE, legend=TRUE, 
  legend_ext=FALSE, legend_pos='topright', legend_ext_pos="topleft",
  cex.axis=14, cex.lab=16, cex.main=20, cex.legend=10, ln.v=TRUE, 
  ln.h=NULL, colvec=NULL, pointsize=12, par_mar=c(5.1,4.1,4.1,2.1),invert=FALSE,...) {
  
  cex.axis   <- cex.axis   / pointsize
  cex.lab    <- cex.lab    / pointsize
  cex.main   <- cex.main   / pointsize
  cex.legend <- cex.legend / pointsize
  
  opar <- par()[c('cex', 'mar')]
  
  mm <- length(INPUTS)
  
  anchor <- INPUTS[[1]]$e
  if (ln.v) { ln.v <- c(0, anchor) } else { ln.v <- NULL }
  
  if (plotScale == 'log2') {
    for (i in 1:mm) { 
      co <- diff(range(INPUTS[[i]]$means));
      INPUTS[[i]]$means <- log2(abs(INPUTS[[i]]$means));
      co <- diff(range(INPUTS[[i]]$means))/co 
      INPUTS[[i]]$stderror <- co*INPUTS[[i]]$stderror
      INPUTS[[i]]$conint <- co*INPUTS[[i]]$conint
    }
  } else if  (plotScale == 'zscore') {
    for (i in 1:mm) { 
      co <- diff(range(INPUTS[[i]]$means))
      INPUTS[[i]]$means <- as.numeric(scale(INPUTS[[i]]$means))
      co <- diff(range(INPUTS[[i]]$means))/co
      INPUTS[[i]]$stderror <- co*INPUTS[[i]]$stderror
      INPUTS[[i]]$conint <- co*INPUTS[[i]]$conint
    }
  }
  
  if (is.null(ylim)) {
    ylim[1] <- min( sapply(INPUTS, function(INPUT){ 
      min(INPUT$means - INPUT$conint, na.rm=TRUE) }) )
    ylim[2] <- max( sapply(INPUTS, function(INPUT){ 
      max(INPUT$means + INPUT$conint, na.rm=TRUE) }) )
  }
  
  err_text <- ''
  if( any(is.na(ylim)) | any(is.infinite(ylim)) ) {
    ylim <- c(0,1)
    xlim <- c(0,1)
    err_text <- 'No data to plot'
  }
  
  if ( !is.null(colvec) ) {
    cols <- colvec
  } else {
    cols <- c("darkblue", "darkgreen", "darkred", "darkmagenta", "darkgray",
              "darkorange", "darkcyan", "black", rainbow(mm-8))
  }
  legendText=NULL
  
  if (type == 'legend') {
    legendText <- sapply(INPUTS, '[[', 'desc')
    plot.new()
    if(legend) {
      legend(
        legend_pos, legendText, col=cols, bg=rainbow(1, alpha=0),  
        bty="n", cex=cex.legend, y.intersp=2, inset=0.0, seg.len=2, 
        title.adj=c(2, 2), lwd=15, pch=0, lty=0
      )
    }
    if(error.estimates && legend_ext) { 
      legend(
        legend_ext_pos, 
        c("Mean\n(solid line)","Std error\n(dark area)", 
          "95% CI\n(light area)"), pch=c(-1, 0, 0), lwd=c(3,15,15), 
        lty=c(1,0,0), col=rgb(0,0,0, c(1,0.5, 0.3)), 
        bg=rainbow(1, alpha=0), bty="n", cex=cex.legend, y.intersp=2, 
        inset=0.0, seg.len=2, title.adj=c(2, 2)) 
    }
  } else {
    
    par(mar=par_mar)
    for (i in 1:mm) {
      INPUT <- INPUTS[[i]]
      xinds <- if (is.null(xlim)) range(INPUT$all_ind) else xlim
      if (i==1) { 
        if(!invert){
          plot(
            INPUT$all_ind, INPUT$means, type="n", main=main, xlab=xlab,
            ylab=ylab, xlim=xlim, ylim=ylim, cex.axis=cex.axis,yaxs='i',xaxs='i', 
            cex.lab=cex.lab, cex.main=cex.main, xaxt="n",yaxt='n',...
          ) 
        } else {
          plot(
            INPUT$means,INPUT$all_ind, type="n", main=main, xlab=ylab,
            ylab=xlab, xlim=ylim, ylim=xlim, cex.axis=cex.axis,yaxs='i',xaxs='i', 
            cex.lab=cex.lab, cex.main=cex.main, xaxt="n",yaxt='n',...
          ) 
          
        }
        if(nchar(err_text)) {
          text(0.5, 0.5, err_text, cex = cex.main)
        }
        if( !any(c( 
          length(list(...)[['xaxt']])>0, 
          list(...)[['axes']]==FALSE )
        )) { 
          if(is.null(anchor)) { 
            axis(
              1, at=c(min(xinds), 0,  max(xinds)), 
              labels=c(num2bp(min(xinds)), '0bp', 
                       num2bp(max(xinds))), cex.axis=cex.axis
            ) 
          } 
        }
      }
      if(!invert){
        if(error.estimates) { 
          dispersion(
            INPUT$all_ind, INPUT$means, INPUT$conint, type="l", 
            fill=adjustcolor(cols[i], 0.5), lty=2
          ) 
        }
        if(error.estimates) { 
          dispersion(
            INPUT$all_ind, INPUT$means, INPUT$stderror, type="l", 
            fill=adjustcolor(cols[i], 0.3), lty=2
          ) 
        }
        lines( INPUT$all_ind, INPUT$means, col=cols[i])
      } else {
        if(error.estimates) { 
          dispersion(
            INPUT$means, INPUT$all_ind, INPUT$conint, type="l", 
            fill=adjustcolor(cols[i], 0.5), lty=2
          ) 
        }
        if(error.estimates) { 
          dispersion(
            INPUT$means,INPUT$all_ind, INPUT$stderror, type="l", 
            fill=adjustcolor(cols[i], 0.3), lty=2
          ) 
        }
        lines(INPUT$means, INPUT$all_ind,  col=cols[i])
      }
      legendText <- c(legendText, INPUT$desc)
    }
    if(legend) {
      legend(
        legend_pos, legendText, col=cols, bg=rainbow(1, alpha=0),  
        bty="n", cex=cex.legend, y.intersp=2, inset=0.0, seg.len=2, 
        title.adj=c(2, 2), lwd=15, pch=0, lty=0
      )
    }
    if(error.estimates && legend_ext) { 
      legend(
        legend_ext_pos, 
        c("Mean\n(solid line)", "Std error\n(dark area)", 
          "95% CI\n(light area)"), pch=c(-1, 0, 0), lwd=c(3,15,15), 
        lty=c(1,0,0), col=rgb(0,0,0, c(1,0.5, 0.3)), 
        bg=rainbow(1, alpha=0), bty="n", cex=cex.legend, y.intersp=2, 
        inset=0.0, seg.len=2, title.adj=c(2, 2)
      ) 
    }
    abline( v=ln.v, h=ln.h, col="gray" )
  }
  par(opar)
  
}


plotMisha <- function(main_f,object,targetGene,outDir='figures/',gene_track=readRDS(paste0('/home/hpc/bonev/annotations/',genome,'/',genome,'_gviz.RDS')),out_f=NULL,upstream=1e5,downstream=1e5,plotClusters=levels(object),chipTracksToExtract = gtrack.ls('scATAC.E14'),conditions='',targetEnh=NULL,plotOrder=list(scores=TRUE,anno=FALSE, domains=FALSE, loops=FALSE, VP=TRUE, arcs=FALSE, rna=FALSE, chip=TRUE,axis=FALSE,scATAC=TRUE, genes=TRUE,ideogram=TRUE),cluster_cols,sample_cells=NULL,arcIntervals=NULL,...){
  suppressWarnings(suppressMessages(source(paste0(main_f,'config.R'))))
  targetCell = cells[1]
  targetRegion = targetGene
  targetEnh <- NULL
  scoreTrackToExtract=conditions
  domainsToPlot = paste0("hic.",targetCell,".ins_250_domains_expanded")
  rnaTracksToExtract = ''
  maxDist <- max(upstream,downstream)
  chipRes = 10
  zoomRegion <- c(50-upstream/maxDist*50,50+downstream/maxDist*50)
  pointCEX = 0.4
  plotScale=FALSE
  cex.axis=1.2
  chipYlim <- ''
  rnaYlim=c(0)
  window_scale=1.1
  img_factor=1
  #plotOrder <- list(scores=FALSE,anno=FALSE, domains=FALSE, loops=FALSE, VP=TRUE, arcs=FALSE, rna=FALSE, chip=TRUE,meth=FALSE,axis=FALSE,scATAC=TRUE, genes=TRUE,ideogram=TRUE)
  plotRatios <- list(unitHeight=120, scores=2.5/window_scale, VP=1.5,UMI4C=2.5, loops=2.2, rna=0.6, chip=0.6,MPRA=0.4,meth=0.6, domains=0.15, genes=0.7, arcs=0.7, axis=0.4,ideogram=0.2,scATAC=4,anno=0.15)
  # if(maxDist>=1e5){plotRatios$genes <- plotRatios$genes*1.5} 
  vTypeScore <- 'max'
  vTypeChip <- 'avg'
  zoomInterval <- zoomRegion
  rnaNames <- c('')
  #  plotClusters <- c('NSC','IPC','PN1','PN2')
  chipNames <- ''
  chipColors <- ''
  methNames <- ''
  methColors <- ''
  rnaColors <- NA
  gene_color <- 'brown'
  figure_mode <- TRUE
  geneNames=TRUE
  radius=25e3
  if (!is.null(sample_cells)){
    scCellsToPlot <- sample(colnames(object),sample_cells)
  } else {
    scCellsToPlot <- NULL
  }
  scATAC_window <- 100
  plot.dgram=FALSE
  binSize <- 2e3
  loopThr <- -101
  arcThr <- 59
  widthCorrFactor <- 15
  arcCol <- c(4,3,2)
  arcType <- c(3,4,2)
  getArcScores <- TRUE
  viewpointOnly <- FALSE
  tssOnly <- FALSE
  exactRegion <- TRUE
  ann2D <- gintervals.ls('domains_expanded')[1]
  annIntervals <- gintervals.ls('domains_expanded')[1]
  leftMargin <- 5
  rightMargin <- 3.35
  suppressWarnings(suppressMessages(source(paste0(main_f,'scripts/temp_functions.R'))))
  geneTableDir <- paste0(main_f,'data/extractedGenes/')
  
  if(exists('chipTracksToExtract')){
    chipYlim <- matrix(,length(chipTracksToExtract),2)
  }
  #########################
  if (length(unlist(strsplit(targetRegion,split = ',')))==1){
    targetGene=targetRegion
  } else if (length(unlist(strsplit(targetRegion,split = ',')))==3){
    targetCoordinates=targetRegion
  } else if (length(unlist(strsplit(targetRegion,split = ',')))==6){
    squareCoordinates=targetRegion
  } else {
    stop('Unknown interval format. Exiting ...')
  }
  ######################################
  if(exists('targetCoordinates') & exactRegion == TRUE){maxDist <- (as.numeric(as.character(unlist(strsplit(targetCoordinates,','))[3]))-as.numeric(as.character(unlist(strsplit(targetCoordinates,','))[2])))/2}
  
  ###########################
  ##################################################################################################################################################################
  ##################################################################################################################################################################
  if (exists("scoreTrackToExtract")){
    if (exists('targetGene')){
      outName <- paste0(targetGene,'_',maxDist,'_',paste0(scoreTrackToExtract,collapse='and'),'_',paste0(zoomRegion,collapse='-'))
    } else if (exists('targetCoordinates')) {
      outName <- paste0(paste0(targetCoordinates,collapse='_'),'_',maxDist,'_',scoreTrackToExtract,'_',paste0(zoomRegion,collapse='-'))
    } else { 
      outName <- paste0(paste0(squareCoordinates,collapse='_'),'_',maxDist,'_',scoreTrackToExtract) 
    }
  } else {
    outName <- paste0(ifelse(exists("targetGene"),targetGene,targetCoordinates),'_',maxDist)
  }
  #message(scoreTrack)
  if (!is.null(out_f)){outName <- out_f}
  ###########################
  availableTables <- list.files(geneTableDir, full.names=F)
  plotMar <- c(leftMargin,rightMargin)
  ###########################
  ### LOAD INTERVALS
  tssCoordinates <- gintervals.load(tss_f)
  tssCoordinates <- tssCoordinates[!duplicated(tssCoordinates[,5]),]  #Fix duplicated names
  rownames(tssCoordinates) <- tssCoordinates$geneName
  geneCoordinates <- gintervals.load(genes_f)
  geneCoordinates <- geneCoordinates[!duplicated(geneCoordinates[,5]),]  #Fix duplicated names
  rownames(geneCoordinates) <- geneCoordinates$geneName
  
  annIntervals <- gintervals.load(annIntervals)
  if(!is.null(arcIntervals)){
    if(!is.data.frame(arcIntervals)){   #Assume is granges
      arcIntervals <- data.frame(chrom1=seqnames(arcIntervals),start1=as.numeric(start(arcIntervals)),end1=as.numeric(end(arcIntervals )),chrom2=arcIntervals$gene_chr,start2=as.numeric(arcIntervals$gene_start),end2=arcIntervals$gene_start+1,Correlation=arcIntervals$Correlation,gene_name=as.character(arcIntervals$gene_name))
    }
    arcIntervals <- arcIntervals[arcIntervals$gene_name==targetGene,]
    arcColor_by='Correlation'
    arcColors <- colorRampPalette(c('darkblue','blue','grey80','red','darkred'))
  }
  #### Evaluate additional arguments #####
  list2env(list(...), envir = environment())
  ##################
  ann2D <- gintervals.load(ann2D)
  ann2D <- ann2D[1,]
  ann2D$cluster=1
  if (!is.null(targetEnh)){
    annIntervals <- rbind(tssCoordinates[targetGene,1:3],targetEnh)
    annIntervals$cluster <- targetGene
    annIntervals$cluster[2:nrow(annIntervals)] <- 'Enhancer'
    annIntervals <- intervals.expand(annIntervals,expansion = 1000)
  } else {
    #plotOrder$anno <- FALSE
  }
  arguments <- list(...)
  paste(arguments)
  ##################################################################################################################################################################
  #radius <- 5e3
  if (!is.null(arcIntervals)){
    loopDetails <- arcIntervals
    # radius <- round(abs(distance_f$dist)/20,-3)
    # loopDetails <- gintervals.2d(distance_f[,1],distance_f[,2],distance_f[,3],distance_f[,8],distance_f[,9],distance_f[,10])
    # if(abs(distance_f$dist)>=1e6){
    #   pointCEX = 0.1
    # } else if(abs(distance_f$dist)>=1e5){
    #   pointCEX = 0.25 
    # } else {
    #   pointCEX = 0.5 
    # }
  } else {
    loopDetails <- gintervals.2d(1,1,2,1,1,2)
  }
  
  genePairs <- c()
  ##################################################################################################################################################################
  
  ##################################################################################################################################################################
  if(!exists('binSize')){binSize=5e+3}
  ##################################################################################################################################################################
  # next
  
  plotSecondary <- 0
  if (exists('squareCoordinates'))
  {
    plotSecondary <- 1
    conditions <- ''
    arcsToDraw <- ''
    zoomRegion <- '0|100'
    targetGene <- paste0(squareCoordinates,collapse='|')
    
    chr <- as.character(unlist(strsplit(squareCoordinates, ","))[1])
    start1 <- as.numeric(as.character(unlist(strsplit(squareCoordinates, ","))[2]))
    stop1 <- as.numeric(as.character(unlist(strsplit(squareCoordinates, ","))[3]))
    start2 <- as.numeric(as.character(unlist(strsplit(squareCoordinates, ","))[5]))
    stop2 <- as.numeric(as.character(unlist(strsplit(squareCoordinates, ","))[6]))
    ##########################################
    fullInterval <- gintervals.2d(chr,start1,stop1,chr,start2,stop2)
    currentCoordinates <- gintervals(chr,start1,stop2)
    currentIntrv <- gintervals(chr,start1,stop2)
  }else if (exists('targetCoordinates'))
  {
    posData <- unlist(strsplit(targetCoordinates, ","))
    mid <- floor(as.numeric(as.character(posData[2]))+(as.numeric(as.character(posData[3]))-as.numeric(as.character(posData[2])))/2)
    currentCoordinates <- gintervals(as.character(posData[1]),mid,mid+1)
    #targetGene <- targetCoordinates
  }else{
    currentCoordinates <- tssCoordinates[which(tssCoordinates$geneName == targetGene),]
    centerTarget=TRUE
    if (nrow(currentCoordinates) == 0){stop(paste0('\n-------------------------\n',targetGene,' NOT FOUND IN DATABASE...\nTRY INTERVAL INPUT... (not available yet, but soon...)'),'\n-------------------------\n');}
  }
  
  ###########################
  # next
  
  plotSecondary <- 0
  
  ###########################
  ### COMPUTE 2D INTERVALS TO EXTRACT
  if(plotSecondary == 0)
  {
    iter2d <- makeIter2D_TSS(currentCoordinates,maxDist,binSize)
    currentIntrv <- gintervals(unique(iter2d$chrom1), min(iter2d$start2), max(iter2d$end2))
    fullInterval <- gintervals.2d(currentIntrv$chrom,min(iter2d$start2),max(iter2d$end2),currentIntrv$chrom,min(iter2d$start2),max(iter2d$end2))
  }
  
  ###########################
  plotLength <- currentIntrv$end-currentIntrv$start
  center <- currentIntrv$start+plotLength/2
  
  ###########################
  
  if(exists('zoomRegion'))
  {
    plotRegionStart <- currentIntrv$start+plotLength*(as.numeric(as.character(zoomInterval[1]))/100)
    plotRegionEnd <- currentIntrv$start+plotLength*(as.numeric(as.character(zoomInterval[2]))/100)
    currentIntrv <- gintervals(currentIntrv$chrom, plotRegionStart, plotRegionEnd)
  }
  
  ############ Extract HiC KNN scores for the full region
  if (exists('scoreTrackToExtract'))
  {
    if (length(scoreTrackToExtract)>=1)
    {
      checkFile <- paste(currentCoordinates$chrom,currentCoordinates$start,currentCoordinates$end,maxDist,paste0(targetGene,collapse='_'),paste0(scoreTrackToExtract,collapse='and'),sep='_')
      if (!(checkFile %in% availableTables))
      { 
        print('EXTRACT TABLE...')
        startTime <- proc.time()
        fullInterval2 <- intervals.2d.expand(fullInterval,(fullInterval$end1-fullInterval$start1)*2,(fullInterval$end2-fullInterval$start2)*2) 
        extractedScores <- list()
        for (scoreTrack in scoreTrackToExtract){
          extractedScores[[scoreTrack]] <- gextract(scoreTrack, fullInterval2,band=c(-2e8,-1000))
        }  
        save(extractedScores, file=paste0(geneTableDir,'/',checkFile))
        stopTime <- proc.time() - startTime
        print(paste(checkFile, 'DONE... in ', stopTime['elapsed'], 'seconds'))
      }else{
        print('USING MEMORY DATA...')			
        extractedScores <- get(load(paste0(geneTableDir,'/',checkFile)))
      }
    } else {plotOrder[['scores']]=FALSE}
  }else{scoreTrackToExtract=''; plotOrder[['scores']]=FALSE}
  
  ############ Extract HiC data for the marginals
  if (exists('conditions'))
  {
    if (nchar(conditions[1]) > 1)
    {
      #   scoreTracks <- gtrack.ls(conditions)
      scoreTracks <- conditions
      vTracks <- paste0("v_",scoreTracks)
      for(set in scoreTracks){gvtrack.create(paste0("v_",set),set, vTypeScore)}
      extReads <- gextract(vTracks, iter2d, iterator=iter2d)
    }else{conditions=''; plotOrder[['VP']]=FALSE; plotOrder[['loops']]=FALSE}
  }else{conditions=''; plotOrder[['VP']]=FALSE; plotOrder[['loops']]=FALSE}
  #######################################################

  ############ Extract MPRA data
  if (exists('mpra_tracks'))
  {
    if (nchar(mpra_tracks[1]) > 1)
    {
      #   scoreTracks <- gtrack.ls(conditions)
      if(!exists('mpra_coords')){
        mpra_iter <- mpra_tracks[1]
      } else {
        mpra_iter <- mpra_coords
      }
      mpraData <- gextract(mpra_tracks, currentIntrv, iterator=mpra_coords)
    }else{mpra_tracks=''; plotOrder[['MPRA']]=FALSE}
  }else{mpra_tracks=''; plotOrder[['MPRA']]=FALSE}
  ########################################################
  
  
  
  ############ Extract HiC data for the UMI4C
  if (exists('umi4c_tracks'))
  {
    if (nchar(umi4c_tracks[1]) > 1)
    {
      library("umi4cPackage")
      p4cLoadConfFiles(conf_dir = umi4c_conf_dir)
      umi4c_list <- list()
      for(set in 1:length(umi4c_tracks)){
        bait_coord <- as.numeric(gtrack.attr.get(umi4c_tracks[set][1], "Bait_coord"))
        umi4c_list[[set]] <- p4cNewProfile(umi4c_tracks[set], scope_5=abs(bait_coord-currentIntrv$start), scope_3=abs(bait_coord-currentIntrv$end))
      }
      if(!exists('umi4c_ylim')){umi4c_ylim=NULL}
      if(!exists('umi4c_main')){umi4c_main=NULL}
    }else{umi4c_tracks=''; plotOrder[['UMI4C']]=FALSE; plotOrder[['loops']]=FALSE}
  }else{umi4c_tracks=''; plotOrder[['UMI4C']]=FALSE; plotOrder[['loops']]=FALSE}
  ########################################################
  
  ############ Extract HiC data for the arcs
  
  ########################################################
  
  ############ Extract ChIPseq data for the full region
  if (exists('chipTracksToExtract'))
  {
    if (nchar(chipTracksToExtract[1]) > 1)
    {
      tracks <- as.vector(unlist(sapply(chipTracksToExtract,function(x){gtrack.ls(x)},simplify=T)))
      #tracks <- gtrack.ls(chipTracksToExtract)
      chipTracks <- paste0("v_",tracks)
      for(set in tracks){gvtrack.create(paste0("v_",set),set, vTypeChip)}
      if (is.null(chipRes)){
        chipRes <- 10		
        if(plotLength > 1e+6 & plotLength < 5e+6){chipRes <- 100
        } else if(plotLength > 5e+6){chipRes <- 1000}
      }
      chipData <- gextract(chipTracks, currentIntrv, iterator=chipRes)
      chipTracksToPlot <- unlist(strsplit(chipTracksToExtract, ","))
      if(plotSecondary ==1)
      {
        chipData_UP <- subset(chipData, chipData$start > fullInterval$start1 & chipData$end < fullInterval$end1)
        chipData_DOWN <- subset(chipData, chipData$start > fullInterval$start2 & chipData$end < fullInterval$end2)
      }
      
    }else{chipTracksToExtract=''; plotOrder[['chip']]=FALSE}
  }else{chipTracksToExtract=''; plotOrder[['chip']]=FALSE}
  ########################################################
  ############ Extract MEth data for the full region
  if (exists('methTracksToExtract'))
  {
    #methTracksToExtract <- methTracksToExtract
    if (nchar(methTracksToExtract[1]) > 1)
    {
      #methTracks <- as.vector(unlist(sapply(methTracksToExtract,function(x){gtrack.ls(x)},simplify=T)))
      methTracks <- methTracksToExtract
      methList <- list()
      for (i in 1:length(methTracks)){
        temp_meth <- gextract(methTracks, currentIntrv, iterator=methTracks[i])
        methList[[i]] <- unique(temp_meth) 
      }
      methData <- do.call("rbind", methList)
      methTracksToPlot <- unlist(strsplit(methTracksToExtract, ","))
      if(plotSecondary ==1)
      {
        methData_UP <- subset(methData, methData$start > fullInterval$start1 & methData$end < fullInterval$end1)
        methData_DOWN <- subset(methData, methData$start > fullInterval$start2 & methData$end < fullInterval$end2)
      }
      
    }else{methTracksToExtract=''; plotOrder[['meth']]=FALSE}
  }else{methTracksToExtract=''; plotOrder[['meth']]=FALSE}
  ########################################################
  ############ Extract RNAseq data for the full region
  if (exists('rnaTracksToExtract'))
  {
    if (nchar(rnaTracksToExtract[1]) >1)
    {
      rTracks <- as.vector(unlist(sapply(rnaTracksToExtract,function(x){gtrack.ls(x)},simplify=T)))
      #rTracks <- gtrack.ls(rnaTracksToExtract)
      rnaTracks <- paste0("v_",rTracks)
      for(set in rTracks){gvtrack.create(paste0("v_",set),set, vTypeChip)}
      ####################
      # chipRes <- 10
      if(plotLength > 1e+6 & plotLength < 5e+6){chipRes <- 100
      } else if(plotLength > 5e+6){chipRes <- 1000}
      ####################
      rnaData <- gextract(rnaTracks, currentIntrv, iterator=chipRes)
      rnaData[is.na(rnaData)] <- 0
      
      if(plotSecondary ==1)
      {
        rnaData_UP <- subset(rnaData, rnaData$start > fullInterval$start1 & rnaData$end < fullInterval$end1)
        rnaData_DOWN <- subset(rnaData, rnaData$start > fullInterval$start2 & rnaData$end < fullInterval$end2)
      }
      
    }else{rnaTracksToExtract=''; plotOrder[['rna']]=FALSE}
  }else{rnaTracksToExtract=''; plotOrder[['rna']]=FALSE}
  ########################################################
  # next
  if(!exists('domainsToPlot')){domainsToPlot=''; plotOrder[['domains']]=FALSE}
  
  ##################################################################################################################################################################
  
  ############ PLOT DATA
  # plotGenes <- subset(geneCoordinates, geneCoordinates$chrom == as.character(currentIntrv$chrom) & geneCoordinates$start > currentIntrv$start & geneCoordinates$end < currentIntrv$end)
  plotGenes <- subset(geneCoordinates, geneCoordinates$chrom == as.character(currentIntrv$chrom) & geneCoordinates$end > currentIntrv$start & geneCoordinates$start < currentIntrv$end)
  plotAnn <- subset(annIntervals, annIntervals$chrom == as.character(currentIntrv$chrom) & annIntervals$end > currentIntrv$start & annIntervals$start < currentIntrv$end)
  
  plot2Dann <- subset(ann2D, ann2D$chrom == as.character(currentIntrv$chrom) & ann2D$end > currentIntrv$start & ann2D$start < currentIntrv$end)
  
  hicNames <- paste0('v_',unlist(strsplit(conditions,"\\|")))
  plotLim <- c(currentIntrv$start,currentIntrv$end)
  
  ############ Set output image parameters
  ######################################################
  scoreSets <- length(conditions)
  chipSets <- length(chipTracksToExtract)
  methSets <- length(methTracksToExtract)
  mpraSets <- length(mpra_tracks)
  rnaSets <- length(rnaTracksToExtract)
  domainSets <- length(domainsToPlot)
  umi4cSets <- ifelse(plot.dgram,length(umi4c_tracks)*2,length(umi4c_tracks))
  ######################################################
  if(plotOrder[['scores']] == TRUE){plotOrder[['scores']] <- length(scoreTrackToExtract)}
  if(plotOrder[['UMI4C']] == TRUE){plotOrder[['UMI4C']] <- umi4cSets}
  if(plotOrder[['chip']] == TRUE){
    plotOrder[['chip']] <- chipSets
    atac_ylim <- c(0,ifelse(ncol(chipData[,grep('ATAC',colnames(chipData),ignore.case = F)])==0,2,max(chipData[,grep('ATAC',colnames(chipData),ignore.case = F)],na.rm=T))*1.1)
  }
  if(plotOrder[['MPRA']] == TRUE){
    plotOrder[['MPRA']] <- mpraSets
    if (!exists('mpra_ylim')){
      mpra_ylim <- c(ifelse(ncol(mpraData[,grep('mpra',colnames(mpraData),ignore.case = T)])==0,-1,min(mpraData[,grep('mpra',colnames(mpraData),ignore.case = T)],na.rm=T))*1.1,ifelse(ncol(mpraData[,grep('mpra',colnames(mpraData),ignore.case = T)])==0,5,max(mpraData[,grep('mpra',colnames(mpraData),ignore.case = T)],na.rm=T))*1.1)
    }
  }
  if(plotOrder[['meth']] == TRUE){plotOrder[['meth']] <- methSets}
  if(plotOrder[['rna']] == TRUE){
    plotOrder[['rna']] <- rnaSets
    rna_ylim <- c(-(max(rnaData[,grep('_Rev',colnames(rnaData))],na.rm=T)*1.1),max(rnaData[,grep('_For',colnames(rnaData))],na.rm=T)*1.1)
  }
  if(plotOrder[['domains']] == TRUE){plotOrder[['domains']] <- domainSets}
  
  ######################################################
  plotSets <- sum(unlist(plotOrder))
  imgGrid <- matrix(1:plotSets, nrow=plotSets)
  ######################################################
  plotHeights <- c()
  for (set in names(plotOrder))
  {
    currentSet <- plotOrder[[set]]
    currentRatio <- plotRatios[[set]]
    if (currentSet == TRUE | currentSet == 1){plotHeights <- append(plotHeights, currentRatio)
    }else if(currentSet > 1){plotHeights <- append(plotHeights, rep(currentRatio,currentSet))}
  }
  
  ############ START THE PLOT
  ########################################################
  imgHeight <- sum(plotHeights)*plotRatios[['unitHeight']]
  imgWidth <- img_factor*10*plotRatios[['unitHeight']]
  #pngName <- paste0(outDir,outName,'.png')
  #if(plotOrder[['scores']] == FALSE)
  #{
  if (!plotOrder[['VP']]){
    pdfName <- paste0(outDir,outName,'.pdf')
  } else {
    pdfName <- paste0(outDir,outName,'_',targetCell,'.pdf')
  }
  imgHeight <- imgHeight/96
  imgWidth <- imgWidth/96
  pdf(file=pdfName,width=imgWidth, height=imgHeight)
  #}else{png(file=pngName,width=imgWidth,height=imgHeight)}
  #png(file=pngName,width=imgWidth, height=imgHeight)
  
  #########################################################
  grid_layout <- grid.layout(nrow=nrow(imgGrid), heights=plotHeights, respect=F)
  plot.new()
  pushViewport(viewport(layout=grid_layout))
  
  i=1
  for (set in names(plotOrder))
  {
    
    if (plotOrder[[set]] == FALSE){next
    } else {
      # print(set)
      if(set == 'genes') {
        pushViewport(viewport(layout.pos.row=i, layout.pos.col=1),plotViewport(c(0.1,plotMar[1]-0.3,0.1,plotMar[2])))
      } else if(set == 'ideogram'){
        pushViewport(viewport(layout.pos.row=i, layout.pos.col=1),plotViewport(c(0,plotMar[1]+1,0,plotMar[2]+1)))
      } else if(set == 'scATAC'){
        pushViewport(viewport(layout.pos.row=i, layout.pos.col=1),plotViewport(c(0,plotMar[1]-0.3,0,plotMar[2]-1.3)))
      } else if(set != 'scores'){
        pushViewport(viewport(layout.pos.row=i, layout.pos.col=1))
      }
      par(fig=gridFIG())
      par(new=TRUE)
      # if(set == 'scores'){
      #   for (scoreSet_counter in 1:length(extractedScores)){
      #     if(scoreSet_counter>1){pushViewport(viewport(layout.pos.row=i, layout.pos.col=1));par(fig=gridFIG());par(new=TRUE)}
      #     plotSingleKNN(extractedScores[[scoreSet_counter]], plotAnn = FALSE, plotMar, plotLim, plot2Dann, loopThr, loopDetails, radius,pointCEX=pointCEX,plotScale=plotScale,cex.axis=cex.axis,window_scale=window_scale)
      #     if(length(extractedScores)>1&scoreSet_counter!=length(extractedScores)){par(new=TRUE);upViewport();i=i+1}
      #   }
      if(set == 'scores'){
        for (scoreSet_counter in 1:length(extractedScores)){
          pushViewport(viewport(layout.pos.row=i, layout.pos.col=1),plotViewport(c(0,plotMar[1],0.05,plotMar[2])))
          p <- plotSingleKNN_ggplot(extractedScores[[scoreSet_counter]], plotAnn = FALSE, plotMar, plotLim, plot2Dann,plotThr = loopThr,loopDetails =  loopDetails,radius =  radius,pointCEX=pointCEX,plotScale=plotScale,cex.axis=cex.axis,window_scale=window_scale)
          print(p, vp=viewport(layout.pos.row = i, layout.pos.col = 1))
          if(scoreSet_counter!=length(extractedScores)){
            popViewport();upViewport();i=i+1
          }
        }
      }else if(set == 'UMI4C'){
        if(plot.dgram){
          for (scoreSet_counter in 1:length(umi4c_tracks)){
            if(scoreSet_counter>1){pushViewport(viewport(layout.pos.row=i, layout.pos.col=1));par(fig=gridFIG());par(new=TRUE)}         #FALSE
            plot(umi4c_list[[scoreSet_counter]],min_win_cov=min_umi4c_covs[scoreSet_counter],plot.trend=T,plot.dgram=F,add_xcoord=F,plotMar=c(0,plotMar[1],0.5,plotMar[2]),plotMar2=c(0,plotMar[1],0.5,plotMar[2]),plotOma=c(0,0,0,0),plot.colorbar=plot.colorbar,plot.misha=T,ylim=umi4c_ylim,main=umi4c_main[scoreSet_counter],cex.axis=cex.axis)
            par(new=TRUE);upViewport();i=i+1      #i=2
            pushViewport(viewport(layout.pos.row=i, layout.pos.col=1));par(fig=gridFIG());par(new=TRUE)
            plot(umi4c_list[[scoreSet_counter]],min_win_cov=min_umi4c_covs[scoreSet_counter],plot.trend=F,plot.dgram=T,add_xcoord=F,plotMar=c(0,plotMar[1],0.5,plotMar[2]),plotMar2=c(0.5,plotMar[1],0,plotMar[2]),plotOma=c(0,0,0,0),plot.colorbar=plot.colorbar,plot.misha=T)
            if(length(umi4c_tracks)>1&scoreSet_counter!=length(umi4c_tracks)){par(new=TRUE);upViewport();i=i+1}          
          }
        }else{
          for (scoreSet_counter in 1:length(umi4c_tracks)){
            if(scoreSet_counter>1){pushViewport(viewport(layout.pos.row=i, layout.pos.col=1));par(fig=gridFIG());par(new=TRUE)}
            plot(umi4c_list[[scoreSet_counter]],min_win_cov=min_umi4c_covs[scoreSet_counter],plot.trend=T,plot.dgram=F,add_xcoord=F,plotMar=c(0.5,plotMar[1],0.5,plotMar[2]),plotMar2=c(0.5,plotMar[1],0,plotMar[2]),plotOma=c(0,0,0,0),plot.colorbar=plot.colorbar,plot.misha=T,ylim=umi4c_ylim,main=umi4c_main[scoreSet_counter],cex.axis=cex.axis)
            if(length(umi4c_tracks)>1&scoreSet_counter!=length(umi4c_tracks)){par(new=TRUE);upViewport();i=i+1}
          }
        }
      }else if(set == 'genes'){plot_genes(genome,gene_track=gene_track[[as.character(currentIntrv$chrom)]], currentIntrv, gene_stacking="squish",gene_color=gene_color,gene_size=0.7,fontsize=20,targetGene=plotGenes[plotGenes$geneName%in%targetGene,],collapseTranscripts='meta',geneNames=geneNames)
      }else if(set == 'rna'){
        for (rnaSet in 1:rnaSets){
          if (rnaSet>1){pushViewport(viewport(layout.pos.row=i, layout.pos.col=1));par(fig=gridFIG());par(new=TRUE)}
          combineRNA(rnaData, rnaTracksToExtract[rnaSet], plotMar, yLim=rna_ylim,cex.axis=cex.axis,setNames=rnaNames[rnaSet],rnaColors=rnaColors[rnaSet],figure_mode=figure_mode)
          if(rnaSets>1&rnaSet!=rnaSets){par(new=TRUE);upViewport();i=i+1}
        } 
      }else if(set == 'chip'){
        for (chipSet in 1:chipSets){
          if (chipSet>1){pushViewport(viewport(layout.pos.row=i, layout.pos.col=1));par(fig=gridFIG());par(new=TRUE)}
          if(grepl('ATAC',chipTracksToPlot[chipSet],ignore.case = F)){
            chipYlim[chipSet,] <- t(as.matrix(atac_ylim))
          } 
          plotChIP(chipData, chipTracksToPlot[chipSet], plotMar, yLim=chipYlim,cex.axis=cex.axis,setNames=chipNames[chipSet],chipColors=chipColors[chipSet],figure_mode=figure_mode,chipIndex=chipSet)
          if(chipSets>1&chipSet!=chipSets){par(new=TRUE);upViewport();i=i+1}
        } 
      }else if(set == 'MPRA'){
        message('MPRA y-axis lims:',mpra_ylim[1],' ',mpra_ylim[2])
        for (mpraSet in 1:mpraSets){
          if (mpraSet>1){pushViewport(viewport(layout.pos.row=i, layout.pos.col=1));par(fig=gridFIG());par(new=TRUE)}
          plotBedGraph_hm(mpraData, track=mpra_tracks[mpraSet], plotMar=plotMar,plotLim=plotLim,yLim=mpra_ylim,cex.axis=cex.axis,setName=mpraNames[mpraSet],cols=mpra_colors,idx=mpraSet,figure_mode=figure_mode)
          if(mpraSets>1&mpraSet!=mpraSets){par(new=TRUE);upViewport();i=i+1}
        } 
      } else if(set == 'meth'){
        for (methSet in 1:methSets){
          if (methSet>1){pushViewport(viewport(layout.pos.row=i, layout.pos.col=1));par(fig=gridFIG());par(new=TRUE)}
          plotMeth(methData, methTracksToPlot[methSet], plotMar, yLim=c(0,100),cex.axis=cex.axis,setNames=methNames[methSet],methColors=methColors[methSet],figure_mode=figure_mode,methIndex=methSet,currentIntrv=currentIntrv)
          if(methSets>1&methSet!=methSets){par(new=TRUE);upViewport();i=i+1}
        }  
      }else if(set == 'VP'){plotViewpoint(hicNames, plotLim, targetGene, extReads, plotMar)
      }else if(set == 'loops'){plotLoops(hicNames, plotLim, targetGene, extReads, plotMar)
      }else if(set == 'domains'){loadDomains(domainsToPlot, plotLim, currentIntrv, plotMar)
      }else if(set == 'anno'){plotSimpleBED(plotAnn,plotMar = plotMar,plotLim = plotLim)
      }else if(set == 'arcs'){
        arcIntervals <- subset(arcIntervals, arcIntervals$start1 > plotLim[1] & arcIntervals$end2 < plotLim[2])
        plotArcs(arcIntervals,plotMar=plotMar,color_by=arcColor_by,cols=arcColors,plotLim=plotLim,flip=F,plotScale=F)
      }else if(set == 'axis'){plotAxis(plotLim, plotMar,currentIntrv$chrom,cex.axis=cex.axis,figure_mode=figure_mode)
      }else if(set == 'ideogram'){plotIdeogram(genome,plotLim,currentIntrv$chrom)
      }else if(set == 'scATAC'){
        p <- scCoveragePlot(object=object,region=currentIntrv,window = scATAC_window,annotation=NULL,idents = plotClusters,cells=scCellsToPlot,fragment.path=fragment.path,plotMar=plotMar,cluster_cols = cluster_cols)
        print(p, vp=viewport(layout.pos.row = i, layout.pos.col = 1))
        # 	}else if(set == 'arcs'){plotArcs(arcReads, arcsToDraw, plotLim, widthCorrection=25, loopThr=50)
      }else{plot(0,xlim=plotLim,xlab='', ylab='', yaxt='n', xaxt='n', frame=T, main=set)}
      par(new=TRUE)
      if(set == 'genes'|set=='ideogram'|set=='scATAC'|set=='scores') {popViewport()}
      upViewport()
      i=i+1
    }  
  }
  
  dev.off()
}

plotSingleKNN_ggplot <- function(set1, plotAnn=F, plotMar=c(3.5,3), plotLim, plot2Dann=F, plotThr=-101, loopDetails=FALSE, radius=25e3,pointCEX=0.05,plotScale=TRUE,cex.axis=2,window_scale=1){
  if(is.null(set1)){return()}
  if(window_scale==1){
    set1 <- subset(set1, set1$start1 > plotLim[1] & set1$start2 < plotLim[2])
  } else {
    set1 <- subset(set1, (set1$start1 + set1$start2)/2 > plotLim[1] & (set1$start1 + set1$start2)/2 < plotLim[2])
  }
  set1 <- subset(set1, set1$start1 < set1$start2)
  set1 <- subset(set1, set1[,7] > plotThr)
  ########################################################
  ########################################################
  winSize <- plotLim[2]-plotLim[1]
  yLim <- c(0,winSize)
  # 	yLim <- c(0, max((set1$start2 - set1$start1)))
  ########################################################
  set1 <- set1[order(set1[,7]),]
  plot_label <- colnames(set1)[7]
  plot_label <- gsub('.score_k200','',plot_label)
  plot_label <- gsub('.score_k100','',plot_label)
  plot_label <- gsub('hic.','',plot_label)
  plot_label <- gsub('E14_','',plot_label)
  plot_label <- gsub('ncx_','',plot_label)
  plot_label <- gsub('_IUE24h','',plot_label)
  plot_label <- gsub('NGN2','Neurog2',plot_label)
  plot_label <- ''
  if(grepl('3DRAM',plot_label)|grepl('Methyl',plot_label)|grepl('Bonev',plot_label)){
    plot_label <- ''
  }
  colnames(set1)[7] <- 'score'
  plot_points <- data.frame(x=(set1$start1+set1$start2)/2,y=(set1$start2-set1$start1)/2,score=set1$score)
  p <- ggplot(plot_points,aes(x=x, y=y, color=score)) + theme(legend.position="none") + geom_point_rast(size=pointCEX,raster.dpi = 72) +
    coord_cartesian(xlim=plotLim,ylim=yLim/window_scale) 
  p <- p + scale_color_gradientn(colors=col.scores[1:201],limits=c(-100,100),breaks=c(-100,0,100),name='') +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_cowplot() + theme_void() + theme(legend.position="none",panel.border = element_rect(colour = "black", fill=NA, size=1)) + annotate(geom='text',label=plot_label,x=plotLim[1]+0.95*(plotLim[2]-plotLim[1]),y=0.9*yLim[2]/window_scale,fontface =2,size=10)
  loopDetails <- intervals.2d.centers(loopDetails)
  plot_anno <- data.frame(x=(loopDetails$start1+loopDetails$start2)/2,y=(loopDetails$start2-loopDetails$start1)/2)
  plot_anno$radius <- radius
  p <- p + ggforce::geom_circle(aes(x0=x,y0=y,r=radius),data=plot_anno,color='black',inherit.aes=F,linetype=2) 
  return(p)
}


plot_DMR<- function(con,files,genome,threshold,difference,min_Cov,qval,high_cov_perc,statistical_method,region,output_path, colours,labels) {
  if (!file.exists(paste0(output_path,'DMR_GpC_',qval,'_levels.txt')) && !file.exists(paste0(output_path,'DMR_CpG_',qval,'_levels.txt'))) {
    CpG_meth=methRead(location=list(paste0('/home/hpc/bonev/projects/ram/data/hg38/RAM_covfiles/replicates/Pax6_rep1.NOMe.GpC.cov')
                                    ,paste0('/home/hpc/bonev/projects/ram/data/hg38/RAM_covfiles/replicates/Pax6_rep2.NOMe.GpC.cov')
                                    ,paste0('/home/hpc/bonev/projects/ram/data/hg38/RAM_covfiles/replicates/Tbr2_rep1.NOMe.GpC.cov')
                                    ,paste0('/home/hpc/bonev/projects/ram/data/hg38/RAM_covfiles/replicates/Tbr2_rep2.NOMe.GpC.cov')),
                      pipeline= "bismarkCoverage",
                      sample.id=list('3DRAM_Pax6_rep1','3DRAM_Pax6_rep2','3DRAM_Tbr2_rep1','3DRAM_Tbr2_rep2'),
                      assembly=genome,
                      treatment=c(0,0,1,1),
                      context='GpC')
    filtered_CpG_meth=filterByCoverage(CpG_meth,lo.count=min_Cov, hi.perc=high_cov_perc)
    united_CpG_meth=unite(filtered_CpG_meth, destrand=F)
    
    #Calculate counts at gNOME peaks and call differential regions
    tiled_CpG_meth=regionCounts(united_CpG_meth,region)
    #perform statistical test
    if (statistical_method=='default') {
      myDiff=calculateDiffMeth(tiled_CpG_meth,mc.cores=16)
    } else if (statistical_method=='over_disp_ftest') 
    {
      myDiff=calculateDiffMeth(tiled_CpG_meth,mc.cores=16,overdispersion="MN")
    } else if (statistical_method=='over_disp_chi') 
    {
      myDiff=calculateDiffMeth(tiled_CpG_meth,mc.cores=16,overdispersion="MN",test="Chisq")
    }
    differential_peaks<-getData(myDiff)
    GpC_counts<-getData(tiled_CpG_meth)
    
    GpC_levels<-merge(differential_peaks,GpC_counts, by=c('chr','start','end'))
    GpC_levels$mean_Pax6_meth<-(((GpC_levels$numCs1/(GpC_levels$numTs1+GpC_levels$numCs1)*100))+((GpC_levels$numCs2/(GpC_levels$numTs2+GpC_levels$numCs2)*100)))/2 
    GpC_levels$mean_Tbr2_meth<-(((GpC_levels$numCs3/(GpC_levels$numTs3+GpC_levels$numCs3)*100))+((GpC_levels$numCs4/(GpC_levels$numTs4+GpC_levels$numCs4)*100)))/2
    GpC_levels$DMR<-'Constant'
    GpC_levels$DMR[GpC_levels$meth.diff >0 &GpC_levels$qvalue <= qval] <- "Increase"
    GpC_levels$DMR[GpC_levels$meth.diff <0 &GpC_levels$qvalue <= qval] <- "Decrease"
    write.table(GpC_levels,file=paste0(output_path,'DMR_GpC_',qval,'_levels.txt'),quote=F, sep='\t', row.names = F, col.names = T)
    
    ### Calculate CpG methylation ####
    CpG_meth=methRead(location=list(paste0('/home/hpc/bonev/projects/ram/data/hg38/RAM_covfiles/replicates/Pax6_rep1.NOMe.CpG.cov')
                                    ,paste0('/home/hpc/bonev/projects/ram/data/hg38/RAM_covfiles/replicates/Pax6_rep2.NOMe.CpG.cov')
                                    ,paste0('/home/hpc/bonev/projects/ram/data/hg38/RAM_covfiles/replicates/Tbr2_rep1.NOMe.CpG.cov')
                                    ,paste0('/home/hpc/bonev/projects/ram/data/hg38/RAM_covfiles/replicates/Tbr2_rep2.NOMe.CpG.cov')),
                      pipeline= "bismarkCoverage",
                      sample.id=list('3DRAM_Pax6_rep1','3DRAM_Pax6_rep2','3DRAM_Tbr2_rep1','3DRAM_Tbr2_rep2'),
                      assembly=genome,
                      treatment=c(0,0,1,1),
                      context='CpG')
    filtered_CpG_meth=filterByCoverage(CpG_meth,lo.count=min_Cov, hi.perc=high_cov_perc)
    united_CpG_meth=unite(filtered_CpG_meth, destrand=F)
    tiled_5mC_meth=regionCounts(united_CpG_meth,region)
    CpG_counts<-getData(tiled_5mC_meth)
    differential_peaks<-getData(myDiff)
    CpG_levels<-merge(differential_peaks,CpG_counts, by=c('chr','start','end'))
    CpG_levels$mean_Pax6_meth<-(((CpG_levels$numCs1/(CpG_levels$numTs1+CpG_levels$numCs1)*100))+((CpG_levels$numCs2/(CpG_levels$numTs2+CpG_levels$numCs2)*100)))/2 
    CpG_levels$mean_Tbr2_meth<-(((CpG_levels$numCs3/(CpG_levels$numTs3+CpG_levels$numCs3)*100))+((CpG_levels$numCs4/(CpG_levels$numTs4+CpG_levels$numCs4)*100)))/2
    CpG_levels$DMR<-'Constant'
    CpG_levels$DMR[CpG_levels$meth.diff >0 & CpG_levels$qvalue <= qval] <- "Increase"
    CpG_levels$DMR[CpG_levels$meth.diff <0 &CpG_levels$qvalue <= qval] <- "Decrease"
    write.table(CpG_levels,file=paste0(output_path,'DMR_CpG_',qval,'_levels.txt'),quote=F, sep='\t', row.names = F, col.names = T)
  }
  DMR<-read.table(file=paste0(output_path,'DMR_',con,'_',qval,'_levels.txt'), header=T)

  p<-ggplot(DMR, aes(x=mean_Pax6_meth ,y=mean_Tbr2_meth))+geom_point_rast(size = 0.8, shape = 21,color='grey',fill='grey',data = DMR[DMR$DMR=='Constant',])+geom_point_rast(size = 0.8, shape = 21,aes(color=DMR,fill=DMR),data = DMR[DMR$DMR!='Constant',])
  p <- p+scale_color_manual(values=colours)+scale_fill_manual(values=colours) 
  #p <- p+annotate("text", x =90, y=10, label = paste0("NSC DAR: ",nrow(DMR[grep("Decrease",DMR$DMR), ])), colour='black',size=3)+annotate("text", x =10, y=95, label = paste0("IPC DAR: ",nrow(DMR[grep("Increase",DMR$DMR), ])), colour='black',size=3)+theme(legend.position = "none") + ylab('IPC GpC Accessibility') + xlab('NSC GpC Accessibility')
  return(p)
}

plot_repeat_scatter <- function(res,metric=c('intraHiC','interHiC','CpG','GpC','ins'),p.cutoff=0.05,diff.cutoff=5,min_n=100,xlim,ylim,cols=c('grey','blue','red'),features=NULL,anno.size=4,point.size=2){
  if(metric=='intraHiC'){
    df <- res[res$N_coveredIntraHiC_pairs>=min_n,c(1:3,ncol(res))]
  } else if (metric=='interHiC'){
    df <- res[res$N_coveredInterHiC_pairs>=min_n,c(5:7,ncol(res))]
  } else if(metric=='CpG'){
    df <- res[res$N_coveredCpG>=min_n,c(9:11,ncol(res))]
  } else if(metric=='GpC'){
    df <- res[res$N_coveredGpC>=min_n,c(13:15,ncol(res))]
  } else if(metric=='ins'){
    df <- res[res$N_coveredINS>=min_n,c(17:19,ncol(res))]
  }
  colnames(df) <- c('NSC','IPC','sign','TE_family')
  df$sign <- p.adjust(df$sign,method = 'BH')
  df$TE_subfamily <- row.names(df)
  p <- ggplot(df,aes(x=NSC,y=IPC,fill=sign)) + geom_point_rast(data = df[df$sign>p.cutoff|abs(df$NSC-df$IPC)<=diff.cutoff,],fill=cols[1],pch = I(21),size = point.size) + geom_point(data = df[df$sign<=p.cutoff&df$NSC-df$IPC>diff.cutoff,],fill=cols[2],pch = I(21),size = point.size)+ geom_point(data = df[df$sign<=p.cutoff&df$IPC-df$NSC>diff.cutoff,],fill=cols[3],pch = I(21),size = point.size)
  if (is.null(features)){
    df1 <- df[df$sign<=p.cutoff&df$NSC-df$IPC>diff.cutoff,]
    df2 <- df[df$sign<=p.cutoff&df$IPC-df$NSC>diff.cutoff,]
    features <- c(head(row.names(df1[order(df1$NSC-df1$IPC,decreasing=T),]),10),head(row.names(df2[order(df2$IPC-df2$NSC,decreasing=T),]),10))
  }
  p <- p + ggrepel::geom_text_repel(
    data = df, size = anno.size,seed = 42,
    box.padding =0.8, min.segment.length = 0,max.iter = 10000,
    aes(x=NSC,y=IPC,color=NULL,label=ifelse(TE_subfamily%in%features, as.character(TE_subfamily), "")),force=10)
  p <- p + theme(legend.position = "none") 
  return(list(p=p,df=df))
}

plotBedGraph_hm <- function(df,track, plotMar=c(3.5,3), plotLim,yLim,cex.axis,setName,cols=c('purple','#D4AF37'),idx,figure_mode=FALSE){
  set <- df[,c('chrom','start','end',track)]
  gen_color <- scales::col_numeric(palette = cols,domain = yLim,na.color = 'white')
  set$col <- gen_color(set[,4])
  
  par(mar = c(0.1,plotMar[1],0.1,plotMar[2]))
  
  #par(mar = c(0, plotMar[1], 0, plotMar[2]))
  
  plot(0 ,ylim=c(0,50), xlim=plotLim, xlab='', ylab='', yaxt='n', xaxt='n',frame=F,yaxs='i',xaxs='i')
  
  for (i in 1:nrow(set))
  {
    rect(set[i,2], 0, set[i,3], 50, col=set[i,'col'], border = NA)
    #text(set[i,4], x=set[i,2]+(set[i,3]-set[i,2]),y=20, pos=4, cex=1.3, las=2)
  }
  mtext(setName, side = 4, outer = F,line=0,at=25,adj=1,las=2,cex=cex.axis)
}


