############################################################################################################################################
##
## plotRegiont.R -- plot a snapshot of a genomic region
## 30 november 2018
## Thomas Faux
## Medical Bioinformatics Centre
############################################################################################################################################

## makeCoverages -- produce the coverage for each replicate

## makeBAMLayout -- create the layout from the BAM info

## getCoverage -- compute the coverage with Rsamtools

## plotAverageBAM -- plot the average of each group

## plotBAM -- plot the bam regions

# produce the coverage for each replicate
#
# @param BAMlayout layout created by the function makeBAMLayout
# @param region GRanges object containing the region to plot
# @return returns a coverages object
############################################################################################################################################
##
## makeCoverages -- produce the coverage for each replicate
##
## BAMlayout   -- layout : layout created by the function makeBAMLayout
##
## region      -- vector : coordiantes of the region, cHromosom, start and end
##
## returns a coverages object
############################################################################################################################################
makeCoverages <- function(region,BAMlayout){
  coverages <- list()
  for(i in 1:length(BAMlayout)){
    coverages[[i]] <- list()
    for(j in 1:length(BAMlayout[[i]])){
      coverages[[i]][[j]] <- getCoverage(region,BAMlayout[[i]][[j]])
    }
  }
  return(coverages)
}

# compute the coverage with Rsamtools
#
# @param bamfile character path to the BAM file
# @param region GRanges object containing the region to plot
# @return returns the coverage of one file
############################################################################################################################################
##
## getCoverage -- compute the coverage with Rsamtools
##
## bamfile  -- character : path to the BAM file
##
## region   -- vector    : coordiantes of the region, chromosom, start and end
##
## returns the coverage of one file
############################################################################################################################################
getCoverage <- function(region, bamfile){

  params <- ScanBamParam(which = region, what = scanBamWhat())
  aln <- scanBam(Rsamtools::BamFile(bamfile), param = params)[[1]]
  coverage <- coverage(IRanges(aln[["pos"]], width=aln[["qwidth"]]))
  coverage <- as.numeric(coverage)[GenomicRanges::start(region):GenomicRanges::end(region)]
  return(coverage)
}

# create the layout from the BAM info
#
# @param BAM list path to the BAM files
# @param region vector of coordiantes of the region, chromosom, start and end
# @return returns a vector with the coverage of one file
############################################################################################################################################
##
## makeBAMLayout -- create the layout from the BAM info
##
## BAM    -- list   : list path to the BAM files
##
## region -- vector : coordiantes of the region, chromosom, start and end
##
## returns a plot
############################################################################################################################################
makeBAMLayout <- function(BAM,region){
  return_object <- list()

  for(i in unique(BAM$group)){
    return_object[[i]] <-BAM[which(BAM$group == i),"files"]
  }
  return(return_object)
}

#  plot the average of each group
#
# @param coverages list of coverages for all the groups
# @param region vector of coordiantes of the region, chromosom, start and end
# @param conditions vector of names of the different conditions
# @param colorPalette vector that contain the colors to be used for each line

############################################################################################################################################
##
## doPlotAverageBAM -- plot the average of each group
##
## coverages    -- coverages   : list of coverages for all the groups
##
## region       -- vector : coordiantes of the region, chromosom, start and end
##
## conditions   -- vector : Names of the different conditions
##
## colorPalette -- vector : Contain the colors to be used for each line
##
## returns a plot
############################################################################################################################################
plotAverageBAM <- function(region,coverages,conditions,colorPalette){
  average_cov <- data.frame(GenomicRanges::start(region):GenomicRanges::end(region))
  for(i in 1:length(coverages)){
    coverage <- cbind(as.data.frame(coverages[[i]]))
    average_cov <- cbind(average_cov,rowMeans(coverage))
  }
  max <- max(as.numeric(names(table(unlist(average_cov[,-1])))))
  colnames(average_cov) <- c("position",sprintf("condition%02d", seq(1,length(coverages))))
  graphics::plot(1,xlim=c(GenomicRanges::start(region),GenomicRanges::end(region)),ylim=c(0,max), main="",ylab="Average coverage",xlab="",xaxt="n")
  for(i in 2:length(average_cov)){
    graphics::lines(average_cov[,"position"],average_cov[,i],col=colorPalette[i],lwd = 2)
  }

}

# plot the bam regions
#
# @param coverages list of coverages for all the groups
# @param region vector of coordiantes of the region, chromosom, start and end
# @param condition character name of the condition plotted
# @param colorPalette vector that contain the colors to be used for each line
# @param max numeric Ordinate maximum for ylim
############################################################################################################################################
##
## doPlotBAM -- plot the bam regions
##
## coverages  -- coverages   : list of coverages for one group
##
## region     -- vector : coordiantes of the region, chromosom, start and end
##
## conditions -- character: Name of the condition
##
## colorPalette -- vector : Contain the colors to be used for each line
##
## returns a plot
############################################################################################################################################
plotBAM <- function(region,coverages, max, condition,colorPalette){
  coverage <- cbind(as.data.frame(coverages),GenomicRanges::start(region):GenomicRanges::end(region))
  colnames(coverage) <- c(sprintf("rep%02d", seq(1,length(coverages))),"position")
  graphics::plot(1,xlim=c(GenomicRanges::start(region),GenomicRanges::end(region)),ylim=c(0,max), main="",ylab=condition,xlab="",xaxt="n")
  for(i in 1:(length(coverage)-1)){
    graphics::lines(coverage[,"position"],coverage[,i],col=colorPalette[i])
  }
}
