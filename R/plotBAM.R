############################################################################################################################################
##
## plotRegiont.R -- plot a snapshot of a genomic region
## 30 november 2018
## Thomas Faux
## Medical Bioinformatics Centre
############################################################################################################################################

## makeCoverages -- produce the coverage for each replicate

## makeBAMLayout -- create the layout from the BAM info

## getCoverage -- create the readerfor rbamtools

## coverageGenomicRegion -- call the rbamtools coverage calculation

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

# create the reader for rbamtools
#
# @param bamfile character path to the BAM file
# @param region GRanges object containing the region to plot
# @return returns the coverage of one file
############################################################################################################################################
##
## getCoverage -- create the readerfor rbamtools
##
## bamfile  -- character : path to the BAM file
##
## region   -- vector    : coordiantes of the region, chromosom, start and end
##
## returns the coverage of one file
############################################################################################################################################
getCoverage <- function(region, bamfile){
  x <- GenomicAlignments::coverage(Rsamtools::BamFile(bamfile))
  coverage <- as.numeric(x[[GenomicRanges::seqnames(region)]][GenomicRanges::ranges(region)])
  return(coverage)
}

# # call the rbamtools coverage calculation
# #
# # @param coord vector of coordiantes of the region, chromosom, start and end
# # @param reader  rbamtools reader
# # @return returns a vector with the coverage of one file
# ############################################################################################################################################
# ##
# ## coverageGenomicRegion -- call the rbamtools coverage calculation
# ##
# ## coord    -- vector : coordiantes of the region, chromosom, start and end
# ##
# ## reader   -- reader : rbamtools reader
# ##
# ## returns a vector with the coverage of one file
# ############################################################################################################################################
# coverageGenomicRegion <- function(coord,reader){
#   vec <- vector()
#   for(i in coord[2]:coord[3]){
#     coords<-rbamtools::getRefCoords(reader,coord[[1]])
#     coords[2]<-i-1
#     coords[3]<-i
#     count<-rbamtools::bamCount(reader,coords)
#     vec <- c(vec,count[10])
#   }
#   return(vec)
# }

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
  max <- max(unlist(average_cov[,-1]))
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
