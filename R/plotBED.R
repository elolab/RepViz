###############################################################################
##
## plotRegiont.R -- plot a snapshot of a genomic region
## 30 november 2018
## Thomas Faux
## Medical Bioinformatics Centre
###############################################################################

## loadPeaks -- load the peaks from the bed files listed in the info

## takePeaksOverlap -- take the peaks that are overlaping with the region to plot

## plotBED -- plot the regions that are contained in BED


# plot the regions that are contained in BED
#
# @param BED  list of dataframes containing the peaks
# @param region GRanges object containing the region to plot
# @param colorPalette Contain vector of colors to be used for each line
# @param verbose Prompt the progress
###############################################################################
##
## plotBED -- plot the regions that are contained in BED
##
## BED    -- list    : list of dataframes containing the peaks
##
## region -- GRanges : GRanges object containing the region to plot
##
## colorPalette -- vector : Contain the colors to be used for each line
##
## returns a plot
###############################################################################

plotBED <- function(BED,region,colorPalette,verbose){

  overlaps <- takePeaksOverlap(BED,region,verbose)
  mysegment <- function(GRanges,y,colorPalette){
    graphics::segments(x0 = GenomicRanges::start(GRanges[[1]]),y0=y-0.5,x1=GenomicRanges::end(GRanges[[1]]),y1 = y-0.5, col = colorPalette[y], lwd = 6)
  }

  graphics::plot(x=1, ylim=c(0,length(overlaps)),xlim=c(GenomicRanges::start(region),GenomicRanges::end(region)), ylab ="",yaxt="n")
  if(dim(GenomicRanges::as.data.frame(overlaps))[1]>0){
    index <- 0
    for(i in names(overlaps)){

      index <- index+1
      if(dim(GenomicRanges::as.data.frame(overlaps[i]))[1]==1){
        mysegment(overlaps[i],index,colorPalette)
      }
      if(dim(GenomicRanges::as.data.frame(overlaps[i]))[1]>1){
        mysegment(overlaps[i][1],index,colorPalette)

      }

    }
  }
  else{
    print(paste0("There is no overlaps with the given region :",region))
  }
}

# take the peaks that are overlaping with the region to plot
#
# @param BED  list of dataframes containing the peaks
# @param region GRanges object containing the region to plot
# @param verbose Prompt the progress
# @return returns a list of dataframes containing the peaks
###############################################################################
##
## takePeaksOverlap -- take the peaks that are overlaping with the region to plot
##
## BED    -- list    : list of dataframes containing the peaks
##
## region -- GRanges : GRanges object containing the region to plot
##
## returns a list of dataframes containing the peaks
###############################################################################
takePeaksOverlap <-function(BED,region,verbose){
  if(verbose == TRUE){
    cat("overlapping peaks with region \n")
  }
  peaks <- loadPeaks(BED = BED,verbose)

  return_object <- GenomicRanges::GRangesList()
  for(software in names(peaks)){
    temp <- peaks[[software]][c(1,2,3)]
    colnames(temp) <- c("seqname","start","end")
    temp <- GenomicRanges::GRanges(temp)
    return_object[[software]] <-temp[S4Vectors::queryHits(GenomicRanges::findOverlaps(temp,region))]
  }

  return(return_object)
}

# load the peaks from the bed files listed in the info
#
# @param BED  list of dataframes containing the peaks
# @param verbose Prompt the progress
# @return returns a list of dataframes containing the peaks
###############################################################################
##
## loadPeaks -- load the peaks from the bed files listed in the info
##
## BED   -- list   : list of path to bed files
##
## returns a list of dataframes containing the peaks
###############################################################################
loadPeaks <- function(BED,verbose){
  if(verbose == TRUE){
    cat("load the peak files \n")
  }
  temp <- lapply(X = BED$files, FUN = function(x) utils::read.table(file=x,sep = "\t"))
  names(temp) <- BED$software
  return(temp)
}

