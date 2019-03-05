############################################################################################################################################
##
## plotRegiont.R -- plot a snapshot of a genomic region
## 30 november 2018
## Thomas Faux
## Medical Bioinformatics Centre
############################################################################################################################################

## createDataObject -- load the data from the files
## unitTest -- check if the files are present, if not returns an error message
## plotRegion -- main function to plot the region
## gg_color_hue -- create a color palette similar to ggplot2
## defineColorPalettes -- create a palette for each plot
## defineLayout -- create a matrix for the layout function

#' @importFrom rbamtools getRefCoords bamReader bamCount
#' @import GenomicRanges
#' @importFrom biomaRt getBM useMart
#' @importFrom IRanges elementNROWS splitAsList
#' @importFrom S4Vectors List subjectHits subjectLength queryHits
#' @importFrom grDevices hcl
#' @importFrom graphics axis legend lines par plot plot.new rect segments text
#' @importFrom utils read.table
NULL

#' Plot a genomic region
#'
#' @param region a GRange object with chr, start, end
#' @param genome a character vector "hg19","hg38" or "mm10"
#' @param BAM a path to the BAM related csv input file
#' @param BED a path to the BED related csv input file
#' @param avgTrack a logical indicating if the average track should be present or not
#' @param geneTrack a logical indicating if the gene track should be present or not
#' @param max a vector of numbers containing the yaxis maximum value of each BAM track
#' @param verbose prompt the progress of the plotting
#'
#' @examples
#' region <- GRanges("chr12:110938000-110940000")
#' setwd(system.file("extdata", package = "RepViz"))
#' RepViz::RepViz(region = region,
#'                genome = "hg19",
#'                BAM = "BAM_input.csv",
#'                BED = "BED_input.csv",
#'                avgTrack = TRUE,
#'                geneTrack = TRUE)
#'
#' @export "RepViz"
RepViz <- function(region,genome,BAM=NULL,BED=NULL,avgTrack=TRUE,geneTrack=TRUE,max=NULL,verbose=TRUE){

  #loading the files
  object_holder <- createDataObject(BAM=BAM,BED=BED)
  plots_list <- list()
  unitTest(object_holder)
  colorsPalette <- defineColorPalettes(object_holder)
  #ploting the coverages
  mat <- defineLayout(object_holder,geneTrack,avgTrack)
  layout(mat,widths=c(5,1))
  graphics::par(mar=c(2, 4, 1, 1))
  if(!is.null(BAM)){
    if(verbose == TRUE){
      cat("plotting the coverages \n")
    }
    layout <- makeBAMLayout(object_holder$BAM,region)
    coord <- c(as.character(GenomicRanges::seqnames(region)),GenomicRanges::start(region),GenomicRanges::end(region))
    coverages <- makeCoverages(coord,layout)

    if(is.null(max)){
      max <- rep(max(unlist(coverages)),max(unique(as.numeric(object_holder$BAM$group))))
    }
    for(group in 1:max(unique(as.numeric(object_holder$BAM$group)))){

      plotBAM(coord,coverages[[group]],max[group],unique(object_holder$BAM$condition)[group],colorsPalette[[1]])
    }
    if(avgTrack == TRUE){
      plotAverageBAM(coord,coverages,unique(object_holder$BAM$condition),colorsPalette[[2]])
    }

  }
  if("BED" %in% names(object_holder)){
    if(verbose == TRUE){
      cat("plotting the BED files \n")
    }
    plotBED(BED = object_holder$BED,region = region,colorsPalette[[3]],verbose)
  }

  if(geneTrack == TRUE){
    if(verbose == TRUE){
      cat("plotting the gene track \n")
    }
    bm <- getBiomaRt(region,genome)
    UTR5 <- findUTR5(region,bm)
    UTR3 <- findUTR3(region,bm)
    gr <- findGenes(region,bm)
    plotGenomicTrack(gr = gr,UTR3 = UTR3, UTR5 = UTR5,region)
  }
  graphics::plot.new()
  legend("center",legend = 1:replicatesNumber(object_holder), col = colorsPalette[[1]],lty=1,lwd = 2,box.col = "white",bg = "white")
  if(avgTrack){
    graphics::plot.new()
    graphics::legend("center",legend = unique(object_holder$BAM$condition), col = colorsPalette[[2]][2:length(colorsPalette[[2]])] ,lty=1,lwd = 2,box.col = "white",bg = "white")
  }
  if(!is.null(BED)){
    graphics::plot.new()
    graphics::legend("center",legend = rev(object_holder$BED$software), col = rev(colorsPalette[[3]][1:length(object_holder$BED$software)]),lty=1,lwd = 2,box.col = "white",bg = "white")
  }

  graphics::par(mfrow=c(1,1))
}


# load the informations from the files
#
# @param BAM a path to the BAM files input file
# @param BED a Path to the BED files input file
# @return returns an object containing the info in the csv files
############################################################################################################################################
##
## createDataObject -- load the informations from the files
##
## BAM   -- path   : Path to the input file of BAM files
## BED   -- path   : Path to the input file of BED files
##
##
## returns an object containing the info in the csv files
############################################################################################################################################
createDataObject <- function(BAM=NULL, BED=NULL){
  #takes into input the pah of the file with the bam and bed
  return_object <- list()

  if(!is.null(BAM)){
    cat(paste0("loading the BAM related data from ",BAM , "\n"))
    files <- utils::read.table(file =BAM ,sep = ";",colClasses = c("character","numeric","character"))
    colnames(files) <- c("files","group","condition")
    return_object$BAM <- files
  }
  if(!is.null(BED)){
    cat(paste0("loading the BED related data from ",BED ,"\n"))
    files <- utils::read.table(file =BED ,sep = ";",colClasses = c("character","character"))
    colnames(files) <- c("files","software")
    return_object$BED <- files
  }

  return(return_object)
}

# create a color palette similar to ggplot2
#
# @param n integer representing number of hues
# @param alpha double for transparancy, ranges from 0.0 to 1.0
# @return returns a list of hues
############################################################################################################################################
##
## gg_color_hue -- create a color palette similar to ggplot2
##
## n      -- integer  : number of hues
## alpha  -- double   : transparancy, ranges from 0.0 to 1.o
##
##
## returns a list of hues
############################################################################################################################################

gg_color_hue <- function(n,alpha) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}

# create a palette for each plot
#
# @param object info contained in the csv input files
# @return returns a colorPalette for each plot
############################################################################################################################################
##
## defineColorPalettes -- create a palette for each plot
##
## object -- object  : info contained in the csv input files
##
## returns a colorPalette for each plot
############################################################################################################################################
defineColorPalettes <- function(object){
  if("BAM" %in% names(object)){
    nbrep <- replicatesNumber(object)
    BAM <- gg_color_hue(nbrep)
    BAMavg <- gg_color_hue(max(unique(object$BAM$group))+1,0.5)
  }
  else{
    BAM<-NULL
    BAMavg<-NULL
  }
  if("BED" %in% names(object)){
    BED <- gg_color_hue(length(object$BED$software)+ceiling(length(object$BAM$condition))/2)
    #BED <- BED[-c(1:ceiling(length(object$BAM$condition))/2)]
  }
  else{
    BED<-NULL
  }

  colorPalettes <- S4Vectors::List(BAM,BAMavg,BED)
  return(colorPalettes)
}

# create a matrix for the layout function
#
# @param object_holder info contained in the csv input files
# @param geneTrack logical if the gene track present or not
# @param avgTrack logical if the average track present or not
# @return returns a matrix to input in the layout function
############################################################################################################################################
##
## defineLayout -- create a matrix for the layout function
##
## object -- object  : info contained in the csv input files
##
## geneTrack -- boolean : is the gene track present or not
##
## returns a matrix to input in the layout function
############################################################################################################################################
defineLayout <- function(object_holder,geneTrack,avgTrack){
  groups <- max(unique(object_holder$BAM$group))
  if(geneTrack){
    mat <- matrix(c(1:(groups+1),groups+2,groups+3,rep(x=groups+4,groups),groups+5,groups+6,groups+7), groups+3, 2, byrow = FALSE)
  }
  if(!avgTrack){
    mat <- matrix(c(1:(groups+1),groups+2,rep(x=groups+3,groups),groups+4,groups+5), groups+2, 2, byrow = FALSE)
  }
  if(!geneTrack){
    mat <- matrix(c(1:(groups+1),groups+2,rep(x=groups+3,groups),groups+4,groups+5), groups+2, 2, byrow = FALSE)
  }
  if(!"BED" %in% names(object_holder)){
    mat <- matrix(c(1:(groups+1),groups+2,rep(x=groups+3,groups),groups+4,groups+5), groups+2, 2, byrow = FALSE)
  }
  return(mat)
}

# check if the files are present, if not returns an error message
#
# @param object_holder info contained in the csv input files

############################################################################################################################################
##
## unitTest -- check if the files are present, if not returns an error message
##
## object -- object  : info contained in the csv input files
##
############################################################################################################################################
unitTest <- function(object_holder){
  for(file in object_holder$BAM$files){
    if(!file.exists(file)){
      cat(paste0("The file : ", file," is not present in the indicated folder \n"))
    }
  }
  for(file in object_holder$BED$files){
    if(!file.exists(file)){
      cat(paste0("The file : ", file," is not present in the indicated folder \n"))
    }
  }
  if(!"BED" %in% names(object_holder)){
    cat(paste0("There is no hits found in the BED files for this region \n"))
  }
}

# calculate the maximum number of replicates
#
# @param object_holder info contained in the csv input files
# @return The maximum number of replicates in the groups

############################################################################################################################################
##
## replicatesNumber -- calculate the maximum number of replicates
##
## object_holder -- object  : info contained in the csv input files
##
############################################################################################################################################
replicatesNumber <- function(object_holder){
  groupnb <- unique(object_holder$BAM$group)

  vec <- c()
  for(i in 1:max(groupnb)){
    len <- dim(object_holder$BAM[which(object_holder$BAM$group == i),])[1]
    vec <- c(vec,len)
  }

  return(max(vec))

}
