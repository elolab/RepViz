###############################################################################
## plotRegiont.R -- plot a snapshot of a genomic region
## 30 november 2018 Thomas Faux
## Medical Bioinformatics Centre
###############################################################################

#' @import Rsamtools
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
#' @param genome a character vector 'hg19','hg38' or 'mm10'
#' @param BAM a path to the BAM related csv input file
#' @param BED a path to the BED related csv input file
#' @param avgTrack a logical indicating if the average track should be included or not
#' @param geneTrack a logical indicating if the gene track should be included or not
#' @param max a numerical vector containing the yaxis maximum value of each BAM track
#' @param verbose a logical indicating whether the progress of the plotting is shown
#' @param cex number indicating the amount by which plotting text and symbols should be scaled relative to the default. 
#' @param col vector of character user can set color of the different BED tracks
#' @return displays the region specified by the user
#'
#'
#' @examples
#' region <- GRanges('chr12:110938000-110940000')
#' setwd(tempdir())
#' #Copying the files to the user working directory
#' file.copy(from = list.files(system.file('extdata', package = 'RepViz'), full.names = TRUE),
#'     to = tempdir())
#' #Generate the visualization of the given region
#' RepViz::RepViz(region = region,
#'     genome = 'hg19',
#'     BAM = 'BAM_input.csv',
#'     BED = 'BED_input.csv',
#'     avgTrack = TRUE,
#'     geneTrack = TRUE)
#'
#' @export 'RepViz'
RepViz <- function(region, genome=c('hg19','hg38','mm10'), BAM = NULL, BED = NULL, avgTrack = TRUE, geneTrack = TRUE,
    max = NULL, verbose = TRUE, cex = 1,col = NULL) {

    # loading the files
    object_holder <- createDataObject(BAM = BAM, BED = BED, verbose = verbose)
    plots_list <- list()
    unitTest(object_holder,region,genome)
    colorsPalette <- defineColorPalettes(object_holder)
    # ploting the coverages
    mat <- defineLayout(object_holder, geneTrack, avgTrack)
    layout(mat, widths = c(5, 1))
    graphics::par(mar = c(2, 6, 1, 0))
    if (!is.null(BAM)) {
        if (verbose == TRUE) {
            message("plotting the coverages \n")
        }
        layout <- makeBAMLayout(object_holder$BAM)
        coverages <- makeCoverages(region, layout)
        if (is.null(max)) {
            max <- rep(max(as.numeric(names(table(unlist(coverages))))), length(unique(object_holder$BAM$group)))
        }

        for (group in seq_len(length(unique(object_holder$BAM$group)))) {
            plotBAM(region, coverages[[group]], max[group], unique(object_holder$BAM$group)[group],
                colorsPalette[[1]], cex)
        }

        if (avgTrack == TRUE) {
            plotAverageBAM(region, coverages, unique(object_holder$BAM$group), colorsPalette[[2]], cex)
        }

    }
    if ("BED" %in% names(object_holder)) {
        if (verbose == TRUE) {
            message("plotting the BED files \n")
        }
      if(!is.null(col)){
        colors <- col
      }else{
        colors <- colorsPalette[[3]]
      }
        plotBED(BED = object_holder$BED, region = region, colors, verbose)
    }

    if (geneTrack == TRUE) {
        if (verbose == TRUE) {
            message("plotting the gene track \n")
        }
        bm <- getBiomaRt(region, genome)
        UTR5 <- findUTR5(region, bm)
        UTR3 <- findUTR3(region, bm)
        gr <- findGenes(region, bm)
        plotGenomicTrack(gr = gr, UTR3 = UTR3, UTR5 = UTR5, region, cex)
    }
    graphics::par(mar = c(1, 1, 1, 1))
    graphics::plot.new()
    plotLegends(object_holder = object_holder, colorsPalette = colorsPalette, avgTrack = avgTrack, BED = BED, cex, col)
    graphics::par(mfrow = c(1, 1))
}

plotLegends <- function(object_holder, colorsPalette, avgTrack, BED, cex, col) {
    legend("left", legend = seq_len(replicatesNumber(object_holder)),
        col = colorsPalette[[1]], lty = 1, lwd = cex, box.col = "white", bg = "white", cex = cex,  text.font  = 2)

    if (avgTrack) {
        graphics::plot.new()
        graphics::legend("left", legend = unique(object_holder$BAM$group),
                        col = colorsPalette[[2]][2:length(colorsPalette[[2]])],
                        lty = 1, lwd = cex, box.col = "white", bg = "white", cex = cex,  text.font  = 2)
    }

    if (!is.null(BED)) {
      if(!is.null(col)){
        colors <- col[length(col):1]
      }else{
        colors <- rev(colorsPalette[[3]][seq_len(length(object_holder$BED$software))])
      }
        graphics::plot.new()
        graphics::legend("left", legend = rev(object_holder$BED$software),
                        col = colors,
                        lty = 1, lwd = cex, box.col = "white", bg = "white", cex = cex,  text.font  = 2)
    }

}
# load the informations from the files
# @param BAM a path to the BAM files input file
# @param BED a Path to the BED files input file
# @return returns an object containing the info in the csv files

createDataObject <- function(BAM = NULL, BED = NULL, verbose = TRUE) {
    # takes into input the pah of the file with the bam and bed
    return_object <- list()

    if (!is.null(BAM)) {
        if (verbose == TRUE)  message("loading the BAM related data from ", BAM, "\n")
      
        if(is.character(BAM)){
            files <- utils::read.table(file = BAM, sep = ",", colClasses = c("character", "character"))
            colnames(files) <- c("files", "group")
            return_object$BAM <- files
        } else if(is.data.frame(BAM)){
            return_object$BAM <- BAM
        } else {
            stop("BAM is not character or data.frame")
        }
      
    }
    if (!is.null(BED)) {
        if (verbose == TRUE) message("loading the BED related data from ", BED, "\n")
      
        if(is.character(BAM)){
            files <- utils::read.table(file = BED, sep = ",", colClasses = c("character", "character"))
            colnames(files) <- c("files", "software")
            return_object$BED <- files
        } else if(is.data.frame(BED)){
            return_object$BED <- BED
        } else {
            stop("BED is not character or data.frame")
        }
    }

    return(return_object)
}

# create a color palette similar to ggplot2
# @param n integer representing number of hues
# @param alpha double for transparancy, ranges from 0.0 to 1.0
# @return returns a list of hues

gg_color_hue <- function(n, alpha) {
    hues = seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100, alpha = alpha)[seq_len(n)]
}

# create a palette for each plot
# @param object info contained in the csv input files
# @return returns a colorPalette for each plot

defineColorPalettes <- function(object) {
    if ("BAM" %in% names(object)) {
        nbrep <- replicatesNumber(object)
        BAM <- gg_color_hue(nbrep)
        BAMavg <- gg_color_hue(length(unique(object$BAM$group)) + 1, 0.5)
    } else {
        BAM <- NULL
        BAMavg <- NULL
    }
    if ("BED" %in% names(object)) {
        BED <- gg_color_hue(length(object$BED$software) + ceiling(length(object$BAM$condition))/2)
    } else {
        BED <- NULL
    }

    colorPalettes <- S4Vectors::List(BAM, BAMavg, BED)
    return(colorPalettes)
}

# create a matrix for the layout function
# @param object_holder info contained in the csv input files
# @param geneTrack logical if the gene track present or not
# @param avgTrack logical if the average track present or not
# @return returns a matrix to input in the layout function

defineLayout <- function(object_holder, geneTrack, avgTrack) {
    groups <- length(unique(object_holder$BAM$group))
    
    A <- as.numeric(geneTrack)
    B <- as.numeric(avgTrack)
    C <- as.numeric("BED" %in% names(object_holder))
    sum <- A+B+C
    
    switch(as.character(sum),
        "0" = {
            mat <- matrix(c(seq_len(groups) , rep(x = groups + 1, groups)), groups, 2, byrow = FALSE)
        },
        "1" = {
            mat <- matrix(c(seq_len(groups), groups + 1 , rep(x = groups + 2, groups),
                          groups + 3 ), groups + 1, 2, byrow = FALSE) 
        },
        "2" = {
            mat <- matrix(c(seq_len(groups), groups + 1 ,groups + 2, rep(x = groups + 3, groups),
                          groups + 4, groups + 5 ), groups + 2, 2, byrow = FALSE)
        },
        "3" = {
            mat <- matrix(c(seq_len(groups), groups + 1 ,groups + 2, groups + 3, rep(x = groups + 4, groups),
                          groups + 5, groups + 6, groups + 7), groups + 3, 2, byrow = FALSE)
        }
      
    )

    return(mat)
}

# check if the files are present, if not returns an error message
# @param object_holder info
# contained in the csv input files

unitTest <- function(object_holder,region,genome) {
    for (file in object_holder$BAM$files) {
        if (!file.exists(file)) {
            stop("The file : ", file, " is not present in the indicated folder \n")
        }
    }
    for (file in object_holder$BED$files) {
        if (!file.exists(file)) {
            stop("The file : ", file, " is not present in the indicated folder \n")
        }
    }
    if (!"BED" %in% names(object_holder)) {
        message("There is no hits found in the BED files for this region \n")
    }
    if (!is(object = region, class2 = "GRanges")){
        message("The given region is not a GRanges object \n")
    }
    stopifnot(genome %in% c('hg19','hg38','mm10'))
}

# calculate the maximum number of replicates
# @param object_holder info contained in the csv input files
# @return The maximum number of replicates in the groups

replicatesNumber <- function(object_holder) {
    groupnb <- unique(object_holder$BAM$group)

    vec <- numeric(length(groupnb))
    for (i in groupnb) {
        len <- dim(object_holder$BAM[which(object_holder$BAM$group == i), ])[1]
        vec[i] <- len
    }

    return(max(vec))

}
