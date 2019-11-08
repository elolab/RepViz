###############################################################################
##  plotRegiont.R -- plot a snapshot of a genomic region
##  30 november 2018 Thomas Faux
##  Medical Bioinformatics Centre
###############################################################################

# produce the coverage for each replicate
# @param BAMlayout layout created by the function makeBAMLayout
# @param region GRanges object containing the region to plot
# @return returns a coverages object

makeCoverages <- function(region, BAMlayout) {
    coverages <- list()
    for (i in seq_len(length(BAMlayout))) {
        coverages[[i]] <- list()
        for (j in seq_len(length(BAMlayout[[i]]))) {
            coverages[[i]][[j]] <- getCoverage(region, BAMlayout[[i]][[j]])
        }
    }
    return(coverages)
}

# compute the coverage with Rsamtools
# @param bamfile character path to the BAM file
# @param region GRanges object containing the region to plot
# @return returns the coverage of one file

getCoverage <- function(region, bamfile) {

    params <- ScanBamParam(which = region, what = scanBamWhat())
    aln <- scanBam(Rsamtools::BamFile(bamfile), param = params)[[1]]
    coverage <- coverage(IRanges::IRanges(aln[["pos"]], width = aln[["qwidth"]]))
    coverage <- as.numeric(coverage)[GenomicRanges::start(region):GenomicRanges::end(region)]
    coverage[is.na(coverage)] <- 0
    return(coverage)
}

# create the layout from the BAM info
# @param BAM list path to the BAM files
# @return returns a vector with

makeBAMLayout <- function(BAM) {
    return_object <- list()

    for (i in unique(BAM$group)) {
        return_object[[i]] <- BAM[which(BAM$group == i), "files"]
    }
    return(return_object)
}

# plot the average of each group
# @param coverages list of coverages for all the groups
# @param region vector of coordiantes of the region, chromosom, start and end
# @param conditions vector of names of the different conditions
# @param colorPalette vector that contain the colors to be used for each line

plotAverageBAM <- function(region, coverages, conditions, colorPalette, cex) {
    average_cov <- data.frame(GenomicRanges::start(region):GenomicRanges::end(region))
    for (i in seq_len(length(coverages))) {
        coverage <- cbind(as.data.frame(coverages[[i]]))
        average_cov <- cbind(average_cov, rowMeans(coverage))
    }
    max <- max(as.numeric(names(table(unlist(average_cov[, -1])))))
    colnames(average_cov) <- c("position", sprintf("condition%02d", seq(1, length(coverages))))
    graphics::plot(1, xlim = c(GenomicRanges::start(region), GenomicRanges::end(region)), ylim = c(0,
        max), main = "", ylab = "", xlab = "", xaxt = "n",font=2,cex = cex)
    graphics::mtext(side=2, line=3, "Average coverage", col="black", font=2, cex=cex)
    for (i in 2:length(average_cov)) {
        graphics::lines(average_cov[, "position"], average_cov[, i], col = colorPalette[i], lwd = 2)
    }

}

# plot the bam regions @param coverages list of coverages for all the groups
# @param region vector of coordiantes of the region, chromosom, start and end
# @param condition character name of the condition plotted
# @param colorPalette vector that contain the colors to be used for each line
# @param max numeric Ordinate maximum for ylim

plotBAM <- function(region, coverages, max, condition, colorPalette, cex) {
    coverage <- cbind(as.data.frame(coverages), GenomicRanges::start(region):GenomicRanges::end(region))
    colnames(coverage) <- c(sprintf("rep%02d", seq(1, length(coverages))), "position")
    graphics::plot(1, xlim = c(GenomicRanges::start(region), GenomicRanges::end(region)), ylim = c(0,
        max), main = "", ylab = "", xlab = "", xaxt = "n",font=2,cex = cex)
    graphics::mtext(side=2, line=3, condition, col="black", font=2, cex=cex)
    for (i in seq_len((length(coverage) - 1))) {
        graphics::lines(coverage[, "position"], coverage[, i], col = colorPalette[i])
    }
}
