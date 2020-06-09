###############################################################################
## plotRegiont.R -- plot a snapshot of a genomic region
## 30 november 2018 Thomas Faux
## Medical Bioinformatics Centre
###############################################################################

# transform the dataframe in GRangesList
# @param df list of dataframes containing the peaks
# @return returns a GRangesList

toGRList <- function(df) {
    grlist <- GenomicRanges::GRangesList()
    for (i in unique(df$external_gene_name)) {
        temp <- GenomicRanges::makeGRangesFromDataFrame(df[which(df$external_gene_name == i),
            ])
        grlist[[i]] <- temp
    }
    names(grlist) <- unique(df$external_gene_name)
    return(grlist)
}

# find genes overlaping with the region for the given database
# @param region GRanges object containing the region to plot
# @param genome Genome used by the db object 'hg19','GRCh38' or 'mm10'
# @param m biomaRt object @return returns GRanges with the Genes hits found in the region

findGenes <- function(region, m) {
    filters <- c("chromosome_name", "start", "end")
    chr <- substr(as.character(seqnames(region)), 4,
        nchar(as.character(seqnames(region))))
    values <- list(chr, GenomicRanges::start(region),
        GenomicRanges::end(region))
    map <- biomaRt::getBM(mart = m, attributes = c("ensembl_gene_id",
                                                    "external_gene_name",
                                                    "ensembl_exon_id",
                                                    "chromosome_name",
                                                    "exon_chrom_start",
                                                    "exon_chrom_end",
                                                    "strand"),
                                                    filters = filters,
                                                    values = values)

    map[which(map$strand == -1), "strand"] <- "-"
    map[which(map$strand == 1), "strand"] <- "+"
    if(dim(map)[1] > 0){
      map$chromosome_name <- as.character(GenomicRanges::seqnames(region))  
    }
    gr <- GenomicRanges::makeGRangesFromDataFrame(map,
                                        keep.extra.columns = TRUE)
    return(gr)
}

# Helper function for the two find UTR functions

f <- function(UTR,region) {
  return_object <- GenomicRanges::GRanges()
  for (gene in unique(UTR$external_gene_name)) {
    temp <- UTR[which(UTR$external_gene_name == gene), ]
    temp <- GenomicRanges::makeGRangesFromDataFrame(temp, keep.extra.columns = TRUE)
    temp <- temp[queryHits(GenomicRanges::findOverlaps(temp, region, ignore.strand = TRUE))]
    if (length(temp) > 0) {
      temp <- temp[which(temp$transcript_length == max(temp$transcript_length))]
    }
    return_object <- c(return_object, temp)
  }
  return(return_object)
}

# find UTR overlaping with the region for the given database
# @param region GRanges object containing the region to plot
# @param m biomaRt object @return returns GRanges with the UTR hits found in the region

findUTR5 <- function(region, m) {
    filters <- c("chromosome_name", "start", "end")
    chr <- substr(as.character(seqnames(region)), 4, nchar(as.character(seqnames(region))))
    values <- list(chr, GenomicRanges::start(region), GenomicRanges::end(region))
    map <- biomaRt::getBM(mart = m, attributes = c("external_gene_name", "chromosome_name", "strand",
        "5_utr_start", "5_utr_end", "transcript_length"), filters = filters, values = values)

    map[which(map$strand == -1), "strand"] <- "-"
    map[which(map$strand == 1), "strand"] <- "+"

    UTR5 <- map[!is.na(map$`5_utr_start`), ]
    if (dim(UTR5)[1] > 0) {
        UTR5$chromosome_name <- as.character(GenomicRanges::seqnames(region))
        UTR5 <- f(UTR5,region)

    } else {
        UTR5 <- GenomicRanges::GRanges()
    }
    return(UTR5)
}


findUTR3 <- function(region, m) {
    filters <- c("chromosome_name", "start", "end")
    chr <- substr(as.character(seqnames(region)), 4, nchar(as.character(seqnames(region))))
    values <- list(chr, GenomicRanges::start(region), GenomicRanges::end(region))
    map <- biomaRt::getBM(mart = m, attributes = c("external_gene_name", "chromosome_name", "strand",
        "3_utr_start", "3_utr_end", "transcript_length"), filters = filters, values = values)

    map[which(map$strand == -1), "strand"] <- "-"
    map[which(map$strand == 1), "strand"] <- "+"

    UTR3 <- map[!is.na(map$`3_utr_start`), ]
    if (dim(UTR3)[1] > 0) {
        UTR3$chromosome_name <- as.character(GenomicRanges::seqnames(region))
        UTR3 <- f(UTR3,region)
    } else {
        UTR3 <- GenomicRanges::GRanges()
    }

    return(UTR3)
}
# get The biomaRt object
# @param region GRanges object containing the region to plot
# @param genome Genome used by the db object 'hg19','GRCh38' or 'mm10'
# @return returns GRanges with the Genes hits found in the region

getBiomaRt <- function(region, genome = c("hg19", "GRCh38", "mm10")) {

  switch(genome, 
         hg19={
             m <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                                 host = "grch37.ensembl.org", path = "/biomart/martservice",
                                 dataset = "hsapiens_gene_ensembl")
         },
         GRCh38={
             m <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                                 host = "grch38.ensembl.org", path = "/biomart/martservice",
                                 dataset = "hsapiens_gene_ensembl")  
         },
         mm10={
             m <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")  
         },
         {
             m <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                                 host = "grch37.ensembl.org", path = "/biomart/martservice",
                                 dataset = "hsapiens_gene_ensembl")
         }
  )

    return(m)
}


# plot the rectangles of the genomic track
# @param rg GRanges list of exons in the gene
# @param gr GRanges object containing the region to plot
# @param y numeric index for color purpose
# @param size value to substract and add around the y value to create a rectangle

plotRectangles <- function(rg, gr, y, size) {
    graphics::rect(GenomicRanges::end(gr), rep(y, length(gr)) - size, GenomicRanges::start(gr),
        rep(y, length(gr)) + size, border = NA, col = "grey")  #y+5)

}

# plot the arrows between the rectangles in the genomic track
# @param x0 list of x start points
# @param y0 list of y start points
# @param x1 list of x end points @param y1 list of y end points
# @param n_arr number of arrows needed
# @param ... graphics::arrows basic parameters

arrowLine <- function(x0, y0, x1, y1, n_arr, ...) {
    x <- seq(x0, x1, length = n_arr + 1)
    y <- seq(y0, y1, length = n_arr + 1)
    graphics::arrows(x[-length(x)], y[-length(y)], x[-1], y[-1], ...)
}

# return the regions contained between the ranges of a GRanges
# @param gr GRanges object
# @return GRanges object

inverseGRanges <- function(gr) {
    temp <- as.data.frame(gr)
    temp <- cbind(as.character(temp[seq_len(dim(temp)[1] - 1), 1]), temp[seq_len(dim(temp)[1] -
        1), 3], temp[2:dim(temp)[1], 2], as.character(temp[seq_len(dim(temp)[1] - 1), 5]))
    temp <- as.data.frame(temp)
    colnames(temp) <- c("seqnames", "start", "end", "strand")
    temp <- GRanges(temp)
    return(temp)
}

# wrapper function to plot the rectangles and arrows of the genomic track
# @param gr GRanges object containing the region to plot
# @param region GRanges object containing the region to plot
# @param y Numeric value used for the y position of the boxes

plotGRanges <- function(gr, region, y, UTR3, UTR5) {
    gr <- GenomicRanges::reduce(gr)

    if (length(gr) > 1) {
        gr_inv <- inverseGRanges(gr)
        for (j in seq_len(length(gr_inv))) {
            value <- (GenomicRanges::end(gr_inv[j]) - GenomicRanges::start(gr_inv[j]))/(GenomicRanges::end(region) -
                GenomicRanges::start(region))
            N <- arrowNumbers(value)
            if (as.character(strand(gr_inv[j])) == "-") {
                arrowLine(GenomicRanges::end(gr_inv[j]), y, GenomicRanges::start(gr_inv[j]), y,
                n_arr = N, length = 0.1)
            }
            if (as.character(strand(gr_inv[j])) == "+") {
                arrowLine(GenomicRanges::start(gr_inv[j]), y, GenomicRanges::end(gr_inv[j]), y,
                n_arr = N, length = 0.1)
            }
        }
    }
    if (length(UTR5) > 0) {
        plotRectangles(region, UTR5, y, size = 0.25)
        if (as.character(unique(strand(gr))) == "-") {
            GenomicRanges::end(
                gr[queryHits(GenomicRanges::findOverlaps(gr, UTR5))]
                ) <- GenomicRanges::start(
                    UTR5[subjectHits(GenomicRanges::findOverlaps(gr,UTR5))])
        }
        if (as.character(unique(strand(gr))) == "+") {
            GenomicRanges::start(
                gr[queryHits(GenomicRanges::findOverlaps(gr, UTR5))]
                ) <- GenomicRanges::end(
                    UTR5[subjectHits(GenomicRanges::findOverlaps(gr,UTR5))])
        }
    }
    if (length(UTR3) > 0) {
        plotRectangles(region, UTR3, y, size = 0.25)
        if (as.character(unique(strand(gr))) == "+") {
            GenomicRanges::end(
                gr[queryHits(GenomicRanges::findOverlaps(gr, UTR3))]
                ) <- GenomicRanges::start(
                    UTR3[subjectHits(GenomicRanges::findOverlaps(gr,UTR3))])
        }
        if (as.character(unique(strand(gr))) == "-") {
            GenomicRanges::start(
                gr[queryHits(GenomicRanges::findOverlaps(gr, UTR3))]
                ) <- GenomicRanges::end(
                    UTR3[subjectHits(GenomicRanges::findOverlaps(gr,UTR3))])
        }
    }
    plotRectangles(region, gr, y, size = 0.5)

}

arrowNumbers <- function(value) {
    N <- 0

    if (value < 0.25) {
        N = 1
    }
    if ((value > 0.25) & (value <= 0.75)) {
        N = 2
    }
    if ((value > 0.75) & (value <= 1)) {
        N = 3
    }
    if ((value > 1) & (value <= 2)) {
        N = 6
    }
    if (value > 2) {
        N = 12
    }
    return(N)
}

# plot the genomic track
# @param gr GRanges object containing the region to plot
# @param region GRanges object containing the region to plot

plotGenomicTrack <- function(gr, UTR3, UTR5, region, cex) {

    if (dim(as.data.frame(gr))[1] > 0) {
        index <- 0
        if (length(unique(gr$external_gene_name)) == 1) {
            graphics::plot(x = 1, ylim = c(0, 2),
                xlim = c(GenomicRanges::start(region),
                GenomicRanges::end(region)),
                xlab = GenomicRanges::seqnames(region), yaxt = "n", ylab = "",
                xaxt = "n")
            graphics::axis(2, at = 1, labels = FALSE)
            graphics::text(x = par()$usr[1] - 0.03 * (par()$usr[2] - par()$usr[1]), 
                           y = 1, labels = unique(gr$external_gene_name),srt = 90, xpd = NA, font=2, cex= cex)
            plotGRanges(gr, region, 1, UTR3, UTR5)
        }
        if (length(unique(gr$external_gene_name)) > 1) {
            graphics::plot(x = 1, ylim = c(0, length(unique(gr$external_gene_name))),
                xlim = c(start(region), end(region)), xlab = seqnames(region),
                yaxt = "n", ylab = "")
            tmin <- 1
            tmax <- length(IRanges::unique(gr$external_gene_name))
            tlab <- seq(tmin, tmax) - 0.5
            lab <- IRanges::unique(gr$external_gene_name)
            graphics::axis(2, at = tlab, labels = FALSE)
            graphics::text(x = par()$usr[1] - 0.01 * (par()$usr[2] - par()$usr[1]), y = tlab,
                labels = lab, srt = 45, xpd = NA, adj = 1, font=2, cex = cex)
            for (i in unique(gr$external_gene_name)) {
                index <- index + 1
                gr2 <- gr[(GenomicRanges::elementMetadata(gr)[, "external_gene_name"] %in% i)]
                UTR3_ind <- UTR3[(GenomicRanges::elementMetadata(UTR3)[, "external_gene_name"] %in%
                    i)]
                UTR5_ind <- UTR5[(GenomicRanges::elementMetadata(UTR5)[, "external_gene_name"] %in%
                    i)]
                plotGRanges(gr2, region, index - 0.5, UTR3 = UTR3_ind, UTR5 = UTR5_ind)
            }

        }
    } else {
        message("There is no genes in the given region :", region)
        graphics::plot(x = 1, ylim = c(0, 2),
                     xlim = c(GenomicRanges::start(region),
                              GenomicRanges::end(region)),
                     xlab = GenomicRanges::seqnames(region), yaxt = "n", ylab = "",
                     xaxt = "n")
    }

}


