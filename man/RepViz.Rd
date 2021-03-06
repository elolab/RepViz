% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotRegion.R
\name{RepViz}
\alias{RepViz}
\title{Plot a genomic region}
\usage{
RepViz(region, genome = c("hg19", "hg38", "mm10"), BAM = NULL,
  BED = NULL, avgTrack = TRUE, geneTrack = TRUE, max = NULL,
  verbose = TRUE, cex = 1)
}
\arguments{
\item{region}{a GRange object with chr, start, end}

\item{genome}{a character vector 'hg19','hg38' or 'mm10'}

\item{BAM}{a path to the BAM related csv input file}

\item{BED}{a path to the BED related csv input file}

\item{avgTrack}{a logical indicating if the average track should be included or not}

\item{geneTrack}{a logical indicating if the gene track should be included or not}

\item{max}{a numerical vector containing the yaxis maximum value of each BAM track}

\item{verbose}{a logical indicating whether the progress of the plotting is shown}

\item{cex}{number indicating the amount by which plotting text and symbols should be scaled relative to the default.}

\item{col}{ vector of character user can set color of the different BED tracks.}
}
\value{
displays the region specified by the user
}
\description{
Plot a genomic region
}
\examples{
region <- GRanges('chr12:110938000-110940000')
setwd(tempdir())
#Copying the files to the user working directory
file.copy(from = list.files(system.file('extdata', package = 'RepViz'), full.names = TRUE),
    to = tempdir())
#Generate the visualization of the given region
RepViz::RepViz(region = region,
    genome = 'hg19',
    BAM = 'BAM_input.csv',
    BED = 'BED_input.csv',
    avgTrack = TRUE,
    geneTrack = TRUE)

}
