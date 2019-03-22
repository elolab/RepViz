## ----setup, include=FALSE------------------------------------------------
suppressMessages(library(GenomicRanges))
library(knitr)
knitr::opts_chunk$set(echo = TRUE)
opts_knit$set(root.dir = system.file("extdata", package = "RepViz"))

## ----download BAM,eval = FALSE-------------------------------------------
#  file.copy(from = list.files(system.file("extdata", package = "RepViz"), full.names = TRUE),to = getwd())

## ----input BAM, echo=FALSE-----------------------------------------------
BAM_table <- read.table(system.file("extdata","BAM_input.csv", package = "RepViz"), sep = ",")
colnames(BAM_table) <- c("bam file","group")
knitr::kable(BAM_table)

## ----input BED, echo=FALSE-----------------------------------------------
BED_table <- read.table(system.file("extdata","BED_input.csv", package = "RepViz"), sep = ",")
colnames(BED_table) <- c("bed file","Legend")
knitr::kable(BED_table)

## ----bam files, echo=FALSE-----------------------------------------------
list.files(system.file("extdata", package = "RepViz"),pattern = ".bam")

## ----plot_code,fig.height=10,fig.width=7,prompt=FALSE,eval = FALSE-------
#  region <- GRanges("chr12:110938000-110940000")
#  RepViz::RepViz(region = region,
#                 genome = "hg19",
#                 BAM = "BAM_input.csv",
#                 BED = "BED_input.csv",
#                 avgTrack = T,
#                 geneTrack = T,
#                 verbose = F)
#  
#  
#  

## ----plot_plot,fig.height=10,fig.width=7,prompt=FALSE, echo=FALSE--------
region <- GRanges("chr12:110938000-110940000")
RepViz::RepViz(region = region,
               genome = "hg19",
               BAM = "BAM_input.csv",
               BED = "BED_input.csv",
               avgTrack = T,
               geneTrack = T,
               verbose = F)




