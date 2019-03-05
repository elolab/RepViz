## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----input BAM-----------------------------------------------------------
BAM_table <- read.table(system.file("extdata","BAM_input.csv", package = "RepViz"), sep = ";")
colnames(BAM_table) <- c("bam file","group","group name")
knitr::kable(BAM_table)

## ----input BED-----------------------------------------------------------
BED_table <- read.table(system.file("extdata","BED_input.csv", package = "RepViz"), sep = ";")
colnames(BED_table) <- c("bed file","Legend")
knitr::kable(BED_table)

## ----bam files, echo=FALSE-----------------------------------------------
list.files(system.file("extdata", package = "RepViz"),pattern = ".bam")

## ----GRanges, echo=FALSE,prompt=FALSE------------------------------------
suppressMessages(library(GenomicRanges))

## ----plot,fig.height=10,fig.width=12-------------------------------------
region <- GRanges("chr12:110938000-110940000")
backup <- getwd()
setwd(system.file("extdata", package = "RepViz"))
RepViz::RepViz(region = region,
               genome = "hg19",
               BAM = "BAM_input.csv",
               BED = "BED_input.csv",
               avgTrack = T,
               geneTrack = T)




## ----backup--------------------------------------------------------------
setwd(backup)

