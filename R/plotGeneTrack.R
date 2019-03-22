############################################################################################################################################
##
## plotRegiont.R -- plot a snapshot of a genomic region
## 30 november 2018
## Thomas Faux
## Medical Bioinformatics Centre
############################################################################################################################################

## geneRanges -- get the genes as GRanges

## splitColumnByOverlap -- helper function

## toGRList -- transform the dataframe in GRangesList

## findGenes -- find genes overlaping with the region for the given database

## plotGenomicTrack -- plot the genomic track

# get the genes as GRanges
#
# @param db object of the type BSgenome.Hsapiens.UCSC.hg19
# @param column character, name of the column you wish to retrive (default ENTREZ_ID)
# @return returns a GRanges of genes
############################################################################################################################################
##
## geneRanges -- get the genes as GRanges
##
## db     -- db object : object of the type BSgenome.Hsapiens.UCSC.hg19
##
## column -- character : name of the column you wish to retrive (default ENTREZ_ID)
##
## returns a GRanges of genes
############################################################################################################################################
geneRanges <- function(db, column="ENTREZID")
  {
    g <- genes(db, columns=column)
    col <- mcols(g)[[column]]
    genes <- granges(g)[rep(seq_along(g), IRanges::elementNROWS(col))]
    mcols(genes)[[column]] <- as.character(unlist(col))
    genes
  }


# helper function to split the genes found
#
# @param query GRanges query object
# @param subject GRanges subject
# @param column character for the name of columns you whish to retrive(default ENTREZ_ID)
# @return returns a GRangesList object
############################################################################################################################################
##
## splitColumnByOverlap -- helper function to split the genes found
##
## query   -- GRanges : query
##
## subject -- GRanges : subject
##
## column -- character : name of the column you wish to retrive (default ENTREZ_ID)
##
## returns a GRangesList object
############################################################################################################################################
splitColumnByOverlap <- function(query, subject, column="ENTREZID")
  {
    olaps <- GenomicRanges::findOverlaps(query, subject)
    f1 <- factor(S4Vectors::subjectHits(olaps),
                 levels=seq_len(S4Vectors::subjectLength(olaps)))
    IRanges::splitAsList(GenomicRanges::mcols(query)[[column]][S4Vectors::queryHits(olaps)], f1)
}


# transform the dataframe in GRangesList
#
# @param df list of dataframes containing the peaks
# @return returns a GRangesList
############################################################################################################################################
## not used ?
## toGRList -- transform the dataframe in GRangesList
##
## df    -- dataframe    : list of dataframes containing the peaks
##
## returns a GRangesList
############################################################################################################################################
toGRList <- function(df){
  grlist <- GenomicRanges::GRangesList()
  for(i in unique(df$external_gene_name)){
    temp <- GenomicRanges::makeGRangesFromDataFrame(df[which(df$external_gene_name == i),])
    grlist[[i]] <- temp
  }
  names(grlist) <- unique(df$external_gene_name)
  return(grlist)
}

# find genes overlaping with the region for the given database
#
# @param region GRanges object containing the region to plot
# @param genome Genome used by the db object "hg19","GRCh38" or "mm10"
# @param m biomaRt object
# @return returns GRanges with the Genes hits found in the region
############################################################################################################################################
##
## findGenes -- find genes overlaping with the region for the given database
##
## region -- GRanges : GRanges object containing the region to plot
##
## m -- biomaRt object
##
## returns genomic range of the genes present in the genomic region
############################################################################################################################################
findGenes <- function(region,m){
  filters <-c("chromosome_name","start","end")
  chr <- substr(as.character(seqnames(region)),4,nchar(as.character(seqnames(region))))
  values <-  list(chr,GenomicRanges::start(region),GenomicRanges::end(region))
  map <- biomaRt::getBM(mart = m,
        attributes = c("ensembl_gene_id", "external_gene_name","ensembl_exon_id","chromosome_name","exon_chrom_start","exon_chrom_end","strand"),
        filters = filters,
        values = values
  )

  map[which(map$strand == -1),"strand"] <- "-"
  map[which(map$strand == 1),"strand"] <- "+"
  map$chromosome_name <- as.character(GenomicRanges::seqnames(region))
  gr <- GenomicRanges::makeGRangesFromDataFrame(map,keep.extra.columns = TRUE)
  return(gr)
}

# find UTR overlaping with the region for the given database
#
# @param region GRanges object containing the region to plot
# @param m biomaRt object
# @return returns GRanges with the UTR hits found in the region
############################################################################################################################################
##
## findUTR -- find genes overlaping with the region for the given database
##
## region -- GRanges : GRanges object containing the region to plot
##
## genome -- character : Genome used by the db object "hg19","GRCh38" or "mm10"
##
## returns genomic rangelist of the UTR regions present in the genomic region per gene
############################################################################################################################################
findUTR5 <- function(region,m){
  f <- function(UTR){
    return_object <- GenomicRanges::GRanges()
    for(gene in unique(UTR$external_gene_name)){
      temp <- UTR[which(UTR$external_gene_name == gene),]
      temp <- GenomicRanges::makeGRangesFromDataFrame(temp,keep.extra.columns = TRUE)
      temp <- temp[queryHits(GenomicRanges::findOverlaps(temp,region,ignore.strand=TRUE))]
      if(length(temp) > 0){
        temp <- temp[which(temp$transcript_length == max(temp$transcript_length))]
      }
      return_object <- c(return_object,temp)
    }
    return(return_object)
  }
  filters <-c("chromosome_name","start","end")
  chr <- substr(as.character(seqnames(region)),4,nchar(as.character(seqnames(region))))
  values <-  list(chr,GenomicRanges::start(region),GenomicRanges::end(region))
  map <- biomaRt::getBM(mart = m,
                        attributes = c("external_gene_name","chromosome_name","strand","5_utr_start","5_utr_end","transcript_length"),
                        filters = filters,
                        values = values
  )

  map[which(map$strand == -1),"strand"] <- "-"
  map[which(map$strand == 1),"strand"] <- "+"

  UTR5 <- map[!is.na(map$`5_utr_start`),]
  if(dim(UTR5)[1] > 0){
    UTR5$chromosome_name <- as.character(GenomicRanges::seqnames(region))
    UTR5 <- f(UTR5)

  }else{
    UTR5 <- GenomicRanges::GRanges()
  }
  return(UTR5)
}


findUTR3 <- function(region,m){
  f <- function(UTR){
    return_object <- GenomicRanges::GRanges()
    for(gene in unique(UTR$external_gene_name)){
      temp <- UTR[which(UTR$external_gene_name == gene),]
      temp <- GenomicRanges::makeGRangesFromDataFrame(temp,keep.extra.columns = TRUE)
      temp <- temp[queryHits(GenomicRanges::findOverlaps(temp,region,ignore.strand=TRUE))]
      if(length(temp) > 0){
        temp <- temp[which(temp$transcript_length == max(temp$transcript_length))]
      }
      return_object <- c(return_object,temp)
    }
    return(return_object)
  }
  filters <-c("chromosome_name","start","end")
  chr <- substr(as.character(seqnames(region)),4,nchar(as.character(seqnames(region))))
  values <-  list(chr,GenomicRanges::start(region),GenomicRanges::end(region))
  map <- biomaRt::getBM(mart = m,
                        attributes = c("external_gene_name","chromosome_name","strand","3_utr_start","3_utr_end","transcript_length"),
                        filters = filters,
                        values = values
  )

  map[which(map$strand == -1),"strand"] <- "-"
  map[which(map$strand == 1),"strand"] <- "+"

  UTR3 <- map[!is.na(map$`3_utr_start`),]
  if(dim(UTR3)[1] > 0){
    UTR3$chromosome_name <- as.character(GenomicRanges::seqnames(region))
    UTR3 <- f(UTR3)
  }else{
    UTR3 <- GenomicRanges::GRanges()
  }

  return(UTR3)
}
# get The biomaRt object
#
# @param region GRanges object containing the region to plot
# @param genome Genome used by the db object "hg19","GRCh38" or "mm10"
# @return returns GRanges with the Genes hits found in the region
############################################################################################################################################
##
## getBiomaRt -- find genes overlaping with the region for the given database
##
## region -- GRanges : GRanges object containing the region to plot
##
## genome -- character : Genome used by the db object "hg19","GRCh38" or "mm10"
##
## returns a biomaRt object
############################################################################################################################################
getBiomaRt <- function(region,genome=c("hg19","GRCh38","mm10")){

  if(genome == "hg19"){
    m <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  }
  if(genome == "GRCh38"){
    m <-  biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch38.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  }
  if(genome == "mm10"){
    m <-  biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  }

  return(m)
}


#  plot the rectangles of the genomic track
#
# @param rg GRanges list of exons in the gene
# @param gr GRanges object containing the region to plot
# @param y numeric index for color purpose
# @param size value to substract and add around the y value to create a rectangle
############################################################################################################################################
##
## plotRectangles -- plot the rectangles of the genomic track
##
## gr    -- GRanges    : list of exons in the gene
##
## region -- GRanges : GRanges object containing the region to plot
##
## y  -- numeric : index for color purpose
##
## size  -- value to substract and add around the y value to create a rectangle
############################################################################################################################################
plotRectangles <- function(rg,gr,y,size){
  graphics::rect(GenomicRanges::end(gr),rep(y,length(gr))-size,GenomicRanges::start(gr),rep(y,length(gr))+size,border=NA, col = "grey")#y+5)

}

#  plot the arrows between the rectangles in the genomic track
#
# @param x0 list of x start points
# @param y0 list of y start points
# @param x1 list of x end points
# @param y1 list of y end points
# @param n_arr number of arrows needed
# @param ... graphics::arrows basic parameters
############################################################################################################################################
##
## arrowLines -- plot the arrows between the rectangles in the genomic track
##
## x  -- vector : list of x points
##
## y  -- vector : list of y points
##
## N  -- vector : number of break points needed
##
############################################################################################################################################
arrowLine <- function(x0,y0, x1,y1,n_arr,...) {
  x<- seq(x0,x1, length=n_arr+1)
  y<-seq(y0,y1, length=n_arr+1)
  graphics::arrows(x[-length(x)],y[-length(y)], x[-1],y[-1],...)
}

#  return the regions contained between the ranges of a GRanges
#
# @param gr GRanges object
# @return GRange object
############################################################################################################################################
##
## inverseGRanges -- return the regions contained between the ranges of a GRanges
##
## gr    -- GRanges    : list of exons in the gene
##
## returns GRange
############################################################################################################################################
inverseGRanges <- function(gr){
  temp <- as.data.frame(gr)
  temp <- cbind(as.character(temp[1:dim(temp)[1]-1,1]),
                temp[1:(dim(temp)[1]-1),3],
                temp[2:dim(temp)[1],2],
                as.character(temp[1:dim(temp)[1]-1,5]))
  temp <- as.data.frame(temp)
  colnames(temp) <- c("seqnames","start","end","strand")
  temp <- GRanges(temp)
  return(temp)
}

#  wrapper function to plot the rectangles and arrows of the genomic track
#
# @param gr GRanges object containing the region to plot
# @param region GRanges object containing the region to plot
# @param y Numeric value used for the y position of the boxes
############################################################################################################################################
##
## plotGRanges-- wrapper function to plot the rectangles and arrows of the genomic track
##
## gr    -- GRanges    : list of exons in the gene
##
## region -- GRanges : GRanges object containing the region to plot
##
## returns a plot
############################################################################################################################################
plotGRanges <- function(gr,region,y,UTR3,UTR5){
  if(length(UTR5) > 0){
    plotRectangles(region,UTR5,y,size=0.25)
    if(as.character(unique(strand(gr))) == "-"){
      GenomicRanges::end(gr[queryHits(GenomicRanges::findOverlaps(gr,UTR5))]) <- GenomicRanges::start(UTR5[subjectHits(GenomicRanges::findOverlaps(gr,UTR5))])
    }
    if(as.character(unique(strand(gr))) == "+"){
      GenomicRanges::start(gr[queryHits(GenomicRanges::findOverlaps(gr,UTR5))]) <- GenomicRanges::end(UTR5[subjectHits(GenomicRanges::findOverlaps(gr,UTR5))])
    }
  }
  if(length(UTR3) > 0){
    plotRectangles(region,UTR3,y,size=0.25)
    if(as.character(unique(strand(gr))) == "+"){
      GenomicRanges::end(gr[queryHits(GenomicRanges::findOverlaps(gr,UTR3))]) <- GenomicRanges::start(UTR3[subjectHits(GenomicRanges::findOverlaps(gr,UTR3))])
    }
    if(as.character(unique(strand(gr))) == "-"){
      GenomicRanges::start(gr[queryHits(GenomicRanges::findOverlaps(gr,UTR3))]) <- GenomicRanges::end(UTR3[subjectHits(GenomicRanges::findOverlaps(gr,UTR3))])
    }
  }
  plotRectangles(region,gr,y,size=0.5)

  if(length(gr) > 1){
    gr <- inverseGRanges(gr)
    for(j in 1:length(gr)){
      N <- 0
      value <- (GenomicRanges::end(gr[j])-GenomicRanges::start(gr[j]))/(GenomicRanges::end(region)-GenomicRanges::start(region))


      if(value < 0.25){
        N=1
      }
      if(value > 0.25){
        N=2
      }
      if(value > 0.75){
        N=3
      }
      if(value > 1){
        N=6
      }
      if(value > 2){
        N=12
      }
      if(as.character(strand(gr[j])) == "-"){
        arrowLine(GenomicRanges::end(gr[j]),y,GenomicRanges::start(gr[j]),y,n_arr=N,length=0.1)
      }
      if(as.character(strand(gr[j])) == "+"){
        arrowLine(GenomicRanges::start(gr[j]),y,GenomicRanges::end(gr[j]),y,n_arr=N,length=0.1)
      }

    }
  }

}

# plot the genomic track
#
# @param gr GRanges object containing the region to plot
# @param region GRanges object containing the region to plot
############################################################################################################################################
##
## plotGenomicTrack -- plot the genomic track
##
## gr    -- GRanges    : list of exons in the gene
##
## region -- GRanges : GRanges object containing the region to plot
##
## returns a plot
############################################################################################################################################
plotGenomicTrack <- function(gr,UTR3,UTR5,region){

  if(dim(as.data.frame(gr))[1]>0){
    index <- 0
    if(length(unique(gr$external_gene_name))==1){
      graphics::plot(x=1, ylim=c(0,2),xlim=c(GenomicRanges::start(region),GenomicRanges::end(region)),xlab=GenomicRanges::seqnames(region),yaxt="n",ylab="",xaxt="n")
      graphics::axis(2, at=1, unique(gr$external_gene_name),srt=45)
      plotGRanges(GenomicRanges::reduce(gr),region,1,UTR3,UTR5)
    }
    if(length(unique(gr$external_gene_name))>1){
      graphics::plot(x=1, ylim=c(0,length(unique(gr$external_gene_name))),xlim=c(start(region),end(region)),xlab=seqnames(region),yaxt="n",ylab="")
      tmin <- 1
      tmax <- length(IRanges::unique(gr$external_gene_name))
      tlab <- seq(tmin, tmax)-0.5
      lab <- IRanges::unique(gr$external_gene_name)
      graphics::axis(2, at=tlab, labels = FALSE)
      graphics::text(x=par()$usr[1]-0.01*(par()$usr[2]-par()$usr[1]),
                     y=tlab,
                     labels=lab,
                     srt=45,
                     xpd=NA,
                     adj=1)
      for(i in unique(gr$external_gene_name)){
        index <- index+1
        gr2 <- gr[(GenomicRanges::elementMetadata(gr)[,"external_gene_name"] %in% i)]
        UTR3_ind <- UTR3[(GenomicRanges::elementMetadata(UTR3)[,"external_gene_name"] %in% i)]
        UTR5_ind <- UTR5[(GenomicRanges::elementMetadata(UTR5)[,"external_gene_name"] %in% i)]
        plotGRanges(GenomicRanges::reduce(gr2),region,index-0.5,UTR3=UTR3_ind,UTR5=UTR5_ind)
      }

    }
  }
  else{
    print(paste0("There is no genes in the given region :",region))
  }

}


