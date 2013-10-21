#' metaseqr pipeline starting from raw BAM or BED files
#'
#' This function is a wrapper for the \code{metaseqr} pipeline, starting from BAM or BED files instead of a read counts file.
#'
#' @param files.list a named list whose members are named vectors. The names of the list correspond to condition names (see the
#' \code{sample.list} argument in the main \code{\link{metaseqr}} function). The names of the vectors are the sample names and the
#' vector elements are full paths to BAM/BED files.
#' @param file.type the type of raw input files. It can be \code{"bed"} for BED files or \code{"bam"} for BAM files. See the same
#' argument in the main \code{\link{metaseqr}} function for the case of auto-guessing.
#' @param annotation see the \code{annotation} argument in the main \code{\link{metaseqr}} function. The \code{"annotation"} parameter
#' here is the result of the same parameter in the main function. See also \code{\link{get.annotation}} and \code{\link{read.annotation}}.
#' @param has.all.fields a logical variable indicating if all annotation fields used by \code{metaseqr} are available (that is apart
#' from the main chromosome, start, end, unique id and strand columns, if also present are the gene name and biotype columns). The
#' default is \code{FALSE}.
#' @return A data frame with counts for each sample, ready to be passed to the main \code{\link{metaseqr}} pipeline.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' # Not yet implenented
bam2metaseqr <- function(files.list,file.type,annotation,has.all.fields=FALSE) {
	if (!require(GenomicRanges))
		stop("The Bioconductor package GenomicRanges is required to proceed!")
	if (file.type=="bed" !require(rtracklayer))
		stop("The Bioconductor package rtracklayer is required to process BED files!")
	if (file.type=="bam" !require(Rsamtools))
		stop("The Bioconductor package Rsamtools is required to process BAM files!")

	# Convert annotation to GRanges
	disp("Converting annotation to GenomicRanges object...")
	if (packageVersion("GenomicRanges")<1.14) { # Classic way
		if (has.all.fields)
			annotation.gr <- GRanges(
				seqnames=Rle(annotation[,1]),
				ranges=IRanges(start=annotation[,2],end=annotation[,3]),
				strand=Rle(annotation[,6]),
				name=as.character(annotation[,4]),
				symbol=as.character(gene.data[,7]),
				biotype=as.character(gene.data[,8])
			)
		else
			annotation.gr <- GRanges(
				seqnames=Rle(annotation[,1]),
				ranges=IRanges(start=annotation[,2],end=annotation[,3]),
				strand=Rle(annotation[,6]),
				name=as.character(annotation[,4])
			)
	}
	else # Use native method in newer versions of GenomicRanges
		annotation.gr <- makeGRangesFromDataFrame(
			df=annotation,
			keep.extra.columns=TRUE,
			seqnames.field="chromosome"
		)

	sample.names <- as.character(sapply(files.list,names))
	sample.files <- unlist(files.list,use.names=FALSE)
	names(sample.files) <- sample.names
	counts <- matrix(0,nrow=length(annotation.gr),ncol=length(sample.names))
	rownames(counts) <- as.character(annotation[,4])
	colnames(counts) <- sample.names
	
	if (file.type=="bed") {
		for (n in sample.names) {
		disp("Reading bed file ",basename(sample.files[n])," for sample with name ",n,". This might take some time...")
		bed <- import.bed(sample.files[n],trackLine=FALSE,asRangedData=FALSE)
		disp("  Counting reads overlapping with given annotation...")
		counts[,n] <- countOverlaps(annotation.gr,bed)
	}
	else if (file.type=="bam") {
		# Requires more reading...
		# 1. They must be indexed
		# 2. Use facilities in GenomicRanges to read the required BAM fields
		# 3. What about paired-end? Probably collapse to single-end...
	}

	return(counts)
}

#' Creates sample list and BAM/BED file list from file
#'
#' Create the main sample list and determine the BAM/BED files for each sample from an external file.
#'
#' @param input a tab-delimited file structured as follows: the first line of the external tab delimited file should contain column 
#' names (names are not important). The first column MUST contain UNIQUE sample names. The second column MUST contain the raw BAM/BED
#' files WITH their full path. Alternatively, the \code{path} argument should be provided (see below). The third column MUST contain 
#' the biological condition where each of the samples in the first column should belong to.
#' @param path an optional path where all the BED/BAM files are placed, to be prepended to the BAM/BED file names in the targets file.
#' @return A named list with three members. The first member is a named list whose names are the conditions of the experiments and its
#' members are the samples belonging to each condition. The second member is like the first, but this time the members are named vectors
#' whose names are the sample names and the vector elements are full path to BAM/BED files. The third member is the guessed type of
#' the input files (BAM or BED). It will be used if not given in the main \code{link\{bam2metaseqr}} function.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' targets <- data.frame(sample=c("C1","C2","T1","T2"),
#'   filename=c("C1_raw.bam","C2_raw.bam","T1_raw.bam","T2_raw.bam"),
#'   condition=c("Control","Control","Treatment","Treatment"))
#' path <- "/home/chakotay/bam"
#' write.table(targets,file="targets.txt",sep="\t",row.names=F,quote="")
#' the.list <- read.targets("targets.txt",path=path)
#' sample.list <- the.list$samples
#' bamfile.list <- the.list$files
#'}
read.targets <- function(input,path=NULL) {
	if (missing(input) || !file.exists(input))
		stop("The targets file should be a valid existing text file!")
	tab <- read.delim(input)
	samples <- as.character(tab[,1])
	conditions <- unique(as.character(tab[,3]))
	rawfiles <- as.character(tab[,2])
	if (!is.null(path)) {
		tmp <- dirname(rawfiles) # Test if there is already a path
		if (any(tmp=="."))
			rawfiles <- file.path(path,basename(rawfiles))
	}
	if (length(samples) != length(unique(samples)))
		stop("Sample names must be unique for each sample!")
	if (length(rawfiles) != length(unique(rawfiles)))
		stop("File names must be unique for each sample!")
	sample.list <- vector("list",length(conditions))
	names(sample.list) <- conditions
	for (n in conditions)
		sample.list[[n]] <- samples[which(as.character(tab[,3]))==n]
	file.list <- vector("list",length(conditions))
	names(sample.list) <- conditions
	for (n in conditions) {
		file.list[[n]] <- rawfiles[which(as.character(tab[,3]))==n]
		names(file.list[[n]]) <- samples[which(as.character(tab[,3]))==n]
	}
	# Guess file type based on only one of them
	tmp <- file.list[[1]][1]
	if (grep("\\.bam$",tmp,ignore.case=TRUE,perl=TRUE)>0)
		type <- "bam"
	else if (grep("\\.bed$",tmp,ignore.case=TRUE,perl=TRUE)>0)
		type <- "bed"
	else
		type <- NULL
	return(list(samples=sample.list,files=file.list,type=type))
}
