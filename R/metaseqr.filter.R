#' Filter gene expression based on exon counts
#'
#' This function performs the gene expression filtering based on exon read counts and a set of exon filter rules. For more details
#' see the main help pages of \code{\link{metaseqr}}.
#'
#' @param the.counts a named list created with the \code{\link{construct.gene.model}} function. See its help page for details.
#' @param gene.data an annotation data frame usually obtained with \code{\link{get.annotation}} containing the unique gene accession
#' identifiers.
#' @param sample.list the list containing condition names and the samples under each condition.
#' @param exon.filters a named list with exon filters and their parameters. See the main help page of \code{\link{metaseqr}} for details.
#' @param restrict.cores in case of parallel execution of several subfunctions, the fraction of the available cores to use. In some 
#' cases if all available cores are used (\code{restrict.cores=1} and the system does not have sufficient RAM, the running machine 
#' might significantly slow down.
#' @return a named list whose names are the exon filter names and its members are the filtered row indices of gene.data.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' exon.counts <- data("hg18.exon.counts",package="metaseqr")
#' gene.data <- get.annotation("hg18","gene")
#' sample.list <- list(CON=c("CON_BR1","CON_BR2"),DOX=c("DOX_BR1","DOX_BR2"))
#' exon.filters <- get.defaults("exon.filter")
#' the.counts <- construct.gene.model(exon.counts,sample.list,gene.data)
#' filter.results <- filter.exons(the.counts,gene.data,sample.list,exon.filters)
#'}
filter.exons <- function(the.counts,gene.data,sample.list,exon.filters,restrict.cores=0.8) {
	multic <- check.parallel(restrict.cores)
	exon.filter.result <- vector("list",length(exon.filters))
	names(exon.filter.result) <- names(exon.filters)
	the.genes <- gene.data$gene_id
	if (!is.null(exon.filters))
	{
		for (xf in names(exon.filters))
		{
			disp("Applying exon filter ",xf,"...")
			switch(xf,
				min.active.exons = {
					pass <- vector("list",length(unlist(sample.list)))
					names(pass) <- names(the.counts)
					for (n in names(pass))
					{
						disp("  Checking read presence in exons for ",n,"...")
						pass[[n]] <- the.genes
						#names(pass[[n]]) <- names(the.gene.counts[[n]]) <- the.genes
						names(pass[[n]]) <- the.genes
						pass[[n]] <- wapply(multic,the.counts[[n]],function(x,f) {
							if (length(x$count) <= f$exons.per.gene)
								if (length(which(x$count!=0)) >= f$min.exons)
									return(TRUE)
								else
									return(FALSE)
							else
								if (length(which(x$count!=0)) >= ceiling(length(x)/f$frac))
									return(TRUE)
								else
									return(FALSE)
							},exon.filters$min.active.exons)
						pass[[n]] <- do.call("c",pass[[n]])
					}
					pass.matrix <- do.call("cbind",pass)
					exon.filter.result[[xf]] <- which(apply(pass.matrix,1,function(x) return(any(x))))
				}
				# More to come...
				# TODO: Write more rules based in exons
			)
		}
	}
	return(exon.filter.result)
}

#' Filter gene expression based on gene counts
#'
#' This function performs the gene expression filtering based on gene read counts and a set of gene filter rules. For more details
#' see the main help pages of \code{\link{metaseqr}}.
#'
#' @param gene.counts a matrix of gene counts, preferably after the normalization procedure.
#' @param gene.data an annotation data frame usually obtained with \code{\link{get.annotation}} containing the unique gene accession
#' identifiers.
#' @param gene.filters a named list with gene filters and their parameters. See the main help page of \code{\link{metaseqr}} for details.
#' @return a named list whose names are the gene filter names and its members are the filtered row indices of gene.data.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' gene.counts <- data("mm9.gene.counts",package="metaseqr")
#' sample.list=list(e15.5=c("e15.5_1","e15.5_2"),P0.5=c("P0.5_1","P0.5_2"),
#'   P60=c("P60_1","P60_2"))
#' gene.counts <- normalize.edger(gene.counts,sample.list)
#' gene.data <- get.annotation("mm9","gene")
#' gene.filters <- get.defaults("gene.filter")
#' filter.results <- filter.genes(gene.counts,gene.data,sample.list,gene.filters)
#'}
filter.genes <- function(gene.counts,gene.data,gene.filters)
{
	gene.filter.result <- vector("list",length(gene.filters))
	names(gene.filter.result) <- names(gene.filters)
	for (gf in names(gene.filters)) {
		disp("Applying gene filter ",gf,"...")
		switch(gf,
			length = { # This is real gene length independently of exons
				gene.filter.result$length <- which(gene.data$end - gene.data$start < gene.filters$length$length)
			},
			avg.reads = {
				avg.mat <- sweep(gene.counts,1,attr(gene.data,"gene.length")/gene.filters$avg.reads$average.per.bp,"/")
				q.t <- max(apply(avg.mat,2,quantile,gene.filters$avg.reads$quantile))
				gene.filter.result$avg.reads <- which(apply(avg.mat,1,filter.low,q.t))
			},
			expression = {
				if (gene.filters$expression$median)
					the.dead.median <- which(apply(gene.counts,1,filter.low,median(gene.counts)))
				else
					the.dead.median <- NULL
				if (gene.filters$expression$mean)
					the.dead.mean <- which(apply(gene.counts,1,filter.low,mean(gene.counts)))
				else
					the.dead.mean <- NULL
				if (!is.na(gene.filters$expression$quantile))
					the.dead.quantile <- which(apply(gene.counts,1,filter.low,quantile(gene.counts,gene.filters$expression$quantile)))
				else
					the.dead.quantile <- NULL
				if (!is.na(gene.filters$expression$known)) {
					bio.cut <- match(gene.filters$expression$known,gene.data$gene_name) # Think about the case of embedded
					bio.cut <- bio.cut[-which(is.na(bio.cut))]
					bio.cut.counts <- as.vector(gene.counts[bio.cut,])
					the.bio.cut <- quantile(bio.cut.counts,0.9)
					the.dead.known <- which(apply(gene.counts,1,filter.low,the.bio.cut))
				}
				else
					the.dead.known <- NULL
				if (!is.na(gene.filters$expression$custom)) {
					# For future use
					the.dead.custom <- NULL
				}
				else
					the.dead.custom <- NULL
				# Derive one common expression filter
				the.dead <- list(the.dead.median,the.dead.mean,the.dead.quantile,the.dead.known,the.dead.custom)
				gene.filter.result$expression <- Reduce("union",the.dead)
			},
			biotype = {
				if (!is.null(gene.filters$biotype)) {
					filter.out <- names(which(unlist(gene.filters$biotype)))
					if (length(grep("three_prime_overlapping_ncrna",filter.out))>0) # Necessary hack because of R naming system
						filter.out <- sub("three_prime_overlapping_ncrna","3prime_overlapping_ncrna",filter.out)
					filter.ind <- vector("list",length(filter.out))
					names(filter.ind) <- filter.out
					for (bt in filter.out)
						filter.ind[[bt]] <- which(gene.data$biotype==bt)
					gene.filter.result$biotype <- Reduce("union",filter.ind)
				}
				else
					gene.filter.result$biotype <- NULL
			}
		)
	}
	return(gene.filter.result)
}

