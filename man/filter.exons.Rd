\name{filter.exons}
\alias{filter.exons}
\title{Filter gene expression based on exon counts}
\usage{
  filter.exons(the.counts, gene.data, sample.list,
    exon.filters, restrict.cores = 0.8)
}
\arguments{
  \item{the.counts}{a named list created with the
  \code{\link{construct.gene.model}} function. See its help
  page for details.}

  \item{gene.data}{an annotation data frame usually
  obtained with \code{\link{get.annotation}} containing the
  unique gene accession identifiers.}

  \item{sample.list}{the list containing condition names
  and the samples under each condition.}

  \item{exon.filters}{a named list with exon filters and
  their parameters. See the main help page of
  \code{\link{metaseqr}} for details.}

  \item{restrict.cores}{in case of parallel execution of
  several subfunctions, the fraction of the available cores
  to use. In some cases if all available cores are used
  (\code{restrict.cores=1} and the system does not have
  sufficient RAM, the running machine might significantly
  slow down.}
}
\value{
  a named list with two members. The first member 
  (\code{result} is a named list whose names are the 
  exon filter names and its members are the filtered 
  rownames of \code{gene.data}. The second member is a 
  matrix of binary flags (0 for non-filtered, 1 for 
  filtered) for each gene. The rownames of the flag 
  matrix correspond to gene ids.
}
\description{
  This function performs the gene expression filtering
  based on exon read counts and a set of exon filter rules.
  For more details see the main help pages of
  \code{\link{metaseqr}}.
}
\examples{
\dontrun{
exon.counts <- data("hg18.exon.counts",package="metaseqr")
gene.data <- get.annotation("hg18","gene")
sample.list <- list(CON=c("CON_BR1","CON_BR2"),DOX=c("DOX_BR1","DOX_BR2"))
exon.filters <- get.defaults("exon.filter")
the.counts <- construct.gene.model(exon.counts,sample.list,gene.data)
filter.results <- filter.exons(the.counts,gene.data,sample.list,exon.filters)
}
}
\author{
  Panagiotis Moulos
}

