\name{filter.genes}
\alias{filter.genes}
\title{Filter gene expression based on gene counts}
\usage{
  filter.genes(gene.counts, gene.data, gene.filters)
}
\arguments{
  \item{gene.counts}{a matrix of gene counts, preferably
  after the normalization procedure.}

  \item{gene.data}{an annotation data frame usually
  obtained with \code{\link{get.annotation}} containing the
  unique gene accession identifiers.}

  \item{gene.filters}{a named list with gene filters and
  their parameters. See the main help page of
  \code{\link{metaseqr}} for details.}
}
\value{
  a named list whose names are the gene filter names and
  its members are the filtered rownames of
  \code{gene.data}.
}
\description{
  This function performs the gene expression filtering
  based on gene read counts and a set of gene filter rules.
  For more details see the main help pages of
  \code{\link{metaseqr}}.
}
\examples{
\dontrun{
gene.counts <- data("mm9.gene.counts",package="metaseqr")
sample.list=list(e15.5=c("e15.5_1","e15.5_2"),P0.5=c("P0.5_1","P0.5_2"),
  P60=c("P60_1","P60_2"))
gene.counts <- normalize.edger(gene.counts,sample.list)
gene.data <- get.annotation("mm9","gene")
gene.filters <- get.defaults("gene.filter")
filter.results <- filter.genes(gene.counts,gene.data,sample.list,gene.filters)
}
}
\author{
  Panagiotis Moulos
}

