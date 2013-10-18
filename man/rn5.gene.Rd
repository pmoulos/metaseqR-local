\docType{data}
\name{rn5.gene}
\alias{rn5.gene}
\title{Rat genome (genes) annotation data, build rn5}
\format{a \code{data.frame} with 8 columns, one row per gene.}
\source{
  Ensembl
}
\description{
  This data set contains predefined annotation for human
  genes, build rn5, if for any reason the user want to use
  fixed annotation for metaseqr, instead of downloading
  (e.g. no internet connection or the biomaRt package is
  not available or for speed). The format is a data frame
  with 8 columns: chromosome, start, end, gene_id,
  gc_content, strand, gene_name, biotype. This data frame
  does not contain information about isoforms. All
  annotation data built with the package can be updated
  with the \code{\link{annotations.update}} function.
}
\author{
  Panagiotis Moulos
}
\keyword{datasets}

