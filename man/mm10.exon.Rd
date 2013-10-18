\docType{data}
\name{mm10.exon}
\alias{mm10.exon}
\title{Mouse genome (exons) annotation data, build mm10}
\format{a \code{data.frame} with 8 columns, one row per exon.}
\source{
  Ensembl
}
\description{
  This data set contains predefined annotation for human
  exons, build mm10, if for any reason the user want to use
  fixed annotation for metaseqr, instead of downloading
  (e.g. no internet connection or the biomaRt package is
  not available or for speed). The format is a data frame
  with 8 columns: chromosome, start, end, exon_id, exon_id,
  strand, exon_name, biotype. This data frame does not
  contain information about isoforms. All annotation data
  built with the package can be updated with the
  \code{\link{annotations.update}} function.
}
\author{
  Panagiotis Moulos
}
\keyword{datasets}

