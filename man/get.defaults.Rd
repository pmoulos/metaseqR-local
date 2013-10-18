\name{get.defaults}
\alias{get.defaults}
\title{Default parameters for several metaseqr functions}
\usage{
  get.defaults(what, method = NULL)
}
\arguments{
  \item{what}{a keyword determining the procedure for which
  to fetch the default settings according to method
  parameter. It can be one of "normalization",
  "statistics", "gene.filter", "exon.filter" or
  "biotype.filter".}

  \item{method}{the supported algorithm included in
  metaseqr for which to fetch the default settings. When
  what is "normalization", method is one of "edaseq",
  "deseq", "edger" or "noiseq". When what is "statistics",
  method is one of "deseq", "edger", "noiseq", "bayseq",
  "limma" or "ebseq". When method is "biotype.filter", what
  is the input organism (see the main
  \code{\link{metaseqr}} help page for a list of supported
  organisms).}
}
\value{
  A list with default setting that can be used directly in
  the call of metaseqr.
}
\description{
  This function returns a list with the default settings
  for each filtering, statistical and normalization
  algorithm included in the metaseqr package. See the
  documentation of the main function and the documentation
  of each statistical and normalization method for details.
}
\examples{
\dontrun{
norm.args.edaseq <- get.defaults("normalization","edaseq")
stat.args.edger <- get.defaults("statistics","edger")
}
}
\author{
  Panagiotis Moulos
}

