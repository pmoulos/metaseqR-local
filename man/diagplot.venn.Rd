\name{diagplot.venn}
\alias{diagplot.venn}
\title{Venn diagrams when performing meta-analysis}
\usage{
  diagplot.venn(pmat, pcut = 0.05,
    nam = as.character(round(1000 * runif(1))),
    output = "x11", path = NULL, alt.names = NULL, ...)
}
\arguments{
  \item{pmat}{a matrix with p-values corresponding to the
  application of each statistical algorithm. The p-value
  matrix must have the colnames attribute and the colnames
  should correspond to the name of the algorithm used to
  fill the specific column (e.g. if
  \code{"statistics"=c("deseq","edger","nbpseq")} then
  \code{colnames(pmat) <-}
  \code{c("deseq","edger","nbpseq")}.}

  \item{pcut}{a p-value cutoff for statistical
  significance. Defaults to \code{0.05}.}

  \item{nam}{a name to be appended to the output graphics
  file (if \code{"output"} is not \code{"x11"}).}

  \item{output}{one or more R plotting device to direct the
  plot result to. Supported mechanisms: \code{"x11"}
  (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
  \code{"pdf"} or \code{"ps"}.}

  \item{path}{the path to create output files.}

  \item{alt.names}{an optional named vector of names, e.g.
  HUGO gene symbols, alternative or complementary to the
  unique gene names which are the rownames of \code{pmat}.
  The names of the vector must be the rownames of
  \code{pmat}.}

  \item{...}{further arguments to be passed to plot
  devices, such as parameter from \code{\link{par}}.}
}
\value{
  The filenames of the plots produced in a named list with
  names the \code{which.plot} argument. If
  output=\code{"x11"}, no output filenames are produced.
}
\description{
  This function uses the R package VennDiagram and plots an
  up to 5-way Venn diagram depicting the common and
  specific to each statistical algorithm genes, for each
  contrast. Mostly for internal use because of its main
  argument which is difficult to construct, but can be used
  independently if the user grasps the logic.
}
\examples{
\dontrun{
# Not yet available...
}
}
\author{
  Panagiotis Moulos
}

