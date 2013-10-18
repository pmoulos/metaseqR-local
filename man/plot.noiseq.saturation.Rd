\name{plot.noiseq.saturation}
\alias{plot.noiseq.saturation}
\title{Simpler implementation of saturation plots inspired from NOISeq package}
\usage{
  plot.noiseq.saturation(x, o, tb, path = NULL)
}
\arguments{
  \item{x}{the count data matrix.}

  \item{o}{one or more R plotting device to direct the plot
  result to. Supported mechanisms: "x11" (default), "png",
  "jpg", "bmp", "pdf" or "ps".}

  \item{tb}{the vector of biotypes, one for each row of x.}

  \item{path}{the path to create output files.}
}
\value{
  The filenames of the plots produced in a named list with
  names the which.plot argument. If output="x11", no output
  filenames are produced.
}
\description{
  Helper function for \code{\link{plot.noiseq}} to plot
  feature detection saturation as presented in the NOISeq
  package vignette. It has two main outputs: a set of
  figures, one for each input sample depicting the
  saturation for each biotype and one single multiplot
  which depicts the saturation of all samples for each
  biotype. It expands the saturation plots of NOISeq by
  allowing more samples to be examined in a simpler way.
  Don't use this function directly. Use either
  \code{\link{plot.metaseqr}} or \code{\link{plot.noiseq}}.
}
\examples{
\dontrun{
require(DESeq)
data.matrix <- counts(makeExampleCountDataSet())
biotype=c(rep("protein_coding",nrow(data.matrix)/2),rep("ncRNA",nrow(data.matrix)/2))
plot.noiseq.saturation(data.matrix,"x11",biotype)
}
}
\author{
  Panagiotis Moulos
}

