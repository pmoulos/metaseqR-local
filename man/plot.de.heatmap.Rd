\name{plot.de.heatmap}
\alias{plot.de.heatmap}
\title{Diagnostic heatmap of differentially expressed genes}
\usage{
  plot.de.heatmap(x, con = NULL, output = "x11",
    path = NULL, ...)
}
\arguments{
  \item{x}{the data matrix to create a heatmap for.}

  \item{con}{an optional string depicting a name (e.g. the
  contrast name) to appear in the title of the volcano
  plot.}

  \item{output}{one or more R plotting device to direct the
  plot result to. Supported mechanisms: "x11" (default),
  "png", "jpg", "bmp", "pdf", "ps".}

  \item{path}{the path to create output files.}

  \item{...}{further arguments to be passed to plot
  devices, such as parameter from \code{\link{par}}.}
}
\value{
  The filenames of the plots produced in a named list with
  names the which.plot argument. If output="x11", no output
  filenames are produced.
}
\description{
  This function plots a heatmap of the differentially
  expressed genes produced by the metaseqr workflow, useful
  for quality control, e.g. whether samples belonging to
  the same group cluster together.
}
\examples{
\dontrun{
require(DESeq)
data.matrix <- counts(makeExampleCountDataSet())
sample.list <- list(A=c("A1","A2"),B=c("B1","B2",B3"))
contrast <- "A_vs_B"
M <- norm.edger(data.matrix,sample.list)
p <- stat.edger(M,sample.list,contrast)
plot.de.heatmap(data.matrix[p[[1]]<0.05])
}
}
\author{
  Panagiotis Moulos
}

