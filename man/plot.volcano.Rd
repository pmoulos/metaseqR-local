\name{plot.volcano}
\alias{plot.volcano}
\title{(Interactive) volcano plots of differentially expressed genes}
\usage{
  plot.volcano(f, p, con = NULL, fcut = 1, pcut = 0.05,
    alt.names = NULL, output = "x11", path = NULL, ...)
}
\arguments{
  \item{f}{the fold changes which are to be plotted on the
  x-axis.}

  \item{p}{the p-values whose -log10 transformation is
  going to be plotted on the y-axis.}

  \item{con}{an optional string depicting a name (e.g. the
  contrast name) to appear in the title of the volcano
  plot.}

  \item{fcut}{a fold change cutoff so as to draw two
  vertical lines indicating the cutoff threshold for
  biological significance.}

  \item{pcut}{a p-value cutoff so as to draw a horizontal
  line indicating the cutoff threshold for statistical
  significance.}

  \item{alt.names}{an optional vector of names, e.g. HUGO
  gene symbols, alternative or complementary to the unique
  names of f or p (one of them must be named!). It is used
  only in JSON output.}

  \item{output}{one or more R plotting device to direct the
  plot result to. Supported mechanisms: "x11" (default),
  "png", "jpg", "bmp", "pdf", "ps" or "json". The latter is
  currently available for the creation of interactive
  volcano plots only when reporting the output, through the
  highcharts javascript library.}

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
  This function plots a volcano plot or returns a JSON
  string which is used to render aninteractive in case of
  HTML reporting.
}
\examples{
\dontrun{
require(DESeq)
data.matrix <- counts(makeExampleCountDataSet())
sample.list <- list(A=c("A1","A2"),B=c("B1","B2",B3"))
contrast <- "A_vs_B"
M <- norm.edger(data.matrix,sample.list)
p <- stat.edger(M,sample.list,contrast)
ma <- apply(M[,sample.list$A],1,mean)
mb <- apply(M[,sample.list$B],1,mean)
f <- log2(mb/ma)
plot.volcano(f,p,con=contrast)
j <- plot.volcano(f,p,con=contrast,output="json")
}
}
\author{
  Panagiotis Moulos
}

