\name{diagplot.avg.ftd}
\alias{diagplot.avg.ftd}
\title{Create average False (or True) Discovery 
  curves}
\usage{
  diagplot.avg.ftd(ftdr.obj, output = "x11",
    path = NULL, draw = TRUE, ...)
}
\arguments{
  \item{ftdr.obj}{a list with outputs from 
  \code{\link{diagplot.ftd}}.}

  \item{output}{one or more R plotting 
  device to direct the plot result to.
  Supported mechanisms: \code{"x11"} (default), 
  \code{"png"}, \code{"jpg"}, \code{"bmp"}, 
  \code{"pdf"} or \code{"ps"}.}

  \item{path}{the path to create output files.}

  \item{draw}{boolean to determine whether 
  to plot the curves or just return the 
  calculated values (in cases where the user 
  wants the output for later averaging
  for example). Defaults to \code{TRUE} (make 
  plots).}

  \item{...}{further arguments to be passed to 
  plot devices, such as parameter from 
  \code{\link{par}}.}
}
\value{
  A named list with two members: the first member 
  (\code{avg.ftdr}) contains a list with the 
  means and the standard deviations of the averaged
  \code{ftdr.obj} and are used to create the plot. 
  The second member (\code{path}) contains the 
  path to the created figure graphic.
}
\description{
  This function creates false (or true) discovery 
  curves using a list containing several outputs 
  from \code{\link{diagplot.ftd}}.
}
\examples{
\dontrun{
# Not yet available
}
}
\author{
  Panagiotis Moulos
}

