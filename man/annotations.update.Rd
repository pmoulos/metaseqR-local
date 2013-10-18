\name{annotations.update}
\alias{annotations.update}
\title{Fixed annotation updater}
\usage{
  annotations.update()
}
\value{
  This function does not return anything. It updates the
  fixed annotation files instead.
}
\description{
  A function to update the fixed annotations contained to
  avoid downloading every time if it's not embedded. It has
  no parameters.
}
\note{
  This function cannot be used by users when the package is
  installed. For this reason it is not exported. If you
  want to maintain a local copy of the package and update
  annotation at will, you can download the package source.
}
\examples{
\dontrun{
library(metaseqr)
annotations.update()
}
}

