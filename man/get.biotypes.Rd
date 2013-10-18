\name{get.biotypes}
\alias{get.biotypes}
\title{Biotype converter}
\usage{
  get.biotypes(a)
}
\arguments{
  \item{a}{the annotation data frame (output of
  \code{\link{get.annotation}}).}
}
\value{
  A character vector of biotypes.
}
\description{
  Returns biotypes as character vector. Internal use.
}
\examples{
\dontrun{
hg19.bt <- get.biotypes(hg19.genes)
}
}
\author{
  Panagiotis Moulos
}

