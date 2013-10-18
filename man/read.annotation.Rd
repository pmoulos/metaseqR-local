\name{read.annotation}
\alias{read.annotation}
\title{Fixed annotation reader}
\usage{
  read.annotation(org, type)
}
\arguments{
  \item{org}{one of the supported organisms.}

  \item{type}{"gene" or "exon".}
}
\value{
  A data frame with the \code{type} annotation for
  \code{org}.
}
\description{
  A function to read fixed annotations from the local
  repository.
}
\examples{
\dontrun{
ann <- read.annotation("hg19","gene")
}
}
\author{
  Panagiotis Moulos
}

