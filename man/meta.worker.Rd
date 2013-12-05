\name{meta.worker}
\alias{meta.worker}
\title{Permutation tests helper}
\usage{
  meta.worker(x,co,sl,cnt,s,r,sa,na,
    ll,el)
}
\arguments{
  \item{x}{a virtual list with the random seed and 
  the permutation index.}

  \item{co}{the counts matrix.}

  \item{sl}{the sample list.}

  \item{s}{the statistical algorithms.}

  \item{sa}{the parameters for each statistical
  algorithm.}

  \item{na}{the parameters for normalization.}

  \item{ll}{a list with library sizes.}

  \item{r}{same as the \code{replace} argument in 
  the \code{\link{sample}} function.}

  \item{weight}{a numeric vector of weights for each
  statistical algorithm}

  \item{el}{min, max or weight.}
}
\value{
  A matrix of p-values.
}
\description{
  This function performs the statistical test for 
  each permutation. Internal use only.
}
\examples{
\dontrun{
# Not yet available
}
}
\author{
  Panagiotis Moulos
}

