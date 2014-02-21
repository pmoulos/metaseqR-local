\name{get.weights}
\alias{get.weights}
\title{Get precalculated statistical test 
  weights}
\usage{
  get.weights(org = c("human", "mouse", "fruitfly", 
    "arabidopsis"), statistics = c("deseq", "edger", 
    "noiseq", "bayseq", "limma", "nbpseq"))
}
\arguments{
  \item{org}{\code{"human"}, \code{"mouse"}, 
  \code{"fruitfly"} or \code{"arabidopsis"}.}

  \item{statistics}{which statistical tests 
  are being used (see \code{link{metaseqr}}
  main help page).}
}
\value{
  A named vector of convex weights.
}
\description{
  This function returns pre-calculated weights 
  for human, mouse and fruitfly, based on the 
  performance of simulated datasets based on 
  real data from the ReCount database 
  (\url{http://bowtie-bio.sourceforge.net/recount/}).
}
\examples{
\dontrun{
wh <- get.weights("human",c("deseq","edger","noiseq"))
}
}
\author{
  Panagiotis Moulos
}

