\name{diagplot.ftd}
\alias{diagplot.ftd}
\title{Create False (or True) Discovery curves}
\usage{
  diagplot.ftd(truth, p, type = "fdc", N = 2000,
    output = "x11", path = NULL, ...)
}
\arguments{
  \item{truth}{the ground truth differential 
  expression vector. It should contain only 
  zero and non-zero elements, with zero denoting
  non-differentially expressed genes and non-zero, 
  differentially expressed genes. Such a vector 
  can be obtained for example by using the 
  \code{\link{make.sim.data.sd}} function, which 
  creates simulated RNA-Seq read counts based on 
  real data. The elements of \code{truth} MUST 
  be named (e.g. each gene's name).}

  \item{p}{a p-value matrix whose rows correspond 
  to each element in the \code{truth} vector. If 
  the matrix has a \code{colnames} attribute, a 
  legend will be added to the plot using these 
  names, else a set of column names will be 
  auto-generated. \code{p} can also be a list or 
  a data frame. The p-values MUST be named (e.g. 
  each gene's name).}

  \item{type}{what to plot, can be \code{"fdc"} 
  for False Discovery Curves (default) or 
  \code{"tdc"} for True Discovery Curves.}

  \item{N}{create the curves based on the 
  top \code{N} ranked genes (default: 2000).}

  \item{output}{one or more R plotting device to 
  direct the plot result to. Supported mechanisms: 
  \code{"x11"} (default), \code{"png"}, \code{"jpg"}, 
  \code{"bmp"}, \code{"pdf"} or \code{"ps"}.}

  \item{path}{the path to create output files.}

  \item{...}{further arguments to be passed to plot
  devices, such as parameter from \code{\link{par}}.}
}
\value{
  A named list with two members: the first member 
  (\code{ftdr}) contains the values used to create 
  the plot. The second member (\code{path}) contains 
  the path to the created figure graphic.
}
\description{
  This function creates false (or true) discovery 
  curves using a matrix of p-values (such a matrix 
  can be derived for example from the result table 
  of \code{\link{metaseqr}} by subsetting the table 
  to get the p-values from several algorithms) 
  given a ground truth vector for differential 
  expression.
}
\examples{
\dontrun{
# Not yet available
}
}
\author{
  Panagiotis Moulos
}

