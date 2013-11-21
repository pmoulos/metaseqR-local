\name{wapply}
\alias{wapply}
\title{List apply helper}
\usage{
  wapply(m, ...)
}
\arguments{
  \item{m}{a logical indicating whether to execute in
  parallel or not.}

  \item{...}{the rest arguments to \code{\link{lapply}} (or
  \code{mclapply})}
}
\description{
  A wrapper around normal and parallel apply
  (\code{\link{mclapply}} or multicore package) to avoid
  excessive coding for control of single or parallel code
  execution. Internal use.
}
\author{
  Panagiotis Moulos
}

