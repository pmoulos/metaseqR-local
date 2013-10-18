\name{check.text.args}
\alias{check.text.args}
\title{Text argument validator}
\usage{
  check.text.args(arg.name, arg.value, arg.list,
    multiarg = FALSE)
}
\arguments{
  \item{arg.name}{the name of the argument that is checked
  (for display purposes).}

  \item{arg.value}{the value(s) of the argument to be
  checked.}

  \item{arg.list}{a vector of valid argument values for
  arg.value to be matched against.}

  \item{multiarg}{a logical scalar indicating whether
  arg.name accepts multiple arguments or not. In that case,
  all of the values in arg.value are checked against
  arg.list.}
}
\description{
  Checks if one or more given textual argument(s) is/are
  member(s) of a list of correct arguments. It's a more
  package-specific function similar to
  \code{\link{match.arg}}. Mostly for internal use.
}
\examples{
\dontrun{
check.text.args("count.type",count.type,c("gene","exon"),multiarg=FALSE)
check.text.args("statistics",statistics,c("deseq","edger","noiseq","bayseq","limma"),
  multiarg=TRUE)
}
}
\author{
  Panagiotis Moulos
}

