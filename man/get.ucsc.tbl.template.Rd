\name{get.ucsc.tbl.template}
\alias{get.ucsc.tbl.template}
\title{Create SQLite UCSC table template defintions}
\usage{
    get.ucsc.tbl(tab, what="queries")
}
\arguments{
    \item{tab}{name of UCSC database table.}
    
    \item{what}{\code{"queries"} for SQLite table 
    definitions or \code{"fields"} for table column names.}
}
\value{
    An SQLite table definition.
}
\description{
    Returns an SQLIte table template defintion, 
    according to  UCSC Genome Browser database 
    table schemas. This functionality is used 
    when the package RMySQL is not available 
    for some reason, e.g. Windows machines. 
    Internal use only.
}
\examples{
\donttest{
db.table.tmpl <- get.ucsc.tbl.template("knownCanonical")
}
}
\author{
    Panagiotis Moulos
}

