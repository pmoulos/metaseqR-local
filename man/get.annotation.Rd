\name{get.annotation}
\alias{get.annotation}
\title{Annotation downloader}
\usage{
    get.annotation(org, type, refdb="ensembl",
        multic=FALSE)
}
\arguments{
    \item{org}{the organism for which to download
    annotation.}

    \item{type}{either \code{"gene"} or \code{"exon"}.}

    \item{refdb}{the online source to use to fetch 
    annotation. It can be \code{"ensembl"} (default), 
    \code{"ucsc"} or \code{"refseq"}. In the later two
    cases, an SQL connection is opened with the UCSC 
    public databases.}

    \item{multic}{a logical value indicating the 
    presence of multiple cores. Defaults to 
    \code{FALSE}. Do not change it if you are not 
    sure whether package parallel has been loaded 
    or not. It is used in the case of 
    \code{type="exon"} to process the return value of 
    the query to the UCSC Genome Browser database.}
}
\value{
    A data frame with the canonical (not isoforms!) genes or
    exons of the requested organism. When
    \code{type="genes"}, the data frame has the following
    columns: chromosome, start, end, gene_id, gc_content,
    strand, gene_name, biotype. When \code{type="exon"} the
    data frame has the following columns: chromosome, start,
    end, exon_id, gene_id, strand, gene_name, biotype. The
    gene_id and exon_id correspond to Ensembl gene and exon
    accessions respectively. The gene_name corresponds to
    HUGO nomenclature gene names.
}
\description{
    This function connects to the EBI's Biomart service using
    the package biomaRt and downloads annotation elements
    (gene co-ordinates, exon co-ordinates, gene
    identifications, biotypes etc.) for each of the supported
    organisms. See the help page of \code{\link{metaseqr}}
    for a list of supported organisms. The function downloads
    annotation for an organism genes or exons.
}
\note{
    The data frame that is returned contains only "canonical"
    chromosomes for each organism. It does not contain
    haplotypes or random locations and does not contain
    chromosome M.
}
\examples{
\donttest{
hg19.genes <- get.annotation("hg19","gene","ensembl")
mm9.exons <- get.annotation("mm9","exon","ucsc")
}
}
\author{
    Panagiotis Moulos
}

