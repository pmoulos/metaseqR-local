#' @docType data
#' @name hg18.exon.counts
#' @title Human RNA-Seq data with two conditions, four samples
#' @description This data set contains RNA-Seq exon read counts for 3 chromosomes. The data are from an experiment studying the
#' effect of a long non-coding RNA related to the ASCL2 gene in WNT signaling and intestinal cancer. It has two conditions (CON, DOX)
#' and four samples (CON_BR1, CON_BR2, DOX_BR1, DOX_BR2). It also contains a predefined \code{sample.list} and \code{libsize.list}
#' named \code{sample.list.hg18} and \code{libsize.list.hg18}.
#' @usage hg18.exon.counts
#' @format a \code{data.frame} with exon read counts and some embedded annotation, one row per exon.
#' @source BSRC Alexander Fleming (http://www.fleming.gr)
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name sample.list.hg18
#' @title Human RNA-Seq data with two conditions, four samples
#' @description The sample list for \code{hg18.exon.counts}. See the data set description.
#' @usage sample.list.hg18
#' @format a named \code{list} with condition and sample names.
#' @source BSRC Alexander Fleming (http://www.fleming.gr)
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name libsize.list.hg18
#' @title Human RNA-Seq data with two conditions, four samples
#' @description The library size list for \code{hg18.exon.counts}. See the data set description.
#' @usage libsize.list.hg18
#' @format a named \code{list} with library sizes.
#' @source BSRC Alexander Fleming (http://www.fleming.gr)
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name mm9.gene.counts
#' @title mouse RNA-Seq data with three conditions, six samples
#' @description This data set contains RNA-Seq gene read counts for 3 chromosomes. The data are from an experiment studying gene
#' expression at several developmental stages of mouse liver. It has three conditions-developmental stages (e15.5, P0.5, P60) and
#' six samples (e15.5_1, e15.5_2, P0.5_1, P0.5_2, P60_1, P60_2). It also contains a predefined \code{sample.list} and \code{libsize.list}
#' named \code{sample.list.mm9} and \code{libsize.list.mm9}.
#' @usage mm9.gene.counts
#' @format a \code{data.frame} with gene read counts and some embedded annotation, one row per gene.
#' @source BSRC Alexander Fleming (http://www.fleming.gr)
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name sample.list.mm9
#' @title Mouse RNA-Seq data with three conditions, six samples
#' @description The sample list for \code{mm9.gene.counts}. See the data set description.
#' @usage sample.list.mm9
#' @format a named \code{list} with condition and sample names.
#' @source BSRC Alexander Fleming (http://www.fleming.gr)
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name libsize.list.mm9
#' @title Mouse RNA-Seq data with three conditions, six samples
#' @description The library size list for \code{mm9.gene.counts}. See the data set description.
#' @usage libsize.list.hg18
#' @format a named \code{list} with library sizes.
#' @source BSRC Alexander Fleming (http://www.fleming.gr)
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name hg18.gene
#' @title Human genome (genes) annotation data, build hg18
#' @description This data set contains predefined annotation for human genes, build hg18, if for any reason the user want to use
#' fixed annotation for metaseqr, instead of downloading (e.g. no internet connection or the biomaRt package is not available or
#' for speed). The format is a data frame with 8 columns: chromosome, start, end, gene_id, gc_content, strand, gene_name, biotype.
#' This data frame does not contain information about isoforms. All annotation data built with the package can be updated with the
#' \code{\link{annotations.update}} function.
#' @usage hg18.gene
#' @format a \code{data.frame} with 8 columns, one row per gene.
#' @source Ensembl
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name hg19.gene
#' @title Human genome (genes) annotation data, build hg19
#' @description This data set contains predefined annotation for human genes, build hg19, if for any reason the user want to use
#' fixed annotation for metaseqr, instead of downloading (e.g. no internet connection or the biomaRt package is not available or
#' for speed). The format is a data frame with 8 columns: chromosome, start, end, gene_id, gc_content, strand, gene_name, biotype.
#' This data frame does not contain information about isoforms. All annotation data built with the package can be updated with the
#' \code{\link{annotations.update}} function.
#' @usage hg19.gene
#' @format a \code{data.frame} with 8 columns, one row per gene.
#' @source Ensembl
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name mm9.gene
#' @title Mouse genome (genes) annotation data, build mm9
#' @description This data set contains predefined annotation for human genes, build mm9, if for any reason the user want to use
#' fixed annotation for metaseqr, instead of downloading (e.g. no internet connection or the biomaRt package is not available or
#' for speed). The format is a data frame with 8 columns: chromosome, start, end, gene_id, gc_content, strand, gene_name, biotype.
#' This data frame does not contain information about isoforms. All annotation data built with the package can be updated with the
#' \code{\link{annotations.update}} function.
#' @usage mm9.gene
#' @format a \code{data.frame} with 8 columns, one row per gene.
#' @source Ensembl
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name mm10.gene
#' @title Mouse genome (genes) annotation data, build mm10
#' @description This data set contains predefined annotation for human genes, build mm10, if for any reason the user want to use
#' fixed annotation for metaseqr, instead of downloading (e.g. no internet connection or the biomaRt package is not available or
#' for speed). The format is a data frame with 8 columns: chromosome, start, end, gene_id, gc_content, strand, gene_name, biotype.
#' This data frame does not contain information about isoforms. All annotation data built with the package can be updated with the
#' \code{\link{annotations.update}} function.
#' @usage mm10.gene
#' @format a \code{data.frame} with 8 columns, one row per gene.
#' @source Ensembl
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name rn5.gene
#' @title Rat genome (genes) annotation data, build rn5
#' @description This data set contains predefined annotation for human genes, build rn5, if for any reason the user want to use
#' fixed annotation for metaseqr, instead of downloading (e.g. no internet connection or the biomaRt package is not available or
#' for speed). The format is a data frame with 8 columns: chromosome, start, end, gene_id, gc_content, strand, gene_name, biotype.
#' This data frame does not contain information about isoforms. All annotation data built with the package can be updated with the
#' \code{\link{annotations.update}} function.
#' @usage rn5.gene
#' @format a \code{data.frame} with 8 columns, one row per gene.
#' @source Ensembl
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name dm3.gene
#' @title Fruitfly genome (genes) annotation data, build dm3
#' @description This data set contains predefined annotation for human genes, build dm3, if for any reason the user want to use
#' fixed annotation for metaseqr, instead of downloading (e.g. no internet connection or the biomaRt package is not available or
#' for speed). The format is a data frame with 8 columns: chromosome, start, end, gene_id, gc_content, strand, gene_name, biotype.
#' This data frame does not contain information about isoforms. All annotation data built with the package can be updated with the
#' \code{\link{annotations.update}} function.
#' @usage dm3.gene
#' @format a \code{data.frame} with 8 columns, one row per gene.
#' @source Ensembl
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name danRer7.gene
#' @title Zebrafish genome (genes) annotation data, build danRer7
#' @description This data set contains predefined annotation for human genes, build danRer7, if for any reason the user want to use
#' fixed annotation for metaseqr, instead of downloading (e.g. no internet connection or the biomaRt package is not available or
#' for speed). The format is a data frame with 8 columns: chromosome, start, end, gene_id, gc_content, strand, gene_name, biotype.
#' This data frame does not contain information about isoforms. All annotation data built with the package can be updated with the
#' \code{\link{annotations.update}} function.
#' @usage danRer7.gene
#' @format a \code{data.frame} with 8 columns, one row per gene.
#' @source Ensembl
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name hg18.exon
#' @title Human genome (exons) annotation data, build hg18
#' @description This data set contains predefined annotation for human exons, build hg18, if for any reason the user want to use
#' fixed annotation for metaseqr, instead of downloading (e.g. no internet connection or the biomaRt package is not available or
#' for speed). The format is a data frame with 8 columns: chromosome, start, end, exon_id, exon_id, strand, exon_name, biotype.
#' This data frame does not contain information about isoforms. All annotation data built with the package can be updated with the
#' \code{\link{annotations.update}} function.
#' @usage hg18.exon
#' @format a \code{data.frame} with 8 columns, one row per exon.
#' @source Ensembl
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name hg19.exon
#' @title Human genome (exons) annotation data, build hg19
#' @description This data set contains predefined annotation for human exons, build hg19, if for any reason the user want to use
#' fixed annotation for metaseqr, instead of downloading (e.g. no internet connection or the biomaRt package is not available or
#' for speed). The format is a data frame with 8 columns: chromosome, start, end, exon_id, exon_id, strand, exon_name, biotype.
#' This data frame does not contain information about isoforms. All annotation data built with the package can be updated with the
#' \code{\link{annotations.update}} function.
#' @usage hg19.exon
#' @format a \code{data.frame} with 8 columns, one row per exon.
#' @source Ensembl
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name mm9.exon
#' @title Mouse genome (exons) annotation data, build mm9
#' @description This data set contains predefined annotation for human exons, build mm9, if for any reason the user want to use
#' fixed annotation for metaseqr, instead of downloading (e.g. no internet connection or the biomaRt package is not available or
#' for speed). The format is a data frame with 8 columns: chromosome, start, end, exon_id, exon_id, strand, exon_name, biotype.
#' This data frame does not contain information about isoforms. All annotation data built with the package can be updated with the
#' \code{\link{annotations.update}} function.
#' @usage mm9.exon
#' @format a \code{data.frame} with 8 columns, one row per exon.
#' @source Ensembl
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name mm10.exon
#' @title Mouse genome (exons) annotation data, build mm10
#' @description This data set contains predefined annotation for human exons, build mm10, if for any reason the user want to use
#' fixed annotation for metaseqr, instead of downloading (e.g. no internet connection or the biomaRt package is not available or
#' for speed). The format is a data frame with 8 columns: chromosome, start, end, exon_id, exon_id, strand, exon_name, biotype.
#' This data frame does not contain information about isoforms. All annotation data built with the package can be updated with the
#' \code{\link{annotations.update}} function.
#' @usage mm10.exon
#' @format a \code{data.frame} with 8 columns, one row per exon.
#' @source Ensembl
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name rn5.exon
#' @title Rat genome (exons) annotation data, build rn5
#' @description This data set contains predefined annotation for human exons, build rn5, if for any reason the user want to use
#' fixed annotation for metaseqr, instead of downloading (e.g. no internet connection or the biomaRt package is not available or
#' for speed). The format is a data frame with 8 columns: chromosome, start, end, exon_id, exon_id, strand, exon_name, biotype.
#' This data frame does not contain information about isoforms. All annotation data built with the package can be updated with the
#' \code{\link{annotations.update}} function.
#' @usage rn5.exon
#' @format a \code{data.frame} with 8 columns, one row per exon.
#' @source Ensembl
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name dm3.exon
#' @title Fruitfly genome (exons) annotation data, build dm3
#' @description This data set contains predefined annotation for human exons, build dm3, if for any reason the user want to use
#' fixed annotation for metaseqr, instead of downloading (e.g. no internet connection or the biomaRt package is not available or
#' for speed). The format is a data frame with 8 columns: chromosome, start, end, exon_id, exon_id, strand, exon_name, biotype.
#' This data frame does not contain information about isoforms. All annotation data built with the package can be updated with the
#' \code{\link{annotations.update}} function.
#' @usage dm3.exon
#' @format a \code{data.frame} with 8 columns, one row per exon.
#' @source Ensembl
#' @author Panagiotis Moulos
NULL
#' @docType data
#' @name danRer7.exon
#' @title Zebrafish genome (exons) annotation data, build danRer7
#' @description This data set contains predefined annotation for human exons, build danRer7, if for any reason the user want to use
#' fixed annotation for metaseqr, instead of downloading (e.g. no internet connection or the biomaRt package is not available or
#' for speed). The format is a data frame with 8 columns: chromosome, start, end, exon_id, exon_id, strand, exon_name, biotype.
#' This data frame does not contain information about isoforms. All annotation data built with the package can be updated with the
#' \code{\link{annotations.update}} function.
#' @usage danRer7.exon
#' @format a \code{data.frame} with 8 columns, one row per exon.
#' @source Ensembl
#' @author Panagiotis Moulos
NULL
