#' The main metaseqr pipeline
#'
#' This function is the main metaseqr workhorse and implements the main metaseqr workflow which performs data read, filtering, normalization 
#' and statistical selection, creates diagnostic plots and exports the results and a report if requested.
#' The metaseqr function is responsible for assembling all the steps of the metaseqr pipeline which i) reads the input gene or exon
#' read count table ii) performs prelimininary filtering of data by removing chrM and other non-essential information for a typical
#' differential gene expression analysis as well as a preliminary expression filtering based on the exon counts, if an exon read
#' count file is provided. iii) performs data normalization with one of currently widely used algorithms, including EDASeq (Risso
#' et al., 2011), DESeq (Anders and Huber, 2010), edgeR (Robinson et al., 2010), NOISeq (Tarazona et al., 2012) or no normalization
#' iv) performs a second stage of filtering based on the normalized gene expression according to several gene filters v) performs
#' statistical testing with one or more of currently widely used algorithms, including DESeq (Anders and Huber, 2010), edgeR (Robinson
#' et al., 2010), NOISeq (Tarazona et al., 2012), limma (Smyth et al., 2005) for RNA-Seq data, baySeq (Hardcastle et al., 2012)
#' vi) in the case of multiple statistical testing algorithms, performs meta-analysis using one of five available methods (see the
#' meta.p argument) vii) exports the resulting differentially expressed gene list in text tab-delimited format viii) creates a set
#' of diagnostic plots either available in the aforementioned packages or metaseqr specific ones and ix) creates a comprehensive HTML
#' report which summarizes the run information, the results and the diagnostic plots. Certain diagnostic plots (e.g. the volcano plot)
#' can be interactive with the use of the external Highcharts (http://www.highcharts.com) JavaScript library for interactive graphs.
#' Although the inputs to the metaseqr workflow are many, in practice, setting only very few of them and accepting the defaults as
#' the rest can result in quite comprehensible results for mainstream organisms like mouse, human, fly and rat.
#'
#' @aliases metaseqr.main
#' @param counts a text tab-delimited file containing gene or exon counts in one of the following formats: i) the first column contains
#' unique gene or exon identifiers and the rest of the columns contain the read counts for each sample. Thus the first cell of each
#' row is a gene or exon accession and the rest are integers representing the counts for that accession. In that case, the \code{annotation}
#' parameter should strictly be \code{"fixed"}, \code{"download"} or an external file in proper format. ii) The first n columns should contain
#' gene or exon annotation elements like chromosomal locations, gene accessions, exon accessions, GC content etc. In that case, the
#' \code{annotation} parameter can also be \code{"embedded"}. The ideal embedded annotation contains 8 columns, chromosome, gene or exon start,
#' gene or exon end, gene or exon accession, GC-content (fraction or percentage), strand, HUGO gene symbol and gene biotype (e.g.
#' "protein_coding" or "ncRNA"). When the \code{annotation} parameter is "embedded", certain of these features are mandatory (co-ordinates
#' and accessions). If they are not present, the pipeline will not run. If additional elements are not present (e.g. GC content or
#' biotypes), certain features of metaseqr will not be available. For example, EDASeq normalization will not be performed based on
#' a GC content covariate but based on gene length which is not what the authors of EDASeq suggest. If biotypes are not present,
#' a lot of diagnostic plots will not be available. If the HUGO gene symbols are missing, the final annotation will contain only
#' gene accessions and thus be less comprehensible. Generally, it's best to set the \code{annotation} parameter to "download" or "fixed"
#' to ensure the most comprehensible results. Finally, counts can be a data frame satisfying the above conditions. It is a data
#' frame by default when \code{read2count} is used.
#' @param sample.list a list containing condition names and the samples under each condition. It should have the format \code{sample.list <-
#' list(ConditionA=c("Sample_A1", "Sample_A2", "Sample_A3"), ConditionB=c("Sample_B1", "Sample_B2"), ConditionC=c("Sample_C1", "Sample_C2"))}.
#' The names of the samples in list members MUST match the column names containing the read counts in the counts file. If they do
#' match, the pipeline will either crash or at best, ignore several of your samples. Alternative, \code{sample.list} can be a small
#' tab-delimited file structured as follows: he first line of the external tab delimited file should contain column names (names are
#' not important). The first column MUST contain UNIQUE sample names and the second column MUST contain the biological condition where
#' each of the samples in the first column should belong to. In this case, the function \code{\link{make.sample.list}} is used. If the
#' \code{counts} argument is missing, the \code{sample.list} argument MUST be a targets text tab-delimited file which contains the 
#' sample names, the BAM/BED file names and the biological conditions/groups for each sample/file. The file should be text tab-delimited
#' and structured as follows: the first line of the external tab delimited file should contain column names (names are not important).
#' The first column MUST contain UNIQUE sample names. The second column MUST contain the raw BAM/BED files WITH their full path.
#' Alternatively, the \code{path} argument should be provided (see below). The third column MUST contain the biological condition where
#' each of the samples in the first column should belong to.
#' @param path an optional path where all the BED/BAM files are placed, to be prepended to the BAM/BED file names in the targets file.
#' If not given and if the files in the second column of the targets file do not contain a path to a directory, the current directory
#' is assumed to be the BAM/BED file container.
#' @param file.type the type of raw input files. It can be \code{"auto"} for auto-guessing, \code{"bed"} for BED files, \code{"sam"}
#' for SAM files or \code{"bam"} for BAM files.
#' @param contrast a character vector of contrasts to be tested in the statistical testing step(s) of the metaseqr pipeline. Each 
#' element of the should STRICTLY have the format "ConditionA_vs_ConditionB_vs_...". A valid example based on the \code{sample.list}
#' above is \code{contrast <- c("ConditionA_vs_ConditionB", "ConditionA_vs_ConditionC", "ConditionA_vs_ConditionB_vs_ConditionC")}. The 
#' first  element of pairwise contrasts (e.g. "ConditionA" above) MUST be the control condition or any reference that ConditionB is 
#' checked  against. metaseqr uses this convention to properly calculate fold changes. If it's NULL, a contrast between the first two  
#' members of the \code{sample.list} will be auto-generated.
#' @param libsize.list an optional named list where names represent samples (MUST be the same as the samples in \code{sample.list}) and members
#' are the library sizes (the sequencing depth) for each sample. For example \code{libsize.list <- list(Sample_A1=32456913, Sample_A2=4346818)}.
#' @param id.col an integer denoting the column number in the file (or data frame) provided with the counts argument, where the unique
#' gene or exon accessions are. Default to 4 which is the standard feature name column in a BED file.
#' @param gc.col an integer denoting the column number in the file (or data frame) provided with the \code{counts} argument, where 
#' each gene's GC content is given. If not provided, GC content normalization provided by EDASeq will not be available.
#' @param name.col an integer denoting the column number in the file (or data frame) provided with the counts argument, where the 
#' HUGO gene symbols are given. If not provided, it will not be available when reporting results. In addition, the "known" gene
#' filter will not be available.
#' @param bt.col an integer denoting the column number in the file (or data frame) provided with the counts argument, where the 
#' gene biotypes are given. If not provided, the "biodetection", "countsbio", "saturation", "filtered" and "biodist" plots will not
#' be available.
#' @param annotation instructs metaseqr where to find the annotation for the given counts file. It can be one of i) "download" 
#' (default) for automatic downloading of the annotation for the organism specified by the org parameter (using biomaRt), ii) "fixed"
#' to retrieve the same annotation data as with "download" but from a fixed location inside the package (the "download" option always
#' download the latest annotation specified by the org parameter), iii) "embedded" if the annotation elements are embedded in the
#' read counts file or iv) a file specified by the user which should be as similar as possible to the "fixed" or "download" case, in
#' terms of column structure.
#' @param org the supported organisms by metaseqr. These can be, for human genomes "hg18" or "hg19", for mouse genomes "mm9", "mm10",
#' for rat genomes "rno5", for drosophila genomes "dm3" and for zebrafish genomes "danRer7".
#' @param count.type the type of reads inside the counts file. It can be one of "gene" or "exon". This is a very important and mandatory
#' parameter as it defines the course of the workflow.
#' @param exon.filters a named list whose names are the names of the supported exon filters and its members the filter parameters.
#' See section "Exon filters" below for details.
#' @param gene.filters a named list whose names are the names of the supported gene filters and its members the filter parameters.
#' See section "Gene filters" below for details.
#' @param normalization the normalization algorithm to be applied on the count data. It can be one of "edaseq" (default) for EDASeq
#' normalization, "deseq" for the normalization algorithm (individual options specified by the \code{norm.args} argument) in the DESeq
#' package, "edger" for the normalization algorithms present in the edgeR package (specified by the \code{norm.args} argument), "noiseq"
#' for the normalization algorithms present in the NOISeq package (specified by the \code{norm.args} argument), "nbpseq" for the
#' normalization algorithms present in the NBPSeq package (specified by the \code{norm.args} argument)  or "none" to not normalize 
#' the data (highly unrecommended).
#' @param norm.args a named list whose names are the names of the normalization algorithm parameters and its members parameter values.
#' See section "Normalization parameters" below for details. Leave NULL for the defaults of \code{normalization}.
#' @param statistics one or more statistical analyses to be performed by the metaseqr pipeline.It can be one or more of "deseq" (default)
#' to conduct statistical test(s) implemented in the DESeq package, "edger" to conduct statistical test(s) implemented in the edgeR
#' package, "limma" to conduct the RNA-Seq version of statistical test(s) implemented in the limma package, "noiseq" to conduct statistical
#' test(s) implemented in the NOISeq package, "bayseq" to conduct statistical test(s) implemented in the baySeq package and "nbpseq" to
#' conduct statistical test(s) implemented in the baySeq package In any case individual algorithm parameters are controlled by the
#' contents of the \code{stat.args} list.
#' @param stat.args a named list whose names are the names of the statistical algorithms used in the pipeline. Each member is another
#' named list whose names are the algorithm parameters and its members are the parameter values. See section "Statistics parameters"
#' below for details. Leave NULL for the defaults of \code{statistics}.
#' @param adjust.method the multiple testing p-value adjustment method. It can be one of \code{\link{p.adjust.methods}} or "qvalue"
#' from the qvalue Bioconductor package. Defaults to "BH" for Benjamini-Hochberg correction.
#' @param meta.p the meta-analysis method to combine p-values from multiple statistical tests. It can be one of "fisher" (default),
#' "perm", "whitlock", "intersection", "union" or "none". For the "fisher" and "perm" methods, see the documentation of the R package
#' MADAM. For the "whitlock" method, see the documentation of the survcomp Bioconductor package. With the "intersection" option, the
#' final p-value is the product of individual p-values derived from each method. However, the product is not used for the statistical
#' cutoff to derive gene lists. In this case, the final gene list is derived from the common differentially expressed genes from all
#' applied methods. Similarly, when meta.p is "union", the final list is derived from the union of individual methods and the final
#' p-values are the sum of individual p-values. The latter can be used as a very lose statistical threshold to aggregate results from
#' all methods regardless of their False Positive Rate.
#' @param pcut a p-value cutoff for exporting differentially genes, default is to export all.
#' @param log.offset an offset to be added to values during logarithmic transformations in order to avoid Infinity (default is 1).
#' @param preset an analysis strictness preset. Not yet implemented but in the end it should be a vector like c("strict","loose",
#' "verystrict","everything") etc.
#' @param qc.plots a set of diagnostic plots to show/create. It can be one or more of "mds", "biodetection", "countsbio", "saturation",
#' "rnacomp", "readnoise", "filtered", "boxplot", "gcbias", "lengthbias", "meandiff", "meanvar", "deheatmap", "volcano", "biodist". The "mds"
#' stands for Mutlti-Dimensional Scaling and it creates a PCA-like plot but using the MDS dimensionality reduction instead. It has
#' been succesfully used for NGS data (e.g. see the package htSeqTools) and it shows how well samples from the same condition cluster
#' together. For "biodetection", "countsbio", "saturation", "rnacomp", "readnoise", "biodist" see the vignette of NOISeq package. The "saturation"
#' case has been rewritten in order to display more samples in a more simple way. See the help page of \code{\link{diagplot.noiseq.saturation}}.
#' In addition, the "readnoise" plots represent an older version or the RNA composition plot included in older versions of NOISeq.
#' For "gcbias", "lengthbias", "meandiff", "meanvar" see the vignette of EDASeq package. "lenghtbias" is similar to "gcbias" but
#' using the gene length instead of the GC content as covariate. The "boxplot" option draws boxplots of log2 transformed gene counts.
#' The "filtered" option draws a 4-panel figure with the filtered genes per chromosome and per biotype, as absolute numbers and as
#' fractions of the genome. See also the help page of \code{\link{diagplot.filtered}}. The "deheatmap" option performs hierarchical
#' clustering and draws a heatmap of differentially expressed genes. In the context of diagnostic plots, it's useful to see if samples
#' from the same groups cluster together after statistical testing. The "volcano" option draws a volcano plot for each contrast and
#' if a report is requested, an interactive volcano plot is presented in the HTML report. Set to \code{NULL} if you don't want any
#' diagnostic plots created.
#' @param fig.format the format of the output diagnostic plots. It can be one or more of "x11" (for direct display), "png", "jpg",
#' "tiff", "bmp", "pdf", "ps".
#' @param out.list a logical controlling whether to export a list with the results in the running environment.
#' @param export.where  an output directory for the project results (report, lists, diagnostic plots etc.)
#' @param export.what the content of the final lists. It can be one or more of "annotation", to bind the annoation elements for each
#' gene, "p.value", to bind the p-values of each method, "adj.p.value", to bind the multiple testing adjusted p-values, "meta.p.value",
#' to bind the combined p-value from the meta-analysis, "adj.meta.p.value", to bind the corrected combined p-value from the meta-analysis,
#' "fold.change", to bind the fold changes of each requested contrast, "stats", to bind several statistics calclulated on raw and
#' normalized counts (see the \code{export.stats} argument), "counts", to bind the raw and normalized counts for each sample.
#' @param export.scale export values from one or more transformations applied to the data. It can be one or more of "natural", "log2",
#' "log10", "vst" (Variance Stabilizing Transormation, see the documentation of DESeq package).
#' @param export.values export raw and normalized counts.
#' @param export.stats calculate and export several statistics on raw and normalized counts, condition-wise. It can be one or more
#' of "mean", "median", "sd", "mad", "cv" for the Coefficient of Variation, "rcv" for a robust version of CV where the median and
#' the MAD are used instead of the mean and the standard deviation. 
#' @param restrict.cores in case of parallel execution of several subfunctions, the fraction of the available cores to use. In some 
#' cases if all available cores are used (\code{restrict.cores=1} and the system does not have sufficient RAM, the running machine 
#' might significantly slow down.
#' @param report a logical value controlling whether to produce a summary report or not. Defaults to TRUE.
#' @param report.template an HTML template to use for the report. Do not change this unless you know what you are doing.
#' @param verbose print informative messages during execution? Defaults to TRUE.
#' @param ... further arguments that may be passed to plotting functions, related to \code{\link{par}}.
#' @return If \code{out.list} is \code{TRUE}, a named list whose length is the same as the number of requested contrasts. Each list
#' member is named according to the corresponding contrast and contains a data frame of differentially expressed genes for that contrast.
#' The contents of the data frame are defined by the \code{export.what, export.scale, export.stats, export.values} parameters. If
#' \code{report} is \code{TRUE}, the output list contains two main elements. The first is described above (the analysis results)
#' and the second contains the same results but in HTML formatted tables.
#' @section Exon filters: The exon filters are a set of filters which are applied after the gene models are assembled from the read
#' counts of individual exons and before the gene expression is summarized from the exons belonging to each gene. These filters can
#' be applied when the input read counts file contains exon reads. It is not applicable when the input file already contains gene
#' counts. Such filters can be for example "accept genes where all the exons contain more than x reads" or "accept genes where there
#' is read presence in at least m/n exons, n being the total exons of the gene". Such filters are NOT meant for detecting differential
#' splicing as also the whole metaseqr pipeline, this they should not be used in that context. The \code{exon.filters} argument is a
#' named list of filters, where the names are the filter names and the members are the filter parameters (named lists with parameter
#' name, parameter value). See the usage of the function for an example of how these lists are structured. The supported exon filters
#' in the current version are: i) \code{min.active.exons} which implements a filter for demanding m out of n exons of a gene to have a 
#' certain read presence with parameters \code{exons.per.gene}, \code{min.exons} and \code{frac}. The filter is described as follows: if
#' a gene has up to \code{exons.per.gene} exons, then read presence is required in at least \code{min.exons} of them, else read presence
#' is required in a \code{frac} fraction of the total exons. With the default values, the filter instructs that if a gene has up to 5
#' exons, read presence is required in at least 2, else in at least 20% of the exons, in order to be accepted. More filters will be
#' implemented in future versions and users are encouraged to propose exon filter ideas to the author by mail. See \code{metaseqr}
#' usage for the defaults. Set exon.filters=NULL to not apply any exon filtering.
#' @section Gene filters: The gene filters are a set of filters applied to gene expression as this is manifested through the read
#' presence on each gene and are preferably applied after normalization. These filters can be applied both when the input file or
#' data frame contains exon read counts and gene read counts. Such filter can be for example "accept all genes above a certain count
#' threshold" or "accept all genes with expression above the median of the normalized counts distribution" or "accept all with length
#' above a certain threshold in kb" or "exclude the 'pseudogene' biotype from further analysis". The supported gene filters in the
#' current version, which have the same structure as the exon filters (named list of lists with filter names, parameter names and
#' parameter arguments)  are: i) \code{length} which implements a length filter where genes are accepted for further analysis if they 
#' are above \code{length} (its parameter) kb. ii) \code{avg.reads} which implements a filter where a gene is accepted for further
#' analysis if it has more average reads than the \code{quantile} of the average count distribution per \code{average.per.bp} base
#' pairs. In summary, the reads of each gene are averaged per \code{average.per.bp} based on each gene's length (in case of exons,
#' input the "gene's length" is the sum of the lengths of exons) and the \code{quantile} quantile of the average counts distribution
#' is calculated for each sample. Genes passing the filter should have an average read count larger than the maximum of the vector
#' of the quantiles calculated above. iii) \code{expression} which implements a filter based on the overall expression of a gene.
#' The parameters of this filter are: \code{median}, where genes below the median of the overall count distribution are not accepted
#' for further analysis (this filter has been used to distinguish between "expressed" and "not expressed" genes in several cases, e.g.
#' (Mokry et al., 2011) with a logical as value, \code{mean} which is the same as \code{median} but using the mean, \code{quantile}
#' which is the same as the previous two but using a specific quantile of the total counts distribution, \code{known}, where in this
#' case, a set of known not-expressed genes in the system under investigation are used to estimate an expression cutoff. This can be
#' quite useful, as the genes are filtered based on a "true biological" cutoff instead of a statistical cutoff. The value of this
#' filter is a character vector of HUGO gene symbols (MUST be contained in the annotation, thus it's better to use \code{annotation=
#' "fixed"} or \code{annotation="download"}) whose counts are used to build a "null" expression distribution. The 90th quantile of
#' this distribution is then the expression cutoff. This filter can be combined with any other filter. Be careful with gene names
#' as they are case sensitive and must match exactly ("Pten" is different than "PTEN"!). iv) \code{biotype} where in this case, genes
#' with a certain biotype (MUST be contained in the annotation, thus it's better to use \code{annotation="fixed"} or \code{annotation=
#' "download"}) are excluded from the analysis. This filter is a named list of logical, where names are the biotypes in each genome
#' and values are \code{TRUE} or \code{FALSE}. If the biotype should be excluded, the value should be \code{TRUE} else \code{FALSE}.
#' See the result of \code{get.defaults("biotype.filter","hg19")} for an example. Finally, in future versions there will be support
#' for user-defined filters in the form of a function.
#' @section Normalization parameters: The normalization parameters are passed again as a named list where the names of the members
#' are the normalization parameter names and the values are the normalization parameter values. You should check the documentation
#' of the packages EDASeq, DESeq, edgeR, NOISeq and NBPSeq for the parameter names and parameter values. There are a few two exceptions
#' in parameter names: in case of \code{normalization="edaseq"} the only parameter names are \code{within.which} and \code{between.which}, 
#' controlling the withing lane/sample and between lanes/samples normalization algorithm. In case of \code{normalization="edger"}, 
#' apart from the rest of the edgeR normalization arguments, there is the argument \code{main.method} which can be either "classic" 
#' or "glm" (see the edgeR's manual for details, briefly these are different algorithms for estimating sample dispersion parameters and
#' normalization factors), and \code{norm.method} which controls the normalization method and replaces the \code{method} parameter in
#' respective edgeR's calls (again see edgeR's manual). For the rest of the algorithms, the parameter names are the same as the
#' names used in the respective packages. For examples, please use the \code{\link{get.defaults}} function.
#' @section Statistics parameters: The statistics parameters as passed to statistical algorithms in metaseqr, exactly with the same
#' way as the normalization parametes above. In this case, there is one more layer in list nesting. Thus, \code{stat.args} is a named
#' list whose names are the names the algorithms used (see the \code{statistics} parameter). Each member is another named list,with
#' parameters to be used for each statistical algorithm. Again, the names of the member lists are parameter names and the values of
#' the member lists are parameter values. You should check the documentations of DESeq, edgeR, NOISeq, baySeq, limma and NBPSeq for 
#' these parameters. There are a few exceptions in parameter names: In case of \code{statistics="edger"}, apart from the rest of the 
#' edgeR statistical testing arguments, there is the argument \code{main.method} which can be either "classic" or "glm", again defining
#' whether the binomial test or GLMs will be used for statistical testing. For examples, please use the \code{\link{get.defaults}}
#' function. When \code{statistics="nbpseq"}, apart from the rest arguments of the NBPSeq functions \code{estimate.disp} and
#' \code{estimate.dispersion}, there is the argument \code{main.method} which can be \code{"nbpseq"} or \code{"nbsmyth"}. This
#' argument determines the parameters to be used by the \code{estimate.dispersion} function or by the \code{estimate.disp} function
#' to estimate RNA-Seq count dispersions. The difference between the two is that they constitute different starting points for two
#' workflows in the package NBPSeq. The first worklfow (with \code{main.method="nbpseq"} and the \code{estimate.dispersion} function
#' is NBPSeq package specific, while the second (with \code{main.method="nbsmyth"} and the \code{estimate.disp} function is similar
#' to the workflow of the edgeR package. For additional information regarding the statistical testing in NBPSeq, please consult
#' the documentation of the NBPSeq package.
#' @note Please note that currently only gene and exon annotation from Ensembl (http://www.ensembl.org) are supported. Thus, the
#' unique gene or exon ids in the counts files should correspond to valid Ensembl gene or exon accessions for the organism of interest.
#' If you are not sure about the source of your counts file or do not know how to produce it, it's better to start from the original
#' BAM files and run the pipeline through the \code{\link{read2count}} wrapper. Keep in mind that in this case the performance will
#' be significantly lower and the overall running time significanlty higher as the R functions which are used to read BAM files to 
#' proper structures (GenomicRanges) and calculate the counts are quite slow. An alternative way is maybe the easyRNASeq package
#' (Delhomme et al, 2012). The read2count function does not use this package but rather makes use of standard Bioconductor functions
#' to handle NGS data. If you wish to work outside R, you can work with other popular read counters such as the HTSeq read counter
#' (http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html). Please also note that in the current version, the members of the
#' gene.filters and exon.filters lists are not checked for validity so be careful to supply with correct names otherwise the pipeline
#' will crash or at the best case scenario, will ignore the filters. Also note that when you are supplying metaseqr wtih an exon counts
#' table, gene annotation is always downloaded or read from a fixed location inside the package. In addition to the above, if you
#' have a multiple core system, be very careful on how you are using the \code{restrict.cores} argument and generally how many cores
#' you are using with scripts purely written in R. The analysis with exon read data can very easily cause memory problems, so unless
#' you have more than 64Gb of RAM available, consider setting restrict.cores to something like 0.2 when working with exon data.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(metaseqr)
#'
#' # An example pipeline with exon counts
#' data("hg18.exon.data",package="metaseqr")
#' metaseqr(
#'  counts=hg18.exon.counts,
#'  sample.list=list(CON=c("CON_BR1","CON_BR2"),DOX=c("DOX_BR1","DOX_BR2")),
#'  contrast=c("CON_vs_DOX"),
#'  libsize.list=list(CON_BR1=17041268,CON_BR2=23579904,DOX_BR1=16018639,DOX_BR2=26294259),
#'  id.col=4,
#'  annotation="download",
#'  org="hg18",
#'  count.type="exon",
#'  normalization="edaseq",
#'  statistics="deseq",
#'  pcut=0.05,
#'  qc.plots=c("mds", "biodetection", "countsbio", "saturation", "rnacomp", "boxplot", "gcbias", "lengthbias", "meandiff",
#'    "meanvar", "readnoise", "deheatmap", "volcano", "biodist", "filtered"),
#'  fig.format=c("png","pdf"),
#'  export.what=c("annotation","p.value","adj.p.value","fold.change","stats","counts"),
#'  export.scale=c("natural","log2","log10","vst"),
#'  export.values=c("raw","normalized"),
#'  export.stats=c("mean","median","sd","mad","cv","rcv"),
#'  restrict.cores=0.8,
#'  gene.filters=list(
#'    length=list(
#'      length=500
#'    ),
#'    avg.reads=list(
#'      average.per.bp=100,
#'      quantile=0.25
#'    ),
#'    expression=list(
#'      median=TRUE,
#'      mean=FALSE
#'    ),
#'    biotype=get.defaults("biotype.filter","hg18")
#'  )
#' )
#'
#' # An example pipeline with gene counts
#' data("mm9.gene.data",package="metaseqr")
#' result <- metaseqr(
#'  counts=mm9.gene.counts,
#'  sample.list=list(e15.5=c("e15.5_1","e15.5_2"), P0.5=c("P0.5_1","P0.5_2"), P60=c("P60_1","P60_2")),
#'  contrast=c("e15.5_vs_P0.5","e15.5_vs_P60","P0.5_vs_P60","e15.5_vs_P0.5_vs_P60"),
#'  libsize.list=list(e15.5_1=46546134, e15.5_2=18439760, P0.5_1=21822789, P0.5_2=40813977,
#'    P60_1=35855191, P60_2=43611778),
#'  annotation="fixed",
#'  org="mm9",
#'  count.type="gene",
#'  normalization="edger",
#'  statistics=c("deseq","edger","noiseq"),
#'  meta.p="fisher",
#'  pcut=0.05,
#'  fig.format=c("png","pdf"),
#'  export.what=c("annotation","p.value","meta.p.value","adj.meta.p.value","fold.change"),
#'  export.scale=c("natural","log2"),
#'  export.values="normalized",
#'  export.stats=c("mean","sd","cv"),
#'  export.where=getwd(),
#'  restrict.cores=0.8,
#'  gene.filters=list(
#'    length=list(
#'      length=500
#'    ),
#'    avg.reads=list(
#' 	    average.per.bp=100,
#' 	    quantile=0.25
#'    ),
#'    expression=list(
#' 	    median=TRUE,
#' 	    mean=FALSE,
#' 	    quantile=NA,
#' 	    known=NA,
#' 	    custom=NA
#'    ),
#'    biotype=get.defaults("biotype.filter","mm9")
#'  ),
#'  out.list=TRUE
#' )
#' head(result$data[["e15.5_vs_P0.5"]])
#' }
metaseqr <- function(
	counts,
	sample.list,
	file.type=c("auto","sam","bam","bed"),
	path=NULL,
	contrast=NULL,
	libsize.list=NULL,
	id.col=4,
	gc.col=NA,
	name.col=NA,
	bt.col=NA,
	annotation=c("download","embedded","fixed"), # It can also be a file at some point, for speed in the final deployment
	org=c("hg18","hg19","mm9","mm10","rno5","dm3","danRer7"),
	count.type=c("gene","exon"),
	exon.filters=list(
		min.active.exons=list(
			exons.per.gene=5,
			min.exons=2,
			frac=1/5
		)
	),
	gene.filters=list(
		length=list(
			length=500
		),
		avg.reads=list(
			average.per.bp=100,
			quantile=0.75
		),
		expression=list(
			median=TRUE,
			mean=FALSE,
			quantile=NA,
			known=NA,
			custom=NA
		),
		biotype=get.defaults("biotype.filter",org[1])
	),
	normalization=c("edaseq","deseq","edger","noiseq","nbpseq","none"),
	norm.args=NULL,
	statistics=c("deseq","edger","noiseq","bayseq","limma","nbpseq"),
	stat.args=NULL,
	adjust.method=sort(c(p.adjust.methods,"qvalue")), # Brings BH first which is the default
	meta.p=if (length(statistics)>1) c("fisher","perm","whitlock","intersection","union","none") else "none",
	pcut=NA, # A p-value cutoff for exporting DE genes, default is to export all
	log.offset=1, # Logarithmic transformation offset to avoid +/-Inf (log2(a+offset/b+offset))
	preset=NULL, # In the end it should be a vector like c("strict","loose","verystrict","everything") etc.
	qc.plots=c(
		"mds","biodetection","countsbio","saturation","readnoise","filtered", # Raw count data
		"boxplot","gcbias","lengthbias","meandiff","meanvar","rnacomp", # Pre and post normalization
		"deheatmap","volcano","biodist" # Post statistical testing
	),
	fig.format=c("x11","png","jpg","tiff","bmp","pdf","ps"),
	out.list=FALSE,
	export.where=NA, # An output directory for the project
	export.what=c("annotation","p.value","adj.p.value","meta.p.value","adj.meta.p.value","fold.change","stats","counts"),
	export.scale=c("natural","log2","log10","vst"),
	export.values=c("raw","normalized"),
	export.stats=c("mean","median","sd","mad","cv","rcv"),
	restrict.cores=0.6,
	report=TRUE,
	report.template="default",
	verbose=TRUE,
	...
)

{
	# Check essential arguments
	from.raw <- FALSE
	if (missing(counts) && (missing(sample.list) || is.list(sample.list)))
		stop(paste("You must provide a file with genomic region (gene, exon, etc.) counts or an input targets file to create input from!",
			"If the counts file is missing, sample.list cannot be missing or a list! It must be a targets file with at least three columns!",
			"See the read.targets function."))
	if (missing(sample.list) || (!is.list(sample.list) && !file.exists(sample.list)))
		stop("You must provide a list with condition names and sample names (same as in the counts file) or an input file to create the sample list from!")
	if (!is.list(sample.list) && file.exists(sample.list) && !missing(counts))
		sample.list <- make.sample.list(sample.list)
	if (!is.list(sample.list) && file.exists(sample.list) && missing(counts))
	{
		counts <- NULL
		the.list <- read.targets(sample.list,path=path)
		sample.list <- the.list$samples
		file.list <- the.list$files
		if (tolower(file.type[1])=="auto")
			file.type <- the.list$type
		if (is.null(file.type))
			stop("The type of the input files could not be recognized! Please specify (BAM or BED)...")
		from.raw <- TRUE
	}

	# Initialize environmental variables
	if (!exists("HOME"))
		init.envar()
	
	# Globalize the project's path and verbosity
	if (from.raw)
		PROJECT.PATH <<- make.project.path(export.where)
	else
		PROJECT.PATH <<- make.project.path(export.where,counts)
	VERBOSE <<- verbose

	file.type <- tolower(file.type[1])
	annotation <- tolower(annotation[1])
	org <- tolower(org[1])
	count.type <- tolower(count.type[1])
	normalization <- tolower(normalization[1])
	adjust.method <- adjust.method[1]
	meta.p <- tolower(meta.p[1])
	statistics <- tolower(statistics)
	fig.format <- tolower(fig.format)
	qc.plots <- tolower(qc.plots)
	export.what <- tolower(export.what)
	export.scale <- tolower(export.scale)
	export.values <- tolower(export.values)
	export.stats <- tolower(export.stats)

	if (!is.data.frame(counts) && !is.null(counts))
	{
		check.file.args("counts",counts)
		counts.name <- basename(counts)
	}
	else
	{
		counts.name <- "imported custom data frame"
	}

	check.text.args("file.type",file.type,c("auto","sam","bam","bed"),multiarg=FALSE)
	check.text.args("annotation",annotation,c("embedded","download","fixed"),multiarg=FALSE)
	check.text.args("org",org,c("hg18","hg19","mm9","mm10","rno5","dm3","danRer7"),multiarg=FALSE)
	check.text.args("count.type",count.type,c("gene","exon"),multiarg=FALSE)
	check.text.args("normalization",normalization,c("edaseq","deseq","edger","noiseq","nbpseq","none"),multiarg=FALSE)
	check.text.args("statistics",statistics,c("deseq","edger","noiseq","bayseq","limma","nbpseq"),multiarg=TRUE)
	check.text.args("meta.p",meta.p,c("fisher","perm","whitlock","intersection","union","none"),multiarg=FALSE)
	check.text.args("fig.format",fig.format,c("x11","png","jpg","tiff","bmp","pdf","ps"),multiarg=TRUE)
	check.text.args("export.what",export.what,c("annotation","p.value","adj.p.value","meta.p.value","adj.meta.p.value","fold.change","stats","counts"),multiarg=TRUE)
	check.text.args("export.scale",export.scale,c("natural","log2","log10","vst"),multiarg=TRUE)
	check.text.args("export.values",export.values,c("raw","normalized"),multiarg=TRUE)
	check.text.args("export.stats",export.stats,c("mean","median","sd","mad","cv","rcv"),multiarg=TRUE)
	if (!is.null(qc.plots))
		check.text.args("qc.plots",qc.plots,c("mds","biodetection","countsbio","saturation","readnoise","boxplot","gcbias","lengthbias",
			"meandiff","meanvar","rnacomp","deheatmap","volcano","biodist","filtered"),multiarg=TRUE)
	if (!is.na(restrict.cores)) check.num.args("restrict.cores",restrict.cores,"numeric",c(0,1),"botheq")
	if (!is.na(pcut)) check.num.args("pcut",pcut,"numeric",c(0,1),"botheq")
	if (!is.na(gc.col)) check.num.args("gc.col",gc.col,"numeric",0,"gt")
	if (!is.na(name.col)) check.num.args("name.col",name.col,"numeric",0,"gt")
	if (!is.na(bt.col)) check.num.args("bt.col",bt.col,"numeric",0,"gt")
	if (!is.na(log.offset)) check.num.args("log.offset",log.offset,"numeric",0,"gt")
	if (!is.null(contrast)) check.contrast.format(contrast,sample.list)
	if ("bayseq" %in% statistics) libsize.list <- check.libsize(libsize.list,sample.list)

	# Check main functionality packages
	check.packages(normalization,statistics,adjust.method,meta.p,export.scale,qc.plots,report)
	# Check if parallel processing is available
	multic <- check.parallel(restrict.cores)
	# Check the case of embedded annotation but not given gc and gene name columns
	if (annotation=="embedded")
	{
		if (is.na(gc.col) && count.type=="gene")
			stop("The column that contains the gene GC content (\"gc.col\") argument is required when \"annotation\" is \"embedded\"!")
		if (is.na(name.col) && !is.na(gene.filters$expression$known))
		{
			warning("The column that contains the HUGO gene symbols (\"bt.col\") is missing with embedded annotation! Gene name expression filter will not be available...",
				call.=FALSE)
			gene.filters$expression$known=NA
			if ("volcano" %in% qc.plots)
				warning("The column that contains the HUGO gene symbols (\"bt.col\") is missing with embedded annotation! Interactive volcano plots will not contain gene names...",
					call.=FALSE)
		}
		if (is.na(bt.col) && count.type=="gene")
		{
			warning("The column that contains the gene biotypes (\"bt.col\") is missing with embedded annotation! Biotype filters and certain plots will not be available...",
				call.=FALSE)
			gene.filters$biotype=NULL
			to.remove <- match(c("biodetection","countsbio","saturation","biodist","filtered"),qc.plots)
			no.match <- which(is.na(to.remove))
			if (length(no.match)>0)
				to.remove <- to.remove[-no.match]
			qc.plots <- qc.plots[-to.remove]
		}
	}
	else if (annotation=="download" || count.type=="exon") # Requires package biomaRt
	{
		if (!require(biomaRt))
			stop("Bioconductor package biomaRt is required when annotation is \"download\" or type argument is \"exon\"!")
	}

	# Check additional input arguments for normalization and statistics
	if (!is.null(norm.args))
	{
		tmp <- norm.args
		norm.args <- get.defaults("normalization",normalization)
		norm.args <- set.arg(norm.args,tmp)
	}
	else
		norm.args <- get.defaults("normalization",normalization)
	for (s in statistics)
	{
		if (!is.null(stat.args[[s]]))
		{
			tmp <- stat.args[[s]]	
			stat.args[[s]] <- get.defaults("statistics",s)
			stat.args[[s]] <- set.arg(stat.args[[s]],tmp)
		}
		else
			stat.args[[s]] <- get.defaults("statistics",s)
	}
	# Override settigs if a preset is given
	if (!is.null(preset))
	{
		# Override filter rules and maybe norm.args and stat.args
	}

	if (report)
	{
		report.messages <- make.report.messages("en")
		if (!is.null(qc.plots) && !("png" %in% fig.format))
		{
			warning("png format is required in order to build a report! Adding to figure output formats...",call.=FALSE)
			fig.format <- c(fig.format,"png")
		}
	}

	# Display initialization report
	disp(strftime(Sys.time()),": Data processing started...\n")
	##############################################################################################################################
	disp("Read counts file: ",counts.name)
	disp("Conditions: ",paste(names(sample.list),collapse=", "))
	disp("Samples: ",paste(unlist(sample.list),collapse=", "))
	disp("Requested contrasts: ",paste(contrast,collapse=", "))
	if (!is.null(libsize.list))
	{
		disp("Library sizes: ")
		for (n in names(libsize.list))
			disp("  ",paste(n,libsize.list[[n]],sep=": "))
	}
	disp("Annotation: ",annotation)
	disp("Organism: ",org)
	disp("Count type: ",count.type)
	if (!is.null(exon.filters))
	{
		disp("Exon filters: ",paste(names(exon.filters),collapse=", "))
		for (ef in names(exon.filters))
		{
			disp("  ",ef,": ")
			for (efp in names(exon.filters[[ef]]))
			{
				if (length(exon.filters[[ef]][[efp]])==1 && is.function(exon.filters[[ef]][[efp]]))
					print(exon.filters[[ef]][[efp]])
				else if (length(exon.filters[[ef]][[efp]])==1)
					disp("    ",paste(efp,exon.filters[[ef]][[efp]],sep=": "))
				else if (length(exon.filters[[ef]][[efp]])>1)
					disp("    ",paste(efp,paste(exon.filters[[ef]][[efp]],collapse=", "),sep=": "))
			}
		}
	}
	if (!is.null(gene.filters))
	{
		disp("Gene filters: ",paste(names(gene.filters),collapse=", "))
		for (gf in names(gene.filters))
		{
			disp("  ",gf,": ")
			for (gfp in names(gene.filters[[gf]]))
			{
				if (length(gene.filters[[gf]][[gfp]])==1 && is.function(gene.filters[[gf]][[gfp]]))
					print(gene.filters[[gf]][[gfp]])
				else if (length(gene.filters[[gf]][[gfp]])==1)
					disp("    ",paste(gfp,gene.filters[[gf]][[gfp]],sep=": "))
				else if (length(gene.filters[[gf]][[gfp]])>1)
					disp("    ",paste(gfp,paste(gene.filters[[gf]][[gfp]],collapse=", "),sep=": "))
			}
		}
	}
	disp("Normalization algorithm: ",normalization)
	if (!is.null(norm.args))
	{
		disp("Normalization arguments: ")
		for (na in names(norm.args))
		{
			if (length(norm.args[[na]])==1 && is.function(norm.args[[na]]))
			{
				disp("  ",na,": ")
				print(norm.args[[na]])
			}
			else if (length(norm.args[[na]])==1)
				disp("  ",paste(na,norm.args[[na]],sep=": "))
			else if (length(norm.args[[na]])>1)
				disp("  ",paste(na,paste(norm.args[[na]],collapse=", "),sep=": "))
		}
	}
	disp("Statistical algorithm: ",paste(statistics,collapse=", "))
	if (!is.null(stat.args))
	{
		disp("Statistical arguments: ")
		for (sa in names(stat.args))
		{
			if (length(stat.args[[sa]])==1 && is.function(stat.args[[sa]]))
			{
				disp("  ",sa,": ")
				print(stat.args[[sa]])
			}
			else if (length(stat.args[[sa]])==1)
				disp("  ",paste(sa,stat.args[[sa]],sep=": "))
			else if (length(stat.args[[sa]])>1)
				disp("  ",paste(sa,paste(stat.args[[sa]],collapse=", "),sep=": "))
		}
	}
	disp("Meta-analysis method: ",meta.p)
	disp("Multiple testing correction: ",adjust.method)
	if (!is.na(pcut)) disp("p-value threshold: ",pcut)
	disp("Logarithmic tranformation offset: ",log.offset)
	if (!is.null(preset)) disp("Analysis preset: ",preset)
	disp("Quality control plots: ",paste(qc.plots,collapse=", "))
	disp("Figure format: ",paste(fig.format,collapse=", "))
	if (!is.na(export.where)) disp("Output directory: ",export.where)
	disp("Output data: ",paste(export.what,collapse=", "))
	disp("Output scale(s): ",paste(export.scale,collapse=", "))
	disp("Output values: ",paste(export.values,collapse=", "))
	disp("Output statistics: ",paste(export.stats,collapse=", "))
	disp("")
	##############################################################################################################################

	if (count.type=="exon")
	{
		if (annotation=="download")
		{
			disp("Downloading gene annotation for ",org,"...")
			gene.data <- get.annotation(org,"gene")
		}
		else
		{
			disp("Reading stored gene annotation for ",org,"...")
			gene.data <- read.annotation(org,"gene")
		}
	
		if (annotation=="download")
		{
			disp("Downloading exon annotation for ",org,"...")
			exon.data <- get.annotation(org,count.type)
		}
		else if (annotation=="fixed")
		{
			disp("Reading stored exon annotation for ",org,"...")
			exon.data <- read.annotation(org,count.type)
		}
		else if (annotation=="embedded") # The following should work if annotation elements are arranged in MeV-like data style
		{
			# Embedded annotation can NEVER occur when receiving data from read2count, so there is no danger here
			if (!is.data.frame(counts))
			{
				disp("Reading counts file ",counts.name,"...")
				exon.counts <- read.delim(counts)
			}
			else
				exon.counts <- counts
			rownames(exon.counts) <- as.character(exon.counts[,id.col])
			all.cols <- 1:ncol(exon.counts)
			sam.cols <- match(unlist(sample.list),colnames(exon.counts))
			sam.cols <- sam.cols[which(!is.na(sam.cols))]
			ann.cols <- all.cols[-sam.cols]
			exon.data <- exon.counts[,ann.cols]
			exon.counts <- exon.counts[,sam.cols]
			colnames(exon.data)[id.col] <- "exon_id"
			if (!is.na(name.col)) colnames(exon.data)[name.col] <- "gene_name"
			if (!is.na(bt.col)) colnames(exon.data)[bt.col] <- "biotype"
			exon.counts <- cbind(exon.data[rownames(exon.counts),c("start","end","exon_id","gene_id")],exon.counts)
		}
		else # Reading from external file, similar to embedded
		{
			disp("Reading external exon annotation for ",org," from ",annotation,"...")
			exon.data <- read.delim(annotation)
			colnames(exon.data)[id.col] <- "exon_id"
		}

		if (annotation!="embedded") # Else everything is provided and done
		{
			if (!is.null(counts)) # Otherwise it's coming ready from read2count
			{
				if (!is.data.frame(counts)) # Else it's already here
				{
					disp("Reading counts file ",counts.name,"...")
					exon.counts <- read.delim(counts)
				}
				else # Already a data frame as input
					exon.counts <- counts
				rownames(exon.counts) <- as.character(exon.counts[,id.col])
				exon.counts <- exon.counts[,unlist(sample.list,use.names=FALSE)]
			}
			else # Coming from read2count
			{
				if (from.raw) # Double check
				{
					r2c <- read2count(file.list,file.type,exon.data)
					exon.counts <- r2c$counts
					if (is.null(libsize.list))
						libsize.list <- r2c$libsize
				}
			}
		}
		exon.counts <- cbind(exon.data[rownames(exon.counts),c("start","end","exon_id","gene_id")],exon.counts[,unlist(sample.list,use.names=FALSE)])

		# Get the exon counts per gene model
		disp("Checking chromosomes in exon counts and gene annotation...")
		gene.data <- reduce.gene.data(exon.data[rownames(exon.counts),],gene.data)
		disp("Processing exons...")
		the.counts <- construct.gene.model(exon.counts,sample.list,gene.data,multic=multic)

		# Apply exon filters
		if (!is.null(exon.filters))
			exon.filter.result <- filter.exons(the.counts,gene.data,sample.list,exon.filters)
		else
			exon.filter.result <- NULL
		
		disp("Summarizing count data...")
		the.gene.counts <- the.exon.lengths <- vector("list",length(unlist(sample.list)))
		names(the.gene.counts) <- names(the.exon.lengths) <- names(the.counts)
		for (n in names(the.gene.counts))
		{
			the.gene.counts[[n]] <- wapply(multic,the.counts[[n]],function(x) return(sum(x$count)))
			the.exon.lengths[[n]] <- wapply(multic,the.counts[[n]],function(x) return(sum(x$length)))
			the.gene.counts[[n]] <- do.call("c",the.gene.counts[[n]])
			the.exon.lengths[[n]] <- do.call("c",the.exon.lengths[[n]])
		}
		gene.counts <- do.call("cbind",the.gene.counts)
		gene.length <- the.exon.lengths[[1]] # Based on the sum of their exon lengths
		
		# In case there are small differences between annotation data and external file, due to e.g. slightly different Ensembl versions
		gene.data <- gene.data[rownames(gene.counts),]
		total.gene.data <- gene.data # We need this for some total stats
	}
	else if (count.type=="gene")
	{
		if (annotation=="download")
		{
			disp("Downloading gene annotation for ",org,"...")
			gene.data <- get.annotation(org,count.type)
		}
		else if (annotation=="fixed")
		{
			disp("Reading stored gene annotation for ",org,"...")
			gene.data <- read.annotation(org,count.type)
		}
		else if (annotation=="embedded") # The following should work if annotation elements are arranged in MeV-like data style
		{
			if (!is.data.frame(counts))
			{
				disp("Reading counts file ",counts.name,"...")
				gene.counts <- read.delim(counts)
			}
			else
				gene.counts <- counts
			rownames(gene.counts) <- as.character(gene.counts[,id.col])
			all.cols <- 1:ncol(gene.counts)
			sam.cols <- match(unlist(sample.list),colnames(gene.counts))
			sam.cols <- sam.cols[which(!is.na(sam.cols))]
			ann.cols <- all.cols[-sam.cols]
			gene.data <- gene.counts[,ann.cols]
			gene.counts <- gene.counts[,sam.cols]
			colnames(gene.data)[id.col] <- "gene_id"
			if (!is.na(gc.col))
			{
				colnames(gene.data)[gc.col] <- "gc_content"
				if (max(gene.data$gc_content<=1)) # Is already divided
					gene.data$gc_content = 100*gene.data$gc_content
			}
			if (!is.na(name.col)) colnames(gene.data)[name.col] <- "gene_name"
			if (!is.na(bt.col)) colnames(gene.data)[bt.col] <- "biotype"
		}
		else # Reading from external file, similar to embedded
		{
			disp("Reading external gene annotation for ",org," from ",annotation,"...")
			gene.data <- read.delim(annotation)
			gene.data <- gene.data[rownames(gene.counts),]
			colnames(gene.data)[id.col] <- "gene_id"
			if (!is.na(gc.col))
			{
				colnames(gene.data)[gc.col] <- "gc_content"
				if (max(gene.data$gc_content<=1)) # Is already divided
					gene.data$gc_content = 100*gene.data$gc_content
			}
			if (!is.na(name.col)) colnames(gene.data)[name.col] <- "gene_name"
			if (!is.na(bt.col)) colnames(gene.data)[bt.col] <- "biotype"
		}
		total.gene.data <- gene.data # We need this for some total stats
		exon.filter.result <- NULL

		if (annotation!="embedded") # Else everything is provided and done
		{
			if (!is.null(counts)) # Otherwise it's coming ready from read2count
			{
				if (!is.data.frame(counts)) # Else it's already here
				{
					disp("Reading counts file ",counts.name,"...")
					gene.counts <- read.delim(counts)
				}
				else # Already a data frame as input
					gene.counts <- counts
				rownames(gene.counts) <- as.character(gene.counts[,id.col])
				gene.counts <- gene.counts[,unlist(sample.list,use.names=FALSE)]
			}
			else # Coming from read2count
			{
				if (from.raw) # Double check
				{
					r2c <- read2count(file.list,file.type,gene.data)
					gene.counts <- r2c$counts
					if (is.null(libsize.list))
						libsize.list <- r2c$libsize
				}
			}
		}

		gene.data <- gene.data[rownames(gene.counts),]
		gene.length <- gene.data$end - gene.data$start # Based on total gene lengths
	}

	# Transform GC-content and biotype
	gene.data$gc_content <- as.numeric(gene.data$gc_content)/100
	if (is.null(gene.data$biotype))
		gene.data$biotype <- rep(NA,nrow(gene.data))
	names(gene.length) <- rownames(gene.counts)
	attr(gene.data,"gene.length") <- gene.length

	# GC bias is NOT alleviated if we do not remove the zeros!!!
	disp("Removing genes with zero counts in all samples...")
	the.zeros <- which(apply(gene.counts,1,filter.low,0))
	if (length(the.zeros)>0)
	{
		gene.counts <- gene.counts[-the.zeros,]
		gene.data <- gene.data[-the.zeros,]
		attr(gene.data,"gene.length") <- gene.length[-the.zeros]
		# Store the filtered, maybe we do some stats
		gene.data.zero <- gene.data[the.zeros,]
		attr(gene.data.zero,"gene.length") <- gene.length[the.zeros]
	}

	disp("Normalizing with: ",normalization)
	switch(normalization,
		edaseq = {
			norm.genes <- normalize.edaseq(gene.counts,sample.list,norm.args,gene.data,output="matrix")
		},
		deseq = {
			norm.genes <- normalize.deseq(gene.counts,sample.list,norm.args,output="native")
		},
		edger = {
			norm.genes <- normalize.edger(gene.counts,sample.list,norm.args,output="native")
		},
		noiseq = {
			norm.genes <- normalize.noiseq(gene.counts,sample.list,norm.args,gene.data,log.offset,output="native")
		},
		nbpseq = {
			norm.genes <- normalize.nbpseq(gene.counts,sample.list,norm.args,libsize.list,output="native")
		},
		none = { # In case some external normalization is applied (e.g. equal read counts from all samples)
			norm.genes <- gene.counts
		}
	)
	
	switch(class(norm.genes),
		CountDataSet = { # Has been normalized with DESeq
			temp.matrix <- round(counts(norm.genes,normalized=TRUE))
		},
		DGEList = { # Has been normalized with edgeR
			if (norm.args$main.method=="classic")
				temp.matrix <- round(norm.genes$pseudo.counts)
			else if (norm.args$main.method=="glm") { # Trick found at http://cgrlucb.wikispaces.com/edgeR+spring2013
				scl <- norm.genes$samples$lib.size * norm.genes$samples$norm.factors
				temp.matrix <- round(t(t(norm.genes$counts)/scl)*mean(scl))
			}
		},
		matrix = { # Has been normalized with EDASeq or NOISeq or nothing
			temp.matrix <- norm.genes
		},
		list = { # Has been normalized with NBPSeq and main method was "nbpseq"
			temp.matrix <- as.matrix(round(sweep(norm.genes$counts,2,norm.genes$norm.factors,"*")))
		},
		nbp = { # Has been normalized with NBPSeq and main method was "nbsmyth"... Jesus...
			 temp.matrix <- as.matrix(round(norm.genes$pseudo.counts))
		}
	)

	# Implement gene filters after normalization
	if (!is.null(gene.filters))
		gene.filter.result <- filter.genes(temp.matrix,gene.data,gene.filters)
	else
		gene.filter.result <- NULL

	# Unify the filters and filter
	the.dead.genes <- list(
		gene.filter.result$expression$median,
		gene.filter.result$expression$mean,
		gene.filter.result$expression$quantile,
		gene.filter.result$expression$known,
		gene.filter.result$expression$custom
	)
	#gene.filter.result$expression <- Reduce("union",the.dead.genes)
	the.dead <- unique(unlist(c(gene.filter.result,exon.filter.result)))
	if (length(the.dead>0)) # All method specific object are row-index subsettable
	{
		switch(class(norm.genes),
			CountDataSet = {
				norm.genes.expr <- norm.genes[-the.dead,]
			},
			DGEList = { # edgeR bug???
				norm.genes.expr <- norm.genes[-the.dead,]
				norm.genes.expr$AveLogCPM <- norm.genes.expr$AveLogCPM[-the.dead]
			},
			matrix = { # Has been normalized with EDASeq or NOISeq
				norm.genes.expr <- norm.genes[-the.dead,]
			},
			list = { # Has been normalized with NBPSeq, main.method="nbpseq"
				norm.genes.expr <- norm.genes
				norm.genes.expr$counts <- as.matrix(norm.genes.expr$counts[-the.dead,])
				norm.genes.expr$rel.frequencies <- norm.genes.expr$rel.frequencies[-the.dead,]
				norm.genes.expr$tags <- as.matrix(norm.genes.expr$tags[-the.dead,])
			},
			nbp = {
				norm.genes.expr <- norm.genes
				norm.genes.expr$counts <- as.matrix(norm.genes.expr$counts[-the.dead,])
				norm.genes.expr$pseudo.counts <- as.matrix(norm.genes.expr$pseudo.counts[-the.dead,])
				norm.genes.expr$pseudo.lib.sizes <- colSums(as.matrix(norm.genes.expr$counts))*rep(1,dim(norm.genes.expr$counts)[2])
			}
		)
		gene.counts.expr <- gene.counts[rownames(norm.genes.expr),]
		gene.data.expr <- gene.data[-the.dead,]
		attr(gene.data.expr,"gene.length") <- attr(gene.data,"gene.length")[-the.dead]
		# Store the filtered, maybe we do some stats
		gene.data.dead <- gene.data[the.dead,]
		attr(gene.data.dead,"gene.length") <- attr(gene.data,"gene.length")[the.dead]
	}
	else
	{
		norm.genes.expr <- norm.genes
		gene.counts.expr <- gene.counts
		gene.data.expr <- gene.data
		gene.data.dead <- NULL
	}
	
	# Store the final filtered, maybe we do some stats
	gene.data.filtered <- rbind(gene.data.zero,gene.data.dead)
	if (nrow(gene.data.filtered)>0)
	{
		if (nrow(gene.data.zero)>0)
			attr(gene.data.filtered,"gene.length") <- c(attr(gene.data.zero,"gene.length"),attr(gene.data.dead,"gene.length"))
		else
			attr(gene.data.filtered,"gene.length") <- attr(gene.data.dead,"gene.length")
	}

	# There is a small case that no genes are left after filtering...
	if(any(dim(norm.genes.expr)==0))
		stop("No genes left after gene and/or exon filtering! Try again with no filtering or less strict filter rules...")

	# Run the statistical test, norm.genes is always a method-specific object, handled in the metaseqr.stat.R stat.* functions
	cp.list <- vector("list",length(contrast))
	names(cp.list) <- contrast
	contrast.list <- make.contrast.list(contrast,sample.list)
	for (n in names(cp.list))
	{
		cp.list[[n]] <- vector("list",length(statistics))
		names(cp.list[[n]]) <- statistics
	}
	for (alg in statistics)
	{
		disp("Running statistical tests with: ",alg)	
		switch(alg,
			deseq = {
				p.list <- stat.deseq(norm.genes.expr,sample.list,contrast.list,norm.args)
			},
			edger = {
				p.list <- stat.edger(norm.genes.expr,sample.list,contrast.list,stat.args[[alg]])
			},
			noiseq = {
				p.list <- stat.noiseq(norm.genes.expr,sample.list,contrast.list,stat.args[[alg]],norm.args,gene.data.expr,log.offset)
			},
			bayseq = {
				p.list <- stat.bayseq(norm.genes.expr,sample.list,contrast.list,stat.args[[alg]],norm.args,libsize.list)
			},
			limma = {
				p.list <- stat.limma(norm.genes.expr,sample.list,contrast.list,stat.args[[alg]])
			},
			nbpseq = {
				p.list <- stat.nbpseq(norm.genes.expr,sample.list,contrast.list,stat.args[[alg]],norm.args,libsize.list)
			}
		)
		for (n in names(p.list))
			cp.list[[n]][[alg]] <- p.list[[n]]
	}
	for (n in names(cp.list))
		cp.list[[n]] <- do.call("cbind",cp.list[[n]])

	# Create the adjusted p-value matrices (if needed)
	if ("adj.p.value" %in% export.what)
	{
		adj.cp.list <- wapply(multic,cp.list,function(x,a) return(apply(x,2,p.adjust,a)),adjust.method)
		for (n in names(cp.list))
		{
			noi <- grep("noiseq",colnames(cp.list[[n]]))
			if (length(noi)>0)
			{
				if (length(strsplit(n,"_vs_")[[1]])==2) # DESeq has not run in this case, FDR cannot be calculated
					adj.cp.list[[n]][,noi] <- rep(NA,nrow(cp.list[[n]]))
			}
		}
	}

	# Calculate meta-statistics, if more than one statistical algorithm has been used
	if (length(statistics)>1)
	{
		disp("Performing meta-analysis with ",meta.p)
		switch(meta.p,
			intersection = {
				sum.p.list <- wapply(multic,cp.list,function(x) return(apply(x,1,prod)))
			},
			union = {
				sum.p.list <- wapply(multic,cp.list,function(x) {
					unp <- apply(x,1,sum)
					unp[unp>1] <- 1
					return(unp)
				})
			},
			fisher = {
				sum.p.list <- wapply(multic,cp.list,function(x) {
					tmp <- fisher.method(x,p.corr="none",zero.sub=1e-32)
					return(tmp$p.value)
				})
			},
			perm = {
				sum.p.list <- wapply(multic,cp.list,function(x) {
					tmp <- fisher.method.perm(x,p.corr="none",zero.sub=1e-32)
					return(tmp$p.value)
				})
			},
			whitlock = {
				sum.p.list <- wapply(multic,cp.list,function(x) return(apply(x,1,combine.test,method="z.transform")))
			},
			none = { # A default value must be there to use with volcanos, we say the one of the first statistic in order of input
				sum.p.list <- wapply(multic,cp.list,function(x) return(x[,1]))
			}
		)
		# ...and the adjusted p-value if requested
	}
	else # We assign the p-values from the only statistic used to sum.p.list in order to use it for stat plots
		sum.p.list <- cp.list
	if ("adj.meta.p.value" %in% export.what) # Useless for one statistics but just for safety
		adj.sum.p.list <- wapply(multic,sum.p.list,function(x,a) return(p.adjust(x,a)),adjust.method)
	
	#disp("Running statistical tests with: ",statistics)
	#contrast.list <- make.contrast.list(contrast,sample.list)
	#switch(statistics,
	#	deseq = {
	#		p.list <- stat.deseq(norm.genes.expr,sample.list,contrast.list,norm.args)
	#	},
	#	edger = {
	#		p.list <- stat.edger(norm.genes.expr,sample.list,contrast.list,stat.args)
	#	},
	#	noiseq = {
	#		p.list <- stat.noiseq(norm.genes.expr,sample.list,contrast.list,stat.args,norm.args,gene.data.expr,log.offset)
	#	},
	#	bayseq = {
	#		p.list <- stat.bayseq(norm.genes.expr,sample.list,contrast.list,stat.args,norm.args,libsize.list)
	#	},
	#	limma = {
	#		p.list <- stat.limma(norm.genes.expr,sample.list,contrast.list,stat.args)
	#	}
	#)
	
	# At this point, all method-specific objects must become a matrices for exporting and plotting
	switch(class(norm.genes.expr),
		CountDataSet = { # Has been processed with DESeq
			norm.genes <- round(counts(norm.genes,normalized=TRUE))
			norm.genes.expr <- round(counts(norm.genes.expr,normalized=TRUE))
		},
		DGEList = { # Has been processed with edgeR
			if (norm.args$main.method=="classic") {
				norm.genes <- round(norm.genes$pseudo.counts)
				norm.genes.expr <- round(norm.genes.expr$pseudo.counts)
			}
			else if (norm.args$main.method=="glm") { # Trick found at http://cgrlucb.wikispaces.com/edgeR+spring2013
				scl.r <- norm.genes$samples$lib.size * norm.genes$samples$norm.factors
				norm.genes <- round(t(t(norm.genes$counts)/scl.r)*mean(scl.r))
				scl.n <- norm.genes.expr$samples$lib.size * norm.genes.expr$samples$norm.factors
				norm.genes.expr <- round(t(t(norm.genes.expr$counts)/scl.n)*mean(scl.n))
			}
		},
		list = {
			norm.genes <- as.matrix(round(sweep(norm.genes$counts,2,norm.genes$norm.factors,"*")))
			norm.genes.expr <- as.matrix(round(sweep(norm.genes.expr$counts,2,norm.genes$norm.factors,"*")))
		},
		nbp = {
			norm.genes <- as.matrix(round(norm.genes$pseudo.counts))
			norm.genes.expr <- as.matrix(round(norm.genes.expr$pseudo.counts))
		}
		# We don't need the matrix case
	)
	
	disp("Building output files...")
	counter <- 1
	if (out.list) out <- make.export.list(contrast) else out <- NULL
	if (report) html <- make.export.list(contrast) else html <- NULL
	if ("normalized" %in% export.values)
		norm.list <- make.transformation(norm.genes.expr,export.scale,log.offset)
	if ("raw" %in% export.values)
		raw.list <- make.transformation(gene.counts.expr,export.scale,log.offset)
	for (cnt in contrast)
	{
		disp("  Contrast: ",cnt)

		export <- data.frame(row.names=rownames(gene.data.expr))
		if (report) export.html <- as.matrix(export)
		the.names <- character(0)
		if ("annotation" %in% export.what)
		{
			disp("    binding annotation...")
			export <- cbind(export,gene.data.expr)
			if (report) export.html <- cbind(export.html,make.html.cells(gene.data.expr,type="text"))
			the.names <- c(the.names,colnames(gene.data.expr))
		}
		if ("p.value" %in% export.what)
		{
			disp("    binding p-values...")
			export <- cbind(export,cp.list[[cnt]])
			if (report) export.html <- cbind(export.html,make.html.cells(cp.list[[cnt]]))
			the.names <- c(the.names,paste("p-value_",colnames(cp.list[[cnt]]),sep=""))
		}
		if ("adj.p.value" %in% export.what)
		{
			disp("    binding FDRs...")
			#if (statistics=="noiseq" && length(strsplit(cnt,"_vs_")[[1]])>2) # DESeq has run instead, FDR can be calculated
			#	export <- cbind(export,wapply(multic,p.list[cnt],wp.adjust,adjust.method))
			#else if (statistics=="noiseq" && length(strsplit(cnt,"_vs_")[[1]])==2) # NOISeq does not return the classical p-value
			#	export <- cbind(export,rep(NA,nrow(export)))
			#else
			#	export <- cbind(export,wapply(multic,p.list[cnt],wp.adjust,adjust.method))
			export <- cbind(export,adj.cp.list[[cnt]])
			#if (report) export.html <- cbind(export.html,make.html.cells(export[,ncol(export)]))
			if (report) export.html <- cbind(export.html,make.html.cells(adj.cp.list[[cnt]]))
			the.names <- c(the.names,paste("FDR_",colnames(adj.cp.list[[cnt]]),sep=""))
		}
		if ("meta.p.value" %in% export.what && length(statistics)>1) # Otherwise it does not exist
		{
			disp("    binding meta p-values...")
			export <- cbind(export,sum.p.list[[cnt]])
			if (report) export.html <- cbind(export.html,make.html.cells(sum.p.list[[cnt]]))
			the.names <- c(the.names,paste("meta_p-value_",cnt,sep=""))
		}
		if ("adj.meta.p.value" %in% export.what && length(statistics)>1)
		{
			disp("    binding adjusted meta p-values...")
			export <- cbind(export,adj.sum.p.list[[cnt]])
			if (report) export.html <- cbind(export.html,make.html.cells(adj.sum.p.list[[cnt]]))
			the.names <- c(the.names,paste("meta_FDR_",cnt,sep=""))
		}
		if ("fold.change" %in% export.what)
		{
			if ("normalized" %in% export.values)
			{
				tmp <- make.fold.change(cnt,sample.list,norm.genes.expr,log.offset)
				if ("natural" %in% export.scale)
				{
					disp("    binding natural normalized fold changes...")
					export <- cbind(export,tmp)
					if (report) export.html <- cbind(export.html,make.html.cells(tmp))
					the.names <- c(the.names,paste("natural_normalized_fold_change_",colnames(tmp),sep=""))
				}
				if ("log2" %in% export.scale)
				{
					disp("    binding log2 normalized fold changes...")
					export <- cbind(export,log2(tmp))
					if (report) export.html <- cbind(export.html,make.html.cells(log2(tmp)))
					the.names <- c(the.names,paste("log2_normalized_fold_change_",colnames(tmp),sep=""))
				}
			}
			if ("raw" %in% export.values)
			{
				tmp <- make.fold.change(cnt,sample.list,gene.counts.expr,log.offset)
				if ("natural" %in% export.scale)
				{
					disp("    binding natural raw fold changes...")
					export <- cbind(export,tmp)
					if (report) export.html <- cbind(export.html,make.html.cells(tmp))
					the.names <- c(the.names,paste("natural_raw_fold_change_",colnames(tmp),sep=""))
				}
				if ("log2" %in% export.scale)
				{
					disp("    binding log2 raw fold changes...")
					export <- cbind(export,log2(tmp))
					if (report) export.html <- cbind(export.html,make.html.cells(log2(tmp)))
					the.names <- c(the.names,paste("log2_raw_fold_change_",colnames(tmp),sep=""))
				}
			}
		}
		if ("stats" %in% export.what)
		{
			conds <- strsplit(cnt,"_vs_")[[1]]
			for (cond in conds)
			{
				if ("normalized" %in% export.values)
				{
					if ("mean" %in% export.stats)
					{
						disp("    binding normalized mean counts...")
						tmp <- make.stat(sample.list[[cond]],norm.list,"mean",export.scale)
						export <- cbind(export,tmp)
						if (report) export.html <- cbind(export.html,make.html.cells(tmp))
						the.names <- c(the.names,paste(colnames(tmp),"_normalized_mean_counts_",cond,sep=""))
					}
					if ("median" %in% export.stats)
					{
						disp("    binding normalized median counts...")
						tmp <- make.stat(sample.list[[cond]],norm.list,"median",export.scale)
						export <- cbind(export,tmp)
						if (report) export.html <- cbind(export.html,make.html.cells(tmp))
						the.names <- c(the.names,paste(colnames(tmp),"_normalized_median_counts_",cond,sep=""))
					}
					if ("sd" %in% export.stats)
					{
						disp("    binding normalized count sds...")
						tmp <- make.stat(sample.list[[cond]],norm.list,"sd",export.scale)
						export <- cbind(export,tmp)
						if (report) export.html <- cbind(export.html,make.html.cells(tmp))
						the.names <- c(the.names,paste(colnames(tmp),"_normalized_sd_counts_",cond,sep=""))
					}
					if ("mad" %in% export.stats)
					{
						disp("    binding normalized count MADs...")
						tmp <- make.stat(sample.list[[cond]],norm.list,"mad",export.scale)
						export <- cbind(export,tmp)
						if (report) export.html <- cbind(export.html,make.html.cells(tmp))
						the.names <- c(the.names,paste(colnames(tmp),"_normalized_mad_counts_",cond,sep=""))
					}
					if ("cv" %in% export.stats)
					{
						disp("    binding normalized count CVs...")
						tmp <- make.stat(sample.list[[cond]],norm.list,"cv",export.scale)
						export <- cbind(export,tmp)
						if (report) export.html <- cbind(export.html,make.html.cells(tmp))
						the.names <- c(the.names,paste(colnames(tmp),"_normalized_cv_counts_",cond,sep=""))
					}
					if ("rcv" %in% export.stats)
					{
						disp("    binding normalized counts RCVs...")
						tmp <- make.stat(sample.list[[cond]],norm.list,"rcv",export.scale)
						export <- cbind(export,tmp)
						if (report) export.html <- cbind(export.html,make.html.cells(tmp))
						the.names <- c(the.names,paste(colnames(tmp),"_normalized_rcv_counts_",cond,sep=""))
					}
				}
				if ("raw" %in% export.values)
				{
					if ("mean" %in% export.stats)
					{
						disp("    binding raw mean counts...")
						tmp <- make.stat(sample.list[[cond]],raw.list,"mean",export.scale)
						export <- cbind(export,tmp)
						if (report) export.html <- cbind(export.html,make.html.cells(tmp))
						the.names <- c(the.names,paste(colnames(tmp),"_raw_mean_counts_",cond,sep=""))
					}
					if ("median" %in% export.stats)
					{
						disp("    binding raw median counts...")
						tmp <- make.stat(sample.list[[cond]],raw.list,"median",export.scale)
						export <- cbind(export,tmp)
						if (report) export.html <- cbind(export.html,make.html.cells(tmp))
						the.names <- c(the.names,paste(colnames(tmp),"_raw_median_counts_",cond,sep=""))
					}
					if ("sd" %in% export.stats)
					{
						disp("    binding raw counts sds...")
						tmp <- make.stat(sample.list[[cond]],raw.list,"sd",export.scale)
						export <- cbind(export,tmp)
						if (report) export.html <- cbind(export.html,make.html.cells(tmp))
						the.names <- c(the.names,paste(colnames(tmp),"_raw_sd_counts_",cond,sep=""))
					}
					if ("mad" %in% export.stats)
					{
						disp("    binding raw counts MADs...")
						tmp <- make.stat(sample.list[[cond]],raw.list,"mad",export.scale)
						export <- cbind(export,tmp)
						if (report) export.html <- cbind(export.html,make.html.cells(tmp))
						the.names <- c(the.names,paste(colnames(tmp),"_raw_mad_counts_",cond,sep=""))
					}
					if ("cv" %in% export.stats)
					{
						disp("    binding raw counts CVs...")
						tmp <- make.stat(sample.list[[cond]],raw.list,"cv",export.scale)
						export <- cbind(export,tmp)
						if (report) export.html <- cbind(export.html,make.html.cells(tmp))
						the.names <- c(the.names,paste(colnames(tmp),"_raw_cv_counts_",cond,sep=""))
					}
					if ("rcv" %in% export.stats)
					{
						disp("    binding raw counts RCVs...")
						tmp <- make.stat(sample.list[[cond]],raw.list,"rcv",export.scale)
						export <- cbind(export,tmp)
						if (report) export.html <- cbind(export.html,make.html.cells(tmp))
						the.names <- c(the.names,paste(colnames(tmp),"_raw_rcv_counts_",cond,sep=""))
					}
				}
			}
		}
		if ("counts" %in% export.what)
		{
			conds <- strsplit(cnt,"_vs_")[[1]]
			for (cond in conds)
			{
				if ("normalized" %in% export.values)
				{					
					disp("    binding all normalized counts for ",cond,"...")
					tmp <- make.matrix(sample.list[[cond]],norm.list,export.scale)
					export <- cbind(export,tmp)
					if (report) export.html <- cbind(export.html,make.html.cells(tmp))
					part.1 <- rep(paste(export.scale,"_normalized_counts_",sep=""),each=length(sample.list[[cond]]))
					part.2 <- paste(part.1,colnames(tmp),sep="")
					the.names <- c(the.names,part.2)
				}
				if ("raw" %in% export.values)
				{
					disp("    binding all raw counts for ",cond,"...")
					tmp <- make.matrix(sample.list[[cond]],raw.list,export.scale)
					export <- cbind(export,tmp)
					if (report) export.html <- cbind(export.html,make.html.cells(tmp))
					part.1 <- rep(paste(export.scale,"_raw_counts_",sep=""),each=length(sample.list[[cond]]))
					part.2 <- paste(part.1,colnames(tmp),sep="")
					the.names <- c(the.names,part.2)
				}
			}
		}
		names(export) <- the.names

		# Adjust the export based on what statistics have been done and a possible p-value cutoff
		if (!is.na(pcut))
		{
			if (length(statistics)>1)
			{
				switch(meta.p,
					intersection = {
						cut.ind <- which(apply(cp.list[[cnt]],1,function(x,p) return(all(x<p)),pcut)) # The correct way
					},
					union = {
						cut.ind <- which(apply(cp.list[[cnt]],1,function(x,p) return(any(x<p)),pcut)) # The correct way
					},
					fisher = {
						cut.ind <- which(sum.p.list[[cnt]]<pcut)
					},
					perm = {
						cut.ind <- which(sum.p.list[[cnt]]<pcut)
					},
					whitlock = {
						cut.ind <- which(sum.p.list[[cnt]]<pcut)
					},
					none = {
						cut.ind <- which(sum.p.list[[cnt]]<pcut)
					}
				)
				export <- export[cut.ind,]
				if (report) export.html <- export.html[cut.ind,]
				pp <- sum.p.list[[cnt]][cut.ind]
			}
			else
			{
				export <- export[sum.p.list[[cnt]]<pcut,]
				if (report) export.html <- export.html[sum.p.list[[cnt]]<pcut,]
				pp <- sum.p.list[[cnt]][sum.p.list[[cnt]]<pcut]
			}
		}
		else pp <- sum.p.list[[cnt]]
		export <- export[order(pp),]
		if (report) export.html <- export.html[order(pp),]
		# Final safety trigger
		na.ind <- grep("NA",rownames(export))
		if (length(na.ind)>0)
		{
			export <- export[-na.ind,]
			export.html <- export.html[-na.ind,]
		}
		res.file <- file.path(PROJECT.PATH[["lists"]],paste("metaseqr_out_",cnt,".txt.gz",sep=""))
		disp("    Writing output...")
		gzfh <- gzfile(res.file,"w")
		write.table(export,gzfh,quote=FALSE,row.names=FALSE,sep="\t")
		close(gzfh)
		if (out.list)
			out[[cnt]] <- export
		if (report)
		{
			the.html.header <- make.html.header(the.names)
			the.html.rows <- make.html.rows(export.html)
			the.html.body <- make.html.body(the.html.rows)
			the.html.table <- make.html.table(the.html.body,the.html.header,id=paste("table_",counter,sep=""))
			html[[cnt]] <- the.html.table
			counter <- counter+1
		}
	}

	if (!is.null(qc.plots))
	{
		disp("Creating quality control graphs...")
		plots <- list(
			raw=c("mds","biodetection","countsbio","saturation","readnoise"),
			norm=c("boxplot","gcbias","lengthbias","meandiff","meanvar","rnacomp"),
			stat=c("deheatmap","volcano","biodist"),
			other=c("filtered")
		)
		fig.raw <- fig.unorm <- fig.norm <- fig.stat <- fig.other <- vector("list",length(fig.format))
		names(fig.raw) <- names(fig.unorm) <- names(fig.norm) <- names(fig.stat) <- names(fig.other) <- fig.format
		for (fig in fig.format)
		{
			disp("Plotting in ",fig," format...")
			fig.raw[[fig]] <- diagplot.metaseqr(gene.counts,sample.list,annotation=gene.data,diagplot.type=intersect(qc.plots,plots$raw),is.norm=FALSE,output=fig,path=PROJECT.PATH$qc) # raw plots
			fig.unorm[[fig]] <- diagplot.metaseqr(gene.counts,sample.list,annotation=gene.data,diagplot.type=intersect(qc.plots,plots$norm),is.norm=FALSE,output=fig,path=PROJECT.PATH$normalization) # un-normalized plots
			fig.norm[[fig]] <- diagplot.metaseqr(norm.genes,sample.list,annotation=gene.data,diagplot.type=intersect(qc.plots,plots$norm),is.norm=TRUE,output=fig,path=PROJECT.PATH$normalization) # normalized plots
			fig.stat[[fig]] <- diagplot.metaseqr(norm.genes.expr,sample.list,annotation=gene.data.expr,contrast.list=contrast.list,p.list=sum.p.list,thresholds=list(p=pcut,f=1),
				diagplot.type=intersect(qc.plots,plots$stat),is.norm=TRUE,output=fig,path=PROJECT.PATH$statistics) # statistical plots
			if (!is.null(gene.data.filtered))
				fig.other[[fig]] <- diagplot.metaseqr(gene.data.filtered,sample.list,annotation=total.gene.data,diagplot.type=intersect(qc.plots,plots$other),is.norm=FALSE,output=fig,path=PROJECT.PATH$qc) # other plots
			else fig.other[[fig]] <- NULL
		}
	}

	if (report)
	{
		disp("Creating HTML report...")
		if (!is.null(qc.plots))
		{
			# First create zip archives of the figures
			disp("Compressing figures...")
			zipfiles <- file.path(PROJECT.PATH$plots,paste("metaseqr_figures_",fig.format,".zip",sep=""))
			names(zipfiles) <- fig.format
			for (f in fig.format)
			{
				files <- c(
					dir(PROJECT.PATH$qc,pattern=paste(".",f,sep=""),full.names=TRUE),
					dir(PROJECT.PATH$normalization,pattern=paste(".",f,sep=""),full.names=TRUE),
					dir(PROJECT.PATH$statistics,pattern=paste(".",f,sep=""),full.names=TRUE)
				)
				zip(zipfiles[f],files)
			}
			# Then create the final figure variables which brew will find...
			fig.raw <- fig.raw[["png"]]
			fig.unorm <- fig.unorm[["png"]]
			fig.norm <- fig.norm[["png"]]
			fig.stat <- fig.stat[["png"]]
			fig.other <- fig.other[["png"]]
		}

		if (tolower(report.template)=="default")
		{
			if (exists("TEMPLATE"))
				report.template=list(
					html=file.path(TEMPLATE,"metaseqr_report.html"),
					css=file.path(TEMPLATE,"styles.css"),
					logo=file.path(TEMPLATE,"logo.png")
				)
			else
				report.template=list(html=NULL,css=NULL,logo=NULL)
		}

		if (!is.null(report.template$html))
		{
			if (file.exists(report.template$html))
			{
				template <- report.template$html
				has.template <- TRUE
			}
			else
			{
				warning(paste("The template file",report.template$html,"was not found! The HTML report will NOT be generated."),
					call.=FALSE)
				has.template <- FALSE
			}
		}
		else
		{
			warning(paste("The report option was enabled but no template file is provided! The HTML report will NOT be generated."),
				call.=FALSE)
			has.template <- FALSE
		}
		if (!is.null(report.template$css))
		{
			if (file.exists(report.template$css))
				file.copy(from=report.template$css,to=PROJECT.PATH$main)
			else
				warning(paste("The stylesheet file",report.template$css,"was not found! The HTML report will NOT be styled."),
					call.=FALSE)
		}
		else
			warning(paste("The report stylesheet file was not provided! The HTML report will NOT be styled."),
				call.=FALSE)
		if (!is.null(report.template$logo))
		{
			if (file.exists(report.template$logo))
				file.copy(from=report.template$logo,to=PROJECT.PATH$main)
			else
				warning(paste("The report logo image",report.template$logo,"was not found!"),
					call.=FALSE)
		}
		else
			warning(paste("The report logo image was not provided!"),
				call.=FALSE)
		
		if (has.template)
		{
			TEMP <<- environment()
			brew(
				file=report.template$html,
				#output=file.path(PROJECT.PATH$main,paste(basename(PROJECT.PATH$main),"html",sep=".")),
				output=file.path(PROJECT.PATH$main,"index.html"),
				envir=TEMP
			)
		}
	}
	
	disp("\n",strftime(Sys.time()),": Data processing finished!\n\n")
	
	if (out.list) return(list(data=out,html=html))
} # End metaseqr

#' Assemble a gene model based on exon counts
#'
#' This function assembles gene models (single genes, not isoforms) based on the input exon read counts file and a gene annotation
#' data frame, either from the fixed annotations included with the package, or with the \code{\link{get.annotation}} function. The
#' gene.data argument should have a specific format and for this reason it's better to use one of the two aforementioned ways to
#' supply it. This function is intended mostly for internal use but can be used if the requirements are met.
#'
#' @param exon.counts the exon counts data frame produced by reading the exon read counts file.
#' @param sample.list the list containing condition names and the samples under each condition.
#' @param gene.data an annotation data frame from the same organism as exon.counts (such the ones produced by \code{get.annotation}).
#' @param multic a logical value indicating the presence of multiple cores. Defaults to FALSE. Do not change it if you are not sure
#' whether package multicore has been loaded or not.
#' @return A named list where names represent samples. Each list member is a also a named list where names correspond to gene ids
#' and members are named vectors. Each vector is named according to the exons corresponding to each gene and contains the read counts
#' for each exon. This structure is used for exon filtering and assembling final gene counts in the metaseqr pipeline.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' #Not yet available... 
construct.gene.model <- function(exon.counts,sample.list,gene.data,multic=FALSE) {
	the.counts <- vector("list",length(unlist(sample.list)))
	names(the.counts) <- unlist(sample.list,use.names=FALSE)
	the.genes <- as.character(unique(gene.data$gene_id))
	#the.exons <- as.character(unique(exon.counts$exon_id))
	#the.exons <- as.character(unique(exon.data[,id.col]))
	for (n in names(the.counts))
	{
		disp("  Separating exons per gene for ",n,"...")
		#the.counts[[n]] <- vector("list",length(the.genes))
		the.counts[[n]] <- the.genes
		names(the.counts[[n]]) <- the.genes
		the.counts[[n]] <- wapply(multic,the.counts[[n]],function(x,d,n) {
			tmp <- d[which(d$gene_id==x),c("start","end","exon_id",n)]
			xx <- tmp[,n]
			yy <- tmp$end - tmp$start
			names(xx) <- names(yy) <- tmp$exon_id
			return(list(count=xx,length=yy))
		},exon.counts,n)
	}
	return(the.counts)
}

#' Reduce the gene annotation in case of not all chromosomes present in counts
#'
#' This function reduces the gene annotation in case of exon reads and when the data to be analyzed do not contain all the standard
#' chromosomes of the genome under investigation. This can greatly reduce processing time in these cases.
#'
#' @param exon.data the exon annotation already reduced to the size of the input exon counts table.
#' @param gene.data an annotation data frame from the same organism as exon.counts (such the ones produced by \code{get.annotation}).
#' @return The \code{gene.data} annotation, reduced to have the same chromosomes as in \code{exon.data}, or the original \code{gene.data}
#' if \code{exon.data} do contain the standard chromosomes.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' #Not yet available...
reduce.gene.data <- function(exon.data,gene.data) {
	exon.chrs <- unique(as.character(exon.data$chromosome))
	gene.chrs <- unique(as.character(gene.data$chromosome))
	if (length(exon.chrs)!=length(gene.chrs)) {
		m <- match(gene.data$chromosome,exon.chrs)
		gene.data <- gene.data[which(!is.na(m)),]
	}
	return(gene.data)
}

#' Initialize environment
#'
#' Initializes metaseqr environmental variables. Internal use only.
#'
#' @author Panagiotis Moulos
init.envar <- function() {
	HOME <<- system.file(package="metaseqr")
	SCRIPT <<- file.path(HOME,"R")
	TEMPLATE <<- HOME
	ANNOTATION <<- file.path(HOME,"data")
}
#init.envar <- function() {
#	#HOME <<- "/media/HD4/Fleming/dev/metaseqr"
#	HOME <<- system.file(package="metaseqr")
#	SCRIPT <<- file.path(HOME,"R")
#	TEMPLATE <<- list(
#		HOME=file.path(HOME,"templates"),
#		HTML=file.path(HOME,"templates","html"),
#		LATEX=NULL,
#		CSS=file.path(HOME,"templates","css"),
#		IMAGE=file.path(HOME,"templates","images")
#	)
#	ANNOTATION <<- list(
#		HOME=file.path(HOME,"annotation"),
#		ENSEMBL=list(
#			GENE=file.path(HOME,"annotation","ensembl","gene"),
#			EXON=file.path(HOME,"annotation","ensembl","exon")
#		)	
#	)
#}
