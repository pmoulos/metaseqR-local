\name{metaseqr}
\alias{metaseqr}
\alias{metaseqr.main}
\title{The main metaseqr pipeline}
\usage{
  metaseqr(counts, sample.list, contrast = NULL,
    libsize.list = NULL, id.col = 4, gc.col = NA,
    name.col = NA, bt.col = NA,
    annotation = c("download", "embedded", "fixed"),
    org = c("hg18", "hg19", "mm9", "mm10", "rno5", "dm3", "danRer7"),
    count.type = c("gene", "exon"),
    exon.filters = list(min.active.exons = list(exons.per.gene = 5, min.exons = 2, frac = 1/5)),
    gene.filters = list(length = list(length = 500), avg.reads = list(average.per.bp = 100, quantile = 0.75), expression = list(median = TRUE, mean = FALSE, quantile = NA, known = NA, custom = NA), biotype = get.defaults("biotype.filter", org[1])),
    normalization = c("edaseq", "deseq", "edger", "noiseq", "none"),
    norm.args = NULL,
    statistics = c("deseq", "edger", "noiseq", "bayseq", "limma"),
    stat.args = NULL,
    adjust.method = sort(c(p.adjust.methods, "qvalue")),
    meta.p = c("fisher", "perm", "sum", "whitlock", "intersection", "union", "none"),
    pcut = NA, log.offset = 1, preset = NULL,
    qc.plots = c("mds", "biodetection", "countsbio", "saturation", "rnacomp", "filtered", "boxplot", "gcbias", "lengthbias", "meandiff", "meanvar", "deheatmap", "volcano", "biodist"),
    fig.format = c("x11", "png", "jpg", "tiff", "bmp", "pdf", "ps"),
    out.list = FALSE, export.where = NA,
    export.what = c("annotation", "p.value", "adj.p.value", "meta.p.value", "adj.meta.p.value", "fold.change", "stats", "counts"),
    export.scale = c("natural", "log2", "log10", "vst"),
    export.values = c("raw", "normalized"),
    export.stats = c("mean", "median", "sd", "mad", "cv", "rcv"),
    restrict.cores = 0.6, report = TRUE,
    report.template = "default", verbose = TRUE, ...)
}
\arguments{
  \item{counts}{a text tab-delimited file containing gene
  or exon counts in one of the following formats: i) the
  first column contains unique gene or exon identifiers and
  the rest of the columns contain the read counts for each
  sample. Thus the first cell of each row is a gene or exon
  accession and the rest are integers representing the
  counts for that accession. In that case, the
  \code{annotation} parameter should strictly be "fixed",
  "download" or an external file in proper format. ii) The
  first n columns should contain gene or exon annotation
  elements like chromosomal locations, gene accessions,
  exon accessions, GC content etc. In that case, the
  \code{annotation} parameter can also be "embedded". The
  ideal embedded annotation contains 8 columns, chromosome,
  gene or exon start, gene or exon end, gene or exon
  accession, GC-content (fraction or percentage), strand,
  HUGO gene symbol and gene biotype (e.g. "protein_coding"
  or "ncRNA"). When the \code{annotation} parameter is
  "embedded", certain of these features are mandatory
  (co-ordinates and accessions). If they are not present,
  the pipeline will not run. If additional elements are not
  present (e.g. GC content or biotypes), certain features
  of metaseqr will not be available. For example, EDASeq
  normalization will not be performed based on a GC content
  covariate but based on gene length which is not what the
  authors of EDASeq suggest. If biotypes are not present, a
  lot of diagnostic plots will not be available. If the
  HUGO gene symbols are missing, the final annotation will
  contain only gene accessions and thus be less
  comprehensible. Generally, it's best to set the
  \code{annotation} parameter to "download" or "fixed" to
  ensure the most comprehensible results. Finally, counts
  can be a data frame satisfying the above conditions. It
  is a data frame by default when \code{bam2metaseqr} is
  used.}

  \item{sample.list}{a list containing condition names and
  the samples under each condition. It should have the
  format \code{sample.list <-
  list(ConditionA=c("Sample_A1", "Sample_A2", "Sample_A3"),
  ConditionB=c("Sample_B1", "Sample_B2"),
  ConditionC=c("Sample_C1", "Sample_C2"))}. The names of
  the samples in list members MUST match the column names
  containing the read counts in the counts file. If they do
  match, the pipeline will either crash or at best, ignore
  several of your samples. Alternative, \code{sample.list}
  can be a small tab-delimited file structured as follows:
  he first line of the external tab delimited file should
  contain column names (names are not important). The first
  column MUST contain UNIQUE sample names and the second
  column MUST contain the biological condition where each
  of the samples in the first column should belong to. In
  this case, the function \code{\link{make.sample.list}} is
  used.}

  \item{contrast}{a character vector of contrasts to be
  tested in the statistical testing step(s) of the metaseqr
  pipeline. Each element of the should STRICTLY have the
  format "ConditionA_vs_ConditionB_vs_...". A valid example
  based on the \code{sample.list} above is \code{contrast
  <- c("ConditionA_vs_ConditionB",
  "ConditionA_vs_ConditionC",
  "ConditionA_vs_ConditionB_vs_ConditionC")}. The first
  element of pairwise contrasts (e.g. "ConditionA" above)
  MUST be the control condition or any reference that
  ConditionB is checked against. metaseqr uses this
  convention to properly calculate fold changes. If it's
  NULL, a contrast between the first two members of the
  \code{sample.list} will be auto-generated.}

  \item{libsize.list}{an optional named list where names
  represent samples (MUST be the same as the samples in
  \code{sample.list}) and members are the library sizes
  (the sequencing depth) for each sample. For example
  \code{libsize.list <- list(Sample_A1=32456913,
  Sample_A2=4346818)}.}

  \item{id.col}{an integer denoting the column number in
  the file (or data frame) provided with the counts
  argument, where the unique gene or exon accessions are.
  Default to 4 which is the standard feature name column in
  a BED file.}

  \item{gc.col}{an integer denoting the column number in
  the file (or data frame) provided with the \code{counts}
  argument, where each gene's GC content is given. If not
  provided, GC content normalization provided by EDASeq
  will not be available.}

  \item{name.col}{an integer denoting the column number in
  the file (or data frame) provided with the counts
  argument, where the HUGO gene symbols are given. If not
  provided, it will not be available when reporting
  results. In addition, the "known" gene filter will not be
  available.}

  \item{bt.col}{an integer denoting the column number in
  the file (or data frame) provided with the counts
  argument, where the gene biotypes are given. If not
  provided, the "biodetection", "countsbio", "saturation",
  "filtered" and "biodist" plots will not be available.}

  \item{annotation}{instructs metaseqr where to find the
  annotation for the given counts file. It can be one of i)
  "download" (default) for automatic downloading of the
  annotation for the organism specified by the org
  parameter (using biomaRt), ii) "fixed" to retrieve the
  same annotation data as with "download" but from a fixed
  location inside the package (the "download" option always
  download the latest annotation specified by the org
  parameter), iii) "embedded" if the annotation elements
  are embedded in the read counts file or iv) a file
  specified by the user which should be as similar as
  possible to the "fixed" or "download" case, in terms of
  column structure.}

  \item{org}{the supported organisms by metaseqr. These can
  be, for human genomes "hg18" or "hg19", for mouse genomes
  "mm9", "mm10", for rat genomes "rno5", for drosophila
  genomes "dm3" and for zebrafish genomes "danRer7".}

  \item{count.type}{the type of reads inside the counts
  file. It can be one of "gene" or "exon". This is a very
  important and mandatory parameter as it defines the
  course of the workflow.}

  \item{exon.filters}{a named list whose names are the
  names of the supported exon filters and its members the
  filter parameters. See section "Exon filters" below for
  details.}

  \item{gene.filters}{a named list whose names are the
  names of the supported gene filters and its members the
  filter parameters. See section "Gene filters" below for
  details.}

  \item{normalization}{the normalization algorithm to be
  applied on the count data. It can be one of "edaseq"
  (default) for EDASeq normalization, "deseq" for the
  normalization algorithm (individual options specified by
  the \code{norm.args} argument) in the DESeq package,
  "edger" for the normalization algorithms present in the
  edgeR package (specified by the \code{norm.args}
  argument), "noiseq" for the normalization algorithms
  present in the NOISeq package (specified by the
  \code{norm.args} argument) or "none" to not normalize the
  data (highly unrecommended).}

  \item{norm.args}{a named list whose names are the names
  of the normalization algorithm parameters and its members
  parameter values. See section "Normalization parameters"
  below for details. Leave NULL for the defaults of
  \code{normalization}.}

  \item{statistics}{one or more statistical analyses to be
  performed by the metaseqr pipeline.It can be one or more
  of "deseq" (default) to conduct statistical test(s)
  implemented in the DESeq package, "edger" to conduct
  statistical test(s) implemented in the edgeR package,
  "limma" to conduct the RNA-Seq version of statistical
  test(s) implemented in the limma package, "noiseq" to
  conduct statistical test(s) implemented in the NOISeq
  package and "baySeq" to conduct statistical test(s)
  implemented in the baySeq package. In any case individual
  algorithm parameters are controlled by the contents of
  the \code{stat.args} list.}

  \item{stat.args}{a named list whose names are the names
  of the statistical algorithms used in the pipeline. Each
  member is another named list whose names are the
  algorithm parameters and its members are the parameter
  values. See section "Statistics parameters" below for
  details. Leave NULL for the defaults of
  \code{statistics}.}

  \item{adjust.method}{the multiple testing p-value
  adjustment method. It can be one of
  \code{\link{p.adjust.methods}} or "qvalue" from the
  qvalue Bioconductor package. Defaults to "BH" for
  Benjamini-Hochberg correction.}

  \item{meta.p}{the meta-analysis method to combine
  p-values from multiple statistical tests. It can be one
  of "fisher" (default), "perm", "whitlock",
  "intersection", "union" or "none". For the "fisher" and
  "perm" methods, see the documentation of the R package
  MADAM. For the "whitlock" method, see the documentation
  of the survcomp Bioconductor package. With the
  "intersection" option, the final p-value is the product
  of individual p-values derived from each method. However,
  the product is not used for the statistical cutoff to
  derive gene lists. In this case, the final gene list is
  derived from the common differentially expressed genes
  from all applied methods. Similarly, when meta.p is
  "union", the final list is derived from the union of
  individual methods and the final p-values are the sum of
  individual p-values. The latter can be used as a very
  lose statistical threshold to aggregate results from all
  methods regardless of their False Positive Rate.}

  \item{pcut}{a p-value cutoff for exporting differentially
  genes, default is to export all.}

  \item{log.offset}{an offset to be added to values during
  logarithmic transformations in order to avoid Infinity
  (default is 1).}

  \item{preset}{an analysis strictness preset. Not yet
  implemented but in the end it should be a vector like
  c("strict","loose", "verystrict","everything") etc.}

  \item{qc.plots}{a set of diagnostic plots to show/create.
  It can be one or more of "mds", "biodetection",
  "countsbio", "saturation", "rnacomp", "filtered",
  "boxplot", "gcbias", "lengthbias", "meandiff", "meanvar",
  "deheatmap", "volcano", "biodist". The "mds" stands for
  Mutlti-Dimensional Scaling and it creates a PCA-like plot
  but using the MDS dimensionality reduction instead. It
  has been succesfully used for NGS data (e.g. see the
  package htSeqTools) and it shows how well samples from
  the same condition cluster together. For "biodetection",
  "countsbio", "saturation", "rnacomp", "biodist" see the
  vignette of NOISeq package. The "saturation" case has
  been rewritten in order to display more samples in a more
  simple way. See the help page of
  \code{\link{plot.noiseq.saturation}}. For "gcbias",
  "lengthbias", "meandiff", "meanvar" see the vignette of
  EDASeq package. "lenghtbias" is similar to "gcbias" but
  using the gene length instead of the GC content as
  covariate. The "boxplot" option draws boxplots of log2
  transformed gene counts. The "filtered" option draws a
  4-panel figure with the filtered genes per chromosome and
  per biotype, as absolute numbers and as fractions of the
  genome. See also the help page of
  \code{\link{plot.filtered}}. The "deheatmap" option
  performs hierarchical clustering and draws a heatmap of
  differentially expressed genes. In the context of
  diagnostic plots, it's useful to see if samples from the
  same groups cluster together after statistical testing.
  The "volcano" option draws a volcano plot for each
  contrast and if a report is requested, an interactive
  volcano plot is presented in the HTML report. Set to
  \code{NULL} if you don't want any diagnostic plots
  created.}

  \item{fig.format}{the format of the output diagnostic
  plots. It can be one or more of "x11" (for direct
  display), "png", "jpg", "tiff", "bmp", "pdf", "ps".}

  \item{out.list}{a logical controlling whether to export a
  list with the results in the running environment.}

  \item{export.where}{an output directory for the project
  results (report, lists, diagnostic plots etc.)}

  \item{export.what}{the content of the final lists. It can
  be one or more of "annotation", to bind the annoation
  elements for each gene, "p.value", to bind the p-values
  of each method, "adj.p.value", to bind the multiple
  testing adjusted p-values, "meta.p.value", to bind the
  combined p-value from the meta-analysis,
  "adj.meta.p.value", to bind the corrected combined
  p-value from the meta-analysis, "fold.change", to bind
  the fold changes of each requested contrast, "stats", to
  bind several statistics calclulated on raw and normalized
  counts (see the \code{export.stats} argument), "counts",
  to bind the raw and normalized counts for each sample.}

  \item{export.scale}{export values from one or more
  transformations applied to the data. It can be one or
  more of "natural", "log2", "log10", "vst" (Variance
  Stabilizing Transormation, see the documentation of DESeq
  package).}

  \item{export.values}{export raw and normalized counts.}

  \item{export.stats}{calculate and export several
  statistics on raw and normalized counts, condition-wise.
  It can be one or more of "mean", "median", "sd", "mad",
  "cv" for the Coefficient of Variation, "rcv" for a robust
  version of CV where the median and the MAD are used
  instead of the mean and the standard deviation.}

  \item{restrict.cores}{in case of parallel execution of
  several subfunctions, the fraction of the available cores
  to use. In some cases if all available cores are used
  (\code{restrict.cores=1} and the system does not have
  sufficient RAM, the running machine might significantly
  slow down.}

  \item{report}{a logical value controlling whether to
  produce a summary report or not. Defaults to TRUE.}

  \item{report.template}{an HTML template to use for the
  report. Do not change this unless you know what you are
  doing.}

  \item{verbose}{print informative messages during
  execution? Defaults to TRUE.}

  \item{...}{further arguments that may be passed to
  plotting functions, related to \code{\link{par}}.}
}
\value{
  If \code{out.list} is \code{TRUE}, a named list whose
  length is the same as the number of requested contrasts.
  Each list member is named according to the corresponding
  contrast and contains a data frame of differentially
  expressed genes for that contrast. The contents of the
  data frame are defined by the \code{export.what,
  export.scale, export.stats, export.values} parameters. If
  \code{report} is \code{TRUE}, the output list contains
  two main elements. The first is described above (the
  analysis results) and the second contains the same
  results but in HTML formatted tables.
}
\description{
  This function is the main metaseqr workhorse and
  implements the main metaseqr workflow which performs data
  read, filtering, normalization and statistical selection,
  creates diagnostic plots and exports the results and a
  report if requested. The metaseqr function is responsible
  for assembling all the steps of the metaseqr pipeline
  which i) reads the input gene or exon read count table
  ii) performs prelimininary filtering of data by removing
  chrM and other non-essential information for a typical
  differential gene expression analysis as well as a
  preliminary expression filtering based on the exon
  counts, if an exon read count file is provided. iii)
  performs data normalization with one of currently widely
  used algorithms, including EDASeq (Risso et al., 2011),
  DESeq (Anders and Huber, 2010), edgeR (Robinson et al.,
  2010), NOISeq (Tarazona et al., 2012) or no normalization
  iv) performs a second stage of filtering based on the
  normalized gene expression according to several gene
  filters v) performs statistical testing with one or more
  of currently widely used algorithms, including DESeq
  (Anders and Huber, 2010), edgeR (Robinson et al., 2010),
  NOISeq (Tarazona et al., 2012), limma (Smyth et al.,
  2005) for RNA-Seq data, baySeq (Hardcastle et al., 2012)
  vi) in the case of multiple statistical testing
  algorithms, performs meta-analysis using one of five
  available methods (see the meta.p argument) vii) exports
  the resulting differentially expressed gene list in text
  tab-delimited format viii) creates a set of diagnostic
  plots either available in the aforementioned packages or
  metaseqr specific ones and ix) creates a comprehensive
  HTML report which summarizes the run information, the
  results and the diagnostic plots. Certain diagnostic
  plots (e.g. the volcano plot) can be interactive with the
  use of the external Highcharts
  (http://www.highcharts.com) JavaScript library for
  interactive graphs. Although the inputs to the metaseqr
  workflow are many, in practice, setting only very few of
  them and accepting the defaults as the rest can result in
  quite comprehensible results for mainstream organisms
  like mouse, human, fly and rat.
}
\note{
  Please note that currently only gene and exon annotation
  from Ensembl (http://www.ensembl.org) are supported.
  Thus, the unique gene or exon ids in the counts files
  should correspond to valid Ensembl gene or exon
  accessions for the organism of interest. If you are not
  sure about the source of your counts file or do not know
  how to produce it, it's better to start from the original
  BAM files and run the pipeline through the
  \code{\link{bam2metaseqr}} wrapper. Keep in mind that in
  this case the performance will be significantly lower and
  the overall running time significanlty higher as the R
  functions which are used to read BAM files to proper
  structures (GenomicRanges) and calculate the counts are
  quite slow. An alternative way is maybe the easyRNASeq
  package (Delhomme et al, 2012). The bam2metaseqr function
  does not use this package but rather makes use of
  standard Bioconductor functions to handle NGS data. If
  you wish to work outside R, you can work with other
  popular read counters such as the HTSeq read counter
  (http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html).
  Please also note that in the current version, the members
  of the gene.filters and exon.filters lists are not
  checked for validity so be careful to supply with correct
  names otherwise the pipeline will crash or at the best
  case scenario, will ignore the filters. Also note that
  when you are supplying metaseqr wtih an exon counts
  table, gene annotation is always downloaded or read from
  a fixed location inside the package. In addition to the
  above, if you have a multiple core system, be very
  careful on how you are using the \code{restrict.cores}
  argument and generally how many cores you are using with
  scripts purely written in R. The analysis with exon read
  data can very easily cause memory problems, so unless you
  have more than 64Gb of RAM available, consider setting
  restrict.cores to something like 0.2 when working with
  exon data.
}
\section{Exon filters}{
  The exon filters are a set of filters which are applied
  after the gene models are assembled from the read counts
  of individual exons and before the gene expression is
  summarized from the exons belonging to each gene. These
  filters can be applied when the input read counts file
  contains exon reads. It is not applicable when the input
  file already contains gene counts. Such filters can be
  for example "accept genes where all the exons contain
  more than x reads" or "accept genes where there is read
  presence in at least m/n exons, n being the total exons
  of the gene". Such filters are NOT meant for detecting
  differential splicing as also the whole metaseqr
  pipeline, this they should not be used in that context.
  The \code{exon.filters} argument is a named list of
  filters, where the names are the filter names and the
  members are the filter parameters (named lists with
  parameter name, parameter value). See the usage of the
  function for an example of how these lists are
  structured. The supported exon filters in the current
  version are: i) \code{min.active.exons} which implements
  a filter for demanding m out of n exons of a gene to have
  a certain read presence with parameters
  \code{exons.per.gene}, \code{min.exons} and \code{frac}.
  The filter is described as follows: if a gene has up to
  \code{exons.per.gene} exons, then read presence is
  required in at least \code{min.exons} of them, else read
  presence is required in a \code{frac} fraction of the
  total exons. With the default values, the filter
  instructs that if a gene has up to 5 exons, read presence
  is required in at least 2, else in at least 20% of the
  exons, in order to be accepted. More filters will be
  implemented in future versions and users are encouraged
  to propose exon filter ideas to the author by mail. See
  \code{metaseqr} usage for the defaults. Set
  exon.filters=NULL to not apply any exon filtering.
}

\section{Gene filters}{
  The gene filters are a set of filters applied to gene
  expression as this is manifested through the read present
  on each gene and are preferably applied after
  normalization. These filters can be applied both when the
  input file or data frame contains exon read counts and
  gene read counts. Such filter can be for example "accept
  all genes above a certain count threshold" or "accept all
  genes with expression above the median of the normalized
  counts distribution" or "accept all with length above a
  certain threshold in kb" or "exclude the 'pseudogene'
  biotype from further analysis". The supported gene
  filters in the current version, which have the same
  structure as the exon filters (named list of lists with
  filter names, parameter names and parameter arguments)
  are: i) \code{length} which implements a length filter
  where genes are accepted for further analysis if they are
  above \code{length} (its parameter) kb. ii)
  \code{avg.reads} which implements a filter where a gene
  is accepted for further analysis if it has more average
  reads than the \code{quantile} of the average count
  distribution per \code{average.per.bp} base pairs. In
  summary, the reads of each gene are averaged per
  \code{average.per.bp} based on each gene's length (in
  case of exons, input the "gene's length" is the sum of
  the lengths of exons) and the \code{quantile} quantile of
  the average counts distribution is calculated for each
  sample. Genes passing the filter should have an average
  read count larger than the maximum of the vector of the
  quantiles calculated above. iii) \code{expression} which
  implements a filter based on the overall expression of a
  gene. The parameters of this filter are: \code{median},
  where genes below the median of the overall count
  distribution are not accepted for further analysis (this
  filter has been used to distinguish between "expressed"
  and "not expressed" genes in several cases, e.g. (Mokry
  et al., 2011) with a logical as value, \code{mean} which
  is the same as \code{median} but using the mean,
  \code{quantile} which is the same as the previous two but
  using a specific quantile of the total counts
  distribution, \code{known}, where in this case, a set of
  known not-expressed genes in the system under
  investigation are used to estimate an expression cutoff.
  This can be quite useful, as the genes are filtered based
  on a "true biological" cutoff instead of a statistical
  cutoff. The value of this filter is a character vector of
  HUGO gene symbols (MUST be contained in the annotation,
  thus it's better to use \code{annotation= "fixed"} or
  \code{annotation="download"}) whose counts are used to
  build a "null" expression distribution. The 90th quantile
  of this distribution is then the expression cutoff. This
  filter can be combined with any other filter. Be careful
  with gene names as they are case sensitive and must match
  exactly ("Pten" is different than "PTEN"!). iv)
  \code{biotype} where in this case, genes with a certain
  biotype (MUST be contained in the annotation, thus it's
  better to use \code{annotation="fixed"} or
  \code{annotation= "download"}) are excluded from the
  analysis. This filter is a named list of logical, where
  names are the biotypes in each genome and values are
  \code{TRUE} or \code{FALSE}. If the biotype should be
  excluded, the value should be \code{TRUE} else
  \code{FALSE}. See the result of
  \code{get.defaults("biotype.filter","hg19")} for an
  example. Finally, in future versions there will be
  support for user-defined filters in the form of a
  function.
}

\section{Normalization parameters}{
  The normalization parameters are passed again as a named
  list where the names of the members are the normalization
  parameter names and the values are the normalization
  parameter values. You should check the documentation of
  the packages EDASeq, DESeq, edgeR and NOISeq for the
  parameter names and parameter values. There are two
  exceptions in parameter names: in case of
  \code{normalization="edaseq"} the only parameter names
  are \code{within.which} and \code{between.which},
  controlling the withing lane/sample and between
  lanes/samples normalization algorithm. In case of
  \code{normalization="edger"}, apart from the rest of the
  edgeR normalization arguments, there is the argument
  \code{main.method} which can be either "classic" or "glm"
  (see the edgeR's manual for details, briefly these are
  different algorithms for estimating sample dispersion
  parameters and normalization factors), and
  \code{norm.method} which controls the normalization
  method and replaces the \code{method} parameter
  respective edgeR's calls (again see edgeR's manual). For
  the rest of the algorithms, the parameter names are the
  same as the names used in the respective packages. For
  examples, please use the \code{\link{get.defaults}}
  function.
}

\section{Statistics parameters}{
  The statistics parameters as passed to statistical
  algorithms in metaseqr, exactly with the same way as the
  normalization parametes above. In this case, there is one
  more layer in list nesting. Thus, \code{stat.args} is a
  named list whose names are the names the algorithms used
  (see the \code{statistics} parameter). Each member is
  another named list,with parameters to be used for each
  statistical algorithm. Again, the names of the member
  lists are parameter names and the values of the member
  lists are parameter values. You should check the
  documentations of DESeq, edgeR, NOISeq, baySeq and limma
  for these parameters. There is one exception in parameter
  names: In case of \code{normalization="edger"}, apart
  from the rest of the edgeR statistical testing arguments,
  there is the argument \code{main.method} which can be
  either "classic" or "glm", again defining whether the
  binomial test or GLMs will be used for statistical
  testing. For examples, please use the
  \code{\link{get.defaults}} function.
}
\examples{
\dontrun{
require(metaseqr)

# An example pipeline with exon counts
data("hg18.exon.data",package="metaseqr")
metaseqr(
 counts=hg18.exon.counts,
 sample.list=list(CON=c("CON_BR1","CON_BR2"),DOX=c("DOX_BR1","DOX_BR2")),
 contrast=c("CON_vs_DOX"),
 libsize.list=list(CON_BR1=17041268,CON_BR2=23579904,DOX_BR1=16018639,DOX_BR2=26294259),
 id.col=4,
 annotation="download",
 org="hg18",
 count.type="exon",
 normalization="edaseq",
 statistics="deseq",
 pcut=0.05,
 qc.plots=c("mds", "biodetection", "countsbio", "saturation", "rnacomp", "boxplot", "gcbias", "lengthbias", "meandiff",
   "meanvar", "deheatmap", "volcano", "biodist", "filtered"),
 fig.format=c("png","pdf"),
 export.what=c("annotation","p.value","adj.p.value","fold.change","stats","counts"),
 export.scale=c("natural","log2","log10","vst"),
 export.values=c("raw","normalized"),
 export.stats=c("mean","median","sd","mad","cv","rcv"),
 restrict.cores=0.8,
 gene.filters=list(
   length=list(
     length=500
   ),
   avg.reads=list(
     average.per.bp=100,
     quantile=0.25
   ),
   expression=list(
     median=TRUE,
     mean=FALSE
   ),
   biotype=get.defaults("biotype.filter","hg18")
 )
)

# An example pipeline with gene counts
data("mm9.gene.data",package="metaseqr")
result <- metaseqr(
 counts=mm9.gene.counts,
 sample.list=list(e15.5=c("e15.5_1","e15.5_2"), P0.5=c("P0.5_1","P0.5_2"), P60=c("P60_1","P60_2")),
 contrast=c("e15.5_vs_P0.5","e15.5_vs_P60","P0.5_vs_P60","e15.5_vs_P0.5_vs_P60"),
 libsize.list=list(e15.5_1=46546134, e15.5_2=18439760, P0.5_1=21822789, P0.5_2=40813977,
   P60_1=35855191, P60_2=43611778),
 annotation="fixed",
 org="mm9",
 count.type="gene",
 normalization="edger",
 statistics=c("deseq","edger","noiseq"),
 meta.p="fisher",
 pcut=0.05,
 fig.format=c("png","pdf"),
 export.what=c("annotation","p.value","meta.p.value","adj.meta.p.value","fold.change"),
 export.scale=c("natural","log2"),
 export.values="normalized",
 export.stats=c("mean","sd","cv"),
 export.where=getwd(),
 restrict.cores=0.8,
 gene.filters=list(
   length=list(
     length=500
   ),
   avg.reads=list(
	    average.per.bp=100,
	    quantile=0.25
   ),
   expression=list(
	    median=TRUE,
	    mean=FALSE,
	    quantile=NA,
	    known=NA,
	    custom=NA
   ),
   biotype=get.defaults("biotype.filter","mm9")
 ),
 out.list=TRUE
)
head(result$data[["e15.5_vs_P0.5"]])
}
}
\author{
  Panagiotis Moulos
}

