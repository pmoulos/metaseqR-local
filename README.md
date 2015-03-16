# metaseqr
An R package for the analysis, meta-analysis and result reporting of RNA-Seq gene expression data

## Summary
The metaseqr workflow performs data read, filtering, normalization and statistical selection, creates diagnostic plots and exports the results and a report if requested. Specifically, the metaseqr pipeline i) reads the input gene or exon, read count table ii) performs prelimininary filtering of data by removing chrM and other non-essential information for a typical differential gene expression analysis as well as a preliminary expression filtering based on the exon counts, if an exon read count file is provided. iii) performs data normalization with one of currently widely used algorithms, including EDASeq (Risso et al., 2011), DESeq (Anders and Huber, 2010), edgeR (Robinson et al., 2010), NOISeq (Tarazona et al., 2012) or no normalization iv) performs a second stage of filtering based on the normalized gene expression according to several gene filters v) performs statistical testing with one or more of currently widely used algorithms, including DESeq (Anders and Huber, 2010), edgeR (Robinson et al., 2010), NOISeq (Tarazona et al., 2012), limma (Smyth et al., 2005) for RNA-Seq data, baySeq (Hardcastle et al., 2012) vi) in the case of multiple statistical testing algorithms, performs meta-analysis using one of five available methods (see the meta.p argument) vii) exports the resulting differentially expressed gene list in text tab-delimited format viii) creates a set of diagnostic plots either available in the aforementioned packages or metaseqr specific ones and ix) creates a comprehensive HTML report which summarizes the run information, the results and the diagnostic plots. Certain diagnostic plots (e.g. the volcano plot) can be interactive with the use of the external Highcharts (http://www.highcharts.com) JavaScript library for interactive graphs. Although the inputs to the metaseqr workflow are many, in practice, setting only very few of them and accepting the defaults as the rest can result in quite comprehensible results for mainstream organisms like mouse, human, fly and rat.

## News
### 01 - 12 - 2014
metaseqR paper published in Nucleid Acids Research!
http://nar.oxfordjournals.org/cgi/content/abstract/gku1273?ijkey=FbJK4tp0FXcSlig&keytype=ref
### 01 - 03 - 2014
metaseqR is now a Bioconductor package! It is included in the latest 2.14 release.

## Installation
metaseqR's installation should take place through Bioconductor's `biocLite()` function in order to ensure the installation of related dependencies. Open an R session and type:
```R
source("http://bioconductor.org/biocLite.R")
biocLite("metaseqR") 
```
Everything else should run smoothly. For usage examples, have a look at the vignette and in the very extensive documentation within the package.
