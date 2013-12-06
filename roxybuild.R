library(roxyPackage)
library(knitr)

VERSION <- "0.9.1"

roxy.package(
	pck.source.dir="/media/HD4/Fleming/dev/metaseqR",
	pck.version=VERSION,
	#R.libs="/usr/local/lib/R/site-library",
	R.libs="/home/panos/R/x86_64-pc-linux-gnu-library/3.0",
	repo.root="/var/www/Rrepos",
	Rcmd.options=c(install="",build="",check="--as-cran",Rd2pdf="--pdf --no-preview"),
	pck.description=data.frame(
		Package="metaseqR",
		Type="Package",
		Title="An R package for the analysis and result reporting of RNA-Seq gene expression data, using multiple statistical algorithms",
		Author="Panagiotis Moulos <moulos@fleming.gr>",
		AuthorsR="c(person(given=\"Panagiotis\", family=\"Moulos\", email=\"moulos@fleming.gr\", role=c(\"aut\", \"cre\")))",
		Maintainer="Panagiotis Moulos <moulos@fleming.gr>",
		Depends="R (>= 2.13.0), Biobase, BiocGenerics, rjson, biomaRt, utils, knitr, EDASeq, DESeq, edgeR, limma, NOISeq, baySeq, NBPSeq, MADAM, survcomp, brew, gplots, corrplot, GenomicRanges, Rsamtools, rtracklayer, Repitools, qvalue, vsn, VennDiagram, log4r",
		Suggests="parallel, TCC",
		Description="Provides an interface to several normalization and statistical testing packages for RNA-Seq gene expression data. Additionally, it creates several diagnostic plots, performs meta-analysis by combinining the results of several statistical tests and reports the results in an interactive way.",
		License="GPL (>= 3)",
		Encoding="UTF-8",
		LazyLoad="yes",
		LazyData="yes",
		URL="http://www.fleming.gr",
		biocViews="Bioinformatics, Preprocessing, HighThroughputSequencing, RNAseq, DifferentialExpression",
		VignetteBuilder="knitr",
		stringsAsFactors=FALSE
	),
	actions=c(
		#"roxy",
		"cite",
		"html",
		"doc",
		#"log",
		"check",
		"package",
		"win",
		"macosx"
	)#,
	#ChangeLog=list(
	#	added=c("First release")
	#)
)
