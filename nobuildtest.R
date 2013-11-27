HOME <<- "/media/HD4/Fleming/dev/metaseqr"
SCRIPT <<- file.path(HOME,"R")
TEMPLATE <<- file.path(HOME,"inst")
ANNOTATION <<- file.path(HOME,"data")

#require(parallel)
require(log4r)
require(GenomicRanges)
require(rtracklayer)
require(Repitools)
require(Rsamtools)
require(EDASeq)
require(DESeq)
require(edgeR)
require(NOISeq)
require(NBPSeq)
require(baySeq)
require(limma)
require(qvalue)
require(MADAM)
require(survcomp)
require(vsn)
require(gplots)
require(corrplot)
require(VennDiagram)
require(brew)
require(biomaRt)
require(utils)
require(rjson)

source(file.path(SCRIPT,"metaseqr.main.R"))
source(file.path(SCRIPT,"metaseqr.argcheck.R"))
source(file.path(SCRIPT,"metaseqr.util.R"))
source(file.path(SCRIPT,"metaseqr.filter.R"))
source(file.path(SCRIPT,"metaseqr.norm.R"))
source(file.path(SCRIPT,"metaseqr.stat.R"))
source(file.path(SCRIPT,"metaseqr.meta.R"))
source(file.path(SCRIPT,"metaseqr.plot.R"))
source(file.path(SCRIPT,"metaseqr.export.R"))

## Redefine function based on the package data and path
#read.annotation <- function(org,type) {
#	if (exists("ANNOTATION")) {
#		load(file.path(ANNOTATION,paste(org,type,"rda",sep=".")))
#		ann <- eval(parse(text=paste(org,type,sep="."))) # Is it loaded?
#		if (type=="gene")
#			rownames(ann) <- ann$gene_id
#		else if (type=="exon")
#			rownames(ann) <- ann$exon_id
#		return(ann)
#	}
#	else
#		stop("metaseqr environmental variables are not properly set up! Annotations cannot be accessed...")
#}
