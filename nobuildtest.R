HOME <<- "/media/HD4/Fleming/dev/metaseqr"
SCRIPT <<- file.path(HOME,"R")
TEMPLATE <<- file.path(HOME,"inst")
ANNOTATION <<- file.path(HOME,"data")

source(file.path(SCRIPT,"metaseqr.main.R"))
source(file.path(SCRIPT,"metaseqr.argcheck.R"))
source(file.path(SCRIPT,"metaseqr.util.R"))
source(file.path(SCRIPT,"metaseqr.filter.R"))
source(file.path(SCRIPT,"metaseqr.norm.R"))
source(file.path(SCRIPT,"metaseqr.stat.R"))
source(file.path(SCRIPT,"metaseqr.plot.R"))

# Redefine function based on the package data and path
read.annotation <- function(org,type) {
	if (exists("ANNOTATION")) {
		load(file.path(ANNOTATION,paste(org,type,"rda",sep=".")))
		ann <- eval(parse(text=paste(org,type,sep="."))) # Is it loaded?
		if (type=="gene")
			rownames(ann) <- ann$gene_id
		else if (type=="exon")
			rownames(ann) <- ann$exon_id
		return(ann)
	}
	else
		stop("metaseqr environmental variables are not properly set up! Annotations cannot be accessed...")
}
