#' Required packages validator
#'
#' Checks if all the required packages are present according to metaseqr input options. Internal use only.
#'
#' @param n normalization method
#' @param d statistics method
#' @param a multiple testing correction method
#' @param m meta-analysis method
#' @param s data transformation method
#' @param p qc plot types
#' @param r reporting option
#' @author Panagiotis Moulos
#' @export
check.packages <- function(n,d,a,m,s,p,r) {
	# Check normalization
	switch(n,
		edaseq = {
			if (!require(EDASeq)) stop("Bioconductor package EDASeq is required for \"edaseq\" normalization!")
		},
		deseq = {
			if (!require(DESeq)) stop("Bioconductor package DESeq is required for \"deseq\" normalization!")
		},
		edger = {
			if (!require(edgeR)) stop("Bioconductor package edgeR is required for \"edger\" normalization!")
		},
		noiseq = {
			if (!require(NOISeq)) stop("Bioconductor package NOISeq is required for \"noiseq\" normalization!")
		},
		nbpseq = {
			if (!require(NBPSeq)) stop("Bioconductor package NBPSeq is required for \"nbpseq\" normalization!")
		},
		none = {}
	)
	# Check differential expression method
	for (dd in d)
	{
		switch(dd,
			deseq = {
				if (!require(DESeq)) stop("Bioconductor package DESeq is required for \"deseq\" differential expression!")
			},
			edger = {
				if (!require(edgeR)) stop("Bioconductor package edgeR is required for \"edger\" differential expression!")
			},
			noiseq = {
				if (!require(NOISeq)) stop("Bioconductor package NOISeq is required for \"noiseq\" differential expression!")
			},
			bayseq = {
				if (!require(baySeq)) stop("Bioconductor package baySeq is required for \"bayseq\" differential expression!")
			},
			limma = {
				if (!require(limma)) stop("Bioconductor package limma is required for \"limma\" differential expression!")
				if (!require(edgeR)) stop("Bioconductor package edgeR is required for \"edger\" differential expression!")
			},
			nbpseq = {
				if (!require(NBPSeq)) stop("Bioconductor package NBPSeq is required for \"nbpseq\" differential expression!")
			},
			none = {}
		)
	}
	# Check multiple testing correction
	if (a=="qvalue" && !require(qvalue))
		stop("Bioconductor package qvalue is required for Storey-Tibshirani q-value correction!")
	# Check meta-analysis packages
	if (m %in% c("fisher","perm","sum") && !require(MADAM))
		stop("R package MADAM is required for \"fisher\", \"perm\" or \"sum\" p-value meta analysis!")
	if (m=="whitlock" && !require(survcomp))
		stop("Bioconductor package survcomp is required for \"whitlock\" p-value meta analysis!")
	# Check VST
	if (("vst" %in% s) && !require(vsn))
		stop("Bioconductor package vsn is required for \"vsn\" count data transformation!")
	# Check plots
	if (any(p %in% c("biodetection","countsbio","saturation","rnacomp","biodist")) && !require(NOISeq))
		stop("Bioconductor package NOISeq is required for some of the selected QC plots!")
	if (any(p %in% c("gcbias","lengthbias","meandiff","meanvar")) && !require(EDASeq))
		stop("Bioconductor package EDASeq is required for some of the selected QC plots!")
	if ("deheatmap" %in% p && !require(gplots))
		stop("R package gplots is required for some of the selected QC plots!")
	if (r && !require(brew))
		stop("R package brew is required to create an HTML report!")
	# Check biomaRt
	if (!require(biomaRt))
		stop("Bioconductor package biomaRt is required!")
	# Check utils
	if (!require(utils))
		stop("R package utils is required!")
	# Check rjson
	if (!require(rjson))
		stop("R package rjson is required!")
}

#' Text argument validator
#'
#' Checks if one or more given textual argument(s) is/are member(s) of a list of correct arguments. It's a more package-specific
#' function similar to \code{\link{match.arg}}. Mostly for internal use.
#' 
#' @param arg.name the name of the argument that is checked (for display purposes).
#' @param arg.value the value(s) of the argument to be checked.
#' @param arg.list a vector of valid argument values for arg.value to be matched against.
#' @param multiarg a logical scalar indicating whether arg.name accepts multiple arguments or not. In that case, all of the values
#' in arg.value are checked against arg.list.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' check.text.args("count.type",count.type,c("gene","exon"),multiarg=FALSE)
#' check.text.args("statistics",statistics,c("deseq","edger","noiseq","bayseq","limma"),
#'   multiarg=TRUE)
#'}
check.text.args <- function(arg.name,arg.value,arg.list,multiarg=FALSE) {
	if (multiarg) {
		arg.value <- tolower(arg.value)
		if (!all(arg.value %in% arg.list))
			stop("\"",arg.name,"\""," parameter must be one or more of ",paste(paste("\"",arg.list,sep=""),collapse="\", "),"\"!")
	}
	else {
		arg.save <- arg.value[1]
		arg.value <- tolower(arg.value[1])
		if (arg.name=="annotation") { # An exception must be added for annotation because it can be an external file too
			if (!(arg.value %in% arg.list) && !file.exists(arg.save))
				stop("\"",arg.name,"\""," parameter must be one of ",paste(paste("\"",arg.list,sep=""),collapse="\", "),"\" or an existing file!")
		}
		else {
			if (!(arg.value %in% arg.list))
				stop("\"",arg.name,"\""," parameter must be one of ",paste(paste("\"",arg.list,sep=""),collapse="\", "),"\"!")
		}
	}
}

#' Numeric argument validator
#'
#' Checks if one or more given numeric argument(s) satisfy several rules concerning numeric arguments, e.g. proper bounds or proper
#' format (e.g. it must be a number and not a character). Mostly for internal use.
#' 
#' @param arg.name the name of the argument that is checked (for display purposes).
#' @param arg.value the value(s) of the argument to be checked.
#' @param arg.type either the string "numeric" to denote generic double-like R numerics or "integer" for integer values.
#' @param arg.bounds a numeric or a vector with 2 elements, restraining arg.value to be within the bounds defined by the input vector
#' or e.g. larger (smaller) than the numeric value. See examples.
#' @param direction a string denoting to which direction the arg.value should be compared with arg.bounds. For example, "both" should
#' be given with a two element vector against which, arg.value will be checked to see whether it is smaller than the low boundary or
#' larger than the higher boundary. In that case, the function will throw an error. The direction parameter can be one of: "both" 
#' (described above), "botheq" (as above, but the arg.val is also checked for equality -closed intervals), "gt" or "gte" (check
#' whether arg.val is smaller or smaller than or equal to the first value of arg.bounds), "lt" or "lte" (check whether arg.val is
#' larger or larger than or equal to the first value of arg.bounds).
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' pcut <- 1.2 # A probability cannot be larger than 1! It will throw an error!
#' check.num.args("pcut",pcut,"numeric",c(0,1),"botheq")
#' pcut <- 0.05 # Pass
#' check.num.args("pcut",pcut,"numeric",c(0,1),"botheq")
#' gc.col <- 3.4 # A column in a file cannot be real! It will throw an error!
#' check.num.args("gc.col",gc.col,"integer",0,"gt")
#' gc.col <- 5 # Pass
#' check.num.args("gc.col",gc.col,"integer",0,"gt")
#'}
check.num.args <- function(arg.name,arg.value,arg.type,arg.bounds,direction) {
	switch(arg.type,
		numeric = {
			if (!is.numeric(arg.value))
				stop("\"",arg.name,"\""," parameter must be a numeric value!")
			if (!missing(arg.bounds)) {
				switch(direction,
					both = {
						if (arg.value<=arg.bounds[1] || arg.value>=arg.bounds[2])
							stop("\"",arg.name,"\""," parameter must be a numeric value larger than or equal to ",arg.bounds[1]," and smaller than or equal to ",arg.bounds[2],"!")
					},
					botheq = {
						if (arg.value<arg.bounds[1] || arg.value>arg.bounds[2])
							stop("\"",arg.name,"\""," parameter must be a numeric value larger than ",arg.bounds[1]," and smaller than ",arg.bounds[2],"!")
					},
					gt = {
						if (arg.value<=arg.bounds[1])
							stop("\"",arg.name,"\""," parameter must be a numeric value greater than ",arg.bounds[1],"!")
					},
					lt = {
						if (arg.value>=arg.bounds[1])
							stop("\"",arg.name,"\""," parameter must be a numeric value lower than ",arg.bounds[1],"!")
					},
					gte = {
						if (arg.value<arg.bounds[1])
							stop("\"",arg.name,"\""," parameter must be a numeric value greater than or equal to ",arg.bounds[1],"!")
					},
					lte = {
						if (arg.value>arg.bounds[1])
							stop("\"",arg.name,"\""," parameter must be a numeric value lower than or equal to ",arg.bounds[1],"!")
					}
				)
			}
		},
		integer = {
			if (!is.integer(arg.value))
				stop("\"",arg.name,"\""," parameter must be an integer!")
			if (!missing(arg.bounds)) {
				switch(direction,
					both = {
						if (arg.value<=arg.bounds[1] || arg.value>=arg.bounds[2])
							stop("\"",arg.name,"\""," parameter must be an integer larger than or equal to ",arg.bounds[1]," and smaller than or equal to ",arg.bounds[2],"!")
					},
					botheq = {
						if (arg.value<arg.bounds[1] || arg.value>arg.bounds[2])
							stop("\"",arg.name,"\""," parameter must be an integer larger than or equal to ",arg.bounds[1]," and smaller than or equal to ",arg.bounds[2],"!")
					},
					gt = {
						if (arg.value<=arg.bounds[1])
							stop("\"",arg.name,"\""," parameter must be an integer greater than ",arg.bounds[1],"!")
					},
					lt = {
						if (arg.value>=arg.bounds[1])
							stop("\"",arg.name,"\""," parameter must be an integer lower than ",arg.bounds[1],"!")
					},
					gte = {
						if (arg.value<arg.bounds[1])
							stop("\"",arg.name,"\""," parameter must be an integer greater than or equal to ",arg.bounds[1],"!")
					},
					lte = {
						if (arg.value>arg.bounds[1])
							stop("\"",arg.name,"\""," parameter must be an integer lower than or equal to ",arg.bounds[1],"!")
					}
				)
			}
		}
	)
}

#' File argument validator
#'
#' Checks if a file exists for specific arguments requiring a file input. Internal use only.
#'
#' @param arg.name argument name to display in a possible error.
#' @param arg.value the filename to check.
#' @author Panagiotis Moulos
#' @export
check.file.args <- function(arg.name,arg.value) {
	if (!file.exists(arg.value))
		stop("\"",arg.name,"\""," parameter must be an existing file!")
}

#' Parallel run validator
#'
#' Checks existence of multiple cores and loads multicore package. Internal use only.
#'
#' @param rc fraction of available cores to use.
#' @author Panagiotis Moulos
#' @export
check.parallel <- function(rc) {
	if (!require(multicore))
			multi <- FALSE
	else {
		multi <- TRUE
		ncores <- multicore:::detectCores()
		if (!missing(rc) || !is.na(rc) || !is.null(rc))
			ncores <- ceiling(rc*ncores)
		options(cores=ncores)
	}
	return(multi)
}

#' Contrast validator
#'
#' Checks if the contrast vector follows the specified format. Internal use only.
#'
#' @param cnt contrasts vector.
#' @param sample.list the input sample list.
#' @author Panagiotis Moulos
#' @export
check.contrast.format <- function(cnt,sample.list) {
	# This function will break cnt and check that all contrast counter parts are members of the names of the sample.list and 
	# contain the string "_vs_" as many times as the names of the sample.list minus 1. If satisfied return TRUE else error.
	cnts <- strsplit(cnt,"_vs_")
	#if (length(unique(unlist(cnts))) != length(names(sample.list)))
	if (!any(unique(unlist(cnts)) %in% names(sample.list)))
		stop("Condition names in sample list and contrast list do not match! Check if the contrasts follow the appropriate format (e.g. \"_vs_\" separating contrasting conditions...")
}

#' Library size validator
#'
#' Checks the names of the supplied library sizes. Internal use only.
#'
#' @param libsize.list the samples-names library size list.
#' @param sample.list the input sample list.
#' @author Panagiotis Moulos
#' @export
check.libsize <- function(libsize.list,sample.list) {
	if (length(intersect(names(libsize.list),unlist(sample.list,use.names=FALSE)))!=length(unlist(sample.list,use.names=FALSE))) {
		warning("Sample names in \"libsize.list\" and \"sample.list\" do not match! Library sized will be estimated from count data...",
			call.=FALSE)
		return(NULL)
	}
	else return(libsize.list)
}
