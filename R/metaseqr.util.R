#' Fixed annotation updater
#'
#' A function to update the fixed annotations contained to avoid downloading every time if it's not embedded. It has no parameters.
#'
#' @return This function does not return anything. It updates the fixed annotation files instead.
#' @note This function cannot be used by users when the package is installed. For this reason it is not exported. If you want to
#' maintain a local copy of the package and update annotation at will, you can download the package source.
# @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' library(metaseqr)
#' annotations.update()
#'}
annotations.update <- function() {
	if(!require(biomaRt))
		stop("Bioconductor package biomaRt is required to update annotations!")
	VERBOSE <<- TRUE
	supported.types <- c("gene","exon")
	supported.orgs <- c("hg18","hg19","mm9","mm10","rn5","dm3","danRer7")
	if (exists("ANNOTATION")) {
		for (type in supported.types) {
			for (org in supported.orgs) {
				disp("Downloading and writing ",type,"s for ",org,"...")
				tryCatch({
					tmp <- get.annotation(org,type)
					var.name <- paste(org,type,sep=".")
					assign(var.name,tmp)
					#if (!file.exists(ANNOTATION$ENSEMBL[[toupper(type)]]))
					#	dir.create(ANNOTATION$ENSEMBL[[toupper(type)]],recursive=TRUE)
					#gzfh <- gzfile(file.path(ANNOTATION$ENSEMBL[[toupper(type)]],paste(org,".txt.gz",sep="")),"w")
					#write.table(tmp,gzfh,sep="\t",row.names=FALSE,quote=FALSE)
					#close(gzfh)},
					save(list=eval(parse(text="var.name")),file=file.path(ANNOTATION,paste(org,type,"rda",sep=".")),compress=TRUE)},
					error=function(e) {
						disp("!!! Probable problem with connection to Biomart...")
					},
					finally=""
				)
			}
		}
		disp("Finished!\n")
	}
	else
		stop("metaseqr environmental variables are not properly set up! Annotations cannot be updated...")
}

#' Fixed annotation reader
#'
#' A function to read fixed annotations from the local repository.
#'
#' @param org one of the supported organisms.
#' @param type "gene" or "exon".
#' @return A data frame with the \code{type} annotation for \code{org}.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' ann <- read.annotation("hg19","gene")
#'}
read.annotation <- function(org,type) {
	data(list=paste(org,type,sep="."))
	ann <- eval(parse(text=paste(org,type,sep=".")))
	if (type=="gene")
		rownames(ann) <- ann$gene_id
	else if (type=="exon")
		rownames(ann) <- ann$exon_id
	return(ann)
}

#' Default parameters for several metaseqr functions
#'
#' This function returns a list with the default settings for each filtering, statistical and normalization algorithm included in
#' the metaseqr package. See the documentation of the main function and the documentation of each statistical and normalization method
#' for details.
#'
#' @param what a keyword determining the procedure for which to fetch the default settings according to method parameter. It can be
#' one of "normalization", "statistics", "gene.filter", "exon.filter" or "biotype.filter".
#' @param method the supported algorithm included in metaseqr for which to fetch the default settings. When what is "normalization",
#' method is one of "edaseq", "deseq", "edger" or "noiseq". When what is "statistics", method is one of "deseq", "edger", "noiseq",
#' "bayseq", "limma" or "ebseq". When method is "biotype.filter", what is the input organism (see the main \code{\link{metaseqr}} 
#' help page for a list of supported organisms).
#' @return A list with default setting that can be used directly in the call of metaseqr.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' norm.args.edaseq <- get.defaults("normalization","edaseq")
#' stat.args.edger <- get.defaults("statistics","edger")
#'}
get.defaults <- function(what,method=NULL) {
	if (what %in% c("normalization","statistics") && is.null(method))
		stop("The method argument must be provided when what is \"normalization\" or \"statistics\"!")
	switch(what,
		normalization = {
			switch(method,
				edaseq = { return(list(within.which="loess",between.which="full")) },
				deseq = { return(list(method="pooled",sharingMode="fit-only",fitType="local")) },
				edger = {
					return(list(
						main.method="classic",									# classic or glm normalization method
						norm.method="TMM",refColumn=NULL,logratioTrim=0.3,
						sumTrim=0.05,doWeighting=TRUE,Acutoff=-1e10,
						p=0.75,rowsum.filter=5,prior.df=10,
						trend="movingave",span=NULL,							# classic estimateCommonDisp arguments
						tag.method="grid",grid.length=11,grid.range=c(-6,6),	# classic estimateTagwiseDisp arguments
						offset=NULL,glm.method="CoxReid",subset=10000,
						AveLogCPM=NULL,trend.method="auto",					# glm estimateGLMCommonDisp and estimateGLMTrendedDisp arguments
						dispersion=NULL										# glm estimateGLMTagwiseDisp arguments
					))
				},
				noiseq = {
					return(list(
						method="rpkm",								# which normalization
						long=1000,lc=1,k=0.5,						# common arguments
						refColumn=1,logratioTrim=0.3,sumTrim=0.05,
						doWeighting=TRUE,Acutoff=-1e+10			# TMM normalization arguments
					))
				}
			)
		},
		statistics = {
			switch(method,
				deseq = { return(list()) },
				edger = {
					return(list(
						main.method="classic",							# classic or glm fit
						dispersion=NULL,offset=NULL,weights=NULL,	# glmFit arguments
						lib.size=NULL,prior.count=0.125,start=NULL,
						method="auto",test="chisq",						# glmLRT arguments
						abundance.trend=TRUE,robust=FALSE,
						winsor.tail.p=c(0.05,0.1)						# glmLFTest arguments
					))
				},
				noiseq = {
					return(list(
						k=0.5,norm="n",replicates="biological",
						factor="class",conditions=NULL,pnr=0.2,
						nss=5,v=0.02,lc=1,						# noiseq general and specific arguments
						nclust=15,r=100,adj=1.5,
						a0per=0.9,filter=0,depth=NULL,		
						cv.cutoff=500,cpm=1						# noiseqbio specific arguments
						
					))
				},
				bayseq = {
					return(list(samplesize=10000,samplingSubset=NULL,equalDispersions=TRUE,estimation="QL",zeroML=FALSE,
						consensus=FALSE,moderate=TRUE,pET="BIC",marginalise=FALSE,subset=NULL,priorSubset=NULL,bootStraps=1,
						conv=1e-4,nullData=FALSE,returnAll=FALSE,returnPD=FALSE,discardSampling=FALSE,cl=NULL))
				},
				limma = { return(list()) }
			)
		},
		gene.filter = {
			return(list(
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
			))
		},
		exon.filter = {
			return(list(
				mnrpx=list(
					exons.per.gene=5,
					min.exons=2,
					frac=1/5
				)
			))
		},
		biotype.filter = {
			switch(method,
				hg18 = {
					return(list(
						unprocessed_pseudogene=TRUE,
						pseudogene=FALSE,
						miRNA=FALSE,
						retrotransposed=FALSE,
						protein_coding=FALSE,
						processed_pseudogene=FALSE,
						snRNA=FALSE,
						snRNA_pseudogene=TRUE,
						Mt_tRNA_pseudogene=TRUE,
						miRNA_pseudogene=TRUE,
						misc_RNA=FALSE,
						tRNA_pseudogene=TRUE,
						snoRNA=FALSE,
						scRNA_pseudogene=TRUE,
						rRNA_pseudogene=TRUE,
						snoRNA_pseudogene=TRUE,
						rRNA=TRUE,
						misc_RNA_pseudogene=TRUE,
						IG_V_gene=FALSE,
						IG_D_gene=FALSE,
						IG_J_gene=FALSE,
						IG_C_gene=FALSE,
						IG_pseudogene=TRUE,
						scRNA=FALSE
					))
				},
				hg19 = {
					return(list(
						pseudogene=FALSE,
						lincRNA=FALSE,
						protein_coding=FALSE,
						antisense=FALSE,
						processed_transcript=FALSE,
						snRNA=FALSE,
						sense_intronic=FALSE,
						miRNA=FALSE,
						misc_RNA=FALSE,
						snoRNA=FALSE,
						rRNA=TRUE,
						polymorphic_pseudogene=FALSE,
						sense_overlapping=FALSE,
						three_prime_overlapping_ncrna=FALSE,
						TR_V_gene=FALSE,
						TR_V_pseudogene=TRUE,
						TR_D_gene=FALSE,
						TR_J_gene=FALSE,
						TR_C_gene=FALSE,
						TR_J_pseudogene=TRUE,
						IG_C_gene=FALSE,
						IG_C_pseudogene=TRUE,
						IG_J_gene=FALSE,
						IG_J_pseudogene=TRUE,
						IG_D_gene=FALSE,
						IG_V_gene=FALSE,
						IG_V_pseudogene=TRUE
					))
				},
				mm9 = {
					return(list(
						pseudogene=FALSE,
						snRNA=FALSE,
						protein_coding=FALSE,
						antisense=FALSE,
						miRNA=FALSE,
						lincRNA=FALSE,
						snoRNA=FALSE,
						processed_transcript=FALSE,
						misc_RNA=FALSE,
						rRNA=TRUE,
						sense_overlapping=FALSE,
						sense_intronic=FALSE,
						polymorphic_pseudogene=FALSE,
						non_coding=FALSE,
						three_prime_overlapping_ncrna=FALSE,
						IG_C_gene=FALSE,
						IG_J_gene=FALSE,
						IG_D_gene=FALSE,
						IG_V_gene=FALSE,
						ncrna_host=FALSE
					))
				},
				mm10 = {
					return(list(
						pseudogene=FALSE,
						snRNA=FALSE,
						protein_coding=FALSE,
						antisense=FALSE,
						miRNA=FALSE,
						snoRNA=FALSE,
						lincRNA=FALSE,
						processed_transcript=FALSE,
						misc_RNA=FALSE,
						rRNA=TRUE,
						sense_intronic=FALSE,
						sense_overlapping=FALSE,
						polymorphic_pseudogene=FALSE,
						IG_C_gene=FALSE,
						IG_J_gene=FALSE,
						IG_D_gene=FALSE,
						IG_LV_gene=FALSE,
						IG_V_gene=FALSE,
						IG_V_pseudogene=TRUE,
						TR_V_gene=FALSE,
						TR_V_pseudogene=TRUE,
						three_prime_overlapping_ncrna=FALSE
					))
				},
				dm3 = {
					return(list(
						protein_coding=FALSE,
						ncRNA=FALSE,
						snoRNA=FALSE,
						pre_miRNA=FALSE,
						pseudogene=FALSE,
						snRNA=FALSE,
						tRNA=FALSE,
						rRNA=TRUE
					))
				},
				rn5 = {
					return(list(
						protein_coding=FALSE,
						pseudogene=FALSE,
						processed_pseudogene=FALSE,
						miRNA=FALSE,
						rRNA=TRUE,
						misc_RNA=FALSE
					))
				},
				danRer7 = {
					return(list(
						antisense=FALSE,
						protein_coding=FALSE,
						miRNA=FALSE,
						snoRNA=FALSE,
						rRNA=TRUE,
						lincRNA=FALSE,
						processed_transcript=FALSE,
						snRNA=FALSE,
						pseudogene=FALSE,
						sense_intronic=FALSE,
						misc_RNA=FALSE,
						polymorphic_pseudogene=FALSE,
						IG_V_pseudogene=TRUE,
						IG_C_pseudogene=TRUE,
						IG_J_pseudogene=TRUE,
						non_coding=FALSE,
						sense_overlapping=FALSE
					))
				}
			)
		}
	)
}

#' Annotation downloader
#'
#' This function connects to the EBI's Biomart service using the package biomaRt and downloads annotation elements (gene co-ordinates,
#' exon co-ordinates, gene identifications, biotypes etc.) for each of the supported organisms. See the help page of \code{\link{metaseqr}}
#' for a list of supported organisms. The function downloads annotation for an organism genes or exons.
#'
#' @param org the organism for which to download annotation.
#' @param type either "gene" or "exon".
#' @return A data frame with the canonical (not isoforms!) genes or exons of the requested organism. When type="genes", the data
#' frame has the following columns: chromosome, start, end, gene_id, gc_content, strand, gene_name, biotype. When type="exon" the
#' data frame has the following columns: chromosome, start, end, exon_id, gene_id, strand, gene_name, biotype. The gene_id and exon_id
#' correspond to Ensembl gene and exon accessions respectively. The gene_name corresponds to HUGO nomenclature gene names.
#' @note The data frame that is returned contains only "canonical" chromosomes for each organism. It does not contain haplotypes or
#' random locations and does not contain chromosome M.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' hg19.genes <- get.annotation("hg19","gene")
#' mm9.exons <- get.annotation("mm9","exon")
#'}
get.annotation <- function(org,type) {
	mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host=get.host(org),dataset=get.dataset(org))
	#mart <- useMart(biomart="ensembl",host=get.host(org),dataset=get.dataset(org))
	chrs.exp <- paste(get.valid.chrs(org),collapse="|")
	if (type=="gene") {
		bm <- getBM(attributes=get.gene.attributes(),mart=mart)
		ann <- data.frame(
			chromosome=paste("chr",bm$chromosome_name,sep=""),
			start=bm$start_position,
			end=bm$end_position,
			gene_id=bm$ensembl_gene_id,
			gc_content=bm$percentage_gc_content,
			strand=ifelse(bm$strand==1,"+","-"),
			gene_name=bm$external_gene_id,
			biotype=bm$gene_biotype
		)
		rownames(ann) <- ann$gene_id
	}
	else if (type=="exon") {
		bm <- getBM(attributes=get.exon.attributes(),mart=mart)
		ann <- data.frame(
			chromosome=paste("chr",bm$chromosome_name,sep=""),
			start=bm$exon_chrom_start,
			end=bm$exon_chrom_end,
			exon_id=bm$ensembl_exon_id,
			gene_id=bm$ensembl_gene_id,
			strand=ifelse(bm$strand==1,"+","-"),
			gene_name=bm$external_gene_id,
			biotype=bm$gene_biotype
		)
		rownames(ann) <- ann$exon_id
	}
	ann <- ann[order(ann$chromosome,ann$start),]
	ann <- ann[grep(chrs.exp,ann$chromosome),]
	ann$chromosome <- as.character(ann$chromosome)
	return(ann)
}

#' Biotype converter
#'
#' Returns biotypes as character vector. Internal use.
#'
#' @param a the annotation data frame (output of \code{\link{get.annotation}}).
#' @return A character vector of biotypes.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' hg19.bt <- get.biotypes(hg19.genes)
#'}
get.biotypes <- function(a) {
	return(as.character(unique(a$biotype)))
}

#' Annotation downloader helper
#'
#' Returns the appropriate Ensembl host address to get different versions of annotation from. Internal use.
#'
#' @param org the organism for which to return the host address.
#' @return A string with the host address.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' mm9.hist <- get.host("mm9")
#'}
get.host <- function(org) {
	switch(org,
		hg18 = { return("may2009.archive.ensembl.org") },
		hg19 = { return("www.ensembl.org") },
		mm9 = { return("may2012.archive.ensembl.org") },
		mm10 = { return("www.ensembl.org") },
		rn5 = { return("www.ensembl.org") },
		dm3 = { return("www.ensembl.org") },
		danRer7 = { return("www.ensembl.org") }
	)
}

#' Annotation downloader helper
#'
#' Returns a dataset (gene or exon) identifier for each organism recognized by the Biomart service for Ensembl. Internal use.
#'
#' @param org the organism for which to return the identifier.
#' @return A string with the dataset identifier.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' dm3.id <- get.dataset("dm3")
#'}
get.dataset <- function(org) {
	switch(org,
		hg18 = { return("hsapiens_gene_ensembl") },
		hg19 = { return("hsapiens_gene_ensembl") },
		mm9 = { return("mmusculus_gene_ensembl") },
		mm10 = { return("mmusculus_gene_ensembl") },
		rn5 = { return("rnorvegicus_gene_ensembl") },
		dm3 = { return("dmelanogaster_gene_ensembl") },
		danRer7 = { return("drerio_gene_ensembl") }
	)
}

#' Annotation downloader helper
#'
#' Returns a vector of chromosomes to maintain after annotation download. Internal use.
#'
#' @param org the organism for which to return the chromosomes. 
#' @return A character vector of chromosomes.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' hg18.chr <- get.valid.chrs("hg18")
#'}
get.valid.chrs <- function(org)
{
	switch(org,
		hg18 = {
			return(c(
				"chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
				"chr2","chr20","chr21","chr22","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
			))
		},
		hg19 = {
			return(c(
				"chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
				"chr2","chr20","chr21","chr22","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
			))
		},
		mm9 = {
			return(c(
				"chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
				"chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
			))
		},
		mm10 = {
			return(c(
				"chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
				"chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
			))
		},
		rn5 = {
			return(c(
				"chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
				"chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX"
			))
		},
		dm3 = {
			return(c(
				"chr2L","chr2LHet","chr2R","chr2RHet","chr3L","chr3LHet","chr3R","chr3RHet",
				"chr4","chrU","chrUextra","chrX","chrXHet","chrYHet"
			))
		},
		danRer7 = {
			return(c(
				"chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr2",
				"chr20","chr21","chr22","chr23","chr24","chr25","chr3","chr4","chr5","chr6","chr7","chr8","chr9"
			))
		}
	)
}

#' Annotation downloader helper
#'
#' Returns a vector of genomic annotation attributes which are used by the biomaRt package in order to fetch the gene annotation for
#' each organism. It has no parameters. Internal use.
#'
#' @return A character vector of Ensembl gene attributes.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' gene.attr <- get.gene.attributes()
#'}
get.gene.attributes <- function() {
	return(c(
		"chromosome_name",
		"start_position",
		"end_position",
		"ensembl_gene_id",
		"percentage_gc_content",
		"strand",
		"external_gene_id",
		"gene_biotype"
	))
}

#' Annotation downloader helper
#'
#' Returns a vector of genomic annotation attributes which are used by the biomaRt package in order to fetch the exon annotation for
#' each organism. It has no parameters. Internal use.
#'
#' @return A character vector of Ensembl exon attributes.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' exon.attr <- get.exon.attributes()
#'}
get.exon.attributes <- function() {
	return(c(
		"chromosome_name",
		"exon_chrom_start",
		"exon_chrom_end",
		"ensembl_exon_id",
		"strand",
		"ensembl_gene_id",
		"external_gene_id",
		"gene_biotype"
	))
}

#' Calculates fold changes
#'
#' Returns a matrix of fold changes based on the requested contrast, the list of all samples and the data matrix which is produced
#' by the metaseqr workflow. For details on the contrast, sample.list and log.offset parameters, see the main usage page of metaseqr.
#' This function is intended mostly for internal use but can also be used independently.
#'
#' @param contrast the vector of requested statistical comparison contrasts.
#' @param sample.list the list containing condition names and the samples under each condition.
#' @param data.matrix a matrix of gene expression data whose column names are the same as the sample names included in the sample
#' list.
#' @param log.offset a number to be added to each element of data matrix in order to avoid Infinity on log type data transformations.
#' @return A matrix of fold change ratios, treatment to control, as these are parsed from contrast.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' data.matrix <- round(1000*matrix(runif(400),100,4))
#' rownames(data.matrix) <- paste("gene_",1:100,sep="")
#' colnames(data.matrix) <- c("C1","C2","T1","T2")
#' fc <- make.fold.change("Control_vs_Treatment",list(Control=c("C1","C2"),
#'   Treatment=c("T1","T2")),data.matrix)
#'}
make.fold.change <- function(contrast,sample.list,data.matrix,log.offset=1) {
	conds <- strsplit(contrast,"_vs_")[[1]]
	fold.mat <- matrix(0,nrow(data.matrix),length(conds)-1)
	for (i in 2:length(conds)) { # First condition is ALWAYS reference
		samples.nom <- sample.list[[conds[i]]]
		samples.denom <- sample.list[[conds[1]]]
		nom <- data.matrix[,match(samples.nom,colnames(data.matrix))]
		denom <- data.matrix[,match(samples.denom,colnames(data.matrix))]
		if (!is.matrix(nom)) nom <- as.matrix(nom) # Cover the case with no replicates...
		if (!is.matrix(denom)) denom <- as.matrix(denom)
		mean.nom <- apply(nom,1,mean)
		mean.denom <- apply(denom,1,mean)
		#mean.nom <- ifelse(mean.nom==0,log.offset,mean.nom)
		if (any(mean.nom==0)) mean.nom <- mean.nom + log.offset
		#mean.denom <- ifelse(mean.denom==0,log.offset,mean.denom)
		if (any(mean.denom==0)) mean.denom <- mean.denom + log.offset
		fold.mat[,i-1] <- mean.nom/mean.denom
	}
	rownames(fold.mat) <- rownames(data.matrix)
	colnames(fold.mat) <- paste(conds[1],"_vs_",conds[2:length(conds)],sep="")
	return(fold.mat)
}

#' HTML report helper
#'
#' Returns a character matrix with html formatted table cells. Essentially, it converts the input data to text and places them in
#' a <td></td> tag set. Internal use.
#'
#' @param mat the data matrix (numeric or character)
#' @param type the type of data in the matrix ("numeric" or "character")
#' @param digits the number of digits on the right of the decimal points to pass to \code{\link{formatC}}. It has meaning when
#' type="numeric".
#' @return A character matrix with html formatted cells.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' data.matrix <- round(1000*matrix(runif(400),100,4))
#' rownames(data.matrix) <- paste("gene_",1:100,sep="")
#' colnames(data.matrix) <- c("C1","C2","T1","T2")
#' the.cells <- make.html.cells(data.matrix)
#'}
make.html.cells <- function(mat,type="numeric",digits=3) {
	if (type=="numeric")
		tmp <- formatC(mat,digits=digits,format="f")
	else
		tmp <- mat
	if (!is.matrix(tmp)) tmp <- as.matrix(tmp)
	tmp <- apply(tmp,c(1,2),function(x) paste("<td>",x,"</td>",sep=""))
	return(tmp)
}

#' HTML report helper
#'
#' Returns a character vector with html formatted rows. Essentially, it collapses every row of a matrix to a single character and
#' puts a <tr></tr> tag set around. It is meant to be applied to the output of \code{\link{make.html.cells}}. Internal use.
#'
#' @param mat the data matrix, usually the output of \code{\link{make.html.cells}} function.
#' @return A character vector with html formatted rows of a matrix.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' data.matrix <- round(1000*matrix(runif(400),100,4))
#' rownames(data.matrix) <- paste("gene_",1:100,sep="")
#' colnames(data.matrix) <- c("C1","C2","T1","T2")
#' the.cells <- make.html.cells(data.matrix)
#' the.rows <- make.html.rows(the.cells)
#'}
make.html.rows <- function(mat) {
	tmp <- apply(mat,1,paste,collapse="")
	tmp <- paste("<tr>",tmp,"</tr>",sep="")
	return(tmp)
}

#' HTML report helper
#'
#' Returns a character vector with an html formatted table head row. Essentially, it collapses the input row to a single character
#' and puts a <th></th> tag set around. It is meant to be applied to the output of \code{\link{make.html.cells}}. Internal use.
#'
#' @param h the colnames of a matrix or data frame, usually as output of \code{\link{make.html.cells}} function.
#' @return A character vector with html formatted header of a matrix.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' data.matrix <- round(1000*matrix(runif(400),100,4))
#' rownames(data.matrix) <- paste("gene_",1:100,sep="")
#' colnames(data.matrix) <- c("C1","C2","T1","T2")
#' the.cells <- make.html.cells(data.matrix)
#' the.header <- make.html.header(the.cells[1,])
#'}
make.html.header <- function(h) {
	tmp <- paste("<th>",h,"</th>",sep="")
	tmp <- paste(tmp,collapse="")
	tmp <- paste("<tr>",tmp,"</tr>",sep="")
	return(tmp)
}

#' HTML report helper
#'
#' Returns a character vector with an html formatted table. Essentially, it collapses the input rows to a single character and puts
#' a <tbody></tbody> tag set around. It is meant to be applied to the output of \code{\link{make.html.rows}}. Internal use.
#'
#' @param mat the character vector produced by \code{\link{make.html.rows}}.
#' @return A character vector with the body of mat formatted in html.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' data.matrix <- round(1000*matrix(runif(400),100,4))
#' rownames(data.matrix) <- paste("gene_",1:100,sep="")
#' colnames(data.matrix) <- c("C1","C2","T1","T2")
#' the.cells <- make.html.cells(data.matrix)
#' the.header <- make.html.header(the.cells[1,])
#' the.rows <- make.html.rows(the.cells)
#' the.body <- make.html.body(the.rows)
#'}
make.html.body <- function(mat) {
	tmp <- paste(mat,collapse="")
	return(tmp)
}

#' HTML report helper
#'
#' Returns a character vector with a fully html formatted table. Essentially, it binds the outputs of \code{\link{make.html.cells}},
#' \code{\link{make.html.rows}}, \code{\link{make.html.header}} and \code{\link{make.html.body}} to the final table and optionally
#' assigns an id attribute. The above functions are meant to format a data table so as it can be rendered by external tools such as
#' DataTables.js during a report creation. It is meant for internal use.
#'
#' @param b the table body as produced by \code{\link{make.html.body}}.
#' @param h the table header as produced by \code{\link{make.html.header}}.
#' @param id the table id attribute.
#' @return A fully formatted html table.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' data.matrix <- round(1000*matrix(runif(400),100,4))
#' rownames(data.matrix) <- paste("gene_",1:100,sep="")
#' colnames(data.matrix) <- c("C1","C2","T1","T2")
#' the.cells <- make.html.cells(data.matrix)
#' the.header <- make.html.header(the.cells[1,])
#' the.rows <- make.html.rows(the.cells)
#' the.body <- make.html.body(the.rows)
#' the.table <- make.html.table(the.body,the.header,id="my_table")
#'}
make.html.table <- function(b,h=NULL,id=NULL) {
	if (!is.null(id))
		html <- paste("<table id=\"",id,"\" class=\"datatable\">",sep="")
	else
		html <- "<table class=\"datatable\">"
	if (!is.null(h))
		html <- paste(html,"<thead>",h,"</thead>",sep="")
	html <- paste(html,"<tbody>",b,"</tbody></table>",sep="")
	return(html)
}	

#' Calculates several transformation of counts
#'
#' Returns a list of transformed (normalized) counts, based on the input count matrix data.matrix. The data transformations are passed
#' from the export.scale parameter and the output list is named accordingly. This function is intended mostly for internal use but can
#' also be used independently.
#'
#' @param data.matrix the raw or normalized counts matrix. Each column represents one input sample.
#' @param export.scale a character vector containing one of the supported data transformations ("natural", "log2","log10","vst").
#' See also the main help page of metaseqr.
#' @param log.offset a number to be added to each element of data.matrix in order to avoid Infinity on log type data transformations.
#' @return A named list whose names are the elements in export.scale. Each list member is the respective transformed data matrix.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' data.matrix <- round(1000*matrix(runif(400),100,4))
#' rownames(data.matrix) <- paste("gene_",1:100,sep="")
#' colnames(data.matrix) <- c("C1","C2","T1","T2")
#' tr <- make.transformation(data.matrix,c("log2","vst"))
#' head(tr$vst)
#'}
make.transformation <- function(data.matrix,export.scale,log.offset=1) {
	mat <- vector("list",length(export.scale))
	names(mat) <- export.scale
	if (!is.matrix(data.matrix)) data.matrix <- as.matrix(data.matrix)
	for (scl in export.scale) {
		switch(scl,
			natural = {
				mat[[scl]] <- data.matrix
			},
			log2 = {
				mat[[scl]] <- nat2log(data.matrix,base=2,off=log.offset)
			},
			log10 = {
				mat[[scl]] <- nat2log(data.matrix,base=10,off=log.offset)
			},
			vst = {
				fit <- vsn2(data.matrix,verbose=FALSE)
				mat[[scl]] <- predict(fit,newdata=data.matrix)
			}
		)
	}
	return(mat)
}

#' Calculates several statistices on read counts
#'
#' Returns a matrix of statistics calculated for a set of given samples. Internal use.
#'
#' @param samples a set of samples from the dataset under processing. They should match sample names from sample.list. See also the
#' main help page of \code{\link{metaseqr}}.
#' @param data.list a list containing natural or transformed data, typically an output from
#' \code{\link{make.transformation}}.
#' @param stat the statistics to calculate. Can be one or more of "mean", "median", "sd", "mad", "cv", "rcv". See also the main help
#' page of metaseqr.
#' @param export.scale the output transformations used as input also to \code{\link{make.transformation}}.
#' @return A matrix of statistics calculated based on the input sample names. The different data transformnations are appended columnwise.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' data.matrix <- round(1000*matrix(runif(400),100,4))
#' rownames(data.matrix) <- paste("gene_",1:100,sep="")
#' colnames(data.matrix) <- c("C1","C2","T1","T2")
#' tr <- make.transformation(data.matrix,c("log2","vst"))
#' st <- make.stat(c("C1","C2"),tr,c("mean","sd"),c("log2","vst"))
#'}
make.stat <- function(samples,data.list,stat,export.scale) {
	stat.result <- vector("list",length(export.scale))
	names(stat.result) <- export.scale
	for (scl in export.scale) {
		stat.data <- data.list[[scl]][,match(samples,colnames(data.list[[scl]]))]
		if (!is.matrix(stat.data)) stat.data <- as.matrix(stat.data)
		switch(stat,
			mean = {
				stat.result[[scl]] <- apply(stat.data,1,function(x,s) {
					if (s=="natural") return(round(mean(x))) else return(mean(x))
				},scl)
			},
			median = {
				stat.result[[scl]] <- apply(stat.data,1,function(x,s) {
					if (s=="natural") return(round(median(x))) else return(median(x))
				},scl)
			},
			sd = {
				stat.result[[scl]] <- apply(stat.data,1,function(x,s) {
					if (s=="natural") return(ceiling(sd(x))) else return(sd(x))
				},scl)
			},
			mad = {
				stat.result[[scl]] <- apply(stat.data,1,function(x,s) {
					if (s=="natural") return(ceiling(mad(x))) else return(mad(x))
				},scl)
			},
			cv = {
				stat.result[[scl]] <- apply(stat.data,1,function(x,s) {
					if (s=="natural") return(ceiling(sd(x))/round(mean(x))) else return(sd(x)/mean(x))
				},scl)
			},
			rcv = {
				stat.result[[scl]] <- apply(stat.data,1,function(x,s) {
					if (s=="natural") return(ceiling(mad(x))/round(median(x))) else return(mad(x)/median(x))
				},scl)
			}
		)
	}
	return(do.call("cbind",stat.result))
}

#' Results output build helper
#'
#' Returns a list of matrices based on the export scales that have been chosen from the main function and a subset of samples based
#' on the sample names provided in the sample.list argument of the main \code{\link{metaseqr}} function. Internal use.
#'
#' @param samples a set of samples from the dataset under processing. They should match sample names from sample.list. See also the
#' main help page of \code{\link{metaseqr}}.
#' @param data.list a list containing natural or transformed data, typically an output from \code{\link{make.transformation}}.
#' @param export.scale the output transformations used as input also to \code{\link{make.transformation}}.
#' @return A named list whose names are the elements in export.scale. Each list member is the respective sample subest data matrix.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' data.matrix <- round(1000*matrix(runif(400),100,4))
#' rownames(data.matrix) <- paste("gene_",1:100,sep="")
#' colnames(data.matrix) <- c("C1","C2","T1","T2")
#' tr <- make.transformation(data.matrix,c("log2","vst"))
#' mm <- make.matrix(c("C1","T1"),tr,"log2")
#' head(tr$vst)
#'}
make.matrix <- function(samples,data.list,export.scale="natural") {
	mat <- vector("list",length(export.scale))
	names(mat) <- export.scale
	for (scl in export.scale) {
		mat.data <- data.list[[scl]][,match(samples,colnames(data.list[[scl]]))]
		if (!is.matrix(mat.data)) mat.data <- as.matrix(mat.data)
		mat[[scl]] <- mat.data
	}
	return(do.call("cbind",mat))
}

#' Create contrast lists from contrast vectors
#'
#' Returns a list, properly structured to be used within the stat.* functions of the metaseqr package. See the main documentation for
#' the structure of this list and the example below. This function is mostly for internal use, as the stat.* functions can be supplied
#' directly with the contrasts vector which is one of the main \code{\link{metaseqr}} arguments.
#'
#' @param contrast a vector of contrasts in the form "ConditionA_vs_ConditionB" or
#' "ConditionA_vs_ConditionB_vs_ConditionC_vs_...".
#' In case of Control vs Treatment designs, the Control condition should ALWAYS be the first.
#' @param sample.list the list of samples in the experiment. See also the main help page of \code{\link{metaseqr}}.
#' @return A named list whose names are the contrasts and its members are named vectors, where the names are the sample names and the
#' actual vector members are the condition names. See the example.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' sample.list <- list(Control=c("C1","C2"),TreatmentA=c("TA1","TA2"),TreatmentB=c("TB1","TB2"))
#' contrast <- c("Control_vs_TreatmentA","Control_vs_TreatmentA_vs_TreatmentB"),
#' cl <- make.contrast.list(contrast,sample.list)
#' cl
#'}
make.contrast.list <- function(contrast,sample.list) {
	# Construction
	contrast.list <- vector("list",length(contrast))
	names(contrast.list) <- contrast
	# First break the contrast vector
	cnts <- strsplit(contrast,"_vs_")
	names(cnts) <- names(contrast.list)
	# Create list members
	for (n in names(contrast.list)) {
		contrast.list[[n]] <- vector("list",length(cnts[[n]]))
		for (i in 1:length(cnts[[n]])) {
			contrast.list[[n]][[i]] <- rep(cnts[[n]][i],length(sample.list[[cnts[[n]][i]]]))
			names(contrast.list[[n]][[i]]) <- sample.list[[cnts[[n]][[i]]]]
		}
	}
	return(contrast.list)
}

#' Creates sample list from file
#'
#' Create the main sample list from an external file.
#'
#' @param input a tab-delimited file structured as follows: the first line of the external tab delimited file should contain column 
#' names (names are not important). The first column MUST contain UNIQUE sample names and the second column MUST contain the biological
#' condition where each of the samples in the first column should belong to.
#' @return A named list whose names are the conditions of the experiments and its members are the samples belonging to each condition.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' targets <- data.frame(sample=c("C1","C2","T1","T2"),
#'   condition=c("Control","Control","Treatment","Treatment"))
#' write.table(targets,file="targets.txt",sep="\t",row.names=F,quote="")
#' sample.list <- make.sample.list("targets.txt")
#'}
make.sample.list <- function(input) {
	if (missing(input) || !file.exists(input))
		stop("File to make sample list from should be a valid existing text file!")
	tab <- read.delim(input)
	samples <- as.character(tab[,1])
	conditions <- unique(as.character(tab[,2]))
	if (length(samples) != length(unique(samples)))
		stop("Sample names must be unique for each sample!")
	sample.list <- vector("list",length(conditions))
	names(sample.list) <- conditions
	for (n in conditions)
		sample.list[[n]] <- samples[which(as.character(tab[,2]))==n]
	return(sample.list)
}

#' Project path constructor
#'
#' Create the main metaseqr project path. Internal use only.
#'
#' @param path The desired project path. Can be NULL for auto-generation.
#' @param f The input counts table file.
#' @return A list with project path elements.
#' @author Panagiotis Moulos
make.project.path <- function(path,f) {
	if (is.na(path) || is.null(path)) {
		if (!is.data.frame(f) && file.exists(f))
			main.path <- file.path(dirname(f),paste("metaseqr_result_",format(Sys.time(),format="%Y%m%d%H%M%S"),sep=""))
		else
			main.path <- file.path(getwd(),paste("metaseqr_result_",format(Sys.time(),format="%Y%m%d%H%M%S"),sep=""))
		project.path <- make.path.struct(main.path)
	}
	else {
		success <- tryCatch(
			if (!file.exists(path)) dir.create(path) else TRUE,
			error=function(e) {
				disp("Cannot create ",path,"! Is it a valid system path? Is there a write permissions problem? Reverting to automatic creation...")
				return(FALSE)
			},
			finally=""
		)
		if (success)
			project.path <- make.path.struct(path)
		else
			project.path <- make.project.path(NA,f)
	}
	return(project.path)
}

#' Project path constructor helper
#'
#' Helper for make.project.path. Internal use only.
#'
#' @param main.path The desired project path.
#' @return A named list whose names are the conditions of the experiments and its members are the samples belonging to each condition.
#' @author Panagiotis Moulos
make.path.struct <- function(main.path) {
	project.path <- list(
		main=main.path,
		lists=file.path(main.path,"lists"),
		plots=file.path(main.path,"plots"),
		qc=file.path(main.path,"plots","qc"),
		normalization=file.path(main.path,"plots","normalization"),
		statistics=file.path(main.path,"plots","statistics")
	)
	for (p in names(project.path))
		if (!file.exists(project.path[[p]]))
			dir.create(project.path[[p]],recursive=TRUE)
	return(project.path)
}

#' Intitialize output list
#'
#' Initializes metaseqr R output. Internal use only.
#'
#' @param con The contrasts.
#' @return An empty named list.
#' @author Panagiotis Moulos
make.export.list <- function(con) {
	f <- vector("list",length(con))
	names(f) <- con
	return(f)
}

#' Optimize rectangular grid plots
#'
#' Returns a vector for an optimized m x m plot grid to be used with e.g. par(mfrow). m x m is as close as possible to the input n.
#' Of course, there will be empty grid positions if n < m x m
#'
#' @param n An integer, denoting the total number of plots to be created.
#' @return A 2-element vector with the dimensions of the grid.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' g1 <- make.grid(16) # Returns c(4,4)
#' g2 <- make.grid(11) # Returns c(4,3)
#'}
make.grid <- function(n) {
	m <- 0
	while (n > m*m)
		m <- m+1
	if (n < m*m) {
		k <- m-1
		if (n > m*k)
			k <- k+1
		else {
			while (n > m*k)
				k=k-1
		}
	}
	else
		k <- m
	return(c(m,k))
}

#' Initializer of report messages
#'
#' Initializes metaseqr report tmeplate messages output. Internal use only.
#'
#' @param lang The language of the report. For now, only english ("en")
#' @return An named list with messages for each input option.
#' @author Panagiotis Moulos
make.report.messages <- function(lang) {
	switch(lang,
		en = {
			messages <- list(
				norm=list(
					edaseq="EDASeq",
					deseq="DESeq",
					edger="edgeR",
					noiseq="NOISeq"
				),
				stat=list(
					deseq="DESeq",
					edger="edgeR",
					noiseq="NOISeq",
					bayseq="baySeq",
					limma="limma"
				),
				meta=list(
					intersection="Intersection of individual resutls",
					union="Union of individual resutls",
					fisher="Fisher's method (R package MADAM)",
					perm="Fisher's method with permutations (R package MADAM)",
					whitlock="Whitlock's Z-transformation method (Bioconductor package survcomp)",
					none="No meta-analysis, p-values from the first supplied statistical algorithm"
				),
				adjust=list(
					holm="Holm FWER",
					hochberg="Hochberg DFR",
					hommel="Hommel FWER",
					bonferroni="Bonferroni FWER",
					bh="Benjamini-Hochberg FDR",
					by="Benjamini-Yekutiely FDR",
					fdr="Benjamini-Hochberg FDR",
					none="No multiple test correction",
					qvalue="Storey-Tibshirani FDR"
				),
				plots=list(
					mds="Multidimensional scaling",
					biodetection="Biotype detection",
					countsbio="Biotype counts",
					saturation="Sample and biotype saturation",
					rnacomp="RNA composition",
					boxplot="Boxplots",
					gcbias="GC-content bias",
					lengthbias="Transcript length bias",
					meandiff="Mean-difference plot",
					meanvar="Mean-variance plot",
					deheatmap="DEG heatmap",
					volcano="Volcano plot",
					biodist="DEG biotype detection",
					filtered="Filtered biotypes",
					none="No plot"
				),
				export=list(
					annotation="Annotation",
					p.value="p-value",
					adj.p.value="Adjusted p-value (FDR)",
					fold.change="Fold change",
					stats="Statistics",
					counts="Read counts",
					natural="Natural scale",
					log2="log2 scale",
					log10="log10 scale",
					vst="Variance stabilization transformation",
					raw="Raw values",
					normalized="Normalized values",
					mean="Mean",
					median="Median",
					sd="Standard deviation",
					mad="Median Absolute Deviation (MAD)",
					cv="Coefficient of Variation",
					rcv="Robust Coefficient of Variation"
				)
			)
		}
	)
	return(messages)
}

#' Interactive volcano plot helper
#'
#' Creates a list which contains the data series of a scatterplot, to be used for serialization with highcharts JavaScript plotting.
#' framework. Internal use only.
#'
#' @param x The x coordinates (should be a named vector!).
#' @param y The y coordinates.
#' @param a Alternative names for each point.
#' @return A list that is later serialized to JSON.
#' @author Panagiotis Moulos
make.highcharts.points <- function(x,y,a) {
	if (length(x)>0) {
		n <- names(x)
		x <- unname(x)
		y <- unname(-log10(y))
		stru <- vector("list",length(x))
		if (is.null(a)) {
			for (i in 1:length(x))
				stru[[i]] <- list(
					x=round(x[i],digits=3),
					y=round(y[i],digits=3),
					name=n[i]
				)
		}
		else {
			for (i in 1:length(x))
				stru[[i]] <- list(
					x=round(x[i],digits=3),
					y=round(y[i],digits=3),
					name=n[i],
					alt_name=a[i]
				)
		}
	}
	else
		stru <- list(x=NULL,y=NULL,name=NULL,alt_name=NULL)
	return(stru)
}

#' Create a class vector
#'
#' Creates a class vector from a sample list. Internal to the stat.* functions. Mostly internal use.
#'
#' @param sample.list the list containing condition names and the samples under each condition
#' @return A vector of condition names.
#' @author Panagiotis Moulos
#' @export
as.class.vector <- function(sample.list) {
	classes <- vector("list",length(sample.list))
	names(classes) <- names(sample.list)
	for (n in names(sample.list))
		classes[[n]] <- rep(n,times=length(sample.list[[n]]))
	classes <- unlist(classes,use.names=FALSE)
	names(classes) <- unlist(sample.list,use.names=FALSE)
	return(classes)
}

#' Argument getter
#'
#' Get argument(s) from a list of arguments, e.g. normalization arguments.
#'
#' @param arg.list the initial list of a method's (e.g. normalization) arguments. Can be created with the \code{\link{get.defaults}}
#' function.
#' @param arg.name the argument name inside the argument list to fetch its value.
#' @return The argument sub-list.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' norm.list <- get.defaults("normalization","egder")
#' a <- get.arg(norm.list,c("main.method","logratioTrim"))
#'}
get.arg <- function(arg.list,arg.name) {
	return(arg.list[arg.name])
}

#' Argument setter
#'
#' Set argument(s) to a list of arguments, e.g. normalization arguments.
#'
#' @param arg.list the initial list of a method's (e.g. normalization) arguments. Can be created with the \code{\link{get.defaults}}
#' function.
#' @param arg.name a named list with names the new arguments to be set, and mebers the values to be set or a vector of argument
#' names. In this case, arg.value must be supplied.
#' @param arg.value when arg.name is a vector of argument names, the values corresponding to these arguments.
#' @return the \code{arg.list} with the changed \code{arg.value} for \code{arg.name}.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' norm.list <- get.defaults("normalization","egder")
#' set.arg(norm.list,list(main.method="glm",logratioTrim=0.4))
#'}
set.arg <- function(arg.list,arg.name,arg.value=NULL) {
	if (is.list(arg.name))
		arg.list[names(arg.name)] <- arg.name
	else if (is.character(arg.name)) {
		tmp <- vector("list",length(arg.name))
		names(tmp) <- arg.name
		i <- 0
		for (n in arg.name) {
			i <- i + 1
			tmp[[n]] <- arg.value[i]
		}
		arg.list[arg.name] <- tmp
	}
	return(arg.list)
}

#' Multiple testing correction helper
#'
#' A wrapper around the \code{\link{p.adjust}} function to include also the qvalue adjustment procedure from the "qvalue" package.
#' Internal use.
#'
#' @param p a vector of p-values.
#' @param m the adjustment method. See the help of \code{\link{p.adjust}}.
#' @export
#' @author Panagiotis Moulos
wp.adjust <- function(p,m) {
	if (m=="qvalue")
		return(qvalue(p))
	else
		return(p.adjust(p,method=m))
}

#' List apply helper
#'
#' A wrapper around normal and parallel apply (\code{\link{mclapply}} or multicore package) to avoid excessive coding for control
#' of single or parallel code execution. Internal use.
#'
#' @param m a logical indicating whether to execute in parallel or not.
#' @param ... the rest arguments to lapply (or mclapply)
#' @export
#' @author Panagiotis Moulos
wapply <- function(m,...) {
	if (m)
		return(mclapply(...))
	else
		return(lapply(...))
}

#' Filtering helper
#'
#' Low score filtering function. Internal use.
#'
#' @param x a data numeric matrix.
#' @param f a threshold.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' filter.low(data,median(data))
#'}
filter.low <- function(x,f) { return(all(x<=f)) }

#' Filtering helper
#'
#' High score filtering function. Internal use.
#'
#' @param x a data numeric matrix.
#' @param f a threshold.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' filter.high(data,median(data))
#'}
filter.high <- function(x,f) { return(all(x>=f)) }

#' Message displayer
#'
#' Displays a message during execution of the several functions. Internal use.
#'
#' @param ... a vector of elements that compose the display message.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' i <- 1
#' disp("Now running iteration ",i,"...")
#'}
disp <- function(...) {
	if (VERBOSE) {
		cat("\n",...,sep="")
		flush.console()
	}
}

## An alternative to compressed text files. It proved that it requires a lot of space.
#update.annotation <- function() {
#	if (!require(RSQLite))
#		stop("R package RSQLite is required to update annotations!")
#	supported.types <- c("gene","exon")
#	supported.orgs <- c("hg18","hg19","mm9","mm10","rn5","dm3","danRer7")
#	if (exists("ANNOTATION")) {
#		db <- dbConnect(SQLite(),dbname=file.path(ANNOTATION$HOME,"annotation.sqlite"))
#		if(dbExistsTable(db,"gene")) dbRemoveTable(db,"gene")
#		if(dbExistsTable(db,"exon")) dbRemoveTable(db,"exon")
#		for (type in supported.types) {
#			for (org in supported.orgs) {
#				disp("Downloading and writing ",type,"s for ",org,"...")
#				tryCatch({
#					tmp <- get.annotation(org,type)
#					tmp <- cbind(tmp,rep(org,nrow(tmp)),rep("ensembl",nrow(tmp)))
#					colnames(tmp)[(ncol(tmp)-1):ncol(tmp)] <- c("organism","name")
#					z <- dbWriteTable(db,name=type,value=tmp,append=TRUE,row.names=FALSE)},
#					error=function(e) {
#						disp("!!! Probable problem with connection to Biomart...")
#					},
#					finally=""
#				)
#			}
#		}
#		dbSendQuery(db,"VACUUM")
#		dbGetQuery(db,"VACUUM")
#		dbDisconnect(db)
#		#dbSendQuery(db,"CREATE INDEX \"primary\" ON gene (gene_id ASC)")
#		#dbSendQuery(db,"CREATE INDEX \"primary\" ON exon (exon_id ASC)")
#		disp("Finished!\n")
#	}
#	else
#		stop("metaseqr environmental variables are not properly set up! Annotations cannot be updated...")
#}
#' Fixed annotation reader
#'
#' A function to read fixed annotations from an SQLite database. It proved that it requieres much more space than compressed files.
#'
#' @param org one of the supported organisms.
#' @param type "gene" or "exon".
#' @return A data frame with the \code{type} annotation for \code{org}.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' ann <- read.annotation("hg19","gene")
#'}
#read.annotation <- function(org,type) {
#	if (!require(RSQLite))
#		stop("R package RSQLite is required to read stored annotations!")
#	if (exists("ANNOTATION")) {
#		db <- dbConnect(SQLite(),dbname=file.path(ANNOTATION,"annotation.sqlite"))
#		if (type=="gene") {
#			query <- "SELECT chromosome, start, end, gene_id, gc_content, strand, gene_name, biotype FROM gene ORDER BY chromosome, start"
#			ann <- dbGetQuery(db,query)
#			rownames(ann) <- ann$gene_id
#		}
#		else if (type=="exon") {
#			query <- "SELECT chromosome, start, end, exon_id, gene_id, strand, gene_name, biotype FROM gene ORDER BY chromosome, start"
#			ann <- dbGetQuery(db,query)
#			rownames(ann) <- ann$exon_id
#		}
#		dbDisconnect(db)
#		return(ann)
#	}
#	else
#		stop("metaseqr environmental variables are not properly set up! Annotations cannot be accessed...")
#}
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
