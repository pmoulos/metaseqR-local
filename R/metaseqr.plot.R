#' Diagnostic plots for the metaseqr package
#'
#' This is the main function for producing sructured quality control and informative graphs base on the results of the various steps
#' of the metaseqr package. The graphs produced span a variety of issues like good sample reproducibility (Multi-Dimensional Scaling
#' plot, biotype detection, heatmaps. plot.metaseqr, apart from implementing certain package-specific plots, is a wrapper around several
#' diagnostic plots present in other RNA-Seq analysis packages such as EDASeq and NOISeq.
#'
#' @param object a matrix or a data frame containing count data derived before or after the normalization procedure, filtered or not
#' by the metaseqr's filters and/or p-value. The object can be fed to any of the plot.metaseqr plotting systems but not every plot is
#' meaningful. For example, it's meaningless to create a "biodist" plot for a count matrix before normalization or statistical testing.
#' @param sample.list the list containing condition names and the samples under each condition.
#' @param annotation a data.frame containing annotation elements for each row in object. Usually, a subset of the annotation obtained
#' by \code{\link{get.annotation}} or a subset of possibly embedded annotation with the input counts table. This parameter is optional
#' and required only when plot.type is any of "biodetection", "countsbio", "saturation","rnacomp", "biodist", "gcbias", "lengthbias" 
#' or "filtered".
#' @param contrast.list a named structured list of contrasts as returned by \code{\link{make.contrast.list}} or just the vector of
#' contrasts as defined in the main help page of \code{\link{metaseqr}}. This parameter is optional and required only when plot.type
#' is any of "deheatmap", "volcano" or "biodist".
#' @param p.list a list of p-values for each contrast as obtained from any of the stat.* methods of the metaseqr package. This parameter
#' is optional and required only when plot.type is any of "deheatmap", "volcano" or "biodist".
#' @param thresholds a list with the elements "p" and "f" which are the p-value and the fold change cutoff when plot.type="volcano".
#' @param plot.type one or more of the diagnostic plots supported in metaseqr package. Many of these plots require the presence of
#' additional package, somethng that is checked while running the main metaseqr function. The supported plots are "mds", "biodetection",
#' "countsbio", "saturation", "rnacomp", "boxplot", "gcbias", "lengthbias", "meandiff", "meanvar", "deheatmap", "volcano", "biodist",
#' "filtered".
#' For a brief description of these plots please see the main \code{\link{metaseqr}} help page.
#' @param is.norm a logical indicating whether object contains raw or normalized data. It is not essential and it serves only plot
#' annotation purposes.
#' @param output one or more R plotting device to direct the plot result to. Supported mechanisms: "x11" (default), "png", "jpg", 
#' "bmp", "pdf", "ps" or "json". The latter is currently available for the creation of interactive volcano plots only when reporting
#' the output, through the highcharts javascript library.
#' @param path the path to create output files.
#' @param ... further arguments to be passed to plot devices, such as parameter from \code{\link{par}}.
#' @return A named list containing the file names of the produced plots. Each list member is names according to the selected plotting
#' device and is also a named list, whose names are the plot types. The final contents are the file names in case the plots are written
#' to a physical location (not meaningful for "x11").
#' @note In order to make the best out of this function, you should generally provide the annotation argument as most and also the
#' most informative plots depend on this. If you don't know what is inside your counts table or how many annotation elements you can
#' provide by embedding it, it's always best to set the annotation parameter of the main metaseqr function to "download" or "fixed"
#' to use predefined annotations that work better with the functions of the whole package.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2",B3"))
#' plot.metaseqr(data.matrix,sample.list,plot.type=c("mds","boxplot"))
#'
#' norm.args <- get.defaults("normalization","deseq")
#' object <- normalize.deseq(data.matrix,sample.list,norm.args)
#' plot.metaseqr(object,sample.list,plot.type="boxplot")
#'
#' p <- stat.deseq(object)
#' plot.metaseqr(object,sample.list,contrast.list=contrast,p.list=p,plot.type="volcano")
#'}
plot.metaseqr <- function(
	object,sample.list,annotation=NULL,contrast.list=NULL,p.list=NULL,thresholds=list(p=0.05,f=1),
	plot.type=c("mds","biodetection","countsbio","saturation","rnacomp","boxplot","gcbias","lengthbias",
		"meandiff","meanvar","deheatmap","volcano","biodist","filtered"),
	is.norm=FALSE,output="x11",path=NULL,...
) {
	# annotation should have the format internally created here... This function can be used outside so it must be checked at some point...
	if (!is.matrix(object) && !is.data.frame(object))
		stop("object argument must be a matrix or data frame!")
	if (is.null(annotation) && any(plot.type %in% c("biodetection","countsbio","saturation","rnacomp","biodist","gcbias","lengthbias","filtered")))
		stop("annotation argument is needed when plot.type is \"biodetection\",\"countsbio\",\"saturation\",\"rnacomp\", \"biodist\", \"gcbias\", \"lengthbias\" or \"filtered\"!")
	if (any(plot.type %in% c("deheatmap","volcano","biodist"))) {
		if (is.null(contrast.list))
			stop("contrast.list argument is needed when plot.type is \"deheatmap\",\"volcano\" or \"biodist\"!")
		if (is.null(p.list))
			stop("The p argument which is a list of p-values for each contrast is needed when plot.type is \"deheatmap\",\"volcano\" or \"biodist\"!")
	}
	if (is.null(path)) path <- getwd()
	if (is.data.frame(object) && !("filtered" %in% plot.type)) object <- as.matrix(object)
	if (any(plot.type %in% c("biodetection","countsbio","saturation","rnacomp","biodist")))
		covars <- list(
			data=object,
			length=annotation$end - annotation$start,
			gc=annotation$gc_content,
			chromosome=annotation[,1:3],
			factors=data.frame(class=as.class.vector(sample.list)),
			biotype=annotation$biotype
		)

	raw.plots <- c("mds","biodetection","countsbio","saturation","rnacomp")
	norm.plots <- c("boxplot","gcbias","lengthbias","meandiff","meanvar")
	stat.plots <- c("deheatmap","volcano","biodist")
	other.plots <- c("filtered")
	files <- list()

	for (p in plot.type) {
		disp("  Plotting ",p,"...")
		if (p %in% raw.plots && !is.norm) {
			switch(p,
				mds = {
					files$mds <- plot.mds(object,sample.list,output=output,path=path)
				},
				biodetection = {
					files$biodetection <- plot.noiseq(object,sample.list,covars,which.plot=p,output=output,path=path,...)
				},
				countsbio = {
					files$countsbio <- plot.noiseq(object,sample.list,covars,which.plot=p,output=output,path=path,...)
				},
				saturation = {
					fil <- plot.noiseq(object,sample.list,covars,which.plot=p,output=output,path=path,...)
					files$saturation$biotype <- fil[["biotype"]]
					files$saturation$sample <- fil[["sample"]]
				},
				rnacomp = {
					files$rnacomp <- plot.noiseq(object,sample.list,covars,which.plot=p,output=output,path=path,...)
				}
			)
		}
		if (p %in% norm.plots) {
			switch(p,
				boxplot = {
					files$boxplot <- plot.boxplot(object,name=sample.list,is.norm=is.norm,output=output,path=path,...)
				},
				gcbias = {
					files$gcbias <- plot.edaseq(object,sample.list,covar=annotation$gc_content,is.norm=is.norm,which.plot=p,output=output,path=path,...)
				},
				lengthbias = {
					files$lengthbias <- plot.edaseq(object,sample.list,covar=annotation$end-annotation$start,is.norm=is.norm,which.plot=p,output=output,path=path,...)
				},
				meandiff = {
					fil <- plot.edaseq(object,sample.list,is.norm=is.norm,which.plot=p,output=output,path=path,...)
					for (n in names(fil))
						files$meandiff[[n]] <- unlist(fil[[n]])
				},
				meanvar = {
					fil <- plot.edaseq(object,sample.list,is.norm=is.norm,which.plot=p,output=output,path=path,...)
					for (n in names(fil))
						files$meanvar[[n]] <- unlist(fil[[n]])
				}
			)
		}
		if (p %in% stat.plots && is.norm) {
			for (cnt in names(contrast.list)) {
			disp("  Contrast: ",cnt)				
				samples <- names(unlist(contrast.list[[cnt]]))
				mat <- as.matrix(object[,match(samples,colnames(object))])
				switch(p,
					deheatmap = {
						files$deheatmap[[cnt]] <- plot.de.heatmap(mat,cnt,output=output,path=path)
					},
					volcano = {
						fc <- log2(make.fold.change(cnt,sample.list,object,1))
						for (contrast in colnames(fc)) {
							files$volcano[[contrast]] <- plot.volcano(fc[,contrast],p.list[[cnt]],contrast,fcut=thresholds$f,pcut=thresholds$p,output=output,path=path)
						}
					},
					biodist = {
						files$biodist[[cnt]] <- plot.noiseq(object,sample.list,covars,which.plot=p,output=output,biodist.opts=list(p=p.list[[cnt]],pcut=thresholds$p,name=cnt),path=path,...)
					}
				)
			}
		}
		if (p %in% other.plots) {
			switch(p,
				filtered = {
					files$filtered <- plot.filtered(object,annotation,output=output,path=path)
				}
			)
		}
	}
	
	return(files)
}

#' Boxplots wrapper for the metaseqr package
#'
#' A wrapper over the general boxplot function, suitable for matrices produced and processed with the metaseqr package. Intended for
#' internal use but can be easily used as stand-alone. It can colors boxes based on group depending on the name argument.
#'
#' @param mat the count data matrix.
#' @param name the names of the samples plotted on the boxplot. If NULL, the function check the column names of mat. If they are also
#' NULL, sample names are autogenerated. If name="none", no sample names are plotted. If name is a list, it should be the sample.list
#' argument provided to the manin metaseqr function. In that case, the boxes are colored per group.
#' @param log.it whether to log transform the values of mat or not. It can be TRUE, FALSE or "auto" for auto-detection. Auto-detection
#' log transforms by default so that the boxplots are smooth and visible.
#' @param y.lim custom y-axis limits. Leave the string "default" for default behavior.
#' @param is.norm a logical indicating whether object contains raw or normalized data. It is not essential and it serves only plot
#' annotation purposes.
#' @param output one or more R plotting device to direct the plot result to. Supported mechanisms: "x11" (default), "png", "jpg", 
#' "bmp", "pdf", "ps" or "json". The latter is currently available for the creation of interactive volcano plots only when reporting
#' the output, through the highcharts javascript library (JSON for boxplots not yet available).
#' @param path the path to create output files.
#' @param ... further arguments to be passed to plot devices, such as parameter from \code{\link{par}}.
#' @return The filename of the boxplot produced if it's a file.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2",B3"))
#' plot.boxplot(data.matrix,sample.list)
#'
#' norm.args <- get.defaults("normalization","deseq")
#' object <- normalize.deseq(data.matrix,sample.list,norm.args)
#' plot.boxplot(object,sample.list)
#'}
plot.boxplot <- function(mat,name=NULL,log.it="auto",y.lim="default",is.norm=FALSE,output="x11",path=NULL,...) {
	if (is.null(path)) path <- getwd()
	if (is.norm)
		status<- "normalized"
	else
		status<- "raw"
	# Need to log?
	if (log.it=="auto") {
		if (diff(range(mat,na.rm=TRUE))>1000)
			mat <- log2disp(mat)
	}
	else if (log.it=="yes")
		mat <- log2disp(mat)
	# Define the axis limits based on user input
	if (!is.numeric(y.lim) && y.lim=="auto") {
		min.y <- floor(min(mat))
		max.y <- ceiling(max(mat))
	}
	else if (is.numeric(y.lim)) {
		min.y <- y.lim[1]
		max.y <- y.lim[2]
	}
	grouped <- FALSE
	if (is.null(name)) {
		if (is.null(colnames(mat)))
			nams <- paste("Sample",1:ncol(mat),sep=" ")
		else
			nams <- colnames(mat)
	}
	else if (length(name)==1 && name=="none")
		nams <- rep("",ncol(mat))
	else if (is.list(name)) { # Is sample.list
		nams <- unlist(name)
		grouped <- TRUE
	}
	cols <- c("red3","green3","blue2","gold","skyblue","orange3","burlywood")
	if (grouped) {
		tmp <- as.numeric(factor(as.class.vector(name)))
		b.cols <- cols[tmp]
	}
	else b.cols <- cols
	mat.list <- list()
	for (i in 1:ncol(mat))
		mat.list[[i]] <- mat[,i]
	fil <- file.path(path,paste("boxplot_",status,".",output,sep=""))
	graphics.open(output,fil)
	if (!is.numeric(y.lim) && y.lim=="default")
		boxplot(mat.list,names=nams,col=b.cols,las=2,main=paste("Boxplot ",status,sep=""),...)
	else
		boxplot(mat.list,names=nams,col=b.cols,ylim=c(min.y,max.y),las=2,main=paste("Boxplot ",status,sep=""),...)
	graphics.close(output)
	return(fil)
}

#' Multi-Dimensinal Scale plots or RNA-Seq samples
#'
#' Creates a Multi-Dimensional Scale plot for the given samples based on the count data matrix. MDS plots are very useful for quality
#' control as you can easily see of samples of the same groups are clustered together based on the whole dataset.
#'
#' @param x the count data matrix.
#' @param sample.list the list containing condition names and the samples under each condition.
#' @param method which correlation method to use. Same as the method parameter in \code{\link{cor}} function.
#' @param log.it whether to log transform the values of x or not.
#' @param output one or more R plotting device to direct the plot result to. Supported mechanisms: "x11" (default), "png", "jpg", 
#' "bmp", "pdf", "ps" or "json". The latter is currently available for the creation of interactive volcano plots only when reporting
#' the output, through the highcharts javascript library.
#' @param path the path to create output files.
#' @param ... further arguments to be passed to plot devices, such as parameter from \code{\link{par}}.
#' @return The filename of the MDS plot produced if it's a file.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2",B3"))
#' plot.mds(data.matrix,sample.list)
#'}
plot.mds <- function(x,sample.list,method="spearman",log.it=TRUE,output="x11",path=NULL,...) {
	if (is.null(path)) path <- getwd()
	classes <- as.factor(as.class.vector(sample.list))
	design <- as.numeric(classes)
	colspace <- c("red","blue","yellowgreen","orange","aquamarine2",
				  "pink2","seagreen4","brown","purple","chocolate")
	pchspace <- c(20,17,15,16,8,3,2,0,1,4)
	if (log.it)
		y <- nat2log(x,base=2)
	else
		y <- x
	d <- as.dist(0.5*(1-cor(y,method=method)))
	mds.obj <- cmdscale(d,eig=TRUE,k=2)
	fil <- file.path(path,paste("mds.",output,sep=""))
	if (output %in% c("pdf","ps","x11"))
		graphics.open(output,fil,width=9,height=7)
	else
		graphics.open(output,fil,width=1024,height=768)
	xr <- diff(range(min(mds.obj$points[,1]),max(mds.obj$points[,1])))
	yr <- diff(range(min(mds.obj$points[,2]),max(mds.obj$points[,2])))
	xlim <- c(min(mds.obj$points[,1])-xr/10,max(mds.obj$points[,1])+xr/10)
	ylim <- c(min(mds.obj$points[,2])-yr/10,max(mds.obj$points[,2])+yr/10)
	plot(mds.obj$points[,1],mds.obj$points[,2],
		 col=colspace[1:length(levels(classes))][design],
		 pch=pchspace[1:length(levels(classes))][design],
		 xlim=xlim,ylim=ylim,
		 main="MDS plot",xlab="MDS 1",ylab="MDS 2",
		 cex=0.9,cex.lab=0.9,cex.axis=0.9,cex.main=0.9)
	text(mds.obj$points[,1],mds.obj$points[,2],labels=colnames(x),pos=3,cex=0.7)
	grid()
	graphics.close(output)
	return(fil)
}

#' Diagnostic plots based on the EDASeq package
#'
#' A wrapper around the plotting functions availale in the EDASeq normalization Bioconductor package. For analytical explanation of
#' each plot please see the vignette of the EDASeq package. It is best to use this function through the main plotting function
#' \code{\link{plot.metaseqr}}.
#'
#' @param x the count data matrix.
#' @param sample.list the list containing condition names and the samples under each condition.
#' @param covar The covariate to plot counts against. Usually "gc" or "length".
#' @param is.norm a logical indicating whether object contains raw or normalized data. It is not essential and it serves only plot
#' annotation purposes.
#' @param which.plot the EDASeq package plot to generate. It can be one or more of "meanvar", "meandiff", "gcbias" or "lengthbias".
#' Please refer to the documentation of the EDASeq package for details on the use of these plots. The which.plot="lengthbias" case
#' is not covered by EDASeq documentation, however it is similar to the GC-bias plot when the covariate is the gene length instead
#' of the GC content.
#' @param output one or more R plotting device to direct the plot result to. Supported mechanisms: "x11" (default), "png", "jpg", 
#' "bmp", "pdf" or "ps".
#' @param path the path to create output files.
#' @param ... further arguments to be passed to plot devices, such as parameter from \code{\link{par}}.
#' @return The filenames of the plot produced in a named list with names the which.plot argument. If output="x11", no output filenames
#' are produced.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2",B3"))
#' gc <- runif(nrow(data.matrix))
#' plot.edaseq(data.matrix,sample.list,covar=gc,which.plot=c("meanvar","meandiff","gcbias"))
#'}
plot.edaseq <- function(x,sample.list,covar=NULL,is.norm=FALSE,which.plot=c("meanvar","meandiff","gcbias","lengthbias"),output="x11",path=NULL,...) {
	if (is.null(path)) path <- getwd()
	check.text.args("which.plot",which.plot,c("meanvar","meandiff","gcbias","lengthbias"),multiarg=TRUE)
	if (is.null(covar) && which.plot %in% c("gcbias","lengthbias"))
		stop("\"covar\" argument is required when \"which.plot\" is \"gcbias\" or \"lengthbias\"!")
	if (is.norm)
		status <- "normalized"
	else
		status <- "raw"
	if (is.null(covar)) covar <- rep(NA,nrow(x))
	s <- newSeqExpressionSet(x,phenoData=AnnotatedDataFrame(data.frame(conditions=as.class.vector(sample.list),row.names=colnames(x))),
		featureData=AnnotatedDataFrame(data.frame(gc=covar,length=covar,row.names=rownames(x))))
	switch(which.plot,
		meandiff = {
			fil <- vector("list",length(sample.list))
			names(fil) <- names(sample.list)
			for (n in names(sample.list)) {
				if (length(sample.list[[n]])==1) {
					warning("Cannot create a mean-difference plot with one sample per condition! Skipping...",call.=FALSE)
					next
				}
				pair.matrix <- combn(1:length(sample.list[[n]]),2)
				fil[[n]] <- vector("list",ncol(pair.matrix))
				for (i in ncol(pair.matrix)) {
					s1 <- sample.list[[n]][pair.matrix[1,i]]
					s2 <- sample.list[[n]][pair.matrix[2,i]]
					fil[[n]][[i]] <- file.path(path,paste(which.plot,"_",status,"_",n,"_",s1,"_",s2,".",output,sep=""))
					names(fil[[n]][i]) <- paste(s1,"vs",s2,sep="_")
					graphics.open(output,fil[[n]][[i]])
					MDPlot(s,y=pair.matrix[,i],main=paste("MD plot for ",n," ",status," samples ",s1," and ",s2,sep=""),cex.main=0.9)
					graphics.close(output)
				}
			}
		},
		meanvar = {
			fil <- vector("list",length(sample.list))
			names(fil) <- names(sample.list)
			for (n in names(sample.list)) {	
				if (length(sample.list[[n]])==1) {
					warning("Cannot create a mean-variance plot with one sample per condition! Skipping...",call.=FALSE)
					next
				}
				pair.matrix <- combn(1:length(sample.list[[n]]),2)
				fil[[n]] <- vector("list",ncol(pair.matrix))
				for (i in ncol(pair.matrix)) {
					s1 <- sample.list[[n]][pair.matrix[1,i]]
					s2 <- sample.list[[n]][pair.matrix[2,i]]
					fil[[n]][[i]] <- file.path(path,paste(which.plot,"_",status,"_",n,"_",s1,"_",s2,".",output,sep=""))
					names(fil[[n]][i]) <- paste(s1,"vs",s2,sep="_")
					graphics.open(output,fil[[n]][[i]])
					meanVarPlot(s,main=paste("MV plot for ",n," ",status," samples ",s1," and ",s2,sep=""),cex.main=0.9)
					graphics.close(output)
				}
			}
		},
		gcbias = {
			fil <- file.path(path,paste(which.plot,"_",status,".",output,sep=""))
			graphics.open(output,fil)
			biasPlot(s,"gc",xlim=c(0.1,0.9),log=TRUE,ylim=c(0,15),main=paste("Expression - GC content ",status,sep=""))
			grid()
			graphics.close(output)
		},
		lengthbias = {
			fil <- file.path(path,paste(which.plot,"_",status,".",output,sep=""))
			graphics.open(output,fil)
			biasPlot(s,"length",log=TRUE,ylim=c(0,10),main=paste("Expression - Gene length ",status,sep=""))
			grid()
			graphics.close(output)
		}
	)
	return(fil)
}

#' Diagnostic plots based on the NOISeq package
#'
#' A wrapper around the plotting functions availale in the NOISeq RNA-Seq analysisBioconductor package. For analytical explanation
#' of each plot please see the vignette of the NOISeq package. It is best to use this function through the main plotting function
#' \code{\link{plot.metaseqr}}.
#'
#' @param x the count data matrix.
#' @param sample.list the list containing condition names and the samples under each condition.
#' @param covars a list (whose annotation elements are ideally a subset of an annotation data frame produced by \code{\link{get.annotation}})
#' with the following members: data (the data matrix), length (gene length), gc (the gene gc_content), chromosome (a data frame with
#' chromosome name and co-ordinates), factors (a factor with the experimental condition names replicated by the number of samples in
#' each experimental condition) and biotype (each gene's biotype as depicted in Ensembl-like annotations).
#' @param which.plot the NOISeq package plot to generate. It can be one or more of "biodetection", "countsbio", "saturation", "rnacomp"
#' or "biodist". Please refer to the documentation of the EDASeq package for details on the use of these plots. The which.plot="saturation"
#' case is modified to be more informative by producing two kinds of plots. See \code{\link{plot.noiseq.saturation}}.
#' @param biodist.opts a list with the following members: p (a vector of p-values, e.g. the p-values of a contrast), pcut (a unique
#' number depicting a p-value cutoff, required for the "biodist" case), name (a name for the "biodist" plot, e.g. the name of the
#' contrast.
#' @param output one or more R plotting device to direct the plot result to. Supported mechanisms: "x11" (default), "png", "jpg", 
#' "bmp", "pdf" or "ps".
#' @param path the path to create output files.
#' @param ... further arguments to be passed to plot devices, such as parameter from \code{\link{par}}.
#' @return The filenames of the plots produced in a named list with names the which.plot argument. If output="x11", no output filenames
#' are produced.
#' @note Please note that in case of "biodist" plots, the behavior of the function is unstable, mostly due to the very specific inputs
#' this plotting function accepts in the NOISeq package. We have tried to predict unstable behavior and avoid exceptions through the
#' use of tryCatch but it's still possible that you might run onto an error.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2",B3"))
#' lengths <- round(1000*runif(nrow(data.matrix)))
#' starts <- round(1000*runif(nrow(data.matrix)))
#' ends <- starts + lengths
#' covars <- list(
#'   data=data.matrix,
#'   length=lengths,
#'   gc=runif(nrow(data.matrix)),
#'   chromosome=data.frame(
#'     chromosome=c(rep("chr1",nrow(data.matrix)/2),rep("chr2",nrow(data.matrix)/2)),
#'     start=starts,
#'     end=ends
#'   ),
#'   factors=data.frame(class=as.class.vector(sample.list)),
#'   biotype=c(rep("protein_coding",nrow(data.matrix)/2),rep("ncRNA",nrow(data.matrix)/2))
#' )
#' p <- runif(nrow(data.matrix))
#' plot.noiseq(data.matrix,sample.list,covars=covars,biodist.opts=list(p=p,pcut=0.1,name="A_vs_B"))
#'}
plot.noiseq <- function(x,sample.list,covars,which.plot=c("biodetection", "countsbio", "saturation", "rnacomp", "biodist"),output="x11",
	biodist.opts=list(p=NULL,pcut=NULL,name=NULL),path=NULL,...) {
	if (is.null(path)) path <- getwd()
	# covars is a list of gc-content, factors, length, biotype, chromosomes, factors, basically copy of the noiseq object
	which.plot <- tolower(which.plot[1])
	check.text.args("which.plot",which.plot,c("biodetection","countsbio","saturation","rnacomp","biodist"),multiarg=FALSE)
	if (missing(covars))
		stop("\"covars\" argument is required with NOISeq specific plots!")
	else {
		covars$biotype <- as.character(covars$biotype)
		names(covars$length) <- names(covars$gc) <- rownames(covars$chromosome) <- names(covars$biotype) <- rownames(x)
	}
	if (which.plot=="biodist") {
		if (is.null(biodist.opts$p))
			stop("A p-value must be provided for the \"biodist\" plot!")
		if (is.null(biodist.opts$pcut) || is.na(biodist.opts$pcut)) biodist.opts$pcut=0.05
	}
	
	# All of these plots are NOISeq specific so we need a local NOISeq object
	if (any(is.na(unique(covars$biotype))))
		covars$biotype=NULL # Otherwise, it will probably crash
	local.obj <- NOISeq::readData(
		data=x,
		length=covars$gene.length,
		gc=covars$gc.content,
		chromosome=covars$chromosome,
		#factors=data.frame(class=covars$factors),
		factors=covars$factors,
		biotype=covars$biotype
	)
	switch(which.plot,
		biodetection = {
			plot.data <- NOISeq::dat(local.obj,type=which.plot)
			samples <- unlist(sample.list)
			fil <- character(length(samples))
			names(fil) <- samples
			for (i in 1:length(samples)) {
				fil[samples[i]] <- file.path(path,paste(which.plot,"_",samples[i],".",output,sep=""))
				if (output %in% c("pdf","ps","x11"))
					graphics.open(output,fil[samples[i]],width=9,height=7)
				else
					graphics.open(output,fil[samples[i]],width=1024,height=768)
				explo.plot(plot.data,samples=i)
				graphics.close(output)
			}
		},
		countsbio = {
			plot.data <- NOISeq::dat(local.obj,type=which.plot,factor=NULL)
			samples <- unlist(sample.list)
			fil <- character(length(samples))
			names(fil) <- samples
			for (i in 1:length(samples)) {
				fil[samples[i]] <- file.path(path,paste(which.plot,"_",samples[i],".",output,sep=""))
				if (output %in% c("pdf","ps","x11"))
					graphics.open(output,fil[samples[i]],width=9,height=7)
				else
					graphics.open(output,fil[samples[i]],width=1024,height=768)
				explo.plot(plot.data,samples=i,plottype="boxplot")
				graphics.close(output)
			}
		},
		saturation = { # Up to 12 samples... Hmmm... I have to improvise
			plot.data <- NOISeq::dat(local.obj,k=0,ndepth=9,type=which.plot) # For 10 saturation points
			d2s <- dat2save(plot.data)
			fil <- plot.noiseq.saturation(d2s,output,covars$biotype,path=path)
		},
		rnacomp = {
			plot.data <- NOISeq::dat(local.obj,type="cd")
			fil <- file.path(path,paste(which.plot,".",output,sep=""))
			graphics.open(output,fil)
			explo.plot(plot.data,main="RNA composition")
			grid()
			graphics.close(output)
		},
		biodist = { # We have to fake a noiseq object
			p <- biodist.opts$p
			if (is.matrix(p)) p <- p[,1]
			dummy <- new("Output",
				comparison=c("Dummy.1","Dummy.2"),
				factor=c("class"),
				k=1,
				lc=1,
				method="n",
				replicates="biological",
				results=list(
					data.frame(
						Dummy.1=rep(1,length(p)),
						Dummy.2=rep(1,length(p)),
						M=rep(1,length(p)),
						D=rep(1,length(p)),
						prob=as.numeric(p),
						ranking=rep(1,length(p)),
						Length=rep(1,length(p)),
						GC=rep(1,length(p)),
						Chrom=as.character(covars$chromosome[,1]),
						GeneStart=covars$chromosome[,2],
						GeneEnd=covars$chromosome[,3],
						Biotype=covars$biotype
					)
				),
				nss=5,
				pnr=0.2,
				v=0.02
			)
			if (!is.null(biodist.opts$name))
				fil <- file.path(path,paste(which.plot,"_",biodist.opts$name,".",output,sep=""))
			else
				fil <- file.path(path,paste(which.plot,".",output,sep=""))
			if (output %in% c("pdf","ps","x11"))
				graphics.open(output,fil,width=10,height=6)
			else
				graphics.open(output,fil,width=1024,height=640)
			tryCatch( # A lot of times, there is a problem with this function
				DE.plot(dummy,chromosomes=NULL,q=biodist.opts$pcut,graphic="distr"),
				error=function(e) {
					disp("      Known problem with NOISeq and external p-values detected! Trying to make a plot with alternative p-values (median of p-value distribution)...")
					tryCatch(
						DE.plot(dummy,chromosomes=NULL,q=quantile(biodist.opts$p,0.5),graphic="distr"),
						error=function(e) {
							disp("      Cannot create DEG biotype plot! This is not related to a problem with the results. Excluding...")
						},
						finally=""
					)
				},
				finally=""
			)
			graphics.close(output)
		}
	)
	return(fil)
}

#' Simpler implementation of saturation plots inspired from NOISeq package
#'
#' Helper function for \code{\link{plot.noiseq}} to plot feature detection saturation as presented in the NOISeq package vignette.
#' It has two main outputs: a set of figures, one for each input sample depicting the saturation for each biotype and one single
#' multiplot which depicts the saturation of all samples for each biotype. It expands the saturation plots of NOISeq by allowing
#' more samples to be examined in a simpler way. Don't use this function directly. Use either \code{\link{plot.metaseqr}} or 
#' \code{\link{plot.noiseq}}.
#'
#' @param x the count data matrix.
#' @param o one or more R plotting device to direct the plot result to. Supported mechanisms: "x11" (default), "png", "jpg", 
#' "bmp", "pdf" or "ps".
#' @param tb the vector of biotypes, one for each row of x.
#' @param path the path to create output files.
#' @return The filenames of the plots produced in a named list with names the which.plot argument. If output="x11", no output filenames
#' are produced.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' biotype=c(rep("protein_coding",nrow(data.matrix)/2),rep("ncRNA",nrow(data.matrix)/2))
#' plot.noiseq.saturation(data.matrix,"x11",biotype)
#'}
plot.noiseq.saturation <- function(x,o,tb,path=NULL) {
	if (is.null(path)) path <- getwd()
	total.biotypes <- table(tb)
	the.biotypes <- names(tb)
	biotypes <- colnames(x[[1]][,2:ncol(x[[1]])])
	colspace <- c("red3","green4","blue2","orange3","burlywood",
				  "lightpink4","gold","skyblue","red2","green2","firebrick3",
				  "orange4","yellow4","skyblue3","tan4","gray40",
				  "brown2","darkgoldenrod","cyan3","coral2","cadetblue",
				  "bisque3","blueviolet","chocolate3","darkkhaki","dodgerblue")
	pchspace <- c(rep(c(15,16,17,18),6),15)

	# Plot all biotypes per sample
	f.sample <- character(length(names(x)))
	names(f.sample) <- names(x)
	for (n in names(x)) {
		f.sample[n] <- file.path(path,paste("saturation_",n,".",o,sep=""))
		if (o %in% c("pdf","ps","x11"))
			graphics.open(o,f.sample[n],width=10,height=7)
		else
			graphics.open(o,f.sample[n],width=1024,height=800)
		y <- x[[n]]
		sep <- match(c("global","protein_coding"),colnames(y))
		yab <- cbind(y[,"depth"],y[,sep])
		ynab <- y[,-sep]
		colnames(yab)[1] <- colnames(ynab)[1] <- "depth"
		xlim <- range(y[,"depth"])
		ylim.ab <- range(yab[,2:ncol(yab)])
		ylim.nab <- range(ynab[,2:ncol(ynab)])
		par(cex.axis=0.9,cex.main=1,cex.lab=0.9,font.lab=2,font.axis=2,pty="m",lty=2,lwd=1.5,mfrow=c(1,2))
		plot.new()
		plot.window(xlim,ylim.nab)
		axis(1,at=pretty(xlim,10),labels=as.character(pretty(xlim,10)/1e+6))
		axis(2,at=pretty(ylim.nab,10))
		title(main="Non abundant biotype detection saturation",xlab="Depth in millions of reads",ylab="Detected features")
		co <- 0
		for (b in biotypes) {
			co <- co + 1
			if (b=="global" || b=="protein_coding") {
				# Silently do nothing
			}
			else {
				lines(ynab[,"depth"],ynab[,b],col=colspace[co])
				points(ynab[,"depth"],ynab[,b],pch=pchspace[co],col=colspace[co],cex=1)
			}
		}
		grid()
		legend(
			x="topleft",legend=colnames(ynab)[2:ncol(ynab)],xjust=1,yjust=0,
			box.lty=0,x.intersp=0.5,cex=0.6,text.font=2,
			col=colspace[1:(ncol(ynab)-1)],pch=pchspace[1:(ncol(ynab)-1)]
		)
		plot.new()
		plot.window(xlim,ylim.ab)
		axis(1,at=pretty(xlim,10),labels=as.character(pretty(xlim,10)/1e+6))
		axis(2,at=pretty(ylim.ab,10))
		title(main="Abundant biotype detection saturation",xlab="Depth in millions of reads",ylab="Detected features")
		co <- 0
		for (b in c("global","protein_coding")) {
			co <- co + 1
			lines(yab[,"depth"],yab[,b],col=colspace[co])
			points(yab[,"depth"],yab[,b],pch=16,col=colspace[co],cex=1.2)
		}
		grid()
		legend(
			x="topleft",legend=c("global","protein_coding"),xjust=1,yjust=0,
			box.lty=0,lty=2,x.intersp=0.5,cex=0.7,text.font=2,
			col=colspace[1:2],pch=pchspace[1:2]
		)
		mtext(n,side=3,line=-1.5,outer=TRUE,font=2,cex=1.3)
		graphics.close(o)
	}

	# Plot all samples per biotype
	g <- make.grid(length(biotypes))
	f.all <- file.path(path,paste("biotype_saturation.",o,sep=""))
	if (o %in% c("pdf","ps"))
		graphics.open(o,f.all,width=12,height=12)
	else
		graphics.open(o,f.all,width=1024,height=1024)
	par(cex.axis=0.8,cex.main=0.9,cex.lab=0.8,pty="m",lty=2,lwd=1.5,mfrow=g,mar=c(3,3,1,1),oma=c(1,1,0,0),mgp=c(2,0.5,0))
	for (b in biotypes) {
		y <- depth <- vector("list",length(x))
		names(y) <- names(depth) <- names(x)
		for (n in names(x)) {
			y[[n]] <- x[[n]][,b]
			depth[[n]] <- x[[n]][,"depth"]
		}
		y <- do.call("cbind",y)
		xlim <- range(do.call("c",depth))
		ylim <- range(y)
		plot.new()
		plot.window(xlim,ylim)
		axis(1,at=pretty(xlim,5),labels=as.character(pretty(xlim,5)/1e+6),line=0.5)
		axis(2,at=pretty(ylim,5),line=0.5)
		title(main=b,xlab="Depth in millions of reads",ylab="Detected features")
		co <- 0
		for (n in colnames(y)) {
			co <- co + 1
			lines(depth[[n]],y[,n],col=colspace[co])
			points(depth[[n]],y[,n],pch=pchspace[co],col=colspace[co])
		}
		grid()
		legend(
			x="bottomright",legend=colnames(y),xjust=1,yjust=0,
			box.lty=0,x.intersp=0.5,
			col=colspace[1:length(colnames(y))],pch=pchspace[1:length(colnames(y))]
		)
	}
	graphics.close(o)

	return(list(sample=f.sample,biotype=f.all))
}

#' (Interactive) volcano plots of differentially expressed genes
#'
#' This function plots a volcano plot or returns a JSON string which is used to render aninteractive in case of HTML reporting.
#'
#' @param f the fold changes which are to be plotted on the x-axis.
#' @param p the p-values whose -log10 transformation is going to be plotted on the y-axis.
#' @param con an optional string depicting a name (e.g. the contrast name) to appear in the title of the volcano plot.
#' @param fcut a fold change cutoff so as to draw two vertical lines indicating the cutoff threshold for biological significance.
#' @param pcut a p-value cutoff so as to draw a horizontal line indicating the cutoff threshold for statistical significance.
#' @param alt.names an optional vector of names, e.g. HUGO gene symbols, alternative or complementary to the unique names of f or p (one of
#' them must be named!). It is used only in JSON output.
#' @param output one or more R plotting device to direct the plot result to. Supported mechanisms: "x11" (default), "png", "jpg", 
#' "bmp", "pdf", "ps" or "json". The latter is currently available for the creation of interactive volcano plots only when reporting
#' the output, through the highcharts javascript library.
#' @param path the path to create output files.
#' @param ... further arguments to be passed to plot devices, such as parameter from \code{\link{par}}.
#' @return The filenames of the plots produced in a named list with names the which.plot argument. If output="x11", no output filenames
#' are produced.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2",B3"))
#' contrast <- "A_vs_B"
#' M <- norm.edger(data.matrix,sample.list)
#' p <- stat.edger(M,sample.list,contrast)
#' ma <- apply(M[,sample.list$A],1,mean)
#' mb <- apply(M[,sample.list$B],1,mean)
#' f <- log2(mb/ma)
#' plot.volcano(f,p,con=contrast)
#' j <- plot.volcano(f,p,con=contrast,output="json")
#'}
plot.volcano <- function(f,p,con=NULL,fcut=1,pcut=0.05,alt.names=NULL,output="x11",path=NULL,...) { # output can be json here...
	if (is.null(path)) path <- getwd()
	if (is.null(con))
		con <- conn <- ""
	else {
		conn <- con
		con <- paste("for ",con)
	}
	fil <- file.path(path,paste("volcano_plot_",conn,".",output,sep=""))
	if (output!="json") {
		if (output %in% c("pdf","ps","x11"))
			graphics.open(output,fil,width=8,height=10)
		else
			graphics.open(output,fil,width=768,height=1024,res=100)
	}
	rem <- which(is.na(p))
	if (length(rem)>0) {
		p <- p[-rem]
		f <- f[-rem]
		if (!is.null(alt.names))
			alt.names <- alt.names[-rem]
	}
	p.zero <- which(p==0)
	if (length(p.zero)>0) # Fix problem with extremely low p-values, only for display purposes though
		p[p.zero] <- runif(length(p.zero),0,1e-256)
	xlim <- c(-max(abs(f)),max(abs(f)))
	ylim <- c(0,ceiling(-log10(min(p))))
	up <- which(f>=fcut & p<pcut)
	down <- which(f<=-fcut & p<pcut)
	u <- union(up,down)
	alt.names.neutral <- NULL
	if (length(u)>0) {
		ff <- f[-u]
		pp <- p[-u]
		if (!is.null(alt.names))
			alt.names.neutral <- alt.names[-u]
	}
	else {
		ff <- f
		pp <- p
		if (!is.null(alt.names))
			alt.names.neutral <- alt.names
	}
	if (output!="json") {
		par(cex.main=1.1,cex.lab=1.1,cex.axis=1.1,font.lab=2,font.axis=2,pty="m",lwd=1.5)
		plot.new()
		plot.window(xlim,ylim)
		axis(1,at=pretty(xlim,10),labels=as.character(pretty(xlim,10)))
		axis(2,at=pretty(ylim,10))
		title(paste(main="Volcano plot",con),xlab="Fold change",ylab="-log10(p-value)")
		points(ff,-log10(pp),pch=20,col="blue2",cex=0.9)
		points(f[down],-log10(p[down]),pch=20,col="green3",cex=0.9)
		points(f[up],-log10(p[up]),pch=20,col="red2",cex=0.9)
		abline(h=-log10(pcut),lty=4)
		abline(v=-fcut,lty=2)
		abline(v=fcut,lty=2)
		grid()
		legend(
			x="topleft",
			legend=c("up-regulated","down-regulated","unregulated","p-value threshold","fold change threshold"),
			col=c("red2","green3","blue1","black","black"),
			pch=c(20,20,20,NA,NA),lty=c(NA,NA,NA,4,2),
			xjust=1,yjust=0,box.lty=0,x.intersp=0.5,cex=0.8,text.font=2
		)
		graphics.close(output)
		return(fil)
	}
	else {
		if (is.null(alt.names))
			point.format="<b>id: </b>{point.name}<br><b>fold change: </b>{point.x}<br><b>significance: </b>{point.y}"
		else
			point.format="<b>name: </b>{point.alt_name}<br><b>id: </b>{point.name}<br><b>fold change: </b>{point.x}<br><b>significance: </b>{point.y}"
		json <- toJSON(
			list(
					chart=list(
					type="scatter",
					zoomType="xy"
				),
				title=list(
					text=paste("Volcano plot for",con)
				),
				xAxis=list(
					title=list(
						enabled=TRUE,
						text="Fold change"
					),
					startOnTick=TRUE,
					endOnTick=TRUE,
					showLastLabel=TRUE,
					gridLineWidth=1,
					min=xlim[1],
					max=xlim[2]
				),
				yAxis=list(
					title=list(
						enabled=TRUE,
						text="-log10(p-value)"
					),
					startOnTick=TRUE,
					endOnTick=TRUE,
					showLastLabel=TRUE,
					gridLineWidth=1,
					min=ylim[1]-2,
					max=ylim[2]
				),
				#legend=list(
				#	layout="vertical",
				#	align="left",
				#	verticalAlign="top",
				#	floating=TRUE,
				#	backgroundColor="#FFFFFF",
				#	borderWidth=1
				#),
				plotOptions=list(
					scatter=list(
						allowPointSelect=TRUE,
						marker=list(
							radius=2,
							states=list(
								hover=list(
									enabled=TRUE,
									lineColor="#333333"
								)
							)
						),
						states=list(
							hover=list(
								marker=list(
									enabled=FALSE
								)
							)
						),
						tooltip=list(
							headerFormat="<span style=\"font-size:12px; color:{series.color}\">{series.name}<br>",
							#pointFormat="<b>name: </b>{point.name}<br><b>fold change: </b>{point.x}<br><b>significance: </b>{point.y}"
							pointFormat=point.format
						),
						turboThreshold=10000
					)
				),
				series=list(
					list(
						name="up-regulated",
						color="#EE0000",
						data=make.highcharts.points(f[up],p[up],alt.names[up])
					),
					list(
						name="down-regulated",
						color="#00CD00",
						data=make.highcharts.points(f[down],p[down],alt.names[down])
					),
					list(
						name="unregulated",
						color="#0000EE",
						data=make.highcharts.points(ff,pp,alt.names.neutral)
					),
					list(
						name="downfold threshold",
						color="#000000",
						type="line",
						dashStyle="Dash",
						marker=list(
							enabled=FALSE
						),
						data=list(c(-fcut,ylim[1]-5),c(-fcut,ylim[2]))
					),
					list(
						name="upfold threshold",
						color="#000000",
						type="line",
						dashStyle="Dash",
						marker=list(
							enabled=FALSE
						),
						data=list(c(fcut,ylim[1]-5),c(fcut,ylim[2]))
					),
					list(
						name="significance threshold",
						color="#000000",
						type="line",
						dashStyle="DashDot",
						marker=list(
							enabled=FALSE
						),
						data=list(c(xlim[1],-log10(pcut)),c(xlim[2],-log10(pcut)))
					)
				)
			)
		)
		return(json)
	}
}

#' Diagnostic heatmap of differentially expressed genes
#'
#' This function plots a heatmap of the differentially expressed genes produced by the metaseqr workflow, useful for quality control,
#' e.g. whether samples belonging to the same group cluster together.
#'
#' @param x the data matrix to create a heatmap for.
#' @param con an optional string depicting a name (e.g. the contrast name) to appear in the title of the volcano plot.
#' @param output one or more R plotting device to direct the plot result to. Supported mechanisms: "x11" (default), "png", "jpg", 
#' "bmp", "pdf", "ps".
#' @param path the path to create output files.
#' @param ... further arguments to be passed to plot devices, such as parameter from \code{\link{par}}.
#' @return The filenames of the plots produced in a named list with names the which.plot argument. If output="x11", no output filenames
#' are produced.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2",B3"))
#' contrast <- "A_vs_B"
#' M <- norm.edger(data.matrix,sample.list)
#' p <- stat.edger(M,sample.list,contrast)
#' plot.de.heatmap(data.matrix[p[[1]]<0.05])
#'}
plot.de.heatmap <- function(x,con=NULL,output="x11",path=NULL,...) {
	if (is.null(path)) path <- getwd()
	if (is.null(con))
		con <- conn <- ""
	else {
		conn <- con
		con <- paste("for ",con)
	}
	y <- nat2log(x,2,1)
	fil <- file.path(path,paste("de_heatmap_",conn,".",output,sep=""))
	if (output %in% c("pdf","ps","x11"))
		graphics.open(output,fil,width=12,height=12)
	else
		graphics.open(output,fil,width=1024,height=1024)
	heatmap.2(y,trace="none",col=bluered(16),labRow="",cexCol=0.9,keysize=1,font.lab=2,main=paste("DEG heatmap",con),cex.main=0.9)
	graphics.close(output)
	return(fil)
}

#' Diagnostic plot for filtered genes
#'
#' This function plots a grid of four graphs depicting: in the first row, the numbers of filtered genes per chromosome in the first
#' column and per biotype in the second column. In the second row, the percentages of filtered genes  per chromosome related to the
#' whole genome in the first columns and per biotype in the second column.
#'
#' @param x an annotation data frame like the ones produced by \code{\link{get.annotation}}. x should be the filtered annotation
#' according to metaseqr's filters.
#' @param y an annotation data frame like the ones produced by \code{\link{get.annotation}}. x should contain the total annotation
#' without the application of any metaseqr filter.
#' @param output one or more R plotting device to direct the plot result to. Supported mechanisms: "x11" (default), "png", "jpg", 
#' "bmp", "pdf" or "ps".
#' @param path the path to create output files.
#' @param ... further arguments to be passed to plot devices, such as parameter from \code{\link{par}}.
#' @return The filenames of the plots produced in a named list with names the which.plot argument. If output="x11", no output filenames
#' are produced.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' y <- get.annotation("mm9","gene")
#' x <- y[-sample(1:nrow(y),10000),]
#' plot.filtered(x,y)
#'}
plot.filtered <- function(x,y,output="x11",path=NULL,...) {
	if (is.null(path)) path <- getwd()
	fil <- file.path(path,paste("filtered_genes.",output,sep=""))
	if (output %in% c("pdf","ps","x11"))
		graphics.open(output,fil,width=12,height=8)
	else
		graphics.open(output,fil,width=1200,height=800,res=100)
	chr <- table(as.character(x$chromosome))
	bt <- table(as.character(x$biotype))
	chr.all <- table(as.character(y$chromosome))
	bt.all <- table(as.character(y$biotype))
	barlab.chr <- as.character(chr)
	barlab.bt <- as.character(bt)
	per.chr <- chr/chr.all[names(chr)]
	per.bt <- bt/bt.all[names(bt)]
	# Some bug...
	per.chr[per.chr>1] <- 1
	per.bt[per.bt>1] <- 1
	#
	suppressWarnings(per.chr.lab <- paste(formatC(100*per.chr,digits=1,format="f"),"%",sep=""))
	suppressWarnings(per.bt.lab <- paste(formatC(100*per.bt,digits=1,format="f"),"%",sep=""))

	par(mfrow=c(2,2),mar=c(1,4,2,1),oma=c(1,1,1,1))

	# Chromosomes
	barx.chr <- barplot(chr,space=0.5,ylim=c(0,max(chr)+ceiling(max(chr)/10)),yaxt="n",xaxt="n",plot=FALSE)
	plot.new()
	plot.window(xlim=c(0,ceiling(max(barx.chr))),ylim=c(0,max(chr)+ceiling(max(chr)/10)),mar=c(1,4,1,1))
	axis(2,at=pretty(0:(max(chr)+ceiling(max(chr)/10))),cex.axis=0.9,padj=1,font=2)
	text(x=barx.chr,y=chr,label=barlab.chr,cex=0.7,font=2,col="green3",adj=c(0.5,-1.3))
	title(main="Filtered genes per chromosome",cex.main=1.1)
	mtext(side=2,text="Number of genes",line=2,cex=0.9,font=2)
	grid()
	barplot(chr,space=0.5,ylim=c(0,max(chr)+ceiling(max(chr)/10)),col="blue3",border="yellow3",yaxt="n",xaxt="n",font=2,add=TRUE)

	# Biotypes
	barx.bt <- barplot(bt,space=0.5,ylim=c(0,max(bt)+ceiling(max(bt)/10)),yaxt="n",xaxt="n",plot=FALSE)
	plot.new()
	plot.window(xlim=c(0,ceiling(max(barx.bt))),ylim=c(0,max(bt)+ceiling(max(bt)/10)),mar=c(1,4,1,1))
	axis(2,at=pretty(0:(max(bt)+ceiling(max(bt)/10))),cex.axis=0.9,padj=1,font=2)
	text(x=barx.bt,y=bt,label=barlab.bt,cex=0.7,font=2,col="blue",adj=c(0.5,-1.3),xpd=TRUE)
	title(main="Filtered genes per biotype",cex.main=1.1)
	mtext(side=2,text="Number of genes",line=2,cex=0.9,font=2)
	grid()
	barplot(bt,space=0.5,ylim=c(0,max(bt)+ceiling(max(bt)/10)),col="red3",border="yellow3",yaxt="n",xaxt="n",font=2,add=TRUE)

	# Chromosome percentage
	barx.per.chr <- barplot(per.chr,space=0.5,ylim=c(0,max(per.chr)),yaxt="n",xaxt="n",plot=FALSE)
	plot.new()
	par(mar=c(9,4,1,1))
	plot.window(xlim=c(0,max(barx.per.chr)),ylim=c(0,max(per.chr)))
	#axis(1,at=barx.per.chr,labels=names(per.chr),cex.axis=0.9,font=2,tcl=-0.3,col="lightgrey",las=2)
	axis(1,at=barx.per.chr,labels=FALSE,tcl=-0.3,col="lightgrey")
	axis(2,at=seq(0,max(per.chr),length.out=5),labels=formatC(seq(0,max(per.chr),length.out=5),digits=2,format="f"),cex.axis=0.9,padj=1,font=2)
	text(barx.per.chr,par("usr")[3]-max(per.chr)/17,labels=names(per.chr),srt=45,adj=c(1,1.1),xpd=TRUE,cex=0.9,font=2)
	text(x=barx.per.chr,y=per.chr,label=per.chr.lab,cex=0.7,font=2,col="green3",adj=c(0.5,-1.3),xpd=TRUE)
	mtext(side=2,text="fraction of total genes",line=2,cex=0.9,font=2)
	grid()
	barplot(per.chr,space=0.5,ylim=c(0,max(per.chr)),col="blue3",border="yellow3",yaxt="n",xaxt="n",font=2,add=TRUE)

	# Biotype percentage
	barx.per.bt <- barplot(per.bt,space=0.5,ylim=c(0,max(per.bt)),yaxt="n",xaxt="n",plot=FALSE)
	plot.new()
	par(mar=c(9,4,1,1))
	plot.window(xlim=c(0,max(barx.per.bt)),ylim=c(0,max(per.bt)))
	#axis(1,at=barx.per.bt,labels=names(per.bt),cex.axis=0.9,font=2,tcl=-0.3,col="lightgrey",las=2)
	axis(1,at=barx.per.bt,labels=FALSE,tcl=-0.3,col="lightgrey")
	axis(2,at=seq(0,max(per.bt),length.out=5),labels=formatC(seq(0,max(per.bt),length.out=5),digits=2,format="f"),cex.axis=0.9,padj=1,font=2)
	text(barx.per.bt,par("usr")[3]-max(per.bt)/17,labels=gsub("prime","'",names(per.bt)),srt=45,adj=c(1,1.1),xpd=TRUE,cex=0.9,font=2)
	text(x=barx.per.bt,y=per.bt,label=per.bt.lab,cex=0.7,font=2,col="blue",adj=c(0.5,-1.3),xpd=TRUE)
	mtext(side=2,text="fraction of total genes",line=2,cex=0.9,font=2)
	grid()
	barplot(per.bt,space=0.5,ylim=c(0,max(per.bt)),col="red3",border="yellow3",yaxt="n",xaxt="n",font=2,add=TRUE)
	
	graphics.close(output)

	return(fil)
}

#' Open plotting device
#'
#' Wrapper function to open a plotting device. Internal use only.
#'
#' @param o the plotting device, see main metaseqr function
#' @param f a filename, if the plotting device requires it (e.g. "pdf")
#' @param ... further arguments to be passed to plot devices, such as parameter from \code{\link{par}}.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' graphics.open("pdf","test.pdf",width=12,height=12)
#'}
graphics.open <- function(o,f,...) {
	if(!check.graphics.type(o))
		stop("Invalid graphics output type!")
	if(check.graphics.file(o) && is.null(f))
		stop("Please specify an output file name for your plot")
	
	switch(o,
		x11 = { x11(...) },
		pdf = { pdf(file=f,pointsize=10,...) },
		ps = { postscript(file=f,pointsize=10,...) },
		png = { png(file=f,pointsize=12,...) },
		jpg = { jpeg(file=f,pointsize=12,quality=100,...) },
		bmp = { bmp(file=f,pointsize=12,...) },
		tiff = { tiff(file=f,pointsize=12,...) }
	)
}

#' Close plotting device
#'
#' Wrapper function to close a plotting device. Internal use only.
#'
#' @param o the plotting device, see main metaseqr function
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' graphics.close("pdf")
#'}
graphics.close <- function(o) {
	if (!is.element(o,c("x11","png","jpg","tiff","bmp","pdf","ps")))
		return(FALSE)
	if (o!="x11")
		dev.off()
}

#' Check plotting device
#'
#' Plotting device checker. Internal use only.
#'
#' @param o the plotting device, see main metaseqr function
#' @author Panagiotis Moulos
check.graphics.type <- function(o) {
	if (!is.element(o,c("x11","png","jpg","tiff","bmp","pdf","ps")))
		return(FALSE)
	else
		return(TRUE)
}

#' Check graphics file
#'
#' Graphics file checker. Internal use only.
#'
#' @param o the plotting device, see main metaseqr function
#' @author Panagiotis Moulos
check.graphics.file <- function(o) {
	if (is.element(o,c("png","jpg","tiff","bmp","pdf","ps")))
		return(TRUE)
	else
		return(FALSE)
}

#' Display value transformation
#'
#' Logarithmic transformation for display purposes. Internal use only.
#'
#' @param mat input data matrix
#' @param base logarithmic base, 2 or 10
#' @author Panagiotis Moulos
log2disp <- function(mat,base=2) {
	mat[mat==0] <- 1
	if (base==10)
		return(log10(mat))
	else
		return(log2(mat))
}

#' General value transformation
#'
#' Logarithmic transformation. Internal use only.
#'
#' @param x input data matrix
#' @param base logarithmic base, 2 or 10
#' @param off offset to avoid Infinity
#' @author Panagiotis Moulos
nat2log <- function(x,base=2,off=1) {
	#x[x==0] <- off
	x <- x + off
	if (base==2)
		return(log2(x))
	else
		return(log10(x))
}
