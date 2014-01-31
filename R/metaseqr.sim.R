#' Estimate AUFC weights
#'
#' This functions automatically estimates weights for the \code{"weight"} and
#' \code{"dperm.weight"} options of metaseqR for combining p-values from multiple
#' statistical tests. It creates simulated dataset based on real data and then
#' performs statistical analysis with metaseqR several times in order to derive
#' False Discovery Curves. Then, the average areas under the false discovery curves
#' are used to construct weights for each algorithm, according to its performance
#' when using simulated data.
#'
#' @param counts the real raw counts table from which the simulation parameters
#' will be estimated. It must not be normalized and must contain only integer
#' counts, without any other annotation elements and unique gene identifiers as
#' the rownames attribute.
#' @param normalization same as \code{normalization} in \code{link{metaseqr}}.
#' @param statistics same as \code{statistics} in \code{link{metaseqr}}.
#' @param nsim the number of simulations to perform to estimate the weights. It
#' default to 10.
#' @param N the number of genes to produce. See \code{link{make.sim.data.sd}}.
#' @param samples a vector with 2 integers, which are the number of samples for
#' each condition (two conditions currently supported).
#' @param ndeg a vector with 2 integers, which are the number of differentially
#' expressed genes to be produced. The first element is the number of up-regulated
#' genes while the second is the number of down-regulated genes.
#' @param fc.basis the minimum fold-change for deregulation.
#' @param top the top \code{top} best ranked (according to p-value) to use, to
#' calculate area under the false discovery curve.
#' @param model.org the organism from which the data are derived. It must be one
#' of \code{\link{metaseqr}} supported organisms.
#' @param seed a list of seed for reproducible simulations. Defaults to \code{NULL}.
#' @param draw.fpc draw the averaged false discovery curves? Default to \code{FALSE}.
#' @param multic whether to run in parallel (if package \code{parallel} is present
#' or not.
#' @value A vector of weights to be used in \link{\code{metaseqr}} with the
#' \code{weights} option.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' # Not yet available
#'}
estimate.aufc.weights <- function(counts,normalization,statistics,nsim=10,
	N=10000,samples=c(3,3),ndeg=c(500,500),top=500,model.org="mm9",seed=NULL,
	draw.fpc=FALSE,multic=FALSE) {
	if (!require(zoo))
		stopwrap("R pacakage zoo is required in order to estimate AUFC weights!")
	if (ncol(counts)<4)
		stopwrap("Cannot estimate AUFC weights with an initial dataset with less ",
			"than 4 samples!")
	else if (ncol(counts)>=4 && ncol(counts)<10) {
		reind <- sample(20,1:ncol(counts),replace=TRUE)
		counts <- counts[,reind]
	}
	if (is.null(seed)) {
		seed.start <- round(100*runif(1))
		seed.end <- seed.start + nsim - 1
		seed <- as.list(seed.start:seed.end)
	}
	par.list <- estimate.sim.params(counts)

	disp("Running simulations... This procedure requires time... Please wait...")
	sim.results <- wapply(seed,function(x,normalization,statistics,N,par.list,
		samples,ndeg,fc.basis,model.org) {
		D <- make.sim.data.sd(N=N,param=par.list,samples=samples,ndeg=ndeg,
			fc.basis=fc.basis,model.org=model.org,seed=x)
		dd <- D$simdata
		
		tmp <- metaseqr(
			counts=dd,
			sample.list=list(G1=paste("G1_rep",1:samples[1],sep=""),
				G2=paste("G2_rep",1:samples[2],sep="")),
			contrast=c("G1_vs_G2"),
			annotation="embedded",
			id.col=4,
			gc.col=5,
			name.col=7,
			bt.col=8,
			count.type="gene",
			normalization=normalization,
			statistics=statistics,
			meta.p="simes",
			fig.format="png",
			preset="all.basic",
			export.where=tempdir(),
			qc.plots=NULL,
			report=FALSE,
			run.log=FALSE,
			out.list=TRUE
		)

		# Retrieve several p-values
		p.list <- vector("list",length(statistics))
		for (s in statistics) {
			field <- paste("p-value",s,sep="_")
			p[[s]] <- tmp$data[[1]][,field]
			names(p[[s]]) <- rownames(tmp$data[[1]])
		}
		p.matrix <- do.call("cbind",p.list)
		return(list(simdata=D,pvalues=p.matrix))
	},normalization,statistics,N,par.list,samples,ndeg,fc.basis,model.org)

	disp("Estimating AUFC weights... Please wait...")
	fpc.obj <- wapply(sim.results,function(x) {
		true.de <- x$simdata$truedeg
		names(true.de) <- rownames(x$simdata$simdata)
		p.matrix <- x$p.matrix
		true.de <- true.de[rownames(p.matrix)]
		fdc <- diagplot.ftd(true.de,p.matrix,type="fpc",draw=FALSE)
	})
	avg.fpc <- diagplot.avg.ftd(fpc.obj,draw=draw.fpc)

	x <- 1:top
	aufc <- apply(avg.fpc$avg.ftdr$means[1:top,],2,function(x,i) {
		return(sum(diff(i)*rollmean(x,2)))
	},x)
	weight.aufc <- (sum(aufc)/aufc)/sum(sum(aufc)/aufc)
	return(weight.aufc)
}

#' Create simulated counts using TCC package
#'
#' This function creates simulated RNA-Seq gene expression datasets using the
#' \code{simulateReadCounts} function from the Bioconductor
#' package TCC and it adds simulated annoation elements. For further information
#' please consult the TCC package documentation. Note that the produced data are
#' based in an Arabidopsis dataset.
#'
#' @param ... parameters to the \code{simulateReadCounts} function.
#' @return A list with the following members: \code{simdata} holding the simulated
#' dataset complying with metaseqr requirements, and \code{simparam} holding the
#' simulation parameters (see TCC documentation).
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' dd <- make.sim.data(Ngene=10000,PDEG=0.2,DEG.assign=c(0.9,0.1),
#'   DEG.foldchange=c(5,5),replicates=c(3,3))
#' head(dd$simdata)
#'}
make.sim.data.tcc <- function(...) {
	if (suppressWarnings(!require(TCC)))
		stopwrap("Bioconductor package TCC is required to create simulated data!")
	#tcc <- simulateReadCounts(Ngene=Ngene,PDEG=PDEG,DEG.assign=DEG.assign,
	#	DEG.foldchange=DEG.foldchange,replicates=replicates)
	tcc <- simulateReadCounts(...)
	n <- nrow(tcc$count)
	# Now we have to simulate annotation
	chromosome <- paste("chr",1+round(20*runif(n)),sep="")
	start <- 1 + round(1e+6*runif(n))
	end <- start + 250 + round(1e+6*runif(n))
	gene_id <- gene_name <- rownames(tcc$count)
	gc_content <- runif(n)
	strand <- sample(c("+","-"),n,replace=TRUE)
	biotype <- sample(paste("biotype",1:10),n,replace=TRUE)
	sim.data <- data.frame(
		chromosome=chromosome,
		start=start,
		end=end,
		gene_id=gene_id,
		gc_content=gc_content,
		strand=strand,
		gene_name=gene_name,
		biotype=biotype
	)
	sim.data <- cbind(sim.data,tcc$count)
	return(list(simdata=sim.data,simparam=tcc$simulation))
}

#' Create simulated counts using the Soneson-Delorenzi method
#'
#' This function creates simulated RNA-Seq gene expression datasets using the
#' method presented in (Soneson and Delorenzi, BMC Bioinformatics, 2013). For the
#' time being, it creates only simulated datasets with two conditions.
#'
#' @param N the number of genes to produce.
#' @param param a named list with negative binomial parameter sets to sample from.
#' The first member is the mean parameter to sample from (\code{mu.hat}} and the
#' second the dispersion (\code{phi.hat}). This list can be created with the
#' \code{\link{estimate.sim.params}} function.
#' @param samples a vector with 2 integers, which are the number of samples for
#' each condition (two conditions currently supported).
#' @param ndeg a vector with 2 integers, which are the number of differentially
#' expressed genes to be produced. The first element is the number of up-regulated
#' genes while the second is the number of down-regulated genes.
#' @param fc.basis the minimum fold-change for deregulation.
#' @param libsize.range a vector with 2 numbers (generally small, see the default),
#' as they are multiplied with \code{libsize.mag}.
#' These numbers control the library sized of the synthetic data to be produced.
#' @param libsize.mag a (big) number to multiply the \code{libsize.range} to
#' produce library sizes.
#' @param model.org the organism from which the real data are derived from. It
#' must be one of the supported organisms (see the main \code{\link{metaseqr}}
#' help page). It is used to sample real values for GC content.
#' @param seed a seed to use with random number generation for reproducibility.
#' @return A named list with two members. The first member (\code{simdata})
#' contains the synthetic dataset 
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' # File "bottomly_read_counts.txt" from the ReCount database
#' N <- 10000
#' par.list <- estimate.sim.params("bottomly_read_counts.txt")
#' sim <- make.sim.data.sd(N,par.list)
#' synth.data <- sim$simdata
#' true.deg <- which(sim$truedeg!=0)
#'}
make.sim.data.sd <- function(N,param,samples=c(5,5),ndeg=rep(round(0.1*N),2),
	fc.basis=1.5,libsize.range=c(0.7,1.4),libsize.mag=1e+7,model.org=NULL,
	seed=NULL) {
	if (!is.null(model.org)) {
		check.text.args("model.org",model.org,c("hg18","hg19","mm9","mm10",
			"rno5","dm3","danRer7"),multiarg=FALSE)
		ann <- get.annotation(model.org,"gene")
		real.gc <- as.numeric(ann$gc_content)
	}
	mu.hat <- param$mu.hat
	phi.hat <- param$phi.hat
	if (!is.null(seed)) set.seed(seed)
	ii <- sample(1:length(mu.hat),N,replace=TRUE)
	s1 <- samples[1]
	s2 <- samples[2]
	if (!is.null(seed)) set.seed(seed)
	L1 <- round(libsize.mag*runif(s1,min=libsize.range[1],max=libsize.range[2]))
	if (!is.null(seed)) set.seed(2*seed)
	L2 <- round(libsize.mag*runif(s2,min=libsize.range[1],max=libsize.range[2]))

	lambda.1 <- do.call("cbind",rep(list(mu.hat[ii]),s1))
	mu.1 <- sweep(lambda.1,2,L1/sum(lambda.1[,1]),"*")
	sim.1 <- matrix(0,N,s1)
	for (j in 1:s1) {
		if (!is.null(seed)) set.seed(seed+j)
		sim.1[,j] <- rnbinom(N,size=1/phi.hat[ii],mu=mu.1[,j])
	}

	v <- numeric(N)
	if (sum(ndeg)>0) {
		if (!is.null(seed)) set.seed(seed)
		i.updown <- sample(1:length(v),sum(ndeg))
		reg.dir <- rep(c(1,-1),c(ndeg[1],ndeg[2]))
		v[i.updown] <- reg.dir
		if (!is.null(seed)) set.seed(seed+19051980)
		lambda.2 <- ((fc.basis + rexp(N))^v)*lambda.1
		mu.2 <- sweep(lambda.2,2,L2/sum(lambda.2[,1]),"*")
		sim.2 <- matrix(0,N,s2)
		for (j in 1:s2)
			sim.2[,j] <- rnbinom(N,size=1/phi.hat[ii],mu=mu.2[,j])
	}
	else {
		if (!is.null(seed)) set.seed(seed+19051980)
		lambda.2 <- lambda.1
		mu.2 <- sweep(lambda.2,2,L2/sum(lambda.2[,1]),"*")
		sim.2 <- matrix(0,N,s2)
		for (j in 1:s2)
			sim.2[,j] <- rnbinom(N,size=1/phi.hat[ii],mu=mu.2[,j])
	}

	# Now we have to simulate annotation
	if (!is.null(seed)) set.seed(seed)
	chromosome <- paste("chr",1+round(20*runif(N)),sep="")
	if (!is.null(seed)) set.seed(seed)
	start <- 1 + round(1e+6*runif(N))
	if (!is.null(seed)) set.seed(seed)
	end <- start + 250 + round(1e+6*runif(N))
	gene_id <- gene_name <- paste("gene",1:N,sep="_")
	if (!is.null(seed)) set.seed(seed)
	if (!is.null(model.org)) {
		if (length(real.gc)<=N)
			gc_content <- sample(real.gc,N)
		else
			gc_content <- sample(real.gc,N,replace=TRUE)
	}
	else
		gc_content <- runif(N)
	if (!is.null(seed)) set.seed(seed)
	strand <- sample(c("+","-"),N,replace=TRUE)
	if (!is.null(seed)) set.seed(seed)
	biotype <- sample(paste("biotype",1:10),N,replace=TRUE)
	sim.data <- data.frame(
		chromosome=chromosome,
		start=start,
		end=end,
		gene_id=gene_id,
		gc_content=gc_content,
		strand=strand,
		gene_name=gene_name,
		biotype=biotype
	)
	colnames(sim.1) <- paste("G1_rep",1:s1,sep="")
	colnames(sim.2) <- paste("G2_rep",1:s2,sep="")
	rownames(sim.1) <- rownames(sim.2) <- names(v) <- gene_id

	return(list(simdata=cbind(sim.data,sim.1,sim.2),truedeg=v))
}

#' Estimate negative binomial parameters from real data
#'
#' This function reads a read counts table containing real RNA-Seq data (preferebly
#' with more than 20 samples so as to get as much accurate as possible estimations)
#' and calculates a population of count means and dispersion parameters which can
#' be used to simulate an RNA-Seq dataset with synthetic genes by drawing from a
#' negative binomial distribution. This function works in the same way as described
#' in (Soneson and Delorenzi, BMC Bioinformatics, 2013) and (Robles et al., BMC
#' Genomics, 2012).
#'
#' @param real.counts a text tab-delimited file with real RNA-Seq data. The file
#' should strictly contain a unique gene name (e.g. Ensembl accession) in the
#' first column and all other columns should contain read counts for each gene.
#' Each column must be named with a unique sample identifier. See examples in the
#' ReCount database \link{http://bowtie-bio.sourceforge.net/recount/}.
#' @param libsize.gt a library size below which samples are excluded from parameter
#' estimation (default: 3000000).
#' @param rowmeans.gt a row means (mean counts over samples for each gene) below
#' which genes are excluded from parameter estimation (default: 5).
#' @param eps the tolerance for the convergence of \code{\link{optimize}} function.
#' Defaults to 1e-11.
#' @param restrict.cores in case of parallel optimization, the fraction of the
#' available cores to use.
#' @param seed a seed to use with random number generation for reproducibility.
#' @return A named list with two members: \code{mu.hat} which contains negative
#' binomial mean estimates and \code{phi.hat} which contains dispersion.
#' estimates
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' par.list <- estimate.sim.params("bottomly_read_counts.txt")
#'}
estimate.sim.params <- function(real.counts,libsize.gt=3e+6,rowmeans.gt=5,
	eps=1e-11,restrict.cores=0.8,seed=42) {
	multic <- check.parallel(restrict.cores)
	if (is.data.frame(real.counts))
		real.data <- real.counts
	else if (file.exists(real.counts))
		real.data <- read.delim(real.counts,row.names=1)
	else
		stopwrap("The input count data must be either a file or a data frame!")
	mat <- as.matrix(real.data)
	low.lib <- which(apply(mat,2,sum)<libsize.gt)
	if (length(low.lib)>0)
		mat <- mat[,-low.lib]
	disp("Downsampling counts...")
	dmat <- downsample.counts(mat,seed)
	low.co <- which(apply(dmat,1,function(x) if (mean(x)<5) TRUE else FALSE))
	if (length(low.co)>0)
		dmat <- dmat[-low.co,]
	mu.hat <- apply(dmat,1,mean)
	disp("Estimating initial dispersion population...")
	phi.est <- apply(dmat,1,function(x) {
		m <- mean(x)
		v <- var(x)
		phi <- (v-m)/m^2
		return(phi)
	})
	phi.ind <- which(phi.est>0)
	phi.est <- phi.est[phi.ind]
	dmat <- dmat[phi.ind,]
	disp("Estimating dispersions using log-likelihood...\n")
	init <- wapply(multic,seq_along(1:nrow(dmat)),function(i,d,p) {
		list(y=d[i,],h=p[i])
	},dmat,phi.est)
	phi.hat <- unlist(wapply(multic,init,function(x,eps) {
		optimize(mlfo,c(x$h-1e-2,x$h+1e-2),y=x$y,tol=eps)$minimum
	},eps))
	return(list(mu.hat=mu.hat[phi.ind],phi.hat=phi.hat))
}

#' Downsample read counts
#'
#' This function downsamples the library sizes of a read counts table to the lowest
#' library size, according to the methdology used in  (Soneson and Delorenzi,
#' BMC Bioinformatics, 2013).
#'
#' @param counts the read counts table which is subjected to downsampling.
#' @param seed random seed for reproducible downsampling.
#' @return The downsampled counts matrix.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' M <- as.matrix(read.delim("bottomly_read_counts.txt",row.names=1))
#' D <- downsample.counts(M)
#'}
downsample.counts <- function(counts,seed=42) {
	lib.sizes <- apply(counts,2,sum)
	target.size <- min(lib.sizes)
	to.remove <- lib.sizes-target.size
	ii <- which(to.remove>0)
	dcounts <- counts
	for (i in ii) {
		tmp <- round(to.remove[i]*(counts[,i]/sum(counts[,i])))
		victim.size <- sum(tmp)
		if (victim.size>to.remove[i]) {
			dif <- victim.size - to.remove[i]
			#victims <- sample(1:length(tmp),dif)
			victims <- sort(tmp,decreasing=TRUE,index.return=TRUE)$ix[1:dif]
			tmp[victims] <- tmp[victims] - 1
		}
		else if (victim.size<to.remove[i]) {
			dif <- to.remove[i] - victim.size
			#victims <- sample(1:length(tmp),dif)
			victims <- sort(tmp,decreasing=TRUE,index.return=TRUE)$ix[1:dif]
			tmp[victims] <- tmp[victims] + 1
		}
		dcounts[,i] <- dcounts[,i] - tmp
	}
	return(dcounts)
}

#' MLE dispersion estimate
#'
#' MLE function used to estimate negative binomial dispersions from real RNA-Seq
#' data, as in (Soneson and Delorenzi, BMC Bioinformatics, 2013) and (Robles et al.,
#' BMC Genomics, 2012). Internal use.
#'
#' @param phi the parameter to be optimized.
#' @param y count samples used to perform the optimization.
#' @return objective function value.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' # Not yet available
#'}
mlfo <- function(phi,y) {
	N <- length(y)
	mu <- mean(y)
	-(sum(lgamma(y+1/phi)) - N*lgamma(1/phi) - sum(lgamma(y+1)) + 
		sum(y*log(mu*phi/(1+mu*phi))) - (N/phi)*log(1+mu*phi))
}

#' Create counts matrix permutations
#'
#' This function creates a permuted read counts matrix based on the \code{contrast}
#' argument (to define new virtual contrasts of the same number) and on the
#' \code{sample.list} to derive the number of samples for each virtual condition.
#' It is a helper for the \code{\link{meta.perm}} function.
#'
#' @param counts the gene read counts matrix.
#' @param sample.list the list containing condition names and the samples under
#' each condition.
#' @param contrast the contrasts vector. See the main \code{\link{metaseqr}} help
#' page.
#' @param repl the same as the replace argument in \code{\link{sample}} function.
#' @return A list with three members: the matrix of permuted per sample read counts,
#' the virtual sample list and the virtual contrast to be used with the \code{stat.*}
#' functions.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' # Not yet available
#'}
make.permutation <- function(counts,sample.list,contrast,repl=FALSE) {
	cnts <- strsplit(contrast,"_vs_")[[1]]
	virtual.contrast <- paste(paste("VirtCond",1:length(cnts),sep=""),
		collapse="_vs_")
	virtual.sample.list <- vector("list",length(sample.list))
	names(virtual.sample.list) <- paste("VirtCond",1:length(sample.list),sep="")
	# Avoid the extreme case of returning a vector with all samples the same
	if (repl) {
		resample <- rep(1,ncol(counts))
		while(length(unique(resample))==1)
			resample <- sample(1:ncol(counts),ncol(counts),replace=repl)
	}
	else
		resample <- sample(1:ncol(counts),ncol(counts),replace=repl)
	virtual.counts <- counts[,resample]
	samples <- paste("VirtSamp",1:ncol(counts),sep="")
	colnames(virtual.counts) <- samples
	nsample <- sapply(sample.list,length)
	virtual.samples <- split(samples,rep(1:length(nsample),nsample))
	names(virtual.samples) <- names(virtual.sample.list)
	for (n in names(virtual.sample.list))
		virtual.sample.list[[n]] <- virtual.samples[[n]]
	return(list(counts=virtual.counts,sample.list=virtual.sample.list,
		contrast=virtual.contrast))
}
