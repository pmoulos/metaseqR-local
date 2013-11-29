#' Meta-analysis using several RNA-Seq statistics
#'
#' This function calculates the combined p-values when multiple statistical algorithms are applied to the input dataset. It is a helper
#' and it requires very specific arguments so it should not be used individually
#'
#' @param cp.list a named list whose names are the contrasts requested from metaseqr. Each member is a p-value matrix whose colnames
#' are the names of the statistical tests applied to the data. See the main \code{\link{metaseqr}} help page.
#' @param meta.p the p-value combination method to use. See the main \code{\link{metaseqr}} help page.
#' @param counts the normalized and possibly filtered read counts matrix. See the main \code{\link{metaseqr}} help page.
#' @param sample.list the list containing condition names and the samples under each condition. See the main \code{\link{metaseqr}}
#' help page.
#' @param statistics the statistical algorithms used in metaseqr. See the main \code{\link{metaseqr}} help page.
#' @param stat.args the parameters for each statistical argument. See the main \code{\link{metaseqr}} help page.
#' @param norm.args the parameters for normalization. See the main \code{\link{metaseqr}} and the \code{stat.*} help pages.
#' @param libsize.list a list with library sizes. See the main \code{\link{metaseqr}} and the \code{stat.*} help pages.
#' @param nperm the number of permutations (Monte Carlo simulations) to perform.
#' @param weight a numeric vector of weights for each statistical algorithm.
#' @param reprod create reproducible permutations. Ideally one would want to create the same set of indices for a given dataset so
#' as to create reproducible p-values. If \code{reprod=TRUE}, a fixed seed is used by \code{meta.perm} for all the datasets analyzed
#' with \code{metaseqr}. If \code{reprod=FALSE}, then the p-values will not be reproducible, although statistical significance is not
#' expected to change for a large number of resambling. Finally, \code{reprod} can be a numeric vector of seeds with the same length as
#' \code{nperm} so that the user can supply his/her own seeds.
#' @param multic use multiple cores to execute the premutations. This is an external parameter and implies the existence of parallel
#' package in the execution environment. See the main \code{\link{metaseqr}} help page.
#' @return A named list with combined p-values. The names are the contrasts and the list members are combined p-value vectors, one
#' for each contrast.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' # Not yet available
#'}
meta.test <- function(cp.list,meta.p=c("simes","fisher","dperm.min","dperm.max","dperm.weight","fperm","whitlock","intersection","union","none"),
	counts,sample.list,statistics,stat.args,norm.args,libsize.list,nperm=10000,weight=rep(1/length(statistics),length(statistics)),
	reprod=TRUE,multic=FALSE) {
	check.text.args("meta.p",meta.p,c("simes","fisher","dperm.min","dperm.max","dperm.weight","fperm","whitlock","intersection","union","none"))
	contrast <- names(cp.list)
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
		fperm = {
			sum.p.list <- wapply(multic,cp.list,function(x) {
				if (multic)
					tmp <- fisher.method.perm(x,p.corr="none",B=nperm,mc.cores=getOption(cores),zero.sub=1e-32)
				else
					tmp <- fisher.method.perm(x,p.corr="none",B=nperm,zero.sub=1e-32)
				return(tmp$p.value)
			})
		},
		whitlock = {
			sum.p.list <- wapply(multic,cp.list,function(x) return(apply(x,1,combine.test,method="z.transform")))
		},
		#hommel = { # Returns a matrix of p-values, not summary
		#	sum.p.list <- wapply(multic,cp.list,function(x) return(apply(x,1,p.adjust,"hommel")))
		#},
		simes = {
			sum.p.list <- wapply(multic,cp.list,function(x) {
				return(apply(x,1,function(p,m) {
					s <- min(m*p)/length(p)
					if (s>1 || is.na(s)) return(1) else return(s)
				},nrow(x)))
			})
		},
		dperm.min = {
			sum.p.list <- vector("list",length(cp.list))
			names(sum.p.list) <- names(cp.list)
			conl <- as.list(contrast)
			names(conl) <- contrast
			temp.p.list <- wapply(multic,conl,meta.perm,
				counts=counts,sample.list=sample.list,
				statistics=statistics,stat.args=stat.args,
				norm.args=norm.args,libsize.list=libsize.list,
				nperm=nperm,weight=weight,
				select="min",reprod=reprod,
				multic=multic)
			original.p.list <- wapply(multic,cp.list,function(x,m,w=NULL) {
				switch(m,
					min = {
						return(apply(x,1,min,na.rm=TRUE))
					},
					max = {
						return(apply(x,1,max,na.rm=TRUE))
					},
					weight = {
						return(apply(x,1,function(p,w) {
							if (any(is.na(p)))
								return(1)
							else
								return(prod(w,p))
						},w))
					}
				)
			},"min")
			for (cc in names(original.p.list))
			{
				pc <- cbind(temp.p.list[[cc]],original.p.list[[cc]])
				ly <- ncol(pc)
				sum.p.list[[cc]] <- apply(pc,1,function(y,m) return(length(which(y[1:(m-1)]<y[m]))/(m-1)),ly)
			}
		},
		dperm.max = {
			sum.p.list <- vector("list",length(cp.list))
			names(sum.p.list) <- names(cp.list)
			conl <- as.list(contrast)
			names(conl) <- contrast
			temp.p.list <- wapply(multic,conl,meta.perm,
				counts=counts,sample.list=sample.list,
				statistics=statistics,stat.args=stat.args,
				norm.args=norm.args,libsize.list=libsize.list,
				nperm=nperm,weight=weight,
				select="max",reprod=reprod,
				multic=multic)
			original.p.list <- wapply(multic,cp.list,function(x,m,w=NULL) {
				switch(m,
					min = {
						return(apply(x,1,min,na.rm=TRUE))
					},
					max = {
						return(apply(x,1,max))
					},
					weight = {
						return(apply(x,1,function(p,w) {
							if (any(is.na(p)))
								return(1)
							else
								return(prod(w,p))
						},w))
					}
				)
			},"max")
			for (cc in names(original.p.list))
			{
				pc <- cbind(temp.p.list[[cc]],original.p.list[[cc]])
				ly <- ncol(pc)
				sum.p.list[[cc]] <- apply(pc,1,function(y,m) return(length(which(y[1:(m-1)]<y[m]))/(m-1)),ly)
			}
		},
		dperm.weight = {
			sum.p.list <- vector("list",length(cp.list))
			names(sum.p.list) <- names(cp.list)
			conl <- as.list(contrast)
			names(conl) <- contrast
			temp.p.list <- wapply(multic,conl,meta.perm,
				counts=counts,sample.list=sample.list,
				statistics=statistics,stat.args=stat.args,
				norm.args=norm.args,libsize.list=libsize.list,
				nperm=nperm,weight=weight,
				select="min",reprod=reprod,
				multic=multic)
			original.p.list <- wapply(multic,cp.list,function(x,m,w=NULL) {
				switch(m,
					min = {
						return(apply(x,1,min))
					},
					max = {
						return(apply(x,1,max))
					},
					weight = {
						return(apply(x,1,function(p,w) {
							if (any(is.na(p)))
								return(1)
							else
								return(prod(w,p))
						},w))
					}
				)
			},"weight",weight)
			for (cc in names(original.p.list))
			{
				pc <- cbind(temp.p.list[[cc]],original.p.list[[cc]])
				ly <- ncol(pc)
				sum.p.list[[cc]] <- apply(pc,1,function(y,m) return(length(which(y[1:(m-1)]<y[m]))/(m-1)),ly)
			}
		},
		none = { # A default value must be there to use with volcanos, we say the one of the first statistic in order of input
			sum.p.list <- wapply(multic,cp.list,function(x) return(x[,1]))
		}
	)
	return(sum.p.list)
}

#' Permutation tests for meta-analysis
#'
#' This function performs permutation tests in order to derive a meta p-value by combining several of the statistical algorithms of
#' metaseqr. This is probably the most accurate way of combining multiple statistical algorithms for RNA-Seq data, as this issue is
#' different from the classic interpretation of the term "meta-analysis" which implies the application of the same statistical test
#' on different datasets treating the same subject/experiment. For other methods, see also the main \code{\link{metaseqr}} help page.
#' You should keep in mind that the permutation procedure can take a long time, even when executed in parallel.
#'
#' @param counts a normalized read counts table, one row for each gene, one column for each sample
#' @param sample.list the list containing condition names and the samples under each condition. See the main \code{\link{metaseqr}}
#' help page.
#' @param contrast the contrasts to be tested by each statistical algorithm. See the main \code{\link{metaseqr}} help page.
#' @param statistics the statistical algorithms used in metaseqr. See the main \code{\link{metaseqr}} help page.
#' @param stat.args the parameters for each statistical algorithm. See the main \code{\link{metaseqr}} help page.
#' @param norm.args the parameters for normalization. See the main \code{\link{metaseqr}} and the \code{stat.*} help pages.
#' @param libsize.list a list with library sizes. See the main \code{\link{metaseqr}} and the \code{stat.*} help pages.
#' @param nperm the number of permutations (Monte Carlo simulations) to perform.
#' @param weight a numeric vector of weights for each statistical algorithm
#' @param select how to select the initial vector of p-values. It can be \code{"min"} to select the minimum p-value for each gene
#' (more conservative), \code{"max"} to select the maximum p-value for each gene (less conservative), \code{"weight"} to apply the
#' weights to the p-value vector for each gene and derive a weighted p-value.
#' @param replace same as the \code{replace} argument in the \code{\link{sample}} function. Implies bootstraping or simple resampling
#' without replacement. It can also be \code{"auto"}, to determine bootstraping or not with the following rule: if \code{ncol(counts)<=4}
#' \code{replace=FALSE else} \code{replace=TRUE}. This protects from the case of having zero variability across resampled conditions.
#' In such cases, most statistical tests would crash.
#' @param multic use multiple cores to execute the premutations. This is an external parameter and implies the existence of parallel
#' package in the execution environment. See the main \code{\link{metaseqr}} help page.
#' @param reprod create reproducible permutations. Ideally one would want to create the same set of indices for a given dataset so
#' as to create reproducible p-values. If \code{reprod=TRUE}, a fixed seed is used by \code{meta.perm} for all the datasets analyzed
#' with \code{metaseqr}. If \code{reprod=FALSE}, then the p-values will not be reproducible, although statistical significance is not
#' expected to change for a large number of resambling. Finally, \code{reprod} can be a numeric vector of seeds with the same length as
#' \code{nperm} so that the user can supply his/her own seeds.
#' @return A vector of meta p-values
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' # Not yet available
#'}
meta.perm <- function(contrast,counts,sample.list,statistics,stat.args,norm.args,libsize.list,nperm=10000,
	weight=rep(1/ncol(counts),ncol(counts)),select=c("min","max","weight"),replace="auto",reprod=TRUE,
	multic=FALSE) {
	check.text.args("select",select,c("min","max","weight"))
	if (replace=="auto") {
		if (ncol(counts)<=4)
			replace=FALSE
		else
			replace=TRUE
	}
	# We will construct relist in a way so that we can assign seeds for random number generation and track progress at the same time
	if (is.logical(reprod)) {
		relist <- vector("list",nperm)
		if (reprod) {
			relist <- wapply(multic,seq_along(relist),function(i) {return(list(seed=i,prog=i))})
		}
		else
			relist <- wapply(multic,seq_along(relist),function(i) {return(list(seed=round(1e+6*runif(1)),prog=i))})
	}
	else if (is.numeric(reprod)) {
		if (length(reprod) != nperm)
			stopwrap("When reprod is numeric, it must have the same length as nperm!")
		relist <- wapply(multic,seq_along(reprod),function(i) {return(list(seed=reprod[i],prog=i))})
	}
	else
		stopwrap("reprod must be either a logical or a numeric vector!")
	disp("  Resampling procedure started...")
	# In this case, we must not use wapply as we want to be able to track progress through mc.preschedule...
	if (multic)
		pp <- mclapply(relist,meta.worker,counts,sample.list,contrast,statistics,replace,
			stat.args,norm.args,libsize.list,select,mc.preschedule=FALSE,mc.cores=getOption("cores"))
	else
		pp <- lapply(relist,meta.worker,counts,sample.list,contrast,statistics,replace,stat.args,norm.args,
			libsize.list,select)
	disp("  Resampling procedure ended...")
	return(do.call("cbind",pp))
}

#' Permutation tests helper
#'
#' This function performs the statistical test for each permutation. Internal use only.
#'
#' @param x a virtual list with the random seed and the permutation index.
#' @param co the counts matrix.
#' @param sl the sample list.
#' @param s the statistical algorithms.
#' @param sa the parameters for each statistical algorithm.
#' @param na the parameters for normalization.
#' @param ll a list with library sizes.
#' @param r same as the \code{replace} argument in the \code{\link{sample}} function.
#' @param el min, max or weight.
#' @return A matrix of p-values.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' # Not yet available
#'}
meta.worker <- function(x,co,sl,cnt,s,r,sa,na,ll,el) {
	set.seed(x$seed)
	disp("    running permutation #",x$prog)
	pl <- make.permutation(co,sl,cnt,r)
	ppmat <- matrix(NA,nrow(co),length(s))
	colnames(ppmat) <- s
	for (alg in s) {
		#disp("      running permutation tests with: ",alg)
		tcl <- make.contrast.list(pl$contrast,pl$sample.list)
		switch(alg,
			deseq = {
				p.list <- suppressMessages(stat.deseq(pl$counts,pl$sample.list,tcl,sa[[alg]],na))
			},
			edger = {
				p.list <- suppressMessages(stat.edger(pl$counts,pl$sample.list,tcl,sa[[alg]]))
			},
			noiseq = {
				p.list <- suppressMessages(stat.noiseq(pl$counts,pl$sample.list,tcl,sa[[alg]],na))
			},
			bayseq = {
				p.list <- suppressMessages(stat.bayseq(pl$counts,pl$sample.list,tcl,sa[[alg]],na,ll))
			},
			limma = {
				p.list <- suppressMessages(stat.limma(pl$counts,pl$sample.list,tcl,sa[[alg]]))
			},
			nbpseq = {
				p.list <- suppressMessages(stat.nbpseq(pl$counts,pl$sample.list,tcl,sa[[alg]],na,ll))
			}
		)
		ppmat[,alg] <- as.numeric(p.list[[1]])
	}
	switch(el,
		min = {
			p.iter <- apply(ppmat,1,min,na.rm=TRUE)
		},
		max = {
			p.iter <- apply(ppmat,1,max,na.rm=TRUE)
		},
		weight = {
			p.iter <- apply(ppmat,1,function(p,w) {
				if (any(is.na(p)))
					return(1)
				else
					return(prod(w,p))
			})
		}
	)
	return(p.iter)
}
