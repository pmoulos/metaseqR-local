\name{diagplot.avg.ftd}
\alias{diagplot.avg.ftd}
\title{Create average False (or True) Discovery 
    curves}
\usage{
    diagplot.avg.ftd(ftdr.obj, output = "x11",
        path = NULL, draw = TRUE, ...)
}
\arguments{
    \item{ftdr.obj}{a list with outputs from 
    \code{\link{diagplot.ftd}}.}

    \item{output}{one or more R plotting 
    device to direct the plot result to.
    Supported mechanisms: \code{"x11"} (default), 
    \code{"png"}, \code{"jpg"}, \code{"bmp"}, 
    \code{"pdf"} or \code{"ps"}.}

    \item{path}{the path to create output files.}

    \item{draw}{boolean to determine whether 
    to plot the curves or just return the 
    calculated values (in cases where the user 
    wants the output for later averaging
    for example). Defaults to \code{TRUE} (make 
    plots).}

    \item{...}{further arguments to be passed to 
    plot devices, such as parameter from 
    \code{\link{par}}.}
}
\value{
    A named list with two members: the first member 
    (\code{avg.ftdr}) contains a list with the 
    means and the standard deviations of the averaged
    \code{ftdr.obj} and are used to create the plot. 
    The second member (\code{path}) contains the 
    path to the created figure graphic.
}
\description{
    This function creates false (or true) discovery 
    curves using a list containing several outputs 
    from \code{\link{diagplot.ftd}}.
}
\examples{
p11 <- 0.001*matrix(runif(300),100,3)
p12 <- matrix(runif(300),100,3)
p21 <- 0.001*matrix(runif(300),100,3)
p22 <- matrix(runif(300),100,3)
p31 <- 0.001*matrix(runif(300),100,3)
p32 <- matrix(runif(300),100,3)
p1 <- rbind(p11,p21)
p2 <- rbind(p12,p22)
p3 <- rbind(p31,p32)
rownames(p1) <- rownames(p2) <- rownames(p3) <-
    paste("gene",1:200,sep="_")
colnames(p1) <- colnames(p2) <- colnames(p3) <-
    paste("method",1:3,sep="_")
truth <- c(rep(1,40),rep(-1,40),rep(0,20),
    rep(1,10),rep(2,10),rep(0,80))
names(truth) <- rownames(p1)
ftd.obj.1 <- diagplot.ftd(truth,p1,N=100,draw=FALSE)
ftd.obj.2 <- diagplot.ftd(truth,p2,N=100,draw=FALSE)
ftd.obj.3 <- diagplot.ftd(truth,p3,N=100,draw=FALSE)
ftd.obj <- list(ftd.obj.1,ftd.obj.2,ftd.obj.3)
avg.ftd.obj <- diagplot.avg.ftd(ftd.obj)
}
\author{
    Panagiotis Moulos
}

