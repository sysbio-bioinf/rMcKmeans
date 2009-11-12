\name{rMcKmeans-package}
\alias{rMcKmeans-package}
\alias{rMcKmeans}
\docType{package}
\title{
Multi-core Kmeans
}
\description{
  Perform a K-means clustering and a cluster number estimation on a data
  matrix. These methods are designed to make use of multi-core hardware.
}
\details{
\tabular{ll}{
Package: \tab rMcKmeans\cr
Type: \tab Package\cr
Version: \tab 0.42\cr
Date: \tab 2009-11-12\cr
License: \tab Artistic license\cr
LazyLoad: \tab yes\cr
}
  K-means clustering is invoked with mckmeans(data, k). Cluster number
  estimation is used with cne(data, nruns, k.max).
}
\author{
Maintainer: Johann M. Kraus <johann.kraus@uni-ulm.de>
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~}

\examples{
  x <- matrix(runif(50), nc=2)
  # uncomment next 2 lines to run examples
  # mckmeans(x, 2)
  # cne(x, 10, 10, plot=T)
}
