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
  This package runs an external Java Virtual Machine, i.e. it
  circumvents memory limitations of R.
}
\details{
\tabular{ll}{
Package: \tab rMcKmeans\cr
Type: \tab Package\cr
Version: \tab 0.42\cr
Date: \tab 2009-11-12\cr
License: \tab Artistic license 2.0\cr
LazyLoad: \tab yes\cr
}
  K-means clustering is invoked with \code{mckmeans(data, k)}. Cluster number
  estimation is used with \code{cne(data, nruns, k.max)}.
}
\author{
Maintainer: Johann M. Kraus <johann.kraus@uni-ulm.de>
}
\keyword{ package }
\seealso{
\code{\link{mckmeans}}, \code{\link{cne}}}

\examples{
  library(rMcKmeans)
  x <- matrix(runif(50), nc=2)
  \dontrun{
  mckmeans(x, 2)
  cne(x, 10, 10, plot=T)
}
}
