\name{mckmeans}
\alias{mckmeans}
\title{
  Multi-core K-means clustering
}
\description{
  Perform k-means clustering on a data matrix.
}
\usage{
mckmeans(x, k, iter.max, Xmx, snp)
}
\arguments{
  \item{x}{
    A numeric matrix of data, or an object that can be coerced to
    such a matrix (such as a numeric vector or a data frame with
    all numeric columns).
  }
  \item{k}{
    The number of clusters. A random set of k (distinct)
    rows in `x' is chosen as the initial centres.
  }
  \item{iter.max}{
    The maximum number of iterations allowed.
  }
  \item{Xmx}{
    This sets the maximum memory heap size for Java to Xmx.
  }
    \item{snp}{
    Set this true, if clustering SNP data (categorical data encoded with 0,1,2).
  }
}
\details{
  This function is a wrapper to call the multi-core implementation of
  Kmeans written in Clojure and compiled to Java bytecode.
}
\value{
  cluster: A vector of integers indicating the cluster to which each
           point is allocated.
  
  centers: A matrix of cluster centres.
}
\references{
%% ~put references to the literature/web site here ~
}
\seealso{
  \code{\link{cne}}
}
\examples{
  # a 2-dimensional example
  x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
             matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
  colnames(x) <- c("x", "y")
  # uncomment next 3 lines to run example
  # res <- mckmeans(x, 2)
  # plot(x, col = res$cluster)
  # points(cl$centers, col = 1:2, pch = 8, cex=2)
}
\keyword{ cluster }
\keyword{ parallel }
