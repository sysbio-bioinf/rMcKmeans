\name{mckmeans}
\alias{mckmeans}
\title{
  Multi-core K-means clustering
}
\description{
  Perform k-means clustering on a data matrix using a parallelized
  k-means algorithm.
}
\usage{
mckmeans(x, k, iter.max, nstart, Xmx, snp, infile)
}
\arguments{
  \item{x}{
    A numeric matrix of data, or an object that can be coerced to
    such a matrix (such as a numeric vector or a data frame with
    numeric columns only). Data is clustered row-wise. Can be omitted,
    if \code{infile is specified}.
  }
  \item{k}{
    The number of clusters. A random set of k (distinct)
    rows in \code{x} is chosen as the initial centres.
  }
  \item{iter.max}{
    The maximum number of iterations allowed.
  }
  \item{nstart}{
    The number of random restarts to be taken.
  }
  \item{Xmx}{
    This sets the maximum memory heap size for Java to Xmx.
  }
  \item{snp}{
    Set this to true if clustering SNP data (categorical data encoded
    with 0,1,2 for 'homozygous reference', 'heterozygous', 'homozygous alternative').
  }
  \item{infile}{
    Path to the data file. If not set, x is taken as input data.
  }
}
\details{
  This function is a wrapper that calls the multi-core implementation of
  K-means written in Clojure and compiled to Java bytecode. Data is
  clustered by rows.
}
\value{
  \item{cluster}{ A vector of integers indicating the cluster to which each
           point is assigned.}
  
  \item{centers}{ A matrix of cluster centres.}
}
\seealso{
  \code{\link{cne}}
}
\examples{
  library(rMcKmeans)
  # a 2-dimensional example
  x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
             matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
  colnames(x) <- c("x", "y")
  \dontrun{
  res <- mckmeans(x, 2)
  plot(x, col = res$cluster+1)
  points(res$centers, col = 1:2, pch = 8, cex=2)
  }
}
\keyword{ cluster }
