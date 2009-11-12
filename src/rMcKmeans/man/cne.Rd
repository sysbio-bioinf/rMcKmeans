\name{cne}
\alias{cne}
\title{
  Cluster number estimation
}
\description{
  Perform a cluster number estimation on a data matrix.
}
\usage{
cne(x, nruns, k.max=10, iter.max=10, plot=F)
}
\arguments{
  \item{x}{
    A numeric matrix of data, or an object that can be coerced to
    such a matrix (such as a numeric vector or a data frame with
    all numeric columns).
  }
  \item{nruns}{
    The number of repeated clusterings. For every k from 2 to k.max
    nruns repeated runs of clustering are performed.
  }
  \item{k.max}{
    The maximal number of clusters. For every k from 2 to k.max a
    random set of k (distinct) rows in `x' is chosen as the initial centres.
  }
  \item{iter.max}{
    The maximum number of iterations allowed.
  }
  \item{plot}{
    Plot the results of the cluster number estimation.
  }
}
\details{
  This is a cluster number estimation via stability analysis. For every
  k from 2 to k.max nruns repeated runs of clustering are
  performed. These clusterings are compared pairwise via the MCA cluster
  validation method. The most stable results with respect to a random
  clustering is reported as the best clustering. The boxplots
  show the collected MCA values for every k from clustering and random clustering.
}
\value{
  cluster: A vector of integers indicating the cluster to which each
           point is allocated. This is the best clustering reported by
	   the cluster number estimation process.

  centers: A matrix of cluster centres.
	   
  mca.cluster: The mean MCA index from the comparison of all cluster
               results along the best clustering k.
  
  mca.base: The mean MCA index from the comparison of all random cluster
            results along the best clustering k.
  
  mca.all: All MCA values collected during the cluster number
           estimation. Two rows for each k give the mca values from clustering
           and random clustering. The entry in the first column is the
	   number of clusters k.
}
\references{
%% ~put references to the literature/web site here ~
}
\seealso{
  \code{\link{mckmeans}}
}
\examples{
  # a 2-dimensional example
  x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
             matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
  colnames(x) <- c("x", "y")
  # uncomment next line to run cluster number estimation
  # res <- cne(x, nruns=10, k.max=10, iter.max=100, plot=T)
}
\keyword{ cluster }
