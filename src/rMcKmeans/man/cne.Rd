\name{cne}
\alias{cne}
\title{
  Cluster number estimation
}
\description{
  Perform a cluster number estimation on a data matrix using a
  parallelized K-means algorithm.
}
\usage{
cne(x, nruns=10, k.max=10, iter.max=10, plot=F, Xmx="512m", snp=F,
infile=NULL)
}
\arguments{
  \item{x}{
    A numeric matrix of data, or an object that can be coerced to
    such a matrix (such as a numeric vector or a data frame with
    numeric columns only). Data is clustered row-wise. Can be omitted,
    if \code{infile} is given.
  }
  \item{nruns}{
    The number of repeated clusterings. For each k from 2 to \code{k.max},
    \code{nruns} repeated runs of clustering are performed.
  }
  \item{k.max}{
    The maximal number of clusters. For each k from 2 to \code{k.max}, a
    random set of k (distinct) rows in \code{x} is chosen as the initial centres.
  }
  \item{iter.max}{
    The maximum number of iterations allowed.
  }
  \item{plot}{
    Plot the results of the cluster number estimation.
  }
  \item{Xmx}{
    This sets the maximum memory heap size for Java to Xmx.
  }  
  \item{snp}{
    Set this to true, if clustering SNP data (categorical data encoded
    with 0,1,2 for 'homozygous reference', 'heterozygous', 'homozygous alternative').
  }
  \item{infile}{
    Path to the data file. If not set, x is taken as input data.
  }
}
\details{
  This is a cluster number estimation via stability analysis. Data is
  clustered row-wise. For each
  k from 2 to \code{k.max} \code{nruns} repeated runs of clustering are
  performed. Pairwise comparisons of these clusterings are performed via the MCA cluster
  validation method. The most stable result with respect to a random
  clustering is reported as the best clustering. If \code{plot} is true,
  boxplots that show the collected MCA values of the clusterings and the
  random clusterings for each k are drawn.
}
\value{
  \item{cluster}{ A vector of integers indicating the cluster to which each
           point is assigned. This is the best clustering reported by
	   the cluster number estimation process.
  }
  \item{centers}{ A matrix of cluster centres.}
	   
  \item{mca.cluster}{ The mean MCA index of the comparison of all cluster
               results along the best clustering k.}
  
  \item{mca.base}{ The mean MCA index of the comparison of all random
    clusterings along the best clustering k.}
  
  \item{mca.all}{ All MCA values collected during the cluster number
           estimation. For each k, two consecutive rows give the MCA values from clustering
           and random clustering.}
}
\seealso{
  \code{\link{mckmeans}}
}
\examples{
  library(rMcKmeans)
  # a 2-dimensional example
  x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
             matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
  colnames(x) <- c("x", "y")
  \dontrun{
  res <- cne(x, nruns=10, k.max=10, iter.max=100, plot=T)
  }
}
\keyword{ cluster }
