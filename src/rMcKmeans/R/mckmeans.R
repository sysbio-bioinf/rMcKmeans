mckmeans <- function(x, k=2, iter.max=10){
  x <- as.matrix(x)
  # write x to file
  infile <- ".rmckmeans_infile.tmp"
  outfile <- ".rmckmeans_outfile.tmp"
  write.table(x, infile, row.names=F, col.names=F, quote=F, sep="\t")
  # run McKmeans
  system(paste("java -Xmx2g -jar", .mckmeansjar, "-i", infile, "-o", outfile, "-k", k, "--maxiter", iter.max))
  # read results from file
  options(warn=-1)
  res <- as.matrix(read.table(outfile, sep=" ", header=F, quote=""))[1,]
  names(res) <- NULL
  options(warn=0) 
  # sweep tmp data
  file.remove(infile)
  file.remove(outfile)
  # return result
  cent <- t(sapply((1:k)-1, function(u) {idx<-which(res==u);colSums(x[idx,])/length(idx)}))
  names(cent) <- NULL
  list(centers=cent, cluster=res)
}

