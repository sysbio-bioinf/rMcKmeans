cne <- function(x, nruns, k.max=10, iter.max=10, plot=F){
  x <- as.matrix(x)
  # write x to file
  infile <- ".rmckmeans_infile.tmp"
  outfile <- ".rmckmeans_outfile.tmp"
  cneoutfile <- ".rmckmeans_cneoutfile.tmp"
  write.table(x, infile, row.names=F, col.names=F, quote=F, sep="\t")
  # run McKmeans
  system(paste("java -Xmx2g -jar", .mckmeansjar, "-i", infile, "-o", outfile, "--maxiter", iter.max, "--cne", "--cnemax", k.max, "--cneruns", nruns, "--cneoutfile", cneoutfile))
  # read results from file
  options(warn=-1)
  res <- as.matrix(read.table(outfile, sep=" ", header=F, quote=""))[1,]
  k <- length(unique(res))
  res.cne <- as.matrix(read.table(cneoutfile, sep=" ", header=F, quote=""))
  names(res) <- NULL
  colnames(res.cne) <- NULL
  options(warn=0) 
  # sweep tmp data
  file.remove(infile)
  file.remove(outfile)
  file.remove(cneoutfile)
  # return result
  cent <- t(sapply((1:k)-1, function(u) {idx<-which(res==u);colSums(x[idx,])/length(idx)}))
  names(cent) <- NULL
  mca.cluster <- mean(res.cne[((k-1)*2)-1,-1])
  mca.base <- mean(res.cne[(k-1)*2,-1])
  if(plot)
    cneplot(res.cne)
  list(centers=cent, cluster=res, mca.cluster=mca.cluster, mca.base=mca.base, mca.all=res.cne)
}

cneplot <- function(mcamatrix){
  ks <- 2:max(mcamatrix[,1])
  tmp <- unlist(apply(mcamatrix[,-1], 1, list), recursive=F)
  boxplot(tmp, ylab="MCA index", names=unlist(lapply(ks, function(u) c(u,u))), col=unlist(lapply(ks,function(u) c(2,4))))
}