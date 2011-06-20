cne <- function(x, nruns=10, k.min = 2, k.max=10, iter.max=10, nstart=10, plot=F, Xmx="512m", snp=F, infile=NULL){
  if(k.min>=k.max)
    stop("k.min must be smaller than k.max")
  datmode <- is.null(infile)
  if(datmode){
    x <- as.matrix(x)
    if(snp & any(x!=0 | x!=1 | x!=2))
      stop("SNP files have to be encoded as 0,1,2 for 'homozygous reference', 'heterozygous', 'homozygous alternative'")
    # write x to file
    if(snp)
      infile <- ".rmckmeans_infile.snp"
    else
      infile <- ".rmckmeans_infile.tmp"
    write.table(x, infile, row.names=F, col.names=F, quote=F, sep=",")
  }
  outfile <- ".rmckmeans_outfile.tmp"
  cneoutfile <- ".rmckmeans_cneoutfile.tmp"
  # run McKmeans
  system(paste("java -Xmx", Xmx, " -jar ", .mckmeansjar, " -i ", infile, " -o ", outfile, " --maxiter ", iter.max, " --cne", " --cnemin ", k.min, " --cnemax ", k.max, " --cneruns ", nruns, " --cnenstart ", nstart, " --cneoutfile ", cneoutfile, sep=""))
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
  if(datmode){
    if(snp)
      cent <- t(sapply((1:k)-1, function(u) apply(x[res==u,,drop=F], 2, function(v) which.max(table(v))-1)))
    else
      cent <- t(sapply((1:k)-1, function(u) {idx<-which(res==u);colSums(x[idx,,drop=F])/length(idx)}))
    names(cent) <- NULL
  }
  else
    cent <- NA
  mca.cluster <- mean(res.cne[((k-1)*2)-1,-1])
  mca.base <- mean(res.cne[(k-1)*2,-1])
  if(plot)
    cneplot(res.cne)
  rownames(res.cne) <- sapply(2:k.max,function(i) c(paste("McKmeans k=",i,sep=""),paste("Random k=",i,sep="")))
  res.cne <- res.cne[,-1]
  list(centers=cent, cluster=res, mca.cluster=mca.cluster, mca.base=mca.base, mca.all=res.cne)
}

cneplot <- function(mcamatrix){
  ks <- 2:max(mcamatrix[,1])
  tmp <- unlist(apply(mcamatrix[,-1], 1, list), recursive=F)
  boxplot(tmp, ylab="MCA index", names=unlist(lapply(ks, function(u) c(u,u))), col=unlist(lapply(ks,function(u) c(2,4))))
}
