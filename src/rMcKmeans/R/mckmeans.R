
mckmeans <- function(x, k=2, iter.max=10, nstart=1, Xmx="512m", snp=F, infile=NULL){
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
  
  # run McKmeans
  system(paste("java -Xmx", Xmx, " -jar ", .mckmeansjar, " -i ", infile, " -o ", outfile, " -k ", k, " --maxiter ", iter.max, " --nstart ", nstart, sep=""))
  # read results from file
  options(warn=-1)
  res <- as.matrix(read.table(outfile, sep=" ", header=F, quote=""))[1,]
  names(res) <- NULL
  options(warn=0) 
  # sweep tmp data
  file.remove(infile)
  file.remove(outfile)
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
  list(centers=cent, cluster=res)
}

