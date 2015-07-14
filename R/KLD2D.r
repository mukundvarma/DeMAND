KLD2D <- function(edge, fgIndex, bgIndex, expData, method=c("integers","bandwidth")[2]){
  #' Calculates the Kullback Leibler divergence (KLD) for the 2D gene expression of X and Y in two conditions 
  #' @param edge a vector holding the names of the two genes for which KLD will be calculated
  #' @param fgIndex a vector of length M holding the indices of the samples in X and Y in condition 1
  #' @param bgIndex a vector of length M holding the indices of the samples in X and Y in condition 2
  #' @param expData a numeric matrix holding the expression data
  #' @docType methods
  #' @return The KL divergence value
  
  
  require(KernSmooth)
  x <- rank(expData[edge[1],c(fgIndex,bgIndex)], ties.method="average")
  y <- rank(expData[edge[2],c(fgIndex,bgIndex)], ties.method="average")
  N <- length(x)
  
  # to calculate the Gaussian kernel smoothing we first need to estimate the bandwidth
  fgWidth <- c(bw.nrd(x[fgIndex]),bw.nrd(y[fgIndex]))
  bgWidth <- c(bw.nrd(x[bgIndex]),bw.nrd(y[bgIndex]))

  gridSize <- switch(method,
                     integers=c(N,N),
                     bandwidth=ceiling(N/c(min(fgWidth[1],bgWidth[1]),min(fgWidth[2],bgWidth[2]))))
    
  ranges <- list(x=c(1,N),y=c(1,N))
  
  fgSmooth <- bkde2D(x=cbind(x[fgIndex],y[fgIndex]),bandwidth=fgWidth,range.x=ranges,gridsize=gridSize)$fhat
  bgSmooth <- bkde2D(x=cbind(x[bgIndex],y[bgIndex]),bandwidth=bgWidth,range.x=ranges,gridsize=gridSize)$fhat
  
  fgSmooth[fgSmooth==0] <- min(fgSmooth[fgSmooth>0])/100
  bgSmooth[bgSmooth==0] <- min(bgSmooth[bgSmooth>0])/100
  
  fgSmooth <- fgSmooth/sum(fgSmooth)
  bgSmooth <- bgSmooth/sum(bgSmooth)
  
  return(sum(fgSmooth*log2(fgSmooth/bgSmooth)) + sum(bgSmooth*log2(bgSmooth/fgSmooth))/2)
}