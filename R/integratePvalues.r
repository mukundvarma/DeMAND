integratePvalues <- function(g, network, expData){
  #' This function integrates the p-values of all the edges surrounding a gene using Fisher's method,
  #' and uses Brown's method to correct for correlations between the p-values.
  #' @param g the name of the gene
  #' @param network a matrix or data frame holding the gene names in the first two columns, followed by the KLD value and followed by the p-balue
  #' @param expData gene expression data with gene names as rownames
  #' @docType methods
  #' @return the corrected p-value for gene g

  # find all the interactions involving g
  neighborEdges<- which(rowSums(matrix(data=network[,1:2] %in% g,nrow=dim(network)[1],ncol=2))>0)
  if(length(neighborEdges)<2)
    return(as.numeric(network[neighborEdges,4]))
  
  # calculate the Fisher method integrated chi square value
  chiVal <- -2*sum(log(as.numeric(network[neighborEdges,4])))
  
  # get the expression of the neighbors to calculate correlations
  neighborGenes <- setdiff(network[neighborEdges,1:2],g)
  if(length(neighborGenes)<2)
    return(as.numeric(network[neighborEdges,4][1]))
  neighborExp <- expData[neighborGenes,]
  gExp <- expData[g,]
    
  # since we are trying to correlate esges, we estimate the residuals from a fit to gene g
  resids <- lm(t(neighborExp)~gExp)$residuals  
  covMat <- cov(resids)
  
  # calculate the correction to Fisher's method
  covMat <- covMat[lower.tri(covMat)]
  I <- covMat<0 # Brown's method requires a different correction for positive and negative values
  covMat[I] <- 3.27 + 0.71*covMat[I]
  covMat[!I] <- 3.25 + 0.75*covMat[!I]
  chiVar <- 4*length(neighborEdges) + sum(covMat)
  dFreedom <- 8*length(neighborEdges)^2/chiVar
  
  return(pchisq(q=chiVal,df=dFreedom,lower.tail=F))
}
