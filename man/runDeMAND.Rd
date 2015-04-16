\name{runDeMAND}
\alias{runDeMAND}
\title{Run DeMAND algorithm}
\usage{
runDeMAND(dobj, fgIndex, bgIndex)
}
\arguments{
\item{dobj}{An instance of class demand which includes: a gene expression data, annotation information, and a molecular interaction network.}
\item{fgIndex}{A numeric vector contains indices of columns which represent case samples (e.g. drug treated). The sample size should be greater than 3.}
\item{bgIndex}{A numeric vector contains indices of columns which represent control samples (e.g. drug treated). The sample size should be greater than 3.}
}
\value{
The moa slot in a DeMAND object will be updated once DeMAND run completes. The prediction is a matrix containing a list of genes, corresponding p-value, and adjusted p-value.
}
\description{
DeMAND.

This function is based on the realization that drugs affect the protein activity of their targets, but not necessarily their mRNA expression levels. In contrast, the change in protein activity directly affects the mRNA expression levels of downstream genes. Based on this hypothesis, DeMAND identifies drug MoA by comparing gene expression profiles following drug perturbation with control samples, and computing the change in the individual interactions within a pre-determined integrated transcriptional and post-translational regulatory model (interactome).
}
\details{
For each edge in the interactome we determine the two-dimensional probability distribution of the gene expression levels both in the control state, and following drug treatment. Any changes in the probability distribution are estimated using the Kullback-Leibler (KL) divergence, from which we determine the statistical significance of the dysregulation of each edge. In the second step of DeMAND, we interrogate each gene independently to determine whether its interactions are enriched in dysregulated ones, suggesting that it is a candidate mechanism of action.
}
\examples{
## Load toy example
data(inputExample)
dobj <- demandClass(exp=bcellExp, anno=bcellAnno, network=bcellNetwork, moa=NULL)
dobj <- runDeMAND(dobj, fgIndex=caseIndex, bgIndex=controlIndex)
## results (head)
printDeMAND(dobj)
## results (all)
print(dobj@moa)
}
\keyword{Action}
\keyword{Drug}
\keyword{Mechanism}
\keyword{of}

