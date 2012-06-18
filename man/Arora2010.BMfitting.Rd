\name{Arora2010.BMfitting}
\alias{Arora2010.BMfitting}
\title{
Fitting a beta-mixture model to association scores
}                                                                        
\description{
The function fits a three-beta mixture model to functional association scores 
computed from rich phenotyping screens in Wang X et al., 2012 for RNAi screening 
data of Ewing's sarcoma in Arora S et al., 2010.
}
\usage{
Arora2010.BMfitting(pheno, model="global", metric="cosine", nPerm=20)
}
\arguments{
	\item{pheno}{
a numeric matrix of z-scores with rows and columns corresponding to samples and 
genes, respectively.
}
	\item{model}{
a character value specifying whether the original (if \code{global}) or extended 
(if \code{stratified}) model is used. 
}
	\item{metric}{
a character value specifying the metric to compute similarity scores. Currently, 
'cosine' and 'correlation' are supported (see \code{\link[PANR:assoScore]{assoScore}} for more details).
}
	\item{nPerm}{
a positive numeric or integer value specifying the number of permutations.
}
}
\details{
The function implements beta-mixture modelling in the application to Ewing's 
sarcoma. 
The function first initiates an object of class \code{\link[PANR:BetaMixture]{BetaMixture}} with
inputted \code{pheno}, \code{model} and \code{metric}.
The `x' component of the beta-mixture model is then fitted to permuted data for 
\code{nPerm} times, and the estimated shape parameters of beta distributions are 
summarized by taking the median values.
After fixing the shape parameters of the `x' component, the function invokes 
\code{\link[PANR:fitBM]{fitBM}} to fit the global or stratified beta-mixture model.
}
\value{
an object of class \code{\link[PANR:BetaMixture]{BetaMixture}} with beta-mixture modelling 
results stored at slot \code{result}
}
\references{
Arora S, Gonzales IM, Hagelstrom RT, et al. (2010). RNAi phenotype profiling of 
kinases identifies potential therapeutic targets in Ewing's sarcoma. Molecular 
Cancer, 9(1), 218.

Wang X, Castro MA, Mulder KW and Markowetz F (2012), Posterior association 
networks and enriched functional gene modules inferred from rich phenotypic 
perturbation screens, PLoS Computational Biology, doi:10.1371/journal.pcbi.1002566.
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}
\seealso{
\code{\link[PANR:BetaMixture]{BetaMixture}}, \code{\link[PANR:fitNULL]{fitNULL}}, \code{\link[PANR:fitBM]{fitBM}}
}
\examples{
\dontrun{
data(Arora2010, package="Mulder2012")
BMfitting(pheno=Arora2010, model="global", metric="cosine", nPerm=20)
}
}












