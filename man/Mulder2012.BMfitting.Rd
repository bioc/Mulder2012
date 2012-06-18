\name{Mulder2012.BMfitting}
\alias{Mulder2012.BMfitting}
\alias{Mulder2012.BMfitting.extended}
\title{
Fitting a beta-mixture model to association scores
}                                                                        
\description{
The function fits a three-beta mixture model to functional association scores 
computed from rich phenotyping screens in Mulder K et al., 2012 and Wang X et 
al., 2012.
}
\usage{
Mulder2012.BMfitting(pheno, model="global", metric="cosine", nPerm=100)
Mulder2012.BMfitting.extended(pheno, model="stratified", metric="cosine", 
nPerm=20, partition) 
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
	\item{partition}{
a numeric of gene partition labels (e.g. c(rep(1, 100), rep(2, 20)) is a valid 
vector of partition labels for a vector of associations falling into two 
categories of interaction types 1 and 2).
}
}
\details{
The two functions implement beta-mixture modelling in the application to epidermal 
stem cells. 
Both functions first initiate an object of class \code{\link[PANR:BetaMixture]{BetaMixture}} with
inputted \code{pheno}, \code{model} and \code{metric} (and \code{partition} for the 
stratified model).
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
Mulder KW, Wang X, Escriu C, et al. (2012). Diverse epigenetic strategies 
interact to control epidermal differentiation. Nature Cell Biology. 
doi:10.1038/ncb2520.

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
data(Mulder2012, package="Mulder2012")
BMfitting(pheno=Mulder2012, model="global", metric="cosine", nPerm=100)
}
}












