\name{Mulder2012.InferPAN}
\alias{Mulder2012.InferPAN}
\title{
Infer a posterior association network given a fitted BetaMixture object
}                                                                        
\description{
The function infers a functional association network given an object of class 
\code{BetaMixture} that has been fitted. 
}
\usage{
Mulder2012.InferPAN(bm, type="SNR", log=TRUE, sign=TRUE, cutoff=log(10), 
filter=FALSE)
}
\arguments{
	\item{bm}{
an object of S4 class \code{BetaMixture}, which models the densities of
first-order similarities between genes (see \code{\link[PANR:BetaMixture]{BetaMixture}}).
}
	\item{type}{
a character value giving the type of edge weights: signal-to-noise ratio
('SNR'), posterior probability ratio ('PPR') or posterior probability ('PP')
}
	\item{log}{
a logical value specifying whether or not to compute logrithms for edge weights.
}
	\item{sign}{
a logical value specifying whether a signed graph should be inferred. It is only 
used when \code{type='SNR'}.
}
	\item{cutoff}{
a numeric value giving the threshold to tell the significance of an edge.
}
	\item{filter}{
a logical value specifying whether or not to filter out genes without any
significant association with all the other genes.
}
}
\details{
This function first initiates an object of class \code{\link[PANR:PAN]{PAN}} with
\code{bm}, which is an object of class \code{\link[PANR:BetaMixture]{BetaMixture}} that is supposed
to have been fitted.
The function \code{\link[PANR:infer]{infer}} is then invoked to infer a posterior
association network given the inputted arguments \code{type}, \code{log},
\code{sign}, \code{cutoff} and \code{filter}.
}
\value{
an object of class \code{\link[PANR:PAN]{PAN}} with the inferred posterior
association network stored at slot \code{graph}
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
\code{\link[PANR:PAN]{PAN}}, \code{\link[PANR:infer]{infer}}
}
\examples{
\dontrun{
data(bm_ext_Mulder2012, package="Mulder2012")
Mulder2012.InferPAN(bm=bm_ext_Mulder2012, type="SNR", log=TRUE, sign=TRUE, 
cutoff=log(10), filter=FALSE) 
}
}












