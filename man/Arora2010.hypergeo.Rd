\name{Arora2010.hypergeo}
\alias{Arora2010.hypergeo}
\title{
Hypergeometric tests for pathway analysis
}                                                                        
\description{
The function is used to perform hypergeometric tests of overrepresentation of 
KEGG pathways in inferred functional modules in the application to Ewing's 
sarcoma. 
}
\usage{
Arora2010.hypergeo(pan, mod.pval.cutoff=0.05, mod.size.cutoff=4, 
avg.degree.cutoff=0.5, filter.effects=TRUE)
}
\arguments{
	\item{pan}{
an object of class \code{PAN}, which includes the inferred posterior association 
network and functional modules.
}
	\item{mod.pval.cutoff}{
a numeric value specifying the p-value cutoff for selecting significant functional 
modules.
}
	\item{mod.size.cutoff}{
an integer value specifying the minimal size of functional modules.
}
	\item{avg.degree.cutoff}{
a numeric value specifying the cutoff of module density, which is the ratio of
predicted significant associations to all possible associations.
}
	\item{filter.effects}{
a logical value specifying whether or not to filter out modules that are not of 
interest. In this application to Ewing's sarcoma, modules associated with 
positive loss-of-function can be filtered out because we are only interested in 
modules related to inhibiting growth and proliferation of cancer cells.
}
}
\details{
The function is used to do pathway analyses in the application of PANs to Ewing's 
sarcoma. Hypergeometric tests are performed to test the overrepresentation of 
KEGG pathways in inferred functional modules. Four arguments can be specified to 
filter modules. More details about module filtering procedures can be found in 
the Methods section of Wang X et al., 2012.
}
\value{
a list including hypergeometric test result for each module. 
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
\code{\link[PANR:PAN]{PAN}}, \code{\link[PANR:infer]{infer}}
}
\examples{
\dontrun{
data(pan_Arora2010, package="Mulder2012")
Arora2010.hypergeo(pan_Arora2010, mod.pval.cutoff=0.05, mod.size.cutoff=4, 
avg.degree.cutoff=0.5, filter.effects=TRUE)
}
}












