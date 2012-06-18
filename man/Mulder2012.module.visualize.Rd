\name{Mulder2012.module.visualize}
\alias{Mulder2012.module.visualize}
\title{
Infer a posterior association network given a fitted BetaMixture object
}                                                                        
\description{
The function represent inferred significant functional modules in \code{RedeR}, which 
is a bioconductor package for network visualization. 
}
\usage{
Mulder2012.module.visualize(rdp, pan, mod.pval.cutoff=0.05, mod.size.cutoff=4, 
avg.degree.cutoff=0.5, edgeWidthLeg=TRUE, filter.effects=TRUE)
}
\arguments{
	\item{rdp}{
an object of class \code{RedPort}, which builds the interface between R and RedeR 
Java server.
}
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
	\item{edgeWidthLeg}{
a logical value specifying whether or not the legend of edge width should be 
created.
}
	\item{filter.effects}{
a logical value specifying whether or not to filter out modules that are not of 
interest. In this application to epidermal stem cells, modules associated with 
negative loss-of-function can be filtered out because we are only interested in 
modules related to stem cell self-renewal.
}
}
\details{
The inferred functional modules are represented in a very compact layout in 
\code{RedeR}, as many of them are nested to each other due to the nature of 
hierarchical clustering. Four arguments can be specified to filter modules. 
More details about module filtering procedures can be found in the Methods 
section of Wang X et al., 2012.
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
\code{\link[RedeR:RedeR]{RedeR}}
}
\examples{
\dontrun{
library(RedeR)
rdp <- RedPort('MyPort')
calld(rdp)
data(pan_ext_Mulder2012, package="Mulder2012")
Mulder2012.module.visualize(rdp, pan_ext_Mulder2012, mod.pval.cutoff=0.05, 
mod.size.cutoff=4, avg.degree.cutoff=0.5)
}
}












