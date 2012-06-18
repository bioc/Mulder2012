\name{Mulder2012.fig}
\alias{Mulder2012.fig}
\title{
A pipeline function to reproduce the main figures of network analyses in Mulder 
K et al 2012 and Wang et al 2012
}                                                                        
\description{
This is a pipeline function for the user to reproduce the main figures regarding 
network analyses in Mulder K et al., 2012 and Wang X et al., 2012.
}
\usage{
Mulder2012.fig(what="ALL")
}
\arguments{
	\item{what}{
a character value specifying what to plot: `NULLfitting', `BMfitting', 
`PPIenrich', `sigMod', `selMod' or `ALL' (see details).
}
}
\details{
The user can run this function with a proper argument to reproduce the main 
figures regarding network analyses in Mulder K et al. 2012 and Wang X et al. 
2012. 

These are the possible options to plot:

`NULLfitting': a figure illustrating the fitting results of a beta distribution 
to permuted RNAi screening data. 

`BMfitting': a figure illustrating the fitting results of a beta-mixture 
model to the real RNAi screening data. 

`PPIenrich': figures describing enrichment results of posterior probabilities 
in protein-protein interactions. 

`sigMod': significant functional modules identified by in \code{\link[pvclust:pvclust]{pvclust}}.

`selMod': a figure of the module selected for further validation using 
combinatorial knock-down experiments.
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
\code{\link[pvclust:pvclust]{pvclust}}, \code{\link[RedeR:RedeR]{RedeR}}
}
\examples{
\dontrun{
Mulder2012.fig("BMfitting")
}
}
