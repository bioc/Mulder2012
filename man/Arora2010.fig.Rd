\name{Arora2010.fig}
\alias{Arora2010.fig}
\title{
A pipeline function to reproduce the main figures of network analyses in Wang 
et al 2012 for Arora S et al 2010
}                                                                        
\description{
This is a pipeline function for the user to reproduce the main figures regarding 
network analyses in Wang X et al., 2012 for RNAi screening data of Ewing's 
sarcoma in Arora S et al., 2010.
}
\usage{
Arora2010.fig(what="ALL")
}
\arguments{
	\item{what}{
a character value specifying what to plot: `NULLfitting', `BMfitting', `sigMod' 
or `ALL' (see details).
}
}
\details{
The user can run this function with a proper argument to reproduce the main 
figures regarding network analyses in Wang X et al., 2012 for RNAi screening 
data of Ewing's sarcoma in Arora S et al., 2010.

These are possible options to plot:

`NULLfitting': a figure illustrating the fitting results of a beta distribution 
to permuted RNAi screening data. 

`BMfitting': a figure illustrating the fitting results of a beta-mixture 
model to the real RNAi screening data. 

`sigMod': significant functional modules identified by in \code{\link[pvclust:pvclust]{pvclust}}.

`pathway': a figure illustrating significantly overrepresented KEGG pathways.
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
\code{\link[pvclust:pvclust]{pvclust}}, \code{\link[RedeR:RedeR]{RedeR}}
}
\examples{
\dontrun{
Arora2010.fig("BMfitting")
}
}
