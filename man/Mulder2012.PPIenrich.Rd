\name{Mulder2012.PPIenrich}
\alias{Mulder2012.PPIenrich}
\title{
Enrichment analyses of posterior probabilities in protein-protein interactions
}                                                                        
\description{
This function performs enrichment analyses of protein-protein interactions in 
posterior probabilities of functional associations belonging to `+', `x', `-' 
mixture component.
}
\usage{
Mulder2012.PPIenrich(pheno, PPI, bm)
}
\arguments{
	\item{pheno}{
a numeric matrix of z-scores with rows and columns corresponding to samples and 
genes, respectively.
}
	\item{PPI}{
a binary adjacency matrix of protein-protein interactions.
}
	\item{bm}{
an object of \code{\link[PANR:BetaMixture]{BetaMixture}}.
}
}
\details{
After fitting a beta-mixture model to the densities of functional gene 
associations, a posterior probability for each pair of genes belonging to 
positive, negative or lack of association component is computed. Using the 
vector of posterior probabilities belonging to each one of the three components, 
an enrichment analysis is performed based on R package \code{\link[HTSanalyzeR:HTSanalyzeR]{HTSanalyzeR}}, 
which does gene set enrichment analysis. The `gene set' we use here is a vector 
of protein-protein interactions, which is extracted from PINdb (Proteins 
Interacting in the Nucleus).
}
\value{
an object of class \code{\link[PANR:BetaMixture]{BetaMixture}} with beta-
mixture modelling results stored at slot \code{result}
}
\references{
Luc PV & Tempst P, PINdb (2004): a database of nuclear protein complexes from 
human and yeast. Bioinformatics 20(9):1413-5.

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
\code{\link[PANR:BetaMixture]{BetaMixture}}, \code{\link[HTSanalyzeR:HTSanalyzeR]{HTSanalyzeR}}
}
\examples{
\dontrun{
data(Mulder2012, package="Mulder2012")
data(bm_Mulder2012, package="Mulder2012")
data(Mulder2012_PPI, package="Mulder2012")
PPIenrich(pheno=Mulder2012, PPI=PPI, bm=bm_Mulder2012)
}
}
