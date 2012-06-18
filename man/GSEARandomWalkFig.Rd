\name{GSEARandomWalkFig}
\alias{GSEARandomWalkFig}
\title{
Plot a figure illustrating enrichment analysis results
}                                                                        
\description{
The function invokes functions in R package \code{HTSanalyzeR} to plot results of 
enrichment analysis of functional interactions for protein-protein interactions.
}
\usage{
GSEARandomWalkFig(pheno, PPI, bm, what="pos")
}
\arguments{
	\item{pheno}{
a numeric matrix of z-scores with rows and columns corresponding to samples and 
genes, respectively.
}
	\item{PPI}{
an adjacency matrix of protein-protein interactions extracted from the PINdb 
(details in \code{\link[Mulder2012:Mulder2012.PPIPre]{Mulder2012.PPIPre}}).
}
	\item{bm}{
an object of S4 class \code{BetaMixture}, which models functional interactions 
between genes (see \code{\link[PANR:BetaMixture]{BetaMixture}}).
}
	\item{what}{
a character value specifying which component should be plotted: `pos' (positive 
association), `neg' (negative association) or `none' (lack of association).
}
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
\code{\link[HTSanalyzeR:gseaScores]{gseaScores}}, \code{\link[HTSanalyzeR:gseaPlots]{gseaPlots}}
}
\examples{
\dontrun{
data(Mulder2012, package="Mulder2012")
data(Mulder2012.PPI, package="Mulder2012")
data(bm_Mulder2012, package="Mulder2012")
GSEARandomWalkFig(Mulder2012, PPI, bm_Mulder2012, what="pos")
}
}












