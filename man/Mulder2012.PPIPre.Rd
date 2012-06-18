\name{Mulder2012.PPIPre}
\alias{Mulder2012.PPIPre}
\title{
Extracting protein-protein interactions from PINdb
}                                                                        
\description{
The function extracts protein-protein interaction information to a binary
adjacency matrix from PINdb (Protein Interacting in the Nucleus). 
}
\usage{
Mulder2012.PPIPre()
}
\value{
a list consisting of a binary adjacency matrix of protein-protein interactions
and a list of protein modules (complexes).
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
\examples{
\dontrun{
ppi<-Mulder2012.PPIPre()
}
}













