\name{Mulder2012.RNAiPre}
\alias{Mulder2012.RNAiPre}
\title{
Preprocessing of the raw RNAi screens
}                                                                        
\description{
The function performs normalization and summarization of RNAi screening data 
for Mulder K et al. 2012 and Wang X et al. 2012.
}
\usage{
Mulder2012.RNAiPre(rawScreenData, rawScreenAnnotation)
}
\arguments{
	\item{rawScreenData}{
a data frame of the raw screening data. 
}
	\item{rawScreenAnnotation}{
a data frame of annotations. 
}
}
\value{
a numeric matrix of z-scores with rows and columns specifying samples and genes, 
respectively.
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
\examples{
\dontrun{
Mulder2012<-Mulder2012.RNAiPre()
}
}













