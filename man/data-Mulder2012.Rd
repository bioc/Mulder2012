\name{data-Mulder2012}
\alias{Mulder2012}
\alias{Mulder2012.rawScreenData}
\alias{rawScreenData}
\alias{Mulder2012.rawScreenAnnotation}
\alias{rawScreenAnnotation}
\alias{Mulder2012.PINdbProteins}
\alias{PINdbProteins}
\alias{Mulder2012.PPI}
\alias{PPI}
\alias{proteinModules}
\alias{Mulder2012.PPIEnrich}
\alias{PPIEnrich_Mulder2012}
\alias{bm_Mulder2012}
\alias{pan_Mulder2012}
\alias{bm_ext_Mulder2012}
\alias{pan_ext_Mulder2012}
\title{
Raw and processed data for the application to epidermal stem cells
}                                                                        
\description{
The data set we use here comes from RNA interference screening on 332 
chromatin modifiers under five different conditions: vehicle, AG1478, 
BMP2/7, AG1478 BMP2/7 and serum stimulation (details in Mulder K et al., 2012).

\code{rawScreenData} and \code{rawScreenAnnotation} are the raw RNAi screening data 
and annotations, respectively.

\code{Mulder2012} is a matrix of z-scores summarized from the raw data and 
annotation.

\code{PINdbProteins} and \code{PPI} are the original PINdb database and extracted 
protein-protein interactions matrix, respectively.
\code{proteinModules} is a list of protein complexes in the PINdb database.

\code{PPIEnrich} is a list of enrichment analyses results including: 
\code{gs}--a matrix of running scores of the real and permuted data in the 
enrichment analyses results. 
\code{pval}--p-values derived from enrichment analyses.

\code{bm_Mulder2012} and \code{bm_ext_Mulder2012} is an object of \code{\link[PANR:BetaMixture]{BetaMixture}} 
including fitting results using the global and stratified BM model, respectively.

\code{pan_Mulder2012} and \code{pan_ext_Mulder2012} is an object of \code{\link[PANR:PAN]{PAN}} including 
the inferred functional network and modules identified by \code{\link[pvclust:pvclust]{pvclust}} 
with and without incorporating prior information, respectively.
}
\usage{
##see example for details
}
\references{
Mulder KW, Wang X, Escriu C, et al. (2012). Diverse epigenetic strategies 
interact to control epidermal differentiation. Nature Cell Biology. 
doi:10.1038/ncb2520.
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}

\examples{
data(Mulder2012, package="Mulder2012")
dim(Mulder2012)
data(bm_Mulder2012, package="Mulder2012")
bm_Mulder2012
data(pan_Mulder2012, package="Mulder2012")
pan_Mulder2012
}










