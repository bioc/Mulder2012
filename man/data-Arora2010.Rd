\name{data-Arora2010}
\alias{Arora2010}
\alias{Arora2010.pathway}
\alias{pw.rslt}
\alias{bm_Arora2010}
\alias{pan_Arora2010}
\title{
Raw and processed data for the application of PANs to Ewing's sarcoma
}                                                                        
\description{
The data used in this application is a matrix (572 x 8) of Z-scores from 
high-throughput RNAi screens run in duplicates targeting 572 human kinases in 
four Ewing's sarcoma cell lines: TC-32, TC-71, SK-ES-1 and RD-ES.
In these phenotyping screens, viability was assessed using a luminescence-based 
cell to quantify each gene's function in cancer cell growth and proliferation.
The screening data was corrected for plate row variations and normalized using 
Z-score method as described in Arora S et al., 2010. 

\code{Arora2010} is a matrix of z-scores summarized from the raw screening data 
in Arora2010.

\code{bm_Arora2010} is an object of \code{\link[PANR:BetaMixture]{BetaMixture}} 
including fitting results using the global BM model.

\code{pan_Arora2010} is an object of \code{\link[PANR:PAN]{PAN}} including 
the inferred functional network and modules identified by \code{\link[pvclust:pvclust]{pvclust}}.

\code{pw.rslt} is a list of hypeogeometric test results for 
overrepresentation of KEGG pathways in identified modules.
}
\usage{
##see example for details
}
\references{
Arora S, Gonzales IM, Hagelstrom RT, et al. (2010). RNAi phenotype profiling of 
kinases identifies potential therapeutic targets in Ewing's sarcoma. Molecular 
Cancer, 9(1), 218.
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}

\examples{
data(Arora2010, package="Mulder2012")
dim(Arora2010)
data(bm_Arora2010, package="Mulder2012")
bm_Arora2010
data(pan_Arora2010, package="Mulder2012")
pan_Arora2010
}










