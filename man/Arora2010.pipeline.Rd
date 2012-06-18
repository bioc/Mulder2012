\name{Arora2010.pipeline}
\alias{Arora2010.pipeline}
\title{
A pipeline to do network analyses in Wang et al 2012 for Arora S et al 2010
}                                                                        
\description{
This pipeline function reproducing results in Wang X et al., 2012 for RNAi 
screening data of Ewing's sarcoma in Arora S et al., 2010.
}
\usage{
Arora2010.pipeline(par4BM=list(model="global", metric="cosine", nPerm=20),
par4PAN=list(type="SNR", log=TRUE, sign=TRUE, cutoff=log(10), filter=FALSE),
par4ModuleSearch=list(nboot=10000, metric="cosine2", hclustMethod="average", 
filter=FALSE))
}
\arguments{
	\item{par4BM}{
a list of parameters for beta-mixture modelling (see \code{\link[Mulder2012:Arora2010.BMfitting]{Arora2010.BMfitting}}) for 
details. 
}
	\item{par4PAN}{
a list of parameters for inferring a posterior association network (see
\code{\link[Mulder2012:Arora2010.InferPAN]{Arora2010.InferPAN}}) for details.
}
	\item{par4ModuleSearch}{
a list of parameters for seaching enriched functional gene modules (see
\code{\link[Mulder2012:Arora2010.ModuleSearchByPvclust]{Arora2010.ModuleSearchByPvclust}}) for details.
}
}
\details{
This function reproduces results of network inference in Wang X et al., 2012 for 
RNAi screening data of Ewing's sarcoma in Arora S et al., 2010.:

(1) Beta-mixture modelling of functional gene associations (see \code{\link[Mulder2012:Arora2010.BMfitting]{Arora2010.BMfitting}}))
for details.

(2) Inferring a posterior association network (see \code{\link[Mulder2012:Arora2010.InferPAN]{Arora2010.InferPAN}}))
for details.

(3) Searching enriched gene modules by pvclust (see \code{\link[Mulder2012:Arora2010.ModuleSearchByPvclust]{Arora2010.ModuleSearchByPvclust}}))
for details.

(4) Pathway analyses (see \code{\link[Mulder2012:Arora2010.hypergeo]{Arora2010.hypergeo}})) for details. 
}
\value{
All results will be saved in the current work directory.
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
\code{\link[Mulder2012:Arora2010.BMfitting]{Arora2010.BMfitting}}, \code{\link[Mulder2012:Arora2010.InferPAN]{Arora2010.InferPAN}}, \code{\link[Mulder2012:Arora2010.ModuleSearchByPvclust]{Arora2010.ModuleSearchByPvclust}}, 
\code{\link[Mulder2012:Arora2010.hypergeo]{Arora2010.hypergeo}}
}
\examples{
\dontrun{
Arora2010.pipeline()
}
}












