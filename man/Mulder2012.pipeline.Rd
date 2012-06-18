\name{Mulder2012.pipeline}
\alias{Mulder2012.pipeline}
\title{
A pipeline to do network analyses for Mulder 2012
}                                                                        
\description{
This pipeline function reproducing results in Mulder 2012.
}
\usage{
Mulder2012.pipeline(par4BM=list(model="global", metric="cosine", nPerm=20),
par4PAN=list(type="SNR", log=TRUE, sign=TRUE, cutoff=log(10), filter=FALSE),
par4ModuleSearch=list(nboot=10000, metric="cosine2", hclustMethod="average", 
filter=FALSE))
}
\arguments{
	\item{par4BM}{
a list of parameters for beta-mixture modelling (see \code{\link[Mulder2012:Mulder2012.BMfitting]{Mulder2012.BMfitting}}) for 
details. 
}
	\item{par4PAN}{
a list of parameters for inferring a posterior association network (see
\code{\link[Mulder2012:Mulder2012.InferPAN]{Mulder2012.InferPAN}}) for details.
}
	\item{par4ModuleSearch}{
a list of parameters for seaching enriched functional gene modules (see
\code{\link[Mulder2012:Mulder2012.ModuleSearchByPvclust]{Mulder2012.ModuleSearchByPvclust}}) for details.
}
}
\details{
This function reproduces results of network inference in Mulder K et al., 2012:

(1) Beta-mixture modelling of functional gene associations (see \code{\link[Mulder2012:Mulder2012.BMfitting]{Mulder2012.BMfitting}}))
for details.

(2) Inferring a functional association network (see \code{\link[Mulder2012:Mulder2012.InferPAN]{Mulder2012.InferPAN}}))
for details.

(3) Searching enriched gene modules by pvclust (see \code{\link[Mulder2012:Mulder2012.ModuleSearchByPvclust]{Mulder2012.ModuleSearchByPvclust}}))
for details.

(4) enrichment analyses of posterior probabilities in protein-protein
interactions (see \code{\link[Mulder2012:Mulder2012.PPIenrich]{Mulder2012.PPIenrich}})) for details. 
}
\value{
All results will be saved in the current work directory.
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
\code{\link[Mulder2012:Mulder2012.BMfitting]{Mulder2012.BMfitting}}, \code{\link[Mulder2012:Mulder2012.InferPAN]{Mulder2012.InferPAN}}, \code{\link[Mulder2012:Mulder2012.ModuleSearchByPvclust]{Mulder2012.ModuleSearchByPvclust}}, 
\code{\link[Mulder2012:Mulder2012.PPIenrich]{Mulder2012.PPIenrich}}
}
\examples{
\dontrun{
Mulder2012.pipeline()
}
}












