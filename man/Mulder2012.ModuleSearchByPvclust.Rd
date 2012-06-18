\name{Mulder2012.ModuleSearchByPvclust}
\alias{Mulder2012.ModuleSearchByPvclust}
\title{
Searching for enriched functional modules by pvclust
}                                                                        
\description{
This function invokes \code{\link[PANR:pvclustModule]{pvclustModule}} to search for enriched 
functional gene modules by hierarchical clustering with multiscale bootstrap 
resampling (using R package \code{pvclust}).
}
\usage{
Mulder2012.ModuleSearchByPvclust(pan, nboot=10000, metric="cosine2", 
hclustMethod="average", filter=TRUE)
}
\arguments{
	\item{pan}{
an object of class \code{\link[PANR:PAN]{PAN}}, in which a posterior association
network has already been inferred.
}
	\item{nboot}{
a numeric value giving the number of bootstraps for \code{pvclust}.
}
	\item{metric}{
a character value specifying which distance metric to use for the hierarchical
clustering: 'cosine' (first-order), 'cosine2' (second-order), 'correlation',
'uncentered', 'abscor' or those allowed by the argument 'method' in \code{\link[stats:dist]{dist}}.
}
	\item{hclustMethod}{
the agglomerative method used in hierarchical clustering: 'average', 'ward',
'single', 'complete', 'mcquitty', 'median' or 'centroid' (see the argument
\code{method} in \code{\link[stats:hclust]{hclust}} for more details).
}
	\item{filter}{
a logical value specifying whether or not to filter out screening data of
genes without significant associations with all the other genes.
}
}
\value{
This function will return an object of class \code{PAN} with inferred
gene modules (\code{modules$clusters}) and corresponding p-values (\code{modules$pval})
updated in slot 'modules'.
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
\code{\link[PANR:PAN]{PAN}}, \code{\link[PANR:pvclustModule]{pvclustModule}}
}
\examples{
\dontrun{
data(pan_ext_Mulder2012, package="Mulder2012")
pan_ext_Mulder2012<-ModuleSearchByPvclust(pan=pan_ext_Mulder2012, nboot=10000, 
metric="cosine2", hclustMethod="average", filter=TRUE)
summarize(pan_ext_Mulder2012, what="module")
}
}












