\name{Arora2010.ModuleSearchByPvclust}
\alias{Arora2010.ModuleSearchByPvclust}
\title{
Searching for enriched functional modules by pvclust
}                                                                        
\description{
This function invokes \code{\link[PANR:pvclustModule]{pvclustModule}} to search for enriched 
functional gene modules by hierarchical clustering with multiscale bootstrap 
resampling (using R package \code{pvclust}).
}
\usage{
Arora2010.ModuleSearchByPvclust(pan, nboot=1000, metric="cosine2", hclustMethod=
"average", filter=TRUE)
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
\code{\link[PANR:PAN]{PAN}}, \code{\link[PANR:pvclustModule]{pvclustModule}}
}
\examples{
\dontrun{
data(pan_Arora2010, package="Mulder2012")
pan_Arora2010<-ModuleSearchByPvclust(pan=pan_Arora2010, nboot=10000, metric=
"cosine2", hclustMethod="average", filter=TRUE)
summarize(pan_Arora2010, what="module")
}
}












