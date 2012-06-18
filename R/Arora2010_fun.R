###############################################################################
#1	Model fitting						      								 
###############################################################################
## Fit a beta-mixture model to functional associations (Global model)
Arora2010.BMfitting<-function(pheno, model="global", metric="cosine", 
	nPerm=20) {
	bm<-new("BetaMixture", pheno=pheno, metric=metric,model=model, order=1)
	bm<-fitNULL(bm, nPerm=nPerm, thetaNULL=c(alphaNULL=4, betaNULL=4),
		sumMethod="median", permMethod="keepRep", verbose=TRUE)
	bm<-fitBM(bm, para=list(zInit=NULL, thetaInit=c(alphaNeg=2,
		betaNeg=4, alphaNULL=bm@result$fitNULL$thetaNULL[["alphaNULL"]],
		betaNULL=bm@result$fitNULL$thetaNULL[["betaNULL"]],
		alphaPos=4, betaPos=2), gamma=NULL),
		ctrl=list(fitNULL=FALSE, tol=1e-1), verbose=TRUE, gradtol=1e-3)
	return(bm)
}
###############################################################################
#2	Inferring a PAN														   
###############################################################################
## Infer a functional association network given a fitted BetaMixture object
Arora2010.InferPAN<-function(bm, type="SNR", log=TRUE, sign=TRUE, 
	cutoff=log(10), filter=FALSE) {
	pan<-new("PAN", bm1=bm)
	pan<-infer(pan, para=list(type=type, log=log, sign=sign, cutoff=cutoff), 
		filter=filter, verbose=TRUE)
	pan<-buildPAN(pan, engine="RedeR", para=list(nodeSumCols=1:2, 
		nodeSumMethod="average", hideNeg=TRUE))
	return(pan)
}
###############################################################################
#3	Searching for functional modules 		
###############################################################################
## Search for functional modules using pvclust
Arora2010.ModuleSearchByPvclust<-function(pan, nboot=1000, metric="cosine2", 
	hclustMethod="average", filter=TRUE) {
#	##!patch pvclust to allow second-order clustering	
#	ns<-getNamespace("pvclust")
#	en<-as.environment("package:pvclust")
	
#	assignInNamespace("dist.pvclust",dist.pvclust4PAN,ns="pvclust", envir=ns)
#	dist.pvclust<-getFromNamespace("dist.pvclust", ns=getNamespace("pvclust"))
#	unlockBinding("parPvclust", ns)
#	assignInNamespace("parPvclust",parPvclust4PAN,ns="pvclust", envir=ns)
#	lockBinding("parPvclust", ns)
	
#	parPvclust<-getFromNamespace("parPvclust", ns)
#	if(is(getOption("cluster"), "cluster") && 
#			"package:snow" %in% search()) {		
#		clusterCall(getOption("cluster"), assignInNamespace, x="dist.pvclust", 
#			value=dist.pvclust4PAN, ns=ns)
#		clusterCall(getOption("cluster"), unlockBinding, sym="parPvclust", env=ns)
#		clusterCall(getOption("cluster"), assignInNamespace, x="parPvclust", 
#			value=parPvclust4PAN, ns=ns)
#		clusterCall(getOption("cluster"), lockBinding, sym="parPvclust", env=ns)
#	}
	pan<-pvclustModule(pan, nboot=nboot, metric=metric, hclustMethod=
					hclustMethod, filter=filter, verbose=TRUE, r=c(6:12/8))
	return(pan)
}

###############################################################################
#4	Pathway analysis	 				  	
###############################################################################
Arora2010.hypergeo<-function(pan, mod.pval.cutoff=0.05, mod.size.cutoff=4, 
	avg.degree.cutoff=0.5, filter.effects=TRUE) {
	univer<-rep(1, length(rownames(pan@bm1@pheno)))
	names(univer)<-toupper(rownames(pan@bm1@pheno))
	
	mod<-Arora2010.PANmoduleSelect(pan, mod.pval.cutoff, mod.size.cutoff, 
		avg.degree.cutoff, filter.effects)
	globalMods<-mod$globalMods
	mod<-mod$mod
	mods<-mod[globalMods]
	Arora2010.pathway<-function(mod) {
		gl<-rep(1, length(mod))
		names(gl)<-mod
		
		gl<-names(annotationConvertor(gl, species="Hs", initialIDs="Symbol", 
			finalIDs="Entrez.gene", keepMultipleMappings=TRUE, verbose=TRUE))
		univer<-names(annotationConvertor(univer, species="Hs", initialIDs="Symbol", 
			finalIDs="Entrez.gene", keepMultipleMappings=TRUE, verbose=TRUE))
		
		kegg<- KeggGeneSets(species="Hs")
		kegg1<-lapply(kegg, function(x) x[x%in%univer])
		names(kegg1)<-names(kegg)
		kegg<-kegg1
		kegg<-kegg[unlist(lapply(kegg, function(x) length(x)>=3))]

		rslt<-lapply(kegg, function(x) {hyperGeoTest(list(gl), univer, x)})
		expectedH<-unlist(lapply(rslt, function(x) x[4]))
		observedH<-unlist(lapply(rslt, function(x) x[5]))
		names(expectedH)<-names(kegg)
		names(observedH)<-names(kegg)
		pval<-unlist(lapply(rslt, function(x) x[6]))
		names(pval)<-names(kegg)
		
		pval<-p.adjust(pval, "BH")
		sel<-names(pval)[which(pval<0.05 & observedH>2)]
		
		pval.sel<-pval[sel]
		expH.sel<-expectedH[sel]
		obsH.sel<-observedH[sel]
		
		rslt<-data.frame(pathway=names(pval.sel), expected=expH.sel, observed=
			obsH.sel, "obs/exp"=obsH.sel/expH.sel,  "p-value"=pval.sel)
		xx <- as.list(KEGGPATHID2NAME)
		keggids<-substr(as.character(rslt[, 1]), 4, 8)
		keggterms<-unlist(xx[keggids])
		rslt<-data.frame(rslt, term=keggterms)
		rslt<-rslt[sort.list(rslt[, 4], decreasing=FALSE), ]
		rslt<-rslt[sort.list(rslt[, 5], decreasing=TRUE), ]	
	}
	rslt<-list()
	for(imod in 1:length(mods)) {
		mod<-mods[[imod]]

		rslt[[imod]]<-Arora2010.pathway(mod)
	}
	return(rslt)
}
###############################################################################
#5	Module visualization 
###############################################################################
##module visualization
Arora2010.module.visualize<-function(rdp, pan, mod.pval.cutoff=0.05, 
	mod.size.cutoff=4, avg.degree.cutoff=0.5, filter.effects=TRUE) {
	
	mod<-Arora2010.PANmoduleSelect(pan, mod.pval.cutoff, mod.size.cutoff, 
		avg.degree.cutoff, filter.effects)
		
	globalMods<-mod$globalMods
	moduleRef<-mod$moduleRef
	isParents<-mod$isParents
	mod<-mod$mod
#	rdp <- RedPort('MyPort')
#	calld(rdp)
	gigRef<-addGraph(rdp, subgraph(pan@iPAN, unique(unlist(mod))), zoom=60)
	iroot<-0
	for(i in 1:length(mod)) {
		ic<-(iroot)%%2+1
		ir<-floor(iroot/2)+1
		ml<-length(mod[[i]])
		if(i>1 && sum(isParents[[i]])>0) {							
			closestParent<-rev(which(isParents[[i]]))[1]
			moduleRef[[i]]<-nestNodes(rdp, nodes=mod[[i]], parent=
				moduleRef[[closestParent]], gcoord=c(50, 50), isAnchor=TRUE, 
				gscale=70, theme="tm1", gatt=list(nestFontSize=22, nestAlias=
				paste("p-value:",as.character(round(pan@modules$pval[
				as.numeric(names(mod)[i])], 3))),sep=""))
		} else {
			moduleRef[[i]]<-nestNodes(rdp, nodes=mod[[i]], parent=gigRef, 
				gcoord=c((ic-1)*40+20, (ir-1)*40+20), isAnchor=TRUE, gscale=40, 
				theme="tm1", gatt=list(nestFontSize=22, nestAlias=paste(
				"p-value:",as.character(round(pan@modules$pval[
				as.numeric(names(mod)[i])], 3))),sep=""))
			iroot<-iroot+1
		}
	}
	zscores<-format(round(pan@legend$nodeColor[, "zscore"], 3), nsmall=3)
	zscore.lab<-c(paste("(",paste(zscores[1:(length(zscores)-1)], 
		zscores[2:length(zscores)],sep=","), "]",sep=""),
		paste("> ", zscores[length(zscores)],sep="")
	  )
	addLegend.color(rdp, colvec=as.character(pan@legend$nodeColor[, "color"]), 
		position="topright",dxtitle=65, labvec=zscore.lab, type="node", 
		title="Loss of function")
	addLegend.size(rdp, sizevec=round(pan@legend$nodeSize[1:10*floor(nrow(
		pan@legend$nodeSize)/10), "diameter"], 2), position="bottomright",
		labvec=pan@legend$nodeSize[1:10*floor(nrow(pan@legend$nodeSize)/10), 
		"degree"], type="node", title="Node degree")
	addLegend.size(rdp, sizevec=round(pan@legend$edgeWidth[, "width"], 2), 
		position="bottomleft", ftsize=15,dxtitle=80, labvec=round(
		pan@legend$edgeWidth[, "SNR"], 2), title="Log(SNR)", type="edge")
	mergeOutEdges(rdp, lb=1, ub=15, nlev=3)
	relax(rdp)
}
Arora2010.PANmoduleSelect<-function(pan, mod.pval.cutoff=0.05, 
	mod.size.cutoff=4, avg.degree.cutoff=0.5, filter.effects=TRUE) {
	
	summ<-list()
	moduleRef<-NULL
	isParents<-NULL
	globalMods<-NULL
	avgdg<-NULL
	##pvclust p-values
	mod<-pan@modules$clusters
	names(mod)<-1:length(mod)
	summ[["All clusters"]]<-length(mod)
	
	mod<-mod[pan@modules$pval<mod.pval.cutoff]
	summ[[paste("p<", mod.pval.cutoff, sep="")]]<-length(mod)
	##module size
	if(length(mod)>0) {
		mod.size<-unlist(lapply(mod, length))
		inds<-which(mod.size>mod.size.cutoff & mod.size<nrow(pan@bm1@pheno)/2)
		if(length(inds)>0) {
			mod<-mod[inds]
			summ[["module size constraint"]]<-length(mod)
			##degree
			if(length(mod)>0) {
				mod.graph<-lapply(mod, function(x) {
					pan@graph[x, x]		
				})
				mod.graph<-lapply(mod.graph, function(x) sign(x))
				avgdg<-unlist(lapply(mod.graph, function(x) {
					mean(rowSums(x))/nrow(x)		
				}))
				inds<-which(avgdg>avg.degree.cutoff)
				if(length(inds)>0) {
					mod<-mod[inds]
					summ[["module density constraint"]]<-length(mod)
					avgdg<-avgdg[inds]
					mod.graph<-mod.graph[inds]
					##positive effects
					if(length(mod)>0) {
						if(filter.effects) {
							mod.g.effects<-lapply(mod.graph, function(x) {
								ids<-rownames(x)
								all(apply(pan@bm1@pheno[ids, ], 1, mean)<0)
								##all(pan@bm1@pheno[ids, ]>0)
							})
							mod.graph<-mod.graph[which(unlist(mod.g.effects))]

							inds<-which(unlist(mod.g.effects))
						} else {
							inds<-1:length(mod.graph)
						}
						
						if(length(inds)>0) {
							mod<-mod[inds]
							avgdg<-avgdg[inds]
							summ[["loss-of-function constraint"]]<-length(mod)
							##root modules
							if(length(mod)>0) {
								x<-NULL
								for(i in 1:length(mod)) {
									y<-sapply(setdiff(1:length(mod), i), 
										function(j) all(mod[[i]]%in%mod[[j]]))
									if(!(any(unlist(y)))) {
										x<-c(x, i)
									}
								}	
								mod<-mod[sort.list(unlist(lapply(mod, length)), 
									decreasing=TRUE)]
								moduleRef<-list()
								isParents<-list()
								globalMods<-c(1)		
								for(i in 1:length(mod)) {
									##find parents
									if(i>1) {
										isParents[[i]]<-sapply(1:(i-1), 
											function(pt) {
												if(all(mod[[i]]%in%mod[[pt]]))
													return(TRUE)
												else
													return(FALSE)
											})
										if(sum(isParents[[i]])==0)
											globalMods<-c(globalMods, i)
									}
								}
								summ[["root modules"]]<-length(globalMods)
							} else {
								summ[["root modules"]]<-0
								
							}
						}
					} else {						
						summ[["module density constraint"]]<-0
					}
				}
			} else {
				summ[["loss-of-function constraint"]]<-0
			}
		}	
	} else {
		summ[["module size constraint"]]<-0
	}
	return(list(mod=mod, moduleRef=moduleRef, isParents=isParents, 
		globalMods=globalMods, summary=summ, avgdg=avgdg))
}
