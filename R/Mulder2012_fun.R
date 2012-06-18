###############################################################################
#1	Data preprocessing						    
###############################################################################
## Extracting protein-protein interactions from PINdb
Mulder2012.PPIPre<-function() {
	data("Mulder2012.PINdbProteins", package="Mulder2012")
	##proteins
	proteins<-unique(as.character(PINdbProteins[, "officialName"]))
	nrP<-length(proteins)
	PPI<-matrix(0, nrP, nrP)
	dimnames(PPI)<-list(sort(proteins), sort(proteins))
	##protein complex
	proteinComps<-unique(as.integer(as.character(PINdbProteins[, "cid"])))
	proteinModules<-list()
	sapply(1:length(proteinComps), function(i) {
		proteinModules[[as.character(i)]]<<-sort(as.character(PINdbProteins[
			PINdbProteins[,"cid"]==proteinComps[i], "officialName"]))
		PPI[proteinModules[[as.character(i)]], proteinModules[[
			as.character(i)]]]<<-1
	})
	return(list(PPI=PPI, proteinModules=proteinModules))
}
## Preprocessing RNAi screens
Mulder2012.RNAiPre<-function(rawScreenData, rawScreenAnnotation) {

	plts<-unique(rawScreenData[, "PLATE"])
	headers<-colnames(rawScreenData)[4:ncol(rawScreenData)]
	ctrWELLs<-c("H02", "H05")
	gsymbols<-as.character(rawScreenAnnotation[, 4])
	gsymbols<-gsymbols[!gsymbols%in%c("siTG1 A","siTG1 B")]
	##
	gsymbols[5]<-"PBRM1_r1"
	gsymbols[27]<-"PBRM1_r2"
	Mulder2012<-matrix(0, length(gsymbols), length(headers))
	dimnames(Mulder2012)<-list(gsymbols, headers)
	##plate-wise analysis
	for(p in 1:length(plts)) {
		rawscn<-rawScreenData[rawScreenData[, "PLATE"]==plts[p], ]
		anno<-rawScreenAnnotation[rawScreenAnnotation[, "PLATE"]==plts[p], ]
		##ctrol signals
		rawctr<-rawscn[rawscn[, "CHANNEL"]==700, ]
		rawctr.mat<-as.matrix(rawctr[, headers])
		rownames(rawctr.mat)<-as.character(rawctr[, "WELL"])
		rawctr.mat<-rawctr.mat[!rownames(rawctr.mat)%in%ctrWELLs, ]
		##screens
		scn<-rawscn[rawscn[, "CHANNEL"]==800, ]
		scn.mat<-as.matrix(scn[, headers])
		rownames(scn.mat)<-as.character(scn[, "WELL"])	
		scn.mat<-scn.mat[!rownames(scn.mat)%in%ctrWELLs, ]
		##
		ctrwells<-rownames(rawctr.mat)
		scnwells<-rownames(scn.mat)
		if(!all(union(ctrwells, scnwells)%in%intersect(ctrwells, scnwells)))
			stop("-Mulder2011 error: screening data is not complete!")
		else {
			siTG1<-as.matrix(rawscn[rawscn[, "WELL"]%in%ctrWELLs & rawscn[, 
				"CHANNEL"]==800, headers])
			avg.siTG1<-as.numeric(sapply(1:(length(headers)/3), function(x) {
				rep(mean(siTG1[, ((x-1)*3+1):(x*3)]), 3)			
			}))
			scn.norm<-t(t(scn.mat)-avg.siTG1)/rawctr.mat
			scn.zscores<-matrix(0, nrow(scn.norm), ncol(scn.norm))
			dimnames(scn.zscores)<-dimnames(scn.norm)
			for(x in 1:length(headers)) {
				scn.zscores[, x]<-(scn.norm[, x]-mean(scn.norm[, x]))/sd(scn.norm[, x])
			}
			thisgsymbols<-as.character(anno[match(scnwells, as.character(
				anno[, "WELL"])), "SYMBOL"])
			if("PBRM1" %in% thisgsymbols) {
				thisgsymbols[thisgsymbols=="PBRM1"]<-c("PBRM1_r1", "PBRM1_r2")
			}
			Mulder2012[thisgsymbols, colnames(scn.zscores)]<-scn.zscores
		}
	}
	Mulder2012<-Mulder2012[sort(rownames(Mulder2012)),]
	colnames(Mulder2012)<-as.character(sapply(colnames(Mulder2012), function(x) 
		strsplit(x,split='_')[[1]][1]))
	Mulder2012
}
###############################################################################
#	2	Model fitting						  
###############################################################################
## Fit a beta-mixture model to functional associations (Global model)
Mulder2012.BMfitting<-function(pheno, model="global", metric="cosine", nPerm=100) {
	bm<-new("BetaMixture", pheno=pheno, metric=metric, order=1, model=model)
	bm<-fitNULL(bm, nPerm=nPerm, thetaNULL=c(alphaNULL=4, betaNULL=4), 
			sumMethod="median", permMethod="keepRep", verbose=TRUE)
	bm<-fitBM(bm, para=list(zInit=NULL, 
					thetaInit=c(alphaNeg=2, betaNeg=4, 
							alphaNULL=bm@result$fitNULL$thetaNULL[["alphaNULL"]], 
							betaNULL=bm@result$fitNULL$thetaNULL[["betaNULL"]], 
							alphaPos=4, betaPos=2), gamma=NULL), 
			ctrl=list(fitNULL=FALSE, tol=1e-3), verbose=TRUE)
	return(bm)
}
## Fit a beta-mixture model to functional associations (Extended model)
Mulder2012.BMfitting.extended<-function(pheno, model="stratified", metric="cosine", 
	nPerm=20, partition) {
	bm<-new("BetaMixture", pheno=pheno, metric=metric, order=1, model=model, 
		partition=partition)
	bm<-fitNULL(bm, nPerm=nPerm, thetaNULL=c(alphaNULL=4, betaNULL=4), 
			sumMethod="median", permMethod="keepRep", verbose=TRUE)
	bm<-fitBM(bm, para=list(zInit=NULL, 
					thetaInit=c(alphaNeg=2, betaNeg=4, 
							alphaNULL=bm@result$fitNULL$thetaNULL[["alphaNULL"]], 
							betaNULL=bm@result$fitNULL$thetaNULL[["betaNULL"]], 
							alphaPos=4, betaPos=2), gamma=NULL), 
			ctrl=list(fitNULL=FALSE, tol=1e-3), verbose=TRUE)
	return(bm)
}
###############################################################################
#	3	Inferring a PAN					
###############################################################################
## Infer a functional association network given a fitted BetaMixture object
Mulder2012.InferPAN<-function(bm, type="SNR", log=TRUE, sign=TRUE, cutoff=log(10), 
	filter=FALSE) {
	pan<-new("PAN", bm1=bm)
	pan<-infer(pan, para=list(type=type, log=log, sign=sign, cutoff=cutoff), 
		filter=filter, verbose=TRUE)
	return(pan)
}
###############################################################################
#4	Enrichment analysis					
###############################################################################
##Enrichment analyses of posterior probabilities in protein-protein interactions
Mulder2012.PPIenrich<-function(pheno, PPI, bm) {
	genes<-rownames(pheno)
	nrg<-nrow(pheno)
	ppi.genes<-rownames(PPI)
	
	PPI.f<-matrix(0, nrg, nrg)
	dimnames(PPI.f)<-list(genes, genes)
	inds<-match(toupper(ppi.genes), toupper(genes))
	inds.ppi<-which(!is.na(inds))
	inds.rnai<-inds[!is.na(inds)]
	PPI.f[inds.rnai, inds.rnai]<-PPI[inds.ppi, inds.ppi]
	
	upper.tri.inds<-which(upper.tri(PPI.f), arr.ind=TRUE)
	PPI.df<-data.frame(gene1=genes[upper.tri.inds[,"row"]], gene2=genes[
		upper.tri.inds[,"col"]], PPI=PPI.f[upper.tri(PPI.f)], 
		neg=bm@result$fitBM$z[, "-"], none=bm@result$fitBM$z[, "x"], 
		pos=bm@result$fitBM$z[, "+"])
	gs<-list()
	pval<-list()
	interList.names<-paste(PPI.df[,"gene1"], PPI.df[,"gene2"],sep="~")
	for(i in c("neg","none","pos")) {
		cat("Running enrichment analysis of PPI in '", i, 
			"' posterior probabilities ...","\n")
		interList<-PPI.df[, i]
		names(interList)<-interList.names
		interSet<-interList.names[which(PPI.df[,"PPI"]==1)]
		gsc<-list(interSet)
		names(gsc)[1]<-"PPI"
		interList<-interList[sort.list(interList, decreasing=TRUE)]
		
		gs[[i]]<-list(Permutation.scores=NULL, Observed.scores=NULL)
		for(j in 1:10) {
			gs.temp<-collectionGsea(gsc, interList, exponent=1, nPermutations=
				1000, minGeneSetSize=1, verbose=FALSE)
			if(j==1) {
				gs[[i]]$Permutation.scores<-gs.temp$Permutation.scores
				gs[[i]]$Observed.scores<-gs.temp$Observed.scores
			} else 
				gs[[i]]$Permutation.scores<-cbind(gs[[i]]$Permutation.scores, 
					gs.temp$Permutation.scores)
		}
		pval[[i]] <- permutationPvalueCollectionGsea(permScores = 
			gs[[i]]$Permutation.scores,dataScores = gs[[i]]$Observed.scores)	
	}
	return(list(gs=gs, pval=pval))
}

###############################################################################
#5	Searching for functional modules 			
###############################################################################
## Search for functional modules using pvclust
Mulder2012.ModuleSearchByPvclust<-function(pan, nboot=10000, 
	metric="cosine2", hclustMethod="average", filter=TRUE) {
	##!patch pvclust to allow second-order clustering	
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
#		clusterCall(getOption("cluster"), unlockBinding, sym="parPvclust", 
#			env=ns)
#		clusterCall(getOption("cluster"), assignInNamespace, x="parPvclust", 
#			value=parPvclust4PAN, ns=ns)
#		clusterCall(getOption("cluster"), lockBinding, sym="parPvclust", 
#			env=ns)
#	}
	pan<-pvclustModule(pan, nboot=nboot, metric=metric, hclustMethod=
		hclustMethod, filter=filter, verbose=TRUE)
	return(pan)
}

###############################################################################
#6	Figures					 		
###############################################################################
GSEARandomWalkFig<-function(pheno, PPI, bm, what="pos") {
	##require(HTSanalyzeR)
	genes<-rownames(pheno)
	nrg<-nrow(pheno)
	ppi.genes<-rownames(PPI)
	
	PPI.f<-matrix(0, nrg, nrg)
	dimnames(PPI.f)<-list(genes, genes)
	inds<-match(toupper(ppi.genes), toupper(genes))
	inds.ppi<-which(!is.na(inds))
	inds.rnai<-inds[!is.na(inds)]
	PPI.f[inds.rnai, inds.rnai]<-PPI[inds.ppi, inds.ppi]
	
	upper.tri.inds<-which(upper.tri(PPI.f), arr.ind=TRUE)
	PPI.df<-data.frame(gene1=genes[upper.tri.inds[,"row"]], gene2=
		genes[upper.tri.inds[,"col"]], PPI=PPI.f[upper.tri(PPI.f)], neg=
		bm@result$fitBM$z[, "-"], none=bm@result$fitBM$z[, "x"], pos=
		bm@result$fitBM$z[, "+"])
	gs<-list()
	pval<-list()
	interList.names<-paste(PPI.df[,"gene1"], PPI.df[,"gene2"],sep="~")
	
	interList<-PPI.df[, what]
	names(interList)<-interList.names
	interSet<-interList.names[which(PPI.df[,"PPI"]==1)]
	gsc<-list(interSet)
	names(gsc)[1]<-"PPI"
	interList<-interList[sort.list(interList, decreasing=TRUE)]
	temp<-gseaScores(interList, interSet, exponent=1, mode="graph")		
	par(mar=c(.5, .5, .5,0))
	gseaPlots(temp$runningScore, temp$enrichmentScore, temp$positions, interList)
}
##filtering modules
Mulder2012.PANmoduleSelect2<-function(pan, mod.pval.cutoff=0.05, 
	mod.size.cutoff=4, avg.degree.cutoff=0.5) {
	##pvclust p-values
	summ<-list()
	mod<-pan@modules$clusters
	names(mod)<-1:length(mod)
	summ[["All clusters"]]<-length(mod)
	
	mod<-mod[pan@modules$pval<mod.pval.cutoff]
	summ[[paste("p<", mod.pval.cutoff, sep="")]]<-length(mod)
	
	##module size
	mod.size<-unlist(lapply(mod, length))
	mod<-mod[mod.size>mod.size.cutoff & mod.size<nrow(pan@bm1@pheno)/2]
	summ[["module size constraint"]]<-length(mod)
	
	##positive effects
	mod.graph<-lapply(mod, function(x) {
				pan@graph[x, x]		
			})
	mod.g.effects<-lapply(mod.graph, function(x) {
				ids<-rownames(x)
				all(apply(pan@bm1@pheno[ids, ], 1, mean)>0)
			})
	mod<-mod[which(unlist(mod.g.effects))]
	summ[["loss-of-function constraint"]]<-length(mod)
	
	avgdg<-NULL
	if(length(mod)>0) {
		##degree
		mod.graph<-mod.graph[which(unlist(mod.g.effects))]
		mod.graph<-lapply(mod.graph, function(x) sign(x))
		avgdg<-unlist(lapply(mod.graph, function(x) {
							mean(rowSums(x))/nrow(x)		
						}))
		inds<-which(avgdg>avg.degree.cutoff)
		mod<-mod[inds]
		avgdg<-avgdg[inds]
		summ[["module density constraint"]]<-length(mod)
	} else {
		summ[["module density constraint"]]<-0
	}
	
	moduleRef<-NULL
	isParents<-NULL
	globalMods<-NULL
	if(length(mod)>0) {
		##
		x<-NULL
		for(i in 1:length(mod)) {
			y<-sapply(setdiff(1:length(mod), i), function(j) 
				all(mod[[i]]%in%mod[[j]]))
			if(!(any(unlist(y)))) {
				x<-c(x, i)
			}
		}
		
		mod<-mod[sort.list(unlist(lapply(mod, length)), decreasing=TRUE)]
		
		moduleRef<-list()
		isParents<-list()
		globalMods<-c(1)		
		for(i in 1:length(mod)) {
			##find parents
			if(i>1) {
				isParents[[i]]<-sapply(1:(i-1), function(pt) {
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
	return(list(mod=mod, moduleRef=moduleRef, isParents=isParents, 
		globalMods=globalMods, summary=summ, avgdg=avgdg))
}
##filtering modules
Mulder2012.PANmoduleSelect<-function(pan, mod.pval.cutoff=0.05, 
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
								all(apply(pan@bm1@pheno[ids, ], 1, mean)>0)
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
						} else {
							summ[["loss-of-function constraint"]]<-0
							summ[["root modules"]]<-0	
						}
					} 
				} else {	
						summ[["module density constraint"]]<-0
						summ[["loss-of-function constraint"]]<-0
						summ[["root modules"]]<-0							
				}
			} 
		} else {
			summ[["module size constraint"]]<-0
			summ[["module density constraint"]]<-0
			summ[["loss-of-function constraint"]]<-0
			summ[["root modules"]]<-0					
		}
	} else {
		summ[[paste("p<", mod.pval.cutoff, sep="")]]<-0
		summ[["module size constraint"]]<-0
		summ[["loss-of-function constraint"]]<-0
		summ[["module density constraint"]]<-0
		summ[["root modules"]]<-0	
	}
	return(list(mod=mod, moduleRef=moduleRef, isParents=isParents, 
		globalMods=globalMods, summary=summ, avgdg=avgdg))
}
Mulder2012.module.visualize<-function(rdp, pan, mod.pval.cutoff=0.05, 
	mod.size.cutoff=4, avg.degree.cutoff=0.5, edgeWidthLeg=TRUE, 
	filter.effects=TRUE) {

	mod<-Mulder2012.PANmoduleSelect(pan, mod.pval.cutoff, mod.size.cutoff, 
		avg.degree.cutoff, filter.effects)
	globalMods<-mod$globalMods
	moduleRef<-mod$moduleRef
	isParents<-mod$isParents
	mod<-mod$mod

	gigRef<-addGraph(rdp, subgraph(pan@iPAN, unique(unlist(mod))), zoom=60)
	iroot<-0
	for(i in 1:length(mod)) {
		ic<-(iroot)%%2+1
		ir<-floor(iroot/2)+1
		ml<-length(mod[[i]])
		if(i>1 && sum(isParents[[i]])>0) {							
			closestParent<-rev(which(isParents[[i]]))[1]
			moduleRef[[i]]<-nestNodes(rdp, nodes=mod[[i]], parent=
				moduleRef[[closestParent]], gcoord=c(50, 50),
				isAnchor=TRUE, gscale=70, theme="tm1", gatt=list(nestFontSize=22, 
				nestAlias=paste("p-value:",as.character(round(pan@modules$pval[
				as.numeric(names(mod)[i])], 3))),sep=""))
		} else {
			moduleRef[[i]]<-nestNodes(rdp, nodes=mod[[i]], parent=gigRef, 
				gcoord=c((ic-1)*40+20, (ir-1)*40+20), isAnchor=TRUE, gscale=40, 
				theme="tm1", gatt=list(nestFontSize=22, nestAlias=paste(
				"p-value:",as.character(round(pan@modules$pval[as.numeric(
				names(mod)[i])], 3))),sep=""))
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
	addLegend.size(rdp, sizevec=round(pan@legend$nodeSize[1:10*floor(
		nrow(pan@legend$nodeSize)/10), "diameter"], 2), position="bottomright",
		labvec=pan@legend$nodeSize[1:10*floor(nrow(pan@legend$nodeSize)/10), 
		"degree"], type="node", title="Node degree")
	if(edgeWidthLeg)
		addLegend.size(rdp, sizevec=round(pan@legend$edgeWidth[, "width"], 2), 
			position="bottomleft", labvec=round(pan@legend$edgeWidth[, "SNR"], 2), 
			title="Log(SNR)", type="edge")
	mergeOutEdges(rdp, lb=1, ub=15, nlev=3)
	relax(rdp)
}
##distance function for pvclust
dist.pvclust4PAN <- function(x, method="euclidean", use.cor="pairwise.complete.obs")
{
	if(!is.na(pmatch(method,"correlation"))){
		res <- as.dist(1 - cor(x, method="pearson", use=use.cor))
		attr(res,"method") <- "correlation"
		return(res)
	}
	else if(!is.na(pmatch(method,"abscor"))){
		res <- as.dist(1 - abs(cor(x,method="pearson",use=use.cor)))
		attr(res,"method") <- "abscor"
		return(res)
	}
	else if(!is.na(pmatch(method,"uncentered"))){
		if(sum(is.na(x)) > 0){
			x <- na.omit(x)
			warning("Rows including NAs were omitted")
		}
		x  <- as.matrix(x)
		P  <- crossprod(x)
		qq <- matrix(diag(P),ncol=ncol(P))
		Q  <- sqrt(crossprod(qq))
		res <- as.dist(1 - P/Q)
		attr(res,"method") <- "uncentered"
		return(res)
	}
	else if(!is.na(pmatch(method, "cosine2"))){
		res<-cosineDist(1-as.matrix(cosineDist(x)))
		attr(res, "method")<-"cosine2"
		return(res)
	}
	else if(!is.na(pmatch(method, "cosine"))){
		res<-cosineDist(x)
		attr(res, "method")<-"cosine"
		return(res)
	}
	else
		dist(t(x),method)
}
##parallel computing for pvclust
parPvclust4PAN <- function(cl, data, method.hclust="average",
		method.dist="correlation", use.cor="pairwise.complete.obs",
		nboot=1000, r=seq(.5,1.4,by=.1), store=FALSE, weight=FALSE, 
		init.rand=TRUE, seed=NULL)
{
	if(!(require(snow))) stop("Package snow is required for parPvclust.")
	
	if((ncl <- length(cl)) < 2 || ncl > nboot) {
		warning("Too small value for nboot: non-parallel version is executed.")
		return(pvclust(data,method.hclust,method.dist,use.cor,nboot,r,store))
	}
	
	if(init.rand) {
		if(is.null(seed))
			seed <- 1:length(cl)
		else if(length(seed) != length(cl))
			stop("seed and cl should have the same length.")
		
		# setting random seeds
		parLapply(cl, as.list(seed), set.seed)
	}
	
	# data: (n,p) matrix, n-samples, p-variables
	n <- nrow(data); p <- ncol(data)
	
	# hclust for original data
	METHODS <- c("ward", "single", "complete", "average", "mcquitty", 
			"median", "centroid", "cosine", "cosine2", "TO")
	method.hclust <- METHODS[pmatch(method.hclust, METHODS)]
	##!patch
	dist.pvclust<-getFromNamespace("dist.pvclust", ns=getNamespace("pvclust"))
	pvclust.node<-getFromNamespace("pvclust.node", ns=getNamespace("pvclust"))
	pvclust.merge<-getFromNamespace("pvclust.merge", ns=getNamespace("pvclust"))
	boot.hclust<-getFromNamespace("boot.hclust", ns=getNamespace("pvclust"))
	
	distance <- dist.pvclust(data, method=method.dist, use.cor=use.cor)
	data.hclust <- hclust(distance, method=method.hclust)
	
	# multiscale bootstrap
	size <- floor(n*r)
	rl <- length(size)
	
	if(rl == 1) {
		if(r != 1.0)
			warning("Relative sample size r is set to 1.0. AU p-values are not 
			calculated\n")
		
		r <- list(1.0)
	}
	else
		r <- as.list(size/n)
	
	nbl <- as.list(rep(nboot %/% ncl,times=ncl))
	
	if((rem <- nboot %% ncl) > 0)
		nbl[1:rem] <- lapply(nbl[1:rem], "+", 1)
	
	cat("Multiscale bootstrap... ")
	
	##!patch
	
	mlist <- parLapply(cl, nbl, pvclust.node,
			r=r, data=data, object.hclust=data.hclust, method.dist=method.dist,
			use.cor=use.cor, method.hclust=method.hclust,
			store=store, weight=weight)
	cat("Done.\n")
	
	mboot <- mlist[[1]]
	
	for(i in 2:ncl) {
		for(j in 1:rl) {
			mboot[[j]]$edges.cnt <- mboot[[j]]$edges.cnt + mlist[[i]][[j]]$edges.cnt
			mboot[[j]]$nboot <- mboot[[j]]$nboot + mlist[[i]][[j]]$nboot
			mboot[[j]]$store <- c(mboot[[j]]$store, mlist[[i]][[j]]$store)
		}
	}
	
	result <- pvclust.merge( data=data, object.hclust=data.hclust, mboot=mboot)
	
	return(result)
}
