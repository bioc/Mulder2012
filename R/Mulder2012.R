library(RedeR)
##pipeline function reproducing results for the Mulder2012 data set
Mulder2012.pipeline<-function(
	par4BM=list(model="global", metric="cosine", nPerm=20),
	par4PAN=list(type="SNR", log=TRUE, sign=TRUE, cutoff=log(10), filter=FALSE),
	par4ModuleSearch=list(nboot=10000, metric="cosine2", hclustMethod="average", 
	filter=FALSE)
) {
	##1. (a) PINdb data preprocessing
	PPIandModules<-Mulder2012.PPIPre()
	PPI<-PPIandModules$PPI
	proteinModules<-PPIandModules$proteinModules
	save(PPI, proteinModules, file=file.path("rslt", "Mulder2012.PPI.RData"))
	##1. (b) RNAi screens preprocessing
	Mulder2012<-Mulder2012.RNAiPre()
	cat("-data preprocessing finished!\n")
	save(Mulder2012, file=file.path("rslt", "Mulder2012.RData")) 
	##2. (a) Beta-mixture modelling
	bm_Mulder2012<-Mulder2012.BMfitting(pheno=Mulder2012, model=par4BM$model, 
			metric=par4BM$metric, nPerm=par4BM$nPerm)
	cat("-beta-mixture (global) model fitting finished!\n")
	save(bm_Mulder2012, file=file.path("rslt", "bm_Mulder2012.RData"))
	##2. (b) extended model
	overlap.genes<-intersect(rownames(PPI), rownames(Mulder2012))
	part<-matrix(0, nrow(Mulder2012), nrow(Mulder2012))
	dimnames(part)<-list(rownames(Mulder2012), rownames(Mulder2012))
	part[overlap.genes, overlap.genes]<-PPI[overlap.genes, overlap.genes]
	part<-part[upper.tri(part)]
	bm_ext_Mulder2012<-Mulder2012.BMfitting.extended(pheno=Mulder2012, 
		model="stratified", metric=par4BM$metric, nPerm=par4BM$nPerm, 
		partition=part)
	cat("-beta-mixture (extended) model fitting finished!\n")
	save(bm_ext_Mulder2012, file=file.path("rslt", "bm_ext_Mulder2012.RData"))
	##3. infer a PAN
	pan_Mulder2012<-Mulder2012.InferPAN(bm=bm_Mulder2012, type=par4PAN$type, 
		log=par4PAN$log, sign=par4PAN$sign, cutoff=par4PAN$cutoff, 
		filter=par4PAN$filter)
	cat("-inferring a PAN (global model) finished!\n")
	pan_Mulder2012<-buildPAN(pan_Mulder2012, engine="RedeR", 
		para=list(nodeSumCols=1:3, nodeSumMethod="average", hideNeg=TRUE), 
		verbose=TRUE)
	##4. search enriched modules
	pan_Mulder2012<-Mulder2012.ModuleSearchByPvclust(pan=pan_Mulder2012, 
		nboot=par4ModuleSearch$nboot, metric=par4ModuleSearch$metric, 
		hclustMethod=par4ModuleSearch$hclustMethod, filter=
		par4ModuleSearch$filter)	
	cat("-module searching (global model) finished!\n")
	save(pan_Mulder2012, file=file.path("rslt", "pan_Mulder2012.RData"))
	##3. infer a PAN
	pan_ext_Mulder2012<-Mulder2012.InferPAN(bm=bm_ext_Mulder2012, type=
		par4PAN$type, log=par4PAN$log, sign=par4PAN$sign, cutoff=par4PAN$cutoff, 
		filter=par4PAN$filter)
	pan_ext_Mulder2012<-buildPAN(pan_ext_Mulder2012, engine="RedeR", 
		para=list(nodeSumCols=1:3, nodeSumMethod="average", hideNeg=TRUE), 
		verbose=TRUE)
	cat("-inferring a PAN (extended model) finished!\n")
	##4. search for modules
	pan_ext_Mulder2012<-Mulder2012.ModuleSearchByPvclust(pan=pan_ext_Mulder2012, 
		nboot=par4ModuleSearch$nboot, metric=par4ModuleSearch$metric, 
		hclustMethod=par4ModuleSearch$hclustMethod, filter=
		par4ModuleSearch$filter)
	cat("-module searching (extended model) finished!\n")
	save(pan_ext_Mulder2012, file=file.path("rslt", "pan_ext_Mulder2012.RData"))
	##5. Enrichment analysis
	PPenrichInPPI<-Mulder2012.PPIenrich(Mulder2012, PPI, bm_Mulder2012)
	cat("-enrichment analysis finished!\n")
	save(PPenrichInPPI, file=file.path("rslt", "PPenrichInPPI.RData"))
}
##pipeline function reproducing figures for the Mulder2012 data set
Mulder2012.fig<-function(what="ALL") {
	##Figure 4. beta-mixture modelling
	##A NULL fitting
	if(what %in% c("NULLfitting", "ALL")) {
		data("bm_Mulder2012", package="Mulder2012")
		##load(file.path("rslt", "bm_Mulder2012.RData"))
		pdf(file=file.path("rslt", "fig4A.pdf"), width=7.874/2, height=7.874/2)
		view(bm_Mulder2012, "fitNULL")
		graphics.off()
		cat("Figure 4(A) generated!\n")
	}
	##B bm fitting
	if(what %in% c("BMfitting", "ALL")) {
		data("bm_Mulder2012", package="Mulder2012")
		##load(file.path("rslt", "bm_Mulder2012.RData"))
		pdf(file=file.path("rslt", "fig4B.pdf"), width=7.874/2, height=7.874/2)
		view(bm_Mulder2012, "fitBM")
		graphics.off()
		cat("Figure 4(B) generated!\n")
	}
	##Figure 5. enrichment analyses
	if(what %in% c("PPIenrich", "ALL")) {
		data("Mulder2012", package="Mulder2012")
		##load(file.path("rslt", "Mulder2012.RData"))
		data("bm_Mulder2012", package="Mulder2012")
		##load(file.path("rslt", "bm_Mulder2012.RData"))
		data("Mulder2012.PPI", package="Mulder2012")
		##load(file.path("rslt", "Mulder2012.PPI.RData"))
		labels<-c("A", "B", "C")
		names(labels)<-c("neg", "none", "pos")
		for(i in c("neg", "none", "pos")) {
			pdf(file=file.path("rslt", paste("fig5", labels[i], ".pdf",sep="")), 
				width=7.874, height=7.874*0.6)
			GSEARandomWalkFig(Mulder2012, PPI, bm_Mulder2012, i)
			graphics.off()
		}
		cat("Figure 5 generated!\n")
	}
	##Figure 6. beta-mixture modelling for extended model
	if(what %in% c("BMfitting.ext", "ALL")) {
		data("bm_ext_Mulder2012", package="Mulder2012")
		##load(file.path("rslt", "bm_ext_Mulder2012.RData"))
		pdf(file=file.path("rslt", "fig6.pdf"), width=7.874, height=7.874/2)
		view(bm_ext_Mulder2012, "fitBM")
		graphics.off()
		cat("Figure 6 generated!\n")
	}
	rdp<-NULL
	##Figure 7. biologically relevent significant modules
	if(what %in% c("sigMod", "ALL")) {
		if(is.null(rdp)) {
			rdp <- RedPort('MyPort')
			calld(rdp)
		} else {
			if(ping(rdp)==0)
				calld(rdp)
			else
				resetd(rdp)
		}
		##require(RedeR)
		system("sleep 5")
		data("pan_ext_Mulder2012", package="Mulder2012")
		##load(file.path("rslt", "pan_ext_Mulder2012.RData"))
		Mulder2012.module.visualize(rdp, pan_ext_Mulder2012, 
			mod.pval.cutoff=0.05, mod.size.cutoff=4, avg.degree.cutoff=0.5)
		cat("Figure 7 generated in RedeR. \nPlease manually improve the 
			layout in RedeR!\n")
		readline(prompt = "Please press any key to reset RedeR!\n")
		resetd(rdp)
	}	
	##Figure 8(A). selected module for validation
	if(what %in% c("selMod", "ALL")) {
		if(is.null(rdp)) {
			rdp <- RedPort('MyPort')
			calld(rdp)
		} else {
			if(ping(rdp)==0)
				calld(rdp)
			else
				resetd(rdp)
		}
		##require(RedeR)
		system("sleep 5")
		data("pan_ext_Mulder2012", package="Mulder2012")
		##load(file.path("rslt", "pan_ext_Mulder2012.RData"))
		addGraph(rdp, subgraph(pan_ext_Mulder2012@iPAN, c("ING5", "BPTF", 
			"BRD1", "UHRF1", "SMARCA5", "EZH2", "SMARCC2", "PRMT1")), zoom=60)
		relax(rdp)
		cat("Figure 8(A) generated in RedeR. \nPlease manually improve the 
			layout in RedeR!\n")
	}	
}
