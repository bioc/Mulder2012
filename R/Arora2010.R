##pipeline function reproducing results for the Arora2010 data set
Arora2010.pipeline<-function(
		par4BM=list(model="global", metric="cosine", nPerm=20),
		par4PAN=list(type="SNR", log=TRUE, sign=TRUE, cutoff=log(10), filter=FALSE),
		par4ModuleSearch=list(nboot=10000, metric="cosine2", hclustMethod="average", filter=FALSE)
		) {
	data("Arora2010", package="Mulder2012")
	##load(file.path("data", "Arora2010.RData"))
	##1. Beta-mixture modelling
	bm_Arora2010<-Arora2010.BMfitting(pheno=Arora2010, model=par4BM$model, 
			metric=par4BM$metric, nPerm=par4BM$nPerm)
	cat("-beta-mixture (global) model fitting finished!\n")
	save(bm_Arora2010, file=file.path("rslt", "bm_Arora2010.RData"))
	##2. infer a PAN
	pan_Arora2010<-Arora2010.InferPAN(bm=bm_Arora2010, type=par4PAN$type, log=par4PAN$log, 
			sign=par4PAN$sign, cutoff=par4PAN$cutoff, filter=par4PAN$filter)
	cat("-inferring a PAN finished!\n")
	##3. search for modules
	pan_Arora2010<-Arora2010.ModuleSearchByPvclust(pan=pan_Arora2010, nboot=par4ModuleSearch$nboot, 
			metric=par4ModuleSearch$metric, hclustMethod=par4ModuleSearch$hclustMethod, 
			filter=par4ModuleSearch$filter)
	save(pan_Arora2010, file=file.path("rslt", "pan_Arora2010.RData"))
	cat("-module searching finished!\n")
	##4. pathway analysis
	pw.rslt<-Arora2010.hypergeo(pan_Arora2010, mod.pval.cutoff=0.05, mod.size.cutoff=4, avg.degree.cutoff=0.5)
	save(pw.rslt, file=file.path("rslt", "Arora2010.pathway.RData"))
	cat("-pathway analysis finished!\n")
}
##pipeline function reproducing figures for the Arora2010 data set
Arora2010.fig<-function(what="ALL") {
	##A. NULL Fitting 
	if(what %in% c("NULLfitting", "ALL")) {
		data("bm_Arora2010", package="Mulder2012")
		##load(file.path("rslt", "bm_Arora2010.RData"))
		pdf(file=file.path("rslt", "fig3A.pdf"), width=7.874/2, height=7.874/2,  
				family=c("font/arial.afm", 
						"font/arialbd.afm", 
						"font/ariali.afm", 
						"font/arialbi.afm"))
#		par(cex.lab=1.7, cex.axis=1.7, mar=c(5, 5, 4, 2))
		view(bm_Arora2010, "fitNULL")
		graphics.off()
		cat("Figure 3(A) generated!\n")
	}
	##B. BM Fitting
	if(what %in% c("BMfitting", "ALL")) {
		data("bm_Arora2010", package="Mulder2012")
		##load(file.path("rslt", "bm_Arora2010.RData"))
		pdf(file=file.path("rslt", "fig3B.pdf"), width=7.874/2, height=7.874/2,  
				family=c("font/arial.afm", 
						"font/arialbd.afm", 
						"font/ariali.afm", 
						"font/arialbi.afm"))
#		par(cex.lab=1.7, cex.axis=1.7, mar=c(5, 5, 4, 2))
		view(bm_Arora2010, "fitBM")
		graphics.off()
		cat("Figure 3(B) generated!\n")
	}
	##C. Module visualization
	if(what %in% c("sigMod", "ALL")) {
		data("pan_Arora2010", package="Mulder2012")
		##load(file.path("rslt", "pan_Arora2010.RData"))
		rdp<-NULL
		if(is.null(rdp)) {
			rdp <- RedPort('MyPort')
			calld(rdp)
		} else {
			if(ping(rdp)==0)
				calld(rdp)
			else
				resetd(rdp)
		}		
		Arora2010.module.visualize(rdp, pan_Arora2010, mod.pval.cutoff=0.05, mod.size.cutoff=4, avg.degree.cutoff=0.5)
		cat("Figure 3(C) generated in RedeR. \nPlease manually improve the layout in RedeR!\n")
	}
	##D. Pathway analysis results
	if(what %in% c("pathway", "ALL")) {
		data("Arora2010.pathway", package="Mulder2012")
		##load(file.path("rslt", "Arora2010.pathway.RData"))
		pdf(file=file.path("rslt", "fig3D.pdf"), width=7.874*0.7, height=7.874*0.5,  
				family=c("font/arial.afm", 
						"font/arialbd.afm", 
						"font/ariali.afm", 
						"font/arialbi.afm"))
		pw.rslt<-pw.rslt[[1]]
		obs.exp<-as.numeric(pw.rslt[, 4])
		names(obs.exp)<-paste(as.character(pw.rslt[, 6]), "  (", format(pw.rslt[, 5], scientific=TRUE, digits=3), ")", sep="")
#		par(mar=c(6, 22, 4, 4))
		par(mar=c(4, 12, 1, 1))
#		barplot((obs.exp), horiz=TRUE, las=2, xlab="Observed/Expected Hits", cex.lab=1.6, cex.axis=1.6)
		barplot((obs.exp), horiz=TRUE, las=2, xlab="Observed/Expected Hits", cex.axis=0.6, cex.lab=0.6, cex=0.6)
		graphics.off()
		cat("Figure 3(D) generated!\n")
	}
}
