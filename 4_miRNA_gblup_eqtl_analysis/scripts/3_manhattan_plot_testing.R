signames <-as.character(unique(rsumtb5$Gene))
signames
signames <- signames[c(7,8,9,11,13,14)]
signames

max(eqtl$gwa.pval[eqtl$gwa.qval < 0.05])
which(eqtl$gwa.qval == (max(eqtl$gwa.qval < 0.05)))
sum(eqtl$gwa.qval >= 0.05)
sum(eqtl$gwa.pval[eqtl$gwa.qval < 0.05])
summary(eqtl$gwa.pval[eqtl$gwa.qval >= 0.05])
sigthresh<-min(eqtl$gwa.pval[eqtl$gwa.qval >= 0.05])
sigthresh

par(mfrow = c(2,3))
par(cex = 0.6)
par(mar = c(0,0,0,0), oma = c(4,4,0.5,0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

for (i in 1:length(signames)){
manhattan_plot(eqtl$gwa.qval[,signames[[i]]], map = final_MSUPRP_miRNA$map, threshold=0.05, ylim = c(0,10))
}

manhattan_plot(eqtl$gwa.pval[,signames[[1]]], map = final_MSUPRP_miRNA$map, threshold=0.05, ylim = c(0,3))