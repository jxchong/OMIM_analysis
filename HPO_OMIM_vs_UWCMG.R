library(RColorBrewer)
colors = brewer.pal(12, "Set3")
# gray, green, blue, orange, purple
bestcolors=colors[c(9,7,6,5,10)]

# Core wrapping function
wrap.it <- function(x, len) { 
  sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

# Call this function with a list or vector
wrap.labels <- function(x, len) {
  if (is.list(x)) {
    lapply(x, wrap.it, len)
	} else {
    wrap.it(x, len)
	}
}


omim.raw <- read.table("OMIM.phenotypecategories.unexplained_HPO_annotated.txt", head=TRUE, sep="\t", comment.char="", quote="\"")
omim <- droplevels(subset(omim.raw, Comments==""|Comments=="to recruit", select=c("morbidmapname","phenomappingkey","isComplex","HPO_category","Comments")))
uwcmg <- read.table("UWCMG_phenotype_categories.txt", head=TRUE, sep="\t", comment.char="")

uwcmg.counts <- table(uwcmg$Single.Feature)
omim.counts <- table(omim$HPO_category)
omim.counts.recruit <- table(subset(omim, select=c(HPO_category,Comments)))

uwcmg.prop <- table(uwcmg$Single.Feature)/length(uwcmg$Single.Feature)
omim.prop <- table(omim$HPO_category)/length(omim$HPO_category)
omim.prop.recruit <- omim.counts.recruit[,"to recruit"]/rowSums(omim.counts.recruit)



write.table(uwcmg.prop, file="uwcmg.prop.txt", quote=FALSE, sep="\t")
write.table(omim.prop, file="omim.prop.txt", quote=FALSE, sep="\t")
write.table(omim.prop.recruit, file="recruit.prop.txt", quote=FALSE, sep="\t")

categories <- c("Abnormality of the abdomen", "Abnormality of blood and blood-forming tissues", "Abnormality of the breast", "Abnormality of the cardiovascular system", "Abnormality of connective tissue", "Abnormality of dentition", "Abnormality of the ear", "Abnormality of the endocrine system", "Abnormality of the eye", "Abnormality of the genitourinary system", "Abnormality of head or neck", "Abnormality of the immune system", "Abnormality of the integument", "Abnormality of limbs", "Abnormality of metabolism/homeostasis", "Abnormality of the musculature", "Neoplasm", "Abnormality of the nervous system", "Abnormality of the respiratory system", "Abnormality of the skeletal system", "Abnormality of the thoracic cavity", "non-syndromic ID/DD", "syndromic ID/DD", "Multiple congenital anomalies")

categories.ordered <- names(sort(omim.counts, decreasing=TRUE))


# check for consistency 
# unique(uwcmg$Single.Feature)[!(unique(uwcmg$Single.Feature) %in% omim$HPO_category)]


merged <- data.frame(matrix(0, nrow=length(categories.ordered), ncol=3), row.names=categories.ordered)
names(merged) <- c("UWCMG", "OMIM", "to recruit")

for (i in 1:length(categories.ordered)) {
	merged[categories.ordered[i],"UWCMG"] <- uwcmg.counts[categories.ordered[i]]
	merged[categories.ordered[i],"OMIM"] <- omim.counts[categories.ordered[i]]
	merged[categories.ordered[i],"to recruit"] <- omim.counts.recruit[categories.ordered[i],2]
}



# side by side barplot
layout(matrix(1:1))
par(xpd=T, mar=c(5.5,4,4,1))
x <- barplot(as.matrix(t(merged)), beside=TRUE, col=bestcolors[1:3], xaxt="n")
text(x=x[2,], y=rep(-10, times=nrow(merged)), labels=gsub("Abnormality of( the )*", "", row.names(merged)), xpd=TRUE, cex=0.7, srt=35, adj=c(1, NA))
legend("topleft", legend=names(merged), fill=bestcolors[1:3], bty="n")



merged.bycategory <- data.frame(matrix(0, nrow=2, ncol=length(categories.ordered)*2), row.names=c("OMIM not recruit","OMIM to recruit"))
names(merged.bycategory) <- rep(categories.ordered, each=2)

for (i in 1:length(categories.ordered)) {
	merged.bycategory[1, (i*2-1)] <- as.matrix(omim.counts.recruit)[categories.ordered[i], 1]
	merged.bycategory[2, (i*2-1)] <- as.matrix(omim.counts.recruit)[categories.ordered[i], "to recruit"]
	merged.bycategory[1, i*2] <- uwcmg.counts[categories.ordered[i]]
	merged.bycategory[2, i*2] <- 0
}
merged.bycategory[is.na(merged.bycategory)] <- 0
names(merged.bycategory)[names(merged.bycategory)=="Abnormality of blood and blood-forming tissues"] <- "Abnormality of blood/blood-forming tissues"
names(merged.bycategory)[names(merged.bycategory)=="Neoplasm"] <- "neoplasm"
names(merged.bycategory)[names(merged.bycategory)=="Multiple congenital anomalies"] <- "multiple congenital anomalies" #paste0("multiple\ncongenital anomalies")


colors.gradient <- colorRampPalette(c(bestcolors[2], "darkgrey"))

# stacked barplot for webpage, without UWCMG phase 1 and without recruiting
pdf("OMIM_unexplained.pdf", width=10, height=4)
layout(matrix(1:1, ncol=1))
par(xpd=T, mar=c(2,3.5,1,10), lwd=0.1)
ncategories.ordered <- length(categories.ordered)
odd <- seq(1,ncategories.ordered*2,2)
even <- seq(2,ncategories.ordered*2,2)
y.bars <- barplot(plot=FALSE, as.matrix(merged.bycategory)[,odd], beside=FALSE, yaxt="n", ylim=c(0,550), xpd=TRUE, xaxt="n")
y.labels <- y.bars
labels.categories <- gsub("Abnormality of (the )*", "", names(merged.bycategory))[odd]
omim.bars <- barplot(colSums(as.matrix(merged.bycategory)[,odd]), beside=FALSE, col=bestcolors[2], ylim=c(0,600), xpd=TRUE, xaxt="n", cex.axis=0.85)
# uwcmg.bars <- barplot(as.matrix(merged.bycategory)[,even], beside=FALSE, col=bestcolors[4], ylim=c(0,550), space=c(0,rep(2, times=ncategories.ordered-1)), xpd=TRUE, xaxt="n", yaxt="n", add=TRUE)
mtext("Phenotype categories", side=1, line=1, cex=0.85)
mtext("Number of unexplained-known phenotypes", side=2, line=2.5, cex=0.85)
text(x=y.labels, y=rep(-12, times=ncategories.ordered), labels=seq(1:length(labels.categories)), xpd=TRUE, cex=0.85)
# text(x=y.labels, y=rep(-10, times=ncategories.ordered), labels=labels.categories, xpd=TRUE, cex=0.7, adj=c(1, NA), srt=60)
legend(y.bars[ncategories.ordered]+.5, y=630, legend=paste(seq(1:length(labels.categories)), " ", labels.categories), bty="n", cex=0.8) 
dev.off()


# stacked barplot, without UWCMG phase 1
pdf("OMIM_unexplained.torecruit.pdf", width=4.5, height=4)
layout(matrix(1:1, ncol=1))
par(xpd=T, mar=c(6.7,3.5,1,0), lwd=0.1)
ncategories.ordered <- length(categories.ordered)
odd <- seq(1,ncategories.ordered*2,2)
even <- seq(2,ncategories.ordered*2,2)
y.bars <- barplot(plot=FALSE, as.matrix(merged.bycategory)[,odd], beside=FALSE, yaxt="n", ylim=c(0,550), xpd=TRUE, xaxt="n")
y.labels <- y.bars
labels.categories <- gsub("Abnormality of (the )*", "", names(merged.bycategory))[odd]
omim.bars <- barplot(as.matrix(merged.bycategory)[,odd], beside=FALSE, col=bestcolors[c(2,3)], ylim=c(0,550), xpd=TRUE, xaxt="n", cex.axis=0.8)
# uwcmg.bars <- barplot(as.matrix(merged.bycategory)[,even], beside=FALSE, col=bestcolors[4], ylim=c(0,550), space=c(0,rep(2, times=ncategories.ordered-1)), xpd=TRUE, xaxt="n", yaxt="n", add=TRUE)
mtext("Phenotype categories", side=1, line=5.5, cex=0.8)
mtext("Number of unexplained-known phenotypes", side=2, line=2.5, cex=0.8)
text(x=y.labels, y=rep(-10, times=ncategories.ordered), labels=labels.categories, xpd=TRUE, cex=0.7, adj=c(1, NA), srt=60)
legend("right", legend=c("OMIM unexplained-known phenotypes", "OMIM unexplained-known phenotypes, high priority"), fill=bestcolors[c(2,3)], bty="n", cex=0.7)
dev.off()


# scatterplot, UWCMG phase 1 vs OMIM
proportion.bycategory <- data.frame(matrix(0, nrow=length(categories.ordered), ncol=2), row.names=categories.ordered)
names(proportion.bycategory) <- c("OMIM", "UWCMG")
for (i in 1:ncategories.ordered) {
	proportion.bycategory[categories.ordered[i], "OMIM"] <- omim.prop[categories.ordered[i]]
	proportion.bycategory[categories.ordered[i], "UWCMG"] <- uwcmg.prop[categories.ordered[i]]
}
proportion.bycategory[is.na(proportion.bycategory)] <- 0
rownames(proportion.bycategory)[rownames(proportion.bycategory)=="Abnormality of blood and blood-forming tissues"] <- "Abnormality of blood/blood-forming tissues"
rownames(proportion.bycategory)[rownames(proportion.bycategory)=="Neoplasm"] <- "neoplasm"
rownames(proportion.bycategory)[rownames(proportion.bycategory)=="Multiple congenital anomalies"] <- "multiple congenital anomalies" #paste0("multiple\ncongenital anomalies")


pdf("OMIM_unexplained_vs_UWCMGphase1.pdf", width=5, height=5, useDingbats=FALSE)
par(mar=c(3.5,3.3,1,1), lwd=0.5)
limits <- c(0, 0.22)
plot(proportion.bycategory$OMIM, proportion.bycategory$UWCMG, xlim=limits, ylim=limits, pch=21, col=bestcolors[4], bg=bestcolors[4], cex=1.3, cex.axis=0.75, xlab="", ylab="")
mtext("Proportion of OMIM unexplained-known phenotypes", side=1, line=2.5, cex=0.75)
mtext("Proportion of UWCMG phase 1 phenotypes", side=2, line=2.5, cex=0.75)
abline(0,1, lty="dashed", lwd=0.5, col="gray")
text(proportion.bycategory$OMIM, proportion.bycategory$UWCMG+0.01, labels=gsub("Abnormality of (the )*", "", rownames(proportion.bycategory)), cex=0.7)
dev.off()



# # stacked barplot, with UWCMG phase 1
# pdf("OMIM_unexplained_with_UWCMGphase1.pdf", width=7.5, height=4)
# layout(matrix(1:1, ncol=1))
# par(xpd=T, mar=c(6.7,3.5,1,0), lwd=0.1)
# ncategories.ordered <- length(categories.ordered)
# odd <- seq(1,ncategories.ordered*2,2)
# even <- seq(2,ncategories.ordered*2,2)
# y.bars <- barplot(plot=FALSE, as.matrix(merged.bycategory)[,odd], beside=FALSE, yaxt="n", ylim=c(0,550), space=c(1,rep(2, times=ncategories.ordered-1)), xpd=TRUE, xaxt="n")
# y.labels <- y.bars-0.5
# labels.categories <- gsub("Abnormality of (the )*", "", names(merged.bycategory))[odd]
# omim.bars <- barplot(as.matrix(merged.bycategory)[,odd], beside=FALSE, col=bestcolors[c(2,3)], ylim=c(0,550), space=c(1,rep(2, times=ncategories.ordered-1)), xpd=TRUE, xaxt="n", cex.axis=0.8)
# uwcmg.bars <- barplot(as.matrix(merged.bycategory)[,even], beside=FALSE, col=bestcolors[4], ylim=c(0,550), space=c(0,rep(2, times=ncategories.ordered-1)), xpd=TRUE, xaxt="n", yaxt="n", add=TRUE)
# mtext("Phenotype categories", side=1, line=5.5, cex=0.8)
# mtext("Number of phenotypes", side=2, line=2.5, cex=0.8)
# text(x=y.labels, y=rep(-10, times=ncategories.ordered), labels=labels.categories, xpd=TRUE, cex=0.7, adj=c(1, NA), srt=55)
# legend("right", legend=c("UWCMG phase 1 phenotypes", "OMIM unexplained-known phenotypes", "OMIM unexplained-known phenotypes with high priority"), fill=bestcolors[c(4,2,3)], bty="n", cex=0.7)
# dev.off()
#

