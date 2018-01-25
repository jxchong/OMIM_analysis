# perl parse_morbidmap.pl --in raw_download_2015-04-02/morbidmap --out 2015-04-02.morbidmap.parsed.txt
# perl parse_omimtxtZ_count_NGS_year.pl --in raw_download_2015-04-02/omim.txt.Z --out 2015-04-02.omimtxtZ.parsed.NGS.year.txt
# perl combine_omimtxtZ_morbidmap_count_NGS_year.pl --morbidmap 2015-04-02.morbidmap.parsed.txt --mim2gene raw_download_2015-04-02/mim2gene.txt --omimtxt 2015-04-02.omimtxtZ.parsed.NGS.year.txt --out 2015-04-02.combinedOMIM.mentionsNGS.year.txt

# colorblind-friendly rainbow palette
#
# 77AADD (119,170,221) # blue
# 4477AA (68,119,170)	# darker
# 77CCCC (119,204,204) # teal
# 44AAAA (68,170,170)
# 88CCAA (136,204,170) # green
# 44AA77 (68,170,119)
# DDDD77 (221,221,119) # yellow
# AAAA44 (170,170,68)
# DDAA77 (221,170,119) # orange-brown
# AA7744 (170,119,68)
# DD7788 (221,119,136) # pink
# AA4455 (170,68,85)
# CC99BB (204,153,187) # mauve
# AA4488 (170,68,136)


reviewred <- "#E72F1C"
reviewblue <- "#05ADEE"


 
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
	stop("Need datestamp argument\n", call.=FALSE)
} else {
	currentdate <- paste0("", args[1])
}


# current download
current <- read.table(paste0(currentdate,".combinedOMIM.mentionsNGS.year.txt"), head=TRUE, sep="\t", stringsAsFactors=FALSE, comment.char="", quote="", strip.white=TRUE)




# library(RColorBrewer)
# colors <- brewer.pal(12, "Set3")
colors <- as.character(c("#77AADD", "#77CCCC", "#88CCAA", "#DDDD77", "#DDAA77", "#DD7788", "#CC99BB"))
colors.darker <- as.character(c("#4477AA", "#44AAAA", "#44AA77", "#AAAA44", "#AA7744", "#AA4455", "#AA4488"))


# read in parsed data
# download February 20, 2015
# feb <- read.table("2015-02-20.combinedOMIM.mentionsNGS.year.txt", head=TRUE, sep="\t", stringsAsFactors=FALSE, comment.char="", quote="", strip.white=TRUE)
# download April 2, 2015
# april <- read.table("2015-04-02.combinedOMIM.mentionsNGS.year.txt", head=TRUE, sep="\t", stringsAsFactors=FALSE, comment.char="", quote="", strip.white=TRUE)


# keep only the entries that are monogenic and have NGS (or synonymous) in the first paragraph of the Molecular Genetics section
# feb.NGS <- subset(feb, select=NGSyear, subset=phenomappingkey==3&isComplex!="yes"&isComplex!="cancer"&mentionsNGS=="yes"&grepl("'1/", NGSparagraph))
# april.NGS <- subset(april, select=NGSyear, subset=phenomappingkey==3&isComplex!="yes"&isComplex!="cancer"&mentionsNGS=="yes"&grepl("'1/", NGSparagraph))
current.NGS <- subset(current, select=NGSyear, subset=phenomappingkey==3&isComplex!="yes"&isComplex!="cancer"&mentionsNGS=="yes"&grepl("'1/", NGSparagraph))

# feb.noNGS <- subset(feb, select=yearDiscovered, subset=phenomappingkey==3&isComplex!="yes"&isComplex!="cancer"&!(mentionsNGS=="yes"&grepl("'1/", NGSparagraph))&yearDiscovered>=2010)
# april.noNGS <- subset(april, select=yearDiscovered, subset=phenomappingkey==3&isComplex!="yes"&isComplex!="cancer"&!(mentionsNGS=="yes"&grepl("'1/", NGSparagraph))&yearDiscovered>=2010)
current.noNGS <- subset(current, select=yearDiscovered, subset=phenomappingkey==3&isComplex!="yes"&isComplex!="cancer"&!(mentionsNGS=="yes"&grepl("'1/", NGSparagraph))&yearDiscovered>=2010)

# feb.noNGS.old <- table(subset(feb, select=yearDiscovered, subset=phenomappingkey==3&isComplex!="yes"&isComplex!="cancer"&!(mentionsNGS=="yes"&grepl("'1/", NGSparagraph))))
# april.noNGS.old <- table(subset(april, select=yearDiscovered, subset=phenomappingkey==3&isComplex!="yes"&isComplex!="cancer"&!(mentionsNGS=="yes"&grepl("'1/", NGSparagraph))))
current.noNGS.old <- table(subset(current, select=yearDiscovered, subset=phenomappingkey==3&isComplex!="yes"&isComplex!="cancer"&!(mentionsNGS=="yes"&grepl("'1/", NGSparagraph))))



# check difference in discoveries per year between Feb 20 and current 2
# table(current.NGS) - table(feb.NGS)
# result: 
# 2010 2011 2012 2013 2014 2015
#    0    0    1    5    7   15
# Note: this is not surprising given OMIM is manually curated, but makes things a little more iffy -- OMIM is adding to entries from previous years


# check counts per year
sum(table(current.NGS))
# 583
table(current.NGS)
# 2010 2011 2012 2013 2014 2015
#   12   55  139  180  171   26
# sum(table(feb.NGS))
# 555
# table(feb.NGS)
# 2010 2011 2012 2013 2014 2015
#   12   55  138  175  164   11


# data summaries for graphing
current.NGS.cumsum <- cumsum(as.vector(table(current.NGS)))
years <- names(table(current.NGS))
years.formatted <- years


# compare discoveries by NGS vs conventional from 2010 onwards
current.NGSvsnoNGS <- matrix(rbind(table(current.NGS)[years], table(current.noNGS)[years]), nrow=2)
# feb.NGSvsnoNGS <- matrix(rbind(table(feb.NGS)[years], table(feb.noNGS)[years]), nrow=2)


# show all years to compare NGS vs no NGS
# feb.NGSvsnoNGS.allyears <- as.matrix(rbind(feb.noNGS.old, table(feb.NGS)[names(feb.noNGS.old)]))
# feb.NGSvsnoNGS.allyears[is.na(feb.NGSvsnoNGS.allyears)] <- 0
# feb.NGSvsnoNGS.allyears <- feb.NGSvsnoNGS.allyears[,colnames(feb.NGSvsnoNGS.allyears)>=1986]
# feb.NGSvsnoNGS.NGSyears <- feb.NGSvsnoNGS.allyears[,feb.NGSvsnoNGS.allyears[2,]!=0]

current.NGSvsnoNGS.allyears <- as.matrix(rbind(current.noNGS.old, table(current.NGS)[names(current.noNGS.old)]))
current.NGSvsnoNGS.allyears[is.na(current.NGSvsnoNGS.allyears)] <- 0
current.NGSvsnoNGS.allyears <- current.NGSvsnoNGS.allyears[,colnames(current.NGSvsnoNGS.allyears)>=1986]
current.NGSvsnoNGS.NGSyears <- current.NGSvsnoNGS.allyears[,current.NGSvsnoNGS.allyears[2,]!=0]




pdf(paste0(currentdate,".NGS_discoveries_per_year.barplot.pdf"), width=5, height=5)
barplot.x <- barplot(table(current.NGS), names.arg=NA, ylim=c(0,200), col=colors.darker[3], xlab="Year", ylab="Number of gene discoveries by NGS", border=NA, cex.main=0.95)
mtext(years.formatted, at=barplot.x, side=1, cex=0.9)
text(x=barplot.x, y=as.vector(table(current.NGS)), labels=as.vector(table(current.NGS)), pos=3, xpd=TRUE, cex=0.9)
dev.off()

# pdf(paste0(currentdate,".NGS_discoveries_per_year.scatter.pdf"), width=5, height=5)
# plot(names(table(current.NGS)), as.vector(table(current.NGS)), ylim=c(0,200), xlab="Year", ylab="Number of gene discoveries by NGS", xaxt="n",col=colors.darker[3], bg=colors.darker[3], pch=21, cex.main=0.95, bty="l")
# axis(labels=years.formatted, at=names(table(current.NGS)), side=1, cex=0.9)
# dev.off()
#
# pdf(paste0(currentdate,".NGS_discoveries_per_year.cumhistline.pdf"), width=5, height=5)
# plot(names(table(current.NGS)), current.NGS.cumsum, xlab="Year", ylab="Cumulative gene discoveries by NGS", xaxt="n", col=colors.darker[3], bg=colors.darker[3], pch=21, cex.main=0.95, cex.lab=0.95, type="b", bty="l")
# axis(labels=years.formatted, at=names(table(current.NGS)), side=1, cex=0.9)
# dev.off()
#
# pdf(paste0(currentdate,".NGS_discoveries_per_year.cumhistbar.pdf"), width=5, height=5)
# cumhistbar.x <- barplot(cumsum(as.vector(table(current.NGS))), names.arg=NA, xlab="Year", ylab="Cumulative gene discoveries by NGS", ylim=c(0,600), col=colors.darker[3], bg=colors.darker[3], pch=21, border=NA, cex.main=0.95)
# mtext(years.formatted, at=cumhistbar.x, side=1, cex=0.9)
# text(x=cumhistbar.x, y=current.NGS.cumsum, labels=current.NGS.cumsum, pos=3, xpd=TRUE, cex=0.9)
# dev.off()




pdf(paste0(currentdate,".NGSvsNoNGS_discoveries_per_year.pdf"), width=6, height=4)
par(mar=c(2.5,4,1,0))
bars.x <- barplot(current.NGSvsnoNGS, beside=TRUE, xaxt="n", ylab="Approximate # of gene discoveries by method", ylim=c(0,200), col=c(reviewred, reviewblue), pch=21, border=NA, cex.main=0.95, cex.axis=0.8)
mtext(years.formatted, at=colMeans(bars.x), side=1, cex=0.9)
mtext("Year", side=1, cex=0.9, line=1.5)
text(x=bars.x, y=current.NGSvsnoNGS, labels=current.NGSvsnoNGS, pos=3, xpd=TRUE, cex=0.9)
legend("topright", fill=c(reviewred, reviewblue), legend=c("WES/WGS", "conventional"), border=NA, bty="n", cex=0.8)
dev.off()


pdf(paste0(currentdate,".NoNGS_discoveries_since1986.pdf"), width=8, height=4)
par(mar=c(3,4,3,0))
bars.x <- barplot(current.NGSvsnoNGS.allyears, beside=FALSE, ylab="Approximate # of gene discoveries by method", ylim=c(0,300), col=c(reviewred, reviewblue), pch=21, border=NA, cex.main=0.95, cex.axis=0.8, xaxt="n")
text(x=bars.x, y=-20, labels=colnames(current.NGSvsnoNGS.allyears), cex=0.8, xpd=TRUE, srt=90)
mtext("Year", side=1, cex=0.8, line=2)
text(x=bars.x, y=colSums(current.NGSvsnoNGS.allyears), labels=colSums(current.NGSvsnoNGS.allyears), pos=3, xpd=TRUE, cex=0.5)
text(x=bars.x[current.NGSvsnoNGS.allyears[2,]!=0], y=current.NGSvsnoNGS.NGSyears[1,]+0.5*current.NGSvsnoNGS.NGSyears[2,], labels=current.NGSvsnoNGS.allyears[2,current.NGSvsnoNGS.allyears[2,]!=0], xpd=TRUE, cex=0.5, col="white")
legend("topleft", fill=c(reviewred, reviewblue), legend=c("conventional", "NGS"), border=NA, bty="n", cex=0.8)
dev.off()




# pdf("2015-02-20.NGSvsNoNGS_discoveries_per_year.pdf", width=6, height=4)
# par(mar=c(2.5,4,1,0))
# bars.x <- barplot(feb.NGSvsnoNGS, beside=TRUE, xaxt="n", ylab="Approximate # of gene discoveries by method", ylim=c(0,200), col=c(reviewred, reviewblue), pch=21, border=NA, cex.main=0.95, cex.axis=0.8)
# mtext(years.formatted, at=colMeans(bars.x), side=1, cex=0.9)
# mtext("Year", side=1, cex=0.9, line=1.5)
# text(x=bars.x, y=feb.NGSvsnoNGS, labels=feb.NGSvsnoNGS, pos=3, xpd=TRUE, cex=0.9)
# legend("topright", fill=c(reviewred, reviewblue), legend=c("WES/WGS", "conventional"), border=NA, bty="n", cex=0.8)
# dev.off()

# pdf("2015-02-20.NGS_discoveries_since1986.pdf", width=8, height=4)
# par(mar=c(3,4,3,0))
# bars.x <- barplot(feb.NGSvsnoNGS.allyears, beside=FALSE, ylab="Approximate # of gene discoveries by method", ylim=c(0,300), col=c(reviewred, reviewblue), pch=21, border=NA, cex.main=0.95, cex.axis=0.8, xaxt="n")
# text(x=bars.x, y=-20, labels=colnames(feb.NGSvsnoNGS.allyears), cex=0.8, xpd=TRUE, srt=90)
# mtext("Year", side=1, cex=0.8, line=2)
# text(x=bars.x, y=colSums(feb.NGSvsnoNGS.allyears), labels=colSums(feb.NGSvsnoNGS.allyears), pos=3, xpd=TRUE, cex=0.5)
# text(x=bars.x[feb.NGSvsnoNGS.allyears[2,]!=0], y=feb.NGSvsnoNGS.NGSyears[1,]+0.5*feb.NGSvsnoNGS.NGSyears[2,], labels=feb.NGSvsnoNGS.allyears[2,feb.NGSvsnoNGS.allyears[2,]!=0], xpd=TRUE, cex=0.5, col="white")
# legend("topleft", fill=c(reviewred, reviewblue), legend=c("conventional", "NGS"), border=NA, bty="n", cex=0.8)
# dev.off()
