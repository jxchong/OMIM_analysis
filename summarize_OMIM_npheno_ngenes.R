# For the file morbidmap, the fields are, in order:
# 1  - Disorder, <disorder MIM no.> (<phene mapping key>)
# 2  - Gene/locus symbols
# 3  - Gene/locus MIM no.
# 4  - cytogenetic location

# Phenotype mapping method - appears in parentheses after a disorder:
# ? - A question mark, "?", before the disease name indicates an unconfirmed or possibly spurious mapping.
# 1 - the disorder is placed on the map based on its association with
# a gene, but the underlying defect is not known.
# 2 - the disorder has been placed on the map by linkage; no mutation has
# been found.
# 3 - the molecular basis for the disorder is known; a mutation has been
# found in the gene.
# 4 - a contiguous gene deletion or duplication syndrome, multiple genes
# are deleted or duplicated causing the phenotype.

# An asterisk (*) before an entry number indicates a gene.
#
# A number symbol (#) before an entry number indicates that it is a descriptive entry, usually of a phenotype, and does not represent a unique locus. The reason for the use of the number symbol is given in the first paragraph of the entry. Discussion of any gene(s) related to the phenotype resides in another entry(ies) as described in the first paragraph.
#
# A plus sign (+) before an entry number indicates that the entry contains the description of a gene of known sequence and a phenotype.
#
# A percent sign (%) before an entry number indicates that the entry describes a confirmed mendelian phenotype or phenotypic locus for which the underlying molecular basis is not known.
#
# No symbol before an entry number generally indicates a description of a phenotype for which the mendelian basis, although suspected, has not been clearly established or that the separateness of this phenotype from that in another entry is unclear.
#
# A caret (^) before an entry number means the entry no longer exists because it was removed from the database or moved to another entry as indicated.
#
# Brackets, "[ ]", indicate "nondiseases," mainly genetic variations that lead to apparently abnormal laboratory test values (e.g., dysalbuminemic euthyroidal hyperthyroxinemia).
#
# Braces, "{ }", indicate mutations that contribute to susceptibility to multifactorial disorders (e.g., diabetes, asthma) or to susceptibility to infection (e.g., malaria).


library(RColorBrewer)
colors <- brewer.pal(12, "Set3")


morbidmap <- read.table("morbidmap.parsed.txt", head=TRUE, sep="\t", quote="")
mim2gene <- read.table("raw_download_2014-02-20/mim2gene.txt", head=TRUE, sep="\t", comment.char="", quote="", col.names=c("MIMnum", "Type", "Gene.IDs", "Approved.Gene.Symbols"))
mim2gene.morbidmap <- merge(subset(mim2gene, subset=Type!="moved/removed"&Type!="gene"&Type!="gene/phenotype", select="MIMnum"), morbidmap, by.x="MIMnum", by.y="phenoMIMnum", all.x=TRUE)
phenonames <- read.table("omimtxtZ.parsed.txt", head=TRUE, sep="\t", comment.char="", quote="")


# How many are non-Mendelian (QTL/somatic/susceptibility/etc)
nrow(subset(morbidmap, subset=isComplex=="yes", select="phenoMIMnum"))
# [1] 1492
sum(is.na(subset(morbidmap, subset=isComplex=="yes", select="phenoMIMnum")))
# [1] 564	# there are many phenotypes with no MIM number for whatever reason, so we need to count these as well; from a manual inspection, we have to assume these are unique
sum(!is.na(subset(morbidmap, subset=isComplex=="yes", select="phenoMIMnum")))
# [1] 928
	length(unique(na.omit(subset(morbidmap, subset=isComplex=="yes")$phenoMIMnum)))
	# [1] 486	# how many unique MIM numbers are there
sum(is.na(subset(morbidmap, subset=isComplex=="yes", select="phenoMIMnum")))+length(unique(na.omit(subset(morbidmap, subset=isComplex=="yes")$phenoMIMnum)))
# [1050]	# total non-Mendelian phenotypes



# How many of these are chromosomal dup/dels (mappingkey=4)
nrow(subset(morbidmap, subset=isComplex=="no"&phenoMappingKey==4))
nrow(subset(mim2gene.morbidmap, subset=isComplex=="no"&phenoMappingKey==4))
# [1] 121


# How many are Mendelian and "monogenic" phenotypes (exclude chrom syndromes and complex) either mapped or unmapped to the genome
# Again, if the phenotype is NOT complex or a chromosomal syndrome and there is no mapping key value, we have to assume that the phenotype IS Mendelian and monogenic because we don't have time to process the raw OMIM text and get the name of every phenotype
mim2gene.morbidmap.nocomplex.nophenokey4 <- subset(mim2gene.morbidmap, subset=(isComplex=="no"&phenoMappingKey!=4)|is.na(phenoMappingKey))
mim2gene.morbidmap.nocomplex.nophenokey4.uniquelocusphenopair.bool <- !duplicated(subset(mim2gene.morbidmap.nocomplex.nophenokey4, select=c(MIMnum,GeneMIMnum)))
mim2gene.morbidmap.nocomplex.nophenokey4.nodupelocusphenopair.bool <- !duplicated(subset(mim2gene.morbidmap.nocomplex.nophenokey4, select=c(MIMnum,GeneMIMnum)), fromLast=TRUE)&!duplicated(subset(mim2gene.morbidmap.nocomplex.nophenokey4, select=c(MIMnum,GeneMIMnum)), fromLast=FALSE)
mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair <- mim2gene.morbidmap.nocomplex.nophenokey4[mim2gene.morbidmap.nocomplex.nophenokey4.uniquelocusphenopair.bool,]
# write.table(mim2gene.morbidmap.nocomplex.nophenokey4[mim2gene.morbidmap.nocomplex.nophenokey4.uniquelocusphenopair.bool,], file="OMIM.mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.txt", quote=FALSE, row.names=FALSE, sep="\t")


# some of these phenotypes don't have names that show that they are actually QTLs; you need to get the name of the phenotype MIM number from omim.txt.Z for it to be obvious that it's a QTL
mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.phenonames <- merge(mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair, phenonames, by.y="MIMnum", by.x="MIMnum", all.x=TRUE, suffixes=c(".morbidmap",".omimtxt"))
mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.phenonames.nocomplex <- subset(mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.phenonames, subset=isComplex.omimtxt=="no")

write.table(mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.phenonames.nocomplex, file="OMIM.phenotypecategories.all.txt", quote=FALSE, row.names=FALSE, sep="\t")

mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.phenonames.nocomplex.mappingkey3 <- subset(mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.phenonames.nocomplex, subset=phenoMappingKey==3)
mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.phenonames.nocomplex.mappingkeynot3 <- subset(mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.phenonames.nocomplex, subset=phenoMappingKey!=3|is.na(phenoMappingKey))

printcolumns <- c("MIMnum","phenoname","phenoMappingKey","LocusSymbols","GeneMIMnum","mimtitle","MIMlink")
write.table(subset(mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.phenonames.nocomplex.mappingkeynot3, select=printcolumns), file="OMIM.phenotypecategories.unexplained.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(subset(mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.phenonames.nocomplex.mappingkey3, select=printcolumns), file="OMIM.phenotypecategories.explained.txt", quote=FALSE, row.names=FALSE, sep="\t")


##################

# How many genes for mapped phenotypes (mapping key == 3)
nrow(unique(subset(mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.phenonames.nocomplex.mappingkey3, select="LocusSymbols")))
# [1] 2822	


# How many are Mendelian and "monogenic" (exclude chrom syndromes and complex) phenotypes are mapped to a general locus (mapping key = 2)
# looks like OMIM usually doesn't assign a phenotype MIM number to phenotypes with mapping key=2 (686 with no phenotype MIM number but OMIM gives them a gene MIM number that corresponds to both the gene and the phenotype)
length(unique(subset(morbidmap, subset=isComplex=="no"&phenoMappingKey==2)$phenoname))
# [1] 696
length(unique(subset(morbidmap, subset=isComplex=="no"&phenoMappingKey==2)$GeneMIMnum))
# [1] 701




# How many are Mendelian and "monogenic" (exclude chrom syndromes and complex) phenotypes are mapped to a specific gene (mapping key = 3)
nrow(subset(mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.phenonames.nocomplex, subset=phenoMappingKey==3))
# write.table(subset(mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.phenonames.nocomplex, subset=phenoMappingKey!=3), file="OMIMpheno.mappingnot3.txt", quote=FALSE, row.names=FALSE, sep="\t")
# [1] 5236


morbidmap.3.nocomplex.phenoMIM <- subset(mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.phenonames.nocomplex, subset=phenoMappingKey==3&!is.na(MIMnum))
# write.table(morbidmap.3.nocomplex.phenoMIM, file="OMIMpheno.mapping3.notcomplex.isnotNA.txt", quote=FALSE, row.names=FALSE, sep="\t")
nrow(morbidmap.3.nocomplex.phenoMIM)
# [1] 4116 have a phenotype MIM number
# write.table(subset(mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.phenonames.nocomplex, subset=isComplex=="yes"|phenoMappingKey!=3|is.na(phenoMIMnum)), file="OMIMpheno.mappingnot3.iscomplex.isNA.txt", quote=FALSE, row.names=FALSE, sep="\t")
# write.table(subset(mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.phenonames.nocomplex, subset=isComplex=="yes"&phenoMappingKey==3), file="OMIMpheno.mapping3.iscomplex.txt", quote=FALSE, row.names=FALSE, sep="\t")
# [1] 138 have no phenotype MIM number; I inspected a few of these and it seems like the phenotype is not "real" (should not be in OMIM to begin with)
# write.table(subset(mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.phenonames.nocomplex, subset=isComplex=="no"&phenoMappingKey==3&is.na(phenoMIMnum)), file="OMIMpheno.mapping3.notcomplex.isNA.txt", quote=FALSE, row.names=FALSE, sep="\t")


# nrow(unique(subset(morbidmap.3.nocomplex.phenoMIM, select=phenoMIMnum)))
# # [1] 3593 have a phenotype MIM number
# morbidmap.3.nocomplex.phenoMIM.uniquepheno.bool <- !duplicated(subset(morbidmap.3.nocomplex.phenoMIM, select=phenoMIMnum))
# morbidmap.3.nocomplex.phenoMIM.dupepheno.bool <- duplicated(subset(morbidmap.3.nocomplex.phenoMIM, select=phenoMIMnum), fromLast=FALSE)|duplicated(subset(morbidmap.3.nocomplex.phenoMIM, select=phenoMIMnum), fromLast=TRUE)
# write.table(morbidmap.3.nocomplex.phenoMIM[morbidmap.3.nocomplex.phenoMIM.uniquepheno.bool,], file="OMIMpheno.mapping3.notcomplex.isnotNA.uniquepheno.txt", quote=FALSE, row.names=FALSE, sep="\t")
# write.table(morbidmap.3.nocomplex.phenoMIM[morbidmap.3.nocomplex.phenoMIM.dupepheno.bool,], file="OMIMpheno.mapping3.notcomplex.isnotNA.notuniquepheno.txt", quote=FALSE, row.names=FALSE, sep="\t")
#

# How many unique phenoMIMnum+GeneMIMnum combinations are there This might be the true number of unique phenotypes
nrow(unique(subset(morbidmap.3.nocomplex.phenoMIM, select=c(MIMnum,GeneMIMnum))))
morbidmap.3.nocomplex.phenoMIM.uniquelocusphenopair.bool <- !duplicated(subset(morbidmap.3.nocomplex.phenoMIM, select=c(MIMnum,GeneMIMnum)))
morbidmap.3.nocomplex.phenoMIM.nodupelocusphenopair.bool <- !duplicated(subset(morbidmap.3.nocomplex.phenoMIM, select=c(MIMnum,GeneMIMnum)), fromLast=TRUE)&!duplicated(subset(morbidmap.3.nocomplex.phenoMIM, select=c(MIMnum,GeneMIMnum)), fromLast=FALSE)
sum(morbidmap.3.nocomplex.phenoMIM.nodupelocusphenopair.bool)
# 3749 completely unique pheno+geneMIM combinations, i.e. values that appear twice are missing completely
# we would need to add one copy of each duplicated value back to get the total number of unique phenotypes
# write.table(morbidmap.3.nocomplex.phenoMIM[morbidmap.3.nocomplex.phenoMIM.uniquelocusphenopair.bool,], file="OMIMpheno.mapping3.notcomplex.isnotNA.uniquelocusphenopair.txt", quote=FALSE, row.names=FALSE, sep="\t")
# [1] 3923
# write.table(morbidmap.3.nocomplex.phenoMIM[!morbidmap.3.nocomplex.phenoMIM.nodupelocusphenopair.bool,], file="OMIMpheno.mapping3.notcomplex.isnotNA.nodupelocusphenopair.txt", quote=FALSE, row.names=FALSE, sep="\t")
# [1] 3749

#
# test1 <- morbidmap.3.nocomplex.phenoMIM[morbidmap.3.nocomplex.phenoMIM.uniquelocusphenopair.bool,]$phenoMIMnum
# test2 <- subset(mim2gene.morbidmap.notcomplex.nophenokey4.uniquelocusphenopair.phenonames.nocomplex, subset=phenoMappingKey==3)$MIMnum
# test1[!(test1 %in% test2)]


#
# # get all duplicated phenoMIMnum+GeneMIMnum combinations including first appearance
# morbidmap.3.nocomplex.phenoMIM.dupelocusphenopair.bool <- duplicated(subset(morbidmap.3.nocomplex.phenoMIM, select=c(MIMnum,GeneMIMnum)), fromLast=TRUE)|duplicated(subset(morbidmap.3.nocomplex.phenoMIM, select=c(MIMnum,GeneMIMnum)), fromLast=FALSE)
# sum(morbidmap.3.nocomplex.phenoMIM.dupelocusphenopair.bool)
# # [1] 367 phenoMIM+GeneMIM pairs appear more than once (could be even more than 2x)
# write.table(morbidmap.3.nocomplex.phenoMIM[morbidmap.3.nocomplex.phenoMIM.dupelocusphenopair.bool,], file="OMIMpheno.mapping3.notcomplex.isnotNA.NOTuniquelocusphenopair.txt", quote=FALSE, row.names=FALSE, sep="\t")
# # [1] 367
# nrow(unique(subset(morbidmap.3.nocomplex.phenoMIM[morbidmap.3.nocomplex.phenoMIM.dupelocusphenopair.bool,], select=c(MIMnum,GeneMIMnum))))
# # [1] 174 unique values within the set of dupes
#









d <- read.table("OMIM.phenotypecategories.explained.txt", head=TRUE, sep="\t", quote="", comment.char="")

nphenopergene <- table(table(d$GeneMIMnum))

pdf("OMIM number of genes vs npheno - barplot.2015-02-22.pdf", width=10, height=5)
x.pos <- barplot(nphenopergene, col=colors[10], xlab="number of Mendelian phenotypes", ylab="number of genes", border=NA)
text(x.pos, nphenopergene, labels=nphenopergene, xpd=TRUE, pos=3, cex=0.8, offset=0.1)
dev.off()


