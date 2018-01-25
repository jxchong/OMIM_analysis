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

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
	stop("Need datestamp argument\n", call.=FALSE)
} else {
	datestamp <- paste0("", args[1])
}

# datestamp <- "../2016-04-11"

library(RColorBrewer)
colors <- brewer.pal(12, "Set3")


combined <- read.table(paste0(datestamp,".combinedOMIM.mentionsNGS.year.txt"), head=TRUE, sep="\t", comment.char="", quote="")
# nrow(combined)
cat(nrow(combined), "phenotypes", "\n")
# [1] 8776 phenotypes

combined.iscomplex <- subset(combined, isComplex=="yes")
# nrow(combined.iscomplex)
cat(nrow(combined.iscomplex), "are complex phenotypes", "\n")
# [1] 1217 are complex phenotypes

combined.iscancer <- subset(combined, isComplex=="cancer")
# nrow(combined.iscancer)
# [1] 121 are cancer phenotypes

combined.isnotcomplex <- subset(combined, isComplex!="yes"&isComplex!="cancer")
# nrow(combined.isnotcomplex)
# [1] 7438 are not complex and not cancer
cat(nrow(combined.isnotcomplex), "are not complex and not cancer", "\n")


combined.ischromosomal <- subset(combined, phenomappingkey==4)
# nrow(combined.ischromosomal)
cat(nrow(combined.ischromosomal), "are chromosomal", "\n")
# [1] 125

combined.nocomplex.nochrom <- subset(combined.isnotcomplex, phenomappingkey!=4)
# nrow(combined.nocomplex.nochrom)
cat(nrow(combined.nocomplex.nochrom), "are not complex or chromosomal", "\n")
# [1] 7315

combined.nocomplex.nochrom.3 <- subset(combined.nocomplex.nochrom, phenomappingkey==3)
# nrow(combined.nocomplex.nochrom.3)
cat(nrow(combined.nocomplex.nochrom.3), "explained phenotypes", "\n")
# [1] 4163 explained phenotypes

combined.nocomplex.nochrom.not3 <- subset(combined.nocomplex.nochrom, phenomappingkey!=3)
# nrow(combined.nocomplex.nochrom.not3)
cat(nrow(combined.nocomplex.nochrom.not3), "unexplained phenotypes", "\n")
# [1] 3152 unexplained phenotypes

combined.nocomplex.nochrom.not3.mapping2 <- subset(combined.nocomplex.nochrom, phenomappingkey==2)
# nrow(combined.nocomplex.nochrom.not3.mapping2)
cat(nrow(combined.nocomplex.nochrom.not3.mapping2), "unexplained phenotypes mapped by linkage/etc", "\n")
# [1] 643 unexplained phenotypes mapped by linkage/etc

# length(unique(combined.nocomplex.nochrom.3$geneMIMnum))
cat(length(unique(combined.nocomplex.nochrom.3$geneMIMnum)), "unique genes explain a Mendelian phenotype", "\n")
# [1] 2937 unique genes explaining phenotypes

write.table(subset(combined.nocomplex.nochrom.3), file=paste0(datestamp,".OMIM.phenotypecategories.explained.txt"), quote=FALSE, row.names=FALSE, sep="\t")
write.table(subset(combined.nocomplex.nochrom.not3), file=paste0(datestamp,".OMIM.phenotypecategories.unexplained.txt"), quote=FALSE, row.names=FALSE, sep="\t")



nphenopergene <- table(table(combined.nocomplex.nochrom.3$geneMIMnum))

pdf(paste0(datestamp,".OMIM number of genes vs npheno - barplot.pdf"), width=10, height=5)
x.pos <- barplot(nphenopergene, col=colors[10], xlab="number of Mendelian phenotypes", ylab="number of genes", border=NA)
text(x.pos, nphenopergene, labels=nphenopergene, xpd=TRUE, pos=3, cex=0.8, offset=0.1)
dev.off()

cat(nphenopergene, "number of phenotypes per gene", "\n")
cat(sum(nphenopergene[2:length(nphenopergene)]), "genes with >=2 phenotypes","\n")


