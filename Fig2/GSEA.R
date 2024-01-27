library(msigdbr)
library(clusterProfiler)
library(enrichplot)
out=read.table(file="infile",header=T,row.names=1,sep="\t")
genelist <- structure(out$rank, names = rownames(out))
#KEGG
genesets = msigdbr(species = "Homo sapiens", category = "C2")
genesetsuse <- subset(genesets, gs_subcat=="CP:KEGG", select = c("gs_name", "gene_symbol"))

res <- GSEA(genelist, TERM2GENE = genesetsuse, eps = 0, pvalueCutoff=1)
res <- data.frame(res)
write.table(res, file="infile.GSEA_KEGG.xls", row.names = F, sep="\t", quote=F)

