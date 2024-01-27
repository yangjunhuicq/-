library(msigdbr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

setwd('/home/yang_junhui/NSCLC_ICI/lung_fib/IPF_COPD_Normal/remove_GSE136831_remove_nonFib/disease_subsclutering_fil/Figures/Fig2.addCOVID19/enricher')
inp1=read.table('all.COVID19vsNormal.filt.noRP_noMT.xls.Up',header=T,row.names=NULL,sep="\t")
inp2=read.table('all.PFvsNormal.filt.noRP_noMT.xls.Up',header=T,row.names=NULL,sep="\t")
inp3=read.table('all.COPDvsNormal.filt.noRP_noMT.xls.Up',header=T,row.names=NULL,sep="\t")
data=list(COVIDvsD=unique(inp1[,1]), PFvsD=unique(inp2[,1]), COPDvsD=unique(inp3[,1]))

genesets = msigdbr(species = "Homo sapiens", category = "C2")
print(unique(genesets$gs_subcat))  # ÓÐ¶à¸öÊý¾Ý¿âÀ´Ô´µÄ»ùÒò¼¯¿ÉÑ¡£¬ÕâÀïÑ¡ÓÃKEGG
genesetsuse1 <- subset(genesets, gs_subcat=="CP:KEGG", select = c("gs_name", "gene_symbol"))
print(length(unique(genesetsuse1$gs_name))) #²é¿´ÓÐ¶àÉÙÌõÍ¨Â·£¨186¸ö£©

genesetsuse2 <- subset(genesets, gs_subcat=="CP:REACTOME", select = c("gs_name", "gene_symbol"))
print(length(unique(genesetsuse2$gs_name))) #²é¿´ÓÐ¶àÉÙÌõÍ¨Â·£¨186¸ö£©

genesetsH = msigdbr(species = "Homo sapiens", category = "H")
genesetsuse3 <- subset(genesetsH, select = c("gs_name", "gene_symbol"))
print(length(unique(genesetsuse3$gs_name))) #²é¿´ÓÐ¶àÉÙÌõÍ¨Â·£¨186¸ö£©


genesetsuse=rbind(genesetsuse1,genesetsuse2,genesetsuse3)
print(length(unique(genesetsuse$gs_name)))

res <- compareCluster(geneClusters=data,fun=enricher, TERM2GENE = genesetsuse, pvalueCutoff=1,qvalueCutoff=0.05,minGSSize=1)

pdf("test.pdf")
dotplot(res,showCategory=15,font.size=6,label_format=100)
dotplot(res,showCategory=15,font.size=6,label_format=100,size="count")
dotplot(res,showCategory=30,font.size=6,label_format=100,size="count")
dotplot(res,showCategory=c("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION","REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES","REACTOME_COLLAGEN_DEGRADATION","REACTOME_COLLAGEN_FORMATION","REACTOME_COLLAGEN_CHAIN_TRIMERIZATION","REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES","KEGG_FOCAL_ADHESION","REACTOME_SIGNALING_BY_PDGF","REACTOME_MET_PROMOTES_CELL_MOTILITY","REACTOME_NCAM1_INTERACTIONS","HALLMARK_MYOGENESIS","REACTOME_LAMININ_INTERACTIONS","HALLMARK_COAGULATION","HALLMARK_TNFA_SIGNALING_VIA_NFKB","REACTOME_RHO_GTPASES_ACTIVATE_IQGAPS","REACTOME_COOPERATION_OF_PREFOLDIN_AND_TRIC_CCT_IN_ACTIN_AND_TUBULIN_FOLDING","HALLMARK_INFLAMMATORY_RESPONSE","KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION","REACTOME_INTERFERON_SIGNALING","KEGG_CELL_ADHESION_MOLECULES_CAMS","REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION"),font.size=6,label_format=100,size="count")
dev.off()

write.table(res,file="test.xls",row.names=F,col.names=T,quote=F,sep="\t")

pdf("test.use.pdf",width=4,height=4)
dotplot(res,showCategory=c("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION","REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES","REACTOME_COLLAGEN_DEGRADATION","REACTOME_COLLAGEN_FORMATION","REACTOME_COLLAGEN_CHAIN_TRIMERIZATION","REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES","KEGG_FOCAL_ADHESION","REACTOME_SIGNALING_BY_PDGF","REACTOME_MET_PROMOTES_CELL_MOTILITY","REACTOME_NCAM1_INTERACTIONS","HALLMARK_MYOGENESIS","REACTOME_LAMININ_INTERACTIONS","HALLMARK_COAGULATION","HALLMARK_TNFA_SIGNALING_VIA_NFKB","REACTOME_RHO_GTPASES_ACTIVATE_IQGAPS","REACTOME_COOPERATION_OF_PREFOLDIN_AND_TRIC_CCT_IN_ACTIN_AND_TUBULIN_FOLDING","HALLMARK_INFLAMMATORY_RESPONSE","KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION","REACTOME_INTERFERON_SIGNALING","KEGG_CELL_ADHESION_MOLECULES_CAMS","REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION"),font.size=6,label_format=40,size="count")+theme(axis.text.x=element_text(angle=45,hjust=1))+ scale_size_area(max_size = 4)
dev.off()




saveRDS(res,file=paste0('test',".compareCluster.enricher.COVIDvsD.PFvsD.COPDvsD.rds"))



Cellular_senescence=read.table('/home/yang_junhui/softwares/msigdb/Cellular_senescence',header=T,row.names=NULL,sep="\t")
genesetsuseadd=rbind(genesetsuse,Cellular_senescence)
print(length(unique(genesetsuseadd$gs_name)))
res2 <- compareCluster(geneClusters=data,fun=enricher, TERM2GENE = genesetsuseadd, pvalueCutoff=1,qvalueCutoff=0.05,minGSSize=1)
write.table(res2,file="test.addCellular_senescence.xls",row.names=F,col.names=T,quote=F,sep="\t")

pdf("test.useadd.pdf",width=4,height=4)
dotplot(res2,showCategory=c("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION","REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES","REACTOME_COLLAGEN_DEGRADATION","REACTOME_COLLAGEN_FORMATION","REACTOME_COLLAGEN_CHAIN_TRIMERIZATION","REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES","KEGG_FOCAL_ADHESION","REACTOME_SIGNALING_BY_PDGF","REACTOME_MET_PROMOTES_CELL_MOTILITY","REACTOME_NCAM1_INTERACTIONS","HALLMARK_MYOGENESIS","REACTOME_LAMININ_INTERACTIONS","HALLMARK_COAGULATION","HALLMARK_TNFA_SIGNALING_VIA_NFKB","REACTOME_RHO_GTPASES_ACTIVATE_IQGAPS","REACTOME_COOPERATION_OF_PREFOLDIN_AND_TRIC_CCT_IN_ACTIN_AND_TUBULIN_FOLDING","HALLMARK_INFLAMMATORY_RESPONSE","KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION","REACTOME_INTERFERON_SIGNALING","KEGG_CELL_ADHESION_MOLECULES_CAMS","REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION","KEGG_CELLULAR_SENESCENCE","REACTOME_DISSOLUTION_OF_FIBRIN_CLOT"),font.size=6,label_format=40,size="count")+theme(axis.text.x=element_text(angle=45,hjust=1))+ scale_size_area(max_size = 4)
dev.off()
saveRDS(res2,file=paste0('test',".compareCluster.enricher.COVIDvsD.PFvsD.COPDvsD.add.rds"))

