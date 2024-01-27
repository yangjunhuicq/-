#/home/yang_junhui/softwares/anaconda/Anaconda3-5.3.1/envs/SCP/bin/R
library(SCP)
library(Seurat)

setwd("/home/yang_junhui/NSCLC_ICI/lung_fib/IPF_COPD_Normal/remove_GSE136831_remove_nonFib/disease_subsclutering_fil/Figures/Fig2.addCOVID19")
inp=readRDS("/home/yang_junhui/NSCLC_ICI/lung_fib/IPF_COVID19_Normal/merge/COVID19_Donor_PF_COPD/remove/remove_peri_meso/COVID19_Donor_PF_COPD.fib.ana.rds")
inp=SetIdent(inp,value="celltype")
table(inp@active.ident)
inp <- subset(inp,idents=c("Adventitial Fibroblast","Alveolar Fibroblast","Fibromyocyte","Lipofibroblast","Myofibroblast"))
table(inp@active.ident)

inp=SetIdent(inp,value="Normal")
table(inp@active.ident)

genes=read.table("COVIDvsD.PFvsD.COPDvsD.merge.new.xls.Up",header=T,row.names=1,sep="\t")

ht <- FeatureHeatmap(srt = inp, group.by = "Normal", features = rownames(genes), slot = "data", species = "Homo_sapiens", db = c("GO_BP", "GO_CC", "GO_MF" , "KEGG", "Reactome"), anno_terms = TRUE, topTerm=15, show_row_names=F, feature_split = NULL, height = 6, width = 3, terms_fontsize = 6, terms_width = grid::unit(5, "cm"), label_size = 0, row_names_side='left', cluster_rows=T, cluster_column_slices=T, max_cells = 1000, n_split=3, split_method = "hclust", ht_params = list(row_names_gp = grid::gpar(fontsize = 6), show_row_dend = T, show_column_dend = FALSE, column_title_gp = grid::gpar(fontsize = 8), row_title_gp = grid::gpar(fontsize = 8), heatmap_legend_param = list(labels_gp =grid::gpar(fontsize=6))))
pdf("COPD_PF_D.DEG_fil.FeatureHeatmap.GO_each.KEGG.Reactome.data.Up.new.pdf",width=20,height=7)
print(ht$plot)
dev.off()

for(i in seq(3,7))
{
ht2 <- FeatureHeatmap(srt = inp, group.by = "Normal", features = rownames(genes), slot = "data", species = "Homo_sapiens", db = c("GO", "KEGG", "Reactome"), anno_terms = TRUE, topTerm=15, show_row_names=F, feature_split = NULL, height = 6, width = 3, terms_fontsize = 6, terms_width = grid::unit(5, "cm"), label_size = 0, row_names_side='left', cluster_rows=T, cluster_column_slices=T, max_cells = 1000, n_split=i, split_method = "hclust", ht_params = list(row_names_gp = grid::gpar(fontsize = 6), show_row_dend = FALSE, show_column_dend = FALSE, column_title_gp = grid::gpar(fontsize = 8), row_title_gp = grid::gpar(fontsize = 8), heatmap_legend_param = list(labels_gp =grid::gpar(fontsize=6))))#feature_annotation = c("TF", "SP"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
pdf(paste0("COPD_PF_D.DEG_fil.FeatureHeatmap.GO.KEGG.Reactome.data.row_",i,".Up.new.pdf"),width=13,height=7)
print(ht2$plot)
dev.off()
}

#clustering_distance_rows="pearson" not work
ht22 <- FeatureHeatmap(srt = inp, group.by = "Normal", features = rownames(genes), slot = "data", species = "Homo_sapiens", db = c("GO", "KEGG", "Reactome"), anno_terms = TRUE, topTerm=15, show_row_names=F, feature_split = NULL, height = 6, width = 3, terms_fontsize = 6, terms_width = grid::unit(5, "cm"), label_size = 0, row_names_side='left', cluster_rows=T, cluster_column_slices=T, max_cells = 1000, n_split=3, split_method = "hclust", ht_params = list(row_names_gp = grid::gpar(fontsize = 6), show_row_dend = FALSE, show_column_dend = FALSE, column_title_gp = grid::gpar(fontsize = 8), row_title_gp = grid::gpar(fontsize = 8), heatmap_legend_param = list(labels_gp =grid::gpar(fontsize=6)), clustering_distance_rows="pearson"))#feature_annotation = c("TF", "SP"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
pdf("COPD_PF_D.DEG_fil.FeatureHeatmap.GO.KEGG.Reactome.data.row_pearson.Up.new.pdf",width=13,height=7)
print(ht22$plot)
dev.off()

