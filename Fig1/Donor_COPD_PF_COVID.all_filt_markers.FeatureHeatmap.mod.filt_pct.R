library(SCP)
library(Seurat)
setwd("/home/yang_junhui/NSCLC_ICI/lung_fib/IPF_COPD_Normal/remove_GSE136831_remove_nonFib/disease_subsclutering_fil/Figures/Fig1.addCOVID19")
inp=readRDS("/home/yang_junhui/NSCLC_ICI/lung_fib/IPF_COVID19_Normal/merge/COVID19_Donor_PF_COPD/remove/remove_peri_meso/COVID19_Donor_PF_COPD.fib.ana.rds")
inp=SetIdent(inp,value="Normal")
table(inp@active.ident)
genes.ori=read.table("Ref_fib_subtypes.marker_genes.xls.classified.mod",header=T,row.names=1,sep="\t")
genes=genes.ori[intersect(rownames(genes.ori),rownames(inp)),seq(1,ncol(genes.ori)),drop=F]
colnames(genes)="Group"
table(genes$Group,useNA='always')

print("filtering pct 10%")
goi_per5=c()
for(g in rownames(genes))
{
 if(length(which(inp[['RNA']]@counts[g,]>0))/length(colnames(inp)) >= 0.1){goi_per5=c(goi_per5,g)}
}
genes=genes[goi_per5,seq(1,ncol(genes)),drop=F]
table(genes$Group,useNA='always')


ht <- FeatureHeatmap(srt = inp, group.by = "celltype", features = rownames(genes), slot = "data", show_row_names=F, feature_split = genes$Group, height = 3, width = 1.5, label_size = 0, row_names_side='left', cluster_row_slices=T, cluster_column_slices=T, max_cells = 50000, ht_params = list(row_names_gp = grid::gpar(fontsize = 6), show_row_dend = FALSE, show_column_dend = FALSE, column_title_gp = grid::gpar(fontsize = 8), row_title_gp = grid::gpar(fontsize = 8), heatmap_legend_param = list(labels_gp =grid::gpar(fontsize=6))))
pdf("Donor_COPD_PF.all_filt_markers.FeatureHeatmap.mod.filt_pct.pdf",width=7,height=4)
print(ht$plot)
dev.off()
htr <- FeatureHeatmap(srt = inp, group.by = "celltype", features = rownames(genes), slot = "data", show_row_names=F, feature_split = genes$Group, height = 3, width = 1.5, label_size = 0, row_names_side='left', cluster_row_slices=T, cluster_column_slices=T, max_cells = 1000, ht_params = list(row_names_gp = grid::gpar(fontsize = 6), show_row_dend = FALSE, show_column_dend = FALSE, column_title_gp = grid::gpar(fontsize = 8), row_title_gp = grid::gpar(fontsize = 8), heatmap_legend_param = list(labels_gp =grid::gpar(fontsize=6))))
pdf("Donor_COPD_PF.all_filt_markers.FeatureHeatmap.random1000.mod.filt_pct.pdf",width=7,height=4)
print(htr$plot)
dev.off()
htrr <- FeatureHeatmap(srt = inp, cell_annotation="celltype", features = rownames(genes), slot = "data", show_row_names=F, feature_split = genes$Group, height = 3, width = 1.5, label_size = 0, row_names_side='left', cluster_row_slices=T, cluster_columns=T, max_cells = 1000, ht_params = list(row_names_gp = grid::gpar(fontsize = 6), show_row_dend = FALSE, show_column_dend = FALSE, column_title_gp = grid::gpar(fontsize = 8), row_title_gp = grid::gpar(fontsize = 8), heatmap_legend_param = list(labels_gp =grid::gpar(fontsize=6))))
pdf("Donor_COPD_PF.all_filt_markers.FeatureHeatmap.random1000.col_nosplit.mod.filt_pct.pdf",width=7,height=4)
print(htrr$plot)
dev.off()



