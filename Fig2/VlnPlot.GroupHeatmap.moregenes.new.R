library(Seurat)
library(SCP)

setwd("/home/yang_junhui/NSCLC_ICI/lung_fib/IPF_COPD_Normal/remove_GSE136831_remove_nonFib/disease_subsclutering_fil/Figures/Fig2.addCOVID19/selected_genes")
lung_fib.my.clean <- readRDS(file = '/home/yang_junhui/NSCLC_ICI/lung_fib/IPF_COVID19_Normal/merge/COVID19_Donor_PF_COPD/remove/remove_peri_meso/COVID19_Donor_PF_COPD.fib.ana.rds')
lung_fib.my.clean

de_top=read.table("selected_genes.addgroup.new",header=T,row.names=NULL,sep="\t")

lung_fib.my.clean$Normal_celltype=paste(lung_fib.my.clean$Normal,lung_fib.my.clean$celltype,sep='_')
table(lung_fib.my.clean$Normal_celltype,useNA='always')
x=table(lung_fib.my.clean$Normal_celltype,useNA='always')
which(as.numeric(x)>=200)
lung_fib.my.clean=SetIdent(lung_fib.my.clean,value='Normal_celltype')
lung_fib.my.clean_fil=subset(lung_fib.my.clean,idents=names(x)[which(as.numeric(x)>=200)])
table(lung_fib.my.clean_fil$Normal_celltype,useNA='always')

library(grid)
pdf("selected.VlnPlot.GroupHeatmap.moregenes.dot.new.pdf",width=5,height=4)
ht <- GroupHeatmap(lung_fib.my.clean_fil,
  features = de_top$gene, feature_split = de_top$group, group.by = "Normal",
  add_dot = TRUE, add_reticle = TRUE, cluster_rows = TRUE, 
  show_row_names = TRUE, ht_params = list(row_names_gp = grid::gpar(fontsize = 6), show_row_dend = FALSE, show_column_dend = FALSE, column_title_gp = grid::gpar(fontsize = 8), row_title_gp = grid::gpar(fontsize = 8), heatmap_legend_param = list(labels_gp =grid::gpar(fontsize=6))), dot_size=unit(4, "mm"), label_size=6
)#heatmap_palette 
print(ht$plot)
dev.off()


