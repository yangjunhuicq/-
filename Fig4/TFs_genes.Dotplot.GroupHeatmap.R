library(Seurat)
library(SCP)

setwd("/home/yang_junhui/NSCLC_ICI/lung_fib/IPF_COPD_Normal/remove_GSE136831_remove_nonFib/disease_subsclutering_fil/Figures/Fig4")
lung_fib.my.clean <- readRDS(file = '/home/yang_junhui/NSCLC_ICI/lung_fib/IPF_COVID19_Normal/merge/COVID19_Donor_PF_COPD/remove/remove_peri_meso/COVID19_Donor_PF_COPD.fib.ana.rds')
lung_fib.my.clean

de_top=read.table("TFs",header=T,row.names=NULL,sep="\t")

lung_fib.my.clean$Normal_celltype=paste(lung_fib.my.clean$Normal,lung_fib.my.clean$celltype,sep='_')
table(lung_fib.my.clean$Normal_celltype,useNA='always')
x=table(lung_fib.my.clean$Normal_celltype,useNA='always')
which(as.numeric(x)>=200)
lung_fib.my.clean=SetIdent(lung_fib.my.clean,value='Normal_celltype')
lung_fib.my.clean_fil=subset(lung_fib.my.clean,idents=names(x)[which(as.numeric(x)>=200)])
table(lung_fib.my.clean_fil$Normal_celltype,useNA='always')

lung_fib.my.clean_fil$Normal_celltype<- factor(lung_fib.my.clean_fil$Normal_celltype,levels=unique(lung_fib.my.clean_fil$Normal_celltype)[c(grep('Adv',unique(lung_fib.my.clean_fil$Normal_celltype)),grep('Alv',unique(lung_fib.my.clean_fil$Normal_celltype)),grep('Myof',unique(lung_fib.my.clean_fil$Normal_celltype)))])
lung_fib.my.clean_fil$Normal_celltype<- factor(lung_fib.my.clean_fil$Normal_celltype,levels=levels(lung_fib.my.clean_fil$Normal_celltype)[c(grep('PF',levels(lung_fib.my.clean_fil$Normal_celltype)),grep('COPD',levels(lung_fib.my.clean_fil$Normal_celltype)),grep('COVID19',levels(lung_fib.my.clean_fil$Normal_celltype)),grep('Donor',levels(lung_fib.my.clean_fil$Normal_celltype)))])
print(unique(lung_fib.my.clean_fil$Normal_celltype))
levels(lung_fib.my.clean_fil)=levels(lung_fib.my.clean_fil$Normal_celltype)
print(levels(lung_fib.my.clean_fil$Normal_celltype))

pdf("TFs_genes.Dotplot.pdf",width=9,height=4)
DotPlot(lung_fib.my.clean_fil, features=de_top$gene, cols = c("grey", "blue"), dot.scale = 6) + RotatedAxis()
DotPlot(lung_fib.my.clean_fil, features=c('NFATC1','RELB','FOSL1'))
DotPlot(lung_fib.my.clean_fil, features=c('NFATC1','RELB','FOSL1',"NFATC2","HMGA1","CARMN",'HLA-DQA1','HLA-DRB5',"HIST2H2BE","CPXM1","SNED1","SEZ6L2",'ARID1B','PIEZO2'), cols = c("grey", "blue"), dot.scale = 6) + RotatedAxis()
DotPlot(lung_fib.my.clean_fil, features=c('NFATC1','RELB','FOSL1',"NFATC2","HMGA1","PLAU",'PLAUR','SERPINE2',"CARMN",'HLA-DQA1','HLA-DRB5',"GJB2","IFI27",'HLA-F',"HIST2H2BE","CPXM1","SNED1","SEZ6L2",'KIAA1324L','MT1F','ARL4C','ARID1B','PIEZO2'), cols = c("grey", "blue"), dot.scale = 6) + RotatedAxis()
dev.off()

library(grid)
pdf("TFs.GroupHeatmap.dot.pdf",width=7,height=7)
ht <- GroupHeatmap(lung_fib.my.clean_fil,
  features = de_top$gene, group.by = "Normal_celltype",
  add_dot = TRUE, add_reticle = TRUE, cluster_rows = TRUE,
  show_row_names = TRUE, ht_params = list(row_names_gp = grid::gpar(fontsize = 6), show_row_dend = FALSE, show_column_dend = FALSE, column_title_gp = grid::gpar(fontsize = 8), row_title_gp = grid::gpar(fontsize = 8), heatmap_legend_param = list(labels_gp =grid::gpar(fontsize=6))), dot_size=unit(4, "mm"), label_size=6
)#heatmap_palette
print(ht$plot)
ht <- GroupHeatmap(lung_fib.my.clean_fil,
  features = c('NFATC1','RELB'), group.by = "Normal_celltype",
  add_dot = TRUE, add_reticle = TRUE, cluster_rows = TRUE,
  show_row_names = TRUE, ht_params = list(row_names_gp = grid::gpar(fontsize = 6), show_row_dend = FALSE, show_column_dend = FALSE, column_title_gp = grid::gpar(fontsize = 8), row_title_gp = grid::gpar(fontsize = 8), heatmap_legend_param = list(labels_gp =grid::gpar(fontsize=6))), dot_size=unit(4, "mm"), label_size=6
)#heatmap_palette
print(ht$plot)
dev.off()


