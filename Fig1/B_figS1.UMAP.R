library(Seurat)
library(ggplot2)
lung_fib.my.clean <- readRDS(file = "/home/yang_junhui/NSCLC_ICI/lung_fib/IPF_COVID19_Normal/merge/COVID19_Donor_PF_COPD/remove/remove_peri_meso/COVID19_Donor_PF_COPD.fib.ana.rds")
lung_fib.my.clean

table(lung_fib.my.clean$Normal)

print("");print("resolution selected:")
lung_fib.my.clean <- SetIdent(lung_fib.my.clean, value="celltype")
table(lung_fib.my.clean@active.ident)


pdf("All_cells.UMAP_plot.pdf",width=3.5,height=2.5)
DimPlot(lung_fib.my.clean, reduction = "umapbyharmony15", label=F)+xlab("UMAP-1")+ylab("UMAP-2")+labs(title="",color = "Fibroblast subtype")+theme(plot.title = element_text(size = 10), axis.title=element_text(colour = "black", size=7, face='bold'), legend.title = element_text(colour = "black", size=7), legend.text=element_text(colour = "black", size=7), axis.text = element_text(size=7))+guides(color=guide_legend(override.aes = list(size=2)))
DimPlot(lung_fib.my.clean, reduction = "umapbyharmony15", label=F, pt.size=0.1,shuffle=T)+xlab("UMAP-1")+ylab("UMAP-2")+labs(title="",color = "Fibroblast subtype")+theme(plot.title = element_text(size = 10), axis.title=element_text(colour = "black", size=7, face='bold'), legend.title = element_text(colour = "black", size=7), legend.text=element_text(colour = "black", size=7), axis.text = element_text(size=7))+guides(color=guide_legend(override.aes = list(size=2)))
DimPlot(lung_fib.my.clean, reduction = "umapbyharmony15", label=T, group.by="Normal")+xlab("UMAP-1")+ylab("UMAP-2")+labs(title="",color = "Disease type")+theme(plot.title = element_text(size = 10), axis.title=element_text(colour = "black", size=7, face='bold'), legend.title = element_text(colour = "black", size=7), legend.text=element_text(colour = "black", size=7), axis.text = element_text(size=7))+guides(color=guide_legend(override.aes = list(size=2)))
dev.off()



print("");print("resolution selected:")
lung_fib.my.clean <- SetIdent(lung_fib.my.clean, value="gse_alias")
table(lung_fib.my.clean@active.ident)


pdf("All_cells.UMAP_plot.gse_alias.pdf",width=3.5,height=2.5)
DimPlot(lung_fib.my.clean, reduction = "umapbyharmony15", label=F, cols=c("GSE132771"="#696969","GSE150247"="#dcdcdc","ERP114453"="#add8e6","GSE128033"="#9932cc","GSE128169"="#008080","GSE135893"="#0000cd","GSE168191"="#ffa07a","GSE173896"="#7b68ee","GSE196341"="#ff69b4","GSE171668"="#8fbc8f", "GSE158127"="#9acd32"), pt.size=0.1,shuffle=T)+xlab("UMAP-1")+ylab("UMAP-2")+labs(title="",color = "Dataset")+theme(plot.title = element_text(size = 10), axis.title=element_text(colour = "black", size=7, face='bold'), legend.title = element_text(colour = "black", size=7), legend.text=element_text(colour = "black", size=7), axis.text = element_text(size=7))+guides(color=guide_legend(override.aes = list(size=2)))
dev.off()

