library(msigdbr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

setwd('/home/yang_junhui/NSCLC_ICI/lung_fib/IPF_COPD_Normal/remove_GSE136831_remove_nonFib/disease_subsclutering_fil/Figures/Fig4')
fsub='Adv'
for(compa in c('HMGA1_target_genes.filt_per15.b2_qval0.001'))
 {
  for(datab in c('GO','KEGG','REACTOME','BIOCARTA','Hallmark'))
  {
   assign(paste0(fsub,'_',compa,'_',datab),readRDS(paste0('enricher/',fsub,'.',compa,'.enricher_',datab,'.rds')))
   print(paste0('enricher/',fsub,'.',compa,'.enricher_',datab,'.rds'))
  }
 }
ls()

res=merge_result(list(Adv_HMGA1_target_genes.filt_per15.b2_qval0.001_GO=Adv_HMGA1_target_genes.filt_per15.b2_qval0.001_GO,Adv_HMGA1_target_genes.filt_per15.b2_qval0.001_KEGG=Adv_HMGA1_target_genes.filt_per15.b2_qval0.001_KEGG,Adv_HMGA1_target_genes.filt_per15.b2_qval0.001_Hallmark=Adv_HMGA1_target_genes.filt_per15.b2_qval0.001_Hallmark,Adv_HMGA1_target_genes.filt_per15.b2_qval0.001_BIOCARTA=Adv_HMGA1_target_genes.filt_per15.b2_qval0.001_BIOCARTA,Adv_HMGA1_target_genes.filt_per15.b2_qval0.001_REACTOME=Adv_HMGA1_target_genes.filt_per15.b2_qval0.001_REACTOME))

selected=read.table('HMGA1.selected_terms',header=T,row.names=1)

dot_df=res@compareClusterResult[which(res@compareClusterResult$ID %in% rownames(selected)),]

write.table(dot_df,file='HMGA1.selected_terms.info.xls',row.names=F,col.names=T,quote=F,sep="\t")


library(forcats)
dot_df2=dot_df
dot_df2$Cluster=gsub('_GO|_KEGG|_REACTOME|_BIOCARTA|_Hallmark','',dot_df2$Cluster,perl=T)
pdf("HMGA1.selected_terms.merge_DB.selected.pdf",width=5,height=2.5)
print(ggplot(dot_df2, aes(x = GeneRatio, y = fct_reorder(Description, qvalue))) +
               geom_point(aes(size = Count, color = -log10(qvalue))) +
               theme_bw(base_size = 6) +
        #scale_colour_gradient(low="white",high='red') +
        ylab(NULL)+theme(axis.text.x = element_text(angle = 45, hjust=1)))
dev.off()

pdf("HMGA1.selected_terms.merge_DB.selected.test.pdf",width=5,height=2.5)
print(ggplot(dot_df2, aes(x = GeneRatio, y = fct_reorder(Description, qvalue))) +
               geom_point(aes(size = Count, color = -log10(qvalue))) +
               theme_bw(base_size = 6) + 
               scale_color_gradient(name = "q-value", limits = c(0, 5), low = "lightblue", high = "darkred") +
        ylab(NULL)+theme(axis.text.x = element_text(angle = 45, hjust=1)))
dev.off()

pdf("HMGA1.selected_terms.merge_DB.selected.test2.pdf",width=5.5,height=2.5)
print(ggplot(dot_df2, aes(x = Count, y = fct_reorder(Description, qvalue))) +
               geom_point(aes(size = -log10(qvalue), color = -log10(qvalue))) +
               theme_bw(base_size = 6) +
               scale_color_gradient(name = "q-value", limits = c(0, 5), low = "lightblue", high = "darkred") +
        ylab(NULL)+theme(axis.text.x = element_text(angle = 45, hjust=1)))
dev.off()

