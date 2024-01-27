library(ggplot2)
library(cowplot)
library(ggpubr)
plots_grid <- function(plots=plots,filename=filename,width=width,height=height)
{
num=length(plots)
x=floor(num/6)-1
y=num-(x+1)*6
print(num)
print(x)
print(y)

pdf(filename, onefile = TRUE, width=width,height=height)
i=-1
if(x>=0)
{
 for(i in c(0:x))
 {
        print(plot_grid(plotlist=plots[c((i*6+1):(i*6+6))],ncol=3,nrow=2,label_size=4))
 }
}
re_list=list()
if(y>0)
{
 for(j in c(1:y))
 {
  re_list[[length(re_list)+1]]=plots[[i*6+6+j]]
 }
 print(plot_grid(plotlist=re_list,ncol=3,nrow=2,label_size=4))
}
if(num==6)
{
 i=0
 print(plot_grid(plotlist=plots[c((i*6+1):(i*6+6))],ncol=3,nrow=2,label_size=4))
}
dev.off()
}

plots=list()
goi=rownames(read.table("selected_genes",header=T,row.names=1))

frac_table=read.table("selected_genes.Patient.nosubtype.cellnum.add_orig.ident.xls",header=T,row.names=NULL,sep="\t")

    plots <- list()
    for (i in 1:length(goi)) {
        inp_use2=frac_table[frac_table$gene == goi[i], ]
        plots[[length(plots)+1]]=ggplot(inp_use2, aes(x=disease, y=expressed_ratio, fill=disease)) + geom_boxplot(size=0.2,outlier.size=0.2)+geom_jitter(size=0.5)+theme_classic()+theme(element_blank(), axis.line.y.left=element_line(color="black"), axis.line.x.bottom=element_line(color="black"), axis.title=element_text(colour = "black", size=7, face='bold'), axis.text = element_text(colour = "black", size=6), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), legend.title = element_text(colour = "black", size=7), legend.text=element_text(colour = "black", size=6), legend.pos='none')+xlab('')+ylab(goi[i])+stat_compare_means(comparisons=list(c('COPD','PF'),c('COPD','COVID19'),c('PF','COVID19'),c('COPD','Donor'),c('PF','Donor'),c('COVID19','Donor')),method = "wilcox.test",size=2)+stat_compare_means()+ylim(c(-0.001,2))
    }
plots_grid(plots=plots,paste0('selected_genes.Patient.nosubtype.',".boxplot.pdf"),width=7,height=5)


    plots <- list()
    for (i in 1:length(goi)) {
        inp_use2=frac_table[which(frac_table$gene == goi[i] & frac_table$cellnum>=50), ]
        plots[[length(plots)+1]]=ggplot(inp_use2, aes(x=disease, y=expressed_ratio, fill=disease)) + geom_boxplot(size=0.2,outlier.size=0.2)+geom_jitter(size=0.5)+theme_classic()+theme(element_blank(), axis.line.y.left=element_line(color="black"), axis.line.x.bottom=element_line(color="black"), axis.title=element_text(colour = "black", size=7, face='bold'), axis.text = element_text(colour = "black", size=6), axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), legend.title = element_text(colour = "black", size=7), legend.text=element_text(colour = "black", size=6), legend.pos='none')+xlab('')+ylab(goi[i])+stat_compare_means(comparisons=list(c('COPD','PF'),c('COPD','COVID19'),c('PF','COVID19'),c('COPD','Donor'),c('PF','Donor'),c('COVID19','Donor')),method = "wilcox.test",size=2)+stat_compare_means()+ylim(c(-0.001,2))
    }
plots_grid(plots=plots,paste0('selected_genes.Patient.nosubtype.',".boxplot.fil_cell50.pdf"),width=7,height=5)

