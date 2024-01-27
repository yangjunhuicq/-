library(ggplot2)
cellfrac=read.table("cellfrac-comp_celltype_by_Normal_by_orig.ident_metadata_cast_pt_melt.txt",header=T,row.names=NULL,sep="\t")
cellfrac$variable=gsub('pt_','',cellfrac$variable)

coluse=c("COVID19"=rgb(51,160,44,maxColorValue=255),"PF"=rgb(178,223,138,maxColorValue=255),"COPD"=rgb(166,206,227,maxColorValue=255),"Donor"=rgb(31,120,180,maxColorValue=255))
coluse=c("#3B4992FF","#EE0000FF","#008B45FF","#631879FF")
coluse=c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00")
#c("Gold3", "CornflowerBlue", "SeaGreen", "MediumPurple")
    pdf(paste0('cellfrac-comp_celltype_by_Normal_by_orig.ident',".size.pdf"),width=3.6,height=3)
    print(ggplot(cellfrac, aes(x=variable, y=value, fill=Normal, color=Normal)) + geom_boxplot(size=0.2,outlier.size=0.2)+xlab('Fibroblast subtypes')+ylab('Cell fraction')+theme_classic()+theme(element_blank(), axis.line.y.left=element_line(color="black"), axis.line.x.bottom=element_line(color="black"), axis.title=element_text(colour = "black", size=7, face='bold'), axis.text = element_text(colour = "black", size=6), axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_text(colour = "black", size=7), legend.text=element_text(colour = "black", size=6))+ylim(c(-0.001,1)) + scale_x_discrete(labels=c("Adv","Alv","FibM","Lipo","MyoF")))#scale_color_manual(values=coluse)+ + scale_fill_manual(values=coluse)
    print(ggplot(cellfrac, aes(x=variable, y=value, fill=Normal, color=Normal)) + geom_boxplot(size=0.2,outlier.size=0.2)+xlab('Fibroblast subtypes')+ylab('Cell fraction')+theme_classic()+theme(element_blank(), axis.line.y.left=element_line(color="black"), axis.line.x.bottom=element_line(color="black"), axis.title=element_text(colour = "black", size=7, face='bold'), axis.text = element_text(colour = "black", size=6), axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_text(colour = "black", size=7), legend.text=element_text(colour = "black", size=6))+ylim(c(-0.001,1)) + scale_color_manual(values=coluse) + scale_fill_manual(values=coluse)+scale_x_discrete(labels=c("Adv","Alv","FibM","Lipo","MyoF")))#
    dev.off()

