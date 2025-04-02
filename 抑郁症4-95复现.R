setwd("your fiels")
options(stringsAsFactors = F)
rm(list=ls()) #清空变量
job <- "DEG" #设定项目名称
options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
#BiocManager::install("edgeR", force = TRUE)
library(DOSE)
library(limma)
library(edgeR)
library(dplyr)
library(topGO)#绘制通路网络图
library(GOplot)#弦图，弦表图，系统聚类图
library(ggplot2) #用于绘制火山图
#library(ggVolcano)
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(pheatmap)
library(KEGGREST)
library(openxlsx)#读取.xlsx文件
library(circlize)#绘制富集分析圈图
library(ggnewscale)
library(enrichplot)#GO,KEGG,GSEA
library(AnnotationDbi)
library(ComplexHeatmap)#绘制图例
library(clusterProfiler)#GO,KEGG,GSEA

#####数据处理#####
expr_data<-read.csv("matrix_treated.csv",header = T,row.names = NULL,sep = ",") 
#输入文件TPM原始值，行名是基因，列名是样本
#expr_data <- expr_data[, c(2, 11:ncol(expr_data))]
names(expr_data)[1] <- "SYMBOL"
df_summarized <- expr_data %>%
  group_by(SYMBOL) %>%
  summarise_all(mean) # 这里使用的是 mean() 函数，你可以根据需要选择其他汇总函数
expr_data <- as.data.frame(df_summarized)
row.names(expr_data) <- expr_data$SYMBOL
expr_data <- expr_data[, !names(expr_data) == "SYMBOL", drop=FALSE]
#expr_data <- select(expr_data, -SYMBOL,)
#expr_data <- expr_data[, !grepl("SYMBOL", names(expr_data))]

#expr_data <- expr_data[which(rowSums(expr_data)!=0),] #删除表达量为0的基因
#expr_data <- expr_data %>% mutate_all(~10^.)
expr_data = log2(expr_data) #log化处理
expr_data[expr_data == -Inf] = 0 #将log化后的负无穷值替换为0
expr_data[expr_data == NA] = 0   #将log化后的负无穷值替换为0
expr_data[expr_data < 0] <- 0
#expr_data[] <- sapply(expr_data, function(x) if(is.numeric(x)) as.integer(x) else x)

write.csv(expr_data, file = "Tmatrix.csv", row.names = TRUE)


###
expr_data <- read.csv("Tmatrix.csv",
                      header = T,row.names = 1,sep = ',')
group<-read.csv("groupinfo.csv",header = T,row.names = 1,sep = ',') #输入文件，样本信息表，包含分组信息

design <- model.matrix(~0+factor(group$Group))
colnames(design) <- levels(factor(group$Group))
rownames(design) <- colnames(expr_data)

contrast.matrix <- makeContrasts(MDDS-MDD, levels = design) #根据实际的样本分组修改，这里对照组CK，处理组HT

fit <- lmFit(expr_data,design) #非线性最小二乘法
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)#用经验贝叶斯调整t-test中方差的部分
DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
DEG <- na.omit(DEG)

DEG$regulate <- ifelse(DEG$P.Value > 0.05, "unchanged",
                       ifelse(DEG$logFC > 0, "up-regulated",
                              ifelse(DEG$logFC < -0, "down-regulated", "unchanged")))

write.table(table(DEG$regulate),file = paste0(job,"_","DEG_result_1_005.txt"),
            sep = "\t",quote = F,row.names = T,col.names = T)
write.table(data.frame(gene_symbol=rownames(DEG),DEG),file = paste0(job,"_","DEG_result.csv"),
            sep = ",",quote = F,row.names = F,col.names = T)


#------------------------------------------------------------------------------#
DEG_genes <- DEG[DEG$P.Value<0.05&abs(DEG$logFC)>0,]
#DEG_genes <- DEG[DEG$P.Value<0.05,]
DEG_genes <- DEG_genes[order(DEG_genes$logFC),]
DEG_gene_expr <- expr_data[rownames(DEG_genes),]
#DEG_gene_expr[is.infinite(DEG_gene_expr)] = 0
#DEG_gene_expr[DEG_gene_expr == -Inf] = 0
tiff(paste0(job,"_MDDvsCTL_","pheatmap.tiff"),res=300,width=2000,height=2000)
pheatmap(DEG_gene_expr,
         color = colorRampPalette(c("blue","white","red"))(100), #颜色
         scale = "row", #归一化的方式
         border_color = NA, #线的颜色
         fontsize = 8, #文字大小
         cluster_rows = T,
         cluster_cols = F,
         show_rownames = F)
dev.off()


#------------------------------------------------------------------------------#
interesting_genes <- read.csv("E:/Experiment_data_backups/4-95/Bio_a_result/2024-4-18/target_genes.csv"
                              , sep = ',')
interesting_genes <- interesting_genes$Genes
#filtered_df <- DEG[rownames(DEG) %in% interesting_genes, ]
#filtered_df <- filtered_df[order(filtered_df$logFC),]
filtered_df <- expr_data[rownames(expr_data) %in% interesting_genes,]
DEG_gene_expr <- expr_data[rownames(filtered_df),]

write.table(data.frame(DEG_gene_expr),file = "TEMP.csv",
            sep = ",",quote = F,row.names = T,col.names = T)

DEG_gene_expr <- read.csv("TEMP.csv", sep = ',' ,row.names = 'X')

tiff(paste0('Target',"_","pheatmap.tiff"),res=300,width=2000,height=1000)
pheatmap(DEG_gene_expr,
         color = colorRampPalette(c("blue","white","red"))(100), #颜色
         scale = "row", #归一化的方式
         border_color = NA, #线的颜色
         cluster_rows = T,
         cluster_cols = F,
         fontsize = 6,) #文字大小)
dev.off()

# 绘图
DEG$Genes <- rownames(DEG)
tiff(paste0(job,"_MDDvsCTL_","volcano.tiff"),width = 7,height = 7,res = 300,units = 'in')
ggplot(DEG,aes(x=logFC,y=-log10(P.Value)))+ #x轴logFC,y轴adj.p.value
  geom_point(alpha=0.5,size=2,aes(color=regulate))+ #点的透明度，大小
  ylab("-log10(P.Value)")+ #y轴的说明
  scale_color_manual(values = c("blue", "grey", "red"))+ #点的颜色
  geom_vline(xintercept = c(-0,0),lty=4,col ="black",lwd=0.8)+ #logFC分界线
  geom_hline(yintercept=-log10(0.05),lty=4,col = "black",lwd=0.8)+ #adj.p.val分界线
  theme_bw()  #火山图绘制
dev.off()
#------------------------------------------------------------------------------#
#载入差异表达数据，只需基因ID(GO,KEGG,GSEA需要)和Log2FoldChange(GSEA需要)即可
info <- read.csv( "DEG_DEG_result.csv", sep = ',',header = T)
#找到上调与下调基因
info_changed <- info[info$regulate == "up-regulated"|info$regulate =='down-regulated',]
info_up <- info[info$regulate == "up-regulated",]
info_down  <- info[info$regulate =='down-regulated',]

#指定富集分析的物种库
GO_database   <- 'org.Mm.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'mmu' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html
#gene ID转换
gene <- bitr(info$gene_symbol, fromType = 'SYMBOL', 
             toType = 'ENTREZID', OrgDb = GO_database)
gene_up <- bitr(info_up$gene_symbol, fromType = 'SYMBOL',
                toType = 'ENTREZID',OrgDb = GO_database)
gene_down <- bitr(info_down$gene_symbol, fromType = 'SYMBOL',
                  toType = 'ENTREZID',OrgDb = GO_database)

GO_up   <-enrichGO(gene_up$ENTREZID,#GO富集分析
                   OrgDb = GO_database,
                   keyType = "ENTREZID",#设定读取的gene ID类型
                   ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                   pvalueCutoff = 0.05,#设定p值阈值
                   qvalueCutoff = 0.05,#设定q值阈值
                   readable = T)
GO_down <-enrichGO(gene_up$ENTREZID,#GO富集分析
                   OrgDb = GO_database,
                   keyType = "ENTREZID",#设定读取的gene ID类型
                   ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                   pvalueCutoff = 0.05,#设定p值阈值
                   qvalueCutoff = 0.05,#设定q值阈值
                   readable = T)

KEGG_up<-enrichKEGG(gene_up$ENTREZID,#KEGG富集分析
                    organism = KEGG_database,
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
KEGG_down<-enrichKEGG(gene_down$ENTREZID,#KEGG富集分析
                      organism = KEGG_database,
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

write.csv(GO_up@result,'GO_up.csv',row.names = F)
write.csv(GO_down@result,'GO_down.csv',row.names = F)
write.csv(KEGG_up@result,'KEGG_up.csv',row.names = F)
write.csv(KEGG_down@result,'KEGG_down.csv',row.names = F)


#------------------------------------------------------------------------------#
#载入差异表达数据，只需基因ID(GO,KEGG,GSEA需要)和Log2FoldChange(GSEAge(GSEA需要)即可
#info <- read.csv( "DEG_DEG_result.txt", sep = '\t',header = T)
info <- distinct(info, gene_symbol,.keep_all = T)
colnames(info)[1] <- "SYMBOL"
gene <- bitr(info$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
#colnames(gene)[1] <- "SYMBOL"
info_merge <- merge(gene,info,by='SYMBOL')#合并转换后的基因ID和Log2FoldChange
info_merge <- info_merge[order(info_merge$logFC),]

GSEA_input <- info_merge$logFC
names(GSEA_input) = info_merge$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)


GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, 
                     pvalueCutoff = 0.05, scoreType='std',
                     minGSSize = 1, maxGSSize = 5000)#GSEA富集分析

GSEA_GO <- gseGO(GSEA_input, OrgDb = GO_database, # 选择相应的物种数据库，这里以人类为例
                 pvalueCutoff = 0.05, scoreType='std',
                 minGSSize = 1, maxGSSize = 5000)

write.csv(GSEA_KEGG@result,'GSEA_KEGG.csv',row.names = F)
write.csv(GSEA_GO@result,'GSEA_GO.csv',row.names = F)


#------------------------------------------------------------------------------#


