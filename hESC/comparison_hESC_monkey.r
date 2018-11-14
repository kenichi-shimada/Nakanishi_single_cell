##
setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/hESC_data")
# hesc.exp.lin <- readRDS(file="uniq.exp.10clust.rds")
hesc.exp.lin <- readRDS(file="uniq.exp.7clust_062818.rds")

##
library(openxlsx)
setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/monkey/")
genes <- read.xlsx("nature19096-s2.xlsx",sheet=1,colNames=TRUE,startRow=2)
hs.genes <- genes[[2]]
table(table(hs.genes)) ## all unique
mf.genes <- genes[[4]]
table(table(mf.genes)) ## all unique
mac2hum <- hs.genes
names(mac2hum) <- mf.genes ##17936

setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/monkey")
mf.exp.lin <- readRDS(file="exprs.14lin_all_062818.rds")
rownames(mf.exp.lin) <- mac2hum[rownames(mf.exp.lin)]
has.gene <- !is.na(rownames(mf.exp.lin))
# has.gene
# FALSE  TRUE 
# 11093 17821 

mf.exp.lin <- mf.exp.lin[has.gene,]
overlapped <- intersect(rownames(hesc.exp.lin),rownames(mf.exp.lin)) ## 15170 genes -> 14893

overlapped.mf.genes <- genes$macFas5_gene_symbol[match(overlapped,genes$hg19_gene_symbol)]
setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/monkey")
if(0)saveRDS(overlapped.mf.genes,file="used.mf.genes_062818.rds")
overlapped.mf.genes <- readRDS("used.mf.genes_062818.rds")

## DEG
deg.tmp <- read.xlsx("GeneList_EDFig5b.xlsx",1)
setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/monkey")
overlapped.mf.genes.1 <- intersect(overlapped.mf.genes,deg.tmp$MF5Symbol) ## 5423
if(0)saveRDS(overlapped.mf.genes.1,file="used.mf.genes.1_062818.rds") ## 5423 -> 5410

## exclusive biomarkers
setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/monkey")
bm.mf.genes <- readRDS(file="biomarkers_kw_wilcox_up_062818.rds")
overlapped.1 <- mac2hum[bm.mf.genes]

## compare Spearman correlation
mf.exp.lin <- mf.exp.lin[overlapped.1,]
hesc.exp.lin <- hesc.exp.lin[overlapped.1,]

cor.mat <- apply(hesc.exp.lin,2,function(hs){
	apply(mf.exp.lin,2,function(mf){
		cors <- cor(hs,mf,use="everything",method="spearman")
	})
})

library(gplots)
library(RColorBrewer)
# display.brewer.all()
cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(20)
# cols <- colorRampPalette(c("black","yellow"))(20)

setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/plots")
png("spearman_correlation_mf_hsESC.png",width=500,height=500)
heatmap.2(cor.mat,trace="none",margins=c(5,10),col=cols)
dev.off()

##
setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/plots")
png("spearman_correlation_mf_hsESC_DEG.png",width=500,height=500)
heatmap.2(cor.mat,trace="none",margins=c(5,10),col=cols)
dev.off()

##
setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/plots")
png("spearman_correlation_mf_hsESC_exclusive_up.png",width=500,height=500)
heatmap.2(cor.mat,trace="none",margins=c(5,10),col=cols)
dev.off()

if(0){
	setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/plots")
	png("spearman_correlation_mf_hsESC_exclusive_up-woEXMC.png",width=500,height=500)
	distfun <- function(x)dist(x,method="euclidean")
	# distfun <- function(x)as.dist(1-cor(t(x),method="pearson"))
	hclustfun <- function(x)hclust(x,method="average")
	heatmap.2(cor.mat[!rownames(cor.mat) %in% "EXMC",],trace="none",margins=c(5,10),col=cols,
		hclustfun=hclustfun,distfun=distfun)
	dev.off()
}else{

}

##
setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/plots")
png("spearman_correlation_mf_hsESC_DEG_exclusive_up.png",width=500,height=500)
heatmap.2(cor.mat,trace="none",margins=c(5,10),col=cols)
dev.off()

setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/plots")
png("spearman_correlation_mf_hsESC_DEG_exclusive_up-woEXMC.png",width=500,height=500)
heatmap.2(cor.mat[!rownames(cor.mat) %in% "EXMC",],trace="none",margins=c(5,10),col=cols)
dev.off()

##
png("spearman_correlation_barplot_hypoblast.png",width=300,height=500)
barplot(cor.mat["Hypoblast",],ylab="Spearman correlation against Hypoblast")
dev.off()

##
barplot(cor.mat["VE/YE",],las=2)

if(0){
	sub.cor.mat <- cor.mat[c("PreE-TE","PreL-TE","Post-paTE","Hypoblast","VE/YE","ICM",
		"Pre-EPI","PostE-EPI","PostL-EPI","Gast1","Gast2a","Gast2b"),
		c("3","1","2","7","4","5","6","10","8","9")]
	colnames(sub.cor.mat) <- 1:10
	rowsidecols <- rep(brewer.pal(5,"Set2"),c(3,2,1,3,3))

	setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/plots")
	png("spearman_correlation_noclust_15lins_1.png",width=500,height=500)
	heatmap.2(sub.cor.mat,trace="none",margins=c(5,10),col=cols,
		hclustfun=hclustfun,distfun=distfun,Rowv=F,Colv=F,dendrogram="none",RowSideColors=rowsidecols)
	dev.off()
}else{
	sub.cor.mat <- cor.mat[c("PreE-TE","PreL-TE","Post-paTE","Hypoblast","VE/YE","ICM",
		"Pre-EPI","PostE-EPI","PostL-EPI","Gast1","Gast2"),]
		# c("3","1","2","7","4","5","6","10","8","9")]
	colnames(sub.cor.mat) <- 1:7
	rowsidecols <- rep(brewer.pal(3,"Set2"),c(3,2,6))

	##
	setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/plots")
	png("spearman_correlation_noclust_14lins.png",width=500,height=500)
	heatmap.2(sub.cor.mat,trace="none",margins=c(5,10),col=cols,
		hclustfun=hclustfun,distfun=distfun,Rowv=F,Colv=F,dendrogram="none")
	dev.off()

	##	
	setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/plots")

	png("spearman_correlation_noclust_14lins_1_062818.png",width=500,height=500)
	heatmap.2(sub.cor.mat,trace="none",margins=c(5,10),col=cols,
		hclustfun=hclustfun,distfun=distfun,Rowv=F,Colv=F,dendrogram="none",RowSideColors=rowsidecols)
	dev.off()
	pdf("spearman_correlation_noclust_14lins_1_062818.pdf",width=5,height=5)
	heatmap.2(sub.cor.mat,trace="none",margins=c(5,10),col=cols,
		hclustfun=hclustfun,distfun=distfun,Rowv=F,Colv=F,dendrogram="none",RowSideColors=rowsidecols)
	dev.off()

}