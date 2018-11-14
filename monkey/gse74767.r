library(SC3)
library(scran)
library(scater)

setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/monkey/")
all.counts <- read.delim(file="GSE74767_SC3seq_Cy_ProcessedData.txt", header=TRUE,
	stringsAsFactors=F)
entrez.id <- all.counts[,1]
gene.symbol <- all.counts[,2]
all.counts <- all.counts[,-1]
all.counts <- all.counts[,-1]
all.counts <- as.matrix(all.counts)
sce <- SingleCellExperiment(assays = list(counts = all.counts, 
	logcounts = log2(all.counts + 1)))
rowData(sce)$entrez_id <- entrez.id
rowData(sce)$feature_symbol <- gene.symbol
rownames(sce) <- gene.symbol
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
norm_exprs(sce) <- logcounts(sce)

##retrieving annotation
library(openxlsx)
setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/monkey/")
annotation <- read.xlsx("nature19096-s1.xlsx",sheet=3,colNames=TRUE,startRow=2)
rownames(annotation) <- annotation$SampleID
annotation <- annotation[-1] ## 581 samples
cellID <- colnames(sce) ## 421 samples

table(cellID %in% rownames(annotation)) ## 421 samples were all in annotation
annotation <- annotation[rownames(annotation) %in% cellID,]

lineage <- annotation$Group[match(cellID,rownames(annotation))]
lineage[lineage %in% c("Gast2a","Gast2b")] <- "Gast2"
embryonic_day <- annotation$Embryonic_Day[match(cellID,rownames(annotation))]

##
table(lineage) ## 15 lineages
cells.lineage <- tapply(cellID,lineage,identity) 

setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/monkey")
saveRDS(names(cells.lineage),file="cell-lineage.rds")

norm.exprs <- norm_exprs(sce)
has.missing <- which(rowSums(is.na(norm.exprs))>0)
norm.exprs <- norm.exprs[-has.missing,]

exprs.lin <- sapply(cells.lineage,function(cells){
	apply(norm.exprs[,cells],1,sum)
})

##
setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/monkey")
if(0)saveRDS(exprs.lin,file="exprs.14lin_all_062818.rds")

## used genes in hESC
setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/monkey")
overlapped.mf.genes <- readRDS(file="used.mf.genes.1_062818.rds")
norm.exprs <- norm.exprs[overlapped.mf.genes,] ## 5423 x 421


## DEG
biomarkers <- mclapply(cells.lineage,function(cells){
	cat("*")
	is.member <- colnames(norm.exprs) %in% cells
	p.vals <- apply(norm.exprs,1,function(x)
		nlp <- -log10(wilcox.test(x[is.member],x[!is.member],alternative="greater")$p.value)
	)
},mc.cores=8)

bm <- do.call(cbind,biomarkers)

##
library(RColorBrewer)
cols <- colorRampPalette(brewer.pal(11,"Spectral"))(15)
plot.ecdf(bm,col=NA)
for(i in seq(biomarkers))plot.ecdf(biomarkers[[i]],col=cols[i],add=T)

bm <- lapply(biomarkers,function(x){
	names(head(sort(x,decreasing=T),50))
})
boxplot(tapply(norm.exprs["LAMP2",],lineage,identity),las=2)

bm.genes <- unique(unlist(bm)) #1447 genes => 663 genes (Kruskal-Wallis => wilcox.test)

setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/monkey")
saveRDS(bm.genes,file="biomarkers_kw_wilcox_up_062818.rds")

boxplot(biomarkers)



##heatmap generation using only PrE cells and provided gene list.
library(pheatmap)
combined3.markers <- scan("combined3_Cy.txt", what="", sep="\t")

set.seed(1)
combined3.markers <- sample(28914,100)
Hypo <-filter(sce, sce$lineage == "Hypoblast")
hypo_combinedgenes3.vals <- norm_exprs(Hypo)[combined3.markers,,drop=FALSE]

png("heatmap_monkeyHypo_plu-selectedPrEgenes3.png", width = 1200, height = 800)
pheatmap(hypo_combinedgenes3.vals, cluster_rows = FALSE, labels_col=NA, cellhight = 15, cellwidth = 2, treeheight_col = 50, fontsize = 15)
dev.off()

##heatmap generation using all cells, provided gene list and lineage annotation (but not sorted)

##
## 1
system.time(d <- dist(t(norm_exprs(sce)),method="euclidean")) ## create dist object
system.time(hc <- hclust(d,method="ward.D2")) ## create hclust object
class(hc)

png("heatmap_monkeyHypo_plu-selectedPrEgenes3.png", width = 1200, height = 800)
pheatmap(all_combinedgenes3.vals, cluster_rows = FALSE, cluster_cols = hc, 
	labels_col=rep("",ncol(sce)), cellhight = 15, cellwidth = 2, treeheight_col = 50, fontsize = 15, 
	annotation_col=data.frame(Cluster=sce$lineage, row.names=colnames(sce)))
dev.off()

## 2
uniq.labs <-  c("cyESCFF","cyESCoF","EXMC","Gast1","Gast2a","Gast2b","Hypoblast","ICM","Post-paTE","PostE-EPI","PostL-EPI","Pre-EPI","PreE-TE","PreL-TE","VE/YE") ## change the order
o <- order(as.numeric(factor(sce$lineage,levels=uniq.labs)),decreasing=F)
all_combinedgenes3.vals <- norm_exprs(sce)[combined3.markers,o,drop=FALSE]

png("heatmap_monkeyHypo_plu-selectedPrEgenes3_no-clust.png", width = 1200, height = 800)
pheatmap(all_combinedgenes3.vals, cluster_rows = FALSE, cluster_cols = FALSE, 
	labels_col=rep("",ncol(sce)), cellhight = 15, cellwidth = 2, treeheight_col = 50, fontsize = 15, 
	annotation_col=data.frame(Cluster=sce$lineage, row.names=colnames(sce)))
dev.off()

