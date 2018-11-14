library(parallel)

setwd(esc.dir <- "/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/hESC_data/")
## annotation
setwd("Annotation")
esc <- read.csv("hESC_grouping_cluster.csv") 
table(esc$cluster) ## 10 classes
#    1    2    3    4    5    6    7    8    9   10 
# 1023 1009  904  773  757  708  602  516  419  353 
if(0){
	bulk.anno <- read.table("sex_annotation_Bulk_H1_H9.txt",sep="\t",header=T)
	ncad.anno <- read.table("sex_annotation_NCADpos_H1_H9.txt",sep="\t",header=T)
	anno <- data.frame(rbind(bulk.anno,ncad.anno),library=rep(c("Bulk","NCADpos"),c(nrow(bulk.anno),nrow(ncad.anno))))
	dim(anno) ## 6914 cells, not 7064

	table(anno$library,anno$Sex)
	table(table(anno$Cell.Barcode)) ## 10 barcodes are duplicated
	dup.barcodes <- names(which(table(anno$Cell.Barcode)==2))
	dup.anno <- data.frame(anno[anno$Cell.Barcode %in% dup.barcodes,])
	dup.anno$Cell.Barcode <- factor(dup.anno$Cell.Barcode)
	table(dup.anno$Cell.Barcode,dup.anno$library) ## at least they are not overlapped between the libraries

	bcds <- paste0(anno$Cell.Barcode,"-",as.numeric(anno$library))
	overlapped <- bcds[bcds %in% rownames(esc)]

	esc <- esc[rownames(esc) %in% overlapped,]
	anno <- anno[bcds %in% overlapped,]
	anno$cluster <- esc$cluster[match(bcds[bcds %in% overlapped],rownames(esc))]
	table(anno$Sex,anno$cluster)
	table(anno$library,anno$cluster)
	table(anno$library,anno$cluster,anno$Sex)
}
##


## 
setwd(esc.dir)
norm.tmp <- readLines("normalized_log_counts.csv")
list.tmp <- mclapply(norm.tmp,function(x){
	gsub("\"","",strsplit(x,",")[[1]])
},mc.cores=9)
genes <- sapply(list.tmp[-1],function(x)x[1])

table(table(genes))
#     1     2 
# 25554    15 
n.genes <- table(genes)
dup.genes <- names(which(n.genes==2))
uniq.genes <- unique(genes) ## 25569 <- 25584

norm.exp <- do.call(rbind,lapply(list.tmp[-1],function(x)as.numeric(x[-1])))
cells <- list.tmp[[1]]
rm(list.tmp)
norm.uniq.exp <- do.call(rbind,mclapply(uniq.genes,function(gn){
	if(n.genes[gn]==1){
		idx <- which(genes==gn)

		return(norm.exp[idx,])
	}else if(n.genes[gn]==2){
		idx <- which(genes==gn)
		return(apply(norm.exp[idx,],2,sum))
	}
},mc.cores=9))
##
colnames(norm.uniq.exp) <- cells
rownames(norm.uniq.exp) <- uniq.genes
##
rm(norm.exp)

##
if(0)saveRDS(norm.uniq.exp,file="uniq.exp.rds")
norm.uniq.exp <- readRDS("uniq.exp.rds")
bcds <- colnames(norm.uniq.exp)
table(bcds %in% rownames(esc)) ## they are all unique and matched.
cells.clust <- tapply(rownames(esc),esc$cluster,identity)
norm.exp.clust <- sapply(cells.clust,function(this.bcs){
	apply(norm.uniq.exp[,this.bcs],1,sum) ## intentionally used sum instead of mean (concerns about small number being eliminated)
})

if(0){
	setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/hESC_data")
	saveRDS(norm.exp.clust,file="uniq.exp.10clust.rds")
}

