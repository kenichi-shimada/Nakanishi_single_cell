library(gdata)

dir <- "/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/"
setwd(dir);setwd("hES/data")
mat <- read.delim("rpkm.txt",header=T,stringsAsFactors=F)
anno <- read.xls("annotation.xls",1,stringsAsFactors=F)

gns <- mat[[1]]
table(table(gns)) ## unique

mat <- mat[-1] ## 26178 x 1529
rownames(mat) <- gns

## sanity check
table(colnames(mat) %in% anno$Source_Name) ## all TRUE - cell names matched

## 
cell.lins <- tapply(anno$Source_Name,anno$Characteristics_inferred_lineage,identity)

## summarized gene expresssion (average)
exprs.lin <- sapply(cell.lins[-2],function(cells){
	apply(mat[,cells],1,sum)
})

## summarized gene expresssion (normalize and average)
mat.1 <- mat/rep(colSums(mat),each=nrow(mat))
exprs.lin.1 <- sapply(cell.lins[-2],function(cells){
	apply(mat.1[,cells],1,sum)
})

dim(exprs.lin)

setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/hES/")
saveRDS(exprs.lin,file="exprs.3lin_062818.rds")
saveRDS(exprs.lin.1,file="exprs.extranorm.3lin_062818.rds")

##
library(parallel)
biomarkers <- mclapply(cell.lins[-2],function(cells){
	cat("*")
	is.member <- colnames(mat) %in% cells
	p.vals <- apply(mat,1,function(x)
		nlp <- -log10(wilcox.test(x[is.member],x[!is.member],alternative="greater")$p.value)
	)
	is.exprs.by.maj <- apply(mat[,cells],1,function(x)sum(x>0)>length(cells)/2)
	p.vals <- p.vals*as.numeric(is.exprs.by.maj)
},mc.cores=3)


## below is after analyze Nakanishi's data.

##
setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/hES/")
overlapped <- readRDS(file="overlapped_062818.rds")
is.overlapped <- gns %in% overlapped
gns.overlapped <- gns[is.overlapped]

library(RColorBrewer)
cols <- colorRampPalette(brewer.pal(11,"Spectral"))(15)

bm <- lapply(biomarkers,function(x){
	gns.overlapped[head(order(x[is.overlapped],decreasing=T),50)]
})

# bm <- lapply(biomarkers,function(x){
# 	gns[head(order(x,decreasing=T),50)]
# })

bm.genes <- unique(unlist(bm)) #1447 genes => 663 genes (Kruskal-Wallis => wilcox.test)

rownames(exprs.lin) <- rownames(exprs.lin.1) <- gns

##
setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/hES/")
saveRDS(bm.genes,file="biomarkers.hs1_062818.rds")

##
head(bm.genes)
