library(Rtsne)
library(RColorBrewer)

dirs <- c("non-normalization","normalized","not_normalized" )

##
setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/normalized")

##
bcs <- readLines("barcodes.tsv")
conditions <- sub(".+-([12])$","\\1",bcs)
table(conditions)

#    1    2 
# 2836 2207 

##
gns <- read.table("genes.tsv",stringsAsFactors=F)
names(gns) <- c("ensembl_id","symbol")
table(table(gns$ensembl)) ## unique 33354
table(table(gns$symbol)) ## non-unique, 1-4 entries per symbol

mt.inds <- grep("hg19_MT-",gns$symbol)

##

norm <- readLines("matrix.mtx")
n <- as.numeric(unlist(strsplit(norm[3]," ")))
mat <- array(0,n[1:2])
stats <- t(sapply(strsplit(norm[-(1:3)]," "),as.numeric))

mat[stats[,1:2]] <- stats[,3]


##
n.reads <- colSums(mat)
n.reads.cds <- tapply(n.reads,conditions,identity)
plot.ecdf(n.reads.cds[[1]],col=1)
plot.ecdf(n.reads.cds[[2]],col=2,add=T)

hist(n.reads,breaks=100)
thres <- 30000
n.reads[n.reads > thres] <- thres
n.levels <- round(n.reads/thres*11)

n.cols <- brewer.pal(11,"Spectral")

##
nmat <- apply(log2(mat+1),2,function(x)x/sum(x))
rownames(nmat) <- gns$ensembl
colnames(nmat) <- bcs


##
system.time(pca_result <- prcomp(t(nmat))) ## 
#    user  system elapsed 
# 979.474  19.533 507.419 

initial_dims <- 10
plot(pca_result$sdev,xlim=c(0,100))
abline(v=10,col=2)
X <- pca_result$x[, 1:initial_dims]
system.time(rt <- Rtsne(X,pca=FALSE,theta=0.5,perplexity=30,max_iter=2000))
#    user  system elapsed 
# 133.964   0.178 134.116 

##
plot(rt$Y,pch=20,col=as.numeric(conditions))

## mt-genes
mt.ratio <- colSums(nmat[mt.inds,])
mt.thres <- max(mt.ratio)
n.levels <- round(mt.ratio/mt.thres*11)
plot(rt$Y,pch=20,col=n.cols[n.levels])

## 
plot(rt$Y,pch=20,col=n.cols[n.levels])
