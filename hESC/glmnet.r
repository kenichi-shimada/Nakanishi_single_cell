setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/monkey")
class <- readRDS(file="cell-lineage.rds")

params <- expand.grid(nrepeat=1:10,
	alpha=seq(0,1,0.1),
	clust=class,
	stringsAsFactors=F) ## 9900

nfolds <- 10
nlambda <- 100

is.mem <- rownames(x) %in% names(thisw)
y <- as.numeric(is.mem)
pos <- sum(y==1)
neg <- sum(y==0)

set.seed(12345) 
foldid <- lapply(seq(nrepeat),function(x){
	ind <- rep(NA,length(y))
	pos.ind <- sample(rep(seq(nfolds),length.out=pos),replace=F)
	neg.ind <- sample(rep(seq(nfolds),length.out=neg),replace=F)
	ind[y==1] <- pos.ind
	ind[y==0] <- neg.ind
	return(ind)
})[[nrepeat]]

## weights
weights <- rep(1,nrow(x))
names(weights) <- rownames(x)
weights[names(thisw)] <- thisw
weights[y==0] <- sum(thisw)/neg

##
set.seed(12345) 
cv <- cv.glmnet(x,y,alpha=alpha,nfolds=nfolds,nlambda=nlambda,foldid=foldid,
	family="binomial",weights=weights,type.measure="auc",parallel=FALSE,keep=TRUE)

setwd("/groups/mitchison/Kenichi/toxicogenomics/tggates/analysis/pathol_classifier_050717")
saveRDS(cv,file=paste0("glmnet_050717_",i,".rds"))