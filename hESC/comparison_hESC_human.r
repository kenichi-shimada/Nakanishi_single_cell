##
setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/hESC_data")
# hesc.exp.lin <- readRDS(file="uniq.exp.10clust.rds")
hesc.exp.lin <- readRDS(file="uniq.exp.7clust_062818.rds")

##
setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/hES/")
hs.exprs.lin <- readRDS(file="exprs.3lin_062818.rds")
hs.exprs.lin.1 <- readRDS(file="exprs.extranorm.3lin_062818.rds")

colnames(hs.exprs.lin) <- colnames(hs.exprs.lin.1) <- c("EPI","PE","TE")
##
if(0){
	overlapped <- intersect(rownames(hesc.exp.lin),rownames(hs.exprs.lin)) ## 18202 -> 17746
	setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/hES/")
	if(0)saveRDS(overlapped,file="overlapped_062818.rds")
	subset <- overlapped
}else{
	setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/hES/")
	bm.genes <- readRDS(file="biomarkers.hs1_062818.rds")
	subset <- 	bm.genes
}

hesc.exp.lin <- hesc.exp.lin[subset,]
hs.exprs.lin <- hs.exprs.lin[subset,]
hs.exprs.lin.1 <- hs.exprs.lin.1[subset,]

##
cor.mat <- apply(hesc.exp.lin,2,function(hs){
	apply(hs.exprs.lin.1,2,function(hs1){
		cors <- cor(hs,hs1,use="everything",method="spearman")
	})
})

library(gplots)
library(RColorBrewer)

cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(20)

sub.cor.mat <- cor.mat[c("EPI","PE","TE"),]
colnames(sub.cor.mat) <- 1:7

setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/hES/plots")
png("spearman_correlation_hs_all_genes_062818.png",width=500,height=300)
heatmap.2(sub.cor.mat,trace="none",margins=c(5,10),col=cols,keysize=2,Rowv=F,Colv=F,dendrogram="none")
dev.off()

pdf("spearman_correlation_hs_biomarkers_only_062818.pdf",width=5,height=5)
heatmap.2(sub.cor.mat,trace="none",margins=c(5,10),col=cols,keysize=1.5,Rowv=F,Colv=F,dendrogram="none")
dev.off()

setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/hES/plots")
png("spearman_correlation_hs_biomarkers_only_062818.png",width=500,height=300)
heatmap.2(sub.cor.mat,trace="none",margins=c(5,10),col=cols,keysize=2,
	Rowv=F,Colv=F,dendrogram="none")
dev.off()


##
if(0){
	setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/hES/plots")
	png("spearman_correlation_mf_hsESC_DEG_062818.png",width=500,height=500)
	heatmap.2(cor.mat,trace="none",margins=c(5,10),col=cols)
	dev.off()


	##
	setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/hES/plots")
	png("spearman_correlation_mf_hsESC_DEG_exclusive_up_062818.png",width=500,height=500)
	heatmap.2(cor.mat,trace="none",margins=c(5,10),col=cols)
	dev.off()

	setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/20180520_hESC/hES/plots")
	png("spearman_correlation_mf_hsESC_DEG_exclusive_up-woEXMC_062818.png",width=500,height=500)
	heatmap.2(cor.mat[!rownames(cor.mat) %in% "EXMC",],trace="none",margins=c(5,10),col=cols)
	dev.off()
}
