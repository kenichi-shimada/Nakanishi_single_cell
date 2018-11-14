library(openxlsx)

setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/non-normalization")
load("non-normalization.rda")

## biomarkers
setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/biomarkers")
bm.list <- read.xlsx("Bhatia_AML_Gene_lists_revised.xlsx",1)
reg <- bm.list$Leukemic.regeneration
stem <- bm.list$Leukemic.stem.cell
stem <- stem[!is.na(stem)]
# stem <- stem[!stem %in% "GPR56"]
lsc17 <- bm.list$LSC
lsc17 <- lsc17[!is.na(lsc17)]

syms <- gns$symbol
syms <- sub("^hg19_","",syms)
bms <- list(regeneration=reg,stem.cell=stem,lsc17=lsc17)

## Are they found in matrix (hg19)?
exists <- lapply(bms,function(l){
	sapply(l,function(x){
		query <- paste0("^",x,"$")
		any(grepl(query,syms))
	})
})
sapply(exists,table)[2:1,]
#       reg stem lsc17
# TRUE  175   38    16
# FALSE   7    3     1

##
bms.exist <- list(regeneration=reg[exists[[1]]],
	stem.cell=stem[exists[[2]]],
	lsc17=lsc17[exists[[3]]])
bms.exist[[2]] <- bms.exist[[2]]

bms.inds <- lapply(bms.exist,function(l){
	sapply(l,function(x){
		query <- paste0("^",x,"$")
		grep(query,syms)
	})
})

library(limma)
bms.all <- unique(unlist(bms.exist))
mat.ind <- cbind(regeneration=bms.all %in% bms.exist$reg,
	stem.cell=bms.all %in% bms.exist$stem,
	lsc17=bms.all %in% bms.exist$lsc17)
vennDiagram(mat.ind,main="biomarker overlap")
intersect(bms.exist$stem,bms.exist$lsc17) ##GPR56 is overlapped

## remove duplicated gene: GPR56
bms.exist$stem,bms.exist$lsc17

##
multi.mapped <- lapply(bms.inds,function(l){
	l[sapply(l,length)>1]
}) ## MUC3A, CCDC177

uniq.inds <- unlist(bms.inds)
syms.inds <- tapply(uniq.inds,syms[uniq.inds],identity)

## overlap
bms.exp <- t(sapply(syms.inds,function(inds){
	apply(mat[inds,,drop=F],2,sum)
}))

bms.nexp <- t(sapply(syms.inds,function(inds){
	apply(nmat[inds,,drop=F],2,sum)
}))

colnames(bms.exp) <- colnames(bms.nexp) <- colnames(mat)

## Number of expressed cells per each biomarker gene per gene set
n.expressed.cells.per.bm.list <- lapply(bms.exist,function(bms){
	t(apply(bms.nexp[bms,],1,function(x){
		tapply(x,conditions,function(y)sum(y>0))
	}))
})

tmp <- do.call(rbind,n.expressed.cells.per.bm.list)
tmp <- cbind(tmp,all=rowSums(tmp))
rownames(tmp)[230] <- "GPR56.1"
table(table(rownames(tmp))) ## all unique entries now

kinds <- factor(rep(names(bms.exist),sapply(bms.exist,length)),
	levels=c("regeneration","stem.cell","lsc17"))
df <- cbind(kind=kinds,tmp)

if(0)write.csv(df,"biomarker_n_expressed_cells.csv")

is.expressed.in.any.cells <- tmp[,"all"]>0
table(kinds,is.expressed.in.any.cells)[,2:1]

##

if(0)write.csv(df[is.expressed.in.any.cells,],"expressed_biomarker_n_expressed_cells.csv")

## Number of expressed genes per cell per condition (1 or 2)
cells <- seq(ncol(mat))
n.expressed.genes.per.bm.per.cell <- tapply(cells,conditions,function(i){
	sapply(bms.exist,function(this.bms){
		apply(bms.nexp[this.bms,i],2,function(y)sum(y>0))
	})
})
tmp.1 <- do.call(rbind,n.expressed.genes.per.bm.per.cell)
tmp.1 <- cbind(tmp.1,all=rowSums(tmp.1))
if(0)
	table(conditions,expressing.any.genes)

n.reads <- colSums(mat)

if(0){
	##
	png("01_total-reads_vs_# of biomarkers expressed.png",width=800,height=400)
	par(mfrow=c(1,2))
	plot(n.reads,tmp.1[,"all"],type="n",main="condition: 1",
		xlab="# total reads/cells",ylab="# biomarker genes observed")
	points(n.reads[conditions=="2"],tmp.1[conditions=="2","all"],col="grey80",pch=20,cex=1)
	points(n.reads[conditions=="1"],tmp.1[conditions=="1","all"],pch=20,col=1)

	plot(n.reads,tmp.1[,"all"],type="n",main="condition: 2",
		xlab="# total reads/cells",ylab="# biomarker genes observed")
	points(n.reads[conditions=="1"],tmp.1[conditions=="1","all"],pch=20,col="grey80")
	points(n.reads[conditions=="2"],tmp.1[conditions=="2","all"],pch=20,col=2)
	dev.off()

	lm1 <- lm(tmp.1[conditions=="1","all"]~n.reads[conditions=="1"])
	lm2 <- lm(tmp.1[conditions=="2","all"]~n.reads[conditions=="2"])

	png("01_smoothscatter_total-reads_vs_# of biomarkers expressed.png",width=800,height=400)
	par(mfrow=c(1,2))
	plot(n.reads,tmp.1[,"all"],type="n",main="condition: 1",
		xlab="# total reads/cells",ylab="# biomarker genes observed")
	smoothScatter(n.reads[conditions=="1"],tmp.1[conditions=="1","all"],pch=20,col=1,add=T,cex=.5)
	abline(lm1,col=2)
	abline(lm2,col="grey60",lty=2)

	plot(n.reads,tmp.1[,"all"],type="n",main="condition: 2",
		xlab="# total reads/cells",ylab="# biomarker genes observed")
	smoothScatter(n.reads[conditions=="2"],tmp.1[conditions=="2","all"],pch=20,col=1,add=T,cex=.5)
	abline(lm1,col="grey60",lty=2)
	abline(lm2,col=2)
	dev.off()

	df.1 <- cbind(condition=conditions,tmp.1)
	write.csv(df.1,"biomarker_n_expressed_genes_per_list_per_cells.csv")

	##

	n.genes <- tmp.1[,"all"]

	thres <- 25
	n.genes[n.genes > thres] <- thres

	hist(n.genes)

	##
	n.cols <- rev(colorRampPalette(brewer.pal(11,"YlGnBu"))(thres+1))
	n.leg.cols <- rev(brewer.pal(5,"YlGnBu"))

	n.levels <- n.genes + 1
	png("02 number of biomarkers expressed on tsne.png",width=400,height=400)
	par(mar=c(4,5,4,3))
	plot(rt$Y,pch=20,col=n.cols[n.levels],
		main="# of biomarkers expressed per cell",
		xlab="tSNE 1st",ylab="tSNE 2nd")
	legend("topright",c("low",rep("",3),"high"),pch=20,col=n.leg.cols,cex=.8)
	dev.off()

	png("02 number of biomarkers expressed on tsne_combined.png",width=800,height=400)
	par(mfrow=c(1,2))
	par(mar=c(4,5,4,3))
	plot(rt$Y,pch=20,col=as.numeric(conditions),
		main="conditions",
		xlab="tSNE 1st",ylab="tSNE 2nd",cex=.7)
	legend("topright",c("control","treated"),pch=20,col=1:2)
	plot(rt$Y,pch=20,col=n.cols[n.levels],
		main="# of biomarkers expressed per cell",
		xlab="tSNE 1st",ylab="tSNE 2nd",cex=.7)
	legend("topright",c("low",rep("",3),"high"),pch=20,col=n.leg.cols,cex=.8)
	dev.off()
}


library(dplyr)
library(ggplot2)
library(tidyr)

ratio.exp.cells <- tmp[,1:2]/rep(table(conditions),each=nrow(tmp))
png("03 ratio.number of cells expressed per gene.png",width=400,height=400)
plot(ratio.exp.cells,pch=20,asp=1,,xlim=c(0,1),ylim=c(0,1),
	main="ratio of # cells expressed per gene",
	xlab="in control (condition:1)",ylab="in experiment (condition:2)")
abline(h=0,v=0,col="grey60")
abline(a=0,b=1,lty=2,col="grey60")
cols1 <- brewer.pal(3,"Set1")
points(ratio.exp.cells[,1:2],pch=20,asp=1,col=cols1[kinds])
legend("bottomright",levels(kinds),col=cols1,pch=20,bg="white")
dev.off()

png("03 ratio.number of cells expressed per gene_split.png",width=1200,height=400)
par(mfrow=c(1,3))
for(i in 1:3){
	plot(ratio.exp.cells,pch=20,asp=1,,xlim=c(0,1),ylim=c(0,1),type="n",
		main="ratio of # cells expressed per gene",
		xlab="in control (condition:1)",ylab="in experiment (condition:2)")
	abline(h=0,v=0,col="grey60")
	abline(a=0,b=1,lty=2,col="grey60")
	cols1 <- rep(NA,3)
	cols1[i] <- brewer.pal(3,"Set1")[i]
	points(ratio.exp.cells[,1:2],pch=20,asp=1,col=cols1[kinds])
	legend("bottomright",levels(kinds)[i],col=cols1[i],pch=20,bg="white")
}
dev.off()

##

df.1 <- data.frame(condition=conditions,tmp.1) %>% gather("class","n.genes",2:4)
df.1$class <- factor(df.1$class,levels=levels(kinds))
levels(df.1$condition) <- c("1 (control)","2 (experiment)")
p.1 <- ggplot(df.1,aes(class,n.genes)) + 
	geom_violin(aes(fill=factor(condition)),scale="count") +
	ggtitle("# of expressed genes per cell")
png("04 ratio.number of genes expressed per cell.png",width=400,height=400)
p.1
dev.off()

boxplot(do.call(c,apply(tmp.1,2,function(x)tapply(x,conditions,identity))))

ratio.n.exp.genes.per.cell <- tmp.1[,1:3]/rep(n.reads,times=3)
hist(ratio.n.exp.genes.per.cell)

png("04 ratio. number of genes expressed per cell.png",width=400,height=400)
plot(ratio.n.exp.genes.per.cells,pch=20,asp=1,,xlim=c(0,1),ylim=c(0,1),
	main="ratio of # cells expressed per gene",
	xlab="in experiment (condition:1)",ylab="in control (condition:2)")
abline(h=0,v=0,col="grey60")
abline(a=0,b=1,lty=2,col="grey60")
cols1 <- brewer.pal(3,"Set1")
points(ratio.exp.cells[,1:2],pch=20,asp=1,col=cols1[kinds])
legend("bottomright",levels(kinds),col=cols1,pch=20,bg="white")
dev.off()

##
if(0){
	is.expressed.in.any.cells <- tmp[,"all"]>0
	expressing.genes.ge5 <- tmp.1[,"all"]>=10
	hist(tmp.1[,"all"])

	library(gplots)
	library(RColorBrewer)
	cols <- brewer.pal(9,"YlGnBu")
	cols1 <- brewer.pal(3,"Set1")
	cols2 <- brewer.pal(3,"Set2")

	sub.bms.nexp <- bms.nexp[is.expressed.in.any.cells,expressing.any.genes] # 96 genes x 4996 cells
	png("heatmap.bms.nexp.png",width=500,height=500)
	heatmap.2(sub.bms.nexp,col=cols,
		RowSideColors=cols1[as.numeric(kinds[is.expressed.in.any.cells])],
		ColSideColors=cols2[as.numeric(conditions[expressing.any.genes])],
		trace="none")
	dev.off()
}

png("05 number of biomarkers expressed on tsne_split.png",width=800,height=800)
par(mfrow=c(2,2))
##
par(mar=c(4,5,4,3))
plot(rt$Y,pch=20,col=as.numeric(conditions),
	main="conditions",
	xlab="tSNE 1st",ylab="tSNE 2nd",cex=.7)
legend("topright",c("control","treated"),pch=20,col=1:2)

for(k in levels(kinds))	{
	n.genes <- tmp.1[,k]
	o <- order(n.genes,decreasing=F)
	if(k==levels(kinds)[1]){
		thres <- 4
	}else if(k==levels(kinds)[2]){
		thres <- quantile(n.genes,.99) ## 13
	}else if(k==levels(kinds)[3]){
		thres <- quantile(n.genes,.99) ## 8
	}
	n.genes[n.genes > thres] <- thres

	##
	n.cols <- c("grey80",rev(colorRampPalette(brewer.pal(9,"YlGnBu"))(13)))
	n.leg.cols <- c("grey80",rev(brewer.pal(4,"YlGnBu")))

	n.levels <- n.genes + 1

	plot(rt$Y[o,],pch=20,col=n.cols[n.levels][o],
		main=paste0("# of biomarkers expressed per cell (",k,")"),
		xlab="tSNE 1st",ylab="tSNE 2nd",cex=.7)
	legend("topright",c(0,"low",rep("",2),"high"),pch=20,col=n.leg.cols,cex=.8)
}
dev.off()

png("05 biomarkers expressed on tsne_1.png",width=400,height=400)

par(mar=c(4,5,4,3))
for(k in levels(kinds)[1])	{
	n.genes <- tmp.1[,k]
	o <- order(n.genes,decreasing=F)
	if(k==levels(kinds)[1]){
		thres <- 4
	}else if(k==levels(kinds)[2]){
		thres <- quantile(n.genes,.99) ## 13
	}else if(k==levels(kinds)[3]){
		thres <- quantile(n.genes,.99) ## 8
	}
	n.genes[n.genes > thres] <- thres

	##
	n.cols <- c("grey80",rev(colorRampPalette(brewer.pal(9,"YlGnBu"))(thres)))
	n.leg.cols <- c("grey80",rev(brewer.pal(4,"YlGnBu")))

	n.levels <- n.genes + 1

	plot(rt$Y[o,],pch=20,col=n.cols[n.levels][o],
		main=paste0("# of biomarkers expressed per cell (",k,")"),
		xlab="tSNE 1st",ylab="tSNE 2nd",cex=.7)
	legend("topright",c(0,"low",rep("",2),"high"),pch=20,col=n.leg.cols,cex=.8)
}
dev.off()

