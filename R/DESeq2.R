library(DESeq2,lib.loc="/lustre/rdi/user/songx/tools/anaconda/lib/R/library")
library(getopt,lib.loc="/lustre/rdi/user/songx/tools/anaconda/lib/R/library")

opt_spec <- matrix(c(
	'help','h', 0, 'logical', 'help manual'
	,'rawcount', 'r', 1, 'character', 'the raw count file'
	,'condition', 'c', 1, 'character', 'the condition file'
	,'compare', 'v', 1, 'character', 'the compare name'
	,'design', 'd',1, 'character', 'the design of experiment'
	,'foldchange', 'f', 1, 'numeric', 'the foldchange value'
	,'pvalue', 'p', 1, 'numeric', 'the p value'
	,'padj','q', 1, 'numeric', 'the p adjust value'
	,'genename','g', 1, 'character', 'the genename file'
	,'outdir', 'o', 1, 'character','the directory of output file'
),byrow=TRUE, ncol=5)

opt = getopt(opt_spec, commandArgs(TRUE))

if (!is.null(opt$help)){
         cat(getopt(opt_spec, usage=TRUE))
         q(save='no', status=1)
}

#stopifnot(file.exists(opt$rawcount,opt$condition,opt$compare,opt$design,opt$foldchange,opt$pvalue,opt$padj,opt$genename,opt$outdir))

setwd(opt$outdir)

if ( !is.na(opt$genename) ) {Genename <- read.table(opt$genename,header=T,row.names=1)}

compares <- strsplit(opt$compare,'vs')[[1]]
condition <- read.delim(opt$condition, header=T)

number <- ncol(condition)

if ( colnames(condition)[number]=='types' ) { ngroup <- number-2}
if ( colnames(condition)[number]!='types' ) { ngroup <- number-1}

if ( ngroup==1 ) {
	groupdata1 <- condition[which(condition$groups1 %in% compares),]
	colnames(groupdata1)[2] <- 'groups'
	groupdata <- groupdata1}
if ( ngroup==2 ) {
	groupdata1 <- condition[which(condition$groups1 %in% compares),-3]
	groupdata2 <- condition[which(condition$groups2 %in% compares),-2]
	colnames(groupdata1)[2] <- 'groups'
	colnames(groupdata2)[2] <- 'groups'
	groupdata <- rbind(groupdata1,groupdata2)
}
if ( ngroup==3 ) {
	groupdata1 <- condition[which(condition$groups1 %in% compares),c(-3,-4)]
	groupdata2 <- condition[which(condition$groups2 %in% compares),c(-2,-4)]
	groupdata3 <- condition[which(condition$groups3 %in% compares),c(-2,-3)]
	colnames(groupdata1)[2] <- 'groups'
	colnames(groupdata2)[2] <- 'groups'
	colnames(groupdata3)[2] <- 'groups'
	groupdata <- rbind(groupdata1,groupdata2,groupdata3)}
if (ngroup==4 ) {
	groupdata1 <- condition[which(condition$groups1 %in% compares),c(-3,-4,-5)]
	groupdata2 <- condition[which(condition$groups2 %in% compares),c(-2,-4,-5)]
	groupdata3 <- condition[which(condition$groups3 %in% compares),c(-2,-3,-5)]
	groupdata4 <- condition[which(condition$groups4 %in% compares),c(-2,-3,-4)]
	colnames(groupdata1)[2] <- 'groups'
	colnames(groupdata2)[2] <- 'groups'
	colnames(groupdata3)[2] <- 'groups'
	colnames(groupdata4)[2] <- 'groups'
	groupdata <- rbind(groupdata1,groupdata2,groupdata3,groupdata4)}
if ( ngroup==5 ) {
	groupdata1 <- condition[which(condition$groups1 %in% compares),c(-3,-4,-5,-6)]
	groupdata2 <- condition[which(condition$groups2 %in% compares),c(-2,-4,-5,-6)]
	groupdata3 <- condition[which(condition$groups3 %in% compares),c(-2,-3,-5,-6)]
	groupdata4 <- condition[which(condition$groups4 %in% compares),c(-2,-3,-4,-6)]
	groupdata5 <- condition[which(condition$groups5 %in% compares),c(-2,-3,-4,-5)]
	colnames(groupdata1)[2] <- 'groups'
	colnames(groupdata2)[2] <- 'groups'
	colnames(groupdata3)[2] <- 'groups'
	colnames(groupdata4)[2] <- 'groups'
	colnames(groupdata5)[2] <- 'groups'
	groupdata <- rbind(groupdata1,groupdata2,groupdata3,groupdata4,groupdata5)}

rownames(groupdata) <- groupdata[,1]
countdata <- read.table(opt$rawcount,header=T,row.names=1)
countdata <- round(countdata)
countdata <- na.omit(countdata)
countdata <- countdata[rowSums(countdata)>1,]
countdata <- subset(countdata,select=rownames(groupdata))

if ( opt$design =='normal' ) {
	dds <- DESeqDataSetFromMatrix(countData=countdata,colData=groupdata,design=~groups)}
if ( opt$design=='pair' ) {
	dds <- DESeqDataSetFromMatrix(countData=countdata,colData=groupdata,design=~type+groups)}

dds <- DESeq(dds)
result <- results(dds,contrast=c('groups',compares[1],compares[2]))
ID <- rownames(result)

NormalizedCount=counts(dds,normalized=TRUE)

if ( is.na(opt$genename) ) {
	result <- cbind(ID,NormalizedCount,as.data.frame(result))
}
if ( !is.na(opt$genename) ) {
	Anno <- Genename[ID,]
	result <- cbind(ID,NormalizedCount,as.data.frame(result),Anno)
}

result$padj[is.na(result$padj)]  <- 1
result <- result[order(result$pvalue),]
result <- subset(result,select=-c(lfcSE,stat,baseMean))
NormalizedCount <- cbind(ID=rownames(NormalizedCount),NormalizedCount)

if ( is.na(opt$pvalue) ) {
	ALL <- subset(result,padj <= opt$padj & abs(log2FoldChange) >= log(opt$foldchange,2))
        UP <- subset(ALL,log2FoldChange > log(opt$foldchange,2))
	DOWN <- subset(ALL,log2FoldChange < -log(opt$foldchange,2))}
if ( is.na(opt$padj) ) {
	ALL <- subset(result,pvalue <= opt$pvalue & abs(log2FoldChange) >= log(opt$foldchange,2))
	UP <- subset(ALL,log2FoldChange > log(opt$foldchange,2))
	DOWN <- subset(ALL,log2FoldChange < -log(opt$foldchange,2))}
 
write.table(result,file=paste(opt$outdir,'/',opt$compare,'/',opt$compare,'_DEG.txt',sep=''),sep='\t',quote=F,row.names=F)
write.table(ALL,file=paste(opt$outdir,'/',opt$compare,'/',opt$compare,'_DEG_ALL.txt',sep=''),sep='\t',quote=F,row.names=F)
write.table(UP,file=paste(opt$outdir,'/',opt$compare,'/',opt$compare,'_DEG_UP.txt',sep=''),sep='\t',quote=F,row.names=F)
write.table(DOWN,file=paste(opt$outdir,'/',opt$compare,'/',opt$compare,'_DEG_DOWN.txt',sep=''),sep='\t',quote=F,row.names=F)
write.table(NormalizedCount,file=paste(opt$outdir,'/','Normalized.txt',sep=''),sep='\t',quote=F,row.names=F)
