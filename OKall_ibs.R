# installing required packages: do this once, then remark with #

setwd('~/Dropbox/pipelines2020') # change this to where your scp'd files are
bams=read.table("bams")[,1] # list of bam files
bams0=bams
# removing leading / trailing filename bits from sample names
bams=sub(".bam","",bams)
bams=sub(".+/","",bams)
length(bams)

ma = as.matrix(read.table("OKall.ibsMat"))
dimnames(ma)=list(bams,bams)

ma[1:6,1:6]

hc=hclust(as.dist(ma),"ave")

plot(hc,cex=0.8) # clustering of samples by IBS (great to detect clones or closely related individuals)

#-------- pruning clonal replicates

cutclones=0.2
abline(h=cutclones,col="red")
grps=cutree(hc,h= cutclones)

# retaining a single representative of each clonal group
pruned=c();i=1
for (i in unique(grps)) {
	pruned[i]=names(grps[grps==i])[1]
}
length(pruned)
ma1=ma[pruned,pruned]
hc=hclust(as.dist(ma1),"ave")
plot(hc,cex=0.8) 

# also removing both K4 and O5 = obviously there is some mixup going on with sample names
pruned=pruned[-which(pruned %in% c("K4","O5"))]
goodbams=bams0[which(bams %in% pruned)]
length(goodbams)

write.table(goodbams,file="bams.nc",quote=F, col.names=F, row.names=F)

# scp bam.nc to TACC, rerun angsd
