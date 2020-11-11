# installing required packages: do this once, then remark with #
#install.packages("vegan")
#install.packages("ggplot2")
# install.packages("RcppCNPy")
library(RcppCNPy)

setwd('~/Dropbox/pipelines2020') # change this to where your scp'd files are
bams=read.table("bams.nc")[,1] # list of bam files
bams0=bams
# removing leading / trailing filename bits from sample names
bams=sub(".bam","",bams)
bams=sub(".+/","",bams)
length(bams)

#------ PCAngsd magic

# kinship (from pcangsd)
kin= npyLoad("pcangsd.kinship.npy") 
dimnames(kin)=list(bams, bams)
pheatmap::pheatmap(kin)
plot(hclust(as.dist(1-kin),method="ave"))

# inbreeding (from pcangsd)
inbr= npyLoad("pcangsd.inbreed.npy") 

# reading pcangsd covariance matrix, converting to correlation-based distances
pcc = as.matrix(read.table("pcangsd.cov"))
dimnames(pcc)=list(bams,bams)
pccd=1-cov2cor(pcc)

# must run admixturePlotting_v5.R first for the following to work
load('pcangsd_clusters.RData')

# reading the ibs matrix
ma = as.matrix(read.table("OK.ibsMat"))
dimnames(ma)=list(bams,bams)

# population designations (first letter of sample names)
pops=substr(bams,0,1)

# ---- sanity checks

# heatmaps of distances
library(pheatmap)
pheatmap(ma)
pheatmap(kin)
pheatmap(pccd)
# which one looks more structured?

# hierarchical clustering trees
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.8) # clustering of samples by IBS 
hc=hclust(as.dist(pccd),"ave")
plot(hc,cex=0.8) # clustering of samples by SNP correlation (from pcangsd) 

#----- vegan stuff (on pcangsd-derived distances)

library(vegan)

# PCoA
ordi=capscale(ma~1)
#ordi=capscale(pccd~1)

# eigenvectors
plot(ordi$CA$eig) 
# how many "interesting" ordination axes do we see?

# plotting vegan way (I really like ordispider and ordiellipse functions)

# extracting "scores" table, to plot
axes2plot=c(1,2) # which MDS to plot
scores=data.frame(scores(ordi,display="sites",choices=axes2plot))

# making color scheme for our groups (colramp trick makes it work with >12 groups)
library(RColorBrewer)
groups2color=cluster.admix
#groups2color=pops
colramp = colorRampPalette(brewer.pal(length(unique(groups2color)),"Set2"))
myColors=colramp(length(unique(groups2color)))
names(myColors) = sort(unique(groups2color))

plot(scores[,1:2],pch=16,col=myColors[groups2color],asp=1)
ordispider(scores[,1:2],group= groups2color,col=myColors)
ordiellipse(scores[,1:2],group= groups2color,col=myColors,draw="polygon",label=T)

# ------- plotting in ggplot2 to visualize if sequencing quality (coverage) affects us

# importing and aligning coverage data (from quality.txt file produced by plotQC.R)
cover=read.table("quality.txt")
cover$V1=sub(".bam","",cover$V1)
row.names(cover)=sub(".+/","",cover$V1)
cover$V1=NULL
cover=cover[bams,]

library(ggplot2)

# colored by coverage
ggplot(scores,aes(scores[,1],scores[,2],color=cover))+scale_color_continuous(type="viridis")+geom_point()+coord_equal()+theme_bw()

# colored by inbreeding
ggplot(scores,aes(scores[,1],scores[,2],color=inbr))+scale_color_continuous(type="viridis")+geom_point()+coord_equal()+theme_bw()

# colored by admixture clusters
ggplot(scores,aes(scores[,1],scores[,2],color=cluster.admix))+geom_point()+coord_equal()+theme_bw()

# colored by sampling site
ggplot(scores,aes(scores[,1],scores[,2],color=pops))+geom_point()+coord_equal()+theme_bw()

#--------- statistical tests (PERMANOVA)

# does coverage and population designation aligns with our "interesting" axes? (use this function to check is any parameter aligns with ordination) 

# assembling data frame of variables to corelate with ordination (cover and pops)
pc=data.frame(cbind(pops,cover),stringsAsFactors=F)
pc$cover=as.numeric(pc$cover)

# fitting them to chosen ordination axes (let's take the first 3)
ef=envfit(ordi,pc,choices=c(1:3))
ef
# note that pop loads on MDS1, cover - on MDS3 (can you make a plot to verify that?). 

# let's add fitted parameters to ourordination
axes2plot=c(1,2) # which PCAs to plot
scores=data.frame(scores(ordi,display="sites",choices=axes2plot))
plot(scores[,1:2],pch=16,col=myColors[groups2color],asp=1)
ordispider(scores[,1:2],group= groups2color,col=myColors)
ordiellipse(scores[,1:2],group= groups2color,col=myColors,draw="polygon",label=T)
plot(ef,choices=axes2plot)
# quantitative variables (cover in this case) plot as vectors, fitted groups (pops) plot as group-centroid locations

# How much *overall* variation is attributable to sampling site and coverage? Are they significant factors?
adonis(ma~cover+pops,pc)
adonis(pccd~cover+pops,pc)
# factor "pops" explains 3% of variation and is significant. 
# factor "cover" is less important but is significant, too.
 
# -------- can we remove effect of coverage from ordination?

# partial ordination - removing effect of coverage
pp1=capscale(pccd~1+Condition(cover))

# let's plot this the vegan way, like before
axes2plot=c(1,2) 
scores=data.frame(scores(pp1,display="sites",choices=axes2plot))
plot(scores[,1:2],pch=16,col=myColors[groups2color],asp=1)
ordispider(scores[,1:2],group= groups2color,col=myColors)
ordiellipse(scores[,1:2],group= groups2color,col=myColors,draw="polygon",label=T)

# how much the ordination actually changed? procrustes function
plot(procrustes(ordi,pp1))
# arrows show where the points moved in pp1 compared to ordi

# are the two ordinations sigificantly similar? Using procrustes test (protest)
# (of course they are but just to show how it is done. Can use this function instead of mantel test, too)
protest(ordi,pp1,permutations=9999)


