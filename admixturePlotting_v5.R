
source('~/Dropbox/pipelines2020/plot_admixture_v5_function.R')

# assembling the input table
dir="~/Dropbox/pipelines2020/" # path to input files
inName="ok_k2.qopt" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run
pops="inds2pops" # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.

#------------

npops=as.numeric(sub(".+(\\d+)\\..+","\\1",inName))
tbl=read.table(paste(dir,inName,sep=""),header=F)
i2p=read.table(paste(dir,pops,sep=""),header=F)

names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
row.names(tbl)=tbl$ind
head(tbl)
# putting populaitons in desired order (edit pop names as needed or skip to plot them alphabetically)
# tbl$pop=factor(tbl$pop,levels=c("O","K"))

ords=plotAdmixture(data=tbl,npops=npops,grouping.method="distance")

# recording cluster affiliations
cluster.admix=apply(tbl[,1:npops],1,function(x) {return(paste(which(x>0.33),collapse=".")) })
save(cluster.admix,file=paste(inName,"_clusters.RData",sep=""))
