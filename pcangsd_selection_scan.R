library(RcppCNPy)
library(ggplot2)

# ------ selection scan (from pcangsd)
# (along all important PCs. In this case there are only two pops inferred, so one important PC)
sel= npyLoad("pcangsd.selection.npy") # Reads results from selection scan

qqchi<-function(x,...){
lambda<-round(median(x)/qchisq(0.5,1),2)
  qqplot(qchisq((1:length(x)-0.5)/(length(x)),1),x,ylab="Observed",xlab="Expected",...);abline(0,1,col=2,lwd=2)
legend("topleft",paste("lambda=",lambda))
}

# qq plot for selection scan - is there inflation of high p-values?
qqchi(sel)

# pvalues - against chisq with df=1
pval=pchisq(sel,1,lower.tail=F)
p.adj=p.adjust(pval,method="BH")

sites=read.table("OK.mafs.gz",header=T)

# ---- making manhattan plot, point size by allele frequency, color by adjusted pval

mh=sites[,c(1,2,5)]
mh$pos.mb=mh$position/1e+6
names(mh)[1:3]=c("chrom","pos","maf")
mh$sel=sel
mh$logp=(-1)*log(pval,10)
mh$logpa=(-1)*log(p.adj,10)
mh$maf3=mh$maf^3

ggplot(mh,aes(pos.mb,logp))+
	geom_point(shape = 21, colour = "grey20", aes(size=maf3,fill=logpa))+
    scale_size_continuous(breaks=c(0.2,0.4,0.6)^3,labels=c(0.2,0.4,0.6))+
	scale_fill_gradient(low="grey80",high="coral")+
	theme_bw() + labs(size = "maf")+theme(axis.text.x=element_text(angle=45, hjust=1))+
	xlab("position,Mb")+facet_grid(~chrom,scale="free_x",space="free_x")

#-------- top selected sites

head(mh[order(mh$logp,decreasing=T),])
