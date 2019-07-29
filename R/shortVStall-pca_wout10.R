## PCA plots of each popn (short vs tall)
### and histograms of k=2 for Short vs tall.

rm(list=ls())
pdf('output/shortVStall-pca_wout10.pdf',width=12,height=12)
par(mfrow=c(3,4))
library(RColorBrewer);library(reshape); library(vegan)

### PCA
meta <- read.csv('data/Spartina_SNP_SiteID.csv')
gprob<- read.table('data/spartinaNov2017.called.subset.mpgl')
ids <- read.csv('data/inds/allpops_individualIDs_subset.csv')[,2]
pop <- substr(ids,1,3)
reg <- meta$State[match(pop,meta$Site_ID)]
gmat <- gprob[,-1]
colnames(gmat) <- ids
sh <- c("SWS","WES","RIS","SBI","SFB","FLS")#,"FJS")
ta <- c("SWT","WET","RIT","TBI","TFB","FLT")#,"SCT")
site <- c("Sweeney, MA","West Marsh, MA","Cole State Park, RI",
          "Bowens Island, SC", "Folly Island, SC","St. Joes Bay, FL")
a <- meta$State[match(ta,meta$Site_ID)]

### Fst for short vs long
nloci = 2735
fst.TS <- c()
filelist <- c("SWS.SWT.pops.fst",
              "WES.WET.pops.fst",
              "RIS.RIT.pops.fst",
              "SBI.TBI.pops.fst",
              "SFB.TFB.pops.fst",
              "FLS.FLT.pops.fst")
# Columns are ordered as: A, AB, f, FST, Pvar; where A is the expectation of genetic variance between populations, 
# AB is the expectation of the total genetic variance, f is the correcting factor for the ratio of expectations, 
# FST is the per-site FST value, Pvar is the probability for the site of being variable.
for (i in 1:length(filelist))
{
tmp <- read.table(paste("data/fst/",as.character(filelist[i]),sep=""))
tmp <- tmp[1:nloci,]
tmp[is.na(tmp[,4]),4] <- 0
tmp[tmp[,4]<0,4] <- 0 ### set all FST <0 to 0
tmp[tmp[,4]>1,4] <- 1 ### convert FST > 1 to 1
fst.tmp <- tmp[,4]
fst.TS <- rbind(fst.TS,data.frame(sites=filelist[i],fst=fst.tmp,type = "TvS"))
}
fst.TS$sites <- factor(fst.TS$sites)
fst.TS$locusnum <- rep(1:2735,length(unique(fst.TS$sites)))
md <- melt(fst.TS,id.vars = c("sites","locusnum","type"),measure.vars = "fst")
### take away 10% outlier
md.90 <- c()
n <- unique(md$sites)

for(i in 1:length(n))
{
  tmp <- md[md$sites==n[i],]
  if(i%in%c(1:6)){
    threshold = sort(tmp$value)[round(length(tmp$value)*.9,0)] # 10% outlier FST
    tmp$value[tmp$value>=threshold] <- NA
  }
  md.90 <- rbind(md.90,data.frame(tmp[,c("locusnum","sites","type","value")]))
}
md.90$sites <- factor(md.90$sites)



out.alldata <- c(); out.90 <- c()
for (i in 1:6)
{

  ### PCA - all data
  tmp <- gmat[,pop==sh[i] | pop==ta[i]]
  pop.tmp <- substr(colnames(tmp),1,3); pop.tmp <- factor(pop.tmp)
  pc.cr <- prcomp(t(tmp))
  anosim.p=anosim(t(tmp),pop.tmp)$signif
  anosim.stat=anosim(t(tmp),pop.tmp)$statistic
  out.alldata <- rbind(out.alldata,data.frame(site[i],PC1=summary(eigenvals(pc.cr))[2,1],PC2=summary(eigenvals(pc.cr))[2,2],anosim.stat,anosim.p))
  col.sub <- c("red","black")[pop.tmp]
  par(mar=c(5,5,5,2))
  plot(pc.cr$x[,1],pc.cr$x[,2],cex=0,xlab="Axis 1",ylab="Axis 2")
  points(pc.cr$x[,1],pc.cr$x[,2],pch=20,cex=2,col=col.sub)
  for (j in 1:2){
    xbar.x <- mean(pc.cr$x[pop.tmp==levels(pop.tmp)[j],1])
    xbar.y <- mean(pc.cr$x[pop.tmp==levels(pop.tmp)[j],2])
    tmp.pc.cr <- pc.cr$x[pop.tmp==levels(pop.tmp)[j],]
    segments(x0 = xbar.x,y0 = xbar.y,x1 = tmp.pc.cr[,1],y1 = tmp.pc.cr[,2],
             col=c("red","black")[j],lwd=.5)
    }
  print(pc.cr$x[1:3,1:3])
  mtext(site[i],cex=1.5,line=2)
  mtext("All loci",line=.5)
  ### PCA - 90% outlier
  tmp <- gmat[,pop==sh[i] | pop==ta[i]]
  pop.tmp <- substr(colnames(tmp),1,3); pop.tmp <- factor(pop.tmp)
  loci.to.use <- md.90$locusnum[md.90$sites==filelist[i] & !is.na(md.90$value)]
  tmp <- tmp[loci.to.use,]  
  pc.cr <- prcomp(t(tmp))
  anosim.p=anosim(t(tmp),pop.tmp)$signif
  anosim.stat=anosim(t(tmp),pop.tmp)$statistic
  out.90 <- rbind(out.90,data.frame(site[i],PC1=summary(eigenvals(pc.cr))[2,1],PC2=summary(eigenvals(pc.cr))[2,2],anosim.stat,anosim.p))
  col.sub <- c("red","black")[pop.tmp]
  par(mar=c(5,5,5,2))
  plot(pc.cr$x[,1],pc.cr$x[,2],cex=0,xlab="Axis 1",ylab="Axis 2")
  points(pc.cr$x[,1],pc.cr$x[,2],pch=20,cex=2,col=col.sub)
  for (j in 1:2){
    xbar.x <- mean(pc.cr$x[pop.tmp==levels(pop.tmp)[j],1])
    xbar.y <- mean(pc.cr$x[pop.tmp==levels(pop.tmp)[j],2])
    tmp.pc.cr <- pc.cr$x[pop.tmp==levels(pop.tmp)[j],]
    segments(x0 = xbar.x,y0 = xbar.y,x1 = tmp.pc.cr[,1],y1 = tmp.pc.cr[,2],
             col=c("red","black")[j],lwd=.5)
  }
  print(pc.cr$x[1:3,1:3])
  mtext(site[i],cex=1.5,line=2)
  mtext("10% outliers removed",line=.5)
}
print(out.alldata)
print(out.90)
dev.off()
