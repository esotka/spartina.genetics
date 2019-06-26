## PCA plots of each popn (short vs tall)
## using big dataset

#rm(list=ls())
pdf('output/shortVStall-pca-allDATA.pdf',width=8,height=12)
par(mfrow=c(3,2))
library(RColorBrewer)

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

for (i in 1:6)
{
  tmp <- gmat[,pop==sh[i] | pop==ta[i]]
  pop.tmp <- substr(colnames(tmp),1,3); pop.tmp <- factor(pop.tmp)
  pc.cr <- prcomp(t(tmp))
  col.sub <- c("red","black")[pop.tmp]
  plot(pc.cr$x[,1],pc.cr$x[,2],cex=0,main=site[i],xlab="PCA axis 1",ylab="PCA Axis 2")
  points(pc.cr$x[,1],pc.cr$x[,2],pch=20,cex=2,col=col.sub)
  for (i in 1:2){
    xbar.x <- mean(pc.cr$x[pop.tmp==levels(pop.tmp)[i],1])
    xbar.y <- mean(pc.cr$x[pop.tmp==levels(pop.tmp)[i],2])
    tmp.pc.cr <- pc.cr$x[pop.tmp==levels(pop.tmp)[i],]
    segments(x0 = xbar.x,y0 = xbar.y,x1 = tmp.pc.cr[,1],y1 = tmp.pc.cr[,2],
             col=c("red","black")[i],lwd=.5)
    }
  #print(stripplot(pc.cr$x[,1]~pop.tmp))
}

dev.off()
