### make histograms of k=2 for Short vs tall.
rm(list=ls())
#quartz(width=6,height=6)
pdf('output/q_Short.vs.Tall.pdf',width=6,height=6)
par(mfrow=c(3,2),mar=c(1,3,3,1))
library(lattice)
k=2
f <- c("k02run01_FLSFLT.qopt",
       "k02run01_SFBTFB.qopt",
       "k02run01_SBITBI.qopt",
       "k02run01_RISRIT.qopt",
       "k02run01_WESWET.qopt",
       "k02run01_SWSSWT.qopt")
f <- f[6:1] # north to south
fid <- c("inds.FL","inds.SFBTFB","inds.SBITBI","inds.RISRIT","inds.WESWET","inds.SWSSWT")
fid <- fid[6:1] # north to south
site <- c("Sweeney, MA","West Marsh, MA","Cole State Park, RI",
          "Bowens Island, SC", "Folly Island, SC","St. Joes Bay, FL")
#tim6equal = c("#00008F", "#005AFF", "#23FFDC", "#ECFF13", "#FF4A00", "#800000")
col.order <- list(c(1,2),c(1,2),c(1,2),c(2,1),c(2,1),c(1,2))
all <- c()
for (i in 1:6)
{
  tmp <- read.delim(paste("data/",f[i],sep=""),sep=" ",header = F)[,-(k+1)]
  tmpids <- read.delim(paste("data/",fid[i],sep=""),sep=" ",header = F)
  site1 <- substr(as.character(tmpids[,1]),1,3)
  all <- rbind(all,data.frame(ind=tmpids,site1,site2=site[i],q1=tmp$V1,q2=tmp$V2))
  barplot(t(tmp),col=c("black","grey")[col.order[[i]]],border=NA,ylim=c(-.1,1.1),
          names.arg = rep("",dim(tmp)[1]),space = 0)
  #box()
  mtext(site[i])
  ### sites
  for(j in 1:2)
  {
    x <- 1:dim(tmp)[1]
    xtmp <- x[site1==unique(site1)[j]]
    segments(max(xtmp),1,max(xtmp),-.05,lwd=3,col="white")
    mtext(side=1,at=mean(xtmp),unique(site1)[j],cex=.75,line=0)
  }
}

dev.off()




#sh <- c("SWS","WES","RIS","SBI","SFB","FLS")#,"FJS")
#ta <- c("SWT","WET","RIT","TBI","TFB","FLT")#,"SCT")
#site <- c("Sweeney, MA","West Marsh, MA","Cole State Park, RI",
#          "Bowens Island, SC", "Folly Island, SC","St. Joes Bay, FL")


