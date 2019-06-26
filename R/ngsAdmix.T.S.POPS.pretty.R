rm(list=ls())
library(fields)
pdf('output/ngsAdmix.T.S.POPS.pretty.pdf',width=8,height=8)
par(mfrow=c(1,16),mar=c(0,0,0,0),xpd = TRUE)
ids <- read.delim('data/indIDS.T.S.POPS.txt',header = F)[,1]
pop <- substr(ids,1,3)
meta <- read.csv('data/Spartina_SNP_SiteID.csv')
state <- meta$State[match(pop,meta$Site_ID)]
siteorder.txt <-     c("FLS","FLT","SFB","TFB","SBI","TBI",
                       "RIS","RIT","SWS","SWT","WES","WET")
siteorderNICE.txt <- c("FLS","FLT","SCFS","SCFT","SCBS","SCBT",
                       "RIS","RIT","MASS","MAST","MAWS","MAWT")
siteorder <- match(pop,siteorder.txt)
pop <- pop[order(siteorder)]
stateorder.txt <- c("FL","SC","RI","MA")
stateorder <- match(state,stateorder.txt)
state <- state[order(siteorder)]

filelist <- c("k2run2.qopt",
              "k3run2.qopt",
              "k4run2.qopt",
              "k5run2.qopt",
              "k6run2.qopt",
              "k7run2.qopt",
              "k8run2.qopt",
              "k9run2.qopt",
              "k10run2.qopt",
              "k11run2.qopt",
              "k12run2.qopt",
              "k13run2.qopt",
              "k14run2.qopt",
              "k15run2.qopt",
              "k16run2.qopt",
              "k17run2.qopt",
              "k18run2.qopt",
              "k19run2.qopt",
              "k20run2.qopt")
ks = 2:16


kcol <- c("black",
          "red",
          "darkgreen",
          "gainsboro",
          "yellow",
          "deepskyblue",
          "brown",
          "dodgerblue4",
          "darkorchid",
          "burlywood",
          "aquamarine",
          "chocolate",
          "khaki",
          "darksalmon",
          "firebrick",
          "forestgreen")
colorder <- list(
  c(1,2),
  c(2,3,1), #ks=3
  c(3,2,1,4), #ks=4
  c(4,2,5,1,3), #ks=5
  c(2,6,3,1,4,5), #ks=6
  c(6,2,7,5,4,1,3), #ks=7
  c(8,1,3,7,4,2,6,5), #ks=8
  c(8,1,5,3,7,2,9,6,4), #ks=9
  c(1,8,4,9,10,5,6,3,2,7), #ks=10
  c(8,3,5,1,11,7,6,2,4,9,10), #ks=11
  c(10,9,5,2,1,4,3,12,11,6,7,8), #ks=12
  c(9,13,3,12,2,11,1,7,8,6,4,5,10), #ks=13
  c(13,14,7,10,3,2,12,4,6,9,11,5,8,1), #ks=14
  c(4,13,7,11,15,14,2,1,3,8,10,9,5,12,6), #ks=15
  c(16,8,7,6,9,11,5,14,12,13,4,3,15,2,1,10) #ks=16
)

for (k in 1:length(ks))
{
  dat <- read.delim(paste("data/admix.runs.TvS.POPS/",filelist[k],sep=""),sep=" ",header = F)
  dat <- dat[order(siteorder),]
  dat <- dat[,-(dim(dat)[2])]
  dat <- dat[,order(colorder[[k]])]
  fig <- barplot(t(dat),col=kcol[1:length(colorder[[k]])],space=0,border=NA,xlab="",ylab="",
                 names.arg = rep("",nrow(dat)),horiz=T)#,main=substr(ks[k],1,3))
  mtext(substr(ks[k],1,3),side=3,line=-1.5,at=.5)
 
  ### lines
  for(i in 1:length(siteorder))
  {
    #pop2 <- pop[order(siteorder)]
    x <- 1:dim(dat)[1]
    xtmp <- x[pop==siteorder.txt[i]]
    segments(1,max(xtmp),-0.05,max(xtmp),col="white")
  }}
fig <- barplot(t(dat),col="white",space=0,border=NA,xlab="",ylab="",
               names.arg = rep("",nrow(dat)),horiz=T)

for(i in 1:length(siteorder))
{
  #pop2 <- pop[order(siteorder)]
  x <- 1:dim(dat)[1]
  xtmp <- x[pop==siteorder.txt[i]]
  text(y=mean(xtmp),x=.5,siteorderNICE.txt[i],cex=1)
}
dev.off()
