rm(list=ls())
library(fields)
pdf('output/ngsAdmix.pretty.pdf',width=8,height=8)
par(mfrow=c(1,16),mar=c(0,0,0,0),xpd = TRUE)
ids <- read.csv('data/allpops_individualIDs_subset.csv')[,2]
pop <- substr(ids,1,3)
meta <- read.csv('data/Spartina_SNP_SiteID.csv')
state <- meta$State[match(pop,meta$Site_ID)]
siteorder.txt <-     c("FLS","FLT","SFB","TFB","SBI","TBI","FJS","HWW","SCT","NCC",
                       "RIS","RIT","SWS","SWT","WES","WET","NHH")
siteorderNICE.txt <- c("FLS","FLT","SCFS","SCFT","SCBS","SCBT","SC1","SC2","SC3","NC1",
                       "RIS","RIT","MASS","MAST","MAWS","MAWT","NH1")
siteorder <- match(pop,siteorder.txt)
pop <- pop[order(siteorder)]
stateorder.txt <- c("FL","SC","NC","RI","MA","NH")
stateorder <- match(state,stateorder.txt)
state <- state[order(siteorder)]

filelist <- c("k03.qopt",
              "k04.qopt",
              "k05.qopt",
              "k06run01.qopt",
              "k07.qopt",
              "k08.qopt",
              "k09.qopt",
              "k10.qopt",
              "k11run01.qopt",
              "k12run01.qopt",
              "k13run01.qopt",
              "k14run01.qopt",
              "k15run01.qopt",
              "k16run01.qopt")
ks = 3:16







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
  c(1,2,3), #ks=3
  c(3,4,1,2), #ks=4
  c(5,1,3,4,2), #ks=5
  c(6,1,2,3,5,4), #ks=6
  c(4,3,5,2,7,1,6), #ks=7
  c(3,4,1,8,6,5,2,7), #ks=8
  c(9,6,2,3,4,7,5,1,8), #ks=9
  c(9,8,2,5,3,10,7,6,1,4), #ks=10
  c(11,4,9,8,1,6,5,2,3,7,10), #ks=11
  c(4,1,9,5,2,10,8,7,11,6,12,3), #ks=12
  c(1,5,4,12,11,2,3,8,7,6,9,10,13), #ks=13
  c(4,12,14,2,13,5,3,10,1,7,9,8,6,11), #ks=14
  c(1,8,11,4,12,3,2,13,10,5,15,7,14,9,6), #ks=15
  c(16,1,12,13,3,5,10,8,7,9,2,4,15,14,11,6) #ks=16
)

for (k in 1:length(ks))
{
  dat <- read.delim(paste("data/",filelist[k],sep=""),sep=" ",header = F)
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
