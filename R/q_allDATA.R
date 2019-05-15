# all pops
# from NGSadmix

pdf('output/q_allDATA.pdf',width=7,height=3)
par(mar=c(.5,2,.5,.5))
rm(list=ls())
library(RColorBrewer)
k=6
dat <- read.delim("data/k06run01.qopt",sep=" ",header = F)[,-(k+1)]
ids <- read.csv('data/allpops_individualIDs_subset.csv')[,2]
pop <- substr(ids,1,3)
meta <- read.csv('data/Spartina_SNP_SiteID.csv')
state <- meta$State[match(pop,meta$Site_ID)]
siteorder.txt <- c("FLS","FLT","SFB","TFB","SBI","TBI","FJS","SCT","HWW","NCC",
               "RIS","RIT","SWS","SWT","WES","WET","NHH")
siteorder <- match(pop,siteorder.txt)
pop <- pop[order(siteorder)]
stateorder.txt <- c("FL","SC","NC","RI","MA","NH")
stateorder <- match(state,stateorder.txt)
state <- state[order(siteorder)]
#tim12equal = c("#00008F", "#0000EA", "#0047FF", "#00A2FF", "#00FEFF", "#5AFFA5", "#B5FF4A", "#FFED00", "#FF9200", "#FF3700", "#DB0000", "#800000")
tim6equal = c("#00008F", "#005AFF", "#23FFDC", "#ECFF13", "#FF4A00", "#800000")
#brewer.pal(8,"Spectral")
dat <- dat[order(siteorder),]
barplot(t(dat),col=tim6equal,border=NA,ylim=c(-.1,1.1),names.arg = rep("",dim(dat)[1]),space = 0)
#box()
segments(0,1.05,0,-0.05)#segments(dim(dat)[1],0,dim(dat)[1],-.05)

### sites
for(i in 1:length(siteorder.txt))
{
  x <- 1:dim(dat)[1]
  xtmp <- x[pop==siteorder.txt[i]]
  segments(max(xtmp),1,max(xtmp),-.05)
  mtext(side=1,at=mean(xtmp),siteorder.txt[i],cex=.6,line=-1.5)
}

### states
for(i in 1:length(stateorder.txt))
{
  x <- 1:dim(dat)[1]
  xtmp <- x[state==stateorder.txt[i]]
  segments(max(xtmp),1.05,max(xtmp),0)
  mtext(side=3,at=mean(xtmp),stateorder.txt[i],cex=1,line=-.75)
}

#text(tapply(1:length(siteorder),siteorder,mean),-.05,siteorder.txt,xpd=T,cex=.6)

#popbylat <- pop[order(lat)]
#popbylat.unique <- unique(popbylat)
#popbylat <- factor(popbylat); popbylat <- factor(popbylat,levels=popbylat.unique)
#text(tapply(1:length(popbylat),popbylat,mean),-0.05,unique(popbylat),xpd=T,cex=.6)
dev.off()



#$param <- as.character(q$param)
#q$group <- substr(q$param,nchar(q$param)-4,nchar(q$param))
#q$ind <- rep(ids,k)
#md <- melt(q[,c("group","ind","mean")])
#md2 <- cast(md,ind~group)
# by ind
#md3 <- data.frame(md2)
#md3$pop <- substr(md3$ind,1,3)
#md3$state <- meta$State[match(md3$pop,meta$Site_ID)]
#md3$lat <- meta$lat[match(md3$pop,meta$Site_ID)]
#md3 <- md3[order(md3$lat),]
#rainbow12equal = c("#BF4D4D", "#BF864D", "#BFBF4D", "#86BF4D", "#4DBF4D", "#4DBF86", "#4DBFBF", "#4D86BF", "#4D4DBF", "#864DBF", "#BF4DBF", "#BF4D86")
#rich12equal = c("#000040", "#000093", "#0020E9", "#0076FF", "#00B8C2", "#04E466", "#49FB25", "#E7FD09", "#FEEA02", "#FFC200", "#FF8500", "#FF3300")
