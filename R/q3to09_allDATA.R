# all pops
# from NGSadmix

pdf('output/q3to09_allDATA.pdf',width=7,height=9)
par(mar=c(.5,2,1.2,3),mfrow=c(6,1),xpd=TRUE)
rm(list=ls())
#library(RColorBrewer)
library(fields)
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
              "k10.qopt")
              #"k09.qopt")
k = 3:8

out.cols <- list(
  c(2,2,3,3,3,3,3,3,3,3,1,1,1,1,1,1,1),
  c(2,2,3,3,3,3,3,3,3,4,1,1,1,1,1,1,1),
  c(5,5,3,3,3,3,3,3,3,4,1,2,2,2,2,2,2),
  c(3,3,4,4,4,4,4,4,4,6,5,2,2,2,2,2,2),
  c(4,4,2,2,2,2,2,2,2,1,3,6,6,7,6,6,6),
  c(7,7,1,1,1,1,1,1,1,2,6,3,3,5,3,3,3),
  c(3,3,4,4,4,4,4,4,4,5,7,8,8,2,8,8,8))

for(i in 1:length(k))
{
#i=3
  dat <- read.delim(paste("data/",filelist[i],sep=""),sep=" ",header = F)[,-(k[i]+1)]
  dat <- dat[order(siteorder),]
  #colors
  #tim12equal = c("#00008F", "#0000EA", "#0047FF", "#00A2FF", "#00FEFF", "#5AFFA5", "#B5FF4A", "#FFED00", "#FF9200", "#FF3700", "#DB0000", "#800000")
  out.cols <- c()
for(ii in 1:length(siteorder.txt))
  {
    tmp.cols <- colSums(dat[pop%in%siteorder.txt[ii],])
    tmp.cols <- (1:ncol(dat))[tmp.cols==max(tmp.cols)]
    out.cols <- c(out.cols,tmp.cols)
  }
  print(out.cols)
  col.to.use <- tim.colors(9)
  col.to.use[out.cols[1]] <- tim.colors(9)[1] ## FLS
  col.to.use[out.cols[6]] <- tim.colors(9)[2] ## TBI
  col.to.use[out.cols[16]] <- tim.colors(9)[3] ## WET
  if(k[i]>=4){col.to.use[out.cols[10]] <- tim.colors(9)[4]} # NCC
  if(k[i]>=5){col.to.use[out.cols[11]] <- tim.colors(9)[5]} # RIS
  if(k[i]>=6){col.to.use[out.cols[15]] <- tim.colors(9)[6]} # WES
  if(k[i]>=7){col.to.use[out.cols[7]] <- tim.colors(9)[7]} # FJS
  if(k[i]>=8){col.to.use[out.cols[14]] <- tim.colors(9)[8]} # SWT
  barplot(t(dat),col=col.to.use,border=NA,ylim=c(-.1,1.1),names.arg = rep("",dim(dat)[1]),space = 0)
  #box()
  segments(0,1.05,0,-0.05)#segments(dim(dat)[1],0,dim(dat)[1],-.05)
  
#tim6equal = 
#brewer.pal(8,"Spectral")

### sites
for(j in 1:length(siteorder.txt))
{
  x <- 1:dim(dat)[1]
  xtmp <- x[pop==siteorder.txt[j]]
  segments(max(xtmp),1,max(xtmp),-.05)
  mtext(side=1,at=mean(xtmp),siteorderNICE.txt[j],cex=.6,line=-1)
}

### states
for(m in 1:length(stateorder.txt))
{
  x <- 1:dim(dat)[1]
  xtmp <- x[state==stateorder.txt[m]]
  segments(max(xtmp),1.1,max(xtmp),0,lwd=2)
  mtext(side=3,at=mean(xtmp),stateorder.txt[m],cex=1,line=-.75)
}
  text(x=320,y=.5,paste("k=",k[i],sep=""),cex=1.5)
}

dev.off()
