## Mean Fst - several categories of distance

library(reshape)
library(ggplot2)
rm(list=ls())
nloci = 2735

### short vs tall comparisons

filelist.tvs <- c("data/fst/FLS.FLT.pops.fst",
                  "data/fst/SBI.TBI.pops.fst",
                  "data/fst/SFB.TFB.pops.fst",
                  "data/fst/RIS.RIT.pops.fst",
                  "data/fst/SWS.SWT.pops.fst",
                  "data/fst/WES.WET.pops.fst"
)
fst.all <- c()
for (i in 1:length(filelist.tvs))
{
  tmp <- read.table(filelist.tvs[i])
  tmp <- tmp[1:2735,]
  tmp[is.na(tmp[,4]),4] <- 0
  tmp[tmp[,4]<0,4] <- 0 ### set all FST <0 to 0
  tmp[tmp[,4]>1,4] <- 1 ### convert FST > 1 to 1
  fst.tmp <- tmp[,4]
  fst.tmp.xbar <- median(fst.tmp,na.rm=T)
  fst.all <- rbind(fst.all,data.frame(sites=filelist.tvs[i],fst=fst.tmp,type = "TvS"))
}
### sites <5km apart
filelist.5 <- c(
  "data/fst/fst.combineTS/HWW.SBITBI.pops.fst",
  "data/fst/HWW.SCT.pops.fst",
  "data/fst/fst.combineTS/HWW.SFBTFB.pops.fst",
  "data/fst/fst.combineTS/SCT.SBITBI.pops.fst",
  "data/fst/fst.combineTS/SCT.SFBTFB.pops.fst",
  "data/fst/fst.combineTS/SFBTFB.SBITBI.pops.fst",
  "data/fst/fst.combineTS/WESWET.SWSSWT.pops.fst")

for (i in 1:length(filelist.5))
{
  tmp <- read.table(filelist.5[i])
  tmp <- tmp[1:2735,]
  tmp[is.na(tmp[,4]),4] <- 0
  tmp[tmp[,4]<0,4] <- 0 ### set all FST <0 to 0
  tmp[tmp[,4]>1,4] <- 1 ### convert FST > 1 to 1
  fst.tmp <- tmp[,4]
  fst.tmp.xbar <- median(fst.tmp,na.rm=T)
  fst.all <- rbind(fst.all,data.frame(sites=filelist.5[i],fst=fst.tmp,type = "marshes.5"))
}

filelist.15 <- c(
  "data/fst/FJS.HWW.pops.fst",
  "data/fst/fst.combineTS/FJS.SBITBI.pops.fst",
  "data/fst/FJS.SCT.pops.fst",
  "data/fst/fst.combineTS/FJS.SFBTFB.pops.fst",
  "data/fst/fst.combineTS/NHH.SWSSWT.pops.fst",
  "data/fst/fst.combineTS/NHH.WESWET.pops.fst")

for (i in 1:length(filelist.15))
{
  tmp <- read.table(filelist.15[i])
  tmp <- tmp[1:2735,]
  tmp[is.na(tmp[,4]),4] <- 0
  tmp[tmp[,4]<0,4] <- 0 ### set all FST <0 to 0
  tmp[tmp[,4]>1,4] <- 1 ### convert FST > 1 to 1
  fst.tmp <- tmp[,4]
  fst.tmp.xbar <- median(fst.tmp,na.rm=T)
  fst.all <- rbind(fst.all,data.frame(sites=filelist.15[i],fst=fst.tmp,type = "marshes.15"))
}



fst.all$sites <- factor(fst.all$sites)
xbar <- tapply(fst.all$fst,fst.all$sites,mean,na.rm=T)
std <- tapply(fst.all$fst,fst.all$sites,sd,na.rm=T)
sterr <- std/dim(tmp)[1]
xbar.out <- data.frame(xbar,std,sterr,type=fst.all$type[match(names(xbar),fst.all$sites)])
xbar.out <- xbar.out[complete.cases(xbar.out),]
m <- lm(xbar~type,data=xbar.out)
print(anova(m))

pdf('output/Fst.morePlots.pdf',width=5,height=4)

xbar.all <- tapply(xbar.out$xbar,xbar.out$type,mean)
std.all <- tapply(xbar.out$xbar,xbar.out$type,sd)
plot(x=1:3,xbar.all,xlim=c(0.5,3.5),ylim=c(0,.3),col="black",
     pch=20,cex=3,xaxt="n",xlab="",ylab="Fst")
arrows(x0=1:3,y0=xbar.all+std.all,x1=1:3,y1=xbar.all-std.all,col=alpha("black",.5),lwd=2,angle = 90,code = 3,length = .1)
x=jitter((1:3)[xbar.out$type],1)
points(y=xbar.out$xbar,x,col="black",cex=.7)
#segments(x,xbar.out$xbar+xbar.out$sterr,x,xbar.out$xbar-xbar.out$sterr)
#c("black","red")[xbar.out$type]
mtext(c("Short vs Tall","Marshes","Marshes"),at=1:3,side=1,line=1)
mtext(c("<200 m","0.5-5 km","5-15 km"),at=1:3,side=1,line=2)
dev.off()
