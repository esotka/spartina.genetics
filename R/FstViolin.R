### make violin density plots of pairwise FST

library(reshape)
library(ggplot2)
rm(list=ls())
nloci = 2735

pdf('output/FstViolin.pdf',width=5,height=4)
### short vs tall comparisons

filelist.tvs <- c("data/SvT/FLS.FLT.pops.fst",
              "data/SvT/SBI.TBI.pops.fst",
              "data/SvT/SFB.TFB.pops.fst",
              "data/SvT/RIS.RIT.pops.fst",
              "data/SvT/SWS.SWT.pops.fst",
              "data/SvT/WES.WET.pops.fst"
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
filelist <- c(
  "data/SvT/HWW.SBITBI.pops.fst",
  "data/SvT/HWW.SCT.pops.fst",
  "data/SvT/HWW.SFBTFB.pops.fst",
  "data/SvT/SCT.SBITBI.pops.fst",
  "data/SvT/SCT.SFBTFB.pops.fst",
  "data/SvT/SFBTFB.SBITBI.pops.fst",
  "data/SvT/WESWET.SWSSWT.pops.fst")


for (i in 1:length(filelist))
{
  tmp <- read.table(filelist[i])
  tmp <- tmp[1:2735,]
  tmp[is.na(tmp[,4]),4] <- 0
  tmp[tmp[,4]<0,4] <- 0 ### set all FST <0 to 0
  tmp[tmp[,4]>1,4] <- 1 ### convert FST > 1 to 1
  fst.tmp <- tmp[,4]
  fst.tmp.xbar <- median(fst.tmp,na.rm=T)
  fst.all <- rbind(fst.all,data.frame(sites=filelist[i],fst=fst.tmp,type = "marshes"))
  }

fst.all$sites <- factor(fst.all$sites)
xbar <- tapply(fst.all$fst,fst.all$sites,mean,na.rm=T)
std <- tapply(fst.all$fst,fst.all$sites,sd,na.rm=T)
sterr <- std/dim(tmp)[1]
xbar.out <- data.frame(xbar,std,sterr,type=fst.all$type[match(names(xbar),fst.all$sites)])
xbar.out <- xbar.out[complete.cases(xbar.out),]
m <- lm(xbar~type,data=xbar.out)
print(anova(m))
fst.all$locusnum <- rep(1:2735,length(unique(fst.all$sites)))
md <- melt(fst.all,id.vars = c("sites","locusnum","type"),measure.vars = "fst")
sites.nice <-  c("FLS vs FLT","SCBS vs SCBT","SCFS vs SCFT","RIS vs RIT","MASS vs MAST",
                 "MAWS vs MAWT","SC2 vs SCB","SC2 vs SC3",  "SC2 vs SCF","SC3 vs SCB","SC3 vs SCF","SCF vs SCB","MAW vs MAS")
md$sites2 <- sites.nice[match(md$sites,unique(md$sites))]
md$sites2 <- factor(md$sites2,levels=sites.nice)
p <- ggplot(md, aes(x=sites2, y=value)) +
  geom_violin(scale="width",aes(fill=type),show.legend = F) + 
  scale_fill_manual(values = c(TvS = "grey", marshes = "black")) +
  labs(x="", y="Fst",cex=3) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12,angle = 90, vjust = .5, hjust = 1),
    axis.text.y = element_text(size=12)
  )

print(p)
dev.off()

pdf('output/FstMeans.pdf',width=4,height=4)

xbar.all <- tapply(xbar.out$xbar,xbar.out$type,mean)
std.all <- tapply(xbar.out$xbar,xbar.out$type,sd)
plot(x=1:2,xbar.all,xlim=c(0.5,2.5),ylim=c(0,.3),col="black",
     pch=20,cex=3,xaxt="n",xlab="",ylab="Fst")
arrows(x0=1:2,y0=xbar.all+std.all,x1=1:2,y1=xbar.all-std.all,col=alpha("black",.5),lwd=2,angle = 90,code = 3,length = .1)
x=jitter((1:2)[xbar.out$type],1)
points(y=xbar.out$xbar,x,col="black",cex=.7)
#segments(x,xbar.out$xbar+xbar.out$sterr,x,xbar.out$xbar-xbar.out$sterr)
#c("black","red")[xbar.out$type]
mtext(c("Short vs Tall","Marshes"),at=1:2,side=1,line=1)
mtext(c("0.02-0.16 km","0.65-4.52 km"),at=1:2,side=1,line=2)
dev.off()
