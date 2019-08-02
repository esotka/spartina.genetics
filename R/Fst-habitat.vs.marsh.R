#### habitat Fst

rm(list=ls())
library(lattice)
library(ggplot2)
fst.all <- read.csv('output/fst.siteXpop.csv')
meta <- read.csv("data/Spartina_SNP_SiteID.csv")

### habitat
st <- c("SWS","WES","RIS","SBI","SFB","FLS",
        "SWT","WET","RIT","TBI","TFB","FLT")

fst.tvs <- fst.all[fst.all$sites%in%c("FLS.FLT.pops.fst",
                  "SBI.TBI.pops.fst",
                  "SFB.TFB.pops.fst",
                  "RIS.RIT.pops.fst",
                  "SWS.SWT.pops.fst",
                  "WES.WET.pops.fst"),]

fst.5km <- fst.all[fst.all$sites%in%
                     c("HWW.SBI.pops.fst",
                     "HWW.TBI.pops.fst",
                       "HWW.SCT.pops.fst",
                       "HWW.SFB.pops.fst",
                      "HWW.TFB.pops.fst",
                       "SBI.SCT.pops.fst",
                     "SCT.TBI.pops.fst",
                       "SCT.SFB.pops.fst",
                     "SCT.TFB.pops.fst",
                       "SBI.SFB.pops.fst",
                       #"SFB.TBI.pops.fst",
                       #"SBI.TFB.pops.fst",
                       "TBI.TFB.pops.fst",
                       "SWS.WES.pops.fst",
                       #"SWT.WES.pops.fst",
                       #"SWS.WET.pops.fst",
                       "SWT.WET.pops.fst"),]

fst.15km  <- fst.all[fst.all$sites%in%
          c("FJS.HWW.pops.fst",
            "FJS.SBI.pops.fst",
            "FJS.TBI.pops.fst",
            "FJS.SCT.pops.fst",
            "FJS.SFBTFB.pops.fst",
            "FJS.TFB.pops.fst",
            "NHH.SWS.pops.fst",
            "NHH.SWT.pops.fst",
            "NHH.WET.pops.fst",
            "NHH.WES.pops.fst"),]



out <- rbind(fst.tvs[,2:3],fst.5km[,2:3],fst.15km[,2:3])
out$type <- factor(c(rep("habitat",dim(fst.tvs)[1]),
              rep("5km",dim(fst.5km)[1]),
              rep("15km",dim(fst.15km)[1])))
out$type <- factor(out$type,levels=levels(out$type)[c(3,2,1)])
print(bwplot(fst~type,data=out,type="p"))
m <- lm(fst~type,data=out)
print(anova(m))  
print(TukeyHSD(aov(m)))

xbar.all <- tapply(out$fst,out$type,mean)
std.all <- tapply(out$fst,out$type,sd)
n.all <- tapply(out$fst,out$type,length)
se.all <- std.all/sqrt(n.all)

pdf('output/Fst-habitat.vs.marsh.pdf',width=5,height=4)
plot(x=1:3,xbar.all,xlim=c(0.5,3.5),ylim=c(0,.3),col="black",
     pch=20,cex=3,xaxt="n",xlab="",ylab="Fst")
arrows(x0=1:3,y0=xbar.all+se.all,x1=1:3,y1=xbar.all-se.all,col=alpha("black",.5),lwd=2,angle = 90,code = 3,length = .1)
x=jitter((1:3)[out$type],1)
points(y=out$fst,x,col="black",cex=.7)
#segments(x,xbar.out$xbar+xbar.out$sterr,x,xbar.out$xbar-xbar.out$sterr)
#c("black","red")[xbar.out$type]
mtext(c("Short vs Tall","Marshes","Marshes"),at=1:3,side=1,line=1)
mtext(c("<200 m","0.5-5 km","5-15 km"),at=1:3,side=1,line=2)
dev.off()
