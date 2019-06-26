## FST generated from ngsFST (see https://github.com/mfumagalli/ngsPopGen)

### convert FST > 1 to 1. Set all FST < 0 equal to zero

### FST by geographic distance (across all populations)

rm(list=ls())
nloci = 2735
filename <- read.table('data/fst/filelist') ## this is all of the pairwise comparison files in this folder
fst.all <- c()
#Columns are ordered as: A, AB, f, FST, Pvar; where A is the expectation of genetic variance between populations, AB is the expectation of the total genetic variance, f is the correcting factor for the ratio of expectations, FST is the per-site FST value, Pvar is the probability for the site of being variable.
for (i in 1:dim(filename)[1])
{tmp <- read.table(paste('data/fst/',as.character(filename[i,1]),sep=""))
tmp <- tmp[1:2735,]
tmp[is.na(tmp[,4]),4] <- 0
tmp[tmp[,4]<0,4] <- 0 ### set all FST <0 to 0
tmp[tmp[,4]>1,4] <- 1 ### convert FST > 1 to 1
#tmp <- tmp[!tmp[,4]>1,] ### remove FST >1
#num.loci <- dim(tmp)[1]
fst.tmp <- tmp[,4]
fst.tmp.xbar <- mean(fst.tmp,na.rm=T)
fst.all <- rbind(fst.all,data.frame(sites=filename[i,1],fst=fst.tmp.xbar))#,num.loci))
#plot(density(fst.tmp),xlim=c(-.1,1),main=filename[i,1])
}
fst.all$pops1 <- substr(fst.all$sites,1,3)
fst.all$pops2 <- substr(fst.all$sites,5,7)
#dev.off()
### get geographic distances
library(geosphere)
meta <- read.csv('data/Spartina_SNP_SiteID.csv')
km <- c()
for (j in 1:dim(fst.all)[1])
{
  km <- c(km,distHaversine(c(meta[match(fst.all$pops1[j],meta$Site_ID),"lat"],meta[match(fst.all$pops1[j],meta$Site_ID),"long"]),c(meta[match(fst.all$pops2[j],meta$Site_ID),"lat"],meta[match(fst.all$pops2[j],meta$Site_ID),"long"]))/1000)
}
fst.all$km <- km
hist(fst.all$fst,breaks=20,col="grey")
fit <- lm(fst~km,data=fst.all)
summary(fit)
visreg(fit)
print(fst.all)
write.csv(fst.all,"output/fst.siteXpop.csv")
