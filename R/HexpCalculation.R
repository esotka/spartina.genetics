### calculate Hexp
library(reshape2)
library(lattice)
library(scales)
rm(list=ls())
meta <- read.csv('data/Spartina_SNP_SiteID.csv')
files <- paste(list.files(path = 'data/allelefreq/',pattern = "."),sep="")
pop.to.use <- substr(files,1,3)

af.all <- c()
for (j in 1:length(files))
{
  tmp <- read.delim(paste('data/allelefreq/',files[j],sep=""),header=T)
  af.all <- rbind(af.all,data.frame(pop=pop.to.use[j],tmp))
  }
af.all$chr_pos <- paste(af.all$chromo,"_",af.all$position,sep="")
md <- melt(af.all,id=c("chr_pos","pop"),measure.vars = "PPmaf")
md$hexp <- 2*(md$value*(1-md$value)) #2pq

stats <- data.frame(xbar=tapply(md$hexp,md$pop,mean,na.rm=T),
                    sd=tapply(md$hexp,md$pop,sd,na.rm=T),
                    n.loci=as.numeric(table(md$pop)))
stats$se <- stats$sd/sqrt(stats$n)
print(stats)
write.csv(stats,"output/Hexp.csv")

