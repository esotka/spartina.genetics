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
siteorder.txt <-     c("FLS","FLT","SFB","TFB","SBI","TBI","FJS","HWW","SCT","NCC",
                       "RIS","RIT","SWS","SWT","WES","WET","NHH")
siteorderNICE.txt <- c("FLS","FLT","SCFS","SCFT","SCBS","SCBT","SCFJ","SC2","SC3","NC1",
                       "RIS","RIT","MASS","MAST","MAWS","MAWT","NH1")

stats$sitenamesNICE <- siteorderNICE.txt[match(rownames(stats),siteorder.txt)]
print(stats)
write.csv(stats,"output/Hexp.csv")

stats$lat <- meta$lat[match(rownames(stats),meta$Site_ID)]
stats$NS <- meta$State[match(rownames(stats),meta$Site_ID)]
stats$NS <- factor(ifelse(stats$NS%in%c("FL","NC","SC"),"S","N"))
stats.noShort <- stats[!stats$sitenamesNICE%in%c("FLS","RIS","SCFS","MASS","MAWS","SCBS"),]
stats.noTall <- stats[!stats$sitenamesNICE%in%c("FLT","RIT","SCFT","MAST","MAWT","SCBT"),]

out <- c()
### with latitude 
# all
out <- rbind(out,data.frame(test="lat",data="all",p=summary(lm(xbar~lat,data=stats))$coefficients[8]))
# NO FL
out <- rbind(out,data.frame(test="lat",data="no FL",p=summary(lm(xbar~lat,data=stats[!rownames(stats)%in%c("FLS","FLT"),]))$coefficients[8]))
# No short all
out <- rbind(out,data.frame(test="lat",data="no short-all",p=summary(lm(xbar~lat,data=stats.noShort))$coefficients[8]))
# No short, no FL
out <- rbind(out,data.frame(test="lat",data="no short-noFL",p=summary(lm(xbar~lat,data=stats.noShort[!rownames(stats.noShort)%in%c("FLS","FLT"),]))$coefficients[8]))
# no tall, all
out <- rbind(out,data.frame(test="lat",data="no tall-all",p=summary(lm(xbar~lat,data=stats.noTall))$coefficients[8]))
# no tall, no FL
out <- rbind(out,data.frame(test="lat",data="no tall-noFL",p=summary(lm(xbar~lat,data=stats.noTall[!rownames(stats.noTall)%in%c("FLS","FLT"),]))$coefficients[8]))

### with region
# all
out <- rbind(out,data.frame(test="NS",data="all",p=summary(lm(xbar~NS,data=stats))$coefficients[8]))
# NO FL
out <- rbind(out,data.frame(test="NS",data="no FL",p=summary(lm(xbar~NS,data=stats[!rownames(stats)%in%c("FLS","FLT"),]))$coefficients[8]))
# No short all
out <- rbind(out,data.frame(test="NS",data="no short-all",p=summary(lm(xbar~NS,data=stats.noShort))$coefficients[8]))
# No short, no FL
out <- rbind(out,data.frame(test="NS",data="no short-noFL",p=summary(lm(xbar~NS,data=stats.noShort[!rownames(stats.noShort)%in%c("FLS","FLT"),]))$coefficients[8]))
# no tall, all
out <- rbind(out,data.frame(test="NS",data="no tall-all",p=summary(lm(xbar~NS,data=stats.noTall))$coefficients[8]))
# no tall, no FL
out <- rbind(out,data.frame(test="NS",data="no tall-noFL",p=summary(lm(xbar~NS,data=stats.noTall[!rownames(stats.noTall)%in%c("FLS","FLT"),]))$coefficients[8]))

### tall vs short

ts <- stats[stats$sitenamesNICE%in%c("FLS","RIS","SCFS","MASS","MAWS","SCBS","FLT","RIT","SCFT","MAST","MAWT","SCBT"),]
ts$site <- factor(meta$Site.Name[match(rownames(ts),meta$Site_ID)])
ts$ts <- substr(ts$sitenamesNICE,nchar(ts$sitenamesNICE),nchar(ts$sitenamesNICE))
library(reshape2)
ts2 <- ts[,c("site","ts","xbar")]
md <- cast(ts2,site~ts,value="xbar")
out <- rbind(out,data.frame(test="TvsS-paired",data="6 marshes",
                            p=t.test(md$S,md$T,paired=T)$p.value))

out$p <- round(out$p,3)
out$sigif <- ifelse(out$p<0.05,"*","")
print(out)

