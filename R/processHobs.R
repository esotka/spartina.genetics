# process heterozygosity estimate from realSFS and angsd
# http://www.popgen.dk/angsd/index.php/Heterozygosity

rm(list=ls())
meta <- read.csv('data/Spartina_SNP_SiteID.csv')
filelist <- list.files('data/Hobs/') # est.ml files. "For diploid single samples the hetereo zygosity is simply second value in the SFS/AFS" 
out <- c()
for (i in 1:length(filelist))
{
  a <- scan(paste('data/Hobs/',filelist[i],sep=""))
  out <- rbind(out,data.frame(file=filelist[i],het=a[2]/sum(a)))
}
out$site <- substr(out$file,1,3)
stats <- data.frame(xbar=tapply(out$het,out$site,mean,na.rm=T),
                 sd=tapply(out$het,out$site,sd,na.rm=T),
                 n=as.numeric(table(out$site)))
stats$se <- stats$sd/sqrt(stats$n)
siteorder.txt <-     c("FLS","FLT","SFB","TFB","SBI","TBI","FJS","HWW","SCT","NCC",
                       "RIS","RIT","SWS","SWT","WES","WET","NHH")
siteorderNICE.txt <- c("FLS","FLT","SCFS","SCFT","SCBS","SCBT","SCFJ","SC2","SC3","NC1",
                       "RIS","RIT","MASS","MAST","MAWS","MAWT","NH1")

stats$sitenamesNICE <- siteorderNICE.txt[match(rownames(stats),siteorder.txt)]
print(stats)
write.csv(stats,"output/Hobs.csv")

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
