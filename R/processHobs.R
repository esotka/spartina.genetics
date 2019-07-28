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
pdf('output/Hobs~lat.pdf')
print(xyplot(xbar~lat,data=stats,type=c("p","r")))
dev.off()

print(summary(lm(xbar~lat,data=stats)))
