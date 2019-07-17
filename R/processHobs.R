# process heterozygosity estimate from realSFS and angsd
# http://www.popgen.dk/angsd/index.php/Heterozygosity

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
print(stats)
write.csv(stats,"output/Hobs.csv")
