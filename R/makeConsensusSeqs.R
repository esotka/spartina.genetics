### make consensus seqs
rm(list=ls())
library(snpStats)
ind <- as.character(read.delim('data/indIDS.T.S.POPS.txt',header=F)[,1])
#snps <- read.delim('data/')
dat <- read.beagle('data/t.s.pops.beagle.gz',header=T,nsnp = 2735,rownames = ind)
out <- c()
for (i in 1:length(ind))
{
  tmp <- as(dat[i,],'character')
  out <- rbind(out,tmp)
}
### convert to As and Ts
out[out=="Uncertain"] <- "N"
out[out=="A/A"] <- "A"
out[out=="A/B"] <- "W"
out[out=="B/B"] <- "T"
print(table(out))
rownames(out) <- paste(">",ind,"\r",sep="")
write.table(out,'output/indIDS.T.S.POPS.fasta',quote=F,row.names = T,col.names=F,sep="")
### make NJ tree
### you need to change the coding to "Windows" within BBEdit
meta <- read.csv("data/Spartina_SNP_SiteID.csv")
dna <- read.dna('output/indIDS.T.S.POPS.fasta',format="fasta")
t.s <- meta$Tall.Short[match(substr(rownames(dna),1,3),meta$Site_ID)]
d <- dist.gene(dna,method="percentage")
pdf('output/indIDS.T.S.POPS.NJ.pdf')
plot(nj(d),type="unrooted",show.tip.label=F)
tiplabels(substr(rownames(dna),1,3),frame="none",cex=.5,col=c("black","red")[t.s])#,show.tip.label=F)
dev.off()


