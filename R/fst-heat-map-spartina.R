### fst heat map

library(lattice);library(visreg);library(RColorBrewer)
rm(list=ls())
fst.all <- read.csv('data/fst.siteXpop.csv')
meta <- read.csv("data/Spartina_SNP_SiteID.csv")

### make site order matrix
#siteorder <- data.frame(pop = unique(c(as.character(fst.all$pops1),as.character(fst.all$pops2))))
#
siteorder <-data.frame(pop= c("FLS","FLT","SFB","TFB","SBI","TBI","FJS","HWW","SCT","NCC","RIS","RIT","SWS","SWT","WES","WET","NHH"),
pop.nice=c("FLS","FLT","SCFS","SCFT","SCBS","SCBT","SCFJ","SC2","SC3","NC1","RIS","RIT","MASS","MAST","MAWS","MAWT","NH1"))
siteorder$state = meta$State[match(siteorder$pop,meta$Site_ID)]
siteorder$pop <- as.character(siteorder$pop)

### make full distance matrix
full <- c()
for (i in 1:length(siteorder$pop))
{
  tmp <- fst.all[fst.all$pops1==siteorder$pop[i] | fst.all$pops2==siteorder$pop[i],c("pops1","pops2","fst")]
  tmp <- rbind(tmp, data.frame(pops1=siteorder$pop[i],pops2=siteorder$pop[i],fst=0))
  tmp$othersite <- ifelse(tmp$pops1==siteorder$pop[i],as.character(tmp$pops2),as.character(tmp$pops1))
  tmp <- tmp[order(match(tmp$othersite,siteorder$pop)),]
  full <- cbind(full,tmp$fst)
}
rownames(full) <- siteorder$pop.nice
colnames(full) <- siteorder$pop.nice
#col.l <- colorRampPalette(c('white',"lightgrey","darkgrey",'black'))(100)
col.l <- colorRampPalette(brewer.pal(11, "RdBu"))
f <- levelplot(full,cex=.3,col.regions=col.l); print(f)
pdf('output/fst-heat-map-spartina.pdf',width=11,height=8)
print(f)
dev.off()

#> head(full)
#hik       shk       waj       hay       cnt       sou       mng
#hik 0.0000000 0.7916595 0.7025675 0.7689090 0.7337842 0.6674295 0.6752981
#shk 0.7916595 0.0000000 0.4443882 0.6267603 0.7232151 0.5934977 0.5998149
#waj 0.7025675 0.4443882 0.0000000 0.4669381 0.6131833 0.4654769 0.4732733

