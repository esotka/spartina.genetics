### map of  sites
rm(list=ls())
pdf('output/map of spartina sites.pdf',width=4.5,height=8)
par(mfrow=c(2,1),mar=c(1,1,0,1))
library(maps);library(mapdata);library(RColorBrewer);library(plot3D)

meta <- read.csv('data/Spartina_SNP_SiteID.csv')
gprob<- read.table('data/spartinaNov2017.called.subset.mpgl')
ids <- read.csv('data/inds/allpops_individualIDs_subset.csv')[,2]
pop <- substr(ids,1,3)
reg <- meta$State[match(pop,meta$Site_ID)]

##map of sites
map("state",xlim=c(-87,-68),ylim=c(29,44),col="gainsboro",fill=TRUE)
#state.meta <- data.frame(states=levels(reg),col=brewer.pal(7,"Dark2")[1:6],
 #                        lon=c(-81.70000,-68.50000,-74.13525,-68.50000,-70.38893,-77.95095),
  #                       lat=c(29.89340,42.21079,34.48443,43.10000,40.80000,32.30089))
  
points(meta$long,meta$lat,col="black",pch=21,cex=1.6)
points(meta$long,meta$lat,col=brewer.pal(7,"Dark2")[meta$State],pch=20,cex=2)
text(x=c(-68.5,-68.5,-70.38893,-74.13525,-77.95095,-81.7),
     y=c(43.1,42.21079,40.8,34.48443,32.30089,29.89340),
     c("NH","MA","RI","NC","SC","FL"),adj=1)
box()

### PCA
gmat <- t(gprob[,-1])
pc.cr <- prcomp(gmat)
col.reg <- brewer.pal(6,"Dark2")[reg]
plot(pc.cr$x[,1],pc.cr$x[,2],col=col.reg,type="n",xlab="",ylab="",xlim=c(-15,18))
#mtext(side=c(1,2),line=2,c("PCA1","PCA2"))
text(pc.cr$x[,1],pc.cr$x[,2],col=col.reg,pop,cex=.5)
text(x=c(-6.147294,-11.427858,-3.110970,8,15),
     y=c(16.592080,-1.396376,-4.144612,-2.80,2.601059),
     c("FL","SC","NC","RI","NH & MA"),cex=.8)

print(summary(eigenvals(pc.cr))[,1:5])
dev.off()
