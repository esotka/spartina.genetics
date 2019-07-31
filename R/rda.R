#### following vegan::RDA from https://popgen.nescent.org/2018-03-27_RDA_GEA.html
#### as suggested by Forester, B. R., J. R. Lasky, H. H. Wagner, and D. L. Urban. 2018. Comparing methods for detecting multilocus adaptation with multivariate genotype–environment associations. Molecular Ecology 27:2215–2233.

### list of candidate loci ###
rm(list=ls())
all.loci <- list()
all.loci.NOTcandidates <- list()
pdf("output/rda.pdf")
library(psych)
library(vegan)
source('R/beagle_calls.fromAllanStrand.R') ### allan's caller. 
loci <- readLines("data/rda/loci.subset2.txt")
Aind <- paste("data/rda/",c("FLT.txt",'TBI.txt',"TFB.txt","RIT.txt","WET.txt","SWT.txt"),sep="")
Afn <- paste("data/rda/",c("FLT.beagle.gz","TBI.beagle.gz",'TFB.beagle.gz','RIT.beagle.gz','WET.beagle.gz',"SWT.beagle.gz"),sep="")
Bind <- paste("data/rda/",c('FLS.txt',"SBI.txt",'SFB.txt',"RIS.txt","WES.txt","SWS.txt"),sep="")
Bfn <- paste("data/rda/",c("FLS.beagle.gz",'SBI.beagle.gz',"SFB.beagle.gz",'RIS.beagle.gz','WES.beagle.gz',"SWS.beagle.gz"),sep="")

for(i in 1:6)
{

### tall datasets
beagle_gl(fn=Afn[i],bamlist = Aind[i],rowstart=0,nrows=-1,support=log(2)) -> Acalls
Acalls.tr <- data.frame(t(Acalls[,-1]))
### short datasets
beagle_gl(fn=Bfn[i],bamlist=Bind[i],rowstart=0,nrows=-1,support=log(2)) -> Bcalls
Bcalls.tr <- data.frame(t(Bcalls[,-1]))

tall.short.inds <- c(readLines(Aind[i]),readLines(Bind[i]))
tall.short <- rbind(Acalls.tr,Bcalls.tr)
tall.short2 <- c()
for(j in 1:ncol(tall.short))
{tall.short2 <- cbind(tall.short2,as.numeric(tall.short[,j]))}
tall.short2 <- data.frame(tall.short2)
rownames(tall.short2) <- c(readLines(Aind[i]),readLines(Bind[i]))
colnames(tall.short2) <- loci
### keep SNPs in which >=50% of individuals have calls
tall.short3 <- tall.short2[,colSums(is.na(tall.short2)) < length(tall.short.inds)/2]

### impute
tall.short.imp <- apply(tall.short3, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
tall.short.imp <- data.frame(tall.short.imp)

############# ENV DATA ##############
env <- data.frame(ind=c(readLines(Aind[i]),readLines(Bind[i])),
                  TvS= c(rep("Tall",length(readLines(Aind[i]))),rep("Short",length(readLines(Bind[i])))))
env$ind <- as.character(env$ind)
##### RUN RDA #######
TvS.rda <- vegan::rda(tall.short.imp~env$TvS,scale=T)
print(RsquareAdj(TvS.rda)) ### print proportion of variance explained
print(summary(eigenvals(TvS.rda, model = "constrained")))
print(signif.full <- anova.cca(TvS.rda, parallel=getOption("mc.cores"))) # default is permutation=999.  The null hypothesis is that no linear relationship exists between the SNP data and the environmental predictors. See ?anova.cca for more details and options

plot(TvS.rda, type="n", scaling=3,main = c(Aind[i],Bind[i]))
points(TvS.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(TvS.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=c("red","black")[env$TvS]) # the ecotypes
text(TvS.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(env$TvS), bty="n", col="gray32", pch=21, cex=1, pt.bg=c("red","black"))

#### candidate loci ####
load.rda <- scores(TvS.rda, choices=1, display="species")  # Species scores for the first constrained axis
#hist(load.rda)

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

#cand1 <- outliers(load.rda[,1],2) # 2.5 SD == outlier
cand.threshold <- sort(abs(load.rda),decreasing = T)[1:20]# top 20 loci
cand1 <- rownames(load.rda)[abs(load.rda)>=min(cand.threshold)]
#all.loci[[i]] <- names(cand1)
#all.loci.NOTcandidates[[i]] <- colnames(tall.short.imp)[!colnames(tall.short.imp)%in%names(cand1)]
all.loci[[i]] <- cand1
all.loci.NOTcandidates[[i]] <- colnames(tall.short.imp)[!colnames(tall.short.imp)%in%cand1]
}
dev.off()
names(all.loci) <- c("FL",'SC-BI',"SC-FB","RI","MA-WE","MA-SW")
names(all.loci.NOTcandidates) <- c("FL",'SC-BI',"SC-FB","RI","MA-WE","MA-SW")
save(all.loci,file = "output/candidateLoci.Rda")
save(all.loci.NOTcandidates,file = "output/NOTcandidateLoci.Rda")
library(gplots)
out <- venn(all.loci,show.plot=F,intersections=F)
write.table(out,"output/candidateLoci.summary.txt",quote=F,row.names=F)
outNO <- venn(all.loci.NOTcandidates,show.plot=F,intersections=F)
write.table(outNO,"output/NOTcandidateLoci.summary.txt",quote=F,row.names=F)


### other
# library(snpStat)
#ind <- as.character(read.delim('RIT.txt',header=F)[,1])
#fn <- read.beagle('RIT.beagle.gz',nsnp=2735,rownames=ind)
