#### analyse output from rda.R
library(reshape)
rm(list=ls())
cand <- load('output/candidateLoci.Rda')
cand.venn <- read.delim('output/candidateLoci.summary.txt',sep=" ")
NOTcand <- load('output/NOTcandidateLoci.Rda')
NOTcand.venn <- read.delim('output/NOTcandidateLoci.summary.txt',sep= " ")

### which loci do we see everywhere?
all <- data.frame()
for (i in 1:6)
{
  all <- rbind(all,data.frame(type="cand",pop=names(all.loci)[i],loci=all.loci[[i]]))
  all <- rbind(all,data.frame(type="NOTcand",pop=names(all.loci.NOTcandidates)[i],loci=all.loci.NOTcandidates[[i]]))
  }
all.md <- melt(all,id=c("pop","loci"))
all2 <- cast(all.md,loci~pop)
all2 <- as.data.frame(all2)
all2 <- all2[complete.cases(all2[,-1]),]
### 904 loci have been analyzed across all 6 marshes
all3 <- data.frame()
for(i in 1:dim(all2)[1])
{
  tmp <- as.numeric(ifelse(all2[i,-1]=="cand","1","0"))  ## 1 = candidate loci
  all3 <- rbind(all3,tmp)
  }
all3 <- data.frame(all2$loci,all3)
colnames(all3) <- colnames(all2)

### how many candidate loci are shared by >1 marsh
### only 13 of 904 (1.4%) are shared by 2 marshes
### NOTE: 4 loci shared within south (FL and SC)
### 5 loci shared within north (MA)
### 4 loci shared between north and south. 
print(cand.shared <- dim(all3[rowSums(all3[,-1])>1,])[1])
print(all3[rowSums(all3[,-1])>1,])
### how many candidate loci are unique to 1 marsh
### 122 of 904 (13%)
print(cand.unique <- dim(all3[rowSums(all3[,-1])==1,])[1])
### how many NON candidate loci are there?
### 769/904 (85%)
print(NOTcand.shared <- dim(all3[rowSums(all3[,-1])==0,])[1])
