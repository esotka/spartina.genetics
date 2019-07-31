#### delta K
rm(list=ls())
ks <- c("k01","k02","k03","k04","k05","k06","k07","k08","k09","k10")
paths <- paste("data/admix.runs.TvS.POPS.byPOP/",
      c("runs.WESWET","runs.FLSFLT","runs.RISRIT","runs.SBITBI","runs.SFBTFB","runs.SWSSWT"),"/",sep="")
pop <- substr(paths,nchar(paths)-6,nchar(paths)-1)
lnl <- c()
for (j in 1:length(paths))
{
filename <- list.files(path=paths[j],pattern=".log")
for (k in 1:length(ks))
{
  filename.k <- filename[grep(pattern=ks[k],filename)]
  for (i in 1:length(filename.k))
  {
    liketmp <- readLines(paste(paths[j],filename.k[i],sep=""))
    liketmp <- liketmp[grep(pattern = "best like=",liketmp)]
    #print(filename.k[i])#;print(liketmp)
    lnl <- rbind(lnl,data.frame(
      pop=pop[j],
      k=ks[k],
      run=i,
      lnl=substr(liketmp,regexpr(pattern="=",liketmp)[1]+1,regexpr(pattern="=",liketmp)[1]+14)))
    
  }
}
}
lnl$lnl <- as.numeric(as.character(lnl$lnl))
#k <- factor(k)
pdf('output/deltaK.TvS.byPop.pdf',width=5,height=8)
par(mfrow=c(2,1),mar=c(4,4,1,1))
for (j in 1:length(pop))
{
  tmp <- lnl[lnl$pop==pop[j],]
  plot(x=tmp$k,y=tmp$lnl,xaxt="n",xlab="k",ylab="lnl",main=pop[j])
  mtext(at=1:10,1:10,side=1,line=1)
  xbar <- tapply(tmp$lnl,tmp$k,mean)
  std <- tapply(tmp$lnl,tmp$k,sd)
  out <- data.frame(pop=pop[j],xbar,std)
# Evanno et al. 2005 ∆K = m(|L(K + 1) − 2 L(K ) + L(K − 1)|)/s[L(K )]
out$L.prime.k <- c(NA,xbar[-1]-xbar[-length(xbar)])
out$L.dblprime.k <- c(NA,out$L.prime.k[-c(1,length(xbar))]-(out$L.prime.k)[-(1:2)],NA)
out$delta <- out$L.dblprime.k/out$std
print(out)
plot(out$delta,xaxt="n",xlab="k",ylab="delta", type="b")
mtext(at=2:9,2:9,side = 1,line=1)

}
dev.off()