#### delta K
rm(list=ls())
ks <- c("k2r","k3r","k4r","k5r","k6r","k7r","k8r","k9r","k10",
        "k11","k12","k13","k14","k15",
        "k16","k17","k18","k19","k20")

#ks <- c("k01","k02","k03","k04","k05","k06","k07","k08","k09")#,"k10")
filename <- list.files(path="data/admix.runs.TvS.POPS/",pattern=".log")
lnl <- c()
for (k in 1:length(ks))
{
  filename.k <- filename[grep(pattern=ks[k],filename)]
  for (i in 1:length(filename.k))
  {
    liketmp <- readLines(paste("data/admix.runs.TvS.POPS/",filename.k[i],sep=""))
    liketmp <- liketmp[grep(pattern = "best like=",liketmp)]
    #print(filename.k[i])#;print(liketmp)
    lnl <- rbind(lnl,
                 data.frame(k=ks[k],run=i,lnl=substr(liketmp,regexpr(pattern="=",liketmp)[1]+1,regexpr(pattern="=",liketmp)[1]+14)))
    
  }
}
lnl$lnl <- as.numeric(as.character(lnl$lnl))
#k <- factor(k)
pdf('output/deltaK.TvS.only.pdf',width=5,height=8)
par(mfrow=c(2,1),mar=c(4,4,1,1))
plot(x=lnl$k,y=lnl$lnl,xaxt="n",xlab="k",ylab="lnl")
mtext(at=1:19,2:20,side=1,line=1)
text(1,max(lnl$lnl)*1.01,"B",cex=2)
#print(y)
xbar <- tapply(lnl$lnl,lnl$k,mean)
std <- tapply(lnl$lnl,lnl$k,sd)
out <- data.frame(xbar,std)
# Evanno et al. 2005 ∆K = m(|L(K + 1) − 2 L(K ) + L(K − 1)|)/s[L(K )]
out$L.prime.k <- c(NA,xbar[-1]-xbar[-length(xbar)])
out$L.dblprime.k <- c(NA,out$L.prime.k[-c(1,length(xbar))]-(out$L.prime.k)[-(1:2)],NA)
out$delta <- out$L.dblprime.k/out$std
print(out)
plot(out$delta,xaxt="n",xlab="k",ylab="delta", type="b")
mtext(at=2:18,3:19,side = 1,line=1)
text(2,max(out$delta,na.rm=T)*.9,"C",cex=2)
#plot(out$L.dblprime.k,xaxt="n",xlab="k",type="b");mtext(at=1:length(xbar),rownames(out),side = 1)

dev.off()