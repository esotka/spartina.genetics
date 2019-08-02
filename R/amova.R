### AMOVA
### method 1 = pegas::vegan on mpgl numbers.
library(pegas)
rm(list=ls())
meta <- read.csv('data/Spartina_SNP_SiteID.csv')
gprob<- read.table('data/spartinaNov2017.called.subset.mpgl')
ids <- read.csv('data/inds/allpops_individualIDs_subset.csv')[,2]
pop <- substr(ids,1,3)
reg <- meta$State[match(pop,meta$Site_ID)]
gmat <- gprob[,-1]
colnames(gmat) <- ids
### tall vs short; 6 marshes
st <- c("SWS","WES","RIS","SBI","SFB","FLS",
            "SWT","WET","RIT","TBI","TFB","FLT")
st.reg <- factor(reg[pop%in%st])
st.subpop <- factor(pop[pop%in%st])
st.SiteName <- factor(meta$Site.Name[match(st.subpop,meta$Site_ID)])
st.trt <- factor(meta$Tall.Short[match(st.subpop,meta$Site_ID)])
gmat2 <- gmat[,pop%in%st]
st.ids <- colnames(gmat2)
st.d <- dist(t(gmat2))## distances between rows
print(m1 <- amova(st.d ~ st.SiteName/st.subpop,nperm=1000))
write.pegas.amova(m1,"output/amova.6marshes.csv")
sig2 <- setNames(m1$varcomp$sigma2,rownames(m1$varcomp))
print(getPhi(sig2))
#print(adonis2(st.d ~ st.SiteName/st.subpop,permutations = 100))

###### tall vs short; within SC ###
st <- c("SBI","SFB",
        "TBI","TFB")
st.reg <- factor(reg[pop%in%st])
st.subpop <- factor(pop[pop%in%st])
st.SiteName <- factor(meta$Site.Name[match(st.subpop,meta$Site_ID)])
st.trt <- factor(meta$Tall.Short[match(st.subpop,meta$Site_ID)])
gmat2 <- gmat[,pop%in%st]
st.ids <- colnames(gmat2)
st.d <- dist(t(gmat2))## distances between rows
print(m1 <- amova(st.d ~ st.SiteName/st.subpop,nperm=1000))
write.pegas.amova(m1,"output/amova.SC.csv")
sig2 <- setNames(m1$varcomp$sigma2,rownames(m1$varcomp))
print(getPhi(sig2))

###### tall vs short; within MA ###
st <- c("SWS","WES",
        "SWT","WET")
st.reg <- factor(reg[pop%in%st])
st.subpop <- factor(pop[pop%in%st])
st.SiteName <- factor(meta$Site.Name[match(st.subpop,meta$Site_ID)])
st.trt <- factor(meta$Tall.Short[match(st.subpop,meta$Site_ID)])
gmat2 <- gmat[,pop%in%st]
st.ids <- colnames(gmat2)
st.d <- dist(t(gmat2))## distances between rows
print(m1 <- amova(st.d ~ st.SiteName/st.subpop,nperm=1000))
write.pegas.amova(m1,"output/amova.MA.csv")
sig2 <- setNames(m1$varcomp$sigma2,rownames(m1$varcomp))
print(getPhi(sig2))

##################################################################
### stampp approach using Nei's genetic distance on imputed calls ###
##################################################################
library(StAMPP)
source('R/beagle_calls.fromAllanStrand.R') ### allan's caller. 
loci <- readLines("data/rda/loci.subset2.txt")
Acalls <- beagle_gl(fn="data/309inds.beagle.gz",
          bamlist = "data/inds/allpops_individualIDs_subset.txt",
          rowstart=0,nrows=-1,support=log(2)) 
Acalls.tr <- data.frame(t(Acalls[,-1]))
### tall vs short; 6 marshes
st <- c("SWS","WES","RIS","SBI","SFB","FLS",
        "SWT","WET","RIT","TBI","TFB","FLT")
st.reg <- factor(reg[pop%in%st])
st.subpop <- factor(pop[pop%in%st])
st.SiteName <- factor(meta$Site.Name[match(st.subpop,meta$Site_ID)])
st.ids <- ids[pop%in%st]
Acalls.tr <- Acalls.tr[pop%in%st,]### subset
## impute
Acalls.tr2 <- c()
for(j in 1:ncol(Acalls.tr))
{Acalls.tr2 <- cbind(Acalls.tr2,as.numeric(Acalls.tr[,j]))}
Acalls.tr2 <- data.frame(Acalls.tr2)
Acalls.impute <- apply(Acalls.tr2, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
Acalls.impute <- data.frame(Acalls.impute)
Acalls.impute[Acalls.impute==1] <- "AA"
Acalls.impute[Acalls.impute==2] <- "AB"
Acalls.impute[Acalls.impute==3] <- "BB"
Acalls.tmp <- data.frame(Inds=st.ids,
                        Pop=st.subpop,
                        Ploidy=2,
                        Format="BiA",
                        Acalls.impute)
colnames(Acalls.tmp)[5:dim(Acalls.tmp)[2]] <- loci

Acalls.stampp <- stamppConvert(Acalls.tmp,"r")
### Nei's genetic distance
d <- stamppNeisD(Acalls.stampp, FALSE)
#pop <- factor(Acalls.stampp$Pop)
print(m2 <- pegas::amova(d~st.SiteName/st.subpop,nperm=1000))
write.pegas.amova(m2,"output/amova.6marshes.V2.csv")
#print(adonis2(d ~ st.SiteName/st.subpop,permutations = 100))

###### tall vs short; within SC ###
st <- c("SBI","SFB",
        "TBI","TFB")
st.reg <- factor(reg[pop%in%st])
st.subpop <- factor(pop[pop%in%st])
st.SiteName <- factor(meta$Site.Name[match(st.subpop,meta$Site_ID)])
SC.tmp <- Acalls.tmp[Acalls.tmp$Pop%in%st,]### subset
SC.stampp <- stamppConvert(SC.tmp,"r")
d <- stamppNeisD(SC.stampp, FALSE)
#pop <- factor(Acalls.stampp$Pop)
print(m2 <- pegas::amova(d~st.SiteName/st.subpop,nperm=1000))
write.pegas.amova(m2,"output/amova.SC.V2.csv")

###### tall vs short; within MA ###
st <- c("SWS","WES",
        "SWT","WET")
st.reg <- factor(reg[pop%in%st])
st.subpop <- factor(pop[pop%in%st])
st.SiteName <- factor(meta$Site.Name[match(st.subpop,meta$Site_ID)])
MA.tmp <- Acalls.tmp[Acalls.tmp$Pop%in%st,]### subset
MA.stampp <- stamppConvert(MA.tmp,"r")
d <- stamppNeisD(MA.stampp, FALSE)
#pop <- factor(Acalls.stampp$Pop)
print(m2 <- pegas::amova(d~st.SiteName/st.subpop,nperm=1000))
write.pegas.amova(m2,"output/amova.MA.V2.csv")
