### AMOVA
### tall vs short; 6 marshes

rm(list=ls())
meta <- read.csv('data/Spartina_SNP_SiteID.csv')
gprob<- read.table('data/spartinaNov2017.called.subset.mpgl')
ids <- read.csv('data/inds/allpops_individualIDs_subset.csv')[,2]
pop <- substr(ids,1,3)
reg <- meta$State[match(pop,meta$Site_ID)]

gmat <- gprob[,-1]
colnames(gmat) <- ids


st <- c("SWS","WES","RIS","SBI","SFB","FLS",
            "SWT","WET","RIT","TBI","TFB","FLT")
st.reg <- factor(reg[pop%in%st])
st.pop <- factor(pop[pop%in%st])
st.trt <- factor(meta$Tall.Short[match(st.pop,meta$Site_ID)])

gmat2 <- gmat[,pop%in%st]
st.d <- dist(t(gmat2))## distances between rows




