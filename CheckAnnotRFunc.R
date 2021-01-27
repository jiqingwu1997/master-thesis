
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("FDb.InfiniumMethylation.hg19")

#RESTART R AFTER INSTALLING IF CANNOT LOAD PACKAGE, WORKED FOR ME
library(FDb.InfiniumMethylation.hg19)

#set your proper working directory...

load("/Users/jiqingwu/Desktop/master thesis/hannum.rda")

probnames = hannum[[3]][[1]] #cpg names from hannum data set
hm450 <- get450k() #written in the manual, gets info of the assay
newprobnames = hm450[probnames] #relates your cpgs to their cpg data base

genes = getNearestGene(newprobnames) #get closest genes (or some are in genes already)

#View new genes!  
View(genes)

# Observe the columns "distance" and "nearestGeneSymbol"
#note that some CpGs have 0 distance from a gene, so from there, we can determine which CpGs are in a gene.  

vec1=c(rep("a",5),rep("b",5))
uni=unique(vec1)

all=data.frame(a=NA,b=NA)
for (i in 1:length(vec1)){
  
  
  dummy=rep(NA,length(uni))
  
  dummy[match(vec1[i],uni)]=1
  
  
  all=rbind(all,dummy)
  
}


load("/Users/jiqingwu/Desktop/master thesis/hannum.rda")
#This code finds the SD for each row:
  rsd <- matrixStats::rowSds(hannum$BVals,na.rm=T)
summary(rsd)
rfilt <- rank(-rsd)                 #rank smallest to largest
fidx <- which( rfilt <= 250000)

library(dplyr)
#This code filters on genes having at least 10 features:
  filt<- genes %>%
  group_by(nearestGeneSymbol) %>%
  summarize(n=n()) %>%
  filter(n>=10) 

Zannot <- genes %>% 
  filter(nearestGeneSymbol %in% filt$nearestGeneSymbol)
