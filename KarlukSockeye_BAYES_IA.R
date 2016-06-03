#### Karluk Smolt Analyses ####
# MSA for 2013 and 2014
# For each year, do 3 time strata of ~253 fish (760 total per year)
# Run all trap fish (760) from each year in BAYES with IA turned on, and also the 3 strata for each year (~253) with IA turned on
# 2014 Fyke net fish from weir are separate, IA
# Performing IA with BAYES and summarizing output at each thinned iteration (accounts for both genotyping and sampling error)
# Mon Apr 06 09:57:06 2015 Kyle Shedd

date()

ls()
rm(list=ls(all=TRUE))
search()
getwd()
setwd("V:/WORK/Sockeye/Kodiak/2013 2014 Karluk Smolt")

## save.image("V:/WORK/Sockeye/Kodiak/2013 2014 Karluk Smolt/KarlukSockeye_2013_2014_IA_RedoAbundances.RData")
## load("V:/WORK/Sockeye/Kodiak/2013 2014 Karluk Smolt/KarlukSockeye_2013_2014_IA_RedoAbundances.RData")

# Create folder directories
sapply(c("BAYES", "Objects", "Output", "Raw genotypes", "Tables", "Genepop"), dir.create)

# This sources all of the new GCL functions to this workspace
source("V:/DATA/R_GEN/GCL Source Scripts/Functions.GCL.r")
# This sources all of Kyle's non-GCL useful functions to this workspace (does not include "tempGCL" functions)
source("V:/WORK/Kyle/R Source Scripts/Functions.GCL_KS.R")


## Get collection SILLYs
KarlukSmolt <- c("SKARL13s","SKARL14s","SKARLW14s")
dput(x=KarlukSmolt,file="Objects/KarlukSmolt.txt")

## Pull all data for each silly code and create .gcl objects for each
ReadLOKI.GCL(sillyvec=KarlukSmolt,markersuite="Sockeye2011_96SNPs"); beep(sound="coin")
objects(pattern="\\.gcl")
unlist(strsplit(x=objects(pattern="\\.gcl"), split="\\.gcl"))

## Create a locus list from Locus Control object
loci96 <- LocusControl$locusnames
dput(x=loci96,file="Objects/loci96.txt")

## Save original LocusControl
dput(x=LocusControl,file="Objects/OriginalLocusControl.txt")

## Save unaltered .gcl's as back-up:
sapply(KarlukSmolt, function(silly) {dput(x=get(paste(silly,".gcl",sep='')),file=paste("Raw genotypes/",x,".txt",sep=''))})

## Old way if doesn't work
#for(Collection in KarlukSmolt){
  #dput(x=get(paste(Collection,".gcl",sep='')),file=paste("Raw genotypes/",Collection,".txt",sep=''))
#}

## Original sample sizes by SILLY
sapply(KarlukSmolt, function(silly) get(paste(silly, ".gcl", sep=""))$n)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata by time for 2013 and 2014 trap fish ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### 2013

## Stratum 1 -  5/16-29

unique(SKARL13s.gcl$attributes$CAPTURE_DATE)
lapply(list(1:8, 9:16, 17:27), function(lst) {sum(table(SKARL13s.gcl$attributes$CAPTURE_DATE)[lst])})

KarlukSmolt2013.1_IDs <- AttributesToIDs.GCL(silly="SKARL13s",attribute="CAPTURE_DATE",matching=unique(SKARL13s.gcl$attributes$CAPTURE_DATE)[1:8])
KarlukSmolt2013.1_IDs <- list(as.numeric(na.omit(KarlukSmolt2013.1_IDs)))
names(KarlukSmolt2013.1_IDs) <- "SKARL13s"

PoolCollections.GCL("SKARL13s",loci=loci96,IDs=KarlukSmolt2013.1_IDs,newname="KarlukSmolt2013.1")
KarlukSmolt2013.1.gcl$n ##  

## Stratum 2 -  5/30-6/10

KarlukSmolt2013.2_IDs <- AttributesToIDs.GCL(silly="SKARL13s",attribute="CAPTURE_DATE",matching=unique(SKARL13s.gcl$attributes$CAPTURE_DATE)[9:16])
KarlukSmolt2013.2_IDs <- list(as.numeric(na.omit(KarlukSmolt2013.2_IDs)))
names(KarlukSmolt2013.2_IDs) <- "SKARL13s"

PoolCollections.GCL("SKARL13s",loci=loci96,IDs=KarlukSmolt2013.2_IDs,newname="KarlukSmolt2013.2")
KarlukSmolt2013.2.gcl$n ##  

## Stratum 3 -  6/11-6/24

KarlukSmolt2013.3_IDs <- AttributesToIDs.GCL(silly="SKARL13s",attribute="CAPTURE_DATE",matching=unique(SKARL13s.gcl$attributes$CAPTURE_DATE)[17:27])
KarlukSmolt2013.3_IDs <- list(as.numeric(na.omit(KarlukSmolt2013.3_IDs)))
names(KarlukSmolt2013.3_IDs) <- "SKARL13s"

PoolCollections.GCL("SKARL13s",loci=loci96,IDs=KarlukSmolt2013.3_IDs,newname="KarlukSmolt2013.3")
KarlukSmolt2013.3.gcl$n ##  



#### 2014

## Stratum 1 =  5/13-30

unique(SKARL14s.gcl$attributes$CAPTURE_DATE)
lapply(list(1:13, 14:25, 26:38), function(lst) {sum(table(SKARL14s.gcl$attributes$CAPTURE_DATE)[lst])})

KarlukSmolt2014.1_IDs <- AttributesToIDs.GCL(silly="SKARL14s",attribute="CAPTURE_DATE",matching=unique(SKARL14s.gcl$attributes$CAPTURE_DATE)[1:13])
KarlukSmolt2014.1_IDs <- list(as.numeric(na.omit(KarlukSmolt2014.1_IDs)))
names(KarlukSmolt2014.1_IDs) <- "SKARL14s"

PoolCollections.GCL("SKARL14s",loci=loci96,IDs=KarlukSmolt2014.1_IDs,newname="KarlukSmolt2014.1")
KarlukSmolt2014.1.gcl$n ##  

## Stratum 2 -  5/31-6/15

KarlukSmolt2014.2_IDs <- AttributesToIDs.GCL(silly="SKARL14s",attribute="CAPTURE_DATE",matching=unique(SKARL14s.gcl$attributes$CAPTURE_DATE)[14:25])
KarlukSmolt2014.2_IDs <- list(as.numeric(na.omit(KarlukSmolt2014.2_IDs)))
names(KarlukSmolt2014.2_IDs) <- "SKARL14s"

PoolCollections.GCL("SKARL14s",loci=loci96,IDs=KarlukSmolt2014.2_IDs,newname="KarlukSmolt2014.2")
KarlukSmolt2014.2.gcl$n ##  

## Stratum 3 -  6/16-7/2

KarlukSmolt2014.3_IDs <- AttributesToIDs.GCL(silly="SKARL14s",attribute="CAPTURE_DATE",matching=unique(SKARL14s.gcl$attributes$CAPTURE_DATE)[26:38])
KarlukSmolt2014.3_IDs <- list(as.numeric(na.omit(KarlukSmolt2014.3_IDs)))
names(KarlukSmolt2014.3_IDs) <- "SKARL14s"

PoolCollections.GCL("SKARL14s",loci=loci96,IDs=KarlukSmolt2014.3_IDs,newname="KarlukSmolt2014.3")
KarlukSmolt2014.3.gcl$n ##  



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### QC ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(xlsx)

KarlukMixtures <- unlist(strsplit(ls(pattern='\\.gcl'),split='\\.gcl'))
dput(x=KarlukMixtures,file="Objects/KarlukMixtures.txt")

KarlukMixtures_SampleSizes <- matrix(data=NA,nrow=length(KarlukMixtures),ncol=5,dimnames=list(KarlukMixtures,c("Genotyped","Alternate","Missing","Duplicate","Final")))


#### Check loci
## Get sample size by locus
Original_KarlukMixtures_SampleSizebyLocus <- SampSizeByLocus.GCL(KarlukMixtures, loci96)
min(Original_KarlukMixtures_SampleSizebyLocus) ## 70
apply(Original_KarlukMixtures_SampleSizebyLocus,1,min)/apply(Original_KarlukMixtures_SampleSizebyLocus,1,max) # Known from QC that SKARLW14s had issues

t(sort(Original_KarlukMixtures_SampleSizebyLocus[KarlukMixtures[9],]/max(Original_KarlukMixtures_SampleSizebyLocus[KarlukMixtures[9],]),decreasing=TRUE)) # One_ZNF.61 is just below the 80% rule for SKARLW14s

#### Check individuals
### Initial
## Get number of individuals per silly before removing missing loci individuals
Original_KarlukMixtures_ColSize <- sapply(paste(KarlukMixtures,".gcl",sep=''), function(x) get(x)$n)
KarlukMixtures_SampleSizes[,"Genotyped"] <- Original_KarlukMixtures_ColSize


### Alternate
## Indentify alternate species individuals
KarlukMixtures_Alternate <- FindAlternateSpecies.GCL(sillyvec=KarlukMixtures, species="sockeye")

## Remove Alternate species individuals
RemoveAlternateSpecies.GCL(AlternateSpeciesReport=KarlukMixtures_Alternate, AlternateCutOff=0.5, FailedCutOff=0.5)

## Get number of individuals per silly after removing alternate species individuals
ColSize_KarlukMixtures_PostAlternate <- sapply(paste(KarlukMixtures,".gcl",sep=''), function(x) get(x)$n)
KarlukMixtures_SampleSizes[,"Alternate"] <- Original_KarlukMixtures_ColSize-ColSize_KarlukMixtures_PostAlternate


### Missing
## Remove individuals with >20% missing data
KarlukMixtures_MissLoci <- RemoveIndMissLoci.GCL(sillyvec=KarlukMixtures,loci=loci96,proportion=0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_KarlukMixtures_PostMissLoci <- sapply(paste(KarlukMixtures,".gcl",sep=''), function(x) get(x)$n)
KarlukMixtures_SampleSizes[,"Missing"] <- ColSize_KarlukMixtures_PostAlternate-ColSize_KarlukMixtures_PostMissLoci


### Duplicate
## Check within collections for duplicate individuals.
KarlukMixtures_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec=KarlukMixtures,loci=loci96,quantile=NULL,minproportion=0.95)
KarlukMixtures_DuplicateCheckReportSummary <- sapply(KarlukMixtures, function(x) KarlukMixtures_DuplicateCheck95MinProportion[[x]]$report)

## Remove duplicate individuals
KarlukMixtures_RemovedDups <- RemoveDups.GCL(KarlukMixtures_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_KarlukMixtures_PostDuplicate <- sapply(paste(KarlukMixtures,".gcl",sep=''), function(x) get(x)$n)
KarlukMixtures_SampleSizes[,"Duplicate"] <- ColSize_KarlukMixtures_PostMissLoci-ColSize_KarlukMixtures_PostDuplicate


### Final
KarlukMixtures_SampleSizes[,"Final"] <- ColSize_KarlukMixtures_PostDuplicate
KarlukMixtures_SampleSizes

write.xlsx(KarlukMixtures_SampleSizes,file="Output/KarlukMixtures_SampleSizes.xlsx")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Combine Loci ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get the final locus set from baseline file to avoid any errors
loci91 <- dget(file="V:/WORK/Sockeye/Kodiak/2014 Karluk Baseline/Objects/loci91.txt")

## Combine loci
CombineLoci.GCL(sillyvec=KarlukMixtures,markerset=c("One_Cytb_26","One_CO1","One_Cytb_17"),update=T)

## NOTE THAT THESE LOCI ARE DROPPED c("One_CO1","One_Cytb_17","One_Cytb_26","One_MHC2_251","One_GPDH","One_Tf_ex3-182")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### MSA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Create folder directories for BAYES
getwd()
invisible(sapply(c("Control","Mixture","Output"), function(x) {dir.create(paste(getwd(),"/BAYES/",x,sep=''))}))
invisible(sapply(c("Control","Mixture"), function(x) {dir.create(paste(getwd(),"/BAYES/",x,"/Done",sep=''))}))
## And Output for each Mixture
invisible(sapply(KarlukMixtures, function(x) {dir.create(paste(getwd(),"/BAYES/Output/",x,sep=''))}))


## Get baseline objects needed for MSA
KarlukGroups2 <- dget(file="V:/WORK/Sockeye/Kodiak/2014 Karluk Baseline/Objects/KarlukGroups2.txt")
Karluk16Pops2FlatPrior <- dget(file="V:/WORK/Sockeye/Kodiak/2014 Karluk Baseline/Objects/Karluk16Pops2FlatPrior.txt")
Karluk16PopsInits <- dget(file="V:/WORK/Sockeye/Kodiak/2014 Karluk Baseline/Objects/Karluk16PopsInits.txt")
Karluk1691Baseline <- dget(file="V:/WORK/Sockeye/Kodiak/2014 Karluk Baseline/Objects/Karluk1691Baseline.txt")
Karluk16Pops <- dget(file="V:/WORK/Sockeye/Kodiak/2014 Karluk Baseline/Objects/Karluk16Pops.txt")
Karluk16GroupVec2 <- dget(file="V:/WORK/Sockeye/Kodiak/2014 Karluk Baseline/Objects/Karluk16GroupVec2.txt")


## Defining the random seeds as the same as WASSIP mixtures for repeatability.
WASSIPSockeyeSeeds <- dget(file="V:/WORK/WASSIP/Sockeye/Mixture/Objects/WASSIPSockeyeSeeds.txt")

## Dump Mixture file format:
KarlukMixtureFormat <- CreateMixture.GCL(sillys=KarlukMixtures[1],loci=loci91,IDs=NULL,mixname=KarlukMixtures[1],dir=paste(getwd(),"/BAYES/Mixture",sep=''),type="BAYES",PT=FALSE)
dput(x=KarlukMixtureFormat,file="V:/WORK/Sockeye/Kodiak/2013 2014 Karluk Smolt/Objects/KarlukMixtureFormat.txt")

## Dump Mixture files
for(Mix in KarlukMixtures){
  CreateMixture.GCL(sillys=Mix,loci=loci91,IDs=NULL,mixname=Mix,dir=paste(getwd(),"/BAYES/Mixture",sep=''),type="BAYES",PT=FALSE)
}


## Dump Control files
for(Mix in KarlukMixtures){
  CreateControlFile.GCL(sillyvec=Karluk16Pops,loci=loci91,mixname=Mix,basename="Karluk16Pops91Markers",suffix="",nreps=40000,nchains=5,
                        groupvec=Karluk16GroupVec2,priorvec=Karluk16Pops2FlatPrior,initmat=Karluk16PopsInits,dir=paste(getwd(),"/BAYES/Control",sep=''),
                        seeds=WASSIPSockeyeSeeds,thin=c(1,1,100),mixfortran=KarlukMixtureFormat,basefortran=Karluk1691Baseline,switches="F T F T T T F")
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Import BAYES Results/Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Summarize mixtures
# By RG
KarlukSmolt_Estimates <- CustomCombineBAYESOutput.GCL(groupvec=1:2,groupnames=KarlukGroups2,maindir="BAYES/Output",mixvec=KarlukMixtures,prior="",ext="RGN",nchains=5,
                                                      burn=0.5,alpha=0.1,PosteriorOutput=FALSE); beep(sound=2)
dput(x=KarlukSmolt_Estimates,file="Objects/KarlukSmolt_Estimates.txt")

## Stratify by outmigration numbers
# 2013
KarlukSmolt_Outmigration_2013 <- as.numeric(readClipboard()) # These are the three outmigration estimates per strata from the "2013Karluk_smolt_pop_est_by_day_age" worksheet, sheet "Samples"
Karluk2013StratifiedEstimates <- StratifiedEstimator.GCL(groupvec=1:2,groupnames=KarlukGroups2,maindir="BAYES/Output",mixvec=KarlukMixtures[1:3],catchvec=KarlukSmolt_Outmigration_2013,
                                                          newname="Karluk2013StratifiedEstimates",priorname="",ext="RGN",nchains=5,burn=0.5,alpha=0.1); beep(sound=2)
Karluk2013StratifiedAbundances <- Karluk2013StratifiedEstimates$Summary*sum(KarlukSmolt_Outmigration_2013)


# 2014
KarlukSmolt_Outmigration_2014 <- as.numeric(readClipboard()) # These are the three outmigration estimates per strata from the "2014Karluk_smolt_pop_est_by_day_age" worksheet, sheet "Samples"
Karluk2014StratifiedEstimates <- StratifiedEstimator.GCL(groupvec=1:2,groupnames=KarlukGroups2,maindir="BAYES/Output",mixvec=KarlukMixtures[1:3],catchvec=KarlukSmolt_Outmigration_2014,
                                                         newname="Karluk2014StratifiedEstimates",priorname="",ext="RGN",nchains=5,burn=0.5,alpha=0.1); beep(sound=2)
Karluk2014StratifiedAbundances <- Karluk2014StratifiedEstimates$Summary*sum(KarlukSmolt_Outmigration_2014)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Exploratory Plots ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

opar=par()
require(gplots)
par(family='serif', mar=c(4.1, 4.1, 1.1, 1.1))
#### 2013 ####
### Proportion
## Both ages through time
barplot2(height=t(matrix(data=c(KarlukSmolt_Estimates$KarlukSmolt2013.1[,3],KarlukSmolt_Estimates$KarlukSmolt2013.2[,3],KarlukSmolt_Estimates$KarlukSmolt2013.3[,3]),nrow=3,ncol=2,byrow=TRUE,
                         dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2))),
         beside=TRUE,plot.ci=TRUE,ylim=c(0,1),ylab="Proportion",xlab="Stratum",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.2,
         ci.l=t(matrix(data=c(KarlukSmolt_Estimates$KarlukSmolt2013.1[,4],KarlukSmolt_Estimates$KarlukSmolt2013.2[,4],KarlukSmolt_Estimates$KarlukSmolt2013.3[,4]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2))),
         ci.u=t(matrix(data=c(KarlukSmolt_Estimates$KarlukSmolt2013.1[,5],KarlukSmolt_Estimates$KarlukSmolt2013.2[,5],KarlukSmolt_Estimates$KarlukSmolt2013.3[,5]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2))))
abline(h=0)
legend(x="topleft",legend=KarlukGroups2,fill=c("blue","red"),bty="n", cex=1.5)
mtext(text="2013", side=3, cex=2, line=-0.5)

## Transposed, periods in legend
barplot2(height=matrix(data=c(KarlukSmolt_Estimates$KarlukSmolt2013.1[,3],KarlukSmolt_Estimates$KarlukSmolt2013.2[,3],KarlukSmolt_Estimates$KarlukSmolt2013.3[,3]),nrow=3,ncol=2,byrow=TRUE,
                         dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2)),
         beside=TRUE,plot.ci=TRUE,ylim=c(0,1),ylab="Proportion",xlab="Stock",col=c("black", "grey", "white"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=matrix(data=c(KarlukSmolt_Estimates$KarlukSmolt2013.1[,4],KarlukSmolt_Estimates$KarlukSmolt2013.2[,4],KarlukSmolt_Estimates$KarlukSmolt2013.3[,4]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2)),
         ci.u=matrix(data=c(KarlukSmolt_Estimates$KarlukSmolt2013.1[,5],KarlukSmolt_Estimates$KarlukSmolt2013.2[,5],KarlukSmolt_Estimates$KarlukSmolt2013.3[,5]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2)))
abline(h=0)
legend(x="topleft",legend=c("May 16-29","May 30-June 10","June 10-June 24"),fill=c("black", "grey", "white"),bty="n", cex=1.5)
mtext(text="2013", side=3, cex=2, line=-0.5)

## Annual total
barplot2(height=Karluk2013StratifiedEstimates$Summary[,4],
         beside=TRUE,plot.ci=TRUE,ylim=c(0,1),ylab="Proportion of annual smolt outmigration",xlab="Stock",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=Karluk2013StratifiedEstimates$Summary[,3],
         ci.u=Karluk2013StratifiedEstimates$Summary[,5])
abline(h=0)
mtext(text="2013", side=3, cex=2, line=-0.5)



### Outmigrants
KarlukSmolt_Outmigrants_2013.1 <- KarlukSmolt_Estimates$KarlukSmolt2013.1 * KarlukSmolt_Outmigration_2013[1]
KarlukSmolt_Outmigrants_2013.2 <- KarlukSmolt_Estimates$KarlukSmolt2013.2 * KarlukSmolt_Outmigration_2013[2]
KarlukSmolt_Outmigrants_2013.3 <- KarlukSmolt_Estimates$KarlukSmolt2013.3 * KarlukSmolt_Outmigration_2013[3]

## Both ages through time
barplot2(height=t(matrix(data=c(KarlukSmolt_Outmigrants_2013.1[,3],KarlukSmolt_Outmigrants_2013.2[,3],KarlukSmolt_Outmigrants_2013.3[,3]),nrow=3,ncol=2,byrow=TRUE,
                         dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2)))/100000,
         beside=TRUE,plot.ci=TRUE,ylim=c(0,2.5),ylab="Outmigrating Smolt (100K)",xlab="Stratum",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.2,
         ci.l=t(matrix(data=c(KarlukSmolt_Outmigrants_2013.1[,4],KarlukSmolt_Outmigrants_2013.2[,4],KarlukSmolt_Outmigrants_2013.3[,4]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2)))/100000,
         ci.u=t(matrix(data=c(KarlukSmolt_Outmigrants_2013.1[,5],KarlukSmolt_Outmigrants_2013.2[,5],KarlukSmolt_Outmigrants_2013.3[,5]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2)))/100000)
abline(h=0)
legend(x="topleft",legend=KarlukGroups2,fill=c("blue","red"),bty="n", cex=1.5)
mtext(text="2013", side=3, cex=2, line=-0.5)

## Transposed, periods in legend
barplot2(height=matrix(data=c(KarlukSmolt_Outmigrants_2013.1[,3],KarlukSmolt_Outmigrants_2013.2[,3],KarlukSmolt_Outmigrants_2013.3[,3]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2))/100000,
         beside=TRUE,plot.ci=TRUE,ylim=c(0,2.5),ylab="Outmigrating Smolt (100K)",xlab="Stock",col=c("black", "grey", "white"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=matrix(data=c(KarlukSmolt_Outmigrants_2013.1[,4],KarlukSmolt_Outmigrants_2013.2[,4],KarlukSmolt_Outmigrants_2013.3[,4]),nrow=3,ncol=2,byrow=TRUE,
                     dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2))/100000,
         ci.u=matrix(data=c(KarlukSmolt_Outmigrants_2013.1[,5],KarlukSmolt_Outmigrants_2013.2[,5],KarlukSmolt_Outmigrants_2013.3[,5]),nrow=3,ncol=2,byrow=TRUE,
                     dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2))/100000)
abline(h=0)
legend(x="topleft",legend=c("May 16-29","May 30-June 10","June 10-June 24"),fill=c("black", "grey", "white"),bty="n", cex=1.5)
mtext(text="2013", side=3, cex=2, line=-0.5)

## Annual total
barplot2(height=Karluk2013StratifiedAbundances[,4]/100000,
         beside=TRUE,plot.ci=TRUE,ylim=c(0,7.0),ylab="Outmigrating Smolt (100K)",xlab="Stock",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=Karluk2013StratifiedAbundances[,3]/100000,
         ci.u=Karluk2013StratifiedAbundances[,5]/100000)
abline(h=0)
mtext(text="2013", side=3, cex=2, line=-0.5)
#### 2014 ####
### Proportion
## Both ages through time
barplot2(height=t(matrix(data=c(KarlukSmolt_Estimates$KarlukSmolt2014.1[,3],KarlukSmolt_Estimates$KarlukSmolt2014.2[,3],KarlukSmolt_Estimates$KarlukSmolt2014.3[,3]),nrow=3,ncol=2,byrow=TRUE,
                         dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2))),
         beside=TRUE,plot.ci=TRUE,ylim=c(0,1),ylab="Proportion",xlab="Stratum",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.2,
         ci.l=t(matrix(data=c(KarlukSmolt_Estimates$KarlukSmolt2014.1[,4],KarlukSmolt_Estimates$KarlukSmolt2014.2[,4],KarlukSmolt_Estimates$KarlukSmolt2014.3[,4]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2))),
         ci.u=t(matrix(data=c(KarlukSmolt_Estimates$KarlukSmolt2014.1[,5],KarlukSmolt_Estimates$KarlukSmolt2014.2[,5],KarlukSmolt_Estimates$KarlukSmolt2014.3[,5]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2))))
abline(h=0)
legend(x="topright",legend=KarlukGroups2,fill=c("blue","red"),bty="n", cex=1.5)
mtext(text="2014", side=3, cex=2, line=-0.5)

## Transposed, periods in legend
barplot2(height=matrix(data=c(KarlukSmolt_Estimates$KarlukSmolt2014.1[,3],KarlukSmolt_Estimates$KarlukSmolt2014.2[,3],KarlukSmolt_Estimates$KarlukSmolt2014.3[,3]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2)),
         beside=TRUE,plot.ci=TRUE,ylim=c(0,1),ylab="Proportion",xlab="Stock",col=c("black", "grey", "white"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=matrix(data=c(KarlukSmolt_Estimates$KarlukSmolt2014.1[,4],KarlukSmolt_Estimates$KarlukSmolt2014.2[,4],KarlukSmolt_Estimates$KarlukSmolt2014.3[,4]),nrow=3,ncol=2,byrow=TRUE,
                     dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2)),
         ci.u=matrix(data=c(KarlukSmolt_Estimates$KarlukSmolt2014.1[,5],KarlukSmolt_Estimates$KarlukSmolt2014.2[,5],KarlukSmolt_Estimates$KarlukSmolt2014.3[,5]),nrow=3,ncol=2,byrow=TRUE,
                     dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2)))
abline(h=0)
legend(x="topleft",legend=c("May 13-30","May 31-June 15","June 16-July 2"),fill=c("black", "grey", "white"),bty="n", cex=1.5)
mtext(text="2014", side=3, cex=2, line=-0.5)

## Annual total
barplot2(height=Karluk2014StratifiedEstimates$Summary[,4],
         beside=TRUE,plot.ci=TRUE,ylim=c(0,1),ylab="Proportion of annual smolt outmigration",xlab="Stock",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=Karluk2014StratifiedEstimates$Summary[,3],
         ci.u=Karluk2014StratifiedEstimates$Summary[,5])
abline(h=0)
mtext(text="2014", side=3, cex=2, line=-0.5)



### Outmigrants
KarlukSmolt_Outmigrants_2014.1 <- KarlukSmolt_Estimates$KarlukSmolt2014.1 * KarlukSmolt_Outmigration_2014[1]
KarlukSmolt_Outmigrants_2014.2 <- KarlukSmolt_Estimates$KarlukSmolt2014.2 * KarlukSmolt_Outmigration_2014[2]
KarlukSmolt_Outmigrants_2014.3 <- KarlukSmolt_Estimates$KarlukSmolt2014.3 * KarlukSmolt_Outmigration_2014[3]

## Both ages through time
barplot2(height=t(matrix(data=c(KarlukSmolt_Outmigrants_2014.1[,3],KarlukSmolt_Outmigrants_2014.2[,3],KarlukSmolt_Outmigrants_2014.3[,3]),nrow=3,ncol=2,byrow=TRUE,
                         dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2)))/100000,
         beside=TRUE,plot.ci=TRUE,ylim=c(0,2.5),ylab="Outmigrating Smolt (100K)",xlab="Stratum",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.2,
         ci.l=t(matrix(data=c(KarlukSmolt_Outmigrants_2014.1[,4],KarlukSmolt_Outmigrants_2014.2[,4],KarlukSmolt_Outmigrants_2014.3[,4]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2)))/100000,
         ci.u=t(matrix(data=c(KarlukSmolt_Outmigrants_2014.1[,5],KarlukSmolt_Outmigrants_2014.2[,5],KarlukSmolt_Outmigrants_2014.3[,5]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2)))/100000)
abline(h=0)
legend(x="topright",legend=KarlukGroups2,fill=c("blue","red"),bty="n", cex=1.5)
mtext(text="2014", side=3, cex=2, line=-0.5)

## Transposed, periods in legend
barplot2(height=matrix(data=c(KarlukSmolt_Outmigrants_2014.1[,3],KarlukSmolt_Outmigrants_2014.2[,3],KarlukSmolt_Outmigrants_2014.3[,3]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2))/100000,
         beside=TRUE,plot.ci=TRUE,ylim=c(0,2.5),ylab="Outmigrating Smolt (100K)",xlab="Stock",col=c("black", "grey", "white"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=matrix(data=c(KarlukSmolt_Outmigrants_2014.1[,4],KarlukSmolt_Outmigrants_2014.2[,4],KarlukSmolt_Outmigrants_2014.3[,4]),nrow=3,ncol=2,byrow=TRUE,
                     dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2))/100000,
         ci.u=matrix(data=c(KarlukSmolt_Outmigrants_2014.1[,5],KarlukSmolt_Outmigrants_2014.2[,5],KarlukSmolt_Outmigrants_2014.3[,5]),nrow=3,ncol=2,byrow=TRUE,
                     dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2))/100000)
abline(h=0)
legend(x="topleft",legend=c("May 13-30","May 31-June 15","June 16-July 2"),fill=c("black", "grey", "white"),bty="n", cex=1.5)
mtext(text="2014", side=3, cex=2, line=-0.5)

## Annual total
barplot2(height=Karluk2014StratifiedAbundances[,4]/100000,
         beside=TRUE,plot.ci=TRUE,ylim=c(0,7.0),ylab="Outmigrating Smolt (100K)",xlab="Stock",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=Karluk2014StratifiedAbundances[,3]/100000,
         ci.u=Karluk2014StratifiedAbundances[,5]/100000)
abline(h=0)
mtext(text="2014", side=3, cex=2, line=-0.5)
#### 2014 Weir ####
barplot2(height=KarlukSmolt_Estimates$SKARLW14s[,3],
         beside=TRUE,plot.ci=TRUE,ylim=c(0,1),ylab="Proportion",xlab="Stock",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=KarlukSmolt_Estimates$SKARLW14s[,4],
         ci.u=KarlukSmolt_Estimates$SKARLW14s[,5])
abline(h=0)
legend(x="topleft",legend=KarlukGroups2,fill=c("blue","red"),bty="n", cex=1.5)
mtext(text="2014 - Weir", side=3, cex=2, line=-0.5)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Incorporate Age Specific Mixtures ####
## These should have been done previously, but were overlooked
## NOTE: the QC step will be skipped as all "problem" fish have already been dealt with in the original SILLYs
## Strata dates will be kept the same as from the original temporal strata with all ages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##### Add Metadata
### 2013
# Metadata inputs
Metadata2013 <- read.table(file="KarlukMetadata2013.csv", as.is=TRUE, header=TRUE, sep=",")
str(Metadata2013)

sapply(KarlukMixtures[c(1:3,7)], function(silly) {assign(x=paste(silly, "_Metadata", sep=""), value=Metadata2013[match(get(paste(silly, ".gcl", sep=""))$attributes$SillySource, Metadata2013$SILLY),], pos=1)})

### 2014
# Metadata inputs
Metadata2014 <- read.table(file="KarlukMetadata2014.csv", as.is=TRUE, header=TRUE, sep=",")
str(Metadata2014)

sapply(KarlukMixtures[c(4:6,8)], function(silly) {assign(x=paste(silly, "_Metadata", sep=""), value=Metadata2014[match(get(paste(silly, ".gcl", sep=""))$attributes$SillySource, Metadata2014$SILLY),], pos=1)})

# Weir
SKARLW14s_Metadata <- read.table(file="KarlukMetadata2014Weir.csv", as.is=TRUE, header=TRUE, sep=",")

## Incorporate Age into Attributes table for SKARL13s and SKARL14s
SKARL13s.gcl$attributes$AGE <- SKARL13s_Metadata$AGE
SKARL14s.gcl$attributes$AGE <- SKARL14s_Metadata$AGE

KarlukSmolt2013.1.gcl$attributes$AGE <- KarlukSmolt2013.1_Metadata$AGE
KarlukSmolt2013.2.gcl$attributes$AGE <- KarlukSmolt2013.2_Metadata$AGE
KarlukSmolt2013.3.gcl$attributes$AGE <- KarlukSmolt2013.3_Metadata$AGE

KarlukSmolt2014.1.gcl$attributes$AGE <- KarlukSmolt2014.1_Metadata$AGE
KarlukSmolt2014.2.gcl$attributes$AGE <- KarlukSmolt2014.2_Metadata$AGE
KarlukSmolt2014.3.gcl$attributes$AGE <- KarlukSmolt2014.3_Metadata$AGE




objects(pattern="\\.gcl")


#### 2013
table(SKARL13s.gcl$attributes$AGE)

## Age 1 - 5/16-6/24 (only 1 Age 1 strata)
KarlukSmolt2013.Age1_IDs <- AttributesToIDs.GCL(silly="SKARL13s",attribute="AGE",matching=1)
KarlukSmolt2013.Age1_IDs <- list(as.numeric(na.omit(KarlukSmolt2013.Age1_IDs)))
names(KarlukSmolt2013.Age1_IDs) <- "SKARL13s"

PoolCollections.GCL("SKARL13s",loci=loci96,IDs=KarlukSmolt2013.Age1_IDs,newname="KarlukSmolt2013.Age1")
KarlukSmolt2013.Age1.gcl$n ## 123


## Age 2 Stratum 1 - 5/16-29
KarlukSmolt2013.Age2.1_IDs <- AttributesToIDs.GCL(silly="KarlukSmolt2013.1",attribute="AGE",matching=2)
KarlukSmolt2013.Age2.1_IDs <- list(as.numeric(na.omit(KarlukSmolt2013.Age2.1_IDs)))
names(KarlukSmolt2013.Age2.1_IDs) <- "SKARL13s"

PoolCollections.GCL("SKARL13s",loci=loci96,IDs=KarlukSmolt2013.Age2.1_IDs,newname="KarlukSmolt2013.Age2.1")
KarlukSmolt2013.Age2.1.gcl$n ## 189

## Age 2 Stratum 2 - 5/30-6/10
KarlukSmolt2013.Age2.2_IDs <- AttributesToIDs.GCL(silly="KarlukSmolt2013.2",attribute="AGE",matching=2)
KarlukSmolt2013.Age2.2_IDs <- list(as.numeric(na.omit(KarlukSmolt2013.Age2.2_IDs)))
names(KarlukSmolt2013.Age2.2_IDs) <- "SKARL13s"

PoolCollections.GCL("SKARL13s",loci=loci96,IDs=KarlukSmolt2013.Age2.2_IDs,newname="KarlukSmolt2013.Age2.2")
KarlukSmolt2013.Age2.2.gcl$n ## 227

## Age 2 Stratum 3 - 6/11-6/24
KarlukSmolt2013.Age2.3_IDs <- AttributesToIDs.GCL(silly="KarlukSmolt2013.3",attribute="AGE",matching=2)
KarlukSmolt2013.Age2.3_IDs <- list(as.numeric(na.omit(KarlukSmolt2013.Age2.3_IDs)))
names(KarlukSmolt2013.Age2.3_IDs) <- "SKARL13s"

PoolCollections.GCL("SKARL13s",loci=loci96,IDs=KarlukSmolt2013.Age2.3_IDs,newname="KarlukSmolt2013.Age2.3")
KarlukSmolt2013.Age2.3.gcl$n ## 141


## Age 3 - 5/16-6/24 (only 1 Age 3 strata)
KarlukSmolt2013.Age3_IDs <- AttributesToIDs.GCL(silly="SKARL13s",attribute="AGE",matching=3)
KarlukSmolt2013.Age3_IDs <- list(as.numeric(na.omit(KarlukSmolt2013.Age3_IDs)))
names(KarlukSmolt2013.Age3_IDs) <- "SKARL13s"

PoolCollections.GCL("SKARL13s",loci=loci96,IDs=KarlukSmolt2013.Age3_IDs,newname="KarlukSmolt2013.Age3")
KarlukSmolt2013.Age3.gcl$n ## 65




#### 2014
table(SKARL14s.gcl$attributes$AGE)

## Age 1 Stratum 12 5/13-6/15 (strata 1 and 2 collapsed for sample size reasons)
KarlukSmolt2014.Age1.12_IDs <- c(AttributesToIDs.GCL(silly="KarlukSmolt2014.1",attribute="AGE",matching=1), AttributesToIDs.GCL(silly="KarlukSmolt2014.2",attribute="AGE",matching=1))
KarlukSmolt2014.Age1.12_IDs <- list(as.numeric(na.omit(KarlukSmolt2014.Age1.12_IDs)))
names(KarlukSmolt2014.Age1.12_IDs) <- "SKARL14s"

PoolCollections.GCL("SKARL14s",loci=loci96,IDs=KarlukSmolt2014.Age1.12_IDs,newname="KarlukSmolt2014.Age1.12")
KarlukSmolt2014.Age1.12.gcl$n ## 87

## Age 1 Stratum 3 6/16-7/2
KarlukSmolt2014.Age1.3_IDs <- AttributesToIDs.GCL(silly="KarlukSmolt2014.3",attribute="AGE",matching=1)
KarlukSmolt2014.Age1.3_IDs <- list(as.numeric(na.omit(KarlukSmolt2014.Age1.3_IDs)))
names(KarlukSmolt2014.Age1.3_IDs) <- "SKARL14s"

PoolCollections.GCL("SKARL14s",loci=loci96,IDs=KarlukSmolt2014.Age1.3_IDs,newname="KarlukSmolt2014.Age1.3")
KarlukSmolt2014.Age1.3.gcl$n ## 163

## Age 2 Stratum 1 - 5/13-30
KarlukSmolt2014.Age2.1_IDs <- AttributesToIDs.GCL(silly="KarlukSmolt2014.1",attribute="AGE",matching=2)
KarlukSmolt2014.Age2.1_IDs <- list(as.numeric(na.omit(KarlukSmolt2014.Age2.1_IDs)))
names(KarlukSmolt2014.Age2.1_IDs) <- "SKARL14s"

PoolCollections.GCL("SKARL14s",loci=loci96,IDs=KarlukSmolt2014.Age2.1_IDs,newname="KarlukSmolt2014.Age2.1")
KarlukSmolt2014.Age2.1.gcl$n ## 229

## Age 2 Stratum 2 - 5/31-6/15
KarlukSmolt2014.Age2.2_IDs <- AttributesToIDs.GCL(silly="KarlukSmolt2014.2",attribute="AGE",matching=2)
KarlukSmolt2014.Age2.2_IDs <- list(as.numeric(na.omit(KarlukSmolt2014.Age2.2_IDs)))
names(KarlukSmolt2014.Age2.2_IDs) <- "SKARL14s"

PoolCollections.GCL("SKARL14s",loci=loci96,IDs=KarlukSmolt2014.Age2.2_IDs,newname="KarlukSmolt2014.Age2.2")
KarlukSmolt2014.Age2.2.gcl$n ## 177

## Age 2 Stratum 3 - 6/16-7/2
KarlukSmolt2014.Age2.3_IDs <- AttributesToIDs.GCL(silly="KarlukSmolt2014.3",attribute="AGE",matching=2)
KarlukSmolt2014.Age2.3_IDs <- list(as.numeric(na.omit(KarlukSmolt2014.Age2.3_IDs)))
names(KarlukSmolt2014.Age2.3_IDs) <- "SKARL14s"

PoolCollections.GCL("SKARL14s",loci=loci96,IDs=KarlukSmolt2014.Age2.3_IDs,newname="KarlukSmolt2014.Age2.3")
KarlukSmolt2014.Age2.3.gcl$n ## 88


# sillyvec for only Age specific mixtures
KarlukMixtures_Age <- unlist(strsplit(ls(pattern='\\.gcl'),split='\\.gcl'))[grep(x=unlist(strsplit(ls(pattern='\\.gcl'),split='\\.gcl')), pattern="Age")]
dput(x=KarlukMixtures_Age,file="Objects/KarlukMixtures_Age.txt")

# Sample sizes for mixtures...
sapply(KarlukMixtures_Age, function(silly) get(paste(silly, ".gcl", sep=""))$n)


## Combine loci
CombineLoci.GCL(sillyvec=KarlukMixtures_Age,markerset=c("One_Cytb_26","One_CO1","One_Cytb_17"),update=T)

# Confirm that loci combined correctly (should be 97 loci now)
sapply(KarlukMixtures_Age, function(mix) {dim(get(paste(mix, ".gcl", sep=""))$counts)}) # cool
## NOTE THAT THESE LOCI ARE DROPPED c("One_CO1","One_Cytb_17","One_Cytb_26","One_MHC2_251","One_GPDH","One_Tf_ex3-182")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### MSA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Create folder directories for BAYES
getwd()
#invisible(sapply(c("Control","Mixture","Output"), function(x) {dir.create(paste(getwd(),"/BAYES/",x,sep=''))}))
invisible(sapply(c("Control","Mixture"), function(x) {dir.create(paste(getwd(),"/BAYES/",x,"/Done",sep=''))}))
## And Output for each Mixture
invisible(sapply(KarlukMixtures_Age, function(x) {dir.create(paste(getwd(),"/BAYES/Output/",x,sep=''))}))


## Dump Mixture files
for(Mix in KarlukMixtures_Age){
  CreateMixture.GCL(sillys=Mix,loci=loci91,IDs=NULL,mixname=Mix,dir=paste(getwd(),"/BAYES/Mixture",sep=''),type="BAYES",PT=FALSE)
}


## Dump Control files
for(Mix in KarlukMixtures_Age){
  CreateControlFile.GCL(sillyvec=Karluk16Pops,loci=loci91,mixname=Mix,basename="Karluk16Pops91Markers",suffix="",nreps=40000,nchains=5,
                        groupvec=Karluk16GroupVec2,priorvec=Karluk16Pops2FlatPrior,initmat=Karluk16PopsInits,dir=paste(getwd(),"/BAYES/Control",sep=''),
                        seeds=WASSIPSockeyeSeeds,thin=c(1,1,100),mixfortran=KarlukMixtureFormat,basefortran=Karluk1691Baseline,switches="F T F T T T F")
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES ... AGAIN!!! ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Import BAYES Results/Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Summarize mixtures
# By RG
KarlukSmolt_Age_Estimates <- CustomCombineBAYESOutput.GCL(groupvec=1:2,groupnames=KarlukGroups2,maindir="BAYES/Output",mixvec=KarlukMixtures_Age,prior="",ext="RGN",nchains=5,
                                                      burn=0.5,alpha=0.1,PosteriorOutput=FALSE); beep(sound=2)
dput(x=KarlukSmolt_Age_Estimates,file="Objects/KarlukSmolt_Age_Estimates.txt")

## Stratify by outmigration numbers
# 2013
# Age 2
KarlukSmolt_Age2_Outmigration_2013 <- as.numeric(readClipboard()) # These are the three outmigration estimates per strata from the "2013Karluk_smolt_pop_est_by_day_age" worksheet, sheet "Samples"
Karluk2013Age2StratifiedEstimates <- StratifiedEstimator.GCL(groupvec=1:2,groupnames=KarlukGroups2,maindir="BAYES/Output",mixvec=KarlukMixtures_Age[2:4],catchvec=KarlukSmolt_Age2_Outmigration_2013,
                                                         newname="Karluk2013Age2StratifiedEstimates",priorname="",ext="RGN",nchains=5,burn=0.5,alpha=0.1); beep(sound=2)
Karluk2013Age2StratifiedAbundances <- Karluk2013Age2StratifiedEstimates$Summary*sum(KarlukSmolt_Age2_Outmigration_2013)


# 2014
# Age 1
KarlukSmolt_Age1_Outmigration_2014 <- as.numeric(readClipboard()) # These are the three outmigration estimates per strata from the "2014Karluk_smolt_pop_est_by_day_age" worksheet, sheet "Samples"
Karluk2014Age1StratifiedEstimates <- StratifiedEstimator.GCL(groupvec=1:2,groupnames=KarlukGroups2,maindir="BAYES/Output",mixvec=KarlukMixtures_Age[6:7],catchvec=KarlukSmolt_Age1_Outmigration_2014,
                                                             newname="Karluk2014Age1StratifiedEstimates",priorname="",ext="RGN",nchains=5,burn=0.5,alpha=0.1); beep(sound=2)
Karluk2014Age1StratifiedAbundances <- Karluk2014Age1StratifiedEstimates$Summary*sum(KarlukSmolt_Age1_Outmigration_2014)

# Age 2
KarlukSmolt_Age2_Outmigration_2014 <- as.numeric(readClipboard()) # These are the three outmigration estimates per strata from the "2014Karluk_smolt_pop_est_by_day_age" worksheet, sheet "Samples"
Karluk2014Age2StratifiedEstimates <- StratifiedEstimator.GCL(groupvec=1:2,groupnames=KarlukGroups2,maindir="BAYES/Output",mixvec=KarlukMixtures_Age[8:10],catchvec=KarlukSmolt_Age2_Outmigration_2014,
                                                             newname="Karluk2014Age2StratifiedEstimates",priorname="",ext="RGN",nchains=5,burn=0.5,alpha=0.1); beep(sound=2)
Karluk2014Age2StratifiedAbundances <- Karluk2014Age2StratifiedEstimates$Summary*sum(KarlukSmolt_Age2_Outmigration_2014)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Exploratory Plots of Age Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

opar=par()
require(gplots)
par(family='serif', mar=c(4.1, 4.1, 1.1, 1.1))
#### 2013 ####
### Proportion
## Age 2
barplot2(height=t(matrix(data=c(KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.1[,3],KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.2[,3],KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.3[,3]),nrow=3,ncol=2,byrow=TRUE,
                         dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2))),
         beside=TRUE,plot.ci=TRUE,ylim=c(0,1),ylab="Proportion",xlab="Stratum",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.2,
         ci.l=t(matrix(data=c(KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.1[,4],KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.2[,4],KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.3[,4]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2))),
         ci.u=t(matrix(data=c(KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.1[,5],KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.2[,5],KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.3[,5]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2))))
abline(h=0)
legend(x="topleft",legend=KarlukGroups2,fill=c("blue","red"),bty="n", cex=1.5)
mtext(text="Age 2 - 2013", side=3, cex=2, line=-0.5)

## Annual total
# Age 1
barplot2(height=KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age1[,3],
         beside=TRUE,plot.ci=TRUE,ylim=c(0,1),ylab="Proportion of annual smolt outmigration",xlab="Stock",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age1[,4],
         ci.u=KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age1[,5])
abline(h=0)
mtext(text="Age 1 - 2013", side=3, cex=2, line=-0.5)

# Age 2
barplot2(height=Karluk2013Age2StratifiedEstimates$Summary[,4],
         beside=TRUE,plot.ci=TRUE,ylim=c(0,1),ylab="Proportion of annual smolt outmigration",xlab="Stock",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=Karluk2013Age2StratifiedEstimates$Summary[,3],
         ci.u=Karluk2013Age2StratifiedEstimates$Summary[,5])
abline(h=0)
mtext(text="Age 2 - 2013", side=3, cex=2, line=-0.5)

# Age 3
barplot2(height=KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age3[,3],
         beside=TRUE,plot.ci=TRUE,ylim=c(0,1),ylab="Proportion of annual smolt outmigration",xlab="Stock",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age3[,4],
         ci.u=KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age3[,5])
abline(h=0)
mtext(text="Age 3 - 2013", side=3, cex=2, line=-0.5)


### Transposed, periods in legend
## Age 2
barplot2(height=matrix(data=c(KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.1[,3],KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.2[,3],KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.3[,3]),nrow=3,ncol=2,byrow=TRUE,
                         dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2)),
         beside=TRUE,plot.ci=TRUE,ylim=c(0,1),ylab="Proportion",xlab="Stratum",col=c("black", "grey", "white"), cex.axis=1.5, cex.lab=1.5, cex.names=1.2,
         ci.l=matrix(data=c(KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.1[,4],KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.2[,4],KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.3[,4]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2)),
         ci.u=matrix(data=c(KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.1[,5],KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.2[,5],KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.3[,5]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2)))
abline(h=0)
legend(x="topleft",legend=c("May 16-29","May 30-June 10","June 10-June 24"),fill=c("black", "grey", "white"),bty="n", cex=1.5)
mtext(text="Age 2 - 2013", side=3, cex=2, line=-0.5)

### Outmigrants
KarlukSmolt_Age2_Outmigrants_2013.1 <- KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.1 * KarlukSmolt_Age2_Outmigration_2013[1]
KarlukSmolt_Age2_Outmigrants_2013.2 <- KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.2 * KarlukSmolt_Age2_Outmigration_2013[2]
KarlukSmolt_Age2_Outmigrants_2013.3 <- KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age2.3 * KarlukSmolt_Age2_Outmigration_2013[3]

## Age 2
barplot2(height=t(matrix(data=c(KarlukSmolt_Age2_Outmigrants_2013.1[,3],KarlukSmolt_Age2_Outmigrants_2013.2[,3],KarlukSmolt_Age2_Outmigrants_2013.3[,3]),nrow=3,ncol=2,byrow=TRUE,
                         dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2)))/100000,
         beside=TRUE,plot.ci=TRUE,ylim=c(0,2.5),ylab="Outmigrating Smolt (100K)",xlab="Stratum",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.2,
         ci.l=t(matrix(data=c(KarlukSmolt_Age2_Outmigrants_2013.1[,4],KarlukSmolt_Age2_Outmigrants_2013.2[,4],KarlukSmolt_Age2_Outmigrants_2013.3[,4]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2)))/100000,
         ci.u=t(matrix(data=c(KarlukSmolt_Age2_Outmigrants_2013.1[,5],KarlukSmolt_Age2_Outmigrants_2013.2[,5],KarlukSmolt_Age2_Outmigrants_2013.3[,5]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2)))/100000)
abline(h=0)
legend(x="topleft",legend=KarlukGroups2,fill=c("blue","red"),bty="n", cex=1.5)
mtext(text="Age 2 - 2013", side=3, cex=2, line=-0.5)

## Annual total
# Age 1
KarlukSmolt_Age1_Outmigration_2013 <- as.numeric(readClipboard())
KarlukSmolt_Age1_Outmigrants_2013 <- KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age1 * KarlukSmolt_Age1_Outmigration_2013

barplot2(height=KarlukSmolt_Age1_Outmigrants_2013[,3]/100000,
         beside=TRUE,plot.ci=TRUE,ylim=c(0,7.0),ylab="Outmigrating Smolt (100K)",xlab="Stock",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=KarlukSmolt_Age1_Outmigrants_2013[,4]/100000,
         ci.u=KarlukSmolt_Age1_Outmigrants_2013[,5]/100000)
abline(h=0)
mtext(text="Age 1 - 2013", side=3, cex=2, line=-0.5)

# Age2
barplot2(height=Karluk2013Age2StratifiedAbundances[,4]/100000,
         beside=TRUE,plot.ci=TRUE,ylim=c(0,7.0),ylab="Outmigrating Smolt (100K)",xlab="Stock",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=Karluk2013Age2StratifiedAbundances[,3]/100000,
         ci.u=Karluk2013Age2StratifiedAbundances[,5]/100000)
abline(h=0)
mtext(text="Age 2 - 2013", side=3, cex=2, line=-0.5)

# Age 3
KarlukSmolt_Age3_Outmigration_2013 <- as.numeric(readClipboard())
KarlukSmolt_Age3_Outmigrants_2013 <- KarlukSmolt_Age_Estimates$KarlukSmolt2013.Age3 * KarlukSmolt_Age3_Outmigration_2013

barplot2(height=KarlukSmolt_Age3_Outmigrants_2013[,3]/100000,
         beside=TRUE,plot.ci=TRUE,ylim=c(0,7.0),ylab="Outmigrating Smolt (100K)",xlab="Stock",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=KarlukSmolt_Age3_Outmigrants_2013[,4]/100000,
         ci.u=KarlukSmolt_Age3_Outmigrants_2013[,5]/100000)
abline(h=0)
mtext(text="Age 3 - 2013", side=3, cex=2, line=-0.5)


### Transposed, periods in legend
## Age 2
barplot2(height=matrix(data=c(KarlukSmolt_Age2_Outmigrants_2013.1[,3],KarlukSmolt_Age2_Outmigrants_2013.2[,3],KarlukSmolt_Age2_Outmigrants_2013.3[,3]),nrow=3,ncol=2,byrow=TRUE,
                         dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2))/100000,
         beside=TRUE,plot.ci=TRUE,ylim=c(0,2.5),ylab="Outmigrating Smolt (100K)",xlab="Stratum",col=c("black", "grey", "white"), cex.axis=1.5, cex.lab=1.5, cex.names=1.2,
         ci.l=matrix(data=c(KarlukSmolt_Age2_Outmigrants_2013.1[,4],KarlukSmolt_Age2_Outmigrants_2013.2[,4],KarlukSmolt_Age2_Outmigrants_2013.3[,4]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2))/100000,
         ci.u=matrix(data=c(KarlukSmolt_Age2_Outmigrants_2013.1[,5],KarlukSmolt_Age2_Outmigrants_2013.2[,5],KarlukSmolt_Age2_Outmigrants_2013.3[,5]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 16-29","May 30-June 10","June 10-June 24"),KarlukGroups2))/100000)
abline(h=0)
legend(x="topleft",legend=c("May 16-29","May 30-June 10","June 10-June 24"),fill=c("black", "grey", "white"),bty="n", cex=1.5)
mtext(text="Age 2 - 2013", side=3, cex=2, line=-0.5)
#### HERE ####
#### 2014 ####
### Proportion
## Age 1
barplot2(height=t(matrix(data=c(KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age1.12[,3],KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age1.3[,3]),nrow=2,ncol=2,byrow=TRUE,
                         dimnames=list(c("May 13-June 15","June 16-July 2"),KarlukGroups2))),
         beside=TRUE,plot.ci=TRUE,ylim=c(0,1),ylab="Proportion",xlab="Stratum",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.2,
         ci.l=t(matrix(data=c(KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age1.12[,4],KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age1.3[,4]),nrow=2,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-June 15","June 16-July 2"),KarlukGroups2))),
         ci.u=t(matrix(data=c(KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age1.12[,5],KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age1.3[,5]),nrow=2,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-June 15","June 16-July 2"),KarlukGroups2))))
abline(h=0)
legend(x="topright",legend=KarlukGroups2,fill=c("blue","red"),bty="n", cex=1.5)
mtext(text="Age 1 - 2014", side=3, cex=2, line=-0.5)

## Age 2
barplot2(height=t(matrix(data=c(KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.1[,3],KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.2[,3],KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.3[,3]),nrow=3,ncol=2,byrow=TRUE,
                         dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2))),
         beside=TRUE,plot.ci=TRUE,ylim=c(0,1),ylab="Proportion",xlab="Stratum",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.2,
         ci.l=t(matrix(data=c(KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.1[,4],KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.2[,4],KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.3[,4]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2))),
         ci.u=t(matrix(data=c(KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.1[,5],KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.2[,5],KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.3[,5]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2))))
abline(h=0)
legend(x="topright",legend=KarlukGroups2,fill=c("blue","red"),bty="n", cex=1.5)
mtext(text="Age 2 - 2014", side=3, cex=2, line=-0.5)

## Annual total
# Age 1
barplot2(height=Karluk2014Age1StratifiedEstimates$Summary[,4],
         beside=TRUE,plot.ci=TRUE,ylim=c(0,1),ylab="Proportion of annual smolt outmigration",xlab="Stock",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=Karluk2014Age1StratifiedEstimates$Summary[,3],
         ci.u=Karluk2014Age1StratifiedEstimates$Summary[,5])
abline(h=0)
mtext(text="Age 1 - 2014", side=3, cex=2, line=-0.5)

# Age 2
barplot2(height=Karluk2014Age2StratifiedEstimates$Summary[,4],
         beside=TRUE,plot.ci=TRUE,ylim=c(0,1),ylab="Proportion of annual smolt outmigration",xlab="Stock",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=Karluk2014Age2StratifiedEstimates$Summary[,3],
         ci.u=Karluk2014Age2StratifiedEstimates$Summary[,5])
abline(h=0)
mtext(text="Age 2 - 2014", side=3, cex=2, line=-0.5)


### Transposed, periods in legend
## Age 1
barplot2(height=matrix(data=c(KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age1.12[,3],KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age1.3[,3]),nrow=2,ncol=2,byrow=TRUE,
                         dimnames=list(c("May 13-June 15","June 16-July 2"),KarlukGroups2)),
         beside=TRUE,plot.ci=TRUE,ylim=c(0,1),ylab="Proportion",xlab="Stratum",col=c("grey37", "white"), cex.axis=1.5, cex.lab=1.5, cex.names=1.2,
         ci.l=matrix(data=c(KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age1.12[,4],KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age1.3[,4]),nrow=2,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-June 15","June 16-July 2"),KarlukGroups2)),
         ci.u=matrix(data=c(KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age1.12[,5],KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age1.3[,5]),nrow=2,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-June 15","June 16-July 2"),KarlukGroups2)))
abline(h=0)
legend(x="topright",legend=c("May 13-June 15","June 16-July 2"),fill=c("grey37", "white"),bty="n", cex=1.5)
mtext(text="Age 1 - 2014", side=3, cex=2, line=-0.5)

## Age 2
barplot2(height=matrix(data=c(KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.1[,3],KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.2[,3],KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.3[,3]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2)),
         beside=TRUE,plot.ci=TRUE,ylim=c(0,1),ylab="Proportion",xlab="Stratum",col=c("black", "grey", "white"), cex.axis=1.5, cex.lab=1.5, cex.names=1.2,
         ci.l=matrix(data=c(KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.1[,4],KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.2[,4],KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.3[,4]),nrow=3,ncol=2,byrow=TRUE,
                     dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2)),
         ci.u=matrix(data=c(KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.1[,5],KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.2[,5],KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.3[,5]),nrow=3,ncol=2,byrow=TRUE,
                     dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2)))
abline(h=0)
legend(x="topleft",legend=c("May 13-30","May 31-June 15","June 16-July 2"),fill=c("black", "grey", "white"),bty="n", cex=1.5)
mtext(text="Age 2 - 2014", side=3, cex=2, line=-0.5)



### Outmigrants
KarlukSmolt_Age1_Outmigrants_2014.12 <- KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age1.12 * KarlukSmolt_Age1_Outmigration_2014[1]
KarlukSmolt_Age1_Outmigrants_2014.3 <- KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age1.3 * KarlukSmolt_Age1_Outmigration_2014[2]

## Age 1
barplot2(height=t(matrix(data=c(KarlukSmolt_Age1_Outmigrants_2014.12[,3],KarlukSmolt_Age1_Outmigrants_2014.3[,3]),nrow=2,ncol=2,byrow=TRUE,
                         dimnames=list(c("May 13-June 15","June 16-July 2"),KarlukGroups2)))/100000,
         beside=TRUE,plot.ci=TRUE,ylim=c(0,2.5),ylab="Outmigrating Smolt (100K)",xlab="Stratum",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.2,
         ci.l=t(matrix(data=c(KarlukSmolt_Age1_Outmigrants_2014.12[,4],KarlukSmolt_Age1_Outmigrants_2014.3[,4]),nrow=2,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-June 15","June 16-July 2"),KarlukGroups2)))/100000,
         ci.u=t(matrix(data=c(KarlukSmolt_Age1_Outmigrants_2014.12[,5],KarlukSmolt_Age1_Outmigrants_2014.3[,5]),nrow=2,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-June 15","June 16-July 2"),KarlukGroups2)))/100000)
abline(h=0)
legend(x="topright",legend=KarlukGroups2,fill=c("blue","red"),bty="n", cex=1.5)
mtext(text="Age 1 - 2014", side=3, cex=2, line=-0.5)


KarlukSmolt_Age2_Outmigrants_2014.1 <- KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.1 * KarlukSmolt_Age2_Outmigration_2014[1]
KarlukSmolt_Age2_Outmigrants_2014.2 <- KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.2 * KarlukSmolt_Age2_Outmigration_2014[2]
KarlukSmolt_Age2_Outmigrants_2014.3 <- KarlukSmolt_Age_Estimates$KarlukSmolt2014.Age2.3 * KarlukSmolt_Age2_Outmigration_2014[3]

## Age 2
barplot2(height=t(matrix(data=c(KarlukSmolt_Age2_Outmigrants_2014.1[,3],KarlukSmolt_Age2_Outmigrants_2014.2[,3],KarlukSmolt_Age2_Outmigrants_2014.3[,3]),nrow=3,ncol=2,byrow=TRUE,
                         dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2)))/100000,
         beside=TRUE,plot.ci=TRUE,ylim=c(0,2.5),ylab="Outmigrating Smolt (100K)",xlab="Stratum",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.2,
         ci.l=t(matrix(data=c(KarlukSmolt_Age2_Outmigrants_2014.1[,4],KarlukSmolt_Age2_Outmigrants_2014.2[,4],KarlukSmolt_Age2_Outmigrants_2014.3[,4]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2)))/100000,
         ci.u=t(matrix(data=c(KarlukSmolt_Age2_Outmigrants_2014.1[,5],KarlukSmolt_Age2_Outmigrants_2014.2[,5],KarlukSmolt_Age2_Outmigrants_2014.3[,5]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2)))/100000)
abline(h=0)
legend(x="topright",legend=KarlukGroups2,fill=c("blue","red"),bty="n", cex=1.5)
mtext(text="Age 2 - 2014", side=3, cex=2, line=-0.5)

## Annual total
# Age 1
barplot2(height=Karluk2014Age1StratifiedAbundances[,4]/100000,
         beside=TRUE,plot.ci=TRUE,ylim=c(0,7.0),ylab="Outmigrating Smolt (100K)",xlab="Stock",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=Karluk2014Age1StratifiedAbundances[,3]/100000,
         ci.u=Karluk2014Age1StratifiedAbundances[,5]/100000)
abline(h=0)
mtext(text="Age 1 - 2014", side=3, cex=2, line=-0.5)

# Age2
barplot2(height=Karluk2014Age2StratifiedAbundances[,4]/100000,
         beside=TRUE,plot.ci=TRUE,ylim=c(0,7.0),ylab="Outmigrating Smolt (100K)",xlab="Stock",col=c("blue","red"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
         ci.l=Karluk2014Age2StratifiedAbundances[,3]/100000,
         ci.u=Karluk2014Age2StratifiedAbundances[,5]/100000)
abline(h=0)
mtext(text="Age 2 - 2014", side=3, cex=2, line=-0.5)


### Transposed, periods in legend
## Age 1
barplot2(height=matrix(data=c(KarlukSmolt_Age1_Outmigrants_2014.12[,3],KarlukSmolt_Age1_Outmigrants_2014.3[,3]),nrow=2,ncol=2,byrow=TRUE,
                         dimnames=list(c("May 13-June 15","June 16-July 2"),KarlukGroups2))/100000,
         beside=TRUE,plot.ci=TRUE,ylim=c(0,2.5),ylab="Outmigrating Smolt (100K)",xlab="Stratum",col=c("grey37","white"), cex.axis=1.5, cex.lab=1.5, cex.names=1.2,
         ci.l=matrix(data=c(KarlukSmolt_Age1_Outmigrants_2014.12[,4],KarlukSmolt_Age1_Outmigrants_2014.3[,4]),nrow=2,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-June 15","June 16-July 2"),KarlukGroups2))/100000,
         ci.u=matrix(data=c(KarlukSmolt_Age1_Outmigrants_2014.12[,5],KarlukSmolt_Age1_Outmigrants_2014.3[,5]),nrow=2,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-June 15","June 16-July 2"),KarlukGroups2))/100000)
abline(h=0)
legend(x="topright",legend=c("May 13-June 15","June 16-July 2"),fill=c("grey37","white"),bty="n", cex=1.5)
mtext(text="Age 1 - 2014", side=3, cex=2, line=-0.5)

## Age 2
barplot2(height=matrix(data=c(KarlukSmolt_Age2_Outmigrants_2014.1[,3],KarlukSmolt_Age2_Outmigrants_2014.2[,3],KarlukSmolt_Age2_Outmigrants_2014.3[,3]),nrow=3,ncol=2,byrow=TRUE,
                       dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2))/100000,
         beside=TRUE,plot.ci=TRUE,ylim=c(0,2.5),ylab="Outmigrating Smolt (100K)",xlab="Stratum",col=c("black", "grey", "white"), cex.axis=1.5, cex.lab=1.5, cex.names=1.2,
         ci.l=matrix(data=c(KarlukSmolt_Age2_Outmigrants_2014.1[,4],KarlukSmolt_Age2_Outmigrants_2014.2[,4],KarlukSmolt_Age2_Outmigrants_2014.3[,4]),nrow=3,ncol=2,byrow=TRUE,
                     dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2))/100000,
         ci.u=matrix(data=c(KarlukSmolt_Age2_Outmigrants_2014.1[,5],KarlukSmolt_Age2_Outmigrants_2014.2[,5],KarlukSmolt_Age2_Outmigrants_2014.3[,5]),nrow=3,ncol=2,byrow=TRUE,
                     dimnames=list(c("May 13-30","May 31-June 15","June 16-July 2"),KarlukGroups2))/100000)
abline(h=0)
legend(x="topleft",legend=c("May 13-30","May 31-June 15","June 16-July 2"),fill=c("black", "grey", "white"),bty="n", cex=1.5)
mtext(text="Age 2 - 2014", side=3, cex=2, line=-0.5)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plots of Strata Dates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### 2013

require(plotrix)

Ymd.format <- "%Y/%m/%d"

Ymd <- function(x){ as.POSIXct(strptime(x, format=Ymd.format))}
gantt.info <- list(
  labels     =c("Stratum 1","Stratum 2","Stratum 3","Age 1","Age 2 Stratum 1","Age 2 Stratum 2","Age 2 Stratum 3","Age 3"),
  starts     =Ymd(c("2013/05/16","2013/05/30","2013/06/11","2013/05/16","2013/05/16","2013/05/30","2013/06/11","2013/05/16")),
  ends       =Ymd(c("2013/05/29","2013/06/10","2013/06/24","2013/06/24","2013/05/29","2013/06/10","2013/06/24","2013/06/24")),
  priorities =c(1,2,3,4,5,6,7,8))

par(family='serif')
gantt.chart(gantt.info, label.cex=1.5, vgridpos=Ymd(c("2013/05/16","2013/05/30","2013/06/11","2013/06/24")), xlim=Ymd(c("2013/05/14","2013/06/26")))
mtext(text="2013 Karluk Smolt Strata", side=3, cex=2, line=2.5)


### 2014

Ymd <- function(x){ as.POSIXct(strptime(x, format=Ymd.format))}
gantt.info <- list(
  labels     =c("Stratum 1","Stratum 2","Stratum 3","Age 1 Stratum 1-2", "Age 1 Stratum 3","Age 2 Stratum 1","Age 2 Stratum 2","Age 2 Stratum 3"),
  starts     =Ymd(c("2014/05/13","2014/05/31","2014/06/16","2014/05/13","2014/06/16","2014/05/13","2014/05/31","2014/06/16")),
  ends       =Ymd(c("2014/05/30","2014/06/15","2014/07/02","2014/06/15","2014/07/02","2014/05/30","2014/06/15","2014/07/02")),
  priorities =c(1,2,3,4,5,6,7,8))

par(family='serif')
gantt.chart(gantt.info, label.cex=1.5, vgridpos=Ymd(c("2014/05/13","2014/05/31","2014/06/16","2014/07/02")), xlim=Ymd(c("2014/05/11","2014/07/05")))
mtext(text="2014 Karluk Smolt Strata", side=3, cex=2, line=2.5)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Individual Assignments with BAYES ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source(file="V:/DATA/R_GEN/tempGCL Source Scripts/IndividualAssignmentSummaryKS.GCL.R")

KarlukSmolt_IndividualAssignments = 
  IndividualAssignmentSummaryKS.GCL(GroupNames=KarlukGroups2,groupvec=Karluk16GroupVec2,mixnames=KarlukMixtures,
                                    BAYESoutputDir="V:/WORK/Sockeye/Kodiak/2013 2014 Karluk Smolt/BAYES/Output",nchains=5,nreps=40000,
                                    burn=1/2,thin=100); beep(sound=2)

str(KarlukSmolt_IndividualAssignments)

### 2013
# Metadata inputs
Metadata2013 <- read.table(file="KarlukMetadata2013.csv", as.is=TRUE, header=TRUE, sep=",")
str(Metadata2013)

sapply(KarlukMixtures[c(1:3,7)], function(silly) {assign(x=paste(silly, "_Metadata", sep=""), value=Metadata2013[match(get(paste(silly, ".gcl", sep=""))$attributes$SillySource, Metadata2013$SILLY),], pos=1)})

### 2014
# Metadata inputs
Metadata2014 <- read.table(file="KarlukMetadata2014.csv", as.is=TRUE, header=TRUE, sep=",")
str(Metadata2014)

sapply(KarlukMixtures[c(4:6,8:9)], function(silly) {assign(x=paste(silly, "_Metadata", sep=""), value=Metadata2014[match(get(paste(silly, ".gcl", sep=""))$attributes$SillySource, Metadata2014$SILLY),], pos=1)})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Calculate straight likelihoods for individuals ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(beepr)
## Get baseline objects
groupnames <- dget(file="V:/WORK/Sockeye/Kodiak/2014 Karluk Baseline/Objects/KarlukGroups2.txt")
groupvec <- dget(file="V:/WORK/Sockeye/Kodiak/2014 Karluk Baseline/Objects/Karluk16GroupVec2.txt")
sillyvec_baseline <- dget(file="V:/WORK/Sockeye/Kodiak/2014 Karluk Baseline/Objects/Karluk16Pops.txt")
loci <- loci91
C=as.vector(table(groupvec))
G=length(C)
popP <- setNames(object=rep(x=1/C/G, C), nm=sillyvec_baseline[order(groupvec)])[sillyvec_baseline] # Use to assure equal probability to each RG, not to each Pop

nsilly=length(sillyvec_baseline)
nloci=length(loci)
alleles=LocusControl$alleles[loci]
nalleles=LocusControl$nalleles[loci]
ploidy=LocusControl$ploidy[loci]


## Read in BAYES baseline .bse file to get baseline allele frequencies
baseline <- scan(file="V:/WORK/Sockeye/Kodiak/2014 Karluk Baseline/BAYES/Baseline/Karluk16Pops91Markers.bse")
str(baseline)

baseline_mat <- matrix(data=baseline, nrow=nsilly*nloci, ncol=1+1+1+max(nalleles), byrow=TRUE)
colnames(baseline_mat) <- c("RG", "loci", "n", paste("Allele ", 1:max(nalleles), sep=""))
str(baseline_mat)
head(baseline_mat)


## Get baseline allele frequencies
q=array(NA,c(nsilly,nloci,max(nalleles)),dimnames=list(sillyvec_baseline,loci,paste("Allele ",1:max(nalleles),sep="")))
for(silly in sillyvec_baseline){
  for(locus in loci){
    q[silly,locus,1:nalleles[locus]]=(baseline_mat[which(baseline_mat[,"loci"]==which(loci==locus) & baseline_mat[,"RG"]==which(sillyvec_baseline==silly)), 4:(3+nalleles[locus])] + 1/nalleles[locus]) / (baseline_mat[which(baseline_mat[,"loci"]==which(loci==locus) & baseline_mat[,"RG"]==which(sillyvec_baseline==silly)), "n"] + 1)
  }
}

## Get mixture genotypes

sillyvec <- KarlukMixtures

my.gcl=lapply(sillyvec,function(silly){get(paste(silly,".gcl",sep=""),pos=1)})
names(my.gcl)=sillyvec  

n=sapply(my.gcl,function(gcl){gcl$n})
names(n)=sillyvec

counts=lapply(my.gcl,function(gcl){gcl$counts[,loci,]})
names(counts)=sillyvec

Y=FreqPop.GCL(sillyvec,loci)

## Likelihood function
genofunc=function(x,theta,NAind){return(prod(sapply(loci[NAind],function(locus){2*dmultinom(x=x[locus,1:nalleles[locus]],size=ploidy[locus],prob=theta[locus,1:nalleles[locus]])})))}

T=lapply(sillyvec_baseline,function(silly1){TT=lapply(sillyvec,function(silly2){rep(NA,n[silly2])});names(TT)=sillyvec;TT});names(T)=sillyvec_baseline
P=lapply(sillyvec_baseline,function(silly1){PP=lapply(sillyvec,function(silly2){rep(NA,n[silly2])});names(PP)=sillyvec;PP});names(P)=sillyvec_baseline

for(silly1 in sillyvec){
  for(m in seq(n[silly1])){
    x=counts[[silly1]][m,1:nloci,1:max(nalleles)]
    NAind=!is.na(x[1:nloci,1])
    theta=q[1:nsilly,1:nloci,1:max(nalleles)]
    for(silly2 in sillyvec_baseline){
      P[[silly2]][[silly1]][m]=popP[[silly2]]*genofunc(x=x,theta=theta[silly2,1:nloci,1:max(nalleles)],NAind)
    }  
    for(silly2 in sillyvec_baseline){
      T[[silly2]][[silly1]][m]=P[[silly2]][[silly1]][m]/sum(sapply(sillyvec_baseline,function(silly){P[[silly]][[silly1]][m]}))
    }  
  }
}


RT=lapply(seq(G),function(g){RTT=lapply(sillyvec,function(silly1){rep(NA,n[silly1])});names(RTT)=sillyvec;RTT});names(RT)=groupnames

for(silly1 in sillyvec){
  for(m in seq(n[silly1])){
    for(g in seq(G)){
      rt=0
      for(silly2 in sillyvec_baseline[which(groupvec==g)]){
        rt=rt+T[[silly2]][[silly1]][m] 
      }
      RT[[groupnames[g]]][[silly1]][m]=rt
    }
  }
}; beep(2)


KarlukSmolt_IndividualAssignments_Likelihood <- setNames(lapply(KarlukMixtures, function(mix) {matrix(data=round(c(unlist(RT[[1]][mix]), unlist(RT[[2]][mix])),3), ncol=2, dimnames=list(seq(length(unlist(RT[[1]][mix]))), groupnames))}), KarlukMixtures)
str(KarlukSmolt_IndividualAssignments_Likelihood)

# Reformat to match output of IndividualAssignmentKS.GCL
KarlukSmolt_IndividualAssignments_Likelihood<- sapply(1:9, function(i) {matrix(data=round(c(unlist(RT[[1]][i]), unlist(RT[[2]][i])),3), ncol=2, dimnames=list(names(unlist(RT[[1]][i])), groupnames))})
str(KarlukSmolt_IndividualAssignments_Likelihood)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Compare BAYES IA to straight Likelihoods ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

str(KarlukSmolt_IndividualAssignments)
str(KarlukSmolt_IndividualAssignments_Likelihood)

cbind(unlist(KarlukSmolt_IndividualAssignments[[1]]), unlist(KarlukSmolt_IndividualAssignments_Likelihood[[1]]))

# Number assigned to each group at threshold
threshold <- 0.9
lapply(KarlukSmolt_IndividualAssignments, function(silly) {c("Early"=length(which(silly[,1]>=threshold)), "Late"=length(which(silly[,2]>=threshold)), "Total"=nrow(silly))})
lapply(KarlukSmolt_IndividualAssignments_Likelihood, function(silly) {c("Early"=length(which(silly[,1]>=threshold)), "Late"=length(which(silly[,2]>=threshold)), "Total"=nrow(silly))})

# Proportions
lapply(lapply(KarlukSmolt_IndividualAssignments, function(silly) {c("Early"=length(which(silly[,1]>=threshold)), "Late"=length(which(silly[,2]>=threshold)), "Total"=nrow(silly))}), function(mix) {mix["Early"]/sum(mix["Early"],mix["Late"])})
lapply(lapply(KarlukSmolt_IndividualAssignments_Likelihood, function(silly) {c("Early"=length(which(silly[,1]>=threshold)), "Late"=length(which(silly[,2]>=threshold)), "Total"=nrow(silly))}), function(mix) {mix["Early"]/sum(mix["Early"],mix["Late"])})

# Likelihood produces proportions that are MUCH closer to the mixture estimates we get from BAYES

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Use Likelihoods to Assign Individuals ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NOTE: we want to provide the "hard assignments" as well as stock specific estimates of the % unclassified and % misclassified

### 2013
## "Relaxed" threshold (0.80)
relax <- 0.8

str(Metadata2013)
Metadata2013$RelaxedAssignment <- rep(as.vector(x=NA, mode="character"), dim(Metadata2013)[1])

# Confirm number of fish
length(as.character(KarlukSmolt2013.1.gcl$attributes$SillySource))
dim(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2013.1)[1]

# Match fish and assign based on likelihood
Metadata2013$RelaxedAssignment[match(as.character(KarlukSmolt2013.1.gcl$attributes$SillySource), Metadata2013$SILLY)] <- ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2013.1[,"KarlukEarly"] >= relax, "Early", 
                                                                                                                               ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2013.1[,"KarlukLate"] >= relax,"Late", NA))
Metadata2013$RelaxedAssignment[match(as.character(KarlukSmolt2013.2.gcl$attributes$SillySource), Metadata2013$SILLY)] <- ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2013.2[,"KarlukEarly"] >= relax, "Early", 
                                                                                                                               ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2013.2[,"KarlukLate"] >= relax,"Late", NA))
Metadata2013$RelaxedAssignment[match(as.character(KarlukSmolt2013.3.gcl$attributes$SillySource), Metadata2013$SILLY)] <- ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2013.3[,"KarlukEarly"] >= relax, "Early", 
                                                                                                                               ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2013.3[,"KarlukLate"] >= relax,"Late", NA))

# Make sure it is okay
head(Metadata2013,20)


## "Strict" threshold (0.95)
strict <- 0.95

Metadata2013$StrictAssignment <- rep(as.vector(x=NA, mode="character"), dim(Metadata2013)[1])

# Match fish and assign based on likelihood
Metadata2013$StrictAssignment[match(as.character(KarlukSmolt2013.1.gcl$attributes$SillySource), Metadata2013$SILLY)] <- ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2013.1[,"KarlukEarly"] >= strict, "Early", 
                                                                                                                               ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2013.1[,"KarlukLate"] >= strict,"Late", NA))
Metadata2013$StrictAssignment[match(as.character(KarlukSmolt2013.2.gcl$attributes$SillySource), Metadata2013$SILLY)] <- ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2013.2[,"KarlukEarly"] >= strict, "Early", 
                                                                                                                               ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2013.2[,"KarlukLate"] >= strict,"Late", NA))
Metadata2013$StrictAssignment[match(as.character(KarlukSmolt2013.3.gcl$attributes$SillySource), Metadata2013$SILLY)] <- ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2013.3[,"KarlukEarly"] >= strict, "Early", 
                                                                                                                               ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2013.3[,"KarlukLate"] >= strict,"Late", NA))

# Make sure it is okay
head(Metadata2013,20)

# How many fish assigned by age
table(Metadata2013$StrictAssignment, Metadata2013$AGE)
table(Metadata2013$RelaxedAssignment, Metadata2013$AGE)

# Add Genetics Number column
Metadata2013$GeneticsNum <- gsub(pattern="SKARL13s_", replacement="", x=Metadata2013$SILLY)
str(Metadata2013)


# Write results
write.csv(x=Metadata2013, file="KarlukSmoltAssignments2013.csv", row.names=FALSE)




### 2014
## "Relaxed" threshold (0.80)
relax <- 0.8

str(Metadata2014)
Metadata2014$RelaxedAssignment <- rep(as.vector(x=NA, mode="character"), dim(Metadata2014)[1])

# Confirm number of fish
length(as.character(KarlukSmolt2014.1.gcl$attributes$SillySource))
dim(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2014.1)[1]

# Match fish and assign based on likelihood
Metadata2014$RelaxedAssignment[match(as.character(KarlukSmolt2014.1.gcl$attributes$SillySource), Metadata2014$SILLY)] <- ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2014.1[,"KarlukEarly"] >= relax, "Early", 
                                                                                                                                ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2014.1[,"KarlukLate"] >= relax,"Late", NA))
Metadata2014$RelaxedAssignment[match(as.character(KarlukSmolt2014.2.gcl$attributes$SillySource), Metadata2014$SILLY)] <- ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2014.2[,"KarlukEarly"] >= relax, "Early", 
                                                                                                                                ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2014.2[,"KarlukLate"] >= relax,"Late", NA))
Metadata2014$RelaxedAssignment[match(as.character(KarlukSmolt2014.3.gcl$attributes$SillySource), Metadata2014$SILLY)] <- ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2014.3[,"KarlukEarly"] >= relax, "Early", 
                                                                                                                                ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2014.3[,"KarlukLate"] >= relax,"Late", NA))

# Make sure it is okay
head(Metadata2014,20)


## "Strict" threshold (0.95)
strict <- 0.95

Metadata2014$StrictAssignment <- rep(as.vector(x=NA, mode="character"), dim(Metadata2014)[1])

# Match fish and assign based on likelihood
Metadata2014$StrictAssignment[match(as.character(KarlukSmolt2014.1.gcl$attributes$SillySource), Metadata2014$SILLY)] <- ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2014.1[,"KarlukEarly"] >= strict, "Early", 
                                                                                                                               ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2014.1[,"KarlukLate"] >= strict,"Late", NA))
Metadata2014$StrictAssignment[match(as.character(KarlukSmolt2014.2.gcl$attributes$SillySource), Metadata2014$SILLY)] <- ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2014.2[,"KarlukEarly"] >= strict, "Early", 
                                                                                                                               ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2014.2[,"KarlukLate"] >= strict,"Late", NA))
Metadata2014$StrictAssignment[match(as.character(KarlukSmolt2014.3.gcl$attributes$SillySource), Metadata2014$SILLY)] <- ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2014.3[,"KarlukEarly"] >= strict, "Early", 
                                                                                                                               ifelse(KarlukSmolt_IndividualAssignments_Likelihood$KarlukSmolt2014.3[,"KarlukLate"] >= strict,"Late", NA))

# Make sure it is okay
head(Metadata2014,20)
tail(Metadata2014,20)

# How many fish assigned by age
table(Metadata2014$StrictAssignment, Metadata2014$AGE)
table(Metadata2014$RelaxedAssignment, Metadata2014$AGE)

# Add Genetics Number column
Metadata2014$GeneticsNum <- gsub(pattern="SKARL14s_", replacement="", x=Metadata2014$SILLY)
str(Metadata2014)


# Write results
write.csv(x=Metadata2014, file="KarlukSmoltAssignments2014.csv", row.names=FALSE)






### 2014 WEIR FISH
## "Relaxed" threshold (0.80)
relax <- 0.8

str(SKARLW14s_Metadata)
SKARLW14s_Metadata$RelaxedAssignment <- rep(as.vector(x=NA, mode="character"), dim(SKARLW14s_Metadata)[1])

# Confirm number of fish
length(as.character(SKARLW14s.gcl$attributes$SillySource))
dim(KarlukSmolt_IndividualAssignments_Likelihood$SKARLW14s)[1]

# Match fish and assign based on likelihood
SKARLW14s_Metadata$RelaxedAssignment[match(as.character(SKARLW14s.gcl$attributes$SillySource), SKARLW14s_Metadata$SILLY)] <- ifelse(KarlukSmolt_IndividualAssignments_Likelihood$SKARLW14s[,"KarlukEarly"] >= relax, "Early", 
                                                                                                                                ifelse(KarlukSmolt_IndividualAssignments_Likelihood$SKARLW14s[,"KarlukLate"] >= relax,"Late", NA))

# Make sure it is okay
head(SKARLW14s_Metadata,20)


## "Strict" threshold (0.95)
strict <- 0.95

SKARLW14s_Metadata$StrictAssignment <- rep(as.vector(x=NA, mode="character"), dim(SKARLW14s_Metadata)[1])

# Match fish and assign based on likelihood
SKARLW14s_Metadata$StrictAssignment[match(as.character(SKARLW14s.gcl$attributes$SillySource), SKARLW14s_Metadata$SILLY)] <- ifelse(KarlukSmolt_IndividualAssignments_Likelihood$SKARLW14s[,"KarlukEarly"] >= strict, "Early", 
                                                                                                                                    ifelse(KarlukSmolt_IndividualAssignments_Likelihood$SKARLW14s[,"KarlukLate"] >= strict,"Late", NA))

# Make sure it is okay
head(SKARLW14s_Metadata,20)
tail(SKARLW14s_Metadata,20)

# How many fish assigned by age
table(SKARLW14s_Metadata$StrictAssignment, SKARLW14s_Metadata$AGE)
table(SKARLW14s_Metadata$RelaxedAssignment, SKARLW14s_Metadata$AGE)

# Add Genetics Number column
SKARLW14s_Metadata$GeneticsNum <- gsub(pattern="SKARLW14s_", replacement="", x=SKARLW14s_Metadata$SILLY)
str(SKARLW14s_Metadata)


# Write results
write.csv(x=SKARLW14s_Metadata, file="KarlukSmoltAssignments2014Weir.csv", row.names=FALSE)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Make Tables for MSA Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ls(pattern="Estimates") 

ls(pattern="Abundances") 

ls(pattern='migrants')

objects(pattern="Outmigration")

objects(pattern="\\.gcl")

# List of abundance estimates (proportion * outmigration)
KarlukAbundances <- setNames(object=lapply(ls(pattern='migrants'), get), nm=gsub(pattern=".gcl", replacement="", x=objects(pattern="\\.gcl"))[c(4,12,13,5:7,14:16,8,1:3,9:11)])
str(KarlukAbundances)

# Vector of outmigration numers
KarlukOutmigration <- setNames(object=unlist(sapply(objects(pattern="Outmigration")[c(6,1,3,5,7,2,4)], get)), nm=gsub(pattern=".gcl", replacement="", x=objects(pattern="\\.gcl"))[-c(17:19)])

# List of abundance estimates (proportion * outmigration)
KarlukFinalSampleSizes <- setNames(object=unlist(lapply(objects(pattern="\\.gcl"), function(silly) {get(silly)$n})), nm=gsub(pattern=".gcl", replacement="", x=objects(pattern="\\.gcl")))
str(KarlukFinalSampleSizes)

# Vector of reporting group names
KarlukGroups <- gsub(pattern="Karluk", replacement="", x=KarlukGroups2)

# Empty vector of CVs
KarlukCVs <- setNames(object=rep(NA, times=length(gsub(pattern=".gcl", replacement="", x=objects(pattern="\\.gcl")))), nm=gsub(pattern=".gcl", replacement="", x=objects(pattern="\\.gcl")))

# Named vecotr of strata dates
KarlukShortDates <- setNames(object=c("5/16-29", "5/30-6/10", "6/11-24", "5/16-6/24", "5/16-29", "5/30-6/10", "6/11-24", "5/16-6/24", 
                                      "5/13-30", "5/31-6/15", "6/16-7/2", "5/13-6/15", "6/16-7/2", "5/13-30", "5/31-6/15", "6/16-7/2"), 
                             nm=gsub(pattern=".gcl", replacement="", x=objects(pattern="\\.gcl"))[-c(17:19)])

AnnualShortDates <- setNames(object=c("5/16-6/24", "5/13-7/2"), nm=2013:2014)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Table 1, 2013 smolt ####
## 17 rows by 15 columns for three strata/year and annual roll up.

KarlukCaption1 <- "Table X.-Estimates of stock composition and stock-specific outmigration for Karluk River sockeye salmon smolt by stratum, 2013.  Reporting group-specific stock composition estimates (%) include median, 90% credibility interval, the probability that the group estimate is equal to zero (P=0), mean and SD.  Stock-specific estimates of outmigration are based upon mark-recapture estimates and variances (CV) of outmigration for each stratum.  See text for details."
Disclaimer <- cbind("Note: Stock composition estimates may not sum to 100% and stock-specific outmigration estimates may not sum to the total outmigration due to rounding error.",'','','','','','','','','','','','','','')

KarlukTable1=cbind(KarlukCaption1,'','','','','','','','','','','','','','')
KarlukTable1=rbind(KarlukTable1,cbind("Stratum","",'','Composition (%)','','','','','','',"Outmigration (number of fish)",'','','',''))
KarlukTable1=rbind(KarlukTable1,cbind('Period','Sample Size','Reporting','',"90% CI",'','','','','','',"90% CI",'','',''))
KarlukTable1=rbind(KarlukTable1,cbind("Dates",'CV','Group',"Median","5%","95%","P=0","Mean","SD",'',"Median","5%","95%","Mean","SD"))
## Early KarlukEarly
KarlukTable1=rbind(KarlukTable1,cbind('Early',paste('n=',KarlukFinalSampleSizes['KarlukSmolt2013.1'],sep=''),KarlukGroups[1],
                                        formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.1']][1,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.1']][1,4]*100,format="f",digits=1),
                                        formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.1']][1,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.1']][1,6],format="f",digits=2),
                                        formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.1']][1,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.1']][1,2]*100,format="f",digits=1),'',
                                        formatC(x=KarlukAbundances[['KarlukSmolt2013.1']][1,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.1']][1,4],format="f",digits=0,big.mark=","),
                                        formatC(x=KarlukAbundances[['KarlukSmolt2013.1']][1,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.1']][1,1],
                                                                                                                                                            format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.1']][1,2],format="f",digits=0,big.mark=",")))
## Early KarlukLate
KarlukTable1=rbind(KarlukTable1,cbind(KarlukShortDates['KarlukSmolt2013.1'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2013.1']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                        formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.1']][2,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.1']][2,4]*100,format="f",digits=1),
                                        formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.1']][2,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.1']][2,6],format="f",digits=2),
                                        formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.1']][2,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.1']][2,2]*100,format="f",digits=1),'',
                                        formatC(x=KarlukAbundances[['KarlukSmolt2013.1']][2,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.1']][2,4],format="f",digits=0,big.mark=","),
                                        formatC(x=KarlukAbundances[['KarlukSmolt2013.1']][2,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.1']][2,1],
                                                                                                                                                            format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.1']][2,2],format="f",digits=0,big.mark=",")))
## Early Total Outmigration row
KarlukTable1=rbind(KarlukTable1,cbind('','','','','','','','','','','','Early Total','',formatC(x=KarlukOutmigration['KarlukSmolt2013.1'],format="f",digits=0,big.mark=","),''))
## Middle KarlukEarly
KarlukTable1=rbind(KarlukTable1,cbind('Middle',paste('n=',KarlukFinalSampleSizes['KarlukSmolt2013.2'],sep=''),KarlukGroups[1],
                                        formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.2']][1,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.2']][1,4]*100,format="f",digits=1),
                                        formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.2']][1,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.2']][1,6],format="f",digits=2),
                                        formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.2']][1,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.2']][1,2]*100,format="f",digits=1),'',
                                        formatC(x=KarlukAbundances[['KarlukSmolt2013.2']][1,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.2']][1,4],format="f",digits=0,big.mark=","),
                                        formatC(x=KarlukAbundances[['KarlukSmolt2013.2']][1,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.2']][1,1],
                                                                                                                                                            format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.2']][1,2],format="f",digits=0,big.mark=",")))
## Middle KarlukLate
KarlukTable1=rbind(KarlukTable1,cbind(KarlukShortDates['KarlukSmolt2013.2'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2013.2']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                        formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.2']][2,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.2']][2,4]*100,format="f",digits=1),
                                        formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.2']][2,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.2']][2,6],format="f",digits=2),
                                        formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.2']][2,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.2']][2,2]*100,format="f",digits=1),'',
                                        formatC(x=KarlukAbundances[['KarlukSmolt2013.2']][2,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.2']][2,4],format="f",digits=0,big.mark=","),
                                        formatC(x=KarlukAbundances[['KarlukSmolt2013.2']][2,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.2']][2,1],
                                                                                                                                                            format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.2']][2,2],format="f",digits=0,big.mark=",")))
## Middle Total Outmigration row
KarlukTable1=rbind(KarlukTable1,cbind('','','','','','','','','','','','Middle Total','',formatC(x=KarlukOutmigration['KarlukSmolt2013.2'],format="f",digits=0,big.mark=","),''))
## Late KarlukEarly
KarlukTable1=rbind(KarlukTable1,cbind('Late',paste('n=',KarlukFinalSampleSizes['KarlukSmolt2013.3'],sep=''),KarlukGroups[1],
                                        formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.3']][1,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.3']][1,4]*100,format="f",digits=1),
                                        formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.3']][1,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.3']][1,6],format="f",digits=2),
                                        formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.3']][1,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.3']][1,2]*100,format="f",digits=1),'',
                                        formatC(x=KarlukAbundances[['KarlukSmolt2013.3']][1,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.3']][1,4],format="f",digits=0,big.mark=","),
                                        formatC(x=KarlukAbundances[['KarlukSmolt2013.3']][1,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.3']][1,1],
                                                                                                                                                            format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.3']][1,2],format="f",digits=0,big.mark=",")))
## Late KarlukLate
KarlukTable1=rbind(KarlukTable1,cbind(KarlukShortDates['KarlukSmolt2013.3'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2013.3']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                        formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.3']][2,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.3']][2,4]*100,format="f",digits=1),
                                        formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.3']][2,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.3']][2,6],format="f",digits=2),
                                        formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.3']][2,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2013.3']][2,2]*100,format="f",digits=1),'',
                                        formatC(x=KarlukAbundances[['KarlukSmolt2013.3']][2,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.3']][2,4],format="f",digits=0,big.mark=","),
                                        formatC(x=KarlukAbundances[['KarlukSmolt2013.3']][2,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.3']][2,1],
                                                                                                                                                            format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.3']][2,2],format="f",digits=0,big.mark=",")))
## Late Total Outmigration row
KarlukTable1=rbind(KarlukTable1,cbind('','','','','','','','','','','','Late Total','',formatC(x=KarlukOutmigration['KarlukSmolt2013.3'],format="f",digits=0,big.mark=","),''))
## 2013 KarlukEarly
KarlukTable1=rbind(KarlukTable1,cbind('2013',paste('n=',sum(KarlukFinalSampleSizes['KarlukSmolt2013.1'],KarlukFinalSampleSizes['KarlukSmolt2013.2'],KarlukFinalSampleSizes['KarlukSmolt2013.3']),sep=''),KarlukGroups[1],
                                        formatC(x=Karluk2013StratifiedEstimates$Summary[1,4]*100,format="f",digits=1),formatC(x=Karluk2013StratifiedEstimates$Summary[1,3]*100,format="f",digits=1),
                                        formatC(x=Karluk2013StratifiedEstimates$Summary[1,5]*100,format="f",digits=1),formatC(x=sum(Karluk2013StratifiedEstimates$Output[Karluk2013StratifiedEstimates$Output[,1]==0])/length(Karluk2013StratifiedEstimates$Output[,1]),format="f",digits=2),
                                        formatC(x=Karluk2013StratifiedEstimates$Summary[1,1]*100,format="f",digits=1),formatC(x=Karluk2013StratifiedEstimates$Summary[1,2]*100,format="f",digits=1),'',
                                        formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.1'],KarlukOutmigration['KarlukSmolt2013.2'],KarlukOutmigration['KarlukSmolt2013.3'])*Karluk2013StratifiedEstimates$Summary[1,4],format="f",digits=0,big.mark=","),
                                        formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.1'],KarlukOutmigration['KarlukSmolt2013.2'],KarlukOutmigration['KarlukSmolt2013.3'])*Karluk2013StratifiedEstimates$Summary[1,3],format="f",digits=0,big.mark=","),
                                        formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.1'],KarlukOutmigration['KarlukSmolt2013.2'],KarlukOutmigration['KarlukSmolt2013.3'])*Karluk2013StratifiedEstimates$Summary[1,5],format="f",digits=0,big.mark=","),
                                        formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.1'],KarlukOutmigration['KarlukSmolt2013.2'],KarlukOutmigration['KarlukSmolt2013.3'])*Karluk2013StratifiedEstimates$Summary[1,1],format="f",digits=0,big.mark=","),
                                        formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.1'],KarlukOutmigration['KarlukSmolt2013.2'],KarlukOutmigration['KarlukSmolt2013.3'])*Karluk2013StratifiedEstimates$Summary[1,2],format="f",digits=0,big.mark=",")))
## 2013 Karluk
KarlukTable1=rbind(KarlukTable1,cbind(AnnualShortDates['2013'],'',KarlukGroups[2],
                                        formatC(x=Karluk2013StratifiedEstimates$Summary[2,4]*100,format="f",digits=1),formatC(x=Karluk2013StratifiedEstimates$Summary[2,3]*100,format="f",digits=1),
                                        formatC(x=Karluk2013StratifiedEstimates$Summary[2,5]*100,format="f",digits=1),formatC(x=sum(Karluk2013StratifiedEstimates$Output[Karluk2013StratifiedEstimates$Output[,2]==0])/length(Karluk2013StratifiedEstimates$Output[,2]),format="f",digits=2),
                                        formatC(x=Karluk2013StratifiedEstimates$Summary[2,1]*100,format="f",digits=1),formatC(x=Karluk2013StratifiedEstimates$Summary[2,2]*100,format="f",digits=1),'',
                                        formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.1'],KarlukOutmigration['KarlukSmolt2013.2'],KarlukOutmigration['KarlukSmolt2013.3'])*Karluk2013StratifiedEstimates$Summary[2,4],format="f",digits=0,big.mark=","),
                                        formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.1'],KarlukOutmigration['KarlukSmolt2013.2'],KarlukOutmigration['KarlukSmolt2013.3'])*Karluk2013StratifiedEstimates$Summary[2,3],format="f",digits=0,big.mark=","),
                                        formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.1'],KarlukOutmigration['KarlukSmolt2013.2'],KarlukOutmigration['KarlukSmolt2013.3'])*Karluk2013StratifiedEstimates$Summary[2,5],format="f",digits=0,big.mark=","),
                                        formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.1'],KarlukOutmigration['KarlukSmolt2013.2'],KarlukOutmigration['KarlukSmolt2013.3'])*Karluk2013StratifiedEstimates$Summary[2,1],format="f",digits=0,big.mark=","),
                                        formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.1'],KarlukOutmigration['KarlukSmolt2013.2'],KarlukOutmigration['KarlukSmolt2013.3'])*Karluk2013StratifiedEstimates$Summary[2,2],format="f",digits=0,big.mark=",")))
## 2013 Total Outmigration row
KarlukTable1=rbind(KarlukTable1,cbind('','','','','','','','','','','','2013 Total','',formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.1'],KarlukOutmigration['KarlukSmolt2013.2'],KarlukOutmigration['KarlukSmolt2013.3']),format="f",digits=0,big.mark=","),''),Disclaimer)

### Write table
write.xlsx(x=as.data.frame(KarlukTable1),file="V:/WORK/Sockeye/Kodiak/2013 2014 Karluk Smolt/Tables/KarlukSmolt2013-2014Tables.xlsx",col.names=F,row.names=F,append=F,sheetName="Karluk 2013 Smolt")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Table 2, 2014 smolt ####
## 17 rows by 15 columns for three strata/year and annual roll up.

KarlukCaption2 <- "Table X.-Estimates of stock composition and stock-specific outmigration for Karluk River sockeye salmon smolt by stratum, 2014.  Reporting group-specific stock composition estimates (%) include median, 90% credibility interval, the probability that the group estimate is equal to zero (P=0), mean and SD.  Stock-specific estimates of outmigration are based upon mark-recapture estimates and variances (CV) of outmigration for each stratum.  See text for details."
Disclaimer <- cbind("Note: Stock composition estimates may not sum to 100% and stock-specific outmigration estimates may not sum to the total outmigration due to rounding error.",'','','','','','','','','','','','','','')

KarlukTable2=cbind(KarlukCaption2,'','','','','','','','','','','','','','')
KarlukTable2=rbind(KarlukTable2,cbind("Stratum","",'','Composition (%)','','','','','','',"Outmigration (number of fish)",'','','',''))
KarlukTable2=rbind(KarlukTable2,cbind('Period','Sample Size','Reporting','',"90% CI",'','','','','','',"90% CI",'','',''))
KarlukTable2=rbind(KarlukTable2,cbind("Dates",'CV','Group',"Median","5%","95%","P=0","Mean","SD",'',"Median","5%","95%","Mean","SD"))
## Early KarlukEarly
KarlukTable2=rbind(KarlukTable2,cbind('Early',paste('n=',KarlukFinalSampleSizes['KarlukSmolt2014.1'],sep=''),KarlukGroups[1],
                                      formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.1']][1,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.1']][1,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.1']][1,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.1']][1,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.1']][1,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.1']][1,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.1']][1,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.1']][1,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.1']][1,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.1']][1,1],
                                                                                                                                       format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.1']][1,2],format="f",digits=0,big.mark=",")))
## Early KarlukLate
KarlukTable2=rbind(KarlukTable2,cbind(KarlukShortDates['KarlukSmolt2014.1'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2014.1']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                      formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.1']][2,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.1']][2,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.1']][2,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.1']][2,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.1']][2,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.1']][2,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.1']][2,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.1']][2,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.1']][2,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.1']][2,1],
                                                                                                                                       format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.1']][2,2],format="f",digits=0,big.mark=",")))
## Early Total Outmigration row
KarlukTable2=rbind(KarlukTable2,cbind('','','','','','','','','','','','Early Total','',formatC(x=KarlukOutmigration['KarlukSmolt2014.1'],format="f",digits=0,big.mark=","),''))
## Middle KarlukEarly
KarlukTable2=rbind(KarlukTable2,cbind('Middle',paste('n=',KarlukFinalSampleSizes['KarlukSmolt2014.2'],sep=''),KarlukGroups[1],
                                      formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.2']][1,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.2']][1,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.2']][1,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.2']][1,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.2']][1,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.2']][1,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.2']][1,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.2']][1,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.2']][1,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.2']][1,1],
                                                                                                                                       format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.2']][1,2],format="f",digits=0,big.mark=",")))
## Middle KarlukLate
KarlukTable2=rbind(KarlukTable2,cbind(KarlukShortDates['KarlukSmolt2014.2'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2014.2']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                      formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.2']][2,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.2']][2,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.2']][2,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.2']][2,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.2']][2,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.2']][2,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.2']][2,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.2']][2,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.2']][2,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.2']][2,1],
                                                                                                                                       format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.2']][2,2],format="f",digits=0,big.mark=",")))
## Middle Total Outmigration row
KarlukTable2=rbind(KarlukTable2,cbind('','','','','','','','','','','','Middle Total','',formatC(x=KarlukOutmigration['KarlukSmolt2014.2'],format="f",digits=0,big.mark=","),''))
## Late KarlukEarly
KarlukTable2=rbind(KarlukTable2,cbind('Late',paste('n=',KarlukFinalSampleSizes['KarlukSmolt2014.3'],sep=''),KarlukGroups[1],
                                      formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.3']][1,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.3']][1,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.3']][1,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.3']][1,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.3']][1,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.3']][1,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.3']][1,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.3']][1,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.3']][1,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.3']][1,1],
                                                                                                                                       format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.3']][1,2],format="f",digits=0,big.mark=",")))
## Late KarlukLate
KarlukTable2=rbind(KarlukTable2,cbind(KarlukShortDates['KarlukSmolt2014.3'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2014.3']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                      formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.3']][2,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.3']][2,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.3']][2,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.3']][2,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.3']][2,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['KarlukSmolt2014.3']][2,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.3']][2,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.3']][2,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.3']][2,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.3']][2,1],
                                                                                                                                       format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.3']][2,2],format="f",digits=0,big.mark=",")))
## Late Total Outmigration row
KarlukTable2=rbind(KarlukTable2,cbind('','','','','','','','','','','','Late Total','',formatC(x=KarlukOutmigration['KarlukSmolt2014.3'],format="f",digits=0,big.mark=","),''))
## 2014 KarlukEarly
KarlukTable2=rbind(KarlukTable2,cbind('2014',paste('n=',sum(KarlukFinalSampleSizes['KarlukSmolt2014.1'],KarlukFinalSampleSizes['KarlukSmolt2014.2'],KarlukFinalSampleSizes['KarlukSmolt2014.3']),sep=''),KarlukGroups[1],
                                      formatC(x=Karluk2014StratifiedEstimates$Summary[1,4]*100,format="f",digits=1),formatC(x=Karluk2014StratifiedEstimates$Summary[1,3]*100,format="f",digits=1),
                                      formatC(x=Karluk2014StratifiedEstimates$Summary[1,5]*100,format="f",digits=1),formatC(x=sum(Karluk2014StratifiedEstimates$Output[Karluk2014StratifiedEstimates$Output[,1]==0])/length(Karluk2014StratifiedEstimates$Output[,1]),format="f",digits=2),
                                      formatC(x=Karluk2014StratifiedEstimates$Summary[1,1]*100,format="f",digits=1),formatC(x=Karluk2014StratifiedEstimates$Summary[1,2]*100,format="f",digits=1),'',
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.1'],KarlukOutmigration['KarlukSmolt2014.2'],KarlukOutmigration['KarlukSmolt2014.3'])*Karluk2014StratifiedEstimates$Summary[1,4],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.1'],KarlukOutmigration['KarlukSmolt2014.2'],KarlukOutmigration['KarlukSmolt2014.3'])*Karluk2014StratifiedEstimates$Summary[1,3],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.1'],KarlukOutmigration['KarlukSmolt2014.2'],KarlukOutmigration['KarlukSmolt2014.3'])*Karluk2014StratifiedEstimates$Summary[1,5],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.1'],KarlukOutmigration['KarlukSmolt2014.2'],KarlukOutmigration['KarlukSmolt2014.3'])*Karluk2014StratifiedEstimates$Summary[1,1],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.1'],KarlukOutmigration['KarlukSmolt2014.2'],KarlukOutmigration['KarlukSmolt2014.3'])*Karluk2014StratifiedEstimates$Summary[1,2],format="f",digits=0,big.mark=",")))
## 2014 Karluk
KarlukTable2=rbind(KarlukTable2,cbind(AnnualShortDates['2014'],'',KarlukGroups[2],
                                      formatC(x=Karluk2014StratifiedEstimates$Summary[2,4]*100,format="f",digits=1),formatC(x=Karluk2014StratifiedEstimates$Summary[2,3]*100,format="f",digits=1),
                                      formatC(x=Karluk2014StratifiedEstimates$Summary[2,5]*100,format="f",digits=1),formatC(x=sum(Karluk2014StratifiedEstimates$Output[Karluk2014StratifiedEstimates$Output[,2]==0])/length(Karluk2014StratifiedEstimates$Output[,2]),format="f",digits=2),
                                      formatC(x=Karluk2014StratifiedEstimates$Summary[2,1]*100,format="f",digits=1),formatC(x=Karluk2014StratifiedEstimates$Summary[2,2]*100,format="f",digits=1),'',
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.1'],KarlukOutmigration['KarlukSmolt2014.2'],KarlukOutmigration['KarlukSmolt2014.3'])*Karluk2014StratifiedEstimates$Summary[2,4],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.1'],KarlukOutmigration['KarlukSmolt2014.2'],KarlukOutmigration['KarlukSmolt2014.3'])*Karluk2014StratifiedEstimates$Summary[2,3],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.1'],KarlukOutmigration['KarlukSmolt2014.2'],KarlukOutmigration['KarlukSmolt2014.3'])*Karluk2014StratifiedEstimates$Summary[2,5],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.1'],KarlukOutmigration['KarlukSmolt2014.2'],KarlukOutmigration['KarlukSmolt2014.3'])*Karluk2014StratifiedEstimates$Summary[2,1],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.1'],KarlukOutmigration['KarlukSmolt2014.2'],KarlukOutmigration['KarlukSmolt2014.3'])*Karluk2014StratifiedEstimates$Summary[2,2],format="f",digits=0,big.mark=",")))
## 2014 Total Outmigration row
KarlukTable2=rbind(KarlukTable2,cbind('','','','','','','','','','','','2014 Total','',formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.1'],KarlukOutmigration['KarlukSmolt2014.2'],KarlukOutmigration['KarlukSmolt2014.3']),format="f",digits=0,big.mark=","),''),Disclaimer)

### Write table
write.xlsx(x=as.data.frame(KarlukTable2),file="V:/WORK/Sockeye/Kodiak/2013 2014 Karluk Smolt/Tables/KarlukSmolt2013-2014Tables.xlsx",col.names=F,row.names=F,append=TRUE,sheetName="Karluk 2014 Smolt")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Table 3, 2014 Weir smolt ####
## 7 rows by 9 columns for three strata/year and annual roll up.

KarlukCaption3 <- "Table X.-Estimates of stock composition for Karluk River sockeye salmon smolt collected at the Weir, 2014.  Reporting group-specific stock composition estimates (%) include median, 90% credibility interval, the probability that the group estimate is equal to zero (P=0), mean and SD.See text for details."
Disclaimer <- cbind("Note: Stock composition estimates may not sum to 100% and stock-specific outmigration estimates may not sum to the total outmigration due to rounding error.",'','','','','','','','')

KarlukTable3=cbind(KarlukCaption3,'','','','','','','','')
KarlukTable3=rbind(KarlukTable3,cbind("Stratum",'','','Composition (%)','','','','',''))
KarlukTable3=rbind(KarlukTable3,cbind('Period','Sample Size','Reporting','',"90% CI",'','','',''))
KarlukTable3=rbind(KarlukTable3,cbind("Dates",'CV','Group',"Median","5%","95%","P=0","Mean","SD"))
## Early KarlukEarly
KarlukTable3=rbind(KarlukTable3,cbind('Early',paste('n=',KarlukFinalSampleSizes['SKARLW14s'],sep=''),KarlukGroups[1],
                                      formatC(x=KarlukSmolt_Estimates[['SKARLW14s']][1,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['SKARLW14s']][1,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Estimates[['SKARLW14s']][1,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['SKARLW14s']][1,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Estimates[['SKARLW14s']][1,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['SKARLW14s']][1,2]*100,format="f",digits=1)))
## Early KarlukLate
KarlukTable3=rbind(KarlukTable3,cbind(KarlukShortDates['SKARLW14s'],paste('CV=',formatC(x=KarlukCVs['SKARLW14s']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                      formatC(x=KarlukSmolt_Estimates[['SKARLW14s']][2,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['SKARLW14s']][2,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Estimates[['SKARLW14s']][2,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['SKARLW14s']][2,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Estimates[['SKARLW14s']][2,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Estimates[['SKARLW14s']][2,2]*100,format="f",digits=1)))
KarlukTable3=rbind(KarlukTable3,Disclaimer)

### Write table
write.xlsx(x=as.data.frame(KarlukTable3),file="V:/WORK/Sockeye/Kodiak/2013 2014 Karluk Smolt/Tables/KarlukSmolt2013-2014Tables.xlsx",col.names=F,row.names=F,append=TRUE,sheetName="Karluk 2014 Weir Smolt")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Table 4, 2013 smolt by age ####
## 15 rows by 15 columns for three strata/year and annual roll up.

KarlukCaption4 <- "Table X.-Estimates of stock composition and stock-specific outmigration for Karluk River sockeye salmon smolt by age, 2013.  Reporting group-specific stock composition estimates (%) include median, 90% credibility interval, the probability that the group estimate is equal to zero (P=0), mean and SD.  Stock-specific estimates of outmigration are based upon mark-recapture estimates and variances (CV) of outmigration for each stratum.  See text for details."
Disclaimer <- cbind("Note: Stock composition estimates may not sum to 100% and stock-specific outmigration estimates may not sum to the total outmigration due to rounding error.",'','','','','','','','','','','','','','')

KarlukTable4=cbind(KarlukCaption4,'','','','','','','','','','','','','','')
KarlukTable4=rbind(KarlukTable4,cbind("Stratum","",'','Composition (%)','','','','','','',"Outmigration (number of fish)",'','','',''))
KarlukTable4=rbind(KarlukTable4,cbind('Age','Sample Size','Reporting','',"90% CI",'','','','','','',"90% CI",'','',''))
KarlukTable4=rbind(KarlukTable4,cbind("Dates",'CV','Group',"Median","5%","95%","P=0","Mean","SD",'',"Median","5%","95%","Mean","SD"))
## Age 1 KarlukEarly
KarlukTable4=rbind(KarlukTable4,cbind('Age 1',paste('n=',KarlukFinalSampleSizes['KarlukSmolt2013.Age1'],sep=''),KarlukGroups[1],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age1']][1,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age1']][1,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age1']][1,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age1']][1,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age1']][1,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age1']][1,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age1']][1,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age1']][1,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age1']][1,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age1']][1,1],
                                                                                                                                       format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age1']][1,2],format="f",digits=0,big.mark=",")))
## Age 1 KarlukLate
KarlukTable4=rbind(KarlukTable4,cbind(KarlukShortDates['KarlukSmolt2013.Age1'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2013.Age1']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age1']][2,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age1']][2,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age1']][2,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age1']][2,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age1']][2,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age1']][2,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age1']][2,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age1']][2,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age1']][2,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age1']][2,1],
                                                                                                                                       format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age1']][2,2],format="f",digits=0,big.mark=",")))
## Age 1 Total Outmigration row
KarlukTable4=rbind(KarlukTable4,cbind('','','','','','','','','','','','Age 1 Total','',formatC(x=KarlukOutmigration['KarlukSmolt2013.Age1'],format="f",digits=0,big.mark=","),''))
## Age 2 KarlukEarly
KarlukTable4=rbind(KarlukTable4,cbind('Age 2',paste('n=',sum(KarlukFinalSampleSizes['KarlukSmolt2013.Age2.1'],KarlukFinalSampleSizes['KarlukSmolt2013.Age2.2'],KarlukFinalSampleSizes['KarlukSmolt2013.Age2.3']),sep=''),KarlukGroups[1],
                                      formatC(x=Karluk2013Age2StratifiedEstimates$Summary[1,4]*100,format="f",digits=1),formatC(x=Karluk2013Age2StratifiedEstimates$Summary[1,3]*100,format="f",digits=1),
                                      formatC(x=Karluk2013Age2StratifiedEstimates$Summary[1,5]*100,format="f",digits=1),formatC(x=Karluk2013Age2StratifiedEstimates$Summary[1,6],format="f",digits=2),
                                      formatC(x=Karluk2013Age2StratifiedEstimates$Summary[1,1]*100,format="f",digits=1),formatC(x=Karluk2013Age2StratifiedEstimates$Summary[1,2]*100,format="f",digits=1),'',
                                      formatC(x=Karluk2013Age2StratifiedAbundances[1,4],format="f",digits=0,big.mark=","),formatC(x=Karluk2013Age2StratifiedAbundances[1,3],format="f",digits=0,big.mark=","),
                                      formatC(x=Karluk2013Age2StratifiedAbundances[1,5],format="f",digits=0,big.mark=","),formatC(x=Karluk2013Age2StratifiedAbundances[1,1],
                                                                                                                                       format="f",digits=0,big.mark=","),formatC(x=Karluk2013Age2StratifiedAbundances[1,2],format="f",digits=0,big.mark=",")))
## Age 2 KarlukLate
KarlukTable4=rbind(KarlukTable4,cbind(KarlukShortDates['KarlukSmolt2013.Age1'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2013.2']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                      formatC(x=Karluk2013Age2StratifiedEstimates$Summary[2,4]*100,format="f",digits=1),formatC(x=Karluk2013Age2StratifiedEstimates$Summary[2,3]*100,format="f",digits=1),
                                      formatC(x=Karluk2013Age2StratifiedEstimates$Summary[2,5]*100,format="f",digits=1),formatC(x=Karluk2013Age2StratifiedEstimates$Summary[2,6],format="f",digits=2),
                                      formatC(x=Karluk2013Age2StratifiedEstimates$Summary[2,1]*100,format="f",digits=1),formatC(x=Karluk2013Age2StratifiedEstimates$Summary[2,2]*100,format="f",digits=1),'',
                                      formatC(x=Karluk2013Age2StratifiedAbundances[2,4],format="f",digits=0,big.mark=","),formatC(x=Karluk2013Age2StratifiedAbundances[2,3],format="f",digits=0,big.mark=","),
                                      formatC(x=Karluk2013Age2StratifiedAbundances[2,5],format="f",digits=0,big.mark=","),formatC(x=Karluk2013Age2StratifiedAbundances[2,1],
                                                                                                                                       format="f",digits=0,big.mark=","),formatC(x=Karluk2013Age2StratifiedAbundances[2,2],format="f",digits=0,big.mark=",")))
## Age 2 Total Outmigration row
KarlukTable4=rbind(KarlukTable4,cbind('','','','','','','','','','','','Age 2 Total','',formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.Age2.1'],KarlukOutmigration['KarlukSmolt2013.Age2.2'],KarlukOutmigration['KarlukSmolt2013.Age2.3']),format="f",digits=0,big.mark=","),''))
## Age 3 KarlukEarly
KarlukTable4=rbind(KarlukTable4,cbind('Age 3',paste('n=',KarlukFinalSampleSizes['KarlukSmolt2013.Age3'],sep=''),KarlukGroups[1],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age3']][1,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age3']][1,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age3']][1,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age3']][1,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age3']][1,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age3']][1,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age3']][1,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age3']][1,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age3']][1,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age3']][1,1],
                                                                                                                                       format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age3']][1,2],format="f",digits=0,big.mark=",")))
## Age 3 KarlukLate
KarlukTable4=rbind(KarlukTable4,cbind(KarlukShortDates['KarlukSmolt2013.Age3'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2013.Age3']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age3']][2,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age3']][2,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age3']][2,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age3']][2,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age3']][2,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age3']][2,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age3']][2,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age3']][2,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age3']][2,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age3']][2,1],
                                                                                                                                       format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age3']][2,2],format="f",digits=0,big.mark=",")))
## Age 3 Total Outmigration row
KarlukTable4=rbind(KarlukTable4,cbind('','','','','','','','','','','','Age 3 Total','',formatC(x=KarlukOutmigration['KarlukSmolt2013.Age3'],format="f",digits=0,big.mark=","),''))

## 2013 Total Outmigration row
KarlukTable4=rbind(KarlukTable4,cbind('','','','','','','','','','','','2013 Total','',formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.Age1'],KarlukOutmigration['KarlukSmolt2013.Age2.1'],KarlukOutmigration['KarlukSmolt2013.Age2.2'],KarlukOutmigration['KarlukSmolt2013.Age2.3'],KarlukOutmigration['KarlukSmolt2013.Age3']),format="f",digits=0,big.mark=","),''),Disclaimer)

### Write table
write.xlsx(x=as.data.frame(KarlukTable4),file="V:/WORK/Sockeye/Kodiak/2013 2014 Karluk Smolt/Tables/KarlukSmolt2013-2014Tables.xlsx",col.names=F,row.names=F,append=TRUE,sheetName="Karluk 2013 Smolt by Age")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Table 5, 2014 smolt by age ####
## 12 rows by 15 columns for three strata/year and annual roll up.

KarlukCaption5 <- "Table X.-Estimates of stock composition and stock-specific outmigration for Karluk River sockeye salmon smolt by age, 2014.  Reporting group-specific stock composition estimates (%) include median, 90% credibility interval, the probability that the group estimate is equal to zero (P=0), mean and SD.  Stock-specific estimates of outmigration are based upon mark-recapture estimates and variances (CV) of outmigration for each stratum.  See text for details."
Disclaimer <- cbind("Note: Stock composition estimates may not sum to 100% and stock-specific outmigration estimates may not sum to the total outmigration due to rounding error.",'','','','','','','','','','','','','','')

KarlukTable5=cbind(KarlukCaption5,'','','','','','','','','','','','','','')
KarlukTable5=rbind(KarlukTable5,cbind("Stratum","",'','Composition (%)','','','','','','',"Outmigration (number of fish)",'','','',''))
KarlukTable5=rbind(KarlukTable5,cbind('Age','Sample Size','Reporting','',"90% CI",'','','','','','',"90% CI",'','',''))
KarlukTable5=rbind(KarlukTable5,cbind("Dates",'CV','Group',"Median","5%","95%","P=0","Mean","SD",'',"Median","5%","95%","Mean","SD"))
## Age 1 KarlukEarly
KarlukTable5=rbind(KarlukTable5,cbind('Age 1',paste('n=',sum(KarlukFinalSampleSizes['KarlukSmolt2014.Age1.12'],KarlukFinalSampleSizes['KarlukSmolt2014.Age1.3']),sep=''),KarlukGroups[1],
                                      formatC(x=Karluk2014Age1StratifiedEstimates$Summary[1,4]*100,format="f",digits=1),formatC(x=Karluk2014Age1StratifiedEstimates$Summary[1,3]*100,format="f",digits=1),
                                      formatC(x=Karluk2014Age1StratifiedEstimates$Summary[1,5]*100,format="f",digits=1),formatC(x=Karluk2014Age1StratifiedEstimates$Summary[1,6],format="f",digits=2),
                                      formatC(x=Karluk2014Age1StratifiedEstimates$Summary[1,1]*100,format="f",digits=1),formatC(x=Karluk2014Age1StratifiedEstimates$Summary[1,2]*100,format="f",digits=1),'',
                                      formatC(x=Karluk2014Age1StratifiedAbundances[1,4],format="f",digits=0,big.mark=","),formatC(x=Karluk2014Age1StratifiedAbundances[1,3],format="f",digits=0,big.mark=","),
                                      formatC(x=Karluk2014Age1StratifiedAbundances[1,5],format="f",digits=0,big.mark=","),formatC(x=Karluk2014Age1StratifiedAbundances[1,1],
                                                                                                                                  format="f",digits=0,big.mark=","),formatC(x=Karluk2014Age1StratifiedAbundances[1,2],format="f",digits=0,big.mark=",")))
## Age 1 KarlukLate
KarlukTable5=rbind(KarlukTable5,cbind(AnnualShortDates['2014'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2014.2']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                      formatC(x=Karluk2014Age1StratifiedEstimates$Summary[2,4]*100,format="f",digits=1),formatC(x=Karluk2014Age1StratifiedEstimates$Summary[2,3]*100,format="f",digits=1),
                                      formatC(x=Karluk2014Age1StratifiedEstimates$Summary[2,5]*100,format="f",digits=1),formatC(x=Karluk2014Age1StratifiedEstimates$Summary[2,6],format="f",digits=2),
                                      formatC(x=Karluk2014Age1StratifiedEstimates$Summary[2,1]*100,format="f",digits=1),formatC(x=Karluk2014Age1StratifiedEstimates$Summary[2,2]*100,format="f",digits=1),'',
                                      formatC(x=Karluk2014Age1StratifiedAbundances[2,4],format="f",digits=0,big.mark=","),formatC(x=Karluk2014Age1StratifiedAbundances[2,3],format="f",digits=0,big.mark=","),
                                      formatC(x=Karluk2014Age1StratifiedAbundances[2,5],format="f",digits=0,big.mark=","),formatC(x=Karluk2014Age1StratifiedAbundances[2,1],
                                                                                                                                  format="f",digits=0,big.mark=","),formatC(x=Karluk2014Age1StratifiedAbundances[2,2],format="f",digits=0,big.mark=",")))
## Age 1 Total Outmigration row
KarlukTable5=rbind(KarlukTable5,cbind('','','','','','','','','','','','Age 1 Total','',formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age1.12'],KarlukOutmigration['KarlukSmolt2014.Age1.3']),format="f",digits=0,big.mark=","),''))
## Age 2 KarlukEarly
KarlukTable5=rbind(KarlukTable5,cbind('Age 2',paste('n=',sum(KarlukFinalSampleSizes['KarlukSmolt2014.Age2.1'],KarlukFinalSampleSizes['KarlukSmolt2014.Age2.2'],KarlukFinalSampleSizes['KarlukSmolt2014.Age2.3']),sep=''),KarlukGroups[1],
                                      formatC(x=Karluk2014Age2StratifiedEstimates$Summary[1,4]*100,format="f",digits=1),formatC(x=Karluk2014Age2StratifiedEstimates$Summary[1,3]*100,format="f",digits=1),
                                      formatC(x=Karluk2014Age2StratifiedEstimates$Summary[1,5]*100,format="f",digits=1),formatC(x=Karluk2014Age2StratifiedEstimates$Summary[1,6],format="f",digits=2),
                                      formatC(x=Karluk2014Age2StratifiedEstimates$Summary[1,1]*100,format="f",digits=1),formatC(x=Karluk2014Age2StratifiedEstimates$Summary[1,2]*100,format="f",digits=1),'',
                                      formatC(x=Karluk2014Age2StratifiedAbundances[1,4],format="f",digits=0,big.mark=","),formatC(x=Karluk2014Age2StratifiedAbundances[1,3],format="f",digits=0,big.mark=","),
                                      formatC(x=Karluk2014Age2StratifiedAbundances[1,5],format="f",digits=0,big.mark=","),formatC(x=Karluk2014Age2StratifiedAbundances[1,1],
                                                                                                                                  format="f",digits=0,big.mark=","),formatC(x=Karluk2014Age2StratifiedAbundances[1,2],format="f",digits=0,big.mark=",")))
## Age 2 KarlukLate
KarlukTable5=rbind(KarlukTable5,cbind(AnnualShortDates['2014'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2014.2']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                      formatC(x=Karluk2014Age2StratifiedEstimates$Summary[2,4]*100,format="f",digits=1),formatC(x=Karluk2014Age2StratifiedEstimates$Summary[2,3]*100,format="f",digits=1),
                                      formatC(x=Karluk2014Age2StratifiedEstimates$Summary[2,5]*100,format="f",digits=1),formatC(x=Karluk2014Age2StratifiedEstimates$Summary[2,6],format="f",digits=2),
                                      formatC(x=Karluk2014Age2StratifiedEstimates$Summary[2,1]*100,format="f",digits=1),formatC(x=Karluk2014Age2StratifiedEstimates$Summary[2,2]*100,format="f",digits=1),'',
                                      formatC(x=Karluk2014Age2StratifiedAbundances[2,4],format="f",digits=0,big.mark=","),formatC(x=Karluk2014Age2StratifiedAbundances[2,3],format="f",digits=0,big.mark=","),
                                      formatC(x=Karluk2014Age2StratifiedAbundances[2,5],format="f",digits=0,big.mark=","),formatC(x=Karluk2014Age2StratifiedAbundances[2,1],
                                                                                                                                  format="f",digits=0,big.mark=","),formatC(x=Karluk2014Age2StratifiedAbundances[2,2],format="f",digits=0,big.mark=",")))
## Age 2 Total Outmigration row
KarlukTable5=rbind(KarlukTable5,cbind('','','','','','','','','','','','Age 2 Total','',formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age2.1'],KarlukOutmigration['KarlukSmolt2014.Age2.2'],KarlukOutmigration['KarlukSmolt2014.Age2.3']),format="f",digits=0,big.mark=","),''))

## 2014 Total Outmigration row
KarlukTable5=rbind(KarlukTable5,cbind('','','','','','','','','','','','2014 Total','',formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age2.1'],KarlukOutmigration['KarlukSmolt2014.Age2.2'],KarlukOutmigration['KarlukSmolt2014.Age2.3'],KarlukOutmigration['KarlukSmolt2014.Age1.12'],KarlukOutmigration['KarlukSmolt2014.Age1.3']),format="f",digits=0,big.mark=","),''),Disclaimer)

### Write table
write.xlsx(x=as.data.frame(KarlukTable5),file="V:/WORK/Sockeye/Kodiak/2013 2014 Karluk Smolt/Tables/KarlukSmolt2013-2014Tables.xlsx",col.names=F,row.names=F,append=TRUE,sheetName="Karluk 2014 Smolt by Age")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Table 6, 2013 smolt age 2 ####
## 17 rows by 15 columns for three strata/year and annual roll up.

KarlukCaption6 <- "Table X.-Estimates of stock composition and stock-specific outmigration for Karluk River Age 2 sockeye salmon smolt by stratum, 2013.  Reporting group-specific stock composition estimates (%) include median, 90% credibility interval, the probability that the group estimate is equal to zero (P=0), mean and SD.  Stock-specific estimates of outmigration are based upon mark-recapture estimates and variances (CV) of outmigration for each stratum.  See text for details."
Disclaimer <- cbind("Note: Stock composition estimates may not sum to 100% and stock-specific outmigration estimates may not sum to the total outmigration due to rounding error.",'','','','','','','','','','','','','','')

KarlukTable6=cbind(KarlukCaption6,'','','','','','','','','','','','','','')
KarlukTable6=rbind(KarlukTable6,cbind("Stratum","",'','Composition (%)','','','','','','',"Outmigration (number of fish)",'','','',''))
KarlukTable6=rbind(KarlukTable6,cbind('Period','Sample Size','Reporting','',"90% CI",'','','','','','',"90% CI",'','',''))
KarlukTable6=rbind(KarlukTable6,cbind("Dates",'CV','Group',"Median","5%","95%","P=0","Mean","SD",'',"Median","5%","95%","Mean","SD"))
## Early KarlukEarly
KarlukTable6=rbind(KarlukTable6,cbind('Early',paste('n=',KarlukFinalSampleSizes['KarlukSmolt2013.Age2.1'],sep=''),KarlukGroups[1],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.1']][1,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.1']][1,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.1']][1,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.1']][1,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.1']][1,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.1']][1,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.1']][1,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.1']][1,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.1']][1,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.1']][1,1],
                                                                                                                                       format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.1']][1,2],format="f",digits=0,big.mark=",")))
## Early KarlukLate
KarlukTable6=rbind(KarlukTable6,cbind(KarlukShortDates['KarlukSmolt2013.Age2.1'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2013.Age2.1']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.1']][2,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.1']][2,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.1']][2,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.1']][2,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.1']][2,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.1']][2,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.1']][2,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.1']][2,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.1']][2,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.1']][2,1],
                                                                                                                                       format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.1']][2,2],format="f",digits=0,big.mark=",")))
## Early Total Outmigration row
KarlukTable6=rbind(KarlukTable6,cbind('','','','','','','','','','','','Early Total','',formatC(x=KarlukOutmigration['KarlukSmolt2013.Age2.1'],format="f",digits=0,big.mark=","),''))
## Middle KarlukEarly
KarlukTable6=rbind(KarlukTable6,cbind('Middle',paste('n=',KarlukFinalSampleSizes['KarlukSmolt2013.Age2.2'],sep=''),KarlukGroups[1],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.2']][1,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.2']][1,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.2']][1,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.2']][1,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.2']][1,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.2']][1,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.2']][1,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.2']][1,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.2']][1,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.2']][1,1],
                                                                                                                                       format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.2']][1,2],format="f",digits=0,big.mark=",")))
## Middle KarlukLate
KarlukTable6=rbind(KarlukTable6,cbind(KarlukShortDates['KarlukSmolt2013.Age2.2'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2013.Age2.2']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.2']][2,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.2']][2,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.2']][2,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.2']][2,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.2']][2,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.2']][2,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.2']][2,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.2']][2,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.2']][2,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.2']][2,1],
                                                                                                                                       format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.2']][2,2],format="f",digits=0,big.mark=",")))
## Middle Total Outmigration row
KarlukTable6=rbind(KarlukTable6,cbind('','','','','','','','','','','','Middle Total','',formatC(x=KarlukOutmigration['KarlukSmolt2013.Age2.2'],format="f",digits=0,big.mark=","),''))
## Late KarlukEarly
KarlukTable6=rbind(KarlukTable6,cbind('Late',paste('n=',KarlukFinalSampleSizes['KarlukSmolt2013.Age2.3'],sep=''),KarlukGroups[1],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.3']][1,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.3']][1,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.3']][1,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.3']][1,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.3']][1,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.3']][1,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.3']][1,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.3']][1,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.3']][1,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.3']][1,1],
                                                                                                                                       format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.3']][1,2],format="f",digits=0,big.mark=",")))
## Late KarlukLate
KarlukTable6=rbind(KarlukTable6,cbind(KarlukShortDates['KarlukSmolt2013.Age2.3'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2013.Age2.3']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.3']][2,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.3']][2,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.3']][2,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.3']][2,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.3']][2,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2013.Age2.3']][2,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.3']][2,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.3']][2,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.3']][2,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.3']][2,1],
                                                                                                                                       format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2013.Age2.3']][2,2],format="f",digits=0,big.mark=",")))
## Late Total Outmigration row
KarlukTable6=rbind(KarlukTable6,cbind('','','','','','','','','','','','Late Total','',formatC(x=KarlukOutmigration['KarlukSmolt2013.Age2.3'],format="f",digits=0,big.mark=","),''))
## 2013 KarlukEarly
KarlukTable6=rbind(KarlukTable6,cbind('2013',paste('n=',sum(KarlukFinalSampleSizes['KarlukSmolt2013.Age2.1'],KarlukFinalSampleSizes['KarlukSmolt2013.Age2.2'],KarlukFinalSampleSizes['KarlukSmolt2013.Age2.3']),sep=''),KarlukGroups[1],
                                      formatC(x=Karluk2013Age2StratifiedEstimates$Summary[1,4]*100,format="f",digits=1),formatC(x=Karluk2013Age2StratifiedEstimates$Summary[1,3]*100,format="f",digits=1),
                                      formatC(x=Karluk2013Age2StratifiedEstimates$Summary[1,5]*100,format="f",digits=1),formatC(x=sum(Karluk2013Age2StratifiedEstimates$Output[Karluk2013Age2StratifiedEstimates$Output[,1]==0])/length(Karluk2013Age2StratifiedEstimates$Output[,1]),format="f",digits=2),
                                      formatC(x=Karluk2013Age2StratifiedEstimates$Summary[1,1]*100,format="f",digits=1),formatC(x=Karluk2013Age2StratifiedEstimates$Summary[1,2]*100,format="f",digits=1),'',
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.Age2.1'],KarlukOutmigration['KarlukSmolt2013.Age2.2'],KarlukOutmigration['KarlukSmolt2013.Age2.3'])*Karluk2013Age2StratifiedEstimates$Summary[1,4],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.Age2.1'],KarlukOutmigration['KarlukSmolt2013.Age2.2'],KarlukOutmigration['KarlukSmolt2013.Age2.3'])*Karluk2013Age2StratifiedEstimates$Summary[1,3],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.Age2.1'],KarlukOutmigration['KarlukSmolt2013.Age2.2'],KarlukOutmigration['KarlukSmolt2013.Age2.3'])*Karluk2013Age2StratifiedEstimates$Summary[1,5],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.Age2.1'],KarlukOutmigration['KarlukSmolt2013.Age2.2'],KarlukOutmigration['KarlukSmolt2013.Age2.3'])*Karluk2013Age2StratifiedEstimates$Summary[1,1],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.Age2.1'],KarlukOutmigration['KarlukSmolt2013.Age2.2'],KarlukOutmigration['KarlukSmolt2013.Age2.3'])*Karluk2013Age2StratifiedEstimates$Summary[1,2],format="f",digits=0,big.mark=",")))
## 2013 Karluk
KarlukTable6=rbind(KarlukTable6,cbind(AnnualShortDates['2013'],'',KarlukGroups[2],
                                      formatC(x=Karluk2013Age2StratifiedEstimates$Summary[2,4]*100,format="f",digits=1),formatC(x=Karluk2013Age2StratifiedEstimates$Summary[2,3]*100,format="f",digits=1),
                                      formatC(x=Karluk2013Age2StratifiedEstimates$Summary[2,5]*100,format="f",digits=1),formatC(x=sum(Karluk2013Age2StratifiedEstimates$Output[Karluk2013Age2StratifiedEstimates$Output[,2]==0])/length(Karluk2013Age2StratifiedEstimates$Output[,2]),format="f",digits=2),
                                      formatC(x=Karluk2013Age2StratifiedEstimates$Summary[2,1]*100,format="f",digits=1),formatC(x=Karluk2013Age2StratifiedEstimates$Summary[2,2]*100,format="f",digits=1),'',
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.Age2.1'],KarlukOutmigration['KarlukSmolt2013.Age2.2'],KarlukOutmigration['KarlukSmolt2013.Age2.3'])*Karluk2013Age2StratifiedEstimates$Summary[2,4],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.Age2.1'],KarlukOutmigration['KarlukSmolt2013.Age2.2'],KarlukOutmigration['KarlukSmolt2013.Age2.3'])*Karluk2013Age2StratifiedEstimates$Summary[2,3],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.Age2.1'],KarlukOutmigration['KarlukSmolt2013.Age2.2'],KarlukOutmigration['KarlukSmolt2013.Age2.3'])*Karluk2013Age2StratifiedEstimates$Summary[2,5],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.Age2.1'],KarlukOutmigration['KarlukSmolt2013.Age2.2'],KarlukOutmigration['KarlukSmolt2013.Age2.3'])*Karluk2013Age2StratifiedEstimates$Summary[2,1],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.Age2.1'],KarlukOutmigration['KarlukSmolt2013.Age2.2'],KarlukOutmigration['KarlukSmolt2013.Age2.3'])*Karluk2013Age2StratifiedEstimates$Summary[2,2],format="f",digits=0,big.mark=",")))
## 2013 Total Outmigration row
KarlukTable6=rbind(KarlukTable6,cbind('','','','','','','','','','','','2013 Total','',formatC(x=sum(KarlukOutmigration['KarlukSmolt2013.Age2.1'],KarlukOutmigration['KarlukSmolt2013.Age2.2'],KarlukOutmigration['KarlukSmolt2013.Age2.3']),format="f",digits=0,big.mark=","),''),Disclaimer)

### Write table
write.xlsx(x=as.data.frame(KarlukTable6),file="V:/WORK/Sockeye/Kodiak/2013 2014 Karluk Smolt/Tables/KarlukSmolt2013-2014Tables.xlsx",col.names=F,row.names=F,append=TRUE,sheetName="Karluk 2013 Smolt Age 2")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Table 7, 2014 smolt age 1 ####
## 17 rows by 15 columns for three strata/year and annual roll up.

KarlukCaption7 <- "Table X.-Estimates of stock composition and stock-specific outmigration for Karluk River Age 1 sockeye salmon smolt by stratum, 2014.  Reporting group-specific stock composition estimates (%) include median, 90% credibility interval, the probability that the group estimate is equal to zero (P=0), mean and SD.  Stock-specific estimates of outmigration are based upon mark-recapture estimates and variances (CV) of outmigration for each stratum.  See text for details."
Disclaimer <- cbind("Note: Stock composition estimates may not sum to 100% and stock-specific outmigration estimates may not sum to the total outmigration due to rounding error.",'','','','','','','','','','','','','','')

KarlukTable7=cbind(KarlukCaption7,'','','','','','','','','','','','','','')
KarlukTable7=rbind(KarlukTable7,cbind("Stratum","",'','Composition (%)','','','','','','',"Outmigration (number of fish)",'','','',''))
KarlukTable7=rbind(KarlukTable7,cbind('Period','Sample Size','Reporting','',"90% CI",'','','','','','',"90% CI",'','',''))
KarlukTable7=rbind(KarlukTable7,cbind("Dates",'CV','Group',"Median","5%","95%","P=0","Mean","SD",'',"Median","5%","95%","Mean","SD"))
## Early KarlukEarly
KarlukTable7=rbind(KarlukTable7,cbind('Early',paste('n=',KarlukFinalSampleSizes['KarlukSmolt2014.Age1.12'],sep=''),KarlukGroups[1],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.12']][1,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.12']][1,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.12']][1,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.12']][1,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.12']][1,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.12']][1,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.12']][1,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.12']][1,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.12']][1,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.12']][1,1],
                                                                                                                                            format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.12']][1,2],format="f",digits=0,big.mark=",")))
## Early KarlukLate
KarlukTable7=rbind(KarlukTable7,cbind(KarlukShortDates['KarlukSmolt2014.Age1.12'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2014.Age1.12']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.12']][2,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.12']][2,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.12']][2,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.12']][2,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.12']][2,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.12']][2,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.12']][2,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.12']][2,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.12']][2,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.12']][2,1],
                                                                                                                                            format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.12']][2,2],format="f",digits=0,big.mark=",")))
## Early Total Outmigration row
KarlukTable7=rbind(KarlukTable7,cbind('','','','','','','','','','','','Early Total','',formatC(x=KarlukOutmigration['KarlukSmolt2014.Age1.12'],format="f",digits=0,big.mark=","),''))
## Late KarlukEarly
KarlukTable7=rbind(KarlukTable7,cbind('Late',paste('n=',KarlukFinalSampleSizes['KarlukSmolt2014.Age1.3'],sep=''),KarlukGroups[1],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.3']][1,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.3']][1,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.3']][1,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.3']][1,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.3']][1,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.3']][1,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.3']][1,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.3']][1,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.3']][1,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.3']][1,1],
                                                                                                                                            format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.3']][1,2],format="f",digits=0,big.mark=",")))
## Late KarlukLate
KarlukTable7=rbind(KarlukTable7,cbind(KarlukShortDates['KarlukSmolt2014.Age1.3'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2014.Age1.3']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.3']][2,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.3']][2,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.3']][2,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.3']][2,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.3']][2,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age1.3']][2,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.3']][2,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.3']][2,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.3']][2,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.3']][2,1],
                                                                                                                                            format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age1.3']][2,2],format="f",digits=0,big.mark=",")))
## Late Total Outmigration row
KarlukTable7=rbind(KarlukTable7,cbind('','','','','','','','','','','','Late Total','',formatC(x=KarlukOutmigration['KarlukSmolt2014.Age1.3'],format="f",digits=0,big.mark=","),''))
## 2014 KarlukEarly
KarlukTable7=rbind(KarlukTable7,cbind('2014',paste('n=',sum(KarlukFinalSampleSizes['KarlukSmolt2014.Age1.12'],KarlukFinalSampleSizes['KarlukSmolt2014.Age1.3']),sep=''),KarlukGroups[1],
                                      formatC(x=Karluk2014Age1StratifiedEstimates$Summary[1,4]*100,format="f",digits=1),formatC(x=Karluk2014Age1StratifiedEstimates$Summary[1,3]*100,format="f",digits=1),
                                      formatC(x=Karluk2014Age1StratifiedEstimates$Summary[1,5]*100,format="f",digits=1),formatC(x=sum(Karluk2014Age1StratifiedEstimates$Output[Karluk2014Age1StratifiedEstimates$Output[,1]==0])/length(Karluk2014Age1StratifiedEstimates$Output[,1]),format="f",digits=2),
                                      formatC(x=Karluk2014Age1StratifiedEstimates$Summary[1,1]*100,format="f",digits=1),formatC(x=Karluk2014Age1StratifiedEstimates$Summary[1,2]*100,format="f",digits=1),'',
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age1.12'],KarlukOutmigration['KarlukSmolt2014.Age1.3'])*Karluk2014Age1StratifiedEstimates$Summary[1,4],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age1.12'],KarlukOutmigration['KarlukSmolt2014.Age1.3'])*Karluk2014Age1StratifiedEstimates$Summary[1,3],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age1.12'],KarlukOutmigration['KarlukSmolt2014.Age1.3'])*Karluk2014Age1StratifiedEstimates$Summary[1,5],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age1.12'],KarlukOutmigration['KarlukSmolt2014.Age1.3'])*Karluk2014Age1StratifiedEstimates$Summary[1,1],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age1.12'],KarlukOutmigration['KarlukSmolt2014.Age1.3'])*Karluk2014Age1StratifiedEstimates$Summary[1,2],format="f",digits=0,big.mark=",")))
## 2014 Karluk
KarlukTable7=rbind(KarlukTable7,cbind(AnnualShortDates['2014'],'',KarlukGroups[2],
                                      formatC(x=Karluk2014Age1StratifiedEstimates$Summary[2,4]*100,format="f",digits=1),formatC(x=Karluk2014Age1StratifiedEstimates$Summary[2,3]*100,format="f",digits=1),
                                      formatC(x=Karluk2014Age1StratifiedEstimates$Summary[2,5]*100,format="f",digits=1),formatC(x=sum(Karluk2014Age1StratifiedEstimates$Output[Karluk2014Age1StratifiedEstimates$Output[,2]==0])/length(Karluk2014Age1StratifiedEstimates$Output[,2]),format="f",digits=2),
                                      formatC(x=Karluk2014Age1StratifiedEstimates$Summary[2,1]*100,format="f",digits=1),formatC(x=Karluk2014Age1StratifiedEstimates$Summary[2,2]*100,format="f",digits=1),'',
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age1.12'],KarlukOutmigration['KarlukSmolt2014.Age1.3'])*Karluk2014Age1StratifiedEstimates$Summary[2,4],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age1.12'],KarlukOutmigration['KarlukSmolt2014.Age1.3'])*Karluk2014Age1StratifiedEstimates$Summary[2,3],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age1.12'],KarlukOutmigration['KarlukSmolt2014.Age1.3'])*Karluk2014Age1StratifiedEstimates$Summary[2,5],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age1.12'],KarlukOutmigration['KarlukSmolt2014.Age1.3'])*Karluk2014Age1StratifiedEstimates$Summary[2,1],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age1.12'],KarlukOutmigration['KarlukSmolt2014.Age1.3'])*Karluk2014Age1StratifiedEstimates$Summary[2,2],format="f",digits=0,big.mark=",")))
## 2014 Total Outmigration row
KarlukTable7=rbind(KarlukTable7,cbind('','','','','','','','','','','','2014 Total','',formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age1.12'],KarlukOutmigration['KarlukSmolt2014.Age1.3']),format="f",digits=0,big.mark=","),''),Disclaimer)

### Write table
write.xlsx(x=as.data.frame(KarlukTable7),file="V:/WORK/Sockeye/Kodiak/2013 2014 Karluk Smolt/Tables/KarlukSmolt2013-2014Tables.xlsx",col.names=F,row.names=F,append=TRUE,sheetName="Karluk 2014 Smolt Age 1")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Table 8, 2014 smolt age 2 ####
## 17 rows by 15 columns for three strata/year and annual roll up.

KarlukCaption8 <- "Table X.-Estimates of stock composition and stock-specific outmigration for Karluk River Age 2 sockeye salmon smolt by stratum, 2014.  Reporting group-specific stock composition estimates (%) include median, 90% credibility interval, the probability that the group estimate is equal to zero (P=0), mean and SD.  Stock-specific estimates of outmigration are based upon mark-recapture estimates and variances (CV) of outmigration for each stratum.  See text for details."
Disclaimer <- cbind("Note: Stock composition estimates may not sum to 100% and stock-specific outmigration estimates may not sum to the total outmigration due to rounding error.",'','','','','','','','','','','','','','')

KarlukTable8=cbind(KarlukCaption8,'','','','','','','','','','','','','','')
KarlukTable8=rbind(KarlukTable8,cbind("Stratum","",'','Composition (%)','','','','','','',"Outmigration (number of fish)",'','','',''))
KarlukTable8=rbind(KarlukTable8,cbind('Period','Sample Size','Reporting','',"90% CI",'','','','','','',"90% CI",'','',''))
KarlukTable8=rbind(KarlukTable8,cbind("Dates",'CV','Group',"Median","5%","95%","P=0","Mean","SD",'',"Median","5%","95%","Mean","SD"))
## Early KarlukEarly
KarlukTable8=rbind(KarlukTable8,cbind('Early',paste('n=',KarlukFinalSampleSizes['KarlukSmolt2014.Age2.1'],sep=''),KarlukGroups[1],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.1']][1,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.1']][1,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.1']][1,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.1']][1,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.1']][1,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.1']][1,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.1']][1,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.1']][1,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.1']][1,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.1']][1,1],
                                                                                                                                            format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.1']][1,2],format="f",digits=0,big.mark=",")))
## Early KarlukLate
KarlukTable8=rbind(KarlukTable8,cbind(KarlukShortDates['KarlukSmolt2014.Age2.1'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2014.Age2.1']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.1']][2,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.1']][2,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.1']][2,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.1']][2,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.1']][2,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.1']][2,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.1']][2,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.1']][2,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.1']][2,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.1']][2,1],
                                                                                                                                            format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.1']][2,2],format="f",digits=0,big.mark=",")))
## Early Total Outmigration row
KarlukTable8=rbind(KarlukTable8,cbind('','','','','','','','','','','','Early Total','',formatC(x=KarlukOutmigration['KarlukSmolt2014.Age2.1'],format="f",digits=0,big.mark=","),''))
## Middle KarlukEarly
KarlukTable8=rbind(KarlukTable8,cbind('Middle',paste('n=',KarlukFinalSampleSizes['KarlukSmolt2014.Age2.2'],sep=''),KarlukGroups[1],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.2']][1,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.2']][1,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.2']][1,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.2']][1,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.2']][1,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.2']][1,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.2']][1,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.2']][1,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.2']][1,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.2']][1,1],
                                                                                                                                            format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.2']][1,2],format="f",digits=0,big.mark=",")))
## Middle KarlukLate
KarlukTable8=rbind(KarlukTable8,cbind(KarlukShortDates['KarlukSmolt2014.Age2.2'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2014.Age2.2']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.2']][2,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.2']][2,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.2']][2,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.2']][2,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.2']][2,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.2']][2,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.2']][2,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.2']][2,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.2']][2,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.2']][2,1],
                                                                                                                                            format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.2']][2,2],format="f",digits=0,big.mark=",")))
## Middle Total Outmigration row
KarlukTable8=rbind(KarlukTable8,cbind('','','','','','','','','','','','Middle Total','',formatC(x=KarlukOutmigration['KarlukSmolt2014.Age2.2'],format="f",digits=0,big.mark=","),''))
## Late KarlukEarly
KarlukTable8=rbind(KarlukTable8,cbind('Late',paste('n=',KarlukFinalSampleSizes['KarlukSmolt2014.Age2.3'],sep=''),KarlukGroups[1],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.3']][1,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.3']][1,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.3']][1,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.3']][1,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.3']][1,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.3']][1,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.3']][1,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.3']][1,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.3']][1,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.3']][1,1],
                                                                                                                                            format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.3']][1,2],format="f",digits=0,big.mark=",")))
## Late KarlukLate
KarlukTable8=rbind(KarlukTable8,cbind(KarlukShortDates['KarlukSmolt2014.Age2.3'],paste('CV=',formatC(x=KarlukCVs['KarlukSmolt2014.Age2.3']*100,format="f",digits=1),"%",sep=''),KarlukGroups[2],
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.3']][2,3]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.3']][2,4]*100,format="f",digits=1),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.3']][2,5]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.3']][2,6],format="f",digits=2),
                                      formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.3']][2,1]*100,format="f",digits=1),formatC(x=KarlukSmolt_Age_Estimates[['KarlukSmolt2014.Age2.3']][2,2]*100,format="f",digits=1),'',
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.3']][2,3],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.3']][2,4],format="f",digits=0,big.mark=","),
                                      formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.3']][2,5],format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.3']][2,1],
                                                                                                                                            format="f",digits=0,big.mark=","),formatC(x=KarlukAbundances[['KarlukSmolt2014.Age2.3']][2,2],format="f",digits=0,big.mark=",")))
## Late Total Outmigration row
KarlukTable8=rbind(KarlukTable8,cbind('','','','','','','','','','','','Late Total','',formatC(x=KarlukOutmigration['KarlukSmolt2014.Age2.3'],format="f",digits=0,big.mark=","),''))
## 2014 KarlukEarly
KarlukTable8=rbind(KarlukTable8,cbind('2014',paste('n=',sum(KarlukFinalSampleSizes['KarlukSmolt2014.Age2.1'],KarlukFinalSampleSizes['KarlukSmolt2014.Age2.2'],KarlukFinalSampleSizes['KarlukSmolt2014.Age2.3']),sep=''),KarlukGroups[1],
                                      formatC(x=Karluk2014Age2StratifiedEstimates$Summary[1,4]*100,format="f",digits=1),formatC(x=Karluk2014Age2StratifiedEstimates$Summary[1,3]*100,format="f",digits=1),
                                      formatC(x=Karluk2014Age2StratifiedEstimates$Summary[1,5]*100,format="f",digits=1),formatC(x=sum(Karluk2014Age2StratifiedEstimates$Output[Karluk2014Age2StratifiedEstimates$Output[,1]==0])/length(Karluk2014Age2StratifiedEstimates$Output[,1]),format="f",digits=2),
                                      formatC(x=Karluk2014Age2StratifiedEstimates$Summary[1,1]*100,format="f",digits=1),formatC(x=Karluk2014Age2StratifiedEstimates$Summary[1,2]*100,format="f",digits=1),'',
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age2.1'],KarlukOutmigration['KarlukSmolt2014.Age2.2'],KarlukOutmigration['KarlukSmolt2014.Age2.3'])*Karluk2014Age2StratifiedEstimates$Summary[1,4],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age2.1'],KarlukOutmigration['KarlukSmolt2014.Age2.2'],KarlukOutmigration['KarlukSmolt2014.Age2.3'])*Karluk2014Age2StratifiedEstimates$Summary[1,3],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age2.1'],KarlukOutmigration['KarlukSmolt2014.Age2.2'],KarlukOutmigration['KarlukSmolt2014.Age2.3'])*Karluk2014Age2StratifiedEstimates$Summary[1,5],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age2.1'],KarlukOutmigration['KarlukSmolt2014.Age2.2'],KarlukOutmigration['KarlukSmolt2014.Age2.3'])*Karluk2014Age2StratifiedEstimates$Summary[1,1],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age2.1'],KarlukOutmigration['KarlukSmolt2014.Age2.2'],KarlukOutmigration['KarlukSmolt2014.Age2.3'])*Karluk2014Age2StratifiedEstimates$Summary[1,2],format="f",digits=0,big.mark=",")))
## 2014 Karluk
KarlukTable8=rbind(KarlukTable8,cbind(AnnualShortDates['2014'],'',KarlukGroups[2],
                                      formatC(x=Karluk2014Age2StratifiedEstimates$Summary[2,4]*100,format="f",digits=1),formatC(x=Karluk2014Age2StratifiedEstimates$Summary[2,3]*100,format="f",digits=1),
                                      formatC(x=Karluk2014Age2StratifiedEstimates$Summary[2,5]*100,format="f",digits=1),formatC(x=sum(Karluk2014Age2StratifiedEstimates$Output[Karluk2014Age2StratifiedEstimates$Output[,2]==0])/length(Karluk2014Age2StratifiedEstimates$Output[,2]),format="f",digits=2),
                                      formatC(x=Karluk2014Age2StratifiedEstimates$Summary[2,1]*100,format="f",digits=1),formatC(x=Karluk2014Age2StratifiedEstimates$Summary[2,2]*100,format="f",digits=1),'',
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age2.1'],KarlukOutmigration['KarlukSmolt2014.Age2.2'],KarlukOutmigration['KarlukSmolt2014.Age2.3'])*Karluk2014Age2StratifiedEstimates$Summary[2,4],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age2.1'],KarlukOutmigration['KarlukSmolt2014.Age2.2'],KarlukOutmigration['KarlukSmolt2014.Age2.3'])*Karluk2014Age2StratifiedEstimates$Summary[2,3],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age2.1'],KarlukOutmigration['KarlukSmolt2014.Age2.2'],KarlukOutmigration['KarlukSmolt2014.Age2.3'])*Karluk2014Age2StratifiedEstimates$Summary[2,5],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age2.1'],KarlukOutmigration['KarlukSmolt2014.Age2.2'],KarlukOutmigration['KarlukSmolt2014.Age2.3'])*Karluk2014Age2StratifiedEstimates$Summary[2,1],format="f",digits=0,big.mark=","),
                                      formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age2.1'],KarlukOutmigration['KarlukSmolt2014.Age2.2'],KarlukOutmigration['KarlukSmolt2014.Age2.3'])*Karluk2014Age2StratifiedEstimates$Summary[2,2],format="f",digits=0,big.mark=",")))
## 2014 Total Outmigration row
KarlukTable8=rbind(KarlukTable8,cbind('','','','','','','','','','','','2014 Total','',formatC(x=sum(KarlukOutmigration['KarlukSmolt2014.Age2.1'],KarlukOutmigration['KarlukSmolt2014.Age2.2'],KarlukOutmigration['KarlukSmolt2014.Age2.3']),format="f",digits=0,big.mark=","),''),Disclaimer)

### Write table
write.xlsx(x=as.data.frame(KarlukTable8),file="V:/WORK/Sockeye/Kodiak/2013 2014 Karluk Smolt/Tables/KarlukSmolt2013-2014Tables.xlsx",col.names=F,row.names=F,append=TRUE,sheetName="Karluk 2014 Smolt Age 2")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Cool Analyses of IA Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

objects(pattern="Metadata")
str(Metadata2013)

# Exploratory plots
require(lattice)

xyplot(FL ~ AGE, data=Metadata2013, groups=StrictAssignment, pch=c(16,17), cex=2, grid=TRUE)
xyplot(FL ~ AGE, data=Metadata2013, groups=StrictAssignment, pch=c(16,17), cex=2, grid=TRUE, panel=panel.smoothScatter)

xyplot(FL ~ AGE | StrictAssignment, data=Metadata2013, pch=16, cex=2, group=StrictAssignment, col=c("blue", "red"))
xyplot(FL ~ AGE | StrictAssignment, data=Metadata2013, pch=16, grid=TRUE, panel=panel.smoothScatter)

xyplot(FL ~ AGE | StrictAssignment, data=Metadata2014, pch=16, cex=2, group=StrictAssignment, col=c("blue", "red"))
xyplot(FL ~ AGE | StrictAssignment, data=Metadata2014, pch=16, grid=TRUE, panel=panel.smoothScatter)


## 2013
# Comparing FL by Age and RG
boxplot(FL ~ StrictAssignment+AGE, data=Metadata2013, col=c("red", "blue"), notch=FALSE, pch=16,  pars=list(boxwex=0.8, staplewex=0.5, outwex=0.5), at=c(1,2,4,5,7,8,10,11,13,14), frame=FALSE, axes=FALSE, ylim=c(90,200), xlim=c(0,14))
axis(side=2, at=seq(from=90, to=200, by=10), lwd=4, cex.axis=1.4)
axis(side=1, pos=90, at=seq(from=1.5, to=13.5, by=3), labels=0:4, cex.axis=1.4, lwd=4); abline(h=90, lwd=4)
mtext(text="Age", side=1, cex=2, line=2)
mtext(text="FL (mm)", side=2, cex=2, line=2.5)
mtext(text="2013 - Strict Assignments", side=3, cex=2, line=1.5)
legend("topleft", legend=c("Early", "Late"), fill=c("red", "blue"), bty="n", cex=2)

# Comparing K by Age and RG
boxplot(K ~ StrictAssignment+AGE, data=Metadata2013, col=c("red", "blue"), notch=FALSE, pch=16,  pars=list(boxwex=0.8, staplewex=0.5, outwex=0.5), at=c(1,2,4,5,7,8,10,11,13,14), frame=FALSE, axes=FALSE, ylim=c(0.6,1.4), xlim=c(0,14))
axis(side=2, at=seq(from=0.6, to=1.4, by=0.2), lwd=4, cex.axis=1.4)
axis(side=1, pos=0.6, at=seq(from=1.5, to=13.5, by=3), labels=0:4, cex.axis=1.4, lwd=4); abline(h=0.6, lwd=4)
mtext(text="Age", side=1, cex=2, line=2.5)
mtext(text="K", side=2, cex=2, line=2.5)
mtext(text="2013 - Strict Assignments", side=3, cex=2, line=1.5)
legend("topleft", legend=c("Early", "Late"), fill=c("red", "blue"), bty="n", cex=2)

# Comparing Age proportions
table(Metadata2013$AGE, Metadata2013$StrictAssignment)
chisq.test(table(Metadata2013$AGE, Metadata2013$StrictAssignment))
table(Metadata2013$AGE, Metadata2013$RelaxedAssignment)
chisq.test(table(Metadata2013$AGE, Metadata2013$RelaxedAssignment))


## 2014
# Clean up
Metadata2014$FL[which(Metadata2014$FL==0)] <- NA
Metadata2014$WEIGHT[which(Metadata2014$WEIGHT==0)] <- NA

# Comparing FL by Age and RG
boxplot(FL ~ StrictAssignment+AGE, data=Metadata2014, col=c("red", "blue"), notch=FALSE, pch=16,  pars=list(boxwex=0.8, staplewex=0.5, outwex=0.5), at=c(1,2,4,5,7,8,10,11), frame=FALSE, axes=FALSE, ylim=c(90,200), xlim=c(0,14))
axis(side=2, at=seq(from=90, to=200, by=10), lwd=4, cex.axis=1.4)
axis(side=1, pos=90, at=seq(from=1.5, to=13.5, by=3), labels=0:4, cex.axis=1.4, lwd=4); abline(h=90, lwd=4)
mtext(text="Age", side=1, cex=2, line=2)
mtext(text="FL (mm)", side=2, cex=2, line=2.5)
mtext(text="2014 - Strict Assignments", side=3, cex=2, line=1.5)
legend("topleft", legend=c("Early", "Late"), fill=c("red", "blue"), bty="n", cex=2)

# Comparing K by Age and RG
boxplot(K ~ StrictAssignment+AGE, data=Metadata2014, col=c("red", "blue"), notch=FALSE, pch=16,  pars=list(boxwex=0.8, staplewex=0.5, outwex=0.5), at=c(1,2,4,5,7,8,10,11), frame=FALSE, axes=FALSE, ylim=c(0.6,1.4), xlim=c(0,14))
axis(side=2, at=seq(from=0.6, to=1.4, by=0.2), lwd=4, cex.axis=1.4)
axis(side=1, pos=0.6, at=seq(from=1.5, to=13.5, by=3), labels=0:4, cex.axis=1.4, lwd=4); abline(h=0.6, lwd=4)
mtext(text="Age", side=1, cex=2, line=2.5)
mtext(text="K", side=2, cex=2, line=2.5)
mtext(text="2014 - Strict Assignments", side=3, cex=2, line=1.5)
legend("topleft", legend=c("Early", "Late"), fill=c("red", "blue"), bty="n", cex=2)


# Comparing Age proportions
table(Metadata2014$AGE, Metadata2014$StrictAssignment)
chisq.test(table(Metadata2014$AGE, Metadata2014$StrictAssignment))
table(Metadata2014$AGE, Metadata2014$RelaxedAssignment)
chisq.test(table(Metadata2014$AGE, Metadata2014$RelaxedAssignment))




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Redo stratified estimates to incorporate uncertainty in OM  estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(beepr)
KarlukMixtures_All <- c(KarlukMixtures, KarlukMixtures_Age)
KarlukSmoltPosteriors <- CustomCombineBAYESOutput.GCL(groupvec = 1:2, groupnames = KarlukGroups2, maindir = "BAYES/Output", mixvec = KarlukMixtures_All, prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE); beep(8)
dir.create("Estimate objects")
dput(x = KarlukSmoltPosteriors, file = "Estimate objects/KarlukSmoltPosteriors.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Tables 1 & 2 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## First, do this with 2013 and 2014 overall since these have easy to get CVs
# Mixtures
KarlukMixturesWithAbundances <- KarlukMixtures_All[1:6]
# OM estimates
Karluk2013Abundances <- setNames(object = as.numeric(readClipboard()), nm = KarlukMixturesWithAbundances[1:3])
Karluk2014Abundances <- setNames(object = as.numeric(readClipboard()), nm = KarlukMixturesWithAbundances[4:6])
# CVs
Karluk2013CVs <- setNames(object = as.numeric(readClipboard()), nm = KarlukMixturesWithAbundances[1:3])
Karluk2014CVs <- setNames(object = as.numeric(readClipboard()), nm = KarlukMixturesWithAbundances[4:6]) # see 2014_popest_karluk.xls or "Daily w CV" tab in 2014Karluk_smolt_pop_est_by_day_age.xlsx


# Define function (took from Tyler's Chignik Smolt script)
IncorporateAbundanceErrors.GCL <- function(posteriorsname,groupnames,mixvec,catchvec,CVvec,alpha=0.1){
  
  results=vector("list",length(mixvec))
  names(results)=mixvec
  Hdata=H=sapply(mixvec,function(mix){NULL},simplify=FALSE)
  
  for(mix in mixvec){
    
    lnvar=log(CVvec[mix]^2+1)
    
    lnmean=log(catchvec[mix])-lnvar/2
    
    H[[mix]]=rlnorm(length(get(posteriorsname)[[2]][[mix]][,1]),lnmean,sqrt(lnvar))
    
    Hdata[[mix]]=H[[mix]]*get(posteriorsname)[[2]][[mix]]
  }
  
  for(mix in mixvec){
    
    results[[mix]]=array(NA,c(length(groupnames),6),dimnames=list(groupnames,c("mean","sd","median",paste(round(alpha/2,3)*100,"%",sep=""),paste(round(1-alpha/2,3)*100,"%",sep=""),"P=0")))
    
    results[[mix]][groupnames,1]=apply(Hdata[[mix]],2,mean)
    
    results[[mix]][groupnames,2]=apply(Hdata[[mix]],2,sd)
    
    results[[mix]][groupnames,3]=apply(Hdata[[mix]],2,median)
    
    results[[mix]][groupnames,4]=apply(Hdata[[mix]],2,function(clm){quantile(clm,alpha/2)})
    
    results[[mix]][groupnames,5]=apply(Hdata[[mix]],2,function(clm){quantile(clm,1-alpha/2)})
    
    results[[mix]][groupnames,6]=apply(Hdata[[mix]],2,function(clm){sum(clm==0)/length(clm)})
  }
  return(list(Stats=results,Output=Hdata))
}

Karluk2014AbundancesWithMRError <- IncorporateAbundanceErrors.GCL(posteriorsname = 'KarlukSmoltPosteriors', groupnames = KarlukGroups2, mixvec = KarlukMixturesWithAbundances[4:6], catchvec = Karluk2014Abundances, CVvec = Karluk2014CVs, alpha = 0.1)
Karluk2013AbundancesWithMRError <- IncorporateAbundanceErrors.GCL(posteriorsname = 'KarlukSmoltPosteriors', groupnames = KarlukGroups2, mixvec = KarlukMixturesWithAbundances[1:3], catchvec = Karluk2013Abundances, CVvec = Karluk2013CVs, alpha = 0.1)

# Stratified Estimates
Karluk2014StratifiedEstimatesWithMRError <- StratifiedEstimator.GCL(groupvec = 1:2, groupnames = KarlukGroups2, maindir = "BAYES/Output", mixvec = KarlukMixturesWithAbundances[4:6], catchvec = Karluk2014Abundances, CVvec = Karluk2014CVs, newname = "Karluk2014StratifiedEstimatesWithMRError", priorname = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1)
Karluk2013StratifiedEstimatesWithMRError <- StratifiedEstimator.GCL(groupvec = 1:2, groupnames = KarlukGroups2, maindir = "BAYES/Output", mixvec = KarlukMixturesWithAbundances[1:3], catchvec = Karluk2013Abundances, CVvec = Karluk2013CVs, newname = "Karluk2013StratifiedEstimatesWithMRError", priorname = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1)

# Stratified Abundances
Karluk2014StratifiedAbundancesWithMRError <- Karluk2014StratifiedEstimatesWithMRError$Summary * sum(Karluk2014Abundances)
Karluk2013StratifiedAbundancesWithMRError <- Karluk2013StratifiedEstimatesWithMRError$Summary * sum(Karluk2013Abundances)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Make tables...again
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~
## Table 1, 2013 Smolt
sheetname = "Karluk 2013 Smolt"
tablenum = 1  
Nrows <- 17
Ncols <- 15
yr = 2013
byAge = FALSE
Ages = FALSE
ages = "Age 2"

dates.vec <- KarlukShortDates[1:3]
samplesizes.vec <- KarlukFinalSampleSizes[1:3]

stats.proportions.list <- KarlukSmoltPosteriors$Stats[1:3]  # List of proportion stats
stratified.proportions <- Karluk2013StratifiedEstimatesWithMRError$Summary  # Stratified porportion stats
stats.abundances.list <- Karluk2013AbundancesWithMRError$Stats[1:3]  # List of abundance stats
stratified.abundances <- Karluk2013StratifiedAbundancesWithMRError  # Stratified abundance stats
abundances.vec <- Karluk2013Abundances  # Vector of abundances
CV.vec <- Karluk2013CVs  # Vector of CVs

#~~~~~~~~
KarlukRedoTable <- matrix(data = "", nrow = Nrows, ncol = Ncols)

KarlukRedoTable[1, 1] <- paste("Table ", tablenum,".-Estimates of stock composition and stock-specific outmigration for Karluk River ", ifelse(Ages, paste(ages, " ", sep = ''), ''),
                               "sockeye salmon smolt ", ifelse(byAge, "by age", "by stratum"), ", ", yr,
                               ".  Reporting group-specific stock composition estimates (%) include median, 90% credibility interval, the probability that the group estimate is equal to zero (P=0), mean, and SD.  Stock-specific estimates of outmigration are based upon mark-recapture estimates and variances (CV) of outmigration for each stratum. See text for details.", sep = "")
KarlukRedoTable[Nrows, 1] <- "Note: Stock composition estimates may not sum to 100% and stock-specific outmigration estimates may not sum to the total outmigration due to rounding error."

KarlukRedoTable[2, ] <- c("Stratum","",'','Composition (%)','','','','','','',"Outmigration (number of fish)",'','','','')
KarlukRedoTable[3, ] <- c('Period','Sample Size','Reporting','',"90% CI",'','','','','','',"90% CI",'','','')
KarlukRedoTable[4, ] <- c("Dates",'CV','Group',"Median","5%","95%","P=0","Mean","SD",'',"Median","5%","95%","Mean","SD")
KarlukRedoTable[c(5, 8, 11, 14), 1] <- c("Early", "Middle", "Late", yr)

KarlukRedoTable[c(6, 9, 12, 15), 1] <- c(dates.vec, paste(unlist(strsplit(x = dates.vec[1], split = "-"))[1],
                                                          ifelse(length(grep(pattern = "/", x = unlist(strsplit(x = dates.vec[3], split = "-"))[2])) == 0, 
                                                                 paste(unlist(strsplit(x = unlist(strsplit(x = dates.vec[3], split = "-"))[1], split = "/"))[1],
                                                                       unlist(strsplit(x = dates.vec[3], split = "-"))[2],
                                                                       sep = "/"), 
                                                                 unlist(strsplit(x = dates.vec[3], split = "-"))[2]),
                                                          sep = '-'))
KarlukRedoTable[c(5:6, 8:9, 11:12, 14:15), 3] <- KarlukGroups

KarlukRedoTable[c(5, 8, 11, 14), 2] <- paste("n=", c(samplesizes.vec, sum(samplesizes.vec)), sep = '')
KarlukRedoTable[c(6, 9, 12), 2] <- paste("CV=", formatC(x = CV.vec * 100, digits = 1, format = "f"), "%", sep = '')

KarlukRedoTable[c(5:6, 8:9, 11:12), c(4:6, 8:9)] <- formatC(x = do.call(what = "rbind", args = stats.proportions.list[1:3])[, c("median", "5%", "95%", "mean", "sd")] * 100, format = "f", digits = 1)
KarlukRedoTable[c(5:6, 8:9, 11:12), 7] <- formatC(x = do.call(what = "rbind", args = stats.proportions.list[1:3])[, c("P=0")] * 100, format = "f", digits = 2)
KarlukRedoTable[c(5:6, 8:9, 11:12), 11:15] <- formatC(x = do.call(what = "rbind", args = stats.abundances.list)[, c("median", "5%", "95%", "mean", "sd")], digits = 0, format = "f", big.mark = ",")

KarlukRedoTable[c(7, 10, 13, 16), 12] <- paste(c("Early", "Middle", "Late", yr), "Total", sep = ' ')
KarlukRedoTable[c(7, 10, 13, 16), 14] <- formatC(x = c(abundances.vec, sum(abundances.vec)), digits = 0, format = "f", big.mark = ",")

KarlukRedoTable[14:15, c(4:6, 8:9)] <- formatC(x = stratified.proportions[, toupper(c("median", "5%", "95%", "mean", "sd"))] * 100, format = "f", digits = 1)
KarlukRedoTable[14:15, 7] <- formatC(x = stratified.proportions[, c("P0")] * 100, format = "f", digits = 2)
KarlukRedoTable[14:15, 11:15] <- formatC(x = stratified.abundances[, toupper(c("median", "5%", "95%", "mean", "sd"))], format = "f", digits = 0, big.mark = ",")


require(xlsx)
write.xlsx(x = KarlukRedoTable, file = "Tables/KarlukSmolt2013-2014TablesNEW.xlsx", sheetName = sheetname, col.names = FALSE, row.names = FALSE, append = TRUE)
#~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~
## Table 2, 2014 Smolt
sheetname = "Karluk 2014 Smolt"
tablenum = 2
Nrows <- 17
Ncols <- 15
yr = 2014
byAge = FALSE
Ages = FALSE
ages = "Age 2"

dates.vec <- KarlukShortDates[9:11]
samplesizes.vec <- KarlukFinalSampleSizes[9:11]

names(KarlukSmoltPosteriors$Stats)
stats.proportions.list <- KarlukSmoltPosteriors$Stats[4:6]  # List of proportion stats
stratified.proportions <- Karluk2014StratifiedEstimatesWithMRError$Summary  # Stratified porportion stats
stats.abundances.list <- Karluk2014AbundancesWithMRError$Stats[1:3]  # List of abundance stats
stratified.abundances <- Karluk2014StratifiedAbundancesWithMRError  # Stratified abundance stats
abundances.vec <- Karluk2014Abundances  # Vector of abundances
CV.vec <- Karluk2014CVs  # Vector of CVs

#~~~~~~~~
KarlukRedoTable <- matrix(data = "", nrow = Nrows, ncol = Ncols)

KarlukRedoTable[1, 1] <- paste("Table ", tablenum,".-Estimates of stock composition and stock-specific outmigration for Karluk River ", ifelse(Ages, paste(ages, " ", sep = ''), ''),
                               "sockeye salmon smolt ", ifelse(byAge, "by age", "by stratum"), ", ", yr,
                               ".  Reporting group-specific stock composition estimates (%) include median, 90% credibility interval, the probability that the group estimate is equal to zero (P=0), mean, and SD.  Stock-specific estimates of outmigration are based upon mark-recapture estimates and variances (CV) of outmigration for each stratum. See text for details.", sep = "")
KarlukRedoTable[Nrows, 1] <- "Note: Stock composition estimates may not sum to 100% and stock-specific outmigration estimates may not sum to the total outmigration due to rounding error."

KarlukRedoTable[2, ] <- c("Stratum","",'','Composition (%)','','','','','','',"Outmigration (number of fish)",'','','','')
KarlukRedoTable[3, ] <- c('Period','Sample Size','Reporting','',"90% CI",'','','','','','',"90% CI",'','','')
KarlukRedoTable[4, ] <- c("Dates",'CV','Group',"Median","5%","95%","P=0","Mean","SD",'',"Median","5%","95%","Mean","SD")
KarlukRedoTable[c(5, 8, 11, 14), 1] <- c("Early", "Middle", "Late", yr)

KarlukRedoTable[c(6, 9, 12, 15), 1] <- c(dates.vec, paste(unlist(strsplit(x = dates.vec[1], split = "-"))[1],
                                                          ifelse(length(grep(pattern = "/", x = unlist(strsplit(x = dates.vec[3], split = "-"))[2])) == 0, 
                                                                 paste(unlist(strsplit(x = unlist(strsplit(x = dates.vec[3], split = "-"))[1], split = "/"))[1],
                                                                       unlist(strsplit(x = dates.vec[3], split = "-"))[2],
                                                                       sep = "/"), 
                                                                 unlist(strsplit(x = dates.vec[3], split = "-"))[2]),
                                                          sep = '-'))
KarlukRedoTable[c(5:6, 8:9, 11:12, 14:15), 3] <- KarlukGroups

KarlukRedoTable[c(5, 8, 11, 14), 2] <- paste("n=", c(samplesizes.vec, sum(samplesizes.vec)), sep = '')
KarlukRedoTable[c(6, 9, 12), 2] <- paste("CV=", formatC(x = CV.vec * 100, digits = 1, format = "f"), "%", sep = '')

KarlukRedoTable[c(5:6, 8:9, 11:12), c(4:6, 8:9)] <- formatC(x = do.call(what = "rbind", args = stats.proportions.list[1:3])[, c("median", "5%", "95%", "mean", "sd")] * 100, format = "f", digits = 1)
KarlukRedoTable[c(5:6, 8:9, 11:12), 7] <- formatC(x = do.call(what = "rbind", args = stats.proportions.list[1:3])[, c("P=0")] * 100, format = "f", digits = 2)
KarlukRedoTable[c(5:6, 8:9, 11:12), 11:15] <- formatC(x = do.call(what = "rbind", args = stats.abundances.list)[, c("median", "5%", "95%", "mean", "sd")], digits = 0, format = "f", big.mark = ",")

KarlukRedoTable[c(7, 10, 13, 16), 12] <- paste(c("Early", "Middle", "Late", yr), "Total", sep = ' ')
KarlukRedoTable[c(7, 10, 13, 16), 14] <- formatC(x = c(abundances.vec, sum(abundances.vec)), digits = 0, format = "f", big.mark = ",")

KarlukRedoTable[14:15, c(4:6, 8:9)] <- formatC(x = stratified.proportions[, toupper(c("median", "5%", "95%", "mean", "sd"))] * 100, format = "f", digits = 1)
KarlukRedoTable[14:15, 7] <- formatC(x = stratified.proportions[, c("P0")] * 100, format = "f", digits = 2)
KarlukRedoTable[14:15, 11:15] <- formatC(x = stratified.abundances[, toupper(c("median", "5%", "95%", "mean", "sd"))], format = "f", digits = 0, big.mark = ",")


require(xlsx)
write.xlsx(x = KarlukRedoTable, file = "Tables/KarlukSmolt2013-2014TablesNEW.xlsx", sheetName = sheetname, col.names = FALSE, row.names = FALSE, append = TRUE)
#~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Table 3 is the same ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Tables 4-8 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Do math 1st, then build tables 2nd

require(xlsx); require(beepr)

# Define function
rdirich <- function(alpha) {
  n = length(alpha)
  alphaprior = alpha + 1/n
  gg = rgamma(n = n, shape = alpha, scale = 1)
  return(gg/sum(gg))
}


#~~~~~~~~~~~~~~~~~~
## 2014
dat <- read.xlsx(file = "2014Karluk_smolt_pop_est_by_day_age.xlsx", sheetName = "AgeSpecificOM")
str(dat)

boots = 100000

dirich <- replicate(n = boots, expr = apply(dat[, c("Age1", "Age2", "Age3", "Age4")], 1, function(day) {rdirich(alpha = day)} )); beep(2)
dimnames(dirich) <- list(Ages = c("Age1", "Age2", "Age3", "Age4"), Day = seq(dim(dirich)[2]), Boot = seq(boots))
str(dirich)

N <- t(apply(dat[, c("DailyOM", "CV")], 1, function(day) {lnvar = log(day[2]^2 + 1)
                                                          lnmean = log(day[1]) - lnvar/2
                                                          rlnorm(n = boots, meanlog = lnmean, sdlog = sqrt(lnvar))} )); beep(2)
dimnames(N) <- list(Day = seq(dim(N)[1]), Boot = seq(boots))
str(N)

Nage <- setNames(object = lapply(seq(dim(dirich)[1]), function(age) {dirich[age, , ] * N} ), nm = c("Age1", "Age2", "Age3", "Age4"))
str(Nage)

# Age 2.1
days <- which(dat$Date %in% seq(from = as.Date("5/13/14", format = "%m/%d/%y"), to = as.Date("5/30/14", format = "%m/%d/%y"), by = "day"))
AbundancePosteriorKarlukSmolt2014.Age2.1 <- apply(Nage$Age2[days, ], 2, function(boot) {sum(boot)} ) * KarlukSmoltPosteriors$Output$KarlukSmolt2014.Age2.1
hist(AbundancePosteriorKarlukSmolt2014.Age2.1[, 1], col = 8, breaks = 100)
c(Median = median(x = AbundancePosteriorKarlukSmolt2014.Age2.1[, 1]), 
  quantile(x = AbundancePosteriorKarlukSmolt2014.Age2.1[, 1], probs = c(0.05, 0.95)),
  Mean = mean(x = AbundancePosteriorKarlukSmolt2014.Age2.1[, 1]),
  SD = sd(x = AbundancePosteriorKarlukSmolt2014.Age2.1[, 1]))

# Age 2.2
days <- which(dat$Date %in% seq(from = as.Date("5/31/14", format = "%m/%d/%y"), to = as.Date("6/15/14", format = "%m/%d/%y"), by = "day"))
AbundancePosteriorKarlukSmolt2014.Age2.2 <- apply(Nage$Age2[days, ], 2, function(boot) {sum(boot)} ) * KarlukSmoltPosteriors$Output$KarlukSmolt2014.Age2.2
hist(AbundancePosteriorKarlukSmolt2014.Age2.2[, 1], col = 8, breaks = 100)
c(Median = median(x = AbundancePosteriorKarlukSmolt2014.Age2.2[, 1]), 
  quantile(x = AbundancePosteriorKarlukSmolt2014.Age2.2[, 1], probs = c(0.05, 0.95)),
  Mean = mean(x = AbundancePosteriorKarlukSmolt2014.Age2.2[, 1]),
  SD = sd(x = AbundancePosteriorKarlukSmolt2014.Age2.2[, 1]))

# Age 2.3
days <- which(dat$Date %in% seq(from = as.Date("6/15/14", format = "%m/%d/%y"), to = as.Date("7/2/14", format = "%m/%d/%y"), by = "day"))
AbundancePosteriorKarlukSmolt2014.Age2.3 <- apply(Nage$Age2[days, ], 2, function(boot) {sum(boot)} ) * KarlukSmoltPosteriors$Output$KarlukSmolt2014.Age2.3
hist(AbundancePosteriorKarlukSmolt2014.Age2.3[, 1], col = 8, breaks = 100)
c(Median = median(x = AbundancePosteriorKarlukSmolt2014.Age2.3[, 1]), 
  quantile(x = AbundancePosteriorKarlukSmolt2014.Age2.3[, 1], probs = c(0.05, 0.95)),
  Mean = mean(x = AbundancePosteriorKarlukSmolt2014.Age2.3[, 1]),
  SD = sd(x = AbundancePosteriorKarlukSmolt2014.Age2.3[, 1]))

# Age 2 Stratified
KarlukSmolt_Age2_Outmigration_2014 # <- as.numeric(readClipboard()) # These are the three outmigration estimates per strata from the "2014Karluk_smolt_pop_est_by_day_age" worksheet, sheet "Samples"

days <- list(which(dat$Date %in% seq(from = as.Date("5/13/14", format = "%m/%d/%y"), to = as.Date("5/30/14", format = "%m/%d/%y"), by = "day")),
             which(dat$Date %in% seq(from = as.Date("5/31/14", format = "%m/%d/%y"), to = as.Date("6/15/14", format = "%m/%d/%y"), by = "day")),
             which(dat$Date %in% seq(from = as.Date("6/16/14", format = "%m/%d/%y"), to = as.Date("7/2/14", format = "%m/%d/%y"), by = "day")))

posterior <- sapply(days, function(Days) {apply(Nage$Age2[Days, ], 2, function(boot) {sum(boot)} )} )
Karluk2014Age2CVs <- setNames(object = apply(posterior, 2, function(strata) {sd(strata) / mean(strata)} ), nm = KarlukMixtures_Age[8:10])

Karluk2014Age2StratifiedEstimatesWithMRError <- StratifiedEstimator.GCL(groupvec = 1:2, groupnames = KarlukGroups2, maindir = "BAYES/Output", mixvec = KarlukMixtures_Age[8:10], catchvec = KarlukSmolt_Age2_Outmigration_2014,
                                                                        CVvec = Karluk2014Age2CVs, newname = "Karluk2014Age2StratifiedEstimatesWithMRError", priorname = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1); beep(2)

days <- which(dat$Date %in% seq(from = as.Date("5/13/14", format = "%m/%d/%y"), to = as.Date("7/2/14", format = "%m/%d/%y"), by = "day"))
posterior <- apply(Nage$Age2[days, ], 2, function(boot) {sum(boot)} )
Karluk2014Age2StratifiedAbundancesWithMRError <- posterior * Karluk2014Age2StratifiedEstimatesWithMRError$Output
Karluk2014Age2overallCV <- setNames(object = sd(posterior) / mean(posterior), nm = unlist(strsplit(x = KarlukMixtures_Age[8], split = "\\.1")))

hist(Karluk2014Age2StratifiedAbundancesWithMRError[, 1], col = 8, breaks = 100)

c(Median = median(x = Karluk2014Age2StratifiedAbundancesWithMRError[, 1]), 
  quantile(x = Karluk2014Age2StratifiedAbundancesWithMRError[, 1], probs = c(0.05, 0.95)),
  Mean = mean(x = Karluk2014Age2StratifiedAbundancesWithMRError[, 1]),
  SD = sd(x = Karluk2014Age2StratifiedAbundancesWithMRError[, 1]))



# Age 1.12
days <- which(dat$Date %in% seq(from = as.Date("5/13/14", format = "%m/%d/%y"), to = as.Date("6/15/14", format = "%m/%d/%y"), by = "day"))
AbundancePosteriorKarlukSmolt2014.Age1.12 <- apply(Nage$Age1[days, ], 2, function(boot) {sum(boot)} ) * KarlukSmoltPosteriors$Output$KarlukSmolt2014.Age1.12
hist(AbundancePosteriorKarlukSmolt2014.Age1.12[, 1], col = 8, breaks = 100)
c(Median = median(x = AbundancePosteriorKarlukSmolt2014.Age1.12[, 1]), 
  quantile(x = AbundancePosteriorKarlukSmolt2014.Age1.12[, 1], probs = c(0.05, 0.95)),
  Mean = mean(x = AbundancePosteriorKarlukSmolt2014.Age1.12[, 1]),
  SD = sd(x = AbundancePosteriorKarlukSmolt2014.Age1.12[, 1]))

# Age 1.3
days <- which(dat$Date %in% seq(from = as.Date("6/15/14", format = "%m/%d/%y"), to = as.Date("7/2/14", format = "%m/%d/%y"), by = "day"))
AbundancePosteriorKarlukSmolt2014.Age1.3 <- apply(Nage$Age1[days, ], 2, function(boot) {sum(boot)} ) * KarlukSmoltPosteriors$Output$KarlukSmolt2014.Age1.3
hist(AbundancePosteriorKarlukSmolt2014.Age1.3[, 1], col = 8, breaks = 100)
c(Median = median(x = AbundancePosteriorKarlukSmolt2014.Age1.3[, 1]), 
  quantile(x = AbundancePosteriorKarlukSmolt2014.Age1.3[, 1], probs = c(0.05, 0.95)),
  Mean = mean(x = AbundancePosteriorKarlukSmolt2014.Age1.3[, 1]),
  SD = sd(x = AbundancePosteriorKarlukSmolt2014.Age1.3[, 1]))

# Age 1 Stratified
KarlukSmolt_Age1_Outmigration_2014 #<- as.numeric(readClipboard()) # These are the three outmigration estimates per strata from the "2014Karluk_smolt_pop_est_by_day_age" worksheet, sheet "Samples"

days <- list(which(dat$Date %in% seq(from = as.Date("5/13/14", format = "%m/%d/%y"), to = as.Date("6/15/14", format = "%m/%d/%y"), by = "day")),
             which(dat$Date %in% seq(from = as.Date("6/16/14", format = "%m/%d/%y"), to = as.Date("7/2/14", format = "%m/%d/%y"), by = "day")))

posterior <- sapply(days, function(Days) {apply(Nage$Age1[Days, ], 2, function(boot) {sum(boot)} )} )
Karluk2014Age1CVs <- setNames(object = apply(posterior, 2, function(strata) {sd(strata) / mean(strata)} ), nm = KarlukMixtures_Age[6:7])

Karluk2014Age1StratifiedEstimatesWithMRError <- StratifiedEstimator.GCL(groupvec = 1:2, groupnames = KarlukGroups2, maindir = "BAYES/Output", mixvec = KarlukMixtures_Age[6:7], catchvec = KarlukSmolt_Age1_Outmigration_2014,
                                                                        CVvec = Karluk2014Age1CVs, newname = "Karluk2014Age1StratifiedEstimatesWithMRError", priorname = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1); beep(2)

days <- which(dat$Date %in% seq(from = as.Date("5/13/14", format = "%m/%d/%y"), to = as.Date("7/2/14", format = "%m/%d/%y"), by = "day"))
posterior <- apply(Nage$Age1[days, ], 2, function(boot) {sum(boot)} )
Karluk2014Age1StratifiedAbundancesWithMRError <- posterior * Karluk2014Age1StratifiedEstimatesWithMRError$Output
Karluk2014Age1overallCV <- setNames(object = sd(posterior) / mean(posterior), nm = unlist(strsplit(x = KarlukMixtures_Age[6], split = "\\.12")))

hist(Karluk2014Age1StratifiedAbundancesWithMRError[, 1], col = 8, breaks = 100)

c(Median = median(x = Karluk2014Age1StratifiedAbundancesWithMRError[, 1]), 
  quantile(x = Karluk2014Age1StratifiedAbundancesWithMRError[, 1], probs = c(0.05, 0.95)),
  Mean = mean(x = Karluk2014Age1StratifiedAbundancesWithMRError[, 1]),
  SD = sd(x = Karluk2014Age1StratifiedAbundancesWithMRError[, 1]))







#~~~~~~~~~~~~~~~~~~
## 2013
dat <- read.xlsx(file = "2013Karluk_smolt_pop_est_by_day_age.xlsx", sheetName = "AgeSpecificOM")
str(dat)

boots = 100000

dirich <- replicate(n = boots, expr = apply(dat[, c("Age1", "Age2", "Age3", "Age4")], 1, function(day) {rdirich(alpha = day)} )); beep(2)
dimnames(dirich) <- list(Ages = c("Age1", "Age2", "Age3", "Age4"), Day = seq(dim(dirich)[2]), Boot = seq(boots))
str(dirich)

N <- t(apply(dat[, c("DailyOM", "CV")], 1, function(day) {lnvar = log(day[2]^2 + 1)
                                                          lnmean = log(day[1]) - lnvar/2
                                                          rlnorm(n = boots, meanlog = lnmean, sdlog = sqrt(lnvar))} )); beep(2)
dimnames(N) <- list(Day = seq(dim(N)[1]), Boot = seq(boots))
str(N)

Nage <- setNames(object = lapply(seq(dim(dirich)[1]), function(age) {dirich[age, , ] * N} ), nm = c("Age1", "Age2", "Age3", "Age4"))
str(Nage)

# Age 2.1
days <- which(dat$Date %in% seq(from = as.Date("5/16/13", format = "%m/%d/%y"), to = as.Date("5/29/13", format = "%m/%d/%y"), by = "day"))
AbundancePosteriorKarlukSmolt2013.Age2.1 <- apply(Nage$Age2[days, ], 2, function(boot) {sum(boot)} ) * KarlukSmoltPosteriors$Output$KarlukSmolt2013.Age2.1
hist(AbundancePosteriorKarlukSmolt2013.Age2.1[, 1], col = 8, breaks = 100)
c(Median = median(x = AbundancePosteriorKarlukSmolt2013.Age2.1[, 1]), 
  quantile(x = AbundancePosteriorKarlukSmolt2013.Age2.1[, 1], probs = c(0.05, 0.95)),
  Mean = mean(x = AbundancePosteriorKarlukSmolt2013.Age2.1[, 1]),
  SD = sd(x = AbundancePosteriorKarlukSmolt2013.Age2.1[, 1]))

# Age 2.2
days <- which(dat$Date %in% seq(from = as.Date("5/30/13", format = "%m/%d/%y"), to = as.Date("6/10/13", format = "%m/%d/%y"), by = "day"))
AbundancePosteriorKarlukSmolt2013.Age2.2 <- apply(Nage$Age2[days, ], 2, function(boot) {sum(boot)} ) * KarlukSmoltPosteriors$Output$KarlukSmolt2013.Age2.2
hist(AbundancePosteriorKarlukSmolt2013.Age2.2[, 1], col = 8, breaks = 100)
c(Median = median(x = AbundancePosteriorKarlukSmolt2013.Age2.2[, 1]), 
  quantile(x = AbundancePosteriorKarlukSmolt2013.Age2.2[, 1], probs = c(0.05, 0.95)),
  Mean = mean(x = AbundancePosteriorKarlukSmolt2013.Age2.2[, 1]),
  SD = sd(x = AbundancePosteriorKarlukSmolt2013.Age2.2[, 1]))

# Age 2.3
days <- which(dat$Date %in% seq(from = as.Date("6/11/13", format = "%m/%d/%y"), to = as.Date("6/24/13", format = "%m/%d/%y"), by = "day"))
AbundancePosteriorKarlukSmolt2013.Age2.3 <- apply(Nage$Age2[days, ], 2, function(boot) {sum(boot)} ) * KarlukSmoltPosteriors$Output$KarlukSmolt2013.Age2.3
hist(AbundancePosteriorKarlukSmolt2013.Age2.3[, 1], col = 8, breaks = 100)
c(Median = median(x = AbundancePosteriorKarlukSmolt2013.Age2.3[, 1]), 
  quantile(x = AbundancePosteriorKarlukSmolt2013.Age2.3[, 1], probs = c(0.05, 0.95)),
  Mean = mean(x = AbundancePosteriorKarlukSmolt2013.Age2.3[, 1]),
  SD = sd(x = AbundancePosteriorKarlukSmolt2013.Age2.3[, 1]))

# Age 2 Stratified
KarlukSmolt_Age2_Outmigration_2013 #<- as.numeric(readClipboard()) # These are the three outmigration estimates per strata from the "2013Karluk_smolt_pop_est_by_day_age" worksheet, sheet "Samples"

days <- list(which(dat$Date %in% seq(from = as.Date("5/16/13", format = "%m/%d/%y"), to = as.Date("5/29/13", format = "%m/%d/%y"), by = "day")),
             which(dat$Date %in% seq(from = as.Date("5/30/13", format = "%m/%d/%y"), to = as.Date("6/10/13", format = "%m/%d/%y"), by = "day")),
             which(dat$Date %in% seq(from = as.Date("6/11/13", format = "%m/%d/%y"), to = as.Date("6/24/13", format = "%m/%d/%y"), by = "day")))

posterior <- sapply(days, function(Days) {apply(Nage$Age2[Days, ], 2, function(boot) {sum(boot)} )} )
Karluk2013Age2CVs <- setNames(object = apply(posterior, 2, function(strata) {sd(strata) / mean(strata)} ), nm = KarlukMixtures_Age[2:4])

Karluk2013Age2StratifiedEstimatesWithMRError <- StratifiedEstimator.GCL(groupvec = 1:2, groupnames = KarlukGroups2, maindir = "BAYES/Output", mixvec = KarlukMixtures_Age[2:4], catchvec = KarlukSmolt_Age2_Outmigration_2013,
                                                                        CVvec = Karluk2013Age2CVs, newname = "Karluk2013Age2StratifiedEstimatesWithMRError", priorname = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1); beep(2)

days <- which(dat$Date %in% seq(from = as.Date("5/16/13", format = "%m/%d/%y"), to = as.Date("6/24/13", format = "%m/%d/%y"), by = "day"))
posterior <- apply(Nage$Age2[days, ], 2, function(boot) {sum(boot)} )
Karluk2013Age2StratifiedAbundancesWithMRError <- posterior * Karluk2013Age2StratifiedEstimatesWithMRError$Output
Karluk2013Age2overallCV <- setNames(object = sd(posterior) / mean(posterior), nm = unlist(strsplit(x = KarlukMixtures_Age[2], split = "\\.1")))

hist(Karluk2013Age2StratifiedAbundancesWithMRError[, 1], col = 8, breaks = 100)

c(Median = median(x = Karluk2013Age2StratifiedAbundancesWithMRError[, 1]), 
  quantile(x = Karluk2013Age2StratifiedAbundancesWithMRError[, 1], probs = c(0.05, 0.95)),
  Mean = mean(x = Karluk2013Age2StratifiedAbundancesWithMRError[, 1]),
  SD = sd(x = Karluk2013Age2StratifiedAbundancesWithMRError[, 1]))


# Age 1
KarlukSmolt_Age1_Outmigration_2013

KarlukSmoltPosteriors$Stats$KarlukSmolt2013.Age1

days <- which(dat$Date %in% seq(from = as.Date("5/16/13", format = "%m/%d/%y"), to = as.Date("6/24/13", format = "%m/%d/%y"), by = "day"))

posterior <- apply(Nage$Age1[days, ], 2, function(boot) {sum(boot)} )
Karluk2013Age1overallCV <- setNames(object = sd(posterior) / mean(posterior), nm = KarlukMixtures_Age[1])

Karluk2013Age1AbundancesWithMRError <- posterior * KarlukSmoltPosteriors$Output$KarlukSmolt2013.Age1
colnames(Karluk2013Age1AbundancesWithMRError) <- KarlukGroups2

hist(Karluk2013Age1AbundancesWithMRError[, 1], col = 8, breaks = 100)

c(Median = median(x = Karluk2013Age1AbundancesWithMRError[, 1]), 
  quantile(x = Karluk2013Age1AbundancesWithMRError[, 1], probs = c(0.05, 0.95)),
  Mean = mean(x = Karluk2013Age1AbundancesWithMRError[, 1]),
  SD = sd(x = Karluk2013Age1AbundancesWithMRError[, 1]))



# Age 3
KarlukSmolt_Age3_Outmigration_2013

KarlukSmoltPosteriors$Stats$KarlukSmolt2013.Age3

days <- which(dat$Date %in% seq(from = as.Date("5/16/13", format = "%m/%d/%y"), to = as.Date("6/24/13", format = "%m/%d/%y"), by = "day"))

posterior <- apply(Nage$Age3[days, ], 2, function(boot) {sum(boot)} )
Karluk2013Age3overallCV <- setNames(object = sd(posterior) / mean(posterior), nm = KarlukMixtures_Age[5])

Karluk2013Age3AbundancesWithMRError <- posterior * KarlukSmoltPosteriors$Output$KarlukSmolt2013.Age3
colnames(Karluk2013Age3AbundancesWithMRError) <- KarlukGroups2

hist(Karluk2013Age3AbundancesWithMRError[, 1], col = 8, breaks = 100)

c(Median = median(x = Karluk2013Age3AbundancesWithMRError[, 1]), 
  quantile(x = Karluk2013Age3AbundancesWithMRError[, 1], probs = c(0.05, 0.95)),
  Mean = mean(x = Karluk2013Age3AbundancesWithMRError[, 1]),
  SD = sd(x = Karluk2013Age3AbundancesWithMRError[, 1]))




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Make Tables!!! ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Table 4, 2013 Smolt by Age
sheetname = "Karluk 2013 Smolt by Age"
tablenum = 4
Nrows <- 15
Ncols <- 15
yr = 2013
byAge = TRUE
Ages = FALSE
ages = "Age 2"

dates.vec <- KarlukShortDates[c(4,4,8)]
samplesizes.vec <- c(KarlukFinalSampleSizes[4], sum(KarlukFinalSampleSizes[5:7]), KarlukFinalSampleSizes[8])

names(KarlukSmoltPosteriors$Stats)
colnames(Karluk2013Age2StratifiedEstimatesWithMRError[[1]]) <- c("mean", "sd", "5%", "median", "95%", "P=0", "GR")
stats.proportions.list <- c(KarlukSmoltPosteriors$Stats[10], list(KarlukSmolt2013.Age2Stratified = Karluk2013Age2StratifiedEstimatesWithMRError[[1]][, colnames(KarlukSmoltPosteriors$Stats[[10]])]), KarlukSmoltPosteriors$Stats[14])  # List of proportion stats

stats.abundances.list <- lapply(list(Karluk2013Age1AbundancesWithMRError, Karluk2013Age2StratifiedAbundancesWithMRError, Karluk2013Age3AbundancesWithMRError), function(Posterior) {
  t(apply(Posterior, 2, function(RG) {c(median = median(x = RG), quantile(x = RG, probs = c(0.05, 0.95)), mean = mean(x = RG), sd = sd(x = RG)) })) })  # List of abundance stats

abundances.vec #<- setNames(object = as.numeric(readClipboard()), nm = paste("Age", 1:3))  # Vector of abundances
CV.vec <- c(Karluk2013Age1overallCV, Karluk2013Age2overallCV, Karluk2013Age3overallCV)  # Vector of CVs

#~~~~~~~~
KarlukRedoTable <- matrix(data = "", nrow = Nrows, ncol = Ncols)

KarlukRedoTable[1, 1] <- paste("Table ", tablenum,".-Estimates of stock composition and stock-specific outmigration for Karluk River ", ifelse(Ages, paste(ages, " ", sep = ''), ''),
                               "sockeye salmon smolt ", ifelse(byAge, "by age", "by stratum"), ", ", yr,
                               ".  Reporting group-specific stock composition estimates (%) include median, 90% credibility interval, the probability that the group estimate is equal to zero (P=0), mean, and SD.  Stock-specific estimates of outmigration are based upon mark-recapture estimates and variances (CV) of outmigration for each stratum. See text for details.", sep = "")
KarlukRedoTable[Nrows, 1] <- "Note: Stock composition estimates may not sum to 100% and stock-specific outmigration estimates may not sum to the total outmigration due to rounding error."

KarlukRedoTable[2, ] <- c("Stratum","",'','Composition (%)','','','','','','',"Outmigration (number of fish)",'','','','')
KarlukRedoTable[3, ] <- c(ifelse(byAge, 'Age', 'Period'),'Sample Size','Reporting','',"90% CI",'','','','','','',"90% CI",'','','')
KarlukRedoTable[4, ] <- c("Dates",'CV','Group',"Median","5%","95%","P=0","Mean","SD",'',"Median","5%","95%","Mean","SD")
KarlukRedoTable[c(5, 8, 11, 14), 1] <- c("Age 1", "Age 2", "Age 3", '')

KarlukRedoTable[c(6, 9, 12), 1] <- dates.vec
KarlukRedoTable[c(5:6, 8:9, 11:12), 3] <- KarlukGroups

KarlukRedoTable[c(5, 8, 11), 2] <- paste("n=", samplesizes.vec, sep = '')
KarlukRedoTable[c(6, 9, 12), 2] <- paste("CV=", formatC(x = CV.vec * 100, digits = 1, format = "f"), "%", sep = '')

KarlukRedoTable[c(5:6, 8:9, 11:12), c(4:6, 8:9)] <- formatC(x = do.call(what = "rbind", args = stats.proportions.list[1:3])[, c("median", "5%", "95%", "mean", "sd")] * 100, format = "f", digits = 1)
KarlukRedoTable[c(5:6, 8:9, 11:12), 7] <- formatC(x = do.call(what = "rbind", args = stats.proportions.list[1:3])[, c("P=0")] * 100, format = "f", digits = 2)
KarlukRedoTable[c(5:6, 8:9, 11:12), 11:15] <- formatC(x = do.call(what = "rbind", args = stats.abundances.list)[, c("median", "5%", "95%", "mean", "sd")], digits = 0, format = "f", big.mark = ",")

KarlukRedoTable[c(7, 10, 13, 14), 12] <- paste(c(paste("Age", 1:3), yr), "Total", sep = ' ')
KarlukRedoTable[c(7, 10, 13, 14), 14] <- formatC(x = c(abundances.vec, sum(abundances.vec)), digits = 0, format = "f", big.mark = ",")

require(xlsx)
write.xlsx(x = KarlukRedoTable, file = "Tables/KarlukSmolt2013-2014TablesNEW.xlsx", sheetName = sheetname, col.names = FALSE, row.names = FALSE, append = TRUE)
#~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~
## Table 5, 2014 Smolt by Age
sheetname = "Karluk 2014 Smolt by Age"
tablenum = 5
Nrows <- 12
Ncols <- 15
yr = 2014
byAge = TRUE
Ages = FALSE
ages = "Age 2"

dates.vec <- rep(paste(unlist(strsplit(x = KarlukShortDates[12], split = "-"))[1],
                       ifelse(length(grep(pattern = "/", x = unlist(strsplit(x = KarlukShortDates[13], split = "-"))[2])) == 0, 
                              paste(unlist(strsplit(x = unlist(strsplit(x = KarlukShortDates[13], split = "-"))[1], split = "/"))[1],
                                    unlist(strsplit(x = KarlukShortDates[13], split = "-"))[2],
                                    sep = "/"), 
                              unlist(strsplit(x = KarlukShortDates[13], split = "-"))[2]),
                       sep = '-'),
                 2)
samplesizes.vec <- c(sum(KarlukFinalSampleSizes[12:13]), sum(KarlukFinalSampleSizes[14:16]))


stats.proportions.list <- list(KarlukSmolt2014.Age1Stratified = Karluk2014Age1StratifiedEstimatesWithMRError$Summary, KarlukSmolt2014.Age2Stratified = Karluk2014Age2StratifiedEstimatesWithMRError$Summary)
colnames(stats.proportions.list[[1]]) <- c("mean", "sd", "5%", "median", "95%", "P=0", "GR")
colnames(stats.proportions.list[[2]]) <- c("mean", "sd", "5%", "median", "95%", "P=0", "GR")

stats.abundances.list <- setNames(lapply(list(Karluk2014Age1StratifiedAbundancesWithMRError, Karluk2014Age2StratifiedAbundancesWithMRError), function(Posterior) {
  t(apply(Posterior, 2, function(RG) {c(median = median(x = RG), quantile(x = RG, probs = c(0.05, 0.95)), mean = mean(x = RG), sd = sd(x = RG)) })) }), nm = c("Age 1", "Age 2"))  # List of abundance stats

abundances.vec #<- setNames(object = as.numeric(readClipboard()), nm = paste("Age", 1:2))  # Vector of abundances
CV.vec <- c(Karluk2014Age1overallCV, Karluk2014Age2overallCV)  # Vector of CVs

#~~~~~~~~
KarlukRedoTable <- matrix(data = "", nrow = Nrows, ncol = Ncols)

KarlukRedoTable[1, 1] <- paste("Table ", tablenum,".-Estimates of stock composition and stock-specific outmigration for Karluk River ", ifelse(Ages, paste(ages, " ", sep = ''), ''),
                               "sockeye salmon smolt ", ifelse(byAge, "by age", "by stratum"), ", ", yr,
                               ".  Reporting group-specific stock composition estimates (%) include median, 90% credibility interval, the probability that the group estimate is equal to zero (P=0), mean, and SD.  Stock-specific estimates of outmigration are based upon mark-recapture estimates and variances (CV) of outmigration for each stratum. See text for details.", sep = "")
KarlukRedoTable[Nrows, 1] <- "Note: Stock composition estimates may not sum to 100% and stock-specific outmigration estimates may not sum to the total outmigration due to rounding error."

KarlukRedoTable[2, ] <- c("Stratum","",'','Composition (%)','','','','','','',"Outmigration (number of fish)",'','','','')
KarlukRedoTable[3, ] <- c(ifelse(byAge, 'Age', 'Period'),'Sample Size','Reporting','',"90% CI",'','','','','','',"90% CI",'','','')
KarlukRedoTable[4, ] <- c("Dates",'CV','Group',"Median","5%","95%","P=0","Mean","SD",'',"Median","5%","95%","Mean","SD")
KarlukRedoTable[c(5, 8), 1] <- c("Age 1", "Age 2")

KarlukRedoTable[c(6, 9), 1] <- dates.vec
KarlukRedoTable[c(5:6, 8:9), 3] <- KarlukGroups

KarlukRedoTable[c(5, 8), 2] <- paste("n=", samplesizes.vec, sep = '')
KarlukRedoTable[c(6, 9), 2] <- paste("CV=", formatC(x = CV.vec * 100, digits = 1, format = "f"), "%", sep = '')

KarlukRedoTable[c(5:6, 8:9), c(4:6, 8:9)] <- formatC(x = do.call(what = "rbind", args = stats.proportions.list[1:2])[, c("median", "5%", "95%", "mean", "sd")] * 100, format = "f", digits = 1)
KarlukRedoTable[c(5:6, 8:9), 7] <- formatC(x = do.call(what = "rbind", args = stats.proportions.list[1:2])[, c("P=0")] * 100, format = "f", digits = 2)
KarlukRedoTable[c(5:6, 8:9), 11:15] <- formatC(x = do.call(what = "rbind", args = stats.abundances.list)[, c("median", "5%", "95%", "mean", "sd")], digits = 0, format = "f", big.mark = ",")

KarlukRedoTable[c(7, 10, 11), 12] <- paste(c(paste("Age", 1:2), yr), "Total", sep = ' ')
KarlukRedoTable[c(7, 10, 11), 14] <- formatC(x = c(abundances.vec, sum(abundances.vec)), digits = 0, format = "f", big.mark = ",")

require(xlsx)
write.xlsx(x = KarlukRedoTable, file = "Tables/KarlukSmolt2013-2014TablesNEW.xlsx", sheetName = sheetname, col.names = FALSE, row.names = FALSE, append = TRUE)
#~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~
## Table 6, 2013 Smolt Age 2
sheetname = "Karluk 2013 Smolt Age 2"
tablenum = 6
Nrows <- 17
Ncols <- 15
yr = 2013
byAge = FALSE
Ages = TRUE
ages = "Age 2"

dates.vec <- KarlukShortDates[5:7]
samplesizes.vec <- KarlukFinalSampleSizes[5:7]

names(KarlukSmoltPosteriors$Stats)
stats.proportions.list <- KarlukSmoltPosteriors$Stats[11:13]  # List of proportion stats
stratified.proportions <- Karluk2013Age2StratifiedEstimatesWithMRError$Summary  # Stratified porportion stats
stats.abundances.list <-  setNames(lapply(list(AbundancePosteriorKarlukSmolt2013.Age2.1, AbundancePosteriorKarlukSmolt2013.Age2.2, AbundancePosteriorKarlukSmolt2013.Age2.3), function(Posterior) {
  t(apply(Posterior, 2, function(RG) {c(median = median(x = RG), quantile(x = RG, probs = c(0.05, 0.95)), mean = mean(x = RG), sd = sd(x = RG)) })) }), nm = paste("KarlukSmolt2013.Age2.", 1:3, sep = ''))  # List of abundance stats
stratified.abundances <- t(apply(Karluk2013Age2StratifiedAbundancesWithMRError, 2, function(RG) {c(median = median(x = RG), quantile(x = RG, probs = c(0.05, 0.95)), mean = mean(x = RG), sd = sd(x = RG)) }))  # Stratified abundance stats
abundances.vec <- KarlukOutmigration[5:7]  # Vector of abundances
CV.vec <- Karluk2013Age2CVs  # Vector of CVs

#~~~~~~~~
KarlukRedoTable <- matrix(data = "", nrow = Nrows, ncol = Ncols)

KarlukRedoTable[1, 1] <- paste("Table ", tablenum,".-Estimates of stock composition and stock-specific outmigration for Karluk River ", ifelse(Ages, paste(ages, " ", sep = ''), ''),
                               "sockeye salmon smolt ", ifelse(byAge, "by age", "by stratum"), ", ", yr,
                               ".  Reporting group-specific stock composition estimates (%) include median, 90% credibility interval, the probability that the group estimate is equal to zero (P=0), mean, and SD.  Stock-specific estimates of outmigration are based upon mark-recapture estimates and variances (CV) of outmigration for each stratum. See text for details.", sep = "")
KarlukRedoTable[Nrows, 1] <- "Note: Stock composition estimates may not sum to 100% and stock-specific outmigration estimates may not sum to the total outmigration due to rounding error."

KarlukRedoTable[2, ] <- c("Stratum","",'','Composition (%)','','','','','','',"Outmigration (number of fish)",'','','','')
KarlukRedoTable[3, ] <- c('Period','Sample Size','Reporting','',"90% CI",'','','','','','',"90% CI",'','','')
KarlukRedoTable[4, ] <- c("Dates",'CV','Group',"Median","5%","95%","P=0","Mean","SD",'',"Median","5%","95%","Mean","SD")
KarlukRedoTable[c(5, 8, 11, 14), 1] <- c("Early", "Middle", "Late", yr)

KarlukRedoTable[c(6, 9, 12, 15), 1] <- c(dates.vec, paste(unlist(strsplit(x = dates.vec[1], split = "-"))[1],
                                                          ifelse(length(grep(pattern = "/", x = unlist(strsplit(x = dates.vec[3], split = "-"))[2])) == 0, 
                                                                 paste(unlist(strsplit(x = unlist(strsplit(x = dates.vec[3], split = "-"))[1], split = "/"))[1],
                                                                       unlist(strsplit(x = dates.vec[3], split = "-"))[2],
                                                                       sep = "/"), 
                                                                 unlist(strsplit(x = dates.vec[3], split = "-"))[2]),
                                                          sep = '-'))
KarlukRedoTable[c(5:6, 8:9, 11:12, 14:15), 3] <- KarlukGroups

KarlukRedoTable[c(5, 8, 11, 14), 2] <- paste("n=", c(samplesizes.vec, sum(samplesizes.vec)), sep = '')
KarlukRedoTable[c(6, 9, 12), 2] <- paste("CV=", formatC(x = CV.vec * 100, digits = 1, format = "f"), "%", sep = '')

KarlukRedoTable[c(5:6, 8:9, 11:12), c(4:6, 8:9)] <- formatC(x = do.call(what = "rbind", args = stats.proportions.list[1:3])[, c("median", "5%", "95%", "mean", "sd")] * 100, format = "f", digits = 1)
KarlukRedoTable[c(5:6, 8:9, 11:12), 7] <- formatC(x = do.call(what = "rbind", args = stats.proportions.list[1:3])[, c("P=0")] * 100, format = "f", digits = 2)
KarlukRedoTable[c(5:6, 8:9, 11:12), 11:15] <- formatC(x = do.call(what = "rbind", args = stats.abundances.list)[, c("median", "5%", "95%", "mean", "sd")], digits = 0, format = "f", big.mark = ",")

KarlukRedoTable[c(7, 10, 13, 16), 12] <- paste(c("Early", "Middle", "Late", yr), "Total", sep = ' ')
KarlukRedoTable[c(7, 10, 13, 16), 14] <- formatC(x = c(abundances.vec, sum(abundances.vec)), digits = 0, format = "f", big.mark = ",")

KarlukRedoTable[14:15, c(4:6, 8:9)] <- formatC(x = stratified.proportions[, c("median", "5%", "95%", "mean", "sd")] * 100, format = "f", digits = 1)
KarlukRedoTable[14:15, 7] <- formatC(x = stratified.proportions[, c("P=0")] * 100, format = "f", digits = 2)
KarlukRedoTable[14:15, 11:15] <- formatC(x = stratified.abundances[, c("median", "5%", "95%", "mean", "sd")], format = "f", digits = 0, big.mark = ",")


require(xlsx)
write.xlsx(x = KarlukRedoTable, file = "Tables/KarlukSmolt2013-2014TablesNEW.xlsx", sheetName = sheetname, col.names = FALSE, row.names = FALSE, append = TRUE)
#~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~
## Table 7, 2014 Smolt Age 1
sheetname = "Karluk 2014 Smolt Age 1"
tablenum = 7
Nrows <- 14
Ncols <- 15
yr = 2014
byAge = FALSE
Ages = TRUE
ages = "Age 1"

dates.vec <- KarlukShortDates[12:13]
samplesizes.vec <- KarlukFinalSampleSizes[12:13]

names(KarlukSmoltPosteriors$Stats)
stats.proportions.list <- KarlukSmoltPosteriors$Stats[15:16]  # List of proportion stats
stratified.proportions <- Karluk2014Age1StratifiedEstimatesWithMRError$Summary  # Stratified porportion stats
colnames(stratified.proportions) <- c("mean", "sd", "5%", "median", "95%", "P=0", "GR")
stats.abundances.list <-  setNames(lapply(list(AbundancePosteriorKarlukSmolt2014.Age1.12, AbundancePosteriorKarlukSmolt2014.Age1.3), function(Posterior) {
  t(apply(Posterior, 2, function(RG) {c(median = median(x = RG), quantile(x = RG, probs = c(0.05, 0.95)), mean = mean(x = RG), sd = sd(x = RG)) })) }), nm = paste("KarlukSmolt2014.Age1.", c("12", "3"), sep = ''))  # List of abundance stats
stratified.abundances <- t(apply(Karluk2014Age1StratifiedAbundancesWithMRError, 2, function(RG) {c(median = median(x = RG), quantile(x = RG, probs = c(0.05, 0.95)), mean = mean(x = RG), sd = sd(x = RG)) }))  # Stratified abundance stats
abundances.vec <- KarlukOutmigration[12:13]  # Vector of abundances
CV.vec <- Karluk2014Age1CVs  # Vector of CVs

#~~~~~~~~
KarlukRedoTable <- matrix(data = "", nrow = Nrows, ncol = Ncols)

KarlukRedoTable[1, 1] <- paste("Table ", tablenum,".-Estimates of stock composition and stock-specific outmigration for Karluk River ", ifelse(Ages, paste(ages, " ", sep = ''), ''),
                               "sockeye salmon smolt ", ifelse(byAge, "by age", "by stratum"), ", ", yr,
                               ".  Reporting group-specific stock composition estimates (%) include median, 90% credibility interval, the probability that the group estimate is equal to zero (P=0), mean, and SD.  Stock-specific estimates of outmigration are based upon mark-recapture estimates and variances (CV) of outmigration for each stratum. See text for details.", sep = "")
KarlukRedoTable[Nrows, 1] <- "Note: Stock composition estimates may not sum to 100% and stock-specific outmigration estimates may not sum to the total outmigration due to rounding error."

KarlukRedoTable[2, ] <- c("Stratum","",'','Composition (%)','','','','','','',"Outmigration (number of fish)",'','','','')
KarlukRedoTable[3, ] <- c('Period','Sample Size','Reporting','',"90% CI",'','','','','','',"90% CI",'','','')
KarlukRedoTable[4, ] <- c("Dates",'CV','Group',"Median","5%","95%","P=0","Mean","SD",'',"Median","5%","95%","Mean","SD")
KarlukRedoTable[c(5, 8, 11), 1] <- c("Early", "Late", yr)

KarlukRedoTable[c(6, 9, 12), 1] <- c(dates.vec, paste(unlist(strsplit(x = dates.vec[1], split = "-"))[1],
                                                          ifelse(length(grep(pattern = "/", x = unlist(strsplit(x = dates.vec[2], split = "-"))[2])) == 0, 
                                                                 paste(unlist(strsplit(x = unlist(strsplit(x = dates.vec[2], split = "-"))[1], split = "/"))[1],
                                                                       unlist(strsplit(x = dates.vec[2], split = "-"))[2],
                                                                       sep = "/"), 
                                                                 unlist(strsplit(x = dates.vec[2], split = "-"))[2]),
                                                          sep = '-'))
KarlukRedoTable[c(5:6, 8:9, 11:12), 3] <- KarlukGroups

KarlukRedoTable[c(5, 8, 11), 2] <- paste("n=", c(samplesizes.vec, sum(samplesizes.vec)), sep = '')
KarlukRedoTable[c(6, 9), 2] <- paste("CV=", formatC(x = CV.vec * 100, digits = 1, format = "f"), "%", sep = '')

KarlukRedoTable[c(5:6, 8:9), c(4:6, 8:9)] <- formatC(x = do.call(what = "rbind", args = stats.proportions.list[1:2])[, c("median", "5%", "95%", "mean", "sd")] * 100, format = "f", digits = 1)
KarlukRedoTable[c(5:6, 8:9), 7] <- formatC(x = do.call(what = "rbind", args = stats.proportions.list[1:2])[, c("P=0")] * 100, format = "f", digits = 2)
KarlukRedoTable[c(5:6, 8:9), 11:15] <- formatC(x = do.call(what = "rbind", args = stats.abundances.list)[, c("median", "5%", "95%", "mean", "sd")], digits = 0, format = "f", big.mark = ",")

KarlukRedoTable[c(7, 10, 13), 12] <- paste(c("Early", "Late", yr), "Total", sep = ' ')
KarlukRedoTable[c(7, 10, 13), 14] <- formatC(x = c(abundances.vec, sum(abundances.vec)), digits = 0, format = "f", big.mark = ",")

KarlukRedoTable[11:12, c(4:6, 8:9)] <- formatC(x = stratified.proportions[, c("median", "5%", "95%", "mean", "sd")] * 100, format = "f", digits = 1)
KarlukRedoTable[11:12, 7] <- formatC(x = stratified.proportions[, c("P=0")] * 100, format = "f", digits = 2)
KarlukRedoTable[11:12, 11:15] <- formatC(x = stratified.abundances[, c("median", "5%", "95%", "mean", "sd")], format = "f", digits = 0, big.mark = ",")


require(xlsx)
write.xlsx(x = KarlukRedoTable, file = "Tables/KarlukSmolt2013-2014TablesNEW.xlsx", sheetName = sheetname, col.names = FALSE, row.names = FALSE, append = TRUE)
#~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~
## Table 8, 2014 Smolt Age 2
sheetname = "Karluk 2014 Smolt Age 2"
tablenum = 8
Nrows <- 17
Ncols <- 15
yr = 2014
byAge = FALSE
Ages = TRUE
ages = "Age 2"

dates.vec <- KarlukShortDates[14:16]
samplesizes.vec <- KarlukFinalSampleSizes[14:16]

names(KarlukSmoltPosteriors$Stats)
stats.proportions.list <- KarlukSmoltPosteriors$Stats[17:19]  # List of proportion stats
stratified.proportions <- Karluk2014Age2StratifiedEstimatesWithMRError$Summary  # Stratified porportion stats
colnames(stratified.proportions) <- c("mean", "sd", "5%", "median", "95%", "P=0", "GR")
stats.abundances.list <-  setNames(lapply(list(AbundancePosteriorKarlukSmolt2014.Age2.1, AbundancePosteriorKarlukSmolt2014.Age2.2, AbundancePosteriorKarlukSmolt2014.Age2.3), function(Posterior) {
  t(apply(Posterior, 2, function(RG) {c(median = median(x = RG), quantile(x = RG, probs = c(0.05, 0.95)), mean = mean(x = RG), sd = sd(x = RG)) })) }), nm = paste("KarlukSmolt2014.Age2.", 1:3, sep = ''))  # List of abundance stats
stratified.abundances <- t(apply(Karluk2014Age2StratifiedAbundancesWithMRError, 2, function(RG) {c(median = median(x = RG), quantile(x = RG, probs = c(0.05, 0.95)), mean = mean(x = RG), sd = sd(x = RG)) }))  # Stratified abundance stats
abundances.vec <- KarlukOutmigration[14:16]  # Vector of abundances
CV.vec <- Karluk2014Age2CVs  # Vector of CVs

#~~~~~~~~
KarlukRedoTable <- matrix(data = "", nrow = Nrows, ncol = Ncols)

KarlukRedoTable[1, 1] <- paste("Table ", tablenum,".-Estimates of stock composition and stock-specific outmigration for Karluk River ", ifelse(Ages, paste(ages, " ", sep = ''), ''),
                               "sockeye salmon smolt ", ifelse(byAge, "by age", "by stratum"), ", ", yr,
                               ".  Reporting group-specific stock composition estimates (%) include median, 90% credibility interval, the probability that the group estimate is equal to zero (P=0), mean, and SD.  Stock-specific estimates of outmigration are based upon mark-recapture estimates and variances (CV) of outmigration for each stratum. See text for details.", sep = "")
KarlukRedoTable[Nrows, 1] <- "Note: Stock composition estimates may not sum to 100% and stock-specific outmigration estimates may not sum to the total outmigration due to rounding error."

KarlukRedoTable[2, ] <- c("Stratum","",'','Composition (%)','','','','','','',"Outmigration (number of fish)",'','','','')
KarlukRedoTable[3, ] <- c('Period','Sample Size','Reporting','',"90% CI",'','','','','','',"90% CI",'','','')
KarlukRedoTable[4, ] <- c("Dates",'CV','Group',"Median","5%","95%","P=0","Mean","SD",'',"Median","5%","95%","Mean","SD")
KarlukRedoTable[c(5, 8, 11, 14), 1] <- c("Early", "Middle", "Late", yr)

KarlukRedoTable[c(6, 9, 12, 15), 1] <- c(dates.vec, paste(unlist(strsplit(x = dates.vec[1], split = "-"))[1],
                                                          ifelse(length(grep(pattern = "/", x = unlist(strsplit(x = dates.vec[3], split = "-"))[2])) == 0, 
                                                                 paste(unlist(strsplit(x = unlist(strsplit(x = dates.vec[3], split = "-"))[1], split = "/"))[1],
                                                                       unlist(strsplit(x = dates.vec[3], split = "-"))[2],
                                                                       sep = "/"), 
                                                                 unlist(strsplit(x = dates.vec[3], split = "-"))[2]),
                                                          sep = '-'))
KarlukRedoTable[c(5:6, 8:9, 11:12, 14:15), 3] <- KarlukGroups

KarlukRedoTable[c(5, 8, 11, 14), 2] <- paste("n=", c(samplesizes.vec, sum(samplesizes.vec)), sep = '')
KarlukRedoTable[c(6, 9, 12), 2] <- paste("CV=", formatC(x = CV.vec * 100, digits = 1, format = "f"), "%", sep = '')

KarlukRedoTable[c(5:6, 8:9, 11:12), c(4:6, 8:9)] <- formatC(x = do.call(what = "rbind", args = stats.proportions.list[1:3])[, c("median", "5%", "95%", "mean", "sd")] * 100, format = "f", digits = 1)
KarlukRedoTable[c(5:6, 8:9, 11:12), 7] <- formatC(x = do.call(what = "rbind", args = stats.proportions.list[1:3])[, c("P=0")] * 100, format = "f", digits = 2)
KarlukRedoTable[c(5:6, 8:9, 11:12), 11:15] <- formatC(x = do.call(what = "rbind", args = stats.abundances.list)[, c("median", "5%", "95%", "mean", "sd")], digits = 0, format = "f", big.mark = ",")

KarlukRedoTable[c(7, 10, 13, 16), 12] <- paste(c("Early", "Middle", "Late", yr), "Total", sep = ' ')
KarlukRedoTable[c(7, 10, 13, 16), 14] <- formatC(x = c(abundances.vec, sum(abundances.vec)), digits = 0, format = "f", big.mark = ",")

KarlukRedoTable[14:15, c(4:6, 8:9)] <- formatC(x = stratified.proportions[, c("median", "5%", "95%", "mean", "sd")] * 100, format = "f", digits = 1)
KarlukRedoTable[14:15, 7] <- formatC(x = stratified.proportions[, c("P=0")] * 100, format = "f", digits = 2)
KarlukRedoTable[14:15, 11:15] <- formatC(x = stratified.abundances[, c("median", "5%", "95%", "mean", "sd")], format = "f", digits = 0, big.mark = ",")


require(xlsx)
write.xlsx(x = KarlukRedoTable, file = "Tables/KarlukSmolt2013-2014TablesNEW.xlsx", sheetName = sheetname, col.names = FALSE, row.names = FALSE, append = TRUE)
#~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Exploratory Plots ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

opar=par()
require(gplots)
par(family='serif', mar=c(4.1, 4.1, 1.1, 1.1))
# Plot size is 630 by 630

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create Functions to ease in plot making
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

StrataProportionsBarplot <- function(posterior, yr, dates, bystock = FALSE, legend.x = "topleft", legend.y = NULL) {

  median.table <- sapply(posterior, function(strata) {apply(strata, 2, median)} )
  dimnames(median.table) <- list(KarlukGroups, dates)
  if(bystock) {median.table <- t(median.table)}
  
  ci.l.table <- sapply(posterior, function(strata) {apply(strata, 2, function(RG) {quantile(x = RG, probs = 0.05)} )} )
  dimnames(ci.l.table) <- list(KarlukGroups, dates)
  if(bystock) {ci.l.table <- t(ci.l.table)}
  
  ci.u.table <- sapply(posterior, function(strata) {apply(strata, 2, function(RG) {quantile(x = RG, probs = 0.95)} )} )
  dimnames(ci.u.table) <- list(KarlukGroups, dates)
  if(bystock) {ci.u.table <- t(ci.u.table)}
  
  if(bystock) {
    if(nrow(median.table) == 3) {
      cols <- c("black", "grey", "white")
    } else {
      cols <- c("grey40", "white")
    }
  } else {
    cols <- c("blue", "red")
  }           
  
  barplot2(height = median.table, beside = TRUE, plot.ci = TRUE, ylim = c(0, 1), ylab = "Proportion", xlab = ifelse(bystock, "Stock", "Stratum"), col = cols, cex.axis = 1.2, cex.lab = 1.2, cex.names = 1.2, ci.l = ci.l.table, ci.u = ci.u.table)
  abline(h = 0)
  
  if(bystock) {
    legend.text <- dates
  } else {
    legend.text <- KarlukGroups
  }
  
  legend(x = legend.x, y = legend.y, legend = legend.text, fill = cols, bty = 'n', cex = 1.2, y.intersp = 1.2)
  mtext(text = yr, side = 3, cex = 1.5, line = -0.5)
}


StrataAbundanceBarplot <- function(posterior, yr, dates, bystock = FALSE, legend.x = "topleft", legend.y = NULL) {
  median.table <- sapply(posterior, function(strata) {apply(strata, 2, median)} ) / 100000
  dimnames(median.table) <- list(KarlukGroups, dates)
  if(bystock) {median.table <- t(median.table)}
  
  ci.l.table <- sapply(posterior, function(strata) {apply(strata, 2, function(RG) {quantile(x = RG, probs = 0.05)} )} ) / 100000
  dimnames(ci.l.table) <- list(KarlukGroups, dates)
  if(bystock) {ci.l.table <- t(ci.l.table)}
  
  ci.u.table <- sapply(posterior, function(strata) {apply(strata, 2, function(RG) {quantile(x = RG, probs = 0.95)} )} ) / 100000
  dimnames(ci.u.table) <- list(KarlukGroups, dates)
  if(bystock) {ci.u.table <- t(ci.u.table)}
  
  if(bystock) {
    if(nrow(median.table) == 3) {
      cols <- c("black", "grey", "white")
    } else {
      cols <- c("grey40", "white")
    }
  } else {
    cols <- c("blue", "red")
  }  
  
  barplot2(height = median.table, beside = TRUE, plot.ci = TRUE, ylim = c(0, 2.5), ylab = "Outmigrating Smolt (100K)", xlab = ifelse(bystock, "Stock", "Stratum"), col = cols, cex.axis = 1.2, cex.lab = 1.2, cex.names = 1.2, ci.l = ci.l.table, ci.u = ci.u.table)
  abline(h = 0)
  
  if(bystock) {
    legend.text <- dates
  } else {
    legend.text <- KarlukGroups
  }
  
  legend(x = legend.x, y = legend.y, legend = legend.text, fill = cols, bty = 'n', cex = 1.2, y.intersp = 1.2)
  mtext(text = yr, side = 3, cex = 1.5, line = -0.5)
}


AnnualProportionBarplot <- function(posterior, yr) {
  median.table <- apply(posterior, 2, median)
  names(median.table) <- KarlukGroups
  barplot2(height = median.table, beside = TRUE, plot.ci = TRUE, ylim = c(0, 1), ylab = "Proportion of annual smolt outmigration", xlab = "Stock", col = c("blue", "red"), cex.axis = 1.2, cex.lab = 1.2, cex.names = 1.2,
           ci.l = apply(posterior, 2, function(RG) {quantile(x = RG, probs = 0.05)} ),
           ci.u = apply(posterior, 2, function(RG) {quantile(x = RG, probs = 0.95)} ))
  abline(h = 0)
  mtext(text = yr, side = 3, cex = 1.5, line = -0.5)
}


AnnualAbundancesBarplot <- function(posterior, yr) {
  
  if(dim(posterior)[1] > 2) {
    median.table <- apply(posterior, 2, median) / 100000
    names(median.table) <- KarlukGroups
    ci.l.table <- apply(posterior, 2, function(RG) {quantile(x = RG, probs = 0.05)} ) / 100000
    ci.u.table <- apply(posterior, 2, function(RG) {quantile(x = RG, probs = 0.95)} ) / 100000
  } else {
    colnames(posterior) <- tolower(colnames(posterior))
    rownames(posterior) <- KarlukGroups
    median.table <- posterior[, "median"] / 100000
    ci.l.table <- posterior[, "5%"] / 100000
    ci.u.table <- posterior[, "95%"] / 100000
  }
    
  barplot2(height = median.table, beside = TRUE, plot.ci = TRUE, ylim = c(0, 7), ylab = "Proportion of annual smolt outmigration", xlab = "Stock", col = c("blue", "red"), cex.axis = 1.2, cex.lab = 1.2, cex.names = 1.2,
           ci.l = ci.l.table, ci.u = ci.u.table)
  abline(h = 0)
  mtext(text = yr, side = 3, cex = 1.5, line = -0.5)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2013 ####
### Proportion
## Strata
StrataProportionsBarplot(posterior = KarlukSmoltPosteriors$Output[1:3], yr = 2013, dates = c("May 16-29", "May 30-June 10", "June 10-24"))
StrataProportionsBarplot(posterior = KarlukSmoltPosteriors$Output[1:3], yr = 2013, dates = c("May 16-29", "May 30-June 10", "June 10-24"), bystock = TRUE)

## Annual total
AnnualProportionBarplot(posterior = Karluk2013StratifiedEstimatesWithMRError$Output, yr = 2013)

#~~~~~~~~~~~~~~~~~~
### Outmigrants
## Strata
StrataAbundanceBarplot(posterior = Karluk2013AbundancesWithMRError$Output, yr = 2013, dates = c("May 16-29", "May 30-June 10", "June 10-24"))
StrataAbundanceBarplot(posterior = Karluk2013AbundancesWithMRError$Output, yr = 2013, dates = c("May 16-29", "May 30-June 10", "June 10-24"), bystock = TRUE)

## Annual total
AnnualAbundancesBarplot(posterior = Karluk2013StratifiedAbundancesWithMRError, yr = 2013)

#~~~~~~~~~~~~~~~~~~
#### 2013 Ages  ####
### Proportion
## Strata
StrataProportionsBarplot(posterior = KarlukSmoltPosteriors$Output[11:13], yr = "Age 2 - 2013", dates = c("May 16-29", "May 30-June 10", "June 10-24"))
StrataProportionsBarplot(posterior = KarlukSmoltPosteriors$Output[11:13], yr = "Age 2 - 2013", dates = c("May 16-29", "May 30-June 10", "June 10-24"), bystock = TRUE)

## Annual total
AnnualProportionBarplot(posterior = KarlukSmoltPosteriors$Output$KarlukSmolt2013.Age1, yr = "Age 1 - 2013")
AnnualProportionBarplot(posterior = Karluk2013Age2StratifiedEstimatesWithMRError$Output, yr = "Age 2 - 2013")
AnnualProportionBarplot(posterior = KarlukSmoltPosteriors$Output$KarlukSmolt2013.Age3, yr = "Age 3 - 2013")

#~~~~~~~~~~~~~~~~~~
### Outmigrants
## Strata
StrataAbundanceBarplot(posterior = list(AbundancePosteriorKarlukSmolt2013.Age2.1,
                                        AbundancePosteriorKarlukSmolt2013.Age2.2,
                                        AbundancePosteriorKarlukSmolt2013.Age2.3),
                       yr = "Age 2 - 2013", dates = c("May 16-29", "May 30-June 10", "June 10-24"))
StrataAbundanceBarplot(posterior = list(AbundancePosteriorKarlukSmolt2013.Age2.1,
                                        AbundancePosteriorKarlukSmolt2013.Age2.2,
                                        AbundancePosteriorKarlukSmolt2013.Age2.3),
                       yr = "Age 2 - 2013", dates = c("May 16-29", "May 30-June 10", "June 10-24"), bystock = TRUE)
## Annual total
AnnualAbundancesBarplot(posterior = Karluk2013Age1AbundancesWithMRError, yr = "Age 1 - 2013")
AnnualAbundancesBarplot(posterior = Karluk2013Age2StratifiedAbundancesWithMRError, yr = "Age 2 - 2013")
AnnualAbundancesBarplot(posterior = Karluk2013Age3AbundancesWithMRError, yr = "Age 3 - 2013")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2014 ####
### Proportion
## Strata
StrataProportionsBarplot(posterior = KarlukSmoltPosteriors$Output[4:6], yr = 2014, dates = c("May 13-30", "May 31-June 15", "June 16-July 2"))
StrataProportionsBarplot(posterior = KarlukSmoltPosteriors$Output[4:6], yr = 2014, dates = c("May 13-30", "May 31-June 15", "June 16-July 2"), bystock = TRUE)

## Annual total
AnnualProportionBarplot(posterior = Karluk2014StratifiedEstimatesWithMRError$Output, yr = 2014)

#~~~~~~~~~~~~~~~~~~
### Outmigrants
## Strata
StrataAbundanceBarplot(posterior = Karluk2014AbundancesWithMRError$Output, yr = 2014, dates = c("May 13-30", "May 31-June 15", "June 16-July 2"))
StrataAbundanceBarplot(posterior = Karluk2014AbundancesWithMRError$Output, yr = 2014, dates = c("May 13-30", "May 31-June 15", "June 16-July 2"), bystock = TRUE)

## Annual total
AnnualAbundancesBarplot(posterior = Karluk2014StratifiedAbundancesWithMRError, yr = 2014)

#~~~~~~~~~~~~~~~~~~
#### 2014 Ages  ####
### Proportion
## Strata
StrataProportionsBarplot(posterior = KarlukSmoltPosteriors$Output[15:16], yr = "Age 1 - 2014", dates = c("May 13-June 15", "June 16-July 2"), legend.x = 5, legend.y = 1)
StrataProportionsBarplot(posterior = KarlukSmoltPosteriors$Output[15:16], yr = "Age 1 - 2014", dates = c("May 13-June 15", "June 16-July 2"), bystock = TRUE, legend.x = 4, legend.y = 1)

StrataProportionsBarplot(posterior = KarlukSmoltPosteriors$Output[17:19], yr = "Age 2 - 2014", dates = c("May 13-30", "May 31-June 15", "June 16-July 2"), legend.x = 6, legend.y = 1)
StrataProportionsBarplot(posterior = KarlukSmoltPosteriors$Output[17:19], yr = "Age 2 - 2014", dates = c("May 13-30", "May 31-June 15", "June 16-July 2"), bystock = TRUE)

## Annual total
AnnualProportionBarplot(posterior = Karluk2014Age1StratifiedEstimatesWithMRError$Output, yr = "Age 1 - 2014")
AnnualProportionBarplot(posterior = Karluk2014Age2StratifiedEstimatesWithMRError$Output, yr = "Age 2 - 2014")


#~~~~~~~~~~~~~~~~~~
### Outmigrants
## Strata
StrataAbundanceBarplot(posterior = list(AbundancePosteriorKarlukSmolt2014.Age1.12,
                                        AbundancePosteriorKarlukSmolt2014.Age1.3),
                       yr = "Age 1 - 2014", dates = c("May 13-June 15", "June 16-July 2"))
StrataAbundanceBarplot(posterior = list(AbundancePosteriorKarlukSmolt2014.Age1.12,
                                        AbundancePosteriorKarlukSmolt2014.Age1.3),
                       yr = "Age 1 - 2014", dates = c("May 13-June 15", "June 16-July 2"), bystock = TRUE)

StrataAbundanceBarplot(posterior = list(AbundancePosteriorKarlukSmolt2014.Age2.1,
                                        AbundancePosteriorKarlukSmolt2014.Age2.2,
                                        AbundancePosteriorKarlukSmolt2014.Age2.3),
                       yr = "Age 2 - 2014", dates = c("May 16-29", "May 30-June 10", "June 10-24"))
StrataAbundanceBarplot(posterior = list(AbundancePosteriorKarlukSmolt2014.Age2.1,
                                        AbundancePosteriorKarlukSmolt2014.Age2.2,
                                        AbundancePosteriorKarlukSmolt2014.Age2.3),
                       yr = "Age 2 - 2014", dates = c("May 16-29", "May 30-June 10", "June 10-24"), bystock = TRUE)
## Annual total
AnnualAbundancesBarplot(posterior = Karluk2014Age1StratifiedAbundancesWithMRError, yr = "Age 1 - 2014")
AnnualAbundancesBarplot(posterior = Karluk2014Age2StratifiedAbundancesWithMRError, yr = "Age 2 - 2014")


### Weir smolt
AnnualProportionBarplot(posterior = KarlukSmoltPosteriors$Output$SKARLW14s, yr = "2014 Weir Grab Sample")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Initial Setup for 2015 Data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/Karluk Smolt 2013-2015/Mixtures")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

username <- "krshedd"
password <- "********"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Pull genotypes from LOKI 2015 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create Locus Control
CreateLocusControl.GCL(markersuite = "Sockeye2011_96SNPs", username = username, password = password)

## Save original LocusControl
loci96 <- LocusControl$locusnames
mito.loci96 <- which(LocusControl$ploidy == 1)

dput(x = LocusControl, file = "Objects/OriginalLocusControl96_2015.txt")
dput(x = loci96, file = "Objects/loci96_2015.txt")
dput(x = mito.loci96, file = "Objects/mito.loci96_2015.txt")

#~~~~~~~~~~~~~~~~~~
## Pull all data for each silly code and create .gcl objects for each
LOKI2R.GCL(sillyvec = "SKARLW15s", username = username, password = password)
rm(username, password)
objects(pattern = "\\.gcl")

## Save unaltered .gcl's as back-up:
dir.create("Raw genotypes/OriginalCollections")
invisible(sapply(c("SKARLW15s"), function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections/" , silly, ".txt", sep = ''))} )); beep(8)

## Original sample sizes by SILLY
SKARLW15s.gcl$n

## Fish IDs
SKARLW15s.gcl$attributes$FK_FISH_ID


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/Karluk Smolt 2013-2015/Mixtures")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

## Get objects
LocusControl <- dget(file = "Objects/OriginalLocusControl96_2015.txt")

## Get original mixtures
SKARLW15s.gcl <- dget(file = "Raw genotypes/OriginalCollections/SKARLW15s.txt")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata by time for 2013 and 2014 trap fish ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### 2013

## Stratum 1 -  5/16-29

unique(SKARL13s.gcl$attributes$CAPTURE_DATE)
lapply(list(1:8, 9:16, 17:27), function(lst) {sum(table(SKARL13s.gcl$attributes$CAPTURE_DATE)[lst])})

KarlukSmolt2013.1_IDs <- AttributesToIDs.GCL(silly="SKARL13s",attribute="CAPTURE_DATE",matching=unique(SKARL13s.gcl$attributes$CAPTURE_DATE)[1:8])
KarlukSmolt2013.1_IDs <- list(as.numeric(na.omit(KarlukSmolt2013.1_IDs)))
names(KarlukSmolt2013.1_IDs) <- "SKARL13s"

PoolCollections.GCL("SKARL13s",loci=loci96,IDs=KarlukSmolt2013.1_IDs,newname="KarlukSmolt2013.1")
KarlukSmolt2013.1.gcl$n ##  

## Stratum 2 -  5/30-6/10

KarlukSmolt2013.2_IDs <- AttributesToIDs.GCL(silly="SKARL13s",attribute="CAPTURE_DATE",matching=unique(SKARL13s.gcl$attributes$CAPTURE_DATE)[9:16])
KarlukSmolt2013.2_IDs <- list(as.numeric(na.omit(KarlukSmolt2013.2_IDs)))
names(KarlukSmolt2013.2_IDs) <- "SKARL13s"

PoolCollections.GCL("SKARL13s",loci=loci96,IDs=KarlukSmolt2013.2_IDs,newname="KarlukSmolt2013.2")
KarlukSmolt2013.2.gcl$n ##  

## Stratum 3 -  6/11-6/24

KarlukSmolt2013.3_IDs <- AttributesToIDs.GCL(silly="SKARL13s",attribute="CAPTURE_DATE",matching=unique(SKARL13s.gcl$attributes$CAPTURE_DATE)[17:27])
KarlukSmolt2013.3_IDs <- list(as.numeric(na.omit(KarlukSmolt2013.3_IDs)))
names(KarlukSmolt2013.3_IDs) <- "SKARL13s"

PoolCollections.GCL("SKARL13s",loci=loci96,IDs=KarlukSmolt2013.3_IDs,newname="KarlukSmolt2013.3")
KarlukSmolt2013.3.gcl$n ##  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### QC ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(xlsx)

KarlukMixtures <- unlist(strsplit(ls(pattern='\\.gcl'),split='\\.gcl'))
dput(x=KarlukMixtures,file="Objects/KarlukMixtures.txt")

KarlukMixtures_SampleSizes <- matrix(data=NA,nrow=length(KarlukMixtures),ncol=5,dimnames=list(KarlukMixtures,c("Genotyped","Alternate","Missing","Duplicate","Final")))


#### Check loci
## Get sample size by locus
Original_KarlukMixtures_SampleSizebyLocus <- SampSizeByLocus.GCL(KarlukMixtures, loci96)
min(Original_KarlukMixtures_SampleSizebyLocus) ## 70
apply(Original_KarlukMixtures_SampleSizebyLocus,1,min)/apply(Original_KarlukMixtures_SampleSizebyLocus,1,max) # Known from QC that SKARLW14s had issues

t(sort(Original_KarlukMixtures_SampleSizebyLocus[KarlukMixtures[9],]/max(Original_KarlukMixtures_SampleSizebyLocus[KarlukMixtures[9],]),decreasing=TRUE)) # One_ZNF.61 is just below the 80% rule for SKARLW14s

#### Check individuals
### Initial
## Get number of individuals per silly before removing missing loci individuals
Original_KarlukMixtures_ColSize <- sapply(paste(KarlukMixtures,".gcl",sep=''), function(x) get(x)$n)
KarlukMixtures_SampleSizes[,"Genotyped"] <- Original_KarlukMixtures_ColSize


### Alternate
## Indentify alternate species individuals
KarlukMixtures_Alternate <- FindAlternateSpecies.GCL(sillyvec=KarlukMixtures, species="sockeye")

## Remove Alternate species individuals
RemoveAlternateSpecies.GCL(AlternateSpeciesReport=KarlukMixtures_Alternate, AlternateCutOff=0.5, FailedCutOff=0.5)

## Get number of individuals per silly after removing alternate species individuals
ColSize_KarlukMixtures_PostAlternate <- sapply(paste(KarlukMixtures,".gcl",sep=''), function(x) get(x)$n)
KarlukMixtures_SampleSizes[,"Alternate"] <- Original_KarlukMixtures_ColSize-ColSize_KarlukMixtures_PostAlternate


### Missing
## Remove individuals with >20% missing data
KarlukMixtures_MissLoci <- RemoveIndMissLoci.GCL(sillyvec=KarlukMixtures,loci=loci96,proportion=0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_KarlukMixtures_PostMissLoci <- sapply(paste(KarlukMixtures,".gcl",sep=''), function(x) get(x)$n)
KarlukMixtures_SampleSizes[,"Missing"] <- ColSize_KarlukMixtures_PostAlternate-ColSize_KarlukMixtures_PostMissLoci


### Duplicate
## Check within collections for duplicate individuals.
KarlukMixtures_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec=KarlukMixtures,loci=loci96,quantile=NULL,minproportion=0.95)
KarlukMixtures_DuplicateCheckReportSummary <- sapply(KarlukMixtures, function(x) KarlukMixtures_DuplicateCheck95MinProportion[[x]]$report)

## Remove duplicate individuals
KarlukMixtures_RemovedDups <- RemoveDups.GCL(KarlukMixtures_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_KarlukMixtures_PostDuplicate <- sapply(paste(KarlukMixtures,".gcl",sep=''), function(x) get(x)$n)
KarlukMixtures_SampleSizes[,"Duplicate"] <- ColSize_KarlukMixtures_PostMissLoci-ColSize_KarlukMixtures_PostDuplicate


### Final
KarlukMixtures_SampleSizes[,"Final"] <- ColSize_KarlukMixtures_PostDuplicate
KarlukMixtures_SampleSizes

write.xlsx(KarlukMixtures_SampleSizes,file="Output/KarlukMixtures_SampleSizes.xlsx")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Combine Loci ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get the final locus set from baseline file to avoid any errors
loci91 <- dget(file="V:/WORK/Sockeye/Kodiak/2014 Karluk Baseline/Objects/loci91.txt")

## Combine loci
CombineLoci.GCL(sillyvec=KarlukMixtures,markerset=c("One_Cytb_26","One_CO1","One_Cytb_17"),update=T)

## NOTE THAT THESE LOCI ARE DROPPED c("One_CO1","One_Cytb_17","One_Cytb_26","One_MHC2_251","One_GPDH","One_Tf_ex3-182")
