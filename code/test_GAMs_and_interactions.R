#################################################
# Malin Pinsky climate velocity analysis
# Test of GAM, including interactions with traits, as well as random forest
# Martin Lindegren, DTU, 20200325
#################################################

rm(list=ls())
setwd("H:/Work/Collaboration ms/Malin climate velocities/Analysis")
datMerge<-read.csv("datMerge.csv",sep = ",",header=T)
head(datMerge)

library(mgcv)
library(randomForest)
library(ggplot2)

# Some exploratory plots
ggplot(aes(y = obslat1, x = factor(DemersPelag3)), data = datMerge)+ geom_boxplot()+theme(axis.text.x = element_text(size=10, angle=45))                                                                                                                                            
ggplot(aes(y = obslat1, x = factor(fished3)), data = datMerge) + geom_boxplot()+theme(axis.text.x = element_text(size=10, angle=45))#+coord_cartesian(ylim = c(0, 1000))                                                                                                                                              
# Assign continuous traits to categories and make box plots by groups
funPlot<-function(B,N){ # B = quantiles, N= variable column)
  datMerge$factorT <- with(datMerge, cut(datMerge[,N],breaks=quantile(datMerge[,N], probs=seq(0,1, by=B), na.rm=TRUE)))
  ggplot(aes(y = obslat1, x = factor(factorT)), data = datMerge) + geom_boxplot()+theme(axis.text.x = element_text(size=10, angle=45))#                                                                                                                                           
}
funPlot(B=0.2,14) # Test 5 groups for MaxLength

# Test distribution of response variables
hist(datMerge$obslat1,breaks=20) # Including direction - normally distributed
hist(datMerge$obslat1abs,breaks=20) # Absolute values - strong left-skew (due to using absolute values) - consdier logging or use appropriate family distribution/link function

# Test GAMs with obslat1 as response and only climate velocity as predictor
fit1<-gam(obslat1~s(gamhlat1,k=3),family = gaussian(link = "identity"),data=datMerge)#
summary(fit1) #
plot(fit1,res=T,rug=F,pch=20,shade=T,shade.col="grey")
hist(residuals(fit1),breaks=20)  
qqnorm(residuals(fit1)) # Hmm - tails and some strong outliers..check when more covariates are added and if remains consider removing the observation

# Explore if the relationship with climate velocity differs between regions and families
fit2<-gam(obslat1~s(gamhlat1,as.factor(region),bs="fs",k=3),family = gaussian(link = "identity"),data=datMerge)#
summary(fit2) #
plot(fit2,res=T,rug=F,pch=20,shade=T,shade.col="grey") # Quite some differences between regions 
# interesting to explore in more detail - do areas with different absolute temperature or range show different responses (i,e., slopes)

fit3<-gam(obslat1~s(gamhlat1,as.factor(family),bs="fs",k=3),family = gaussian(link = "identity"),data=datMerge)#
summary(fit3) #
plot(fit3,res=T,rug=F,pch=20,shade=T,shade.col="grey") # Highly variable responses to climate velocities between families 
# interesting to explore in more detail, particularly if families with particular traits show different responses (i.e., slopes)
# As an additional test we could for instance extract fitted intercepts and slopes and test against traits?

# Test adding family as random factor while adding tratis as fixed effects: use simple formulation of random effects (bs="re"): https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/random.effects.html
# Explore if the relationship with climate velocity differs between demersal/pelagic and fishes/unfished
fit4<-gam(obslat1~s(gamhlat1,as.factor(DemersPelag3),bs="fs",k=3)+s(family,bs="re"),family = gaussian(link = "identity"),data=datMerge)#
summary(fit4) #
plot(fit4,sel=1,res=T,rug=F,pch=20,shade=T,shade.col="grey") # Hmm, slightly different responses between groups - worth looking into

fit5<-gam(obslat1~s(gamhlat1,as.factor(fished3),bs="fs",k=3)+s(family,bs="re"),family = gaussian(link = "identity"),data=datMerge)#
summary(fit5) #
plot(fit5,sel=1, res=T,rug=F,pch=20,shade=T,shade.col="grey") # No, rather similar responses in both groups - hence fishing does not seem to be a strong factor

# Test if the relationship with climate velocity differs between species with different trait combinations
# Create a categorical (factor) of slow_fast instead of continous (here 4 levels) 
datMerge$slowfastFac <- with(datMerge, cut(slow_fast,breaks=quantile(slow_fast, probs=seq(0,1, by=0.25), na.rm=TRUE)))#
table(datMerge$slowfastFac)

fit6<-gam(obslat1~s(gamhlat1,slowfastFac,bs="fs",k=3)+s(family,bs="re"),family = gaussian(link = "identity"),data=datMerge)#
summary(fit6) #
plot(fit6,sel=1,res=T,rug=F,pch=20,shade=T,shade.col="grey",ylim=c(-0.10,0.15)) 
# Hmm - potentially really interesting. Some differences in the responses to climate velocity between "slow" and "fast" species. Worth looking into in more detail

# Let?s explore this potential trait dependent responses for other continuous traits -re-classified to categorical factors
funTest<-function(B,N){ # B = number of breaks (i.e., categories, N= variable column)
  datMerge$factorT <- with(datMerge, cut(datMerge[,N],breaks=quantile(datMerge[,N], probs=seq(0,1, by=B), na.rm=TRUE)))#,labels=c("a","b","c","d")))# 
  fit<-gam(obslat1~s(gamhlat1,as.factor(factorT),bs="fs",k=3)+s(family,bs="re"),family = gaussian(link = "identity"),data=datMerge)#+s(biomass.b.lm,k=3)+s(tsdur,k=3)
  plot(fit,sel=1,res=T,rug=T)
  return(summary(fit)) #
}
#names(datMerge)
#[1] "sciname"       "region"        "family"        "genus"         "species"       "obslat1"       "gamhlat1"     
#[8] "obsdepth1"     "gamhdepth1"    "fishinvert"    "DemersPelag3"  "biomass.b.lm"  "fished3"       "MaxLength"    
#[15] "Troph"         "K"             "range"         "tsdur"         "extentlat"     "slow_fast"     "Equilibrium"  
#[22] "Periodic"      "Opportunistic" "obslat1abs"    "gamhlat1abs"   "slowfastFac"
funTest(B=0.25,N=14) #14 is max length - seem to show rather similar responses between size groups
funTest(B=0.25,N=15) #15 is Trophic level - similar responses
funTest(B=0.25,N=16) #16 is K - rather different responses! Potentially rather interesting - Could indicate that slow growing species (black) respond more rapidly to small changes in climate velocities and vice versa for fast growers.
funTest(B=0.5,N=16) # K: Even more pronounced if we just cut K into two groups (slow - fast)
funTest(B=0.25,N=21) #21 is Equilibrium - a bit more varied responses (i.e., slopes)
funTest(B=0.25,N=22) #22 is Periodic - rather similar responses
funTest(B=0.25,N=23) #23 is Opportunistic - rather similar responses

datMerge$Kcategory <- with(datMerge, cut(K,breaks=quantile(K, probs=seq(0,1, by=0.5), na.rm=TRUE)))#
table(datMerge$Kcategory)

# Test adding more covariates besides K - range tsdur etc etc - but should be using a formal model selection instead!
fit7<-gam(obslat1~s(gamhlat1,Kcategory,bs="fs",k=3)+s(range,k=3)+s(tsdur,k=3)+DemersPelag3+fished3+s(family,bs="re"),data=datMerge,family = gaussian(link = "identity"))#
summary(fit7) # 
plot(fit7,sel=1,res=T,rug=F,pch=20,shade=T,shade.col="grey",ylim=c(-0.06,0.1)) 

# Test some diagnostics
gam.check(fit7)
hist(residuals(fit7),breaks=20)  
qqnorm(residuals(fit7)) # Hmm - ok BUT two strong outliers..check and remove

which(residuals(fit7)>0.15) # Observation 10


############################################
# Test random forest
datMerge2<-datMerge[complete.cases(datMerge), ] #RF does not handle NAs
fitR <- randomForest(obslat1 ~ gamhlat1+DemersPelag3+fished3+MaxLength+Troph+K+slow_fast+range+tsdur,  data=datMerge2, ntree=1000, mtry=3)
print(fitR) # Well - in terms of r2 it is much worse than the above mentioned GAMs where explicit interactions were added
importance(fitR, tytype = 2) # Relative importance scores
varImpPlot(fitR) # As expected climate velocity has the highest importance

plot(fitR)
par(mfrow=c(4,2),mar=c(4,4,0.5,0.5))
partialPlot(fitR, datMerge2,gamhlat1,rug=F,main="")
partialPlot(fitR, datMerge2,MaxLength,rug=F,main="")
partialPlot(fitR, datMerge2,Troph,rug=F,main="")
partialPlot(fitR, datMerge2,K,rug=F,main="")
partialPlot(fitR, datMerge2,slow_fast,rug=F,main="")
partialPlot(fitR, datMerge2,range,rug=F,main="")
partialPlot(fitR, datMerge2,tsdur,rug=F,main="")
par(mfrow=c(1,1),mar=c(5,4,4,2))


