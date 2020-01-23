###########################################################################################
## Multiple Linear Regression of traits vs. observed shifts in species geographic ranges ##
## Only fish
## Incorporating new trait data from Beukhof et al. 2019 
###########################################################################################
require(visreg) # for partial residuals
require(AICcmodavg) # for AICcs
require(dplyr)
source('code/lm.aic.all.r')
directionlh = read.csv('data/directionlh.csv', row.names=1) # data from Pinsky et al. 2013 (shifts and some traits)
traits = read.csv('data/Beaukhof_traits.csv',header = T, sep = ";", row.names = NULL) # data from Beukhof et al. 2019 (new traits)

###############################################
# Add more trait data for fishes
###############################################
# First remove invertebrates
directionlh<-directionlh[-c(which(directionlh$fishinvert==FALSE)),]
dim(directionlh)

# Add LME and FAO columns and recode survey regions into LME/FAO numbers (as in the trait data file)
directionlh$LME<-directionlh$region
unique(directionlh$LME)
directionlh$LME2 <- as.factor(recode(directionlh$LME,`AFSC_Aleutians`=65,`AFSC_EBS`=1,`AFSC_GOA`=2,`DFO_Newfoundland_Fall`=9,
                                  `DFO_ScotianShelf`=8,`DFO_SoGulf`=8,`NEFSC_Spring`=7,`SEFSC_GOMex`=5,`WestCoast_Tri`=3))
unique(directionlh$LME2) # looks good (no NAs)
#directionlh[is.na(directionlh$LME2), c('LME', 'LME2')] # check NAs if they exist

# Add FAO areas
directionlh$FAO <- directionlh$LME2
directionlh$FAO <- recode(directionlh$FAO,`65`=67,`1`=61,`2`=67,`9`=21,`8`=21,`7`=21,`5`=31,`3`=77)
unique(directionlh$FAO) #

# Fix some species names in directionlh to match names in traits
directionlh$species <- as.character(directionlh$species) # first change to character so that we can edit
directionlh$genus <- as.character(directionlh$genus)
inds <- directionlh$genus == 'Theragra' & directionlh$species == 'chalcogramma'
directionlh$genus[inds] <- 'Gadus'
directionlh$species[inds] <- 'chalcogrammus'
inds <- directionlh$genus == 'Clupea' & directionlh$species == 'pallasii'
directionlh$species[inds] <- 'pallasii pallasii'


# Create unique identifier FAO /LME /Genus /species
directionlh$SpecUQ <- paste(directionlh$FAO, directionlh$LME2, directionlh$genus, directionlh$species)

# Replace empty cells with NA
traits <- traits %>% mutate_all(na_if,"")
# Create unique identifier same as before
traits$SpecUQ <- paste(traits$FAO, traits$LME, traits$genus, traits$species)

# Check for matching records
MatchRec <- match(unique(sort(directionlh$SpecUQ)), unique(sort(traits$SpecUQ)))
unique(sort(directionlh$SpecUQ))[is.na(MatchRec)] # 
# 18 fish species are not matching. This seem to be due to
# trait records for the particular LME that are for some reason missing (e.g., Icelus spatula and Sebastes mentella in LME 9). 
# These need to be checked and filled out manually or by assigning trait values from neighbouring LMEs.

# Once the above issue is fully fixed - merge the data files by unique id (here excluding columns with references)
datMerge<-merge(directionlh, traits[,c(42,9,11,13,16,17,18,21,24,26,29,32,33,36,39)], by="SpecUQ")#,all.x=T)
dim(datMerge) # 241 matching species
unique(datMerge$SpecUQ) # but only 35 unique species entries 
which(duplicated(datMerge$SpecUQ)==TRUE) # Hmm - unsure why these are duplicated - check
datMerge[datMerge$SpecUQ=='21 8 Amblyraja radiata',] # an example. They are duplicated because there are two surveys (and therefore two measurements of a population shift) in this LME
#datMerge<-datMerge[-which(duplicated(datMerge$SpecUQ)==TRUE),] # keep the duplicates, so don't run this line

# Compare some old vs new trait values that exist in both data sets - seem ok
plot(datMerge[,15], datMerge[,36], ylab="Lmax Esther", xlab="Lmax old"); abline(0, 1)
plot(datMerge[,17],datMerge[,35],ylab="K Esther",xlab="K old"); abline(0, 1)
plot(datMerge[,16],datMerge[,26],ylab="TL Esther",xlab="TL old"); abline(0, 1)
names(datMerge)
summary(datMerge[,c(1,24:37)]) # Note that some trait values are still missing, especially aspect ratio (AR) that cannot really be calculated for some species. Note that Aurore mentioned that she has updated some missing values)



#############################################
# Statistical models
###########################################

# Add log variables and absolute rate of shift
# use new variables from 
datMerge$loggrowth.coefficient = log(datMerge$growth.coefficient)
datMerge$loglength.max = log(datMerge$length.max)
datMerge$loglength.infinity = log(datMerge$length.infinity)
datMerge$logage.max = log(datMerge$age.max)
datMerge$logoffspring.size = log(datMerge$offspring.size)
datMerge$logfecundity = log(datMerge$fecundity)

datMerge$obslat1abs = abs(datMerge$obslat1)
datMerge$gamhlat1abs = datMerge$gamhlat1*sign(datMerge$obslat1) # predicted shift rate, signed by whether in the direction of observed shift or not

# Sample size by region
table(datMerge$region)

# Sample size by factor (ordered by decreasing size)
sum(!is.na(datMerge$obslat1)) # 241 (up from 199 fishes in Pinsky et al. 2013)
with(datMerge, sum(!is.na(obslat1) & !is.na(gamhlat1abs))) # 241 climate velocity
with(datMerge, sum(!is.na(obslat1) & !is.na(habitat))) # 241 from Beukhof
with(datMerge, sum(!is.na(obslat1) & !is.na(biomass.b.lm))) # 241 biomass trend from survey data
with(datMerge, sum(!is.na(obslat1) & !is.na(fished3))) # 241 commercial fishing or not
with(datMerge, sum(!is.na(obslat1) & !is.na(loglength.max))) # , 241 # max length from Beukhof
with(datMerge, sum(!is.na(obslat1) & !is.na(tl))) # , 241 trophic level from Beukhof
with(datMerge, sum(!is.na(obslat1) & !is.na(body.shape))) # 241 from Beukhof
with(datMerge, sum(!is.na(obslat1) & !is.na(fin.shape))) # 241 from Beukhof
with(datMerge, sum(!is.na(obslat1) & !is.na(logoffspring.size))) # 239 from Beukhof
with(datMerge, sum(!is.na(obslat1) & !is.na(spawning.type))) # 237 from Beukhof
with(datMerge, sum(!is.na(obslat1) & !is.na(logfecundity))) # 237 from Beukhof
with(datMerge, sum(!is.na(obslat1) & !is.na(age.maturity))) # 234 from Beukhof
with(datMerge, sum(!is.na(obslat1) & !is.na(loggrowth.coefficient))) # 232 VB K from Beukhof
with(datMerge, sum(!is.na(obslat1) & !is.na(logage.max))) # 232 from Beukhof
with(datMerge, sum(!is.na(obslat1) & !is.na(loglength.infinity))) # 232 from Beukhof
with(datMerge, sum(!is.na(obslat1) & !is.na(AR))) # 214 aspect ratio? from Beukhof
with(datMerge, sum(!is.na(obslat1) & !is.na(range))) # 196 latitudinal range size from fishbase
# tsdur: time series duration
# extentlat: survey region latitudinal extent
# 

# Choose which variables to include
# Pick ONE of these options and then continue on with "Trim to complete cases"
vars = c('gamhlat1abs', 'fished3', 'habitat', 'biomass.b.lm', 'loglength.max', 'tl', 'body.shape', 'fin.shape', 
         'tsdur', 'extentlat'); yvar = 'logabs' # only complete vars for log(obslat1abs)
vars = c('gamhlat1abs', 'fished3', 'habitat', 'biomass.b.lm', 'loglength.max', 'tl', 'loggrowth.coefficient', 'body.shape', 'fin.shape', 
         'logoffspring.size', 'spawning.type', 'logfecundity', 'age.maturity', 'logage.max', 'loglength.infinity', 'AR', 'range', 'tsdur', 
         'extentlat'); yvar = 'logabs' # all vars for log(obslat1abs). This is very slow and memory intensive to fit!

vars = c('gamhlat1', 'fished3', 'habitat', 'biomass.b.lm', 'loglength.max', 'tl', 'body.shape', 'fin.shape', 
         'tsdur', 'extentlat'); yvar = 'raw' # only complete vars for obslat1
vars = c('gamhlat1', 'fished3', 'habitat', 'biomass.b.lm', 'loglength.max', 'tl', 'loggrowth.coefficient', 'body.shape', 'fin.shape', 
         'logoffspring.size', 'spawning.type', 'logfecundity', 'age.maturity', 'logage.max', 'loglength.infinity', 'AR', 'range', 'tsdur', 
         'extentlat'); yvar = 'logabs' # all vars for obslat1. This is probably slow to fit!


# Trim to complete cases
i = complete.cases(datMerge[,vars[!grepl(':', vars)]]); sum(i)


# Colinearity?
r2 = matrix(NA, nrow=length(vars), ncol=length(vars))
colnames(r2) = vars; rownames(r2) = vars
for(k in 1:(length(vars)-1)){
  for(j in (k+1):length(vars)){
    r2[k,j] = round(cor(as.numeric(datMerge[[vars[k]]][i]), as.numeric(datMerge[[vars[j]]][i]))^2,2)
  }
}
r2	

# All subsets lm		
if(yvar=='logabs') mod = lm.aic.all(y=log(datMerge$obslat1abs[i]), x=datMerge[i,vars], models='all') # test all subsets with log(abs shift)
if(yvar=='raw') mod = lm.aic.all(y=datMerge$obslat1[i], x=datMerge[i,vars], models='all') # test all subsets with raw shift value
aic = numeric(length(mod))
for(i in 1:length(mod)) aic[i] = AIC(mod[[i]])
aicc = numeric(length(mod)) # small-sample version
for(i in 1:length(mod)) aicc[i] = AICc(mod[[i]])
r2 = numeric(length(mod))
for(i in 1:length(mod)) r2[i] = summary(mod[[i]])$r.squared

cor.test(aicc, aic)

modord = order(aic) # order from best to worst

# AIC model weights
aicd = aic - min(aic) # delta aic
aicw = exp(-aicd/2)/sum(exp(-aicd/2), na.rm=TRUE) # AIC weights: Burnham & Anderson p. 75 section 2.9.1

# term weights
for(i in 1:length(vars)){
  cat(paste(format(vars[i], width=17), signif(sum(aicw[grep(vars[i], names(mod))], na.rm=TRUE), 3), '\n'))
}

# full model
summary(mod[[1]])

# best model
summary(mod[[modord[1]]])
aic[modord[1]]
aicw[modord[1]]
#par(mfrow=c(1,3)); visreg(mod[[modord[1]]], trans=exp, type='conditional')

# plot survey duration and climate velocity
#plot(datMerge$tsdur, datMerge$gamhlat1abs)
#summary(lm(datMerge$gamhlat1abs ~ datMerge$tsdur))

# plot survey duration and lags
#plot(datMerge$tsdur, datMerge$obslat1abs - datMerge$gamhlat1abs)
#summary(lm(I(datMerge$obslat1abs - datMerge$gamhlat1abs) ~ datMerge$tsdur))
#i = = which(datMerge$obslat1abs - datMerge$gamhlat1abs < 0.2) # remove an outlier
#summary(lm(I(datMerge$obslat1abs[i] - datMerge$gamhlat1abs[i]) ~ datMerge$tsdur[i]))

# Multi-model inference
options(warn=2) # turn warnings to errors
vars2 = c('Intercept', vars)
param = numeric(length(vars2))
se = numeric(length(vars2))
for(i in 1:length(vars2)){
  print(vars2[i])
  subs = grep(vars2[i], names(mod))
  if(vars2[i] == 'Intercept') subs = 1:length(mod)
  for(j in subs){
    if(all(!is.na(mod[[j]]))){ # skip any models that couldn't be fit
      c = summary(mod[[j]])$coefficients
      k = grep(vars2[i], rownames(c))
      if(vars2[i] == 'Intercept') k=1
      if(length(k)>1) stop('matched >1 row')
      param[i] = param[i] + aicw[j] * c[k,1] # Eq. 4.1 in Burnham & Anderson (p. 150)
      se[i] = se[i] + aicw[j] * sqrt(c[k,2]^2 + (c[k,1]- param[i])^2) # Eq. 4.9 in Burnham & Anderson (p. 162)
    }
  }
}
params = data.frame(var = vars2, param=param, se=se)
print(params)



##########################################################################
# write out results for top 6 models and
# only climate velocity, only survey duration, and only survey extent
##########################################################################
if(yvar=='logabs') outvars = c('family', 'gamhlat1abs', 'maxlatshift', 'extentlat', 'tsdur', 'fishinvert', 'logK', 'fished3', 'DemersPelag3', 'range', 'biomass.b.lm', 'Troph', 'logMaxLength',  'logLongevity') # so I can set the order in the table
if(yvar=='raw') outvars = c('family', 'gamhlat1', 'maxlatshift', 'extentlat', 'tsdur', 'fishinvert', 'logK', 'fished3', 'DemersPelag3', 'range', 'biomass.b.lm', 'Troph', 'logMaxLength',  'logLongevity') # so I can set the order in the table
outvars = outvars[outvars %in% vars]
outvarslong = outvars # for matching to coefficient names
outvarslong[outvarslong=='DemersPelag3'] = 'DemersPelag3pelagic'
outvarslong[outvarslong=='fished3'] = 'fished3no'
outvarslong[outvarslong=='fishinvert'] = 'fishinvertTRUE' # is it a fish?
n = length(outvars)+3
maxmods = 6
out = data.frame(Taxa = tax, Factor = c(outvarslong, 'delAIC', 'r2', 'w'), Weights = character(n), M1 = character(n), M2 = character(n), M3 = character(n), M4 = character(n), M5 = character(n), M6 = character(n), Full = character(n), Ave = character(n), CV = character(n), SD = character(n), SE = character(n), SDSE = character(n), SDSECV = character(n), SppChar = character(n), stringsAsFactors=FALSE)
if(maxmods >6){
  for(i in 7:maxmods){
    nm = paste('M', i, sep='')
    out[[nm]] = character(n)
  }
}
modord = order(aic)
for(i in 1:min(maxmods,length(mod))){
  v = names(coef(mod[[modord[i]]])) # variable names in this model
  s2<- signif(coef(mod[[modord[i]]]),2)
  nm = paste('M', i, sep='')
  for(j in 1:length(outvars)){
    if(any(grepl(outvarslong[j], v))){
      out[[nm]][j] = s2[grep(outvarslong[j], v)]
    }
  }
  out[[nm]][j+1] = signif(aic[modord[i]] - min(aic),3)
  out[[nm]][j+2] = signif(r2[modord[i]],2)
  out[[nm]][j+3] = signif(aicw[modord[i]],2)
}
# Full model
v = names(coef(mod[[1]])) # variable names in this model
s2<- signif(coef(mod[[1]]),2)
for(j in 1:length(outvars)){
  out$Full[j] = s2[grep(outvarslong[j], v)]
}
out$Full[j+1] = signif(aic[1] - min(aic),3)
out$Full[j+2] = signif(r2[1],2)	
out$Full[j+3] = signif(aicw[1],2)
# Averaged model
v = names(coef(mod[[1]])) # variable names in this model (same as full model)
for(j in 1:length(outvars)){
  out$Ave[j] = signif(params$param[grep(outvars[j], params$var)],2)
}
# special models to include as well
if(yvar=='logabs') specialmods = c('y~gamhlat1abs', 'y~tsdur', 'y~extentlat', 'y~tsdur+extentlat', 'y~gamhlat1abs+tsdur+extentlat')
if(yvar=='raw') specialmods = c('y~gamhlat1', 'y~tsdur', 'y~extentlat', 'y~tsdur+extentlat', 'y~gamhlat1+tsdur+extentlat')
if(tax=='All') specialmods = c(specialmods, 'y~fishinvert+fished3+DemersPelag3+biomass.b.lm')
if(tax=='Fish') specialmods = c(specialmods, 'y~fished3+DemersPelag3+biomass.b.lm+logMaxLength+Troph+logK+range')
specialnames = c('CV', 'SD', 'SE', 'SDSE', 'SDSECV', 'SppChar')
for(j in 1:length(specialmods)){
  i = which(names(mod) == specialmods[j])
  v = names(coef(mod[[i]])) # variable names in this model
  s2<- signif(coef(mod[[i]]),2)
  nm = specialnames[j]
  for(j in 1:length(outvars)){
    if(any(grepl(outvarslong[j], v))){
      out[[nm]][j] = s2[grep(outvarslong[j], v)]
    }
  }
  out[[nm]][j+1] = signif(aic[i] - min(aic),3)
  out[[nm]][j+2] = signif(r2[i],2)
  out[[nm]][j+3] = signif(aicw[i],2)
}
# Single-variable models
for(j in 1:length(outvars)){ # add columns for single-variable models
  nm = outvars[j]
  out[[nm]] = character(n) # add the column
  i = which(names(mod) == paste('y~', outvars[j], sep='')) # find the right model
  v = names(coef(mod[[i]])) # variable names in this model
  s2<- signif(coef(mod[[i]]),2)
  out[[nm]][j] = s2[2]
  out[[nm]][length(outvars)+1] = signif(aic[i] - min(aic),3)
  out[[nm]][length(outvars)+2] = signif(r2[i],2)
  out[[nm]][length(outvars)+3] = signif(aicw[i],2)
}
# Term weights
for(i in 1:length(outvars)){
  out$Weights[i] = signif(sum(aicw[grep(outvars[i], names(mod))], na.rm=TRUE), 3)
}

out			

# Write out table of models
write.csv(out, file=paste('output/traitmods_', tax, '_', yvar, '.csv', sep=''), row.names=FALSE)


