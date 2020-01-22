################################
## Multiple Linear Regression ##
################################
setwd('/Users/mpinsky/Documents/Princeton/Trawl Data/')
require(visreg) # for partial residuals
require(AICcmodavg) # for AICcs
source('lm.aic.all 2012-09-07.r')
directionlh = read.csv('Output/directionlh_allyrsonlytemp_2013-07-03.csv', row.names=1)

# Add log variables and absolute rate of shift
	directionlh$logK = log(directionlh$K)
	directionlh$logMaxLength = log(directionlh$MaxLength)
	directionlh$logLongevity = log(directionlh$LongevityWild)

	directionlh$obslat1abs = abs(directionlh$obslat1)
	directionlh$gamhlat1abs = directionlh$gamhlat1*sign(directionlh$obslat1) # predicted shift rate, 

# Sample size by region
	table(directionlh$region)

# Sample size by factor for all taxa
	nrow(directionlh) # 325
	with(directionlh, sum(!is.na(obslat2) & !is.na(DemersPelag3))) # , 325
	with(directionlh, sum(!is.na(obslat2) & !is.na(fished))) # , 325

	with(directionlh, sum(!is.na(obslat2) & !is.na(Troph))) # , 255
	with(directionlh, sum(!is.na(obslat2) & !is.na(fecundity))) # , 180

# Sample size by factor for fish
	sum(directionlh$fishinvert) # , 259 (w/out family, w/out genus)
	with(directionlh, sum(!is.na(obslat2) & fishinvert & !is.na(DemersPelag3))) # , 259
	with(directionlh, sum(!is.na(obslat2) & fishinvert  & !is.na(fished))) # , 259
	with(directionlh, sum(!is.na(obslat2) & fishinvert  & !is.na(MaxLength))) # , 241
	with(directionlh, sum(!is.na(obslat2) & fishinvert  & !is.na(Troph))) # , 241
	with(directionlh, sum(!is.na(obslat2) & fishinvert  & !is.na(K))) # , 231
	with(directionlh, sum(!is.na(obslat2) & fishinvert  & !is.na(range))) # 203

# Choose which variables to include
	vars = c('gamhlat1abs', 'fishinvert', 'fished3', 'DemersPelag3', 'biomass.b.lm', 'tsdur', 'extentlat'); tax='All'; yvar = 'logabs' # for inverts + fish log(obslat1abs) 
	
	vars = c('gamhlat1abs', 'fished3', 'DemersPelag3', 'biomass.b.lm', 'logMaxLength', 'Troph', 'logK', 'tsdur', 'extentlat', 'range'); tax='Fish'; yvar = 'logabs' # for fish log(obslat1abs)

	vars = c('gamhlat1', 'fishinvert', 'fished3', 'DemersPelag3', 'biomass.b.lm', 'tsdur', 'extentlat'); tax='All'; yvar = 'raw' # for inverts + fish obslat1

	vars = c('gamhlat1', 'fished3', 'DemersPelag3', 'biomass.b.lm', 'logMaxLength', 'Troph', 'logK', 'tsdur', 'extentlat', 'range'); tax='Fish'; yvar = 'raw' # for fish obslat1 (no longevity)

# Trim to complete cases
	if(tax=='All'){ i = complete.cases(directionlh[,vars[!grepl(':', vars)]]); tax = 'All'; sum(i)} # with all data
	if(tax=='Fish'){	i = complete.cases(directionlh[,vars[!grepl(':', vars)]]) & directionlh$fishinvert; tax = 'Fish'; sum(i)} # with all data among fishes


# Colinearity?
	r2 = matrix(NA, nrow=length(vars), ncol=length(vars))
		colnames(r2) = vars; rownames(r2) = vars
	for(k in 1:(length(vars)-1)){
		for(j in (k+1):length(vars)){
			r2[k,j] = round(cor(as.numeric(directionlh[[vars[k]]][i]), as.numeric(directionlh[[vars[j]]][i]))^2,2)
		}
	}
	r2	

# All subsets lm		
	if(yvar=='logabs') mod2 = lm.aic.all(y=log(directionlh$obslat1abs[i]), x=directionlh[i,vars], models='all') # test all subsets with log(abs shift)
	if(yvar=='raw') mod2 = lm.aic.all(y=directionlh$obslat1[i], x=directionlh[i,vars], models='all') # test all subsets with raw shift value
		aic = numeric(length(mod2))
		for(i in 1:length(mod2)) aic[i] = AIC(mod2[[i]])
		aicc = numeric(length(mod2)) # small-sample version
		for(i in 1:length(mod2)) aicc[i] = AICc(mod2[[i]])
		r2 = numeric(length(mod2))
		for(i in 1:length(mod2)) r2[i] = summary(mod2[[i]])$r.squared

		cor.test(aicc, aic)

		modord = order(aic) # order from best to worst
		
	# AIC model weights
	aicd = aic - min(aic) # delta aic
	aicw = exp(-aicd/2)/sum(exp(-aicd/2), na.rm=TRUE) # AIC weights: Burnham & Anderson p. 75 section 2.9.1

	# term weights
	for(i in 1:length(vars)){
		cat(paste(format(vars[i], width=17), signif(sum(aicw[grep(vars[i], names(mod2))], na.rm=TRUE), 3), '\n'))
	}

	# full model
	summary(mod2[[1]])
	
	# best model
	summary(mod2[[modord[1]]])
	aic[modord[1]]
	aicw[modord[1]]
	#par(mfrow=c(1,3)); visreg(mod2[[modord[1]]], trans=exp, type='conditional')

	# plot survey duration and climate velocity
		#plot(directionlh$tsdur, directionlh$gamhlat1abs)
		#summary(lm(directionlh$gamhlat1abs ~ directionlh$tsdur))

	# plot survey duration and lags
		#plot(directionlh$tsdur, directionlh$obslat1abs - directionlh$gamhlat1abs)
		#summary(lm(I(directionlh$obslat1abs - directionlh$gamhlat1abs) ~ directionlh$tsdur))
		#i = = which(directionlh$obslat1abs - directionlh$gamhlat1abs < 0.2) # remove an outlier
		#summary(lm(I(directionlh$obslat1abs[i] - directionlh$gamhlat1abs[i]) ~ directionlh$tsdur[i]))

	# Multi-model inference
	options(warn=2) # turn warnings to errors
	vars2 = c('Intercept', vars)
	param = numeric(length(vars2))
	se = numeric(length(vars2))
	for(i in 1:length(vars2)){
		print(vars2[i])
		subs = grep(vars2[i], names(mod2))
		if(vars2[i] == 'Intercept') subs = 1:length(mod2)
		for(j in subs){
			if(all(!is.na(mod2[[j]]))){ # skip any models that couldn't be fit
				c = summary(mod2[[j]])$coefficients
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
	
			
	# write out results for top 6 models, only climate velocity, only survey duration, and only survey extent (for paper)
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
	for(i in 1:min(maxmods,length(mod2))){
		v = names(coef(mod2[[modord[i]]])) # variable names in this model
		s2<- signif(coef(mod2[[modord[i]]]),2)
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
		v = names(coef(mod2[[1]])) # variable names in this model
		s2<- signif(coef(mod2[[1]]),2)
		for(j in 1:length(outvars)){
			out$Full[j] = s2[grep(outvarslong[j], v)]
		}
		out$Full[j+1] = signif(aic[1] - min(aic),3)
		out$Full[j+2] = signif(r2[1],2)	
		out$Full[j+3] = signif(aicw[1],2)
	# Averaged model
		v = names(coef(mod2[[1]])) # variable names in this model (same as full model)
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
		i = which(names(mod2) == specialmods[j])
		v = names(coef(mod2[[i]])) # variable names in this model
		s2<- signif(coef(mod2[[i]]),2)
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
		i = which(names(mod2) == paste('y~', outvars[j], sep='')) # find the right model
		v = names(coef(mod2[[i]])) # variable names in this model
		s2<- signif(coef(mod2[[i]]),2)
		out[[nm]][j] = s2[2]
		out[[nm]][length(outvars)+1] = signif(aic[i] - min(aic),3)
		out[[nm]][length(outvars)+2] = signif(r2[i],2)
		out[[nm]][length(outvars)+3] = signif(aicw[i],2)
	}
	# Term weights
	for(i in 1:length(outvars)){
		out$Weights[i] = signif(sum(aicw[grep(outvars[i], names(mod2))], na.rm=TRUE), 3)
	}

	out			

# Write out table of models (Tables 1, S4, and S5)
	write.csv(out, file=paste('Tables/traitmods_', tax, '_', yvar, '_', Sys.Date(), '.csv', sep=''), row.names=FALSE)


