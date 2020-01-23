# Tests all possible multivariate regressions of y on x
# DOES account for differences in sample size, so models are comparable
# x must be a dataframe with column names
# log is TRUE or FALSE to indicate which columns in x should be log-transformed
# Does not examine interactions

lm.aic.all = function(y, x, log=FALSE, models = 'best'){
	vars = names(x)
	if(length(log) != length(vars)){
		log = rep(FALSE, length(vars))
	}
	
	# Prep log-transform for formula
	varsmod = as.character(vars)
	for(i in 1:length(varsmod)){
		if(log[i]) varsmod[i] = paste('log(', vars[i], ')', sep='')
	}
	
	# Make the formulas in every possible length
	forms = character(0)
	nterms = character(0)
	for(i in length(varsmod):1){
		m = combn(varsmod, i)
		for(j in 1:ncol(m)){
			forms = c(forms, paste('y~', paste(m[,j], collapse='+'), sep=''))
			nterms = c(nterms, i)
		}
	}
	
	# Add the null model
	forms = c(forms, 'y~1')

	# Set the data trimming
	inds = complete.cases(x)
	y = y[inds]

	# Fit the models
	mods = list(0)
	aic = numeric(0)
	for(j in 1:length(forms)){
		temp = lm(formula(forms[j]), data=x[inds,])
		mods[[j]] = temp
		aic = c(aic, AIC(temp))
	}
	names(mods) = forms
	best = which.min(aic)	
	print(paste('Tried', length(forms), 'models'))

	if(models == 'best') return(mods[[best]])
	if(models == 'all') return(mods)
}