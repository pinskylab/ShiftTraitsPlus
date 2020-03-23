Life History Strategies vs. observed shifts in species geographic ranges
================

``` r
require(visreg) # for partial residuals
```

    ## Loading required package: visreg

``` r
require(AICcmodavg) # for AICcs
```

    ## Loading required package: AICcmodavg

``` r
require(dplyr)
```

    ## Loading required package: dplyr

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
source('lm.aic.all.r')

knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) # tell RStudio to use project root directory as the root for this notebook. Needed since we are storing code in a separate directory.
```

# Multiple Linear Regression of Life History Strategies vs. observed shifts in species geographic ranges

### Only fish

Load
data

``` r
directionlh = read.csv('data/directionlh.csv', row.names=1) # data from Pinsky et al. 2013 (shifts and some traits)
lhs = read.csv('data/sp_LHS.csv',header = T, sep = ";", row.names = NULL) # Life History Strategy axes from Laurene Pecuchet
```

Merge data for fishes

``` r
# First remove invertebrates
directionlh<-directionlh[-c(which(directionlh$fishinvert==FALSE)),]

# Fix some species names in directionlh to match names in lhs
directionlh$species <- as.character(directionlh$species) # first change to character so that we can edit
directionlh$genus <- as.character(directionlh$genus)
directionlh$sciname <- as.character(directionlh$sciname)
inds <- directionlh$genus == 'Theragra' & directionlh$species == 'chalcogramma'
directionlh$genus[inds] <- 'Gadus'
directionlh$species[inds] <- 'chalcogrammus'
directionlh$sciname[inds] <- 'Gadus chalcogrammus'

inds <- directionlh$genus == 'Clupea' & directionlh$species == 'pallasii'
directionlh$species[inds] <- 'pallasii pallasii'
directionlh$sciname[inds] <- 'Clupea pallasii pallasii'


# check LHS scientific names before merging
MatchRecLHS <- match(unique(sort(directionlh$sciname)), unique(sort(lhs$taxon)))
print('Records that do not match')
```

    ## [1] "Records that do not match"

``` r
unique(sort(directionlh$sciname))[is.na(MatchRecLHS)] # 10 non-matching records, most as Genus level
```

    ##  [1] "Artediellus spp."     "Atheresthes spp."     "Bathyraja spp."      
    ##  [4] "Cynoscion spp."       "Gaidropsarus spp."    "Hippoglossoides spp."
    ##  [7] "Lepidopsetta spp."    "Myoxocephalus spp."   "Sebastes spp."       
    ## [10] "Triglops spp."

``` r
# create LHS genus-level axes
lhs$genus <- vapply(strsplit(as.character(lhs$taxon), " "), `[`, 1, FUN.VALUE=character(1)) # genus name in lower case

lhsgen <- aggregate(list(slow_fast = lhs$slow_fast, Equilibrium = lhs$Equilibrium, Periodic = lhs$Periodic, Opportunistic = lhs$Opportunistic), 
                    by = list(taxon = paste0(lhs$genus, ' spp.')), FUN = mean, na.rm= TRUE) # calc average by genus
lhs <- rbind(lhs[, names(lhsgen)], lhsgen) # append to species values
MatchRecLHS <- match(unique(sort(directionlh$sciname)), unique(sort(lhs$taxon))) # recheck matching names
unique(sort(directionlh$sciname))[is.na(MatchRecLHS)] # all match now
```

    ## character(0)

``` r
# merge LHS with shift data
datMerge<-merge(directionlh, lhs, by.x = 'sciname', by.y = 'taxon')
nrow(directionlh) # 259
```

    ## [1] 259

``` r
nrow(datMerge) # 259
```

    ## [1] 259

``` r
datMerge[is.na(datMerge$Equilibrium), c('sciname', 'Equilibrium')] # 10 species without LHS data
```

    ##                      sciname Equilibrium
    ## 74     Gymnocanthus galeatus          NA
    ## 75  Halieutichthys aculeatus          NA
    ## 80  Hemitripterus americanus          NA
    ## 81  Hemitripterus americanus          NA
    ## 82  Hemitripterus americanus          NA
    ## 83      Hemitripterus bolini          NA
    ## 84      Hemitripterus bolini          NA
    ## 85      Hemitripterus bolini          NA
    ## 104  Lagocephalus laevigatus          NA
    ## 166           Peprilus burti          NA
    ## 167     Peprilus triacanthus          NA
    ## 239   Synaphobranchus kaupii          NA
    ## 247      Trichodon trichodon          NA
    ## 257          Zaprora silenus          NA
    ## 258          Zaprora silenus          NA

# Statistical models

``` r
# Add absolute rate of shift
datMerge$obslat1abs = abs(datMerge$obslat1)
datMerge$gamhlat1abs = datMerge$gamhlat1*sign(datMerge$obslat1) # predicted shift rate, signed by whether in the direction of observed shift or not
```

Sample size by region

``` r
table(datMerge$region)
```

    ## 
    ##        AFSC_Aleutians              AFSC_EBS              AFSC_GOA 
    ##                    28                    24                    49 
    ## DFO_Newfoundland_Fall      DFO_ScotianShelf            DFO_SoGulf 
    ##                    41                    16                     6 
    ##          NEFSC_Spring           SEFSC_GOMex         WestCoast_Tri 
    ##                    24                    31                    40

Sample size by factor (ordered by decreasing
size)

``` r
sum(!is.na(datMerge$obslat1)) # 259 (up from 199 fishes in Pinsky et al. 2013)
```

    ## [1] 259

``` r
with(datMerge, sum(!is.na(obslat1) & !is.na(slow_fast))) # 259
```

    ## [1] 259

``` r
with(datMerge, sum(!is.na(obslat1) & !is.na(Equilibrium))) # 244. Missing LHS data for some species
```

    ## [1] 244

``` r
with(datMerge, sum(!is.na(obslat1) & !is.na(Periodic))) # 244
```

    ## [1] 244

``` r
with(datMerge, sum(!is.na(obslat1) & !is.na(Opportunistic))) # 244
```

    ## [1] 244

Set up which variables to
include

``` r
vars <- list(logabs_slowfast = c('gamhlat1abs', 'slow_fast', 'tsdur', 'extentlat'), # slow_fast for log absolute rate of shift
    logabs_lhs = c('gamhlat1abs', 'Equilibrium', 'Opportunistic', 'Periodic', 'tsdur', 'extentlat'), # LHS for log abs
    raw_slowfast = c('gamhlat1', 'slow_fast', 'tsdur', 'extentlat'), # slow_fast for raw rate & direction of shift
    raw_lhs = c('gamhlat1', 'Equilibrium', 'Opportunistic', 'Periodic', 'tsdur', 'extentlat')) # LHS for raw

yvars <- c(logabs_slowfast = 'logabs', logabs_lhs = 'logabs', raw_slowfast = 'raw', raw_lhs = 'raw')
```

Cycle through the possible models and print the results

``` r
for(iter in 1:length(yvars)){
    cat('\n########################\n') # a crude way to delineate the breaks between the four options we run
    cat(paste('Model type:', names(vars)[iter], '\n'))
        
    # Trim to complete cases
    i = complete.cases(datMerge[,vars[[iter]][!grepl(':', vars[[iter]])]])
    cat(paste('Complete cases:', sum(i), '\n'))
    
    
    cat('Colinearity?\n')
    r2 = matrix(NA, nrow=length(vars[[iter]]), ncol=length(vars[[iter]]))
    colnames(r2) = vars[[iter]]; rownames(r2) = vars[[iter]]
    for(k in 1:(length(vars[[iter]])-1)){
        for(j in (k+1):length(vars[[iter]])){
            r2[k,j] = round(cor(as.numeric(datMerge[[vars[[iter]][k]]][i]), as.numeric(datMerge[[vars[[iter]][j]]][i]))^2,2)
        }
    }
    print(r2)   
    
    # All subsets lm        
    if(yvars[iter]=='logabs') mod = lm.aic.all(y=log(datMerge$obslat1abs[i]), x=datMerge[i,vars[[iter]]], models='all') # test all subsets with log(abs shift)
    if(yvars[iter]=='raw') mod = lm.aic.all(y=datMerge$obslat1[i], x=datMerge[i,vars[[iter]]], models='all') # test all subsets with raw shift value
    aic = numeric(length(mod))
    for(i in 1:length(mod)) aic[i] = AIC(mod[[i]])
    aicc = numeric(length(mod)) # small-sample version
    for(i in 1:length(mod)) aicc[i] = AICc(mod[[i]])
    r2 = numeric(length(mod))
    for(i in 1:length(mod)) r2[i] = summary(mod[[i]])$r.squared
    
    modord = order(aic) # order from best to worst
    
    # AIC model weights
    aicd = aic - min(aic) # delta aic
    aicw = exp(-aicd/2)/sum(exp(-aicd/2), na.rm=TRUE) # AIC weights: Burnham & Anderson p. 75 section 2.9.1
    
    # term weights
    cat('\n# Relative Variable Importance (RVI)\n')
    for(i in 1:length(vars[[iter]])){
        cat(paste(format(vars[[iter]][i], width=17), signif(sum(aicw[grep(vars[[iter]][i], names(mod))], na.rm=TRUE), 3), '\n'))
    }
    
    # full model
    cat('\n# Full model\n')
    print(summary(mod[[1]]))
    
    # best model
    cat('\n# Best model\n')
    print(summary(mod[[modord[1]]]))
    cat(paste('\nAIC:', aic[modord[1]]))
    cat(paste('\nAIC weight:', aicw[modord[1]], '\n\n'))
    #par(mfrow=c(1,3)); visreg(mod[[modord[1]]], trans=exp, type='conditional')
    
}
```

    ## 
    ## ########################
    ## Model type: logabs_slowfast 
    ## Complete cases: 259 
    ## Colinearity?
    ##             gamhlat1abs slow_fast tsdur extentlat
    ## gamhlat1abs          NA         0  0.01      0.11
    ## slow_fast            NA        NA  0.00      0.01
    ## tsdur                NA        NA    NA      0.27
    ## extentlat            NA        NA    NA        NA
    ## [1] "Tried 16 models"
    ## 
    ## # Relative Variable Importance (RVI)
    ## gamhlat1abs       1 
    ## slow_fast         0.337 
    ## tsdur             0.952 
    ## extentlat         1 
    ## 
    ## # Full model
    ## 
    ## Call:
    ## lm(formula = formula(forms[j]), data = x[inds, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -4.3497 -0.4220  0.2100  0.6512  2.2628 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -4.46857    0.41339 -10.810  < 2e-16 ***
    ## gamhlat1abs 13.51909    2.21236   6.111 3.71e-09 ***
    ## slow_fast   -0.03593    0.04516  -0.796  0.42702    
    ## tsdur       -0.02899    0.01030  -2.814  0.00528 ** 
    ## extentlat    0.12038    0.02140   5.626 4.85e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.014 on 254 degrees of freedom
    ## Multiple R-squared:  0.3869, Adjusted R-squared:  0.3773 
    ## F-statistic: 40.08 on 4 and 254 DF,  p-value: < 2.2e-16
    ## 
    ## 
    ## # Best model
    ## 
    ## Call:
    ## lm(formula = formula(forms[j]), data = x[inds, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -4.3618 -0.4214  0.2012  0.6466  2.2929 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -4.44937    0.41239 -10.789  < 2e-16 ***
    ## gamhlat1abs 13.63221    2.20620   6.179 2.54e-09 ***
    ## tsdur       -0.02915    0.01029  -2.832    0.005 ** 
    ## extentlat    0.11837    0.02123   5.575 6.28e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.013 on 255 degrees of freedom
    ## Multiple R-squared:  0.3854, Adjusted R-squared:  0.3782 
    ## F-statistic:  53.3 on 3 and 255 DF,  p-value: < 2.2e-16
    ## 
    ## 
    ## AIC: 747.891012411695
    ## AIC weight: 0.63170098008523 
    ## 
    ## 
    ## ########################
    ## Model type: logabs_lhs 
    ## Complete cases: 244 
    ## Colinearity?
    ##               gamhlat1abs Equilibrium Opportunistic Periodic tsdur
    ## gamhlat1abs            NA        0.03          0.02     0.00  0.01
    ## Equilibrium            NA          NA          0.17     0.28  0.00
    ## Opportunistic          NA          NA            NA     0.31  0.10
    ## Periodic               NA          NA            NA       NA  0.07
    ## tsdur                  NA          NA            NA       NA    NA
    ## extentlat              NA          NA            NA       NA    NA
    ##               extentlat
    ## gamhlat1abs        0.11
    ## Equilibrium        0.03
    ## Opportunistic      0.00
    ## Periodic           0.03
    ## tsdur              0.28
    ## extentlat            NA
    ## [1] "Tried 64 models"
    ## 
    ## # Relative Variable Importance (RVI)
    ## gamhlat1abs       1 
    ## Equilibrium       0.404 
    ## Opportunistic     0.302 
    ## Periodic          0.327 
    ## tsdur             0.924 
    ## extentlat         1 
    ## 
    ## # Full model
    ## 
    ## Call:
    ## lm(formula = formula(forms[j]), data = x[inds, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -4.2808 -0.4352  0.1751  0.6480  2.3878 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   -9.274e+01  1.358e+04  -0.007   0.9946    
    ## gamhlat1abs    1.391e+01  2.295e+00   6.061 5.28e-09 ***
    ## Equilibrium    8.790e+01  1.358e+04   0.006   0.9948    
    ## Opportunistic  8.828e+01  1.358e+04   0.007   0.9948    
    ## Periodic       8.838e+01  1.358e+04   0.007   0.9948    
    ## tsdur         -2.890e-02  1.158e-02  -2.496   0.0132 *  
    ## extentlat      1.223e-01  2.282e-02   5.359 1.98e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.027 on 237 degrees of freedom
    ## Multiple R-squared:  0.389,  Adjusted R-squared:  0.3736 
    ## F-statistic: 25.15 on 6 and 237 DF,  p-value: < 2.2e-16
    ## 
    ## 
    ## # Best model
    ## 
    ## Call:
    ## lm(formula = formula(forms[j]), data = x[inds, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -4.3445 -0.4253  0.2049  0.6515  2.3097 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -4.41223    0.43405 -10.165  < 2e-16 ***
    ## gamhlat1abs 13.70227    2.25131   6.086 4.54e-09 ***
    ## tsdur       -0.03024    0.01080  -2.800  0.00552 ** 
    ## extentlat    0.11708    0.02218   5.280 2.89e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.024 on 240 degrees of freedom
    ## Multiple R-squared:  0.3853, Adjusted R-squared:  0.3776 
    ## F-statistic: 50.15 on 3 and 240 DF,  p-value: < 2.2e-16
    ## 
    ## 
    ## AIC: 709.964971889704
    ## AIC weight: 0.255298792024711 
    ## 
    ## 
    ## ########################
    ## Model type: raw_slowfast 
    ## Complete cases: 259 
    ## Colinearity?
    ##           gamhlat1 slow_fast tsdur extentlat
    ## gamhlat1        NA         0  0.01      0.21
    ## slow_fast       NA        NA  0.00      0.01
    ## tsdur           NA        NA    NA      0.27
    ## extentlat       NA        NA    NA        NA
    ## [1] "Tried 16 models"
    ## 
    ## # Relative Variable Importance (RVI)
    ## gamhlat1          1 
    ## slow_fast         0.349 
    ## tsdur             0.294 
    ## extentlat         0.942 
    ## 
    ## # Full model
    ## 
    ## Call:
    ## lm(formula = formula(forms[j]), data = x[inds, ])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.153700 -0.011757  0.000719  0.013139  0.230769 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -0.0170881  0.0148494  -1.151  0.25091    
    ## gamhlat1     0.8026235  0.0795630  10.088  < 2e-16 ***
    ## slow_fast   -0.0013732  0.0015652  -0.877  0.38114    
    ## tsdur        0.0001518  0.0003626   0.419  0.67571    
    ## extentlat    0.0021114  0.0007960   2.653  0.00849 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.03522 on 254 degrees of freedom
    ## Multiple R-squared:  0.4136, Adjusted R-squared:  0.4043 
    ## F-statistic: 44.78 on 4 and 254 DF,  p-value: < 2.2e-16
    ## 
    ## 
    ## # Best model
    ## 
    ## Call:
    ## lm(formula = formula(forms[j]), data = x[inds, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.15304 -0.01226  0.00059  0.01297  0.22750 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -0.0108510  0.0055831  -1.944  0.05305 .  
    ## gamhlat1     0.8094776  0.0779906  10.379  < 2e-16 ***
    ## extentlat    0.0018716  0.0006672   2.805  0.00541 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.03514 on 256 degrees of freedom
    ## Multiple R-squared:  0.4114, Adjusted R-squared:  0.4068 
    ## F-statistic: 89.47 on 2 and 256 DF,  p-value: < 2.2e-16
    ## 
    ## 
    ## AIC: -994.428417865364
    ## AIC weight: 0.436319873656814 
    ## 
    ## 
    ## ########################
    ## Model type: raw_lhs 
    ## Complete cases: 244 
    ## Colinearity?
    ##               gamhlat1 Equilibrium Opportunistic Periodic tsdur extentlat
    ## gamhlat1            NA        0.02          0.03     0.00  0.01      0.21
    ## Equilibrium         NA          NA          0.17     0.28  0.00      0.03
    ## Opportunistic       NA          NA            NA     0.31  0.10      0.00
    ## Periodic            NA          NA            NA       NA  0.07      0.03
    ## tsdur               NA          NA            NA       NA    NA      0.28
    ## extentlat           NA          NA            NA       NA    NA        NA
    ## [1] "Tried 64 models"
    ## 
    ## # Relative Variable Importance (RVI)
    ## gamhlat1          1 
    ## Equilibrium       0.412 
    ## Opportunistic     0.596 
    ## Periodic          0.527 
    ## tsdur             0.298 
    ## extentlat         0.948 
    ## 
    ## # Full model
    ## 
    ## Call:
    ## lm(formula = formula(forms[j]), data = x[inds, ])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.154014 -0.013695  0.000821  0.014068  0.224294 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   -7.599e+02  4.698e+02  -1.617   0.1071    
    ## gamhlat1       7.896e-01  8.257e-02   9.563   <2e-16 ***
    ## Equilibrium    7.598e+02  4.698e+02   1.617   0.1071    
    ## Opportunistic  7.598e+02  4.698e+02   1.617   0.1071    
    ## Periodic       7.602e+02  4.700e+02   1.617   0.1071    
    ## tsdur         -1.130e-04  4.061e-04  -0.278   0.7810    
    ## extentlat      2.162e-03  8.520e-04   2.537   0.0118 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.03568 on 237 degrees of freedom
    ## Multiple R-squared:  0.4281, Adjusted R-squared:  0.4136 
    ## F-statistic: 29.57 on 6 and 237 DF,  p-value: < 2.2e-16
    ## 
    ## 
    ## # Best model
    ## 
    ## Call:
    ## lm(formula = formula(forms[j]), data = x[inds, ])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.154792 -0.014694  0.000175  0.013750  0.225136 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   -0.0041544  0.0071007  -0.585  0.55905    
    ## gamhlat1       0.7758449  0.0814398   9.527  < 2e-16 ***
    ## Opportunistic -0.0216587  0.0116125  -1.865  0.06338 .  
    ## extentlat      0.0021139  0.0007005   3.018  0.00282 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0357 on 240 degrees of freedom
    ## Multiple R-squared:   0.42,  Adjusted R-squared:  0.4128 
    ## F-statistic: 57.93 on 3 and 240 DF,  p-value: < 2.2e-16
    ## 
    ## 
    ## AIC: -927.856897003257
    ## AIC weight: 0.148892343669801
