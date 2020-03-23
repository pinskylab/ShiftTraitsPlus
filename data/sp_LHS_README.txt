From: Aurore Maureaud <auma@aqua.dtu.dk>
Subject: FW: LHS scores
Date: February 25, 2020 at 7:35:38 AM GMT+1
To: Malin Pinsky <malin.pinsky@rutgers.edu>, Martin Lindegren <mli@aqua.dtu.dk>

Hi Malin and Martin,
 
Sorry it took so long, Laurene received all my emails in her junk box & I got a little behind with the thesis submission coming up.
Attached the life-history strategy scores from Laurene, and she also gave me the scores from the fast-slow continuum (based on a simple PCA)
The analyses concern >2000 species.
In the file, OS is offspring size, PC is parental care (and the 1-3 scores correspond to the 3 modalities we include in the trait data)
She also shared her code to re-run the archetypal analysis (see below)
I told her we would invite her for our next skype meeting, hopefully with some interesting results!
 
Best,
Aurore
 
library(archetypes)
 
aa<-stepArchetypes(My_trait_dataset, k=1:10, nrep=10,archetypesFamily(which = c("robust"))) # Trying from 1 to 10 archetypes, with 10 runs for each
 
#### Look at the Goodness of fit -> to determine how many Archetypes are needed to summarize the data
rss<-as.data.frame(rss(aa))
rss$mean<-rowMeans(rss, na.rm=T)
 
rss$se<-apply(rss[,1:10],1,function(x){sd(x, na.rm=T)})
rss$k<-factor(1:10)
 
tiff("RSS_AA.tiff", height = 4, width = 6, units = 'cm', compression = "lzw", res = 400)
ggplot(rss, aes(x=k, y=mean)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.4) +
  geom_line(aes(group = 1), size=0.3) +
  geom_point(size=0.4) +
  ylab("Residuals Sum of Squares (RSS)") +
  xlab("Numbers of Strategies (k)") +
  theme_bw() +
  theme(axis.text.x = element_text(size=6),axis.text.y = element_text(size=6),axis.title.x=element_text(size=6),axis.title.y=element_text(size=6))
dev.off()
 
#### 3 Archetypes 
aa_3<-bestModel(aa[[3]]) # Select the run with 3 archetypes and lowest RSS
barplot(aa_3,SptrA,below.compressed.height = 2) # Trait values of the archetypes
 
coef(aa_3) ## The species LHS scores
 
#### Plot the 3 groups Ternary plot
library(Ternary)
 
tiff("LHS_AA_Ternary.tiff", height = 8, width = 8, units = 'cm', compression = "lzw", res = 600)
par(mar=c(0,0,0,0))
TernaryPlot(alab="Equilibrium \u2192", blab=expression("Opportunistic "*symbol('\256')), clab=expression(symbol('\254')*" Periodic"))
ColourTernary(TernaryDensity(coef(aa_3), resolution=10L))
TernaryPoints(coef(aa_3), col='red', pch=16, cex=0.2)
#TernaryDensityContour(coef(aa_3), resolution=30L)
dev.off()
