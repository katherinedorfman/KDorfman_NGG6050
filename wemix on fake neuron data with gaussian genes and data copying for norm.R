rm(list=ls())

library(lme4) # for the analysis
library(haven) # to load the SPSS .sav file
library(tidyverse) # needed for data manipulation.
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(lmerTest)# to get p-value estimations that are not part of the standard lme4 packages
library(sjstats)
library(WeMix)

#### create fake cell dataset for 1000 nuclei
# 4 samples: 400 nucs from sample 1, 300 from sample 2, 200 from sample 3, 100 from sample 4
# 5 cell types: 500 excitatory neurons, 200 inhibitory neurons, 150 astrocytes, 100 oligs, 50 microglia
# 3 genes: geneA, geneB, and geneC, with expression values ranging from 1 to 5 
# geneA: high in excitatory neurons, geneB: high is astrocytes, geneC: random
cellnumber <- c(1:1000)

cellsample <- vector(length=1000)
cellsample[1:400] <- 1
cellsample[401:700] <- 2
cellsample[701:900] <- 3
cellsample[901:1000] <- 4

celltype <- vector(length=1000)
celltype[1:500] <- "excitatory neuron"
celltype[501:700] <- "inhibitory neuron"
celltype[701:850] <- "astrocyte"
celltype[851:950] <- "oligo"
celltype[951:1000] <- "microglia"
celltype <- sample(celltype)

geneA <- vector(length=1000)
geneA <- rnorm(geneA, mean=2.5, sd=1)
geneA[which(celltype=='excitatory neuron')] <- 10
geneB <- vector(length=1000)
geneB <- rnorm(geneB, mean=2.5, sd=1)
geneB[which(celltype=='astrocyte')] <-4
geneC <- vector(length=1000)
geneC <- rnorm(geneC, mean=2.5, sd=1)

fakeneurondata <- data.frame(cellnumber,cellsample,celltype,geneA,geneB,geneC)


#### wemix
# create weights
fakeneurondata$cellsampleweights <- vector(length=1000)
fakeneurondata$cellsampleweights[which(fakeneurondata$cellsample==1)] <- 1/(sum(fakeneurondata$cellsample==1)/length(fakeneurondata$cellsample))
fakeneurondata$cellsampleweights[which(fakeneurondata$cellsample==2)] <- 1/(sum(fakeneurondata$cellsample==2)/length(fakeneurondata$cellsample))
fakeneurondata$cellsampleweights[which(fakeneurondata$cellsample==3)] <- 1/(sum(fakeneurondata$cellsample==3)/length(fakeneurondata$cellsample))
fakeneurondata$cellsampleweights[which(fakeneurondata$cellsample==4)] <- 1/(sum(fakeneurondata$cellsample==4)/length(fakeneurondata$cellsample))

fakeneurondata$celltypeweights <- vector(length=1000)
fakeneurondata$celltypeweights[which(fakeneurondata$celltype=='excitatory neuron')] <- 1/(sum(fakeneurondata$celltype=='excitatory neuron')/length(fakeneurondata$celltype))
fakeneurondata$celltypeweights[which(fakeneurondata$celltype=='inhibitory neuron')] <- 1/(sum(fakeneurondata$celltype=='inhibitory neuron')/length(fakeneurondata$celltype))
fakeneurondata$celltypeweights[which(fakeneurondata$celltype=='astrocyte')] <- 1/(sum(fakeneurondata$celltype=='astrocyte')/length(fakeneurondata$celltype))
fakeneurondata$celltypeweights[which(fakeneurondata$celltype=='oligo')] <- 1/(sum(fakeneurondata$celltype=='oligo')/length(fakeneurondata$celltype))
fakeneurondata$celltypeweights[which(fakeneurondata$celltype=='microglia')] <- 1/(sum(fakeneurondata$celltype=='microglia')/length(fakeneurondata$celltype))


## Plotting
# plot neuron type vs geneA
ggplot(data  = fakeneurondata,
       aes(x = celltype,
           y = geneA),
       title = 'geneA')+
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  geom_smooth(method = lm,
              se     = FALSE, 
              col    = "black",
              size   = .5, 
              alpha  = .8)+ # to add regression line
  ylim(0,10)+
  theme_minimal()
# plot neuron type vs geneB
ggplot(data  = fakeneurondata,
       aes(x = celltype,
           y = geneB),
       title = 'geneB')+
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  geom_smooth(method = lm,
              se     = FALSE, 
              col    = "black",
              size   = .5, 
              alpha  = .8)+ # to add regression line
  ylim(0,10)+
  theme_minimal()
# plot neuron type vs geneC
ggplot(data  = fakeneurondata,
       aes(x = celltype,
           y = geneC),
       title = 'geneC')+
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  geom_smooth(method = lm,
              se     = FALSE, 
              col    = "black",
              size   = .5, 
              alpha  = .8)+ # to add regression line
  ylim(0,10)+
  theme_minimal()

data <- data.frame(celltype,geneA,geneB,geneC)
boxplot(geneA ~ celltype)

library(MASS) 
library(reshape2) 
library(reshape) 
data <- melt(data, id="celltype")

ggplot(data, aes(x = celltype, y = value, color = variable)) +  # ggplot function
  geom_boxplot()

# means
mean(fakeneurondata$geneA)
weighted.mean(fakeneurondata$geneA,fakeneurondata$cellsampleweights)
weighted.mean(fakeneurondata$geneA,fakeneurondata$celltypeweights)
mean(subset(fakeneurondata$geneA, celltype %in% 'astrocyte'))
mean(subset(fakeneurondata$geneA, celltype %in% 'excitatory neuron'))
mean(subset(fakeneurondata$geneA, celltype %in% 'inhibitory neuron'))
mean(subset(fakeneurondata$geneA, celltype %in% 'microglia'))
mean(subset(fakeneurondata$geneA, celltype %in% 'oligo'))

mean(fakeneurondata$geneB)
weighted.mean(fakeneurondata$geneB,fakeneurondata$cellsampleweights)
weighted.mean(fakeneurondata$geneB,fakeneurondata$celltypeweights)
mean(fakeneurondata$geneC)
weighted.mean(fakeneurondata$geneC,fakeneurondata$cellsampleweights)
weighted.mean(fakeneurondata$geneC,fakeneurondata$celltypeweights)


### Copy data
data_copy <- fakeneurondata
data_copy$celltype[1:length(data_copy$celltype)] <- "aaa"
fakeneurondata <- rbind(fakeneurondata, data_copy)

# model with one random effect
# For the weights argument, the weights are specified as a vector with the
# first element giving the name of the weights for level 1 units, the second element giving the name of the weights for
# the level 2 groups, and—when fitting a model with three levels—the third element giving the name of the weights
# for the level 3 groups.
# sort in ascending weights order: fakeneurondata <- fakeneurondata[order(fakeneurondata$celltypeweights), ]

## LMER: write code to run through all genes and save outputs
#number_genes = 3
#lmermodel1 <- list()

#for (i in 1:number_genes){
#  f <- colnames(fakeneurondata)[(number_genes+i)]
#  formula <- paste(f,'~ celltype + (1|cellsample)')
#  formula <- formula(formula)
#  lmermodel1[[i]] <- lmer(formula = formula, 
#                          data = fakeneurondata)
#  
#} 



### LMER MODELS

lmermodel1A <- lmer(formula = geneA ~ celltype + (1|cellsample), 
                    data = fakeneurondata)
summary(lmermodel1A)
a<-fixef(lmermodel1A) #extract fixed effects
which(as.vector(a) > 0)
b<-ranef(lmermodel1A) #extract random effects

lmermodel2A <- lmer(formula = geneB ~ celltype + (1|celltype), 
                    data = fakeneurondata)
summary(lmermodel2A)
lmermodel1Aw <- lmer(formula = geneA ~ celltype + (1|cellsample),
                     data = fakeneurondata,
                     weights = fakeneurondata$celltypeweights)
summary(lmermodel1Aw)
lmermodel2Aw <- lmer(formula = geneB ~ celltype + (1|cellsample),
                     data = fakeneurondata,
                     weights = fakeneurondata$celltypeweights)
summary(lmermodel2Aw)
lmermodel3Aw <- lmer(formula = geneC ~ celltype + (1|cellsample),
                     data = fakeneurondata,
                     weights = fakeneurondata$celltypeweights)
summary(lmermodel3Aw)


# wemix
mixmodel1A <- mix(geneA ~ celltype + (1|cellsample), data=fakeneurondata, 
                  weights = c("celltypeweights","cellsampleweights"))
summary(mixmodel1A)
mixmodel2A <- mix(geneB ~ celltype + (1|cellsample), data=fakeneurondata, 
                  weights = c("celltypeweights","cellsampleweights"))
summary(mixmodel2A)
mixmodel3A <- mix(geneC ~ celltype + (1|cellsample), data=fakeneurondata, 
                  weights = c("celltypeweights","cellsampleweights"))
summary(mixmodel3A)



fakeneurondata_alph <- fakeneurondata$
  
  
  mixmodel1B <- mix(geneA ~ celltype + (1|cellsample), data=fakeneurondata, 
                    weights = c("celltypeweights","cellsampleweights"))
summary(mixmodel1B)
mixmodel1C <- mix(geneC ~ celltype + (1|cellsample), data=fakeneurondata, 
                  weights = c("celltypeweights","cellsampleweights"))
summary(mixmodel1C)

mixmodel2A <- mix(geneA ~ celltype + (1|cellsample), data=fakeneurondata, 
                  weights = c('celltypeweights',"W2"))
summary(mixmodel2A)
mixmodel2B <- mix(geneB ~ celltype + (1|cellsample), data=fakeneurondata, 
                  weights = c('celltypeweights',"W2"))
summary(mixmodel2B)
