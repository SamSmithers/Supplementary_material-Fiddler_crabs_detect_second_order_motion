# This is a copy of the R script used to plot the results and perform the statistical analysis for Smithers et al.
# "Fiddler crabs detect second-order motion in both intensity and polarization". For an explanation of the statistical 
# analysis conducted, please refer to the statistical analysis section of the methods in the main manuscript.

# Script written by Dr Samuel P. Smithers, Northeastern University, 2022
# Last edited August 2023

# Corresponding authors SPS (s.smithers@northeastern.edu) and NWR (nicholas.roberts@bristol.ac.uk)

# Included in this script:
# - Code used to generate Figure 3 from the main manuscript.
# - Code used for the statistical analysis of data from the intensity and polarization experiments. 

# Session info:
# R version 4.2.2 (2022-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] Hmisc_5.1-0           cowplot_1.1.1         ggplot2_3.4.0         bannerCommenter_1.0.0 DHARMa_0.4.6         
# [6] dplyr_1.0.10          lme4_1.1-31           Matrix_1.5-1         
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.9         lattice_0.20-45    digest_0.6.30      foreach_1.5.2      utf8_1.2.2        
# [6] mime_0.12          R6_2.5.1           plyr_1.8.8         backports_1.4.1    gap.datasets_0.0.5
# [11] evaluate_0.18      pillar_1.8.1       rlang_1.0.6        data.table_1.14.4  minqa_1.2.5       
# [16] rstudioapi_0.15.0  nloptr_2.0.3       rpart_4.1.19       gap_1.3-1          checkmate_2.2.0   
# [21] rmarkdown_2.23     splines_4.2.2      stringr_1.4.1      foreign_0.8-83     htmlwidgets_1.5.4 
# [26] munsell_0.5.0      shiny_1.7.3        compiler_4.2.2     httpuv_1.6.6       xfun_0.39         
# [31] pkgconfig_2.0.3    base64enc_0.1-3    mgcv_1.8-41        htmltools_0.5.3    nnet_7.3-18       
# [36] tidyselect_1.2.0   tibble_3.1.8       gridExtra_2.3      htmlTable_2.4.1    codetools_0.2-18  
# [41] fansi_1.0.3        crayon_1.5.2       withr_2.5.0        later_1.3.0        MASS_7.3-58.1     
# [46] grid_4.2.2         nlme_3.1-160       xtable_1.8-4       gtable_0.3.1       lifecycle_1.0.3   
# [51] magrittr_2.0.3     scales_1.2.1       cli_3.4.1          stringi_1.7.8      farver_2.1.1      
# [56] promises_1.2.0.1   doParallel_1.0.17  ellipsis_0.3.2     generics_0.1.3     vctrs_0.5.0       
# [61] boot_1.3-28        qgam_1.3.4         Formula_1.2-5      iterators_1.0.14   tools_4.2.2       
# [66] glue_1.6.2         parallel_4.2.2     fastmap_1.1.0      colorspace_2.0-3   cluster_2.1.4     
# [71] knitr_1.41        


##Required dependencies/packages. 
#For stats
library(lme4)
library(dplyr)
library(DHARMa) # For residual diagnostics.
#Other
library(bannerCommenter)

rm(list=ls(all=TRUE))

setwd("***Pathway to folder containing this R script and 'Data from intensity experiment.csv' and 'Data from polarization experiment.csv' raw data files***")

###########################################################################
###########################################################################
###                                                                     ###
###                           INTENSITY EXPERIMENT                      ###
###                                                                     ###
###########################################################################
###########################################################################

banner("Intensity experiment:", emph = TRUE)

#Load in data
I_data.with.NA<-read.csv("Data from intensity experiment.csv",header=T)

#Remove rows containing NAs (i.e. remove trials that were rejected)
I_data <-na.omit(I_data.with.NA)

#Make binary response variable a factor
I_data$Response <-factor(I_data$BinaryRes)

#Check data
str(I_data)

##================================================================
##    Mixed effects binary logistic regression for intensity     =
##================================================================

boxup("Mixed effects binary logistic regression for intensity", bandChar = "=")

#Set up full model for intensity experiment
m1 <- glmer(Response ~ Stimulus + Sex + Order + Size_mm + (1|CrabID) , family="binomial", data=I_data)
#Model fails to converge properly so we try a different optimizer. 
ss <- getME(m1,c("theta","fixef"))
m1_1 <- update(m1,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
#The model now converges properly. 

# Plot model residual and check residual diagnostics.
simulationOutput <- simulateResiduals(fittedModel = m1_1, plot = T)
# Both plots look good and diagnostic tests are all non significant.  

# Test significant of the fixed effects using a likelihood ratio rest to compare the full model 
# with the same model but with the fixed effect of interest removed. 
Stim<-update(m1,~.-Stimulus)
anova(m1_1, Stim, test = "LRT") #Chisq(4)= 62.944, p<0.001**

Sx <-update(m1,~.-Sex)
anova(m1_1, Sx, test = "LRT") #Chisq(1)= 2.9489, p= 0.08594

Sz <-update(m1,~.-Size_mm)
anova(m1_1, Sz, test = "LRT") #Chisq(1)= 3.2608, p= 0.07095

Ord <-update(m1,~.-Order)
anova(m1_1, Ord, test = "LRT") #Chisq(1)= 0.0442, p=0.8335

##===============================================================
##            Pairwise McNemar tests for intensity              =
##===============================================================

boxup("Pairwise McNemar tests for intensity ", bandChar = "=")

### Control vs first order ###
# Get numbers to put in table needed to run the McNemar's test
pairwiseData <- subset(I_data, Stimulus == "Control" | Stimulus == 'First-order stimulus')

#remove all crabs that do not contribute data to both stimuli
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 

#Extract data to build table
Stim_1 <- subset(pairwiseData, Stimulus == "Control")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "First-order stimulus")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID

#Check each response in a row are both from the same individual
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE

#Make data frame with response to each stimulus as a separate column 
ResponseTable <- data.frame(Control = Stim_1_Res, FO = Stim_2_Res)

#Convert ResponseTable into a contingency table
contingencyTable <- table(ResponseTable) 
contingencyTable

#Perform McNemar test
mcnemar.test(contingencyTable, y = NULL, correct = TRUE) # A continuity correction (correct = TRUE) is applied if any of the counts are <5.
#McNemar's chi-squared = 19.531, df = 1, p-value = 9.897e-06

#Use the contingency table to calculate odds ratio. 
# Considering a 2 x 2 table, with a and d being the concordant cells and b and c being the discordant cells, the odds ratio is simply the greater of (b/c) or (c/b) (https://rcompanion.org/handbook/H_05.html)
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3])) 
# OR = 9.666667

### Control vs mean grey ###
pairwiseData <- subset(I_data, Stimulus == "Control" | Stimulus == 'Mean grey')
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 
Stim_1 <- subset(pairwiseData, Stimulus == "Control")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "Mean grey")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE
ResponseTable <- data.frame(Control = Stim_1_Res, MeanGrey = Stim_2_Res)
contingencyTable <- table(ResponseTable) 
contingencyTable
mcnemar.test(contingencyTable, y = NULL, correct = TRUE)
#McNemar's chi-squared = 15.613, df = 1, p-value = 7.772e-05
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3]))
# OR = 6.75

#Control vs flicker stim
pairwiseData <- subset(I_data, Stimulus == "Control" | Stimulus == 'Second-order flicker stimulus')
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 
Stim_1 <- subset(pairwiseData, Stimulus == "Control")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "Second-order flicker stimulus")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE
ResponseTable <- data.frame(Control = Stim_1_Res, Flicker = Stim_2_Res)
contingencyTable <- table(ResponseTable) 
contingencyTable
mcnemar.test(contingencyTable, y = NULL, correct = TRUE)
#McNemar's chi-squared = 31.61, df = 1, p-value = 1.885e-08
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3]))
# OR = 19.5

#Control vs flicker control
pairwiseData <- subset(I_data, Stimulus == "Control" | Stimulus == 'Flicker control')
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 
Stim_1 <- subset(pairwiseData, Stimulus == "Control")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "Flicker control")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE
ResponseTable <- data.frame(Control = Stim_1_Res,  Flicker_control= Stim_2_Res)
contingencyTable <- table(ResponseTable) 
contingencyTable
mcnemar.test(contingencyTable, y = NULL, correct = TRUE)
# McNemar's chi-squared = 9.3333, df = 1, p-value = 0.00225
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3]))
# OR = 6

#First order vs mean grey
pairwiseData <- subset(I_data, Stimulus == "First-order stimulus" | Stimulus == 'Mean grey')
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 
Stim_1 <- subset(pairwiseData, Stimulus == "First-order stimulus")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "Mean grey")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE
ResponseTable <- data.frame(FO = Stim_1_Res,  midgrey= Stim_2_Res)
contingencyTable <- table(ResponseTable) 
contingencyTable
mcnemar.test(contingencyTable, y = NULL, correct = FALSE)
#McNemar's chi-squared = 0.36, df = 1, p-value = 0.5485 WITHOUT continuity correction (since none of the counts are below 5)
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3]))
# OR = 1.272727

#First order vs flicker stim
pairwiseData <- subset(I_data, Stimulus == "First-order stimulus" | Stimulus == 'Second-order flicker stimulus')
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 
Stim_1 <- subset(pairwiseData, Stimulus == "First-order stimulus")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "Second-order flicker stimulus")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE
ResponseTable <- data.frame(FO = Stim_1_Res,  Flicker= Stim_2_Res)
contingencyTable <- table(ResponseTable) 
contingencyTable
mcnemar.test(contingencyTable, y = NULL, correct = TRUE)
#McNemar's chi-squared = 4.5, df = 1, p-value = 0.03389
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3]))
# OR = 3.5

#First order vs flicker control
pairwiseData <- subset(I_data, Stimulus == "First-order stimulus" | Stimulus == 'Flicker control')
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 
Stim_1 <- subset(pairwiseData, Stimulus == "First-order stimulus")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "Flicker control")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE
ResponseTable <- data.frame(FO = Stim_1_Res,  Flicker_control= Stim_2_Res)
contingencyTable <- table(ResponseTable) 
contingencyTable
mcnemar.test(contingencyTable, y = NULL, correct = FALSE)
#McNemar's chi-squared = 5.2609, df = 1, p-value = 0.02181 WITHOUT continuity correction (since none of the counts are below 5)
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3]))
# OR = 2.833333

#mean grey vs flicker stim
pairwiseData <- subset(I_data, Stimulus == "Mean grey" | Stimulus == 'Second-order flicker stimulus')
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 
Stim_1 <- subset(pairwiseData, Stimulus == "Mean grey")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "Second-order flicker stimulus")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE
ResponseTable <- data.frame(Midgrey = Stim_1_Res,  Flicker = Stim_2_Res)
contingencyTable <- table(ResponseTable) 
contingencyTable
mcnemar.test(contingencyTable, y = NULL, correct = TRUE)
#McNemar's chi-squared = 6.8571, df = 1, p-value = 0.008829
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3]))
# OR = 4.25

#mean grey vs flicker control
pairwiseData <- subset(I_data, Stimulus == "Mean grey" | Stimulus == 'Flicker control')
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 
Stim_1 <- subset(pairwiseData, Stimulus == "Mean grey")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "Flicker control")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE
ResponseTable <- data.frame(Midgrey = Stim_1_Res,  Flicker_control= Stim_2_Res)
contingencyTable <- table(ResponseTable) 
contingencyTable
mcnemar.test(contingencyTable, y = NULL, correct = FALSE)
# McNemar's chi-squared = 2.6667, df = 1, p-value = 0.1025 WITHOUT continuity correction (since none of the counts are below 5)
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3]))
# OR = 2

#flicker stim vs flicker control
pairwiseData <- subset(I_data, Stimulus == "Second-order flicker stimulus" | Stimulus == 'Flicker control')
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 
Stim_1 <- subset(pairwiseData, Stimulus == "Second-order flicker stimulus")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "Flicker control")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE
ResponseTable <- data.frame(Flicker = Stim_1_Res,  Flicker_control= Stim_2_Res)
contingencyTable <- table(ResponseTable) 
contingencyTable
mcnemar.test(contingencyTable, y = NULL, correct = TRUE)
# McNemar's chi-squared = 15.75, df = 1, p-value = 7.229e-05
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3]))
# OR = 8.333333

##### INTESNITY: APPLY BONFERRONI CORRECTION TO P VALUES FROM PAIRWISE MCNEMAR'S TESTS ####
# List all of the p-values
P<-c(9.897e-06, 7.772e-05, 1.885e-08, 0.00225, 0.5485, 0.03389, 0.02181, 0.008829, 0.1025, 7.229e-05)
P <- sort(P, decreasing = FALSE)

# Apply Bonferroni correct
P_adjusted <- p.adjust(P, method = "bonferroni", n = length(P))
P_adjusted

# Copy of test results from above showing adjusted P values. 
#Control vs first order: McNemar's chi-squared = 19.531, df = 1, p-value = 9.897e-06, p-adjusted = 9.897e-05*** 
#Control vs mean grey: McNemar's chi-squared = 15.613, df = 1, p-value = 7.772e-05, p-adjusted = 0.0007772***
#Control vs flicker stim: McNemar's chi-squared = 31.61, df = 1, p-value = 1.885e-08, p-adjusted = 1.885e-07***
#Control vs flicker control: McNemar's chi-squared = 9.3333, df = 1, p-value = 0.00225, p-adjusted = 0.0225*
#First order vs mean grey: McNemar's chi-squared = 0.36, df = 1, p-value = 0.5485, p-adjusted = 1 
#First order vs flicker stim: McNemar's chi-squared = 4.5, df = 1, p-value = 0.03389, p-adjusted = 0.3389 
#First order vs flicker control: McNemar's chi-squared = 5.2609, df = 1, p-value = 0.02181, p-adjusted = 0.2181
#mean grey vs flicker stim: McNemar's chi-squared = 6.8571, df = 1, p-value = 0.008829, p-adjusted = 0.08829 
#mean grey vs flicker control: McNemar's chi-squared = 2.6667, df = 1, p-value = 0.1025, p-adjusted = 1
#flicker stim vs flicker control: McNemar's chi-squared = 15.75, df = 1, p-value = 7.229e-05, p-adjusted = 0.0007229** 

############################################################################
############################################################################
###                                                                      ###
###                       POLARIZATION EXPERIMENT:                       ###
###                                                                      ###
############################################################################
############################################################################

banner("Polarization experiment:", emph = TRUE)

rm(list=ls(all=TRUE))
#Load in data
P_data.with.NA<-read.csv("Data from polarization experiment.csv",header=T)

#Remove rows containing NAs (i.e. remove trials that were rejected)
P_data <-na.omit(P_data.with.NA)

#Make binary response variable a factor
P_data$Response <-factor(P_data$BinaryRes)

#Check data
str(P_data)

##===============================================================
##  Mixed effects binary logistic regression for polarization   =
##===============================================================

boxup("Mixed effects binary logistic regression for polarization", bandChar = "=")

m1 <- glmer(Response ~ Stimulus + Sex + Size_mm + Order + (1|CrabID) , family="binomial", data=P_data)
#Mode fails to converge properly so we try a different optimizer.
ss <- getME(m1,c("theta","fixef"))
m1_1 <- update(m1,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
#The model now converges properly. 

#Plot model residual and check residual diagnostics.
simulationOutput <- simulateResiduals(fittedModel = m1_1, plot = T)
# Both plots look good and diagnostic tests are all non significant. 

# Test significant of the fixed effects using a likelihood ratio rest to compare the full model 
# with the same model but with the fixed effect of interest removed. 
Stim <-update(m1,~.-Stimulus)
anova(m1_1, Stim, test = "LRT") #Chisq(4)= 45.729, p<0.001**

Sx <-update(m1,~.-Sex)
anova(m1_1, Sx, test = "LRT") #Chisq(1)= 2.5371, p=0.1112

Sz <-update(m1,~.-Size_mm)
anova(m1_1, Sz, test = "LRT") #Chisq(1)= 0.5303, p=0.4665

Ord <-update(m1,~.-Order)
anova(m1_1, Ord, test = "LRT") #Chisq(1)= 0.7312, p=0.3925

##================================================================
##            Pairwise McNemar tests for polarization            =
##================================================================

boxup("Pairwise McNemar tests for polarization ", bandChar = "=")

### Control vs first order
# Get numbers to put in table needed to run the McNemar's test
pairwiseData <- subset(P_data, Stimulus == "Control" | Stimulus == 'First-order stimulus')

#remove all crabs that do not contribute data to both stimuli
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 

#Extract data to build table
Stim_1 <- subset(pairwiseData, Stimulus == "Control")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "First-order stimulus")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID

#Check each response in a row are both from the same individual
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE

#Make data frame with response to each stimulus as a separate column 
ResponseTable <- data.frame(Control = Stim_1_Res, FO = Stim_2_Res)

#Convert ResponseTable into a contingency table
contingencyTable <- table(ResponseTable) 
contingencyTable

#Perform McNemar test
mcnemar.test(contingencyTable, y = NULL, correct = TRUE) # A continuity correction (correct = TRUE) is applied if any of the counts are <5.
#McNemar's chi-squared = 25.037, df = 1, p-value = 5.624e-07

#Use the contingency table to calculate odds ratio. 
# Considering a 2 x 2 table, with a and d being the concordant cells and b and c being the discordant cells, the odds ratio is simply the greater of (b/c) or (c/b) (https://rcompanion.org/handbook/H_05.html)
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3])) 
# OR = Inf

#Control vs mean grey'
pairwiseData <- subset(P_data, Stimulus == "Control" | Stimulus == 'Mean grey')
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 
Stim_1 <- subset(pairwiseData, Stimulus == "Control")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "Mean grey")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE
ResponseTable <- data.frame(Control = Stim_1_Res, MeanGrey = Stim_2_Res)
contingencyTable <- table(ResponseTable) 
contingencyTable
mcnemar.test(contingencyTable, y = NULL, correct = TRUE)
#McNemar's chi-squared = 4.3478, df = 1, p-value = 0.03706
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3]))
# OR = 2.833333

#Control vs flicker stim
pairwiseData <- subset(P_data, Stimulus == "Control" | Stimulus == 'Second-order flicker stimulus')
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 
Stim_1 <- subset(pairwiseData, Stimulus == "Control")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "Second-order flicker stimulus")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE
ResponseTable <- data.frame(Control = Stim_1_Res, Flicker = Stim_2_Res)
contingencyTable <- table(ResponseTable) 
contingencyTable
mcnemar.test(contingencyTable, y = NULL, correct = TRUE)
#McNemar's chi-squared = 16, df = 1, p-value = 6.334e-05
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3]))
# OR = 11.5

#Control vs flicker control
pairwiseData <- subset(P_data, Stimulus == "Control" | Stimulus == 'Flicker control')
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 
Stim_1 <- subset(pairwiseData, Stimulus == "Control")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "Flicker control")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE
ResponseTable <- data.frame(Control = Stim_1_Res,  Flicker_control= Stim_2_Res)
contingencyTable <- table(ResponseTable) 
contingencyTable
mcnemar.test(contingencyTable, y = NULL, correct = TRUE)
#McNemar's chi-squared = 0.083333, df = 1, p-value = 0.7728
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3]))
# OR = 1.4

#First order vs mean grey
pairwiseData <- subset(P_data, Stimulus == "First-order stimulus" | Stimulus == 'Mean grey')
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 
Stim_1 <- subset(pairwiseData, Stimulus == "First-order stimulus")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "Mean grey")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE
ResponseTable <- data.frame(FO = Stim_1_Res,  midgrey= Stim_2_Res)
contingencyTable <- table(ResponseTable) 
contingencyTable
mcnemar.test(contingencyTable, y = NULL, correct = FALSE)
#McNemar's chi-squared = 8, df = 1, p-value = 0.004678 WITHOUT continuity correction (since none of the counts are below 5)
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3]))
# OR = 3

#First order vs flicker stim
pairwiseData <- subset(P_data, Stimulus == "First-order stimulus" | Stimulus == 'Second-order flicker stimulus')
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 
Stim_1 <- subset(pairwiseData, Stimulus == "First-order stimulus")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "Second-order flicker stimulus")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE
ResponseTable <- data.frame(FO = Stim_1_Res,  Flicker= Stim_2_Res)
contingencyTable <- table(ResponseTable) 
contingencyTable
mcnemar.test(contingencyTable, y = NULL, correct = FALSE)
#McNemar's chi-squared = 1.6364, df = 1, p-value = 0.2008 WITHOUT continuity correction (since none of the counts are below 5)
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3]))
# OR = 1.75

#First order vs flicker control
pairwiseData <- subset(P_data, Stimulus == "First-order stimulus" | Stimulus == 'Flicker control')
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 
Stim_1 <- subset(pairwiseData, Stimulus == "First-order stimulus")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "Flicker control")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE
ResponseTable <- data.frame(FO = Stim_1_Res,  Flicker_control= Stim_2_Res)
contingencyTable <- table(ResponseTable) 
contingencyTable
mcnemar.test(contingencyTable, y = NULL, correct = TRUE)
#McNemar's chi-squared = 17.455, df = 1, p-value = 2.943e-05
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3]))
# OR = 7.25

#mean grey vs flicker stim
pairwiseData <- subset(P_data, Stimulus == "Mean grey" | Stimulus == 'Second-order flicker stimulus')
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 
Stim_1 <- subset(pairwiseData, Stimulus == "Mean grey")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "Second-order flicker stimulus")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE
ResponseTable <- data.frame(Midgrey = Stim_1_Res,  Flicker = Stim_2_Res)
contingencyTable <- table(ResponseTable) 
contingencyTable
mcnemar.test(contingencyTable, y = NULL, correct = FALSE)
#McNemar's chi-squared = 3.125, df = 1, p-value = 0.0771 WITHOUT continuity correction (since none of the counts are below 5)
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3]))
# OR = 1.909091

#mean grey vs flicker control
pairwiseData <- subset(P_data, Stimulus == "Mean grey" | Stimulus == 'Flicker control')
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 
Stim_1 <- subset(pairwiseData, Stimulus == "Mean grey")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "Flicker control")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE
ResponseTable <- data.frame(Midgrey = Stim_1_Res,  Flicker_control= Stim_2_Res)
contingencyTable <- table(ResponseTable) 
contingencyTable
mcnemar.test(contingencyTable, y = NULL, correct = TRUE)
#McNemar's chi-squared = 3.3684, df = 1, p-value = 0.06646
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3]))
# OR = 2.8

#flicker stim vs flicker control
pairwiseData <- subset(P_data, Stimulus == "Second-order flicker stimulus" | Stimulus == 'Flicker control')
pairwiseData <- pairwiseData %>% group_by(CrabID) %>% filter(n()>1) 
Stim_1 <- subset(pairwiseData, Stimulus == "Second-order flicker stimulus")
Stim_1_Res <- Stim_1$BinaryRes
Stim_1_CrabID <- Stim_1$CrabID
Stim_2 <- subset(pairwiseData, Stimulus == "Flicker control")
Stim_2_Res <- Stim_2$BinaryRes
Stim_2_CrabID <- Stim_2$CrabID
table(Stim_1_CrabID == Stim_2_CrabID) #Check that this is TRUE
ResponseTable <- data.frame(Flicker = Stim_1_Res,  Flicker_control= Stim_2_Res)
contingencyTable <- table(ResponseTable) 
contingencyTable
mcnemar.test(contingencyTable, y = NULL, correct = TRUE)
#McNemar's chi-squared = 12, df = 1, p-value = 0.000532
max((contingencyTable[3]/contingencyTable[2]), (contingencyTable[2]/contingencyTable[3]))
# OR = 5.75

##### POLARIZATION APPLY BONFERRONI CORRECTION TO P VALUES FROM PAIRWISE MCNEMAR'S TESTS ####
#List all of the p-values
P<-c(5.624e-07, 0.03706, 6.334e-05, 0.7728, 0.004678, 0.2008, 2.943e-05, 0.0771, 0.06646, 0.000532)
P <- sort(P, decreasing = FALSE)

# Apply Bonferroni correct
P_adjusted <- p.adjust(P, method = "bonferroni", n = length(P))
P_adjusted

# Copy of test results from above showing adjusted P values. 
#Control vs first order: McNemar's chi-squared = 25.037, df = 1, p-value = 5.624e-07, p-adjusted = 5.624e-06
#Control vs mean grey: McNemar's chi-squared = 4.3478, df = 1, p-value = 0.03706, p-adjusted = 0.3706
#Control vs flicker stim: McNemar's chi-squared = 16, df = 1, p-value = 6.334e-05, p-adjusted = 0.0006334
#Control vs flicker control: McNemar's chi-squared = 0.083333, df = 1, p-value = 0.7728, p-adjusted = 1 
#First order vs mean grey: McNemar's chi-squared = 8, df = 1, p-value = 0.004678, p-adjusted = 0.04678
#First order vs flicker stim: McNemar's chi-squared = 1.6364, df = 1, p-value = 0.2008, p-adjusted = 1
#First order vs flicker control: McNemar's chi-squared = 17.455, df = 1, p-value = 2.943e-05, p-adjusted = 0.0002943
#mean grey vs flicker stim: McNemar's chi-squared = 3.125, df = 1, p-value = 0.0771, p-adjusted = 0.771
#mean grey vs flicker control: McNemar's chi-squared = 3.3684, df = 1, p-value = 0.06646, p-adjusted = 0.6646
#flicker stim vs flicker control: McNemar's chi-squared = 12, df = 1, p-value = 0.000532, p-adjusted = 0.00532

###########################################################################
###########################################################################
###                                                                     ###
###               PLOTTING RESULTS FROM BOTH EXPERIMENTS:               ###
###                                                                     ###
###########################################################################
###########################################################################

banner("Plotting results from both experiments:", emph = TRUE)

##Required dependencies/packages. 
#For plots
library(Hmisc) #For calculating confidence intervals
library(ggplot2)

rm(list=ls(all=TRUE))

# Set up graph theme
Graph.theme<-theme_bw()+
  theme(axis.text=element_text(size=12,colour="black"),
        axis.title=element_text(size=15),
        axis.title.y=element_text(vjust=1),axis.title.x=element_text(vjust=-1),
        #plot.margin = unit(c(1,1,1,1), "cm"),#This ajusts the margines of the plot so that bits don't get cut off.
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.key.size=unit(1, "cm"),
        legend.background=element_rect(),
        legend.key = element_blank(), #this gets ride of the annoying box around each of
        panel.border = element_rect(colour = "black", linewidth= 1),
        #panel.border=element_blank(), axis.line=element_line(colour = "black", size= 1),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.spacing = unit(0, "lines")) + #set the size of the space between subplots when using facet_grid()
  theme(strip.background = element_blank(),strip.text.x = element_blank())#This removes the labels that appear above each facet

##================================================================
##              Raw plots for intensity experiment               =
##================================================================

boxup("Raw plots for intensity experiment", bandChar = "=")

#Load data
I_data.with.NA<-read.csv("Data from intensity experiment.csv",header=T)
I_data <-na.omit(I_data.with.NA)

#Prepare data
XNP_I = matrix(nrow=5, ncol=4)
j <- 1
while (j<6) {
  mm2<-subset(I_data, StimNumber==j)
  Replication_check <-as.data.frame(table((duplicated(mm2$CrabID, incomparables = FALSE))))
  XNP_I[j,3]<- Replication_check[2,2] #number of pseudoreplicates (NA=good)
  XNP_I[j,1]<- sum(mm2$BinaryRes) #total number of responses 
  XNP_I[j,2]<- nrow(mm2) #n 
  XNP_I[j,4]<- mm2[1,8] #Stimulus
  j<-j+1
}
print(XNP_I) ###CHECK THAT THE THIRD COL IS ALL NA

#Calculate response probably and Wilson score intervals. 
I_CI = matrix(nrow=5, ncol=6)
I_CI[,1:5] <- binconf(XNP_I[,1],XNP_I[,2], alpha=0.05,
                       method=c("wilson"),
                       include.x=TRUE, include.n=TRUE, return.df=FALSE)
I_CI[,6] <- XNP_I[,4]
colnames(I_CI) <- c("X","N","PointEst","Lower","Upper","Stimulus")
I_G <- data.frame(I_CI)
I_G["Stimuli"] <- c("0Control0","1FO1","2SO12","3SO23","4FC4")
I_G

Fig_I <-ggplot(I_G, aes(x=Stimuli,y=PointEst))  + geom_point(color="black",size=3)+
  geom_errorbar(aes(ymin=Lower, ymax=Upper),width=.15) + ylab("Response probability
                                                              ") + xlab("Stimulus
                                                                        ")+ 
  scale_y_continuous(breaks= seq(0,1,.2), limits=c(0, 1),expand = c(0, 0)) +
  Graph.theme+ background_grid(major = "y", minor = "none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4))

Fig_I
#N= 60-59

##===============================================================
##            Raw plots for polarization experiment             =
##===============================================================

boxup("Raw plots for polarization experiment", bandChar = "=")

P_data.with.NA<-read.csv("Data from polarization experiment.csv",header=T)
P_data <-na.omit(P_data.with.NA)

XNP_P = matrix(nrow=5, ncol=4)
j <- 1
while (j<6) {
  mm3<-subset(P_data, StimNumber==j)
  Replication_check <-as.data.frame(table((duplicated(mm3$CrabID, incomparables = FALSE))))
  XNP_P[j,3]<- Replication_check[2,2] #number of pseudoreplicates (NA=good)
  XNP_P[j,1]<- sum(mm3$BinaryRes) #total number of responses 
  XNP_P[j,2]<- nrow(mm3) #n 
  XNP_P[j,4]<- mm3[1,8] #Stimulus
  j<-j+1
}
print(XNP_P) ###CHECK THAT THE THIRD COL IS ALL NAs.

P_CI = matrix(nrow=5, ncol=6)
P_CI[,1:5] <- binconf(XNP_P[,1],XNP_P[,2], alpha=0.05,
                       method=c("wilson"),
                       include.x=TRUE, include.n=TRUE, return.df=FALSE)
P_CI[,6] <- XNP_P[,4]
colnames(P_CI) <- c("X","N","PointEst","Lower","Upper","Stimulus")
P_G <- data.frame(P_CI)
P_G["Stimuli"] <- c("0Control0","1FO1","2SO12","3SO23","4FC4")
P_G

Fig_P <-ggplot(P_G, aes(x=Stimuli,y=PointEst))  + geom_point(color="black",size=3)+
  geom_errorbar(aes(ymin=Lower, ymax=Upper),width=.15) + ylab("Response probability
                                                              ") + xlab("Stimulus
                                                                        ")+ 
  scale_y_continuous(breaks= seq(0,1,.2), limits=c(0, 1),expand = c(0, 0)) +
  Graph.theme+ background_grid(major = "y", minor = "none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4))

Fig_P
#N=58-56

#############################################################################################################
### For the paper these plots were saved as a .svg file to enable them to be 
### corrected formatted and labeled outside R. ### 
#############################################################################################################

