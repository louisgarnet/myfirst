####################################################################################################
##                                                                                                ##
##  This script runs a generalized linear mixed model (GLMM) to test differences between          ##
##     soil CO2 concentration s[CO2] in soil affected or not by leaf cutter ant presence. It also ##
##     creates tables with correlation between s[CO2] and depth for different precipitation       ##
##     patterns. It runs a GLMM to test for differences between soil CO2 emissions and soil       ##
##     water content. Leaf cutter ant species studied is Atta cephalotes.                         ##
##                                                                                                ##
##            Author: Angel S. Fernandez Bou - afernandezbou@ucmerced.edu                         ##
##                                                                                                ##
####################################################################################################

####################################################################################################
#-----     SCRIPT HISTORY LOG     -----------------------------------------------------------------#
#   Published on March 03, 2018                                                                    #
#     Created by Angel S. Fernandez Bou (afernandezbou@ucmerced.edu) and Diego Dierick             #
#                                                                                                  #
#                                                                                                  #
####################################################################################################

####################################################################################################
#-----     IMPORTANT NOTES     --------------------------------------------------------------------#
#  - Open RStudio or equivalent using the file containing this script; this way it will detect     #
#       the folder to define the directory location and that will make running this code automatic #
#  - This code must be run from a folder containing 2 other folders                                #
#      "Data": with the csv file with the data set                                                 #
#      "Results": where the results are placed after running the code                              #
####################################################################################################

####################################################################################################
#-----     LIBRARIES USED IN THIS SCRIPT     ------------------------------------------------------#
#   If you do not have these libraries, uncomment removing the "#" symbols first time you run it   #
#                                                                                                  #
#install.packages("MASS") ## For glmmPQL                                                           #
#install.packages("car")    ## For qqp                                                             #
#install.packages(fitdistrplus) ## For plotting models                                             #
#install.packages("GGally")                                                                        #
#install.packages("zoo")                                                                           #
####################################################################################################

####################################################################################################
#-----     SUMMARY     ----------------------------------------------------------------------------#
#                                                                                                  #
# PART 1: soil CO2 concentration. Comparison of Control vs Abandoned nest vs Nest considering 3    #
#           different period as a proxy of the possible legacy effect of the abandonment.          #
#           Compares seasonality in different patterns of Dry vs Wet (accumulated precipitation)   #
# PART 2: Soil CO2 efflux and moisture analysis                                                    #
####################################################################################################

####################################################################################################
##########----------          RUN THE WHOLE SCRIPT HERE                         ----------##########
# {   # This bracket closes at the end of the script                                                 #
####################################################################################################



####################################################################################################
#####------------------------------------------------------------------------------------------#####
##### PART 1 --- Control vs Abandoned Nest vs Nest --- START                                   #####
#####------------------------------------------------------------------------------------------#####
##### Comparison of Control vs Abandoned Nest vs Nest with actual data                         #####
##### Since some nests were abandoned while the experiment was performed, we did in addition   #####
##### the same analyses explained before after rearranging the data into 3 plot types:         ##### 
##### control, nest and abandoned nest.                                                        #####
##### We did it both for actual data and transformed data. Analysis (3) was performed only for ##### 
##### primary canopy, since there was not enough data to perform it for the secondary canopy,  ##### 
##### and analysis (5) was not performed due to lack of data (too few repetitions/locations).  #####
#####------------------------------------------------------------------------------------------#####
####################################################################################################

####################################################################################################
#-----     DESCRIPTION OF THE DATASET     ---------------------------------------------------------#
##                                                                                                ##
## 	TagNr: ID of each gas well (3 per plot, one at each depth)                                    ##
## 	Plot_ID_original: code given to each plot (E.g., ALPN1)                                       ##
## 	SiteID: ALP, REP, RES                                                                         ##
## 	plot_type: Control and Nest (disregarding abandoned nests)                                    ##
## 	Depth: 20 cm, 60 cm, 100 cm                                                                   ##
## 	Date: day of the data collection. Base for the campaign                                       ##
## 	sCO2_raw: soil CO2 concentration including #N/A                                               ## 
## 	sCO2: soil CO2 concentration where NA replace #N/A                                            ##
## 	ID: plot ID with depth (not used in this script)                                              ##
## 	PlotID: plot ID considering abandoned nests                                                   ##
## 	Soil: AL is alluvial and RE is residual                                                       ##
## 	Canopy: P is primary and S is secondary                                                       ##
## 	Campaign: sequential number of each time data were collected (once a month, normally)         ##
##        Included to capture the variability of time in the GLMM                                 ##
## 	Type: Control (C), Abandoned Nest (L), Nest (N)                                               ##
## 	Type2: Control (C), Abandoned Nest (A), Nest (N)                                              ##
##                                                                                                ##
####################################################################################################

# { # RUN PART 1
  { # RUN Start PART 1
rm(list=ls(all=TRUE))                   ## This deletes all variables

## Define directory location                          
CurrentDirectory = getwd()              ## R script folder
ReadDirectory = file.path(getwd(),"")   ## data folder
setwd(ReadDirectory)                    ## This sets "ReadDirectory" as the directory to retrieve files

## Reading the file with all the gas well data
gaswells=read.csv('Data/data.dat',header=TRUE,stringsAsFactors = FALSE)

rm(list=setdiff(ls(),c("gaswells","CurrentDirectory","ReadDirectory"))) # this removes all the variables but these ones listed

## Removing rows with NA
GW <- na.omit(gaswells)     # for initial analysis with sCO2

GW$Date <- as.POSIXct(strptime(GW$Date, format="%Y-%m-%d", tz="UTC"))
### (UTC to avoid daylight saving issues, but values are local times)

MeteoData = read.csv('Data/meteodata.dat',header=TRUE,stringsAsFactors = FALSE)# fill rows of unequal length
# colnames(MeteoData)[1] <- "Date"
MeteoData$Date = as.POSIXct(MeteoData$Date, tz="UTC")                          # make date string to POSIXct
MeteoData$Date[ which(is.na(MeteoData$Rain_Tot)) ]                             # a few days have NA precipitation data
MeteoData$Rain_Tot[ which(is.na(MeteoData$Rain_Tot)) ] = 0                     # we will set it to 0 mm
    # For the MeteoData data frame, Rain_Tot and Date are the only variables used

Historic <- read.csv('Data/historic.dat',header=TRUE,stringsAsFactors = FALSE) # Data from 1985 to 2015
# In this meteorological data set, NA data are 0.
# using the average of the day before and after the NA (when possible), hardly changes the value (~0.4%)
colnames(Historic) <- list('Date','Rain')
Historic$Date <- as.POSIXct(strptime(Historic$Date, format="%Y-%m-%d", tz="UTC"))

## Declaring libraries
library("MASS") ## For glmmPQL
library("car")    ## For qqp, to decide the best density distribution
library("fitdistrplus") ## For plotting models
library("GGally") ## for ggpairs and plot the correlations between variables (used to decide the seasonality patterns)
library("zoo")  
} # END Start PART 1


##################################################
#####----------------------------------------#####
#####---- This script generates the GLMM ----#####
#####---- results for differences among  ----#####
#####---- plot types                     ----#####
#####----------------------------------------#####
##################################################
{
MeteoData = MeteoData[which(MeteoData$Date > as.POSIXct("2014-03-01 00:00", tz="UTC")),]        # gaswell data only startes March 2015
    ## MeteoData is used to calculate the moving average of the precipitation from 1 to 360 prior to each sampling campaign
TagNr = unique(GW$TagNr)
result = NULL
resultline <- NULL

#####----- This loop calculates the significance of the differences among plot types using a GLMM
for (window in 1:360)
{
  runAvgRain = rollapply(MeteoData$Rain_Tot, window, mean, align="right", fill=("extend"))         # average rain for window length
  thresholdmean <- mean(Historic[(Historic$Date >= "1986-01-01"),]$Rain)
  subSet = GW[-c(15:41,50:59)]
  avgRain = NULL 
  for (date in 1:nrow(subSet) )                                                                   # for each line dataset
  {
    avgRain = c(avgRain,runAvgRain[which(MeteoData$Date==subSet$Date[date])])                     # look up avg. rain for that date
  }
  subSet$avgRain <- avgRain    # this creates the variable Season and Season2, which define if a period is wet or dry when compared with the historic daily mean precipitation
  subSet$Season[subSet$avgRain >= thresholdmean]="Wet"
  subSet$Season[subSet$avgRain < thresholdmean]="Dry"
  subSet$Season2[subSet$avgRain >= thresholdmean]="Wet"
  subSet$Season2[subSet$avgRain < thresholdmean]="zDry"
  subSet<-na.omit(subSet)
  subSet$Depth <- as.factor(subSet$Depth)
  
  summD <- summary(glmmPQL(sCO2 ~ Type + Depth + Season #fix factors 
                  +Soil + Canopy #Soil and Canopy not significant, but they are kept since they are low-level factors that bring information
                  # +plot*depth+ # this interaction is not significant and seems to be affected by the seasonality # not significant 91.7% of the times # 70% of the times p-value > 0.3
                  +Type*Season + Depth*Season #interactions
                  ,c(~1|TagNr), # random factors
                  family = gaussian(link="log"), # density distribution
                  data=subSet, # data set ALL
                  verbose=FALSE))
  
  summW <- summary(glmmPQL(sCO2 ~Type+Depth+Season2 #fix factors; Season2 is like Season but with different sorting
                  +Soil + Canopy #Soil and Canopy not significant, can be removed to simplify
                  # +plot*depth+ # this interaction is not significant and seems to be affected by the seasonality # not significant 91.7% of the times
                  +Type*Season2 + Depth*Season2 #interactions
                  ,c(~1|TagNr), # random factors
                  family = gaussian(link="log"), # density distribution
                  data=subSet, # data set ALL
                  verbose=FALSE))
  
  summ2D <- summary(glmmPQL(sCO2 ~Type2+Depth+Season #fix factors; Type2 is like Type but with different sorting
                  +Soil +Canopy #Soil and Canopy not significant, can be removed to simplify
                  # +plot*depth+ # this interaction is not significant and seems to be affected by the seasonality # not significant 91.7% of the times
                  +Type2*Season + Depth*Season #interactions
                  ,c(~1|TagNr), # random factors
                  family = gaussian(link="log"), # density distribution
                  data=subSet, # data set ALL
                  verbose=FALSE))

  summ2W <- summary(glmmPQL(sCO2 ~Type2+Depth+Season2 #fix factors 
                 +Soil +Canopy #Soil and Canopy not significant, can be removed to simplify
                 # +plot*depth+ # this interaction is not significant and seems to be affected by the seasonality # not significant 91.7% of the times
                 +Type2*Season2 + Depth*Season2 #interactions
                 ,c(~1|TagNr), # random factors
                 family = gaussian(link="log"), # density distribution
                 data=subSet, # data set ALL
                 verbose=FALSE))
  
    CvsN.D <- round(summD$tTable[3,5],4)
    CvsA.D <- round(summD$tTable[2,5],4)
    NvsA.D <- round(summ2D$tTable[3,5],4)
    CvsN.W <- round(summW$tTable[3,5],4)
    CvsA.W <- round(summW$tTable[2,5],4)
    NvsA.W <- round(summ2W$tTable[3,5],4)
    resultline <- c(resultline,CvsN.D,CvsA.D,NvsA.D,CvsN.W,CvsA.W,NvsA.W)
  
  result <- rbind(result,resultline)
  resultline <- NULL#matrix(,nrow=0,ncol=0)
}
colnames(result) <- c('CvsN Dry','CvsA Dry','AvsN Dry','CvsN Wet','CvsA Wet','AvsN Wet')
rownames(result) <- seq(1,360)
write.csv(result,'Results/p.csv')
      # p.csv is saved on the folder Results (previously created). p presents the p-values of the differences among plot types for different dry/wet periods
}
##################################################
#####-- End of p.csv                       --#####
##################################################


##################################################
#####----------------------------------------#####
#####----- Table with sCO2 monthly means-----#####
#####----------------------------------------#####
##################################################
{
MeansRow <- NULL
MeansMatrix <- NULL
for (camp in 1:35)
 {
  MeansRow <- round(c(mean(GW[(GW$Campaign==camp),]$sCO2),mean(GW[(GW$Campaign==camp & GW$Type=="C"),]$sCO2),mean(GW[(GW$Campaign==camp & GW$Type=="L"),]$sCO2),mean(GW[(GW$Campaign==camp & GW$Type=="N"),]$sCO2)),0)
  MeansMatrix <- rbind(MeansMatrix,MeansRow)
  MeansRow <- NULL 
 }
colnames(MeansMatrix) <- list("All","Control", "Abandoned","Nest")
rownames(MeansMatrix) <- unique(GW$Campaign)
write.csv(MeansMatrix,'Results/sCO2means.csv')
}
##################################################
#####-- End of sCO2means.csv               --#####
##################################################


##################################################
#####----------------------------------------#####
#####----- Table with mean, SD, median, -----#####
#####----- estima and p-value           -----#####
#####----- for P90d Dry and Wet         -----#####
#####----------------------------------------#####
##################################################
{
result = NULL
resultline = NULL
for(seasDW in unique(GW$P90dseason))
{
  for (profundidad in c(20,60,100))
  {
    for (PlotType in unique(GW$Type)) 
    {
      resultline <- c(resultline,PlotType,profundidad,seasDW,
                      round(c(mean(GW[(GW$Type==PlotType & GW$Depth==profundidad & GW$P90dseason==seasDW),]$sCO2),
                              sd(GW[(GW$Type==PlotType & GW$Depth==profundidad & GW$P90dseason==seasDW),]$sCO2),
                              median(GW[(GW$Type==PlotType & GW$Depth==profundidad & GW$P90dseason==seasDW),]$sCO2)
                      ),0))
    }
    result <- rbind(result,resultline)
    resultline=NULL
  }
}
colnames(result) <- c('Type','Depth', 'Season', 'mean','sd','median','Type','Depth', 'Season', 'mean','sd','median','Type','Depth', 'Season', 'mean','sd','median')

# for the estimates 
estimates.P90d = NULL
estimates.P90d <- rbind(estimates.P90d,round(c(exp(Nd20W90[1,1]),exp(Cd20W90[1,1]),exp(Ad20W90[1,1])),0))
estimates.P90d <- rbind(estimates.P90d,round(c(exp(Nd60W90[1,1]),exp(Cd60W90[1,1]),exp(Ad60W90[1,1])),0))
estimates.P90d <- rbind(estimates.P90d,round(c(exp(Nd100W90[1,1]),exp(Cd100W90[1,1]),exp(Ad100W90[1,1])),0))
estimates.P90d <- rbind(estimates.P90d,round(c(exp(Nd20D90[1,1]),exp(Cd20D90[1,1]),exp(Ad20D90[1,1])),0))
estimates.P90d <- rbind(estimates.P90d,round(c(exp(Nd60D90[1,1]),exp(Cd60D90[1,1]),exp(Ad60D90[1,1])),0))
estimates.P90d <- rbind(estimates.P90d,round(c(exp(Nd100D90[1,1]),exp(Cd100D90[1,1]),exp(Ad100D90[1,1])),0))
colnames(estimates.P90d) <- c("Estimate N","C","A")
row.names(estimates.P90d) <- list('Wet 20 cm','60 cm','100 cm','Dry 20 cm','60 cm','100 cm')

table.results <- cbind(result,estimates.P90d)
write.csv(table.results,'Results/mean sd median estimate.csv')
}
##################################################
#####-- End of mean sd median estimate.csv --#####
##################################################


##################################################
#####----------------------------------------#####
#####-----        CORRELATIONS          -----#####
#####----------------------------------------#####
#####--this loop calculates the correlation--#####
#####--between the daily precipitation     --#####
#####--moving average from 1 to 360 days   --#####
#####--and the sCO2, and the same by depth --#####
#####----------------------------------------#####
##################################################
{
MeteoData = MeteoData[which(MeteoData$Date > as.POSIXct("2014-03-01 00:00", tz="UTC")),]        ### gaswell data only startes March 2015

gaswell = unique(GW$TagNr)
PlotType <- unique(GW$Type)
result = NULL                                             ### prepare variable to hold results
resultline <- NULL

#####----- This loop calculates the correlation for control, abandoned and nest (it does not split by depth)
for (window in 1:360)
{
runAvgRain = rollapply(MeteoData$Rain_Tot, window, mean, align="right", fill=("extend"))           ### average rain for window length
for (PlotType in unique(GW$Type))
{
  subSet = GW[which(GW$Type == PlotType),]
  avgRain = NULL                                                                                   ### init array to hold avg. rain on matched dates
  for (date in 1:nrow(subSet) )                                                                    ### for each line dataset
  {
    avgRain = c(avgRain,runAvgRain[which(MeteoData$Date==subSet$Date[date])])                   ### look up avg. rain for that date
  }
  subSet$avgRain <- avgRain
  subSet<-na.omit(subSet)
  correl <- cor(subSet$sCO2,subSet$avgRain)
  cor.p <- cor.test(subSet$sCO2,subSet$avgRain)$p.value
  
  resultline <- c(resultline,correl,round(cor.p,4))
}
result <- rbind(result,resultline)
resultline <- NULL
}
colnames(result) <- c('Abandoned','A p','Control','C p','Nest','N p')
rownames(result) <- seq(1,360)
write.csv(result,'Results/correl.csv')
rm(subSet,result)            ### cleaning





#####----- This loop calculates the correlation by depth 
Deep <- unique(GW$Depth)
result = NULL                                             ### prepare variable to hold results
resultline <- NULL

#####----- This loop calculates the correlation for control, abandoned and nest (it does not split by depth)
for (Deep in unique(GW$Depth)) 
{
for (PlotType in unique(GW$Type))
{
  subSet = GW[which(GW$Type == PlotType & GW$Depth == Deep),]
  
  for (window in 1:360)
  {
    runAvgRain = rollapply(MeteoData$Rain_Tot, window, mean, align="right", fill=("extend"))           ### average rain for window length
    
    avgRain = NULL                                                                                   ### init array to hold avg. rain on matched dates
    
    for (date in 1:nrow(subSet) )                                                                    ### for each line dataset
    {
      avgRain = c(avgRain,runAvgRain[which(MeteoData$Date==subSet$Date[date])])                   ### look up avg. rain for that date
    }
    subSet$avgRain <- avgRain
    subSet<-na.omit(subSet)
    correl <- cor(subSet$sCO2,subSet$avgRain)
    
    resultline <- c(resultline,correl)
  }
  resultline <- c(PlotType,Deep,resultline)
  result <- rbind(result,resultline)
  resultline <- NULL
}
}
colnames(result) <- c('Type','Depth',1:360)
write.csv(result,'Results/correl by depth.csv',row.names = FALSE)
rm(subSet,result)            ### cleaning
}
##################################################
#####-----     END TABLE correlations   -----#####
##################################################
  
  
# } # END PART 1
####################################################################################################
#####------------------------------------------------------------------------------------------#####
##### END --- PART 1 --- Control vs Abandoned Nest vs Nest                                     #####
#####------------------------------------------------------------------------------------------#####
####################################################################################################




####################################################################################################
#####------------------------------------------------------------------------------------------#####
##### PART 2 --- Soil CO2 efflux and moisture analysis                                         #####
#####------------------------------------------------------------------------------------------#####
#####------------------------------------------------------------------------------------------#####
####################################################################################################
{

####################################################################################################
#-----     DESCRIPTION OF THE DATASET     ---------------------------------------------------------#
##                                                                                                ##
## 	Date: day of the data collection.                                                             ##
## 	ID: code given to each plot                                                                   ##
## 	VWC: volumetric water content                                                                 ##
## 	Flux: efflux measured with our low-cost CO2 detection chamber                                 ##
## 	Soil: AL is alluvial and RE is residual                                                       ##
## 	CANOPY: P is primary and M is modified (similar to secondary)                                 ##
## 	Local: identification number given to each pair (C,N or C,A) or triplet (C,N,A) of plots      ##
## 	CAMPAIGN: sequential number of each time data were collected                                  ##
##        Included to capture the variability of time in the GLMM                                 ##
## 	Type: Control (C), Abandoned Nest (L), Nest (N)                                               ##
## 	Type2: Control (C), Abandoned Nest (A), Nest (N)                                              ##
##                                                                                                ##
####################################################################################################


rm(list=ls(all=TRUE))                                                 ## This deletes all variables

## Define directory location                                                                     ##
CurrentDirectory = getwd()                                            ## R script folder
ReadDirectory = file.path(getwd(),"")                 ## data folder
setwd(ReadDirectory)                                                  ## This sets "ReadDirectory" as the directory to retrieve files

# Reading the file
efflux.raw=read.csv('Data/efflux.dat',header=TRUE, stringsAsFactors = FALSE)

# # Remove rows with NA
# efflux <- na.omit(efflux.raw)

# # Format Date
# efflux$Date <- strptime(as.character(efflux$Date), "%m/%d/%y")

## Declaring libraries
library("MASS") ## For glmmPQL
library("car")    ## For qqp
library("fitdistrplus") ## For plotting models
library("GGally")
library("nortest") ## to test normality ad.test(x)
}

##################################################
#####----------------------------------------#####
#####-----     TABLE FLUX                ----#####
#####----------------------------------------#####
#####----------------------------------------#####
##################################################
{
efflux <- efflux.raw[!is.na(efflux.raw$Flux),]

# Means and standard deviation
mean(efflux$Flux)
sd(efflux$Flux)
mean(efflux[(efflux$Type=="C"),]$Flux)
sd(efflux[(efflux$Type=="C"),]$Flux)
mean(efflux[(efflux$Type=="L"),]$Flux)
sd(efflux[(efflux$Type=="L"),]$Flux)
mean(efflux[(efflux$Type=="N"),]$Flux)
sd(efflux[(efflux$Type=="N"),]$Flux)

### model 
CvsLvsN<- glmmPQL(Flux~Type+CANOPY,~1|Local/CAMPAIGN,family = gaussian(link = "identity"), efflux)
summary(CvsLvsN)

LvsN<- glmmPQL(Flux~Type2+CANOPY,~1|Local/CAMPAIGN,family = gaussian(link = "identity"), efflux)
summary(LvsN)

### plot
par(mfrow=c(1,1),oma=c(0,0.2,0,0))
boxplot(Flux~Type,efflux,col=c("green","blue","red"),names=c("Control","Abandoned","Nest"),ylab="Flux (µmol / m² s)",las=1,cex.lab=1.2)#,xlab="Type")
title(main=expression("Soil CO"[2]*" emissions"))
title(main=expression("                                           a"),cex.main = 1.5)
text(1,12,"5.0 ± 2.0")
text(2,12,"5.0 ± 2.1")
text(3,12,"5.3 ± 2.2")
}

##################################################
#####----------------------------------------#####
#####-----     TABLE SOIL MOISTURE       ----#####
#####----------------------------------------#####
#####----------------------------------------#####
##################################################
{
effluxVWC <- efflux.raw[!is.na(efflux.raw$VWC),]

# Means and standard deviation
mean(effluxVWC$VWC)
sd(effluxVWC$VWC)
mean(effluxVWC[(effluxVWC$Type=="C"),]$VWC)
sd(effluxVWC[(effluxVWC$Type=="C"),]$VWC)
mean(effluxVWC[(effluxVWC$Type=="L"),]$VWC)
sd(effluxVWC[(effluxVWC$Type=="L"),]$VWC)
mean(effluxVWC[(effluxVWC$Type=="N"),]$VWC)
sd(effluxVWC[(effluxVWC$Type=="N"),]$VWC)

#### VWC efflux's data ####
### model 
VWCCvsLvsN<- glmmPQL(VWC~Type+CANOPY,~1|Local/CAMPAIGN,family = gaussian(link = "identity"), effluxVWC)
summary(VWCCvsLvsN)

VWCLvsN<- glmmPQL(VWC~Type2+CANOPY,~1|Local/CAMPAIGN,family = gaussian(link = "identity"), effluxVWC)
summary(VWCLvsN)

## plot
par(mfrow=c(1,1),oma=c(0,0.2,0,0))
boxplot(VWC~Type,efflux,col=c("green","blue","red"),names=c("Control","Abandoned","Nest"),ylab="VWC (m³ / m³)",ylim=c(0.35,0.7),las=1,cex.lab=1.2)#,xlab="Plot type")
title(main=expression("Soil Moisture"))
title(main=expression("                                          b"),cex.main = 1.5)
# text(1,0.68,"0.519 ±0.045
# A")
# text(2,0.68,"0.494 ±0.047
# B")
# text(3,0.68,"0.503 ±0.053
# B")
text(1,0.68,bquote(.(round(mean(effluxVWC[(effluxVWC$Type=="C"),]$VWC),2)) ~"±"
                   ~ .(round(sd(effluxVWC[(effluxVWC$Type=="C"),]$VWC),2))))
text(1,0.65,"A")
text(2,0.68,bquote(.(round(mean(effluxVWC[(effluxVWC$Type=="L"),]$VWC),2)) ~"±"
                   ~ .(round(sd(effluxVWC[(effluxVWC$Type=="L"),]$VWC),2))))
text(2,0.65,"B")
text(3,0.68,bquote(.(round(mean(effluxVWC[(effluxVWC$Type=="N"),]$VWC),2)) ~"±"
                   ~ .(round(sd(effluxVWC[(effluxVWC$Type=="N"),]$VWC),2))))
text(3,0.65,"B")

median(effluxVWC[(effluxVWC$Type=="C"),]$VWC)
median(effluxVWC[(effluxVWC$Type=="L"),]$VWC)
median(effluxVWC[(effluxVWC$Type=="N"),]$VWC)

quantile(effluxVWC[(effluxVWC$Type=="C"),]$VWC,0.71) # VWC=0.55 for that quantile
quantile(effluxVWC[(effluxVWC$Type=="L"),]$VWC,0.907)
quantile(effluxVWC[(effluxVWC$Type=="N"),]$VWC,0.8216)
## Quantiles show that 29%, 9% and 18% of the time Control, abandoned and nest soils were as wet or wetter than 0.55
}
####################################################################################################
#####------------------------------------------------------------------------------------------#####
##### END --- PART 2 --- Soil CO2 efflux and moisture analysis                                 #####
#####------------------------------------------------------------------------------------------#####
####################################################################################################

## } # END OF SCRIPT
####################################################################################################
#####     END OF SCRIPT                                                                        #####
####################################################################################################


