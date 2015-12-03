rm(list=ls())
library(reshape2)
library(ggplot2)

load("~/Desktop/Marmot ews /Data/inds_notcorrected.Rda")
#source('~/Desktop/Marmot ews /R code/Generalised code - composite ews.R', chdir = TRUE)
source('~/Desktop/Marmot ews /R code/Generalised code - GAM composite ews .R', chdir = TRUE)
head(inds)

########################################################################################
##pop sizes
abun<-melt(tapply(inds$id, list(inds$year), length));names(abun)<-c("Year", "Abundance")
abun<-subset(abun, Year>=1976)

##plot the abundances
ggplot(data=abun, aes(x=Year,  y=Abundance))+geom_line()+theme_classic()

##ews for abundance only data
ews<-composite_ews(abun, c("cv", "acf", "ar1", "dr", "rr"), threshold=2, T, knots=-1)

ggplot(data=abun, aes(x=Year,  y=Abundance))+geom_line()+theme_classic()+geom_point(data=ews, aes(x=time, y= GAM.singals+100))

########################################################################################
##calcualte mean body size
size<-melt(tapply(inds$mass, list(inds$year), mean, na.rm=T));names(size)<-c("Year", "Mass")
size<-subset(size, Year>=1976)
##plot the sizes
ggplot(data=size, aes(x=Year,  y=Mass))+geom_line()+theme_classic()
########################################################################################


