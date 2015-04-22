#=============================================================================================================#
# Script created by Meghan Duffy, duffymeg@umich.edu
# Script created in version R 3.1.2 
# This script is for analyzing data related to the Duffy et al. paper re-describing
# the Daphnia brood parasite Blastulidium paedophthorum (Bp)
#=============================================================================================================#

# Set working directory
#use line below if working on pc
setwd("C:/Users/duffymeg/Box Sync/Parasites/BroodParasite/CodeforDuffyetalBPpaper")
#use line below if working on mac
# setwd("~/Box Sync/Parasites/BroodParasite/CodeforDuffyetalBPpaper")

# load all data
# North Fall 2013 life table
brooddataNorth2013 <- read.table("NorthBPLifeTable2013.txt", header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
# Now moving on to Fall 2014 data 
# North Lake:
brooddataNorth2014 <- read.csv("NorthBPLifeTable2014.csv", header=TRUE, na.strings="?", dec=".", strip.white=TRUE)
brooddataNorth2014dentifera <- subset(brooddataNorth2014,Treatment=="Infected Dent" | Treatment=="Uninfected Dent")
brooddataNorth2014retrocurva <- subset(brooddataNorth2014,Treatment=="Infected Retro" | Treatment=="Uninfected Retro")
# Cedar Lake:
brooddataCedar2014 <- read.csv("CedarBPLifeTable2014.csv", header=TRUE, na.strings="?", dec=".", strip.white=TRUE)
brooddataCedar2014dentifera <- subset(brooddataCedar2014,Treatment=="Infected Dent" | Treatment=="Uninfected Dent")
brooddataCedar2014retrocurva <- subset(brooddataCedar2014,Treatment=="Infected Retro" | Treatment=="Uninfected Retro")

# Survival analyses
library(survival)
#Note "Censored" column is whether or not the animal died before experimental takedown 1= yes, it died during the experiment, 0 = no, it survived until the end

# North 2013 data:
model1 <- survfit(Surv(brooddataNorth2013$Day.Death, brooddataNorth2013$Censored) ~ brooddataNorth2013$Treatment)
plot(model1, lty=c(5,1),lwd=3)

tapply(brooddataNorth2013$Day.Death,brooddataNorth2013$Treatment,mean,na.rm=TRUE)

# Cox's proportional hazard model 
coxfit1 <- coxph(Surv(brooddataNorth2013$Day.Death, brooddataNorth2013$Censored) ~ brooddataNorth2013$Treatment)  
summary(coxfit1)   

# North 2014 dentifera data:
model2 <- survfit(Surv(brooddataNorth2014dentifera$Day.Death, brooddataNorth2014dentifera$Censored) ~ brooddataNorth2014dentifera$Treatment)
plot(model2, lty=c(5,1),lwd=3)

tapply(brooddataNorth2014dentifera$Day.Death,brooddataNorth2014dentifera$Treatment,mean,na.rm=TRUE)

coxfit2 <- coxph(Surv(brooddataNorth2014dentifera$Day.Death, brooddataNorth2014dentifera$Censored) ~ brooddataNorth2014dentifera$Treatment)  
summary(coxfit2)   

# North 2014 retrocurva data:
model3 <- survfit(Surv(brooddataNorth2014retrocurva$Day.Death, brooddataNorth2014retrocurva$Censored) ~ brooddataNorth2014retrocurva$Treatment)
plot(model3, lty=c(5,1),lwd=3)

tapply(brooddataNorth2014retrocurva$Day.Death,brooddataNorth2014retrocurva$Treatment,mean,na.rm=TRUE)

coxfit3 <- coxph(Surv(brooddataNorth2014retrocurva$Day.Death, brooddataNorth2014retrocurva$Censored) ~ brooddataNorth2014retrocurva$Treatment)  
summary(coxfit3)   

# Cedar 2014 dentifera data:
model4 <- survfit(Surv(brooddataCedar2014dentifera$Day.Death, brooddataCedar2014dentifera$Censored) ~ brooddataCedar2014dentifera$Treatment)
plot(model4, lty=c(1,5),lwd=3)

tapply(brooddataCedar2014dentifera$Day.Death,brooddataCedar2014dentifera$Treatment,mean,na.rm=TRUE)

coxfit4 <- coxph(Surv(brooddataCedar2014dentifera$Day.Death, brooddataCedar2014dentifera$Censored) ~ brooddataCedar2014dentifera$Treatment)  
summary(coxfit4)   

# Cedar 2014 retrocurva data:
model5 <- survfit(Surv(brooddataCedar2014retrocurva$Day.Death, brooddataCedar2014retrocurva$Censored) ~ brooddataCedar2014retrocurva$Treatment)
plot(model5, lty=c(5,1),lwd=3)

tapply(brooddataCedar2014retrocurva$Day.Death,brooddataCedar2014retrocurva$Treatment,mean,na.rm=TRUE)

coxfit5 <- coxph(Surv(brooddataCedar2014retrocurva$Day.Death, brooddataCedar2014retrocurva$Censored) ~ brooddataCedar2014retrocurva$Treatment)  
summary(coxfit5)  


# One big figure with all the survival plots:
par(mfrow=c(3,2))
plot(model1, lty=c(1,1),lwd=3, col=c("gray","black"))
plot(model1, lty=c(1,1),lwd=3, col=c("gray","black")) # placeholder; will delete when export
plot(model2, lty=c(1,1),lwd=3, col=c("gray","black"))
plot(model3, lty=c(1,1),lwd=3, col=c("gray","black"))
plot(model4, lty=c(1,1),lwd=3, col=c("gray","black"))
plot(model5, lty=c(1,1),lwd=3, col=c("gray","black"))

## Now combining all the data into one big dataset to plot reproduction:
brooddataNorth2013$LakeYear <- "North2013"
brooddataNorth2013$Host <- "dentifera"
brooddataNorth2013$Infection <- brooddataNorth2013$Treatment

library(stringr)
library(plyr)

list <- str_split(brooddataNorth2014$Treatment, " ", n = 2)
df <- ldply(list)
colnames(df) <- c("Infection", "Host")

brooddataNorth2014$Infection <- df$Infection
brooddataNorth2014$Host <- df$Host

list <- str_split(brooddataCedar2014$Treatment, " ", n = 2)
df <- ldply(list)
colnames(df) <- c("Infection", "Host")

brooddataCedar2014$Infection <- df$Infection
brooddataCedar2014$Host <- df$Host

brooddataNorth2014$LakeYear <- "North2014"
brooddataCedar2014$LakeYear <- "Cedar2014"

brooddataNorth2013$Other.Infect <- NULL

allbrooddata <- rbind(brooddataNorth2013,brooddataNorth2014,brooddataCedar2014)

allbrooddata$Host <- str_replace(allbrooddata$Host,"dentifera","Dent")

allbrooddata$Dummy <- "Q"
allbrooddata$Dummy[1:25] <- "A"
allbrooddata$Dummy[26:58] <- "B"
allbrooddata$Dummy[59:71] <- "C"
allbrooddata$Dummy[72:85] <- "D"
allbrooddata$Dummy[86:100] <- "G"
allbrooddata$Dummy[101:115] <- "H"
allbrooddata$Dummy[116:123] <- "E"
allbrooddata$Dummy[124:133] <- "F"
allbrooddata$Dummy[134:137] <- "I"
allbrooddata$Dummy[138:147] <- "J"

library(ggplot2)
ggplot(allbrooddata, aes(x=Dummy, y=Sum.Offspring)) + 
  geom_boxplot() + theme_bw() +
  labs(x=" ", y=expression(paste("Total offspring"))) + theme(
    axis.title.y = element_text(size=18,vjust=0.9),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14)
  ) 

# now analyzing 
#For North 2013
#get mean reproduction of the two treatments
tapply(brooddataNorth2013$Sum.Offspring,brooddataNorth2013$Treatment,mean,na.rm=TRUE)

# non-parametric permutation test using the coin package:
library(coin)

# North 2013:
oneway_test(Sum.Offspring ~ Treatment, data=brooddataNorth2013,
            distribution=approximate(B=9999))

#For North 2014
#get mean reproduction of the two treatments
tapply(brooddataNorth2014$Sum.Offspring,brooddataNorth2014$Treatment,mean,na.rm=TRUE)

#test North 2014 dentifera using the coin package:
oneway_test(Sum.Offspring ~ Treatment, data=brooddataNorth2014dentifera,
            distribution=approximate(B=9999))

#North 2014 retrocurva
oneway_test(Sum.Offspring ~ Treatment, data=brooddataNorth2014retrocurva,
            distribution=approximate(B=9999))

#For Cedar 2014
#get mean reproduction of the treatments
tapply(brooddataCedar2014$Sum.Offspring,brooddataCedar2014$Treatment,mean,na.rm=TRUE)

#test Cedar 2014 dentifera using the coin package:
oneway_test(Sum.Offspring ~ Treatment, data=brooddataCedar2014dentifera,
            distribution=approximate(B=9999))

# test Cedar 2014 retrocurva:
oneway_test(Sum.Offspring ~ Treatment, data=brooddataCedar2014retrocurva,
            distribution=approximate(B=9999))


# looking at recovery of infected hosts
tapply(allbrooddata$Did.Clear,list(allbrooddata$LakeYear,allbrooddata$Host,allbrooddata$Infection),mean,na.rm=TRUE)
tapply(allbrooddata$Did.Clear,list(allbrooddata$LakeYear,allbrooddata$Host,allbrooddata$Infection),sum,na.rm=TRUE)
tapply(allbrooddata$Sum.Offspring,list(allbrooddata$LakeYear,allbrooddata$Host,allbrooddata$Infection,allbrooddata$Did.Clear),mean,na.rm=TRUE)

#######################CODE FOR FIELD DATA ANALYSES#####################################
library(reshape2)
library(scales)

parasitedata <- read.csv("2014ParasiteSurveyJustBrood.csv")

colnames(parasitedata)

#remove species-lake-date combinations where scanned fewer than 20 individuals of a given species
parasitedata <- subset(parasitedata,Total>19)

#remove ambigua and mendotae
shortdata <- subset(parasitedata,!Species%in%c("Ambigua","Mendotae"))

#format dates correctly
shortdata <- mutate(shortdata,Date=as.Date(Date,format="%d-%b-%y"))

shortdata$Lake <- str_replace(shortdata$Lake,"Little Appleton","LittleAppleton")
shortdata$Lake <- factor(shortdata$Lake)
shortdata$Lake2 <- shortdata$Lake
shortdata$Lake2 <- mapvalues(shortdata$Lake2, from = c("Appleton", "Bishop", "Bruin", "Cedar", "CrookedP", "CrookedW", "Gosling", "LittleAppleton", "Mill", "North", "Pickerel", "Sullivan", "Walsh", "Whitmore", "Woodland"), to = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"))
shortdata$Species <- factor(shortdata$Species)
shortdata$Species2 <- shortdata$Species
shortdata$Species2 <- mapvalues(shortdata$Species2, from = c("Ceriodaphnia", "Dentifera", "Dubia", "Parvula", "Pulicaria", "Retrocurva"), to = c("Cdu", "Dde", "Ddu", "Dpa", "Dpu", "Dre"))

#Looking at proportion of dates where brood was present
sum(shortdata$PropBrood!=0)/length(shortdata$PropBrood)

library(dplyr)
proppresent<-shortdata %>% group_by(Species2, Lake2) %>% summarize(prop = mean(PropBrood > 0), count = n())
enoughdates <- subset(proppresent,count>4)
enoughdates$perc <- enoughdates$prop*100
summarypresent <- enoughdates %>% group_by(Species2) %>% summarize (minprop = min(prop),maxprop = max(prop))
summarypresentbylake <- enoughdates %>% group_by(Lake2) %>% summarize (minprop = min(prop),maxprop = max(prop))

summarypresent
summarypresentbylake

#Plotting presence/absence data:
give.n <- function(x){
  return(c(y = 95, label = length(x)))
}

panela <- ggplot(enoughdates, aes(x=Species2, y=perc)) + 
  geom_boxplot() + theme_bw() + 
  stat_summary(fun.data = give.n, geom = "text") + 
  labs(y="Percent of sampling\ndates present") + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=16,vjust=-1.9),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=14)
  ) 

panelb <- ggplot(enoughdates, aes(x=Lake2, y=perc)) + 
  geom_boxplot() + theme_bw() + 
  stat_summary(fun.data = give.n, geom = "text") + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  ) 

#Plotting maxima
MaxAdBrood <- ddply(shortdata,.(Species2,Lake2),summarise,max.adbrood= max(PropAdBrood, na.rm = TRUE))
MaxAsexAdBrood <- ddply(shortdata,.(Species2,Lake2),summarise,max.asexadbrood= max(PropAsexAdBrood, na.rm = TRUE))
MaxBrood <- ddply(shortdata,.(Species2,Lake2),summarise,max.brood= max(PropBrood, na.rm = TRUE))
MaxAdBrood$max.percadbrood <- MaxAdBrood$max.adbrood*100
MaxAsexAdBrood$max.percasexadbrood <- MaxAsexAdBrood$max.asexadbrood*100
MaxBrood$max.percbrood <- MaxBrood$max.brood*100

#updating give.n function so that the output is at the right point on the next set of panels
give.n <- function(x){
  return(c(y = 9.5, label = length(x)))
}

panelc <- ggplot(MaxBrood, aes(x=Species2, y=max.percbrood)) + 
  geom_boxplot() + theme_bw() + 
  stat_summary(fun.data = give.n, geom = "text") +
  labs(y=expression(paste("Maximum\nprevalence (%)"))) + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=16),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=14)
  ) 

paneld <- ggplot(MaxBrood, aes(x=Lake2, y=max.percbrood)) + 
  geom_boxplot() + theme_bw() + 
  stat_summary(fun.data = give.n, geom = "text") + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  ) 

#updating give.n function so that the output is at the right point on the next set of panels
give.n <- function(x){
  return(c(y = 21.5, label = length(x)))
}

panele <- ggplot(MaxAsexAdBrood, aes(x=Species2, y=max.percasexadbrood)) + 
  geom_boxplot() + theme_bw() +
  stat_summary(fun.data = give.n, geom = "text") +
  labs(x="Host species", y=expression(paste("Maximum prevalence in\nasexual adults (%)"))) + theme(
    axis.title.x = element_text(size=16,vjust=-0.35),
    axis.title.y = element_text(size=16,vjust=0.9),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14)
  ) 

panelf <- ggplot(MaxAsexAdBrood, aes(x=Lake2, y=max.percasexadbrood)) + 
  geom_boxplot() + theme_bw() +
  stat_summary(fun.data = give.n, geom = "text") +
  labs(x="Lake") + theme(
    axis.title.x = element_text(size=16,vjust=-0.35),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_blank()
  ) 

library(grid)
grid.draw(rbind(ggplotGrob(panela),ggplotGrob(panelc),ggplotGrob(panele),size="last"))
grid.draw(rbind(ggplotGrob(panelb),ggplotGrob(paneld),ggplotGrob(panelf),size="last"))


tapply(MaxBrood$max.brood,list(MaxBrood$Species),max,na.rm=TRUE)
tapply(MaxAsexAdBrood$max.asexadbrood,list(MaxBrood$Species),max,na.rm=TRUE)

tapply(MaxBrood$max.brood,list(MaxBrood$Lake),max,na.rm=TRUE)
tapply(MaxAsexAdBrood$max.asexadbrood,list(MaxBrood$Lake),max,na.rm=TRUE)