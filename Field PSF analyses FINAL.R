library(car)
library(lme4)
library(ggplot2)
library(MuMIn)
library(reshape2)
library(nlme)
library(visreg)
library(lsmeans)
library(multcomp)
library(plyr)

#Longform data with all censuses
data <- read.csv("allPSFdataLong 170221.csv")
#Get rid of dead plants
data<-subset(data,Surv1==1)
data<-subset(data,is.na(Height)==F)

#Assign "block" to plots------------------------------------
data$Block<-c()
for(i in 1:length(data[,1])){
  if (data$Plot[i]==19|data$Plot[i]==23) data$Block[i]<-1
  if (data$Plot[i]==34|data$Plot[i]==36) data$Block[i]<-2
  if (data$Plot[i]==49|data$Plot[i]==50) data$Block[i]<-3
  if (data$Plot[i]==52|data$Plot[i]==53) data$Block[i]<-4
  if (data$Plot[i]==55|data$Plot[i]==56) data$Block[i]<-5
  if (data$Plot[i]==59|data$Plot[i]==61) data$Block[i]<-6
  if (data$Plot[i]==67|data$Plot[i]==69) data$Block[i]<-7
  if (data$Plot[i]==77|data$Plot[i]==79) data$Block[i]<-8
  if (data$Plot[i]==89|data$Plot[i]==90) data$Block[i]<-9
  if (data$Plot[i]==93|data$Plot[i]==95) data$Block[i]<-10
}

#Subset data --------------------------------------
dataLV<-subset(data,Type=="Live")
BOGRLV<-subset(dataLV,Transplant=="BOGR")
BOERLV<-subset(dataLV,Transplant=="BOER")

#Soil effects on plant biomass----------------
#BOGR
BOGR.bio<-lmer(logMass.harv~PSF*State+Census+(1|Block)+(1|ID),data=BOGRLV)
summary(BOGR.bio)
Anova(BOGR.bio,test="F")#census super sig, state 0.03, interaction 0.01
#Pairwise contrasts
lsmeans(BOGR.bio,pairwise~PSF|State)#D contrast p=0.23; S contrast p=0.01
#Check assumptions
qqPlot(resid(BOGR.bio))#not terrible
plot(fitted(BOGR.bio),resid(BOGR.bio))#fine

#BOER
BOER.bio<-lmer(logMass.harv~PSF*State+Census+(1|Block)+(1|ID),data=BOERLV)
summary(BOER.bio)
Anova(BOER.bio,test="F")#Census is super sig
#Pairwise contrasts
lsmeans(BOER.bio,pairwise~PSF|State)#D contrast p=0.76; S contrast p=0.04
#Check assumptions
qqPlot(resid(BOER.bio))#fine
plot(fitted(BOER.bio),resid(BOER.bio))#fine

#Soil effects on Height----------------
#BOGR
BOGR.bioH<-lmer(Height~PSF*State+Census+(1|Block)+(1|ID),data=BOGRLV)
Anova(BOGR.bioH,test="F")#census super sig, state 0.042 interaction not sig
lsmeans(BOGR.bioH,pairwise~PSF|State)#No sig contrasts
#Trend is the same though, taller plants in away than home under static.

#BOER
BOER.bioH<-lmer(Height~PSF*State+Census+(1|Block)+(1|ID),data=BOERLV)
Anova(BOER.bioH,test="F")#census super sig
lsmeans(BOER.bioH,pairwise~PSF|State)#S contrast 0.03

#Soil effects on Tillers----------------
#BOGR
BOGR.bioT<-lmer(Tiller~PSF*State+Census+(1|Block)+(1|ID),data=BOGRLV)
Anova(BOGR.bioT,test="F")#census super sig, State 0.05, PSF:state 0.02
lsmeans(BOGR.bioT,pairwise~PSF|State)#S contrast 0.02

#BOER
BOER.bioT<-lmer(Tiller~PSF*State+Census+(1|Block)+(1|ID),data=BOERLV)
Anova(BOER.bioT,test="F")#census super sig, PSF 0.08
lsmeans(BOER.bioT,pairwise~PSF|State)#S contrast 0.02

#Calculate PSF ln ratio using all censuses-------------------
#make unique plot-census identifier
BOERLV$plot.census<-paste(BOERLV$Plot,BOERLV$Census,sep=".")
BOGRLV$plot.census<-paste(BOGRLV$Plot,BOGRLV$Census,sep=".")

#Calculate PSF ln ratios
source("funcPSF.R")
PSF.BOGR<-func.PSF(BOGRLV)#183 rows
PSF.BOER<-func.PSF(BOERLV)#69 rows only (more BOER mortality)
#double-check that function ran properly (TRUE)
nrow(PSF.BOGR)==length(which(table(BOGRLV$plot.census)==2))
nrow(PSF.BOER)==length(which(table(BOERLV$plot.census)==2))

#Visualize
hist(PSF.BOGR$lnRatio)
hist(PSF.BOER$lnRatio)

#Analyze BOGR PSF
psfBOGR<-lmer(lnRatio~factor(State)+Census+(1|Block)+(1|Plot),data=PSF.BOGR)#Plot is the repeated measures, and Block is for spatial
qqPlot(resid(psfBOGR))#good enough
plot(psfBOGR)#not terrible but a hint of hump shape
Anova(psfBOGR)#State significant p=0.02
summary(psfBOGR)
lsmeans(psfBOGR,pairwise~State)#overlaps zero for D and does not for S
#Use means parameterization to figure out difference from zero
psfBOGR1<-lme(lnRatio~factor(State)-1,random=~1|Plot,data=PSF.BOGR)#use lme to get p values for means parameterization
summary(psfBOGR1)#S sig D insig
#Repeated measures but test against zero using means parameterization

#Analyze BOER PSF
psfBOER<-lmer(lnRatio~factor(State)+Census+(1|Block)+(1|Plot),data=PSF.BOER)
qqPlot(resid(psfBOER))#good
plot(psfBOER)#good
Anova(psfBOER)#State p=0.07, Census super sig
summary(psfBOER)#census and block both contribute to variance as random effects
#Use means parameterization to figure out difference from zero
psfBOER1<-lme(lnRatio~factor(State)-1,random=~1|Plot,data=PSF.BOER)
summary(psfBOER1)#both insig

#Calculate and analyze Bever's Is------------------------------------
#Combine the two output datasets from previous function
IsData<-rbind.data.frame(PSF.BOER,PSF.BOGR)
IsData$Plant<-c(rep("BOER",nrow(PSF.BOER)),rep("BOGR",nrow(PSF.BOGR)))

#Calculate Is
Is<-IsData[1,-c(1,7,8)]#initialize dataframe
Is$Is<-0
for (i in 1: length(unique(IsData$plot.census))){
  subdata<-subset(IsData,plot.census==unique(IsData$plot.census)[i])
  if (nrow(subdata)!=2) next #don't bother calculating Is if no living quartet
  else {
    #somehow this clunky thing is the only way the data can be added
    result<-c(subdata$plot.census[1],subdata$Plot[1],subdata$Census[1],subdata$Block[1],subdata$State[1],sum(subdata$lnRatio))
    Is<-rbind(Is,result)
  }
}
#Get rid of dummy first row
Is<-Is[-1,]#69 rows
#double-check that function ran properly (TRUE)
nrow(Is)==length(which(table(IsData$plot.census)==2))

#Analyze Is
hist(as.numeric(Is$Is))#Somehow all the columns in Is are "characters" now
#Good distribution
m.Is<-lmer(as.numeric(Is)~factor(State)+as.numeric(Census)+(1|Block)+(1|Plot),data=Is)#again Plot for repeated measures
qqPlot(resid(m.Is))#good
plot(m.Is)#there's a funnel...
aggregate(as.numeric(Is$Is),list(factor(Is$State)),var)#within 4 fold
aggregate(as.numeric(Is$Is),list(as.numeric(Is$Census)),var)#large variation here
Anova(m.Is)#state sig
summary(m.Is)
lsmeans(m.Is,pairwise~State)#CL overlaps zero for D but not for S
#Use means parameterization to figure out difference from zero
m1.Is<-lme(as.numeric(Is)~factor(State)-1,random=~1|Plot,data=Is)#use lme to get p values for means parameterization
summary(m1.Is)#S sig D insig

#Static patches conspecific freq using rare/common ----------------------------
#Soil effects
#subset data by static only
BOGRstat<-subset(BOGRLV,State=="S")
BOERstat<-subset(BOERLV,State=="S")
#BOGR
BOGR.S<-lmer(logMass.harv~PSF*Dom_SPP+Census+(1|Block)+(1|ID),data=BOGRstat)
plot(BOGR.S)#fine
Anova(BOGR.S,test="F")#PSF 0.055, census 0.07, no sig interaction
lsmeans(BOGR.S,pairwise~PSF|Dom_SPP)#But significant contrast (p=0.02) when BOGR rare but not when common
#BOER
BOER.S<-lmer(logMass.harv~PSF*Dom_SPP+Census+(1|Block)+(1|ID),data=BOERstat)
plot(BOER.S)#fine
Anova(BOER.S,test="F")#PSF 0.09, census<0.001, interaction 0.04
lsmeans(BOER.S,pairwise~PSF|Dom_SPP)#significant contrast (p=0.006) when BOER rare but not when common

#PSF response
#subset data by static only
PSF.BOGRS<-subset(PSF.BOGR,State=="S")
PSF.BOERS<-subset(PSF.BOER,State=="S")
#Analyze BOGR PSF
psfBOGRS<-lmer(lnRatio~factor(Dom_spp)+Census+(1|Plot),data=PSF.BOGRS)#Plot is the repeated measures, and Block is for spatial
qqPlot(resid(psfBOGRS))#good enough
plot(psfBOGRS)#good
Anova(psfBOGRS)#only census sig
#Use means parameterization to figure out difference from zero
psfBOGRS1<-lme(lnRatio~factor(Dom_spp)-1,random=~1|Plot,data=PSF.BOGRS)#use lme to get p values for means parameterization
summary(psfBOGRS1)#rare p=0.055 common p=0.66
#Analyze BOER PSF
psfBOERS<-lmer(lnRatio~factor(Dom_spp)+Census+(1|Plot),data=PSF.BOERS)#Plot is the repeated measures, and Block is for spatial
qqPlot(resid(psfBOERS))#good enough
plot(psfBOERS)#good
Anova(psfBOERS)#only census sig
#Use means parameterization to figure out difference from zero
psfBOERS1<-lme(lnRatio~factor(Dom_spp)-1,random=~1|Plot,data=PSF.BOERS)#use lme to get p values for means parameterization
summary(psfBOERS1)#rare p=0.16 common p=0.66

#Look at the trends
aggregate(PSF.BOGRS$lnRatio~PSF.BOGRS$Dom_spp,FUN="mean")
aggregate(PSF.BOERS$lnRatio~PSF.BOERS$Dom_spp,FUN="mean")

#Survival ---------------------------------------
#Who survived all years? 
survID<-data$ID[which(data$Census==15&data$Surv2==1)]
surv<-subset(data,Census==1)
for(i in 1:length(surv[,1])){
  if (surv$ID[i]%in%survID) surv$Surv15[i]<-1
  else surv$Surv15[i]<-0
}

#Survival with species separate 179224
BOERsurv<-subset(surv,Transplant=="BOER")
BOGRsurv<-subset(surv,Transplant=="BOGR")

#BOER
surv.boer<-glmer(Surv15~State*PSF+(1|Block),family="binomial",data=BOERsurv)
#warning message (this is because for D-Away, all values=0)
#Diagnositic plots are kinda terrible but there's not much to do for binomial data
plot(surv.boer)
qqPlot(resid(surv.boer))
#Try for results anyway
summary(surv.boer)#All p's are 0.95
Anova(surv.boer)#All insig

#Pairwise contrasts
head(model.matrix(surv.boer))# (Intercept) StateS PSFHome PSFSterile StateS:PSFHome StateS:PSFSterile
contrast<-rbind("AwayD-HomeD"=c(1,0,-1,0,0,0),
                "AwayD-SterileD"=c(1,0,0,-1,0,0),
                "HomeD-SterileD"=c(1,0,1,-1,0,0),
                "AwayS-HomeS"=c(1,1,0,0,-1,0),
                "AwayS-SterileS"=c(1,1,0,0,0,-1),
                "HomeS-SterileS"=c(1,0,0,0,1,-1))
multiStat<-glht(surv.boer,contrast)
summary(multiStat)

#ggplot
plotdata<-aggregate(BOERsurv$Surv15~BOERsurv$State*BOERsurv$PSF,FUN="mean")
sd<-aggregate(BOERsurv$Surv15~BOERsurv$State*BOERsurv$PSF,FUN="sd")[,3]
df<-aggregate(BOERsurv$Surv15~BOERsurv$State*BOERsurv$PSF,FUN="length")[,3]
se<-sd/sqrt(df)
plotdata$se<-se
colnames(plotdata)<-c("State","PSF","Survival","se")
theme_set(theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
dodge <- position_dodge(width=0.5)
gp <- ggplot(plotdata, aes(x=PSF, y=Survival,group=State, colour=State))
gp  + geom_point(size=7,position=dodge)+
  geom_errorbar(aes(ymin=Survival-se, ymax=Survival+se), width=.2,position=dodge)+                    # Width of the error bars
  theme(text=element_text(size=20))+ 
  xlab("\nFeedback environment") + ylab("Mean BOER survival\n")+
  ylim(-0.01,0.8)

#BOGR
surv.bogr<-glmer(Surv15~State*PSF+(1|Block),family="binomial",data=BOGRsurv)
#Diagnositic plots are kinda terrible but there's not much to do for binomial data
plot(surv.bogr)
qqPlot(resid(surv.bogr))
#Try for results anyway
summary(surv.bogr)#all insig
Anova(surv.boer)#All insig
#multi comparisons
head(model.matrix(surv.bogr))# (Intercept) StateS PSFHome PSFSterile StateS:PSFHome StateS:PSFSterile
contrast<-rbind("AwayD-HomeD"=c(1,0,-1,0,0,0),
                "AwayD-SterileD"=c(1,0,0,-1,0,0),
                "HomeD-SterileD"=c(1,0,1,-1,0,0),
                "AwayS-HomeS"=c(1,1,0,0,-1,0),
                "AwayS-SterileS"=c(1,1,0,0,0,-1),
                "HomeS-SterileS"=c(1,0,0,0,1,-1))
multiStat<-glht(surv.bogr,contrast)
summary(multiStat)

#ggplot
plotdata<-aggregate(BOGRsurv$Surv15~BOGRsurv$State*BOGRsurv$PSF,FUN="mean")
sd<-aggregate(BOGRsurv$Surv15~BOGRsurv$State*BOGRsurv$PSF,FUN="sd")[,3]
df<-aggregate(BOGRsurv$Surv15~BOGRsurv$State*BOGRsurv$PSF,FUN="length")[,3]
se<-sd/sqrt(df)
plotdata$se<-se
colnames(plotdata)<-c("State","PSF","Survival","se")
theme_set(theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
dodge <- position_dodge(width=0.5)
gp <- ggplot(plotdata, aes(x=PSF, y=Survival,group=State, colour=State))
gp  + geom_point(size=7,position=dodge)+
  geom_errorbar(aes(ymin=Survival-se, ymax=Survival+se), width=.2,position=dodge)+                    # Width of the error bars
  theme(text=element_text(size=20))+ 
  xlab("\nFeedback environment") + ylab("Mean BOGR survival\n")+
  ylim(0,0.8)

#Field germination--------------------------------------
field <- read.csv("Field germ with survival 170511.csv")

#BOER only 
boer<-subset(field,Transplant=="BOER")
m1<-glmer(maxgerm~State*PSF+(1|Block),data=boer,family=poisson)
qqPlot(resid(m1))#looks decent; warnings about NAs...not sure why
plot(m1)#not bad
summary(m1)#StatS:Sterile 0.02 (StateS p=0.09)
#Block is tiny part of the variance
Anova(m1)#can't do F test
leveneTest(glm(maxgerm~State*PSF,data=boer,family=poisson))#actully fine
#Conlusion: BOER germ is higher in Sterile than Away, possibly also Home, but only in static plots
#Pairwise contrasts
head(model.matrix(m1))# (Intercept) StateS PSFHome PSFSterile StateS:PSFHome StateS:PSFSterile
contrast<-rbind("D.Away-D.Sterile"=c(1,0,0,-1,0,0),
                "D.Home-D.Sterile"=c(1,0,1,-1,0,0),
                "D.Away-D.Home"=c(1,0,-1,0,0,0),
                "S.Away-S.Sterile"=c(1,1,0,0,0,-1),
                "S.Home-S.Sterile"=c(1,0,0,0,1,-1),
                "S.Away-S.Home"=c(1,1,0,0,-1,0))
multi<-glht(m1,contrast)
summary(multi)#with errors

#Pretty plot
boer$levels<-paste(boer$State,boer$PSF)
df<-aggregate(boer$maxgerm,list(boer$levels), mean)
df$se <- with(boer,aggregate(boer$maxgerm, list(boer$levels), 
                             function(x) sd(x)/sqrt(length(x))))[,2]
theme_set(theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
Xgroup<-c("Away","Home","Sterile","Away","Home","Sterile")
State<-c("Dynamic","Dynamic","Dynamic","Static","Static","Static")
dodge <- position_dodge(width=0.5)
gp <- ggplot(df, aes(x=Xgroup, y=x,group=State, colour=State))
gp  + geom_point(size=7,position=dodge)+
  geom_errorbar(aes(ymin=x-se, ymax=x+se), width=.2,position=dodge)+                    # Width of the error bars
  theme(text=element_text(size=20))+ 
  xlab("Feedback environment") + ylab("Average BOER seeds germinated")+
  #ggtitle("B. eriopda")+
  ylim(0,4.2)

#BOGR only 
bogr<-subset(data,Transplant=="BOGR")
m2<-glmer(maxgerm~State*PSF+(1|Block),data=bogr,family=poisson)
qqPlot(resid(m2))#looks decent
plot(m2)#not bad
summary(m2)#Still no sig effects (StateS p=0.09)
#Block is pretty big part of the variance
Anova(m2)#State is 0.07 here
leveneTest(glm(maxgerm~State*PSF,data=bogr,family=poisson))#actully fine
#Conlusion: BOGR germ is trending lower in S vs. D
#Pairwise contrasts
head(model.matrix(m2))# (Intercept) StateS PSFHome PSFSterile StateS:PSFHome StateS:PSFSterile
contrast<-rbind("D.Away-D.Sterile"=c(1,0,0,-1,0,0),
                "D.Home-D.Sterile"=c(1,0,1,-1,0,0),
                "D.Away-D.Home"=c(1,0,-1,0,0,0),
                "S.Away-S.Sterile"=c(1,1,0,0,0,-1),
                "S.Home-S.Sterile"=c(1,0,0,0,1,-1),
                "S.Away-S.Home"=c(1,1,0,0,-1,0))
multi<-glht(m2,contrast)
summary(multi)#The D.Away and D.Home comparisons against sterile are sig; with errors
#lsmeans
lsmeans(m2,pairwise~PSF|State)#won't run right because model is saturated 
#Try running separate models 
m2.s<-glmer(maxgerm~PSF+(1|Block),data=subset(bogr,State=="S"),family=poisson)
m2.d<-glmer(maxgerm~PSF+(1|Block),data=subset(bogr,State=="D"),family=poisson)
lsmeans(m2.s,pairwise~PSF)

#Pretty plot
bogr$levels<-paste(bogr$State,bogr$PSF)
df<-aggregate(bogr$maxgerm,list(bogr$levels), mean)
df$se <- with(bogr,aggregate(bogr$maxgerm, list(bogr$levels), 
                             function(x) sd(x)/sqrt(length(x))))[,2]
theme_set(theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
Xgroup<-c("Away","Home","Sterile","Away","Home","Sterile")
State<-c("Dynamic","Dynamic","Dynamic","Static","Static","Static")
dodge <- position_dodge(width=0.5)
gp <- ggplot(df, aes(x=Xgroup, y=x,group=State, colour=State))
gp  + geom_point(size=7,position=dodge)+
  geom_errorbar(aes(ymin=x-se, ymax=x+se), width=.2,position=dodge)+  
  theme(text=element_text(size=20))+ 
  xlab("Feedback environment") + ylab("Average BOGR seeds germinated")+
  #ggtitle("B. gracilis")+
  ylim(0,4.2)

#Pot germination experiment------------------------------
pot <- read.csv("Pot germination 160922.csv")

#Add extra variable: total germination at end of experiment
pot$BOERall<-pot$BOERgerm+pot$BOERdeath
pot$BOGRall<-pot$BOGRgerm+pot$BOGRdeath

#analyze final germination
final<-subset(pot,Day==28)

#Check distribution of pot
hist(final$BOERall)#bad
hist(final$BOGRall)#decent
hist(log(final$BOERall+0.05))#better
hist(log(final$BOGRall+0.05))#little weird

#BOER
boer1<-lme(log(BOERall+0.05)~State*InocSPP,random=~1|Block,data=final)
qqPlot(boer1$resid) #fine
summary(boer1)#nothing
boer2<-lme(BOERall~State*InocSPP,random=~1|Block,data=final)#report this
qqPlot(boer2$resid) #fine
summary(boer2)#still nothing

#BOGR
bogr1<-lme(log(BOGRall+0.05)~State*InocSPP,random=~1|Block,data=final)
qqPlot(bogr1$resid) #fine
summary(bogr1)#state is sig...lower germination from Static soils
bogr2<-lme(BOGRall~State*InocSPP,random=~1|Block,data=final)#report this
qqPlot(bogr2$resid) #better
summary(bogr2)#state is sig (p=0.004)...lower germination from Static soils



