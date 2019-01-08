library(ggplot2)
library(grid)
library(gridExtra)
library(lme4)
library(emmeans)
library(nlme)

#Data processing-------------
#Longform data with all censuses
data <- read.csv("allPSFdataLong 170221.csv")
#Get rid of dead plants
data<-subset(data,Surv1==1)
data<-subset(data,is.na(Height)==F)

#Assign blocks
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

#Subset data
dataLV<-subset(data,Type=="Live")
BOGRLV<-subset(dataLV,Transplant=="BOGR")
BOERLV<-subset(dataLV,Transplant=="BOER")
dataS<-subset(dataLV,State=="S")
dataD<-subset(dataLV,State=="D")

#make unique plot-census identifier
BOERLV$plot.census<-paste(BOERLV$Plot,BOERLV$Census,sep=".")
BOGRLV$plot.census<-paste(BOGRLV$Plot,BOGRLV$Census,sep=".")

#Calculate PSFs
source("funcPSF.R")
PSF.BOGR<-func.PSF(BOGRLV)#183 rows
PSF.BOER<-func.PSF(BOERLV)#69 rows only (more BOER mortality)

#Combine the two output datasets from previous function
IsData<-rbind.data.frame(PSF.BOER,PSF.BOGR)
IsData$Plant<-c(rep("BOER",nrow(PSF.BOER)),rep("BOGR",nrow(PSF.BOGR)))

#Separate out based on Static v Dynamic
PSF.S<-subset(IsData,State=="S")
PSF.D<-subset(IsData,State=="D")

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

#Static only
BOGRstat<-subset(BOGRLV,State=="S")
BOERstat<-subset(BOERLV,State=="S")
PSF.BOGRS<-subset(PSF.BOGR,State=="S")
PSF.BOERS<-subset(PSF.BOER,State=="S")
#Separate out based on rare/common
rareData<-rbind.data.frame(subset(BOGRstat,Dom_SPP=="BOER"),subset(BOERstat,Dom_SPP=="BOGR"))
commonData<-rbind.data.frame(subset(BOGRstat,Dom_SPP=="BOGR"),subset(BOERstat,Dom_SPP=="BOER"))
PSF.R<-rbind.data.frame(subset(PSF.BOGRS,Dom_spp=="BOER"),subset(PSF.BOERS,Dom_spp=="BOGR"))
PSF.R$Plant<-c(rep("BOGR",nrow(subset(PSF.BOGRS,Dom_spp=="BOER"))),rep("BOER",nrow(subset(PSF.BOERS,Dom_spp=="BOGR"))))
PSF.C<-rbind.data.frame(subset(PSF.BOGRS,Dom_spp=="BOGR"),subset(PSF.BOERS,Dom_spp=="BOER"))
PSF.C$Plant<-c(rep("BOGR",nrow(subset(PSF.BOGRS,Dom_spp=="BOGR"))),rep("BOER",nrow(subset(PSF.BOERS,Dom_spp=="BOER"))))

#Growth responses ------------------------------
#Fit biomass models for fitted SEs
mod.biomassBOER<-lmer(exp(logMass.harv)~PSF*State+Census+(1|Block)+(1|ID),data=BOERLV)
lsm.biomassBOER<-lsmeans(mod.biomassBOER,pairwise~PSF|State)
s.biomassBOER<-summary(lsm.biomassBOER$lsmeans)
mod.biomassBOGR<-lmer(exp(logMass.harv)~PSF*State+Census+(1|Block)+(1|ID),data=BOGRLV)
lsm.biomassBOGR<-lsmeans(mod.biomassBOGR,pairwise~PSF|State)
s.biomassBOGR<-summary(lsm.biomassBOGR$lsmeans)
#Fit PSF models for fitted SEs
#Means parameterization to figure out difference from zero
psfBOGR1<-lme(lnRatio~factor(State)-1,random=~1|Plot,data=PSF.BOGR)#use lme to get p values for means parameterization
s.psfBOGR1<-summary(psfBOGR1)#S sig D insig
psfBOER1<-lme(lnRatio~factor(State)-1,random=~1|Plot,data=PSF.BOER)#use lme to get p values for means parameterization
s.psfBOER1<-summary(psfBOER1)#S sig D insig
#Fit Is model
#Use means parameterization to figure out difference from zero
m1.Is<-lme(as.numeric(Is)~factor(State)-1,random=~1|Plot,data=Is)#use lme to get p values for means parameterization
s.Is<-summary(m1.Is)#S sig D insig

#Dynamic growth barchart panel
#Make data frame
dataD$levels<-paste(dataD$Transplant,dataD$PSF)
df.D<-aggregate(exp(dataD$logMass.harv),list(dataD$levels), mean)
df.D$fittedSE<-c(s.biomassBOER$SE[1:2],s.biomassBOGR$SE[1:2])
theme_set(theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
df.D$PSF<-c("Heterospecific","Conspecific","Heterospecific","Conspecific")
df.D$Transplant<-c("Black grama","Black grama","Blue grama","Blue grama")
dodge <- position_dodge(width=0.6)
#Plot
bar.D <- ggplot(df.D, aes(x=Transplant, y=x,group=PSF,fill=PSF))+
  geom_bar(stat="identity",position=dodge,width=0.6,colour="black")+
  geom_errorbar(aes(ymin=x, ymax=x+fittedSE), colour="black",width=.3,position=dodge)+   
  theme(text=element_text(size=15),legend.position = c(0.7, 0.8))+ 
  xlab("Plant species") + ylab("Aboveground biomass (g)\n")+
  scale_fill_manual(values=c("grey","white"),name="Soil type"
                    ,labels=c("Conspecific","Heterospecific")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.2))+
  annotate("text",x=0.55,y=2.05,label="E",size=8)

#Static growth barchart panel
#Make data frame
dataS$levels<-paste(dataS$Transplant,dataS$PSF)
df.S<-aggregate(exp(dataS$logMass.harv),list(dataS$levels), mean)
df.S$fittedSE<-c(s.biomassBOER$SE[3:4],s.biomassBOGR$SE[3:4])
theme_set(theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
df.S$PSF<-c("Heterospecific","Conspecific","Heterospecific","Conspecific")
df.S$Transplant<-c("Black grama","Black grama","Blue grama","Blue grama")
dodge <- position_dodge(width=0.6)
bar.S <- ggplot(df.S, aes(x=Transplant, y=x,group=PSF,fill=PSF))+
  geom_bar(stat="identity",position=dodge,width=0.6,colour="black")+
  geom_errorbar(aes(ymin=x, ymax=x+fittedSE), colour="black",width=.3,position=dodge)+   
  theme(text=element_text(size=15),legend.position = "none")+ 
  xlab("Plant species")  + ylab("Aboveground biomass (g)\n")+
  scale_fill_manual(values=c("grey","white")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.2))+
  annotate("text",x=0.55,y=2.05,label="F",size=8)+
  annotate("text",x=1,y=1,label="*",size=10)+
  annotate("text",x=2,y=2.1,label="*",size=10)

#Dynamic PSF points panel
psf.D<-aggregate(PSF.D$lnRatio~PSF.D$Plant,FUN="mean")
psf.D$se <- c(s.psfBOER1$tTable[1,2],s.psfBOGR1$tTable[1,2])
colnames(psf.D)<-c("Plant","lsmean","se")
theme_set(theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
points.D.ytitle <- expression(paste(" PSF (",italic("ln"),"-ratio)"))
points.D <- ggplot(psf.D, aes(x=Plant, y=lsmean))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=lsmean-se, ymax=lsmean+se), width=.2)+    
  theme(text=element_text(size=15))+ 
  geom_hline(linetype="dashed",yintercept=0)+
  xlab(" ") + ylab(points.D.ytitle)+
  scale_x_discrete(labels=c("Black grama", "Blue grama"))+
  scale_y_continuous(expand = c(0, 0), limits = c(-0.7, 0.5))+
  annotate("text",x=0.53,y=0.3,label="C",size=8)

#Static PSF points panel
psf.S<-aggregate(PSF.S$lnRatio~PSF.S$Plant,FUN="mean")
psf.S$se <- c(s.psfBOER1$tTable[2,2],s.psfBOGR1$tTable[2,2])
colnames(psf.S)<-c("Plant","lsmean","se")
theme_set(theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
points.S.ytitle <- expression(paste(" PSF (",italic("ln"),"-ratio)"))
points.S <- ggplot(psf.S, aes(x=Plant, y=lsmean))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=lsmean-se, ymax=lsmean+se), width=.2)+    
  theme(text=element_text(size=15))+ 
  geom_hline(linetype="dashed",yintercept=0)+
  xlab(" ") + ylab(points.S.ytitle)+
  scale_x_discrete(labels=c("Black grama", "Blue grama"))+
  scale_y_continuous(expand = c(0, 0), limits = c(-0.8, 0.5))+
  annotate("text",x=0.53,y=0.3,label="D",size=8)+
  annotate("text",x=2,y=-0.3,label="*",size=10)

#Is plot data
df.IS<-aggregate(as.numeric(Is$Is)~factor(Is$State),FUN="mean")
df.IS$se <- s.Is$tTable[,2]
colnames(df.IS)<-c("State","Ismean","se")
Is.ytitle<-expression(italic("I")[s])

#Is panel separated
D.Is <- ggplot(subset(df.IS,State=="D"), aes(x=State, y=as.numeric(Ismean)))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=as.numeric(Ismean)-as.numeric(se), ymax=as.numeric(Ismean)+as.numeric(se)), width=.1)+    
  theme(text=element_text(size=15),plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(margin=margin(t=0,r=12,b=0,l=0,unit="pt")))+ 
  geom_hline(linetype="dashed",yintercept=0)+
  xlab("")+  ylab(Is.ytitle)+ggtitle("Dynamic patches")+
  scale_y_continuous(expand = c(0, 0), limits = c(-1, 0.7))+
  scale_x_discrete(labels= c("Dynamic"))+
  annotate("text",x=0.56,y=0.37,label="A",size=8)

S.Is <- ggplot(subset(df.IS,State=="S"), aes(x=State, y=as.numeric(Ismean)))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=as.numeric(Ismean)-as.numeric(se), ymax=as.numeric(Ismean)+as.numeric(se)), width=.1)+    
  theme(text=element_text(size=15),plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(margin=margin(t=0,r=12,b=0,l=0,unit="pt")))+ 
  geom_hline(linetype="dashed",yintercept=0)+
  xlab("") +  ylab(Is.ytitle)+ggtitle("Static patches")+
  scale_y_continuous(expand = c(0, 0), limits = c(-1.2, 0.7))+
  scale_x_discrete(labels= c("Static"))+
  annotate("text",x=0.56,y=0.37,label="B",size=8)+
  annotate("text",x=1,y=-0.5,label="*",size=10)

#Combine panels
lay2<-rbind(c(1,2),
            c(3,4),
            c(5,6),
            c(5,6))
grid.arrange(D.Is,S.Is,points.D,points.S,bar.D,bar.S,layout_matrix = lay2)

#Make biomass panel with sterile controls
BOER<-subset(data,Transplant=="BOER")
BOGR<-subset(data,Transplant=="BOGR")
alldatD<-subset(data,State=="D")
alldatS<-subset(data,State=="S")

#Fit biomass models for fitted SEs
mod.BOER<-lmer(exp(logMass.harv)~PSF*State+Census+(1|Block)+(1|ID),data=BOER)
lsm.BOER<-lsmeans(mod.BOER,pairwise~PSF|State)
s.BOER<-summary(lsm.BOER$lsmeans)
mod.BOGR<-lmer(exp(logMass.harv)~PSF*State+Census+(1|Block)+(1|ID),data=BOGR)
lsm.BOGR<-lsmeans(mod.BOGR,pairwise~PSF|State)
s.BOGR<-summary(lsm.BOGR$lsmeans)

#Dynamic growth barchart panel
#Make data frame
alldatD$levels<-paste(alldatD$Transplant,alldatD$PSF)
df.D<-aggregate(exp(alldatD$logMass.harv),list(alldatD$levels), mean)
df.D$fittedSE<-c(s.BOER$SE[1:3],s.BOGR$SE[1:3])
theme_set(theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
df.D$PSF<-c("Heterospecific","Conspecific","No fungi","Heterospecific","Conspecific","No fungi")
df.D$Transplant<-c("Black grama","Black grama","Black grama","Blue grama","Blue grama","Blue grama")
dodge <- position_dodge(width=0.6)
#Plot
bar.D <- ggplot(df.D, aes(x=Transplant, y=x,group=PSF,fill=PSF))+
  geom_bar(stat="identity",position=dodge,width=0.6,colour="black")+
  geom_errorbar(aes(ymin=x, ymax=x+fittedSE), colour="black",width=.3,position=dodge)+   
  theme(text=element_text(size=15),legend.position = c(0.7, 0.8))+ 
  xlab("Plant species") + ylab("Aboveground biomass (g)\n")+
  scale_fill_manual(values=c("black","grey","white"),name="Feedback env."
                    ,labels=c("Conspecific","Heterospecific","No fungi")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.2))+
  annotate("text",x=0.55,y=2.05,label="A",size=8)

#Static growth barchart panel
#Make data frame
alldatS$levels<-paste(alldatS$Transplant,alldatS$PSF)
df.S<-aggregate(exp(alldatS$logMass.harv),list(alldatS$levels), mean)
df.S$fittedSE<-c(s.BOER$SE[4:6],s.BOGR$SE[4:6])
theme_set(theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
df.S$PSF<-c("Heterospecific","Conspecific","No fungi","Heterospecific","Conspecific","No fungi")
df.S$Transplant<-c("Black grama","Black grama","Black grama","Blue grama","Blue grama","Blue grama")
dodge <- position_dodge(width=0.6)
#Plot
bar.S <- ggplot(df.S, aes(x=Transplant, y=x,group=PSF,fill=PSF))+
  geom_bar(stat="identity",position=dodge,width=0.6,colour="black")+
  geom_errorbar(aes(ymin=x, ymax=x+fittedSE), colour="black",width=.3,position=dodge)+   
  theme(text=element_text(size=15),legend.position = "none")+ 
  xlab("Plant species") + ylab("Aboveground biomass (g)\n")+
  scale_fill_manual(values=c("black","grey","white")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.2))+
  annotate("text",x=0.55,y=2.05,label="B",size=8)

#Combine panels
lay3<-rbind(c(1,2))
grid.arrange(bar.D,bar.S,layout_matrix = lay3)


#Rare/common----------------------------
#exponentiated biomass models for fitted SEs
BOGR.S<-lmer(exp(logMass.harv)~PSF*Dom_SPP+Census+(1|Block)+(1|ID),data=BOGRstat)
lsm.BOGR.S<-lsmeans(BOGR.S,pairwise~PSF|Dom_SPP)#But significant contrast (p=0.02) when BOGR rare but not when common
s.BOGR.S<-summary(lsm.BOGR.S$lsmeans)
BOER.S<-lmer(exp(logMass.harv)~PSF*Dom_SPP+Census+(1|Block)+(1|ID),data=BOERstat)
lsm.BOER.S<-lsmeans(BOER.S,pairwise~PSF|Dom_SPP)#significant contrast (p=0.006) when BOER rare but not when common
s.BOER.S<-summary(lsm.BOER.S$lsmeans)
#models for ln ratio PSF fitted SEs
psfBOGRS1<-lme(lnRatio~factor(Dom_spp)-1,random=~1|Plot,data=PSF.BOGRS)#use lme to get p values for means parameterization
s.psfBOGRS1<-summary(psfBOGRS1)#rare p=0.055 common p=0.66
psfBOERS1<-lme(lnRatio~factor(Dom_spp)-1,random=~1|Plot,data=PSF.BOERS)#use lme to get p values for means parameterization
s.psfBOERS1<-summary(psfBOERS1)#rare p=0.16 common p=0.66

#Rare growth barchart panel
rareData$levels<-paste(rareData$Transplant,rareData$PSF)
df.R<-aggregate(exp(rareData$logMass.harv),list(rareData$levels), mean)
df.R$se <- c(s.BOER.S$SE[3:4],s.BOGR.S$SE[1:2])
theme_set(theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
df.R$PSF<-c("Heterospecific","Conspecific","Heterospecific","Conspecific")
df.R$Transplant<-c("Black grama","Black grama","Blue grama","Blue grama")
dodge <- position_dodge(width=0.6)
bar.R <- ggplot(df.R, aes(x=Transplant, y=x,group=PSF,fill=PSF))+
  geom_bar(stat="identity",position=dodge,width=0.6,colour="black")+
  geom_errorbar(aes(ymin=x, ymax=x+se), colour="black",width=.3,position=dodge)+   
  theme(text=element_text(size=15),legend.position = "none",
        axis.title.y = element_text(margin=margin(t=0,r=23,b=0,l=0,unit="pt")))+ 
  xlab("Plant species") + ylab("Aboveground biomass (g)")+
  scale_fill_manual(values=c("grey","white"),name="Type"
                    ,labels=c("Conspecific","Heterospecific")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3))+
  annotate("text",x=0.55,y=2.8,label="C",size=8)+
  annotate("text",x=1,y=1,label="*",size=10)+
  annotate("text",x=2,y=2.8,label="*",size=10)

#Common growth barchart panel
commonData$levels<-paste(commonData$Transplant,commonData$PSF)
df.C<-aggregate(exp(commonData$logMass.harv),list(commonData$levels), mean)
df.C$se <- c(s.BOER.S$SE[1:2],s.BOGR.S$SE[3:4])
theme_set(theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
df.C$PSF<-c("Heterospecific","Conspecific","Heterospecific","Conspecific")
df.C$Transplant<-c("Black grama","Black grama","Blue grama","Blue grama")
dodge <- position_dodge(width=0.6)
bar.C <- ggplot(df.C, aes(x=Transplant, y=x,group=PSF,fill=PSF))+
  geom_bar(stat="identity",position=dodge,width=0.6,colour="black")+
  geom_errorbar(aes(ymin=x, ymax=x+se), colour="black",width=.3,position=dodge)+   
  theme(text=element_text(size=15),legend.position = c(0.7,0.87),
        axis.title.y = element_text(margin=margin(t=0,r=23,b=0,l=0,unit="pt")))+ 
  xlab("Plant species") + ylab("Aboveground biomass (g)")+
  scale_fill_manual(values=c("grey","white"),name="Soil type"
                    ,labels=c("Conspecific","Heterospecific")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3))+
  annotate("text",x=0.55,y=2.8,label="D",size=8)

#Rare PSF points panel
psf.R<-aggregate(PSF.R$lnRatio~PSF.R$Plant,FUN="mean")
psf.R$se <- c(s.psfBOERS1$tTable[2,2],s.psfBOGRS1$tTable[1,2])
colnames(psf.R)<-c("Plant","lsmean","se")
theme_set(theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
points.R.ytitle <- expression(paste(" PSF (",italic("ln"),"-ratio)"))
points.R <- ggplot(psf.R, aes(x=Plant, y=lsmean))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=lsmean-se, ymax=lsmean+se), width=.2)+    
  theme(text=element_text(size=15),plot.title = element_text(hjust = 0.5))+ 
  geom_hline(linetype="dashed",yintercept=0)+
  xlab(" ") + ylab(points.R.ytitle)+ggtitle("Rare")+
  scale_x_discrete(labels=c("Black grama", "Blue grama"))+
  scale_y_continuous(expand = c(0, 0), limits = c(-1.5, 0.5))+
  annotate("text",x=0.55,y=0.25,label="A",size=8)

#Common PSF points panel
psf.C<-aggregate(PSF.C$lnRatio~PSF.C$Plant,FUN="mean")
psf.C$se <- c(s.psfBOERS1$tTable[1,2],s.psfBOGRS1$tTable[2,2])
colnames(psf.C)<-c("Plant","lsmean","se")
theme_set(theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
points.C.ytitle <- expression(paste(" PSF (",italic("ln"),"-ratio)"))
points.C <- ggplot(psf.C, aes(x=Plant, y=lsmean))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=lsmean-se, ymax=lsmean+se), width=.2)+    
  theme(text=element_text(size=15),plot.title = element_text(hjust = 0.5))+ 
  geom_hline(linetype="dashed",yintercept=0)+
  xlab(" ") + ylab(points.C.ytitle)+ggtitle("Common")+
  scale_x_discrete(labels=c("Black grama", "Blue grama"))+
  scale_y_continuous(expand = c(0, 0), limits = c(-1.5, 0.5))+
  annotate("text",x=0.55,y=0.25,label="B",size=8)

#Combine panels
lay3<-rbind(c(1,2),
            c(3,4),
            c(3,4))
grid.arrange(points.R,points.C,bar.R,bar.C,layout_matrix = lay3)

