

library(car)
library(MuMIn)
library(nlme)
library(visreg)

# Glom 4m sections --------------------------------------------------------
rawdata <- read.csv("sev004_transects_11292012_1.txt")
dwdata<-subset(rawdata,rawdata$location=="DW") #subset deep well data
#get rid of data error
dwdata<-subset(dwdata, dwdata$start<400)
dwdata<-subset(dwdata, dwdata$stop<401)
#Remake interval data to include total biomass
newdata<-dwdata[,c(1,4,5,6,7)] #use only relevant columns
newdata$interval<-rep(0,length(newdata$year))
extra<-newdata[1,] #for the new split rows
for (i in 1:length(newdata$year)) {
  if (floor(newdata$start[i]/4)== floor(newdata$stop[i]/4))
    newdata$interval[i]<-floor(newdata$start[i]/4)+1
  else{
    newdata<-rbind(newdata,extra)
    newdata$year[length(newdata$year)]<-newdata$year[i]
    newdata$species_code[length(newdata$year)]<-newdata$species_code[i]
    newdata$period[length(newdata$year)]<-newdata$period[i]
    newdata$start[length(newdata$year)]<-floor(newdata$stop[i]/4)*4
    newdata$stop[length(newdata$year)]<-newdata$stop[i]
    newdata$interval[length(newdata$year)]<-floor(newdata$start[length(newdata$year)]/4)+1
    newdata$stop[i]<-floor(newdata$stop[i]/4)*4
    newdata$interval[i]<-floor(newdata$start[i]/4)+1
  }  
}
##Takes 2+ hr on high performance power to run
##Test to see if it is done right
#How many new rows are supposed to be added on?
olddata<-dwdata[,c(1,5,6,7)]
mark<-rep(0,150000)
k=1
for (i in 1:length(olddata$year)) {
  if (floor(olddata$start[i]/4)!= floor(olddata$stop[i]/4))
  {mark[k]<-1
  k=k+1}
  else
  {k=k}
}
sum(mark)==length(newdata$year)-length(dwdata$year) #ok!

library(reshape)

allplants$cover<-allplants$stop-allplants$start
#take out the lines where cover=0 due to the cut-off code
#take out LITT and SOIL since we're not interested
allplants<-subset(allplants,allplants$cover>0)
bareIndex<-which(allplants$species_code=="LITT"|allplants$species_code=="SOIL"|allplants$species_code=="SOIL ")
allplants<-allplants[-bareIndex, ]
allplants[which(allplants$species_code=="BOER4 "),3]<-"BOER4"
#split spring and fall data (period= 1 vs. 3)
sp_alldata<-subset(allplants,allplants$period==1)#33293 obs
fl_alldata<-subset(allplants,allplants$period==3)#28119 obs


# Reshape spring and fall data -----------------------------------------------------
#Reshape spring data
mdata<-melt(sp_alldata[,c(1,2,3,6,7)], id.vars=c("year", "species_code","interval","period"))
sp_comp<-cast(mdata,interval+year~species_code+variable, fun.aggregate=sum)
for(i in 1:length(sp_comp$year)){
  sp_comp$tot_cover[i]<-sum(sp_comp[i,3:106])
}#takes 3 minutes to run on high performance
#104 spp in the spring

#Reshape fall data
mdata<-melt(fl_alldata[,c(1,2,3,6,7)], id.vars=c("year", "species_code","interval","period"))
fl_comp<-cast(mdata,interval+year~species_code+variable, fun.aggregate=sum)
for(i in 1:length(fl_comp$year)){
  fl_comp$tot_cover[i]<-sum(fl_comp[i,3:108])
} 
#107 spp in the fall

#Output the datasets
write.table(sp_comp,file="sp_comp.csv",sep=",",row.names=F)
write.table(fl_comp,file="fl_comp.csv",sep=",",row.names=F)


# Curate SPRING and FALL data  ---------------------
#SPRING
sp_comp <- read.csv("sp_comp.csv")
sp_bout<-sp_comp[,c(1,2,20,21,107)]#just the boutlouas and total
sp_bout$frac_BOER<-sp_bout$BOER4_cover/sp_bout$tot_cover 
sp_bout$frac_BOGR<-sp_bout$BOGR2_cover/sp_bout$tot_cover
#Only use sections with >35% mean cover
avg.tot.cover<-rep(0,100)
for (i in 1:100){
  temp<-subset(sp_bout,sp_bout$interval==i)
  avg.tot.cover[i]<-mean(temp$tot_cover)
}
hist(avg.tot.cover)#pretty nicely distributed
mean(avg.tot.cover)#1.65, maybe use 1.4 as cutoff (35%)
length(which(avg.tot.cover>1.4))#this gives 86 sections
covered.sections<-which(avg.tot.cover>1.4)
sp_bout_covered<-subset(sp_bout,sp_bout$interval %in% covered.sections)
#Only use sections where mean(BOER.frac+BOGR.frac)>50% 
avg.bout.cover<-cbind(covered.sections,rep(0,86))
for (i in 1:86){
  temp<-subset(sp_bout_covered,sp_bout_covered$interval==covered.sections[i])
  temp$bout<-temp$frac_BOER+temp$frac_BOGR
  avg.bout.cover[i,2]<-mean(temp$bout)
}
hist(avg.bout.cover[,2])
length(which(avg.bout.cover[,2]>0.5))#leaves 77 sections
bout.sections<-avg.bout.cover[which(avg.bout.cover[,2]>0.5),1]
sp_bout_good<-subset(sp_bout_covered,sp_bout_covered$interval %in% bout.sections )
#output dataset
write.table(sp_bout_good,file="sp_bout_good.csv",sep=",",row.names=F)

#FALL
fl_comp <- read.csv("fl_comp.csv")
fl_bout<-fl_comp[,c(1,2,20,21,109)]#just the boutlouas and total
fl_bout$frac_BOER<-fl_bout$BOER4_cover/fl_bout$tot_cover 
fl_bout$frac_BOGR<-fl_bout$BOGR2_cover/fl_bout$tot_cover
#Only use sections with >35% mean cover
avg.tot.cover<-rep(0,100)
for (i in 1:100){
  temp<-subset(fl_bout,fl_bout$interval==i)
  avg.tot.cover[i]<-mean(temp$tot_cover)
}
hist(avg.tot.cover)#pretty nicely distributed
mean(avg.tot.cover)#2.03, still use 1.4 as cutoff (35%)
length(which(avg.tot.cover>1.4))#this gives 100 sections (ALL sections have avg cover above 1.4)
covered.sections<-which(avg.tot.cover>1.4)
fl_bout_covered<-subset(fl_bout,fl_bout$interval %in% covered.sections)
#Only use sections where mean(BOER.frac+BOGR.frac)>50% 
avg.bout.cover<-cbind(covered.sections,rep(0,100))
for (i in 1:100){
  temp<-subset(fl_bout_covered,fl_bout_covered$interval==covered.sections[i])
  temp$bout<-temp$frac_BOER+temp$frac_BOGR
  avg.bout.cover[i,2]<-mean(temp$bout)
}
hist(avg.bout.cover[,2])
length(which(avg.bout.cover[,2]>0.5))#leaves 85 sections
bout.sections<-avg.bout.cover[which(avg.bout.cover[,2]>0.5),1]
fl_bout_good<-subset(fl_bout_covered,fl_bout_covered$interval %in% bout.sections )
#output dataset
write.table(fl_bout_good,file="fl_bout_good.csv",sep=",",row.names=F)


#Calculate Euclidean distances in spring data --------------------------------------------
#SPRING Euclidian distances between fractions of BOER and BOGR
interval<-bout.sections#make sure this is the spring version (n=77)
sum.dist.frac<-rep(0,77)
cv.dist.frac<-rep(0,77)
sum.dist<-rep(0,77)
cv.dist<-rep(0,77)
sum.dist.all<-rep(0,77)
cv.dist.all<-rep(0,77)
mean.all<-rep(0,77)
mean.x<-rep(0,77)#where x=%BOER-%BOGR "BOER-ness"
bc.dist.frac<-rep(0,77)#bray-curtis (city block) distance
bout.dist<-data.frame(interval,sum.dist,cv.dist,sum.dist.frac,cv.dist.frac,sum.dist.all,cv.dist.all,mean.all,mean.x,bc.dist.frac) #make new dataframe 
for (i in 1:77){
  subdata<-subset(sp_bout_good,sp_bout_good$interval==bout.dist$interval[i])
  temp.x<-rep(0,length(subdata$year))
  temp.frac<-rep(0,length(subdata$year))
  temp<-rep(0,length(subdata$year))
  temp.all<-rep(0,length(subdata$year))
  temp.bc<-rep(0,length(subdata$year))
  for(j in 1:(length(subdata$year)-1)){
    d.frac<-sqrt((subdata$frac_BOGR[j]-subdata$frac_BOGR[j+1])^2+(subdata$frac_BOER[j]-subdata$frac_BOER[j+1])^2)
    d<-sqrt((subdata$BOGR2_cover[j]-subdata$BOGR2_cover[j+1])^2+(subdata$BOER4_cover[j]-subdata$BOER4_cover[j+1])^2)
    d.all<-abs(subdata$tot_cover[j]-subdata$tot_cover[j+1])
    d.bc<-abs(subdata$BOGR2_cover[j]-subdata$BOGR2_cover[j+1])+abs(subdata$BOER4_cover[j]-subdata$BOER4_cover[j+1])
    temp.frac[j]<-d.frac
    temp[j]<-d
    temp.all[j]<-d.all
    temp.bc[j]<-d.bc
  }
  for (k in 1:length(subdata)){
    x<-subdata$frac_BOER[k]-subdata$frac_BOGR[k]
    temp.x[k]<-x
  }
  bout.dist$sum.dist[i]<-sum(temp)
  bout.dist$cv.dist[i]<-(sd(temp)/mean(temp))
  bout.dist$sum.dist.frac[i]<-sum(temp.frac)
  bout.dist$cv.dist.frac[i]<-(sd(temp.frac)/mean(temp.frac))
  bout.dist$sum.dist.all[i]<-sum(temp.all)
  bout.dist$cv.dist.all[i]<-(sd(temp.all)/mean(temp.all))
  bout.dist$mean.all[i]<-mean(subdata$tot_cover)
  bout.dist$mean.x[i]<-mean(temp.x)
  bout.dist$bc.dist.frac[i]<-sum(temp.bc)
}#loop for bout.dist
#Output bout.dist
write.table(bout.dist,file="bout.dist.csv",sep=",",row.names=F)


# Calculate dominance flips in spring data----------------------
#Load relevant data files
bout.dist <- read.csv("bout.dist.csv")
#Euclidean distances for spring curated data
sp_bout_good <- read.csv("sp_bout_good.csv")
#Cover data by year for spring curated data

sp_bout_good$delta<-sp_bout_good$frac_BOER-sp_bout_good$frac_BOGR

flips<-rep(0,77)#77 intervals in the good bout data
BOGR_to_BOER<-c()
BOER_to_BOGR<-c()
for(i in 1:77){
  subdata<-subset(sp_bout_good,sp_bout_good$interval==bout.dist$interval[i])
  temp<-cbind(rep(0,length(subdata$interval)-1),rep(0,length(subdata$interval)-1),rep(0,length(subdata$interval)-1))
  for (j in 1:(length(subdata$interval)-1)){
    if ((subdata$delta[j]*subdata$delta[j+1])<0 & abs(subdata$delta[j])>0.05 & abs(subdata$delta[j+1])>0.05) {
      temp[j,1]<-1
      if (subdata$delta[j]<0) temp[j,2]<-1
      else temp[j,3]<-1
    }
    #count as flip if there is switch in dominance AND the difference in percent cover in both cases is >5%
    else temp[j,1]<-0
  }
  flips[i]<-sum(temp[,1])
  BOGR_to_BOER[i]<-sum(temp[,2])
  BOER_to_BOGR[i]<-sum(temp[,3])
}

hist(flips)
hist(BOGR_to_BOER)
hist(BOER_to_BOGR)
bout.dist$flips<-flips
bout.dist$BOGR_to_BOER<-BOGR_to_BOER
bout.dist$BOER_to_BOGR<-BOER_to_BOGR

#Do the number of flip types differ?
mean(bout.dist$BOGR_to_BOER)#1.558442
quantile(bout.dist$BOGR_to_BOER)#0    1    1    2    7
mean(bout.dist$BOER_to_BOGR)#1.779221
quantile(bout.dist$BOER_to_BOGR)#0    1    2    3    8 
#BOER to BOGR slightly more frequent but CI's overlap

#A quick measure of the distribution of deltas for flip situations
index<-c()
for (j in 1:(length(sp_bout_good$delta)-1)){
  if ((sp_bout_good$delta[j]*sp_bout_good$delta[j+1])<0) {
    index<-c(index,j)  }
}
summary(sp_bout_good$delta[index])#mean is -0.01419
qt(0.975,df=length(index)-1)*sd(abs(sp_bout_good$delta[index]))/sqrt(length(index))#95%CI based on t distribution
#mean +/- 0.02408976

#How do #flips and sum.dist.frac compare as measures of dynamism?
plot(bout.dist$sum.dist.frac,bout.dist$flips)
#phew, there is still a positive relationsip!
#test
cor.test(bout.dist$sum.dist.frac,bout.dist$flips)#p=9.306e-11, r=0.656
#what happens when the big outlier point is removed?
sub<-subset(bout.dist,flips!=15)
cor.test(sub$sum.dist.frac,sub$flips)
#still significant relationship p=1.623e-07, r=0.558


#Combine dominance flips and euclidean dists as one metric-----------------
#Standardize
hist(bout.dist$flips)#poisson
hist(bout.dist$sum.dist.frac)#skewed normal
qqPlot(bout.dist$sum.dist.frac)#only slightly skewed at high values
#Transform both to standard normal distr
hist(log(bout.dist$flips+1))#+1 to everything won't make a huge difference and helps with distribution tp get rid of zeros
qqPlot(log(bout.dist$flips+1))#not great, but acceptable?
#convert to z scores (# of stdevs)
bout.dist$std.sumdistfrac<-(bout.dist$sum.dist.frac-mean(bout.dist$sum.dist.frac))/sd(bout.dist$sum.dist.frac)
hist(bout.dist$std.sumdistfrac)
bout.dist$std.flips<-(log(bout.dist$flips+1)-mean(log(bout.dist$flips+1)))/sd(log(bout.dist$flips+1))
hist(bout.dist$std.flips)#weird bimodal?
plot(bout.dist$std.flips,bout.dist$std.sumdistfrac)#still looks correlated
bout.dist$dyn.score<-bout.dist$std.sumdistfrac+bout.dist$std.flips#dynamic score
hist(bout.dist$dyn.score)#huh pretty decent

#output a new bout.dist file with new metrics
write.table(bout.dist,file="bout.dist.csv",sep=",",row.names=F)


#Candidate patch selection------------------------------
bout.dist <- read.csv("bout.dist.csv")
sp_bout_good <- read.csv("sp_bout_good.csv")
#Consider differences between euclid and bray curtis dist
plot(bout.dist$sum.dist.frac,bout.dist$bc.dist.frac)#pretty straight-forward positive relationship

#look at the distribution of "BOER-ness" (x)
hist(bout.dist$mean.x) #peaks at 0.95-0.1 (just a little more dominated by BOER)

#STATIC BOER CANDIDATES
BOER<-subset(bout.dist,bout.dist$mean.x>0)
#look at the distribution of euclid distances
hist(BOER$sum.dist.frac)
#want to select ones with high x and min sum.dist.frac
plot(BOER$sum.dist.frac,BOER$mean.x) #go for plots in top left quadrant (mean.x>0.1, and sum.dist.frac<5.5)
#candidate plots
index<-which(BOER$sum.dist.frac<5.5&BOER$mean.x>0.1)#20 candidates
BOER.cand<-BOER[index,]
#check to see if bc dist gets similar results
plot(BOER$bc.dist.frac,BOER$mean.x)#not quite the same
bc.index<-which(BOER$bc.dist.frac<15&BOER$mean.x>0.1)#gets 18
match(BOER$interval[index],BOER$interval[bc.index])#pretty much getting all the same ones

#STATIC BOGR CANDIDATES
BOGR<-subset(bout.dist,bout.dist$mean.x<0)
#look at the distribution of euclid distances
hist(BOGR$sum.dist.frac)
#want to select ones with low x and min sum.dist.frac
plot(BOGR$sum.dist.frac,BOGR$mean.x) #not many plots in low left quad
#candidate plots
index<-which(BOGR$sum.dist.frac<6&BOGR$mean.x<(-0.075))
BOGR.cand<-BOGR[index,]
#check to see if bc dist gets similar results
plot(BOGR$bc.dist.frac,BOGR$mean.x)#practically nothing in lower left quad
bc.index<-which(BOGR$bc.dist.frac<15&BOGR$mean.x<(-0.075))#gets 9
match(BOGR$interval[index],BOGR$interval[bc.index])#pretty much getting all the same ones


#DYNAMIC CANDIDATES
dynamic<-subset(bout.dist,bout.dist$mean.x>(-0.075)&bout.dist$mean.x<0.075)
#look at the distribution of euclid distances
hist(dynamic$sum.dist.frac)
#want ones with high dists
dynamic.cand<-subset(dynamic,dynamic$sum.dist.frac>5)#23 candidates
#These need to be hand-curated to check
#check to see if bc dist gets similar results
hist(dynamic$bc.dist.frac)
bc.index<-which(dynamic$bc.dist.frac>15)#gets 12
index<-which(dynamic$sum.dist.frac>5)
match(dynamic$interval[index],dynamic$interval[bc.index])#pretty much getting all the same ones

#Dynamic plot hand curation--------------------------------------------------------------
int19<-subset(sp_bout_good,sp_bout_good$interval==19)
plot(int19$year,int19$tot_cover,col="red",ylim=c(0, 4),main="plot 19 cover")
points(int19$year,int19$BOER4_cover,pch=16) 
points(int19$year,int19$BOGR2_cover,col="blue")
#good

int36<-subset(sp_bout_good,sp_bout_good$interval==36)
plot(int36$year,int36$tot_cover,col="red",ylim=c(0, 4),main="plot 36 cover")
points(int36$year,int36$BOER4_cover,pch=16) 
points(int36$year,int36$BOGR2_cover,col="blue")
#maybe

int40<-subset(sp_bout_good,sp_bout_good$interval==40)
plot(int40$year,int40$tot_cover,col="red",ylim=c(0, 4),main="plot 40 cover")
points(int40$year,int40$BOER4_cover,pch=16) 
points(int40$year,int40$BOGR2_cover,col="blue")

int45<-subset(sp_bout_good,sp_bout_good$interval==45)
plot(int45$year,int45$tot_cover,col="red",ylim=c(0, 4),main="plot 45 cover")
points(int45$year,int45$BOER4_cover,pch=16) 
points(int45$year,int45$BOGR2_cover,col="blue")
#nope

int47<-subset(sp_bout_good,sp_bout_good$interval==47)
plot(int47$year,int47$tot_cover,col="red",ylim=c(0, 4),main="plot 47 cover")
points(int47$year,int47$BOER4_cover,pch=16) 
points(int47$year,int47$BOGR2_cover,col="blue")
#mmaybe

int48<-subset(sp_bout_good,sp_bout_good$interval==48)
plot(int48$year,int48$tot_cover,col="red",ylim=c(0, 4),main="plot 48 cover")
points(int48$year,int48$BOER4_cover,pch=16) 
points(int48$year,int48$BOGR2_cover,col="blue")
#maybe

int<-subset(sp_bout_good,sp_bout_good$interval==52)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 52 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")
#good

int<-subset(sp_bout_good,sp_bout_good$interval==56)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 56 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")
#nope

int<-subset(sp_bout_good,sp_bout_good$interval==57)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 57 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")
#good

int<-subset(sp_bout_good,sp_bout_good$interval==58)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 58 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")
#good

int<-subset(sp_bout_good,sp_bout_good$interval==59)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 59 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")
#maybe

int<-subset(sp_bout_good,sp_bout_good$interval==60)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 60 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")
#maybe

int<-subset(sp_bout_good,sp_bout_good$interval==62)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 62 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")

int<-subset(sp_bout_good,sp_bout_good$interval==73)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 73 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")

int<-subset(sp_bout_good,sp_bout_good$interval==78)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 78 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")
#good

int<-subset(sp_bout_good,sp_bout_good$interval==86)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 86 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")
#maybe

int<-subset(sp_bout_good,sp_bout_good$interval==92)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 92 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")
#maybe

#sum.dist.frac 5-5.5 potentials
int<-subset(sp_bout_good,sp_bout_good$interval==40)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 40 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")
#maybe

int<-subset(sp_bout_good,sp_bout_good$interval==61)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 61 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")
#yes

int<-subset(sp_bout_good,sp_bout_good$interval==62)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 62 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")
#yes

int<-subset(sp_bout_good,sp_bout_good$interval==73)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 73 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")
#maybe

int<-subset(sp_bout_good,sp_bout_good$interval==94)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 94 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")
#maybe

int49<-subset(sp_bout_good,sp_bout_good$interval==49)
plot(int49$year,int49$tot_cover,col="red",ylim=c(0, 4),main="plot 49 cover")
points(int49$year,int49$BOER4_cover,pch=16) 
points(int49$year,int49$BOGR2_cover,col="blue")
#maybe

int67<-subset(sp_bout_good,sp_bout_good$interval==67)
plot(int67$year,int67$tot_cover,col="red",ylim=c(0, 4),main="plot 67 cover")
points(int67$year,int67$BOER4_cover,pch=16) 
points(int67$year,int67$BOGR2_cover,col="blue")
#maybe

int77<-subset(sp_bout_good,sp_bout_good$interval==77)
plot(int67$year,int77$tot_cover,col="red",ylim=c(0, 4),main="plot 77 cover")
points(int77$year,int77$BOER4_cover,pch=16) 
points(int77$year,int77$BOGR2_cover,col="blue")
#maybe

#STATIC CANDIDATES
int<-subset(sp_bout_good,sp_bout_good$interval==23)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 23 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")

int<-subset(sp_bout_good,sp_bout_good$interval==34)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 34 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")

int<-subset(sp_bout_good,sp_bout_good$interval==39)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 39 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")

int<-subset(sp_bout_good,sp_bout_good$interval==53)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 53 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")

int<-subset(sp_bout_good,sp_bout_good$interval==55)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 55 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")

int<-subset(sp_bout_good,sp_bout_good$interval==64)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 64 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")

int<-subset(sp_bout_good,sp_bout_good$interval==69)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 69 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")

int<-subset(sp_bout_good,sp_bout_good$interval==79)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 79 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")

int<-subset(sp_bout_good,sp_bout_good$interval==83)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 83 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")

int<-subset(sp_bout_good,sp_bout_good$interval==90)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 90 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")

int<-subset(sp_bout_good,sp_bout_good$interval==95)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 95 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")

int<-subset(sp_bout_good,sp_bout_good$interval==29)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 29 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")

int<-subset(sp_bout_good,sp_bout_good$interval==81)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 81 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")

int<-subset(sp_bout_good,sp_bout_good$interval==50)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 50 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")

int<-subset(sp_bout_good,sp_bout_good$interval==61)
plot(int$year,int$tot_cover,col="red",ylim=c(0, 4),main="plot 61 cover")
points(int$year,int$BOER4_cover,pch=16) 
points(int$year,int$BOGR2_cover,col="blue")