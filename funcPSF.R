#Function to calculate PSF ln ratio
func.PSF<-function(data){
  #Make new dataframe to store
  PSF<-data.frame(plot.census=character(),
                  Plot=integer(), 
                  Census=integer(),
                  State=character(), 
                  Dom_spp=character(),
                  Block=integer(),
                  lnRatio=double(),
                  stringsAsFactors=FALSE) 
  for (i in 1: length(unique(data$plot.census))){
    subdata<-subset(data,plot.census==unique(data$plot.census)[i])
    if (nrow(subdata)!=2) next #don't bother calculating PSF if no living pair
    else {
      options(stringsAsFactors=FALSE)
      PSF<-rbind(PSF, list(lnRatio = subdata$logMass.harv[which(subdata$PSF=="Home")]-subdata$logMass.harv[which(subdata$PSF=="Away")],
                           plot.census=subdata$plot.census[1],
                           Plot=subdata$Plot[1],
                           Census=subdata$Census[1],
                           Block=subdata$Block[1],
                           State=subdata$State[1],
                           Dom_spp=subdata$Dom_SPP[1]))
    }
  }
  return(PSF)
}