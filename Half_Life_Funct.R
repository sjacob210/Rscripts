#!/Users/bin/env

#@Author: Samson Jacob 1.2015

calcHalflife <- function(x,tcc) {
  newdf = data.frame(row.names=row.names(x))##sets the dimensions for the newdf output
  for(i in 1:nrow(x)){ ##lnR+1 ratio
    theseratios = as.numeric(x[i,2:6]) ##correct the index 
    time = as.numeric(c(15,30,60,120,180)) ##correct the timepoints
    df = as.data.frame(t(rbind(time,theseratios)))
    fit = lm(theseratios~time-1,data=df)## generate linear model of the data
    newdf$Coef[i] = coef(fit) # attain the coefficient
    newdf$kDEG[i]=  coef(fit)-(log(2,base=exp(1)/tcc)) # calculate degradation
    newdf$HalfLifeMin[i]=  log(2,base=exp(1))/newdf$kDEG[i] # calc half life
    newdf$R2[i]= summary(fit)$r.squared # view the fit
  }
  return(newdf)
}
