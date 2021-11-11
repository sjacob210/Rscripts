#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#@Author: Samson Jacob
#Purpose: Peptide LEVEL analysis COMPLETE

#write out info
##writeLines(capture.output(sessionInfo()), "SAMSON_sessionInfo.txt")

#Functions
sigvals = function(x,dec){
    return(round(as.numeric(x),2))
}

#helper function to break dataframe columns
get_subrowsfromdataframe=function(df,val){
  subfr = subset(df,df$Sample==val)
  vals = as.numeric(unlist(subfr[1]))
}

tblgrob = function(output){
  #grob function for plotting with a ggplot figure
  require(gridExtra)
  themeSJ= gridExtra::ttheme_default(
    core = list(fg_params=list(cex = .7)),
    colhead = list(fg_params=list(cex = .7)),
    rowhead = list(fg_params=list(cex = .7)))
  v=as.data.frame(output[1])
  return(tableGrob(v,rows=NULL,theme = themeSJ))
}

PEPT_Incorporation_frame = function(df,dfcolumn){
  ## Read in a maxquant peptide file generate a list of 2 dataframes
  ## ARG:: df = dataframe
  ## ARG:: dfcolumn = string representing sample of interest
  require(dplyr) # full join function
  rowvals = c('Arg_peptides','Lys_peptides','FILTERED')
  colvals = c('Mean','Median','Standard_Deviation','Sample_Name','Subpopulation','Peptides{ROWS}')
  val = as.character(dfcolumn)
  #filteredpeptides = filter(df, othercol == 'Median') #make sure the other column is also the median
  filterpeptides=filter(df,df$Ratio.H.L=='Median')
  filteredpeptides = df[is.na(df[[dfcolumn]]) !=TRUE,] # drop the missing values from the particular sample
  filteredpeptides = subset(filteredpeptides,filteredpeptides$Potential.contaminant!='*' & filteredpeptides$Reverse!='*' & filteredpeptides$Missed.cleavages==0) #additional filters # drop missed cleavages
  filteredpeptides2 = subset(filteredpeptides,filteredpeptides$PEP<.0001) # add posterior error probability filter < .0001 to larger population
  filteredpeptides2 = subset(filteredpeptides2,filteredpeptides2$Last.amino.acid%in%c('R','K')) # filter on the Posterior Error Probability
  Rpeptides = filteredpeptides2[filteredpeptides2$R.Count>0,] # validate missed cleavage, have the opposite AA in the cterminus (R in Missed cleavages)
  Rpeptides = Rpeptides[Rpeptides$K.Count==0,]
  Rpeptides = subset(Rpeptides,Rpeptides$Last.amino.acid=='R') #ends in R
  Kpeptides = filteredpeptides2[filteredpeptides2$K.Count>0,] # no more than 1 K
  Kpeptides = Kpeptides[Kpeptides$R.Count==0,]
  Kpeptides = subset(Kpeptides,Kpeptides$Last.amino.acid=='K') # ends in K
  rvals= row.names(Rpeptides)
  kvals = row.names(Kpeptides)
  incl = append(rvals,kvals)
  ALL_sample = as.data.frame(100*filteredpeptides2[[dfcolumn]]/(filteredpeptides2[[dfcolumn]]+1)) # convert to score
  ALL_sample$Sample = rep('Filtered',nrow(ALL_sample))
  K = as.data.frame(100 * Kpeptides[[dfcolumn]]/(Kpeptides[[dfcolumn]]+1)) # K
  K$Sample =rep('Lys',nrow(K))
  R = as.data.frame(100 * Rpeptides[[dfcolumn]]/(Rpeptides[[dfcolumn]]+1)) # R
  R$Sample = rep('Arg',nrow(R))
  alls= full_join(R,K)
  alls = as.data.frame(full_join(alls,ALL_sample))
  colnames(alls)=c('Arginine_Score','Sample','Lysine_Score','Filtered_Score')
  ag=as.data.frame(aggregate(alls,by=list(alls$Sample),FUN = median,na.rm=TRUE))
  outfilename= paste(val,'agginfo.txt',sep='_')
  write.table(ag,file=outfilename,sep='\t')
  KVALS=get_subrowsfromdataframe(K,'Lys')
  RVALS = get_subrowsfromdataframe(R,'Arg')
  FILVALS = get_subrowsfromdataframe(ALL_sample,'Filtered')
  filstats= summary(FILVALS)
  Rstats = summary(RVALS)
  Kstats = summary(KVALS)
  bound_sd = sigvals(c(sd(RVALS),sd(KVALS),sd(FILVALS)),3)
  bound_length=sigvals(c(length(RVALS),length(KVALS),length(FILVALS)),3)
  bound_mean = sigvals(c(as.numeric(Rstats[4]),as.numeric(Kstats[4]),as.numeric(filstats[4])),3)
  bound_median=sigvals(c(as.numeric(Rstats[3]),as.numeric(Kstats[3]),as.numeric(filstats[3])),3)
  sampletypcol = rep(dfcolumn,3)
  outframe = as.data.frame(cbind(bound_mean,bound_median,bound_sd,sampletypcol,rowvals,bound_length))
  colnames(outframe)=colvals
  outfile=paste(dfcolumn,'_PEPTIDE_Subpopulation.txt')
  write.table(outframe,file = outfile,sep='\t',quote =FALSE)
  return(list(outframe,alls))
}

bplot_subs=function(output,nam){
  # barplot of the subpopulations
  ##::ARG::output= list of dataframe from PEPT_Incorporation_frame(function) result
  outdf = as.data.frame(output[2])
  tdf = as.data.frame(table(outdf$Sample))
  vec = as.numeric(tdf[,2])
  ylim <- c(0, 1.1*max(vec))
  val= as.character(nam)
  pdfname=paste(val,'sub.pdf',sep='_')
  pdf(pdfname)
  xbar=barplot(table(outdf$Sample),col='navy',main=val)
  text(x=xbar,y=vec-100,labels=as.character(tdf[,2]),col='orange')
  dev.off()
}

plot_Scores=function(output,title){
  #returns a graph of all the distributions with a table of stats
  #::ARG::output:: a list from PEPT_Incorporation_from(function)
  require(grid)
  require(gridExtra)
  require(ggplot2)
  require(reshape2)
  picname =paste(title,'subpopulationHIST.png',sep='_')
  b = tblgrob(output)
  m = melt(output[2])
  graph2= ggplot(m, aes(x=value, fill=variable)) + geom_histogram(alpha=0.35, position="identity",bins=200)+ggtitle(title)+theme(legend.title=element_blank())
  plt=grid.arrange(graph2,b,nrow=2,heights=c(4,1),as.table=TRUE)
  ggsave(picname,plot=plt,width = 20, height = 20, units = "cm")
}

# plotpep = function(testf,val){
#   ##plot the PEP distribution
#   filtdf = subset(testf,testf$Reverse!='+' & testf$Potential.contaminant!='+')
#   h <- hist(testf$PEP,breaks=10000,plot=FALSE)
#   axis(side=1,at=seq(0,0.01,.001))
#   cuts <- cut(h$breaks, c(-Inf,.0001))
#   pdf("Export_Plots.pdf", width = 16 , height = 10, title = "EDA Plots for data")
#   for(i in 1:2){
#   plot(h, col=cuts,xlim=c(0,0.0005),main='PEP values:: UNFILTERED PEPTIDES \n range(0-0.04)',xlab='black is < .0001)')
#   plot(h2,col=cuts,xlim=c(0,0.0005),main='PEP values:: FILTERED REV/CONT \n range(0-0.04)',xlab='black is <.0001')
#   }
#   dev.off()
#   }
# 
# plotpep1= function(testf,val){
#   pname = paste(val,'Unfiltered.pdf',sep='_')
#   pdf(pname)
#   h = hist(testf$PEP,breaks=10000,plot=FALSE)
#   axis(side=1,at=seq(0,0.01,.001))
#   cuts =cut(h$breaks, c(-Inf,.0001))
#   plot(h)
#   dev.off()
# }


# mainfunct= function(peptidefile,samplefile){
#   df = read.csv(samplefile,header=TRUE,sep='\t')
#   samps = as.vector(df[,1])
#   pept = read.table(peptidefile,sep='\t',header=TRUE)
#   for(s in samps){
#     h = PEPT_Incorporation_frame(pept,s)
#     plot_Scores(h,s)
#     bplot_subs(h,s)
#  }
# }


##COMMANDLINE
#commandline
if (length(args)>0) {
  pept = read.table(args[1], header=TRUE,sep='\t')
  df = read.table(args[2], header=TRUE,sep='\t')
  tim= format(Sys.time(), "%X%b%d%Y")
  s=gsub(":","_",tim)
  mainDir <- getwd()
  subDir <- paste("PEP_RES",s,sep='_')
  dir.create(file.path(mainDir,subDir),showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))
  samps = df$samples
  for(s in samps){
    h = PEPT_Incorporation_frame(pept,s)
    plot_Scores(h,s)
    bplot_subs(h,s)
  }
}

