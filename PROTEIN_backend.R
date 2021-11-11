#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#@Author:Samson Jacob
#Purpose: Protein and Peptide level analysis

#FUNCTIONS

sigvals = function(x,dec){
    return(round(as.numeric(x),2))
}

#normalized values
get_normedsub= function(peptidedf,proteindf,samplename){
  require(dplyr)
  filt_prot=subset(proteindf,proteindf$Potential.contaminant!='+' & proteindf$Reverse!='+')
  filt_prot = filt_prot[is.nan(filt_prot[[samplename]]) !=TRUE,]
  filt_pep=subset(peptidedf,peptidedf$Potential.contaminant!='+' & peptidedf$Reverse!='+') #drop reverse and contaminants (peptide)
  pep2 = subset(filt_pep,filt_pep$PEP<.0001 & filt_pep$Last.amino.acid %in% c('R','K')) # only choose peptides that end in R and K
  pep2 = pep2[is.nan(pep2[[samplename]]) !=TRUE,]#drop rows with missing values for that sample (drop missing values at the peptide level)
  pep3 = subset(pep2,pep2$Missed.cleavages== 0)
  normpept=select(pep3,contains("normalized")) 
  normprot = select(filt_prot,contains('normalized'))
  valnew = gsub('Ratio.H.L.','',samplename)
  pept1=select(normpept,contains(valnew))
  prot1=select(normprot,contains(valnew))
  pept1=as.data.frame(apply(pept1[1],1,function(x){100 * as.numeric(x)/sum(as.numeric(x),1)}))
  colnames(pept1)=paste('NORMpept',samplename,sep='_')
  prot1=as.data.frame(apply(prot1[1],1,function(x){100* as.numeric(x)/sum(as.numeric(x),1)}))
  colnames(prot1)=paste('NORMprot',samplename,sep='_')
  return(list(prot1,pept1)) #protein first peptide second
}


get_bestPeps=function(peptidedf,protdf,sample){
  colimp=c(sample,'Protein.group.IDs') # columns to keep
  filteredproteins = subset(protdf,protdf$Potential.contaminant!='+' & protdf$Reverse!='+') # drop reverse and contaminants(protein)
  filt_pep=subset(peptidedf,peptidedf$Potential.contaminant!='+' & peptidedf$Reverse!='+') #drop reverse and contaminants (peptide)
  pep2 = subset(filt_pep,filt_pep$PEP<.0001 & filt_pep$Last.amino.acid %in% c('R','K')) # only choose peptides that end in R and K
  pep2 = pep2[is.nan(pep2[[sample]]) !=TRUE,]#drop rows with missing values for that sample (drop missing values at the peptide level)
  pep3 = subset(pep2,pep2$Missed.cleavages== 0) # chose peptides with 0 missed cleavages
  pep3=as.data.frame(pep3[, (names(pep3) %in% colimp)]) # subset the columns of interest
  s = strsplit(as.character(pep3$Protein.group.IDs), split = ";")# explode the protein IDs of the best peptides(peptide level)
  out=data.frame(V1 = rep(pep3[[sample]],sapply(s,length)),V2=unlist(s)) # create a dummy dataframe to store the exploded values; retaining mapping
  out$PEPscore = apply(out[1],1,function(x){100 * as.numeric(x)/sum(as.numeric(x),1)})
  colnames(out)=c(sample,'Protein.group.IDs','Score') #correct headers
  prt = subset(filteredproteins,filteredproteins$id %in% out$Protein.group.IDs) # get the Protein Groups from exploded list of Protein_Group IDs
  exactsample = as.data.frame(prt[[sample]]) # get the Ratio H/L values from the Protein Level for that sample (after Proteins are filtered above)
  prt$Score=apply(exactsample[1],1,function(x){100 * as.numeric(x)/sum(as.numeric(x),1)}) # apply the score
  return(list(prt,out))
}

getprotgrps_pepstats=function(output,sampname){
  #from the output of the peptide/proteing groups analysis return descriptive statistic values
  best_protgroups= as.data.frame(output[1]) #proteinlevel
  best_peps = as.data.frame(output[2]) #peptide level
  orig= c("BEST_Peptide","QUERIED_PROTEIN_GROUPS")
  colls = c('Sample','Mean','Median','Variance','SD','Column')
  #outfr= data.frame(row.names=c('Best_Proteins','Best_Peptides'))
  pepsum = summary(best_peps$Score)
  protsum = summary(best_protgroups$Score)
  medvals= sigvals(c(pepsum[3],protsum[3]),3)
  meanvals=sigvals(c(pepsum[4],protsum[4]),3)
  sdvals = sigvals(c(sd(best_peps$Score,na.rm = TRUE),sd(best_protgroups$Score,na.rm=TRUE)),3)
  varvals = sigvals(c(var(best_peps$Score,na.rm=TRUE),sd(best_protgroups$Score,na.rm=TRUE)),3)
  out=as.data.frame(cbind(meanvals,medvals,sdvals,varvals,orig))
  out$Sample = rep(sampname,nrow(out))
  row.names(out)=seq(1,nrow(out))
  return(out)
  oufile = paste(sampname,'ProteinLevel.txt',sep='_')
  write.table(out,file=oufile,sep='\t')
}

plot_protgrps_pep=function(output,sampname){
  require(data.table)
  require(ggplot2)
  require(grid)
  require(gridExtra)
  picname=paste(sampname,'_DBplot.png',sep='')
  themeSJ=ttheme_default(
    core = list(fg_params=list(cex = .6)),
    colhead = list(fg_params=list(cex = .6)),
    rowhead = list(fg_params=list(cex = .6)))
  restab=getprotgrps_pepstats(output,sampname)
  grb=tableGrob(restab,rows=NULL,theme = themeSJ)
  best_protgroups= as.data.frame(output[1])
  best_peps = as.data.frame(output[2])
  f = as.data.frame(best_peps[['Score']])
  f$Origin = rep('Peptide',nrow(f))
  g = as.data.frame(best_protgroups[['Score']])
  g$Origin = rep('Protein',nrow(g))
  tog=rbindlist(list(g,f))
  colnames(tog)=c('Score','Origin')
  grph1=ggplot(tog, aes(x=tog$Score, fill=tog$Origin)) + geom_histogram(alpha=0.35, position="identity",bins=200)+ggtitle('Protein_Level and Peptide')+
  xlab(sampname)+guides(fill=guide_legend(title=NULL))+theme(axis.title.x = element_text(colour = "red"))
  plt=grid.arrange(grph1,grb,nrow=2,heights=c(3,1),as.table=TRUE)
  ggsave(picname,plot=plt,width = 20, height = 20, units = "cm")
}

#MAIN

### COMMANDLINE
if (length(args)>1) {
  pept = read.table(args[1],sep='\t',header=TRUE)
  prot = read.table(args[2],sep='\t',header=TRUE)
  df = read.table(args[3], header=TRUE,sep='\t')
  samps = df$samples
  tim= format(Sys.time(), "%X%b%d%Y")
  s=gsub(":","_",tim)
  mainDir <- getwd()
  subDir <- paste("PROTEIN_RES",s,sep='_')
  dir.create(file.path(mainDir,subDir),showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))
  for(c in samps){
    ou= get_bestPeps(pept,prot,c)
    plot_protgrps_pep(ou,c)
  }
}

