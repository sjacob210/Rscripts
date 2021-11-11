#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Purpose: Conver the distributions for a specific sample into 50% cwentering

library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(data.table)
library(plyr)
library(dplyr)



##FUNCTIONS
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

#for raw
scale_to_spec = function(df,dfcol,val){
  comps = df[is.nan(df[[dfcol]]) !=TRUE,] # drop missing
  s = as.numeric(comps[[dfcol]]) # take the raw H/L value
  me = median(s) # calculate the median
  mx = max(s)# calculate the max
  mi= min(s)# calculate the minimum
  diff1 = mx-val # current maximum - desired val(median)
  diff2 = val-mi #  desired val(same) -current minimum
  subr = as.data.frame(s-me) # subtract all values by the median
  below = subset(subr,s<0) # all negative values
  above = subset(subr,s>0)# all positive values
  y = ((above/mx) *diff1) # divide all positive values from current max; multiply this by difference of desired max  desired median
  x = ((below/abs(mi))*diff2) # divide all negative values from abs-current min; multiply by difference of desired median and curren min
  outf = as.data.frame(rbind(y,x))
  outf= outf+val # add the desired median
  outf$Score=apply(outf[1],1,function(x){100 * as.numeric(x)/sum(as.numeric(x),1)}) # apply the score
  outf$Original = s
  outf=outf[,c(3,1,2)]
  colnames(outf)=c(paste('RAW',dfcol,sep='_'),paste('transformed_RAW',dfcol,sep='_'),'Transformed_Score')
  return(outf)
}

fix_genename = function(df1col){
  # function takes a dataframe column (df1col) as input and corrects the column name to grep GENENAME only #
  m = as.character(df1col) #no reason for 2 steps
  return(regmatches(m,regexpr("([0-9,A-Z])\\w+_HUMAN",m,perl=TRUE))) #apply regex to get the Uniprot ID only (careful not to use with columns with multiple IDs)
}

sigvals = function(x,dec){
    return(round(as.numeric(x),2))
}

calc_sigma=function(dfcol){
  x = as.numeric(dfcol)
  me = sigvals(median(x,na.rm=TRUE),2)
  s = sigvals(sd(x,na.rm=TRUE),2)
  return(c(me,s))
}

scale_to_score = function(df,dfcol,val){
  ## val is designed to be 50 for score
  comps = df[is.nan(df[[dfcol]]) !=TRUE,] # drop missing
  nims = fix_genename(comps$Protein.IDs)
  s = as.numeric(comps$Score) # take the score
  me = median(s) # calculate the median of the score
  x = val/me
  outf= data.frame(row.names=seq(1:length(s)))
  outf$Proteins = nims
  outf$Score = s
  outf$Converted =  outf$Score * x
  colnames(outf)=c('Proteins',paste('Score',dfcol,sep='_'),paste('transformed_SCORE',dfcol,sep='_'))
  return(outf)
}

graph_vals=function(outdfcol){
  vals= calc_sigma(outdfcol) # get the sd and the median
  medval = as.numeric(vals[1]) # median
  sig = as.numeric(vals[2]) # sd
  twosig = sigvals((2 * sig)+medval,2) # 2 sigmas above the median
  threesig= sigvals((3 * sig)+medval,2) # 3 sigmas
  hh = c(medval,sig,twosig,threesig)
  return(hh)
}

plot_dens=function(datafrcol,lim1,lim2,titlename){
  lbl1=as.character(expression(paste('2',sigma)))
  lbl2=as.character(expression(paste('3',sigma)))
  vv=as.data.frame(datafrcol)
  vv$x = seq(1:length(nrow(vv)))
  colnames(vv)=c('x','y')
  dat = with(density(vv$x),data.frame(x,y))
  a=c(0,25,50,lim1,lim2,max(dat$x))
  brks=sprintf("%0.2f", round(a, digits = 2))
  p=ggplot(data=dat,mapping=aes(x=x,y=y))+
  geom_line(color='black',size=1, alpha=0.4)+
  geom_area(mapping= aes(x= ifelse(x<50,x,0)),fill='grey')+
  geom_area(mapping = aes(x = ifelse(x>50 & x< lim1 , x, 0)), fill = "#5079ff")+
  geom_area(mapping =  aes(x=  ifelse(x>lim1 & x<lim2, x,0)),fill='#79ff50')+ 
  geom_area(mapping = aes(x= ifelse(x>lim2,x,0)),fill='#ff033e')+
    xlim(0,100)+ylim(0,.12)+xlab('Transformed_Score')+
    ylab('')+ggtitle(titlename)+
    geom_vline(xintercept=c(50,lim1,lim2),linetype='dotted',color='black')+
    scale_x_continuous(breaks=as.numeric(brks))+
    annotate('text', x = c(lim1,lim2), y = c(.110,.110), label = c(lbl1,lbl2),parse = TRUE,size=5)
  return(p+theme_bw() + theme(panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))
}


#both distributions
#outputcol = output from scale toscore
plot_both = function(outputcol,val){
  require(reshape2)
  require(ggplot2)  
  mdf = melt(outputcol[,2:3])
  p=ggplot(mdf,aes(x=value,fill=variable))+geom_density(alpha=0.1)+
    ggtitle(val)+
    theme_bw() + theme(legend.position="bottom",panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  return(p)
}

# plot_both(ou$transformed_SCORE_Ratio.H.L.9,'test')
# mdf = melt(ou[,2:3])

#TEST
# plot_dens(ou$transformed_SCORE_Ratio.H.L.1,v[3],v[4],'Ratio_HL_1')
if (length(args)>1) {
  pept = read.table(args[1],sep='\t',header=TRUE)
  prot = read.table(args[2],sep='\t',header=TRUE)
  df = read.table(args[3], header=TRUE,sep='\t')
  samples = df$Samples
  tim= format(Sys.time(), "%X%b%d%Y")
  s=gsub(":","_",tim)
  mainDir <- getwd()
  subDir <- paste("CENTERED_RES",s,sep='_')
  dir.create(file.path(mainDir,subDir),showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))
  rws=c()
  pltlist=list()
  combolist=list()
  for(x in 1:length(samples)){
    nm= as.character(samples[x]) # name of the sample
    oufile= paste(nm,'Centered.txt',sep='_') # ouput
    clname = paste('transformed_SCORE',nm,sep='_') #transformed score name
    output = get_bestPeps(pept,prot,nm)# test Sample1
    outputdf = as.data.frame(output[1]) # sample1df
    ou = scale_to_score(outputdf,nm,50)
    hh=plot_both(ou,nm)
    v = graph_vals(ou[[clname]]) # get median,sd, 2sigma,3sigma
    ou$Significant = ifelse(ou[[clname]]>v[3],1,0)
    rws=append(rws,as.numeric(nrow(ou)))
    write.table(ou,file = oufile,sep='\t',quote=FALSE,row.names=FALSE)
    p=plot_dens(ou[3],v[3],v[4],nm)
    pltlist[[x]]=p
    combolist[[x]]=hh
  }
}

ml <- marrangeGrob(pltlist, nrow=3, ncol=1)
ggsave('Densities_Tranformed.pdf',ml)
ml2 = marrangeGrob(combolist,nrow=3,ncol=1)
ggsave('Transormations_dens.pdf',ml2)

# calculate totals
aggs = list.files(pattern='Centered.txt') # get all the outputs
y <- lapply(aggs, read.table, header=TRUE,sep='\t',colClasses=c(NA,'NULL','NULL',NA))

## HACK
test= do.call(rbind.fill, y)
#rws=c(372,360,198,234,207,178,267,153,238,396,248,517)

l=c()
for(s in 1:length(aggs)){
  file=as.character(aggs[s])
  reps = as.numeric(rws[s])
  l=append(l,rep(file,reps))
}
f = data.frame(SampleName=cbind(l))
bb = cbind(test,f)
colnames(bb)=c('Protein','Significant','Sample')

xx=aggregate(bb$Significant, by=list(PROTEIN=bb$Protein), FUN=sum)
colnames(xx)=c('PROTEIN_NAME','Aggregated_Counts')
write.table(xx,file='AggregatedCounts.txt',sep='\t',quote=FALSE,row.names=FALSE)
write.table(bb,file='SigmaResults.txt',sep='\t',quote = FALSE,row.names=FALSE)
yy = aggregate(bb$Significant,by=list(PROTEIN=bb$Protein,SAMP=bb$Sample),FUN=sum)
y2 = subset(yy,yy$x>0)

#sigs=subset(xx,xx$Aggregated_Counts>0)
#barplot(table(sigs$Aggregated_Counts))


