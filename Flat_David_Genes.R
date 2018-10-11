#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# usage for this script.
# Arg1 Path of the results from the DAVID web 
# Arg2: A value to supply for coding the results 

#ie:
# Rscript --vanilla Flat_David_Genes.R ./Down_DavidResults.txt Down


#@Author: Samson Jacob 1.14.2018
#Purpose: Take the results of the David Outfile and flatten the genes associated with an ontology
#Source for the Coded file is from Cronstein2.R

######################################################################### Attain the Up and Down regulated Genes

#find the coded output
#file1 = list.files(".",pattern = "DESEQ_CODED_mouse.txt")

#load the coded file 
#df_coded = read.table(file1[1],sep='\t',header=TRUE)

#UP regulated genes
#Up = subset(df_coded,df_coded$VAL=='UP')
#Dw = subset(df_coded,df_coded$VAL=='DOWN')

#write_out
#write.table(Dw,'Down_Mouse.txt',sep = '\t',row.names = TRUE,quote = FALSE)
#write.table(Up,'UP_Mouse.txt',sep='\t',row.names=TRUE,quote = FALSE)

#rm(list=ls())

############################################################################# Explode DAVID Output

#load the path of the results from DAVID analysis
#downDAVIDpth = './Downloads/Down_DavidResults.txt'
#downKEGGpth = './Downloads/KEGG_DOWN_mouse.txt'

#column headings to keep 
keepcols = c('Term','Count','PValue','Genes','FDR')

#function attains the unique biological processes
get_BPcats = function(f){
  hh = read.table(f,sep='\t',header=TRUE)
  hh = hh[,keepcols]
  bb = data.frame(do.call('rbind', strsplit(as.character(hh$Term),'~',fixed=TRUE))) #splt the Term column into 2 columns
  colnames(bb)=c('GO_ID','Biological_Process')
  return(as.character(unique(bb$Biological_Process)))
}


#function to convert output from DAVID to explode genes
flat_BPdavid_out = function(df1,subset1){
  require(splitstackshape)
  require(org.Mm.eg.db) 
  h = df1[,(names(df1)%in% keepcols)]
  ext = cSplit(h,'Genes',',','long') # perform the explosion on the column Genes
  bb = data.frame(do.call('rbind', strsplit(as.character(ext$Term),'~',fixed=TRUE))) #split the Term column into 2 columns
  bb$genes = ext$Genes
  colnames(bb)=c('GO_ID','Biological_Process','Genes')
  sub = subset(bb,bb$Biological_Process==subset1) #subset by biological process name
  names1 = as.character(sub$Genes) # convert ENTREZID into a character
  eg = AnnotationDbi::select(org.Mm.eg.db,keys=names1,columns='SYMBOL',keytype='ENTREZID') #convert from ENTREZID to GENESYMBOL
  return(merge(eg,sub,by.x='ENTREZID',by.y='Genes')) #merge the converted names with the subset
}


#read in the file
#df = read.table(downDAVIDpth,sep='\t',header=TRUE)

#attain all of the unique Biological Processes
#process = get_cats(downDAVIDpth)

# #loop through each of the BPs to generate a converted file for each
# for(c in process){
#   m=gsub(" ", "_", c, fixed = TRUE)#remove white space
#   oufile = paste(m,'.txt',sep='')
#   nu=flat_BPdavid_out(df,c)
#   write.table(nu,file=oufile,sep='\t',row.names = FALSE)
# }

#main function
conv_explode = function(pth,val){
  df = read.table(pth,sep='\t',header=TRUE) # read in this file (DAVID output)
  process = get_BPcats(pth) #acquire all of the unique
  for(c in process){ # loop through all the biological processes
    m=gsub(" ", "_", c, fixed = TRUE)#remove white space
    oufile = paste(m,'.txt',sep='') # add the text suffix
    oufile = paste(val,oufile,sep='') # add the value prefix ie up/down for organization
    nu=flat_BPdavid_out(df,c) # run the function with the subset c
    write.table(nu,file=oufile,sep='\t',row.names = FALSE,quote = FALSE) # write out the subset 
  }
}

# commandline arguments 

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("2 arguments must be supplied ", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  stop("Please provide both the path to the file and a Value(up/down)",call.=FALSE)
}else if (length(args)==2){
  conv_explode(args[1],args[2])
}


