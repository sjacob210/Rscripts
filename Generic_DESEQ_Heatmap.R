#!/Users/bin/env/RScript

args <- commandArgs(trailingOnly=TRUE)

# install packages if not installed
required_base_Packages = c('dplyr','ggplot2','gdata','reshape2','pheatmap','grid')
for(p in required_base_Packages){
  if(!require(p,character.only = TRUE)) install.packages(p)
}

# load libraries after install
library(ggplot2)
library(dplyr)
library(gdata)
library(reshape2)
library(pheatmap)
library(grid)

#Generic Heatmap for DESEQ output as .xlsx file
plot_Heatmap <-function(xlsxfile,colname){
  df <- read.xls (xlsxfile,sheet=1,header=TRUE)
  df2 <- subset(df,df$padj<.05)
  df2 <-df2[c(1,3,6,7,10:17)]
  df2$RowSums<- apply(df2[,5:12],1,function(x){sum(as.numeric(x))})
  logoutdf <- data.frame(row.names = 1:length(df2[[colname]]))
  logoutdf[,1:8] <- log2(df2[,5:12]/df2[,13]) # normalized to the RowSums hardwired
  logoutdf[[colname]]<-df2[[colname]]
  vals<-colorRampPalette(c('#2e287a','#453f92','navy','#f2d054','yellow'))(100)
  purplevals <- colorRamp(c('#3f0854','#5B0C79','#FBC33A'))(100)
  outmap<-pheatmap(as.matrix(logoutdf[,1:8]),
                   col =vals,
                   distfun=function(x) as.dist((1-cor(t(x)))/2),
                   treeheight_row=0,
                   show_rownames = F,
                   main='Comparison\n')
  pdf('Heatmap.pdf',width=7,height=7)
  grid.newpage()
  grid.draw(outmap$gtable)
  dev.off()
  }

  ## comandline args
if (length(args) > 0) {
  xlsxfile = args[1]
  tim<- format(Sys.time(), "%X%b%d%Y") # get the current time
  s<-gsub(":","_",tim) # replace : with _ for the current time
  mainDir <- getwd() # get the current working directory
  subDir <- paste("Heatmap_Output",s,sep='_') # create the name of the subdir called HeatMap output with the time
  dir.create(file.path(mainDir,subDir),showWarnings = FALSE) # create the subdir in the current working directory
  setwd(file.path(mainDir, subDir)) # set the working directory as the newly created subdir
  plot_Heatmap(xlsxfile)
}