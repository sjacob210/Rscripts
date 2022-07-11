#require(devtools)
#install_version("arsenal", version = "3.6.3", repos = "http://cran.us.r-project.org")

library(arsenal)
library(pheatmap)
library(tidyverse)


find_col<-function(df,value){
  # Find the column index for a character value
  if(value%in%names(df)){
    return(grep(value, colnames(df)))
  }else{
    print("Column Not found")
  }
}

return_rows<-function(df,vecofixes,column){
  return(df[vecofixes,]%>%
           select(column)) # non-specific
}

com_df<-function(prevdf,curdf){
  # find differences in curdf
     require(arsenal)
  # run compardf function, return differences for only the curdf
     curdf_fail<-select(as.data.frame(summary(comparedf(x=prevdf,y=curdf))["diffs.table"]),contains("y"))
     colnames(curdf_fail)<-c("COLUMN","DIFFERENT","ROW") # rename columns
     return(curdf_fail)
   }

create_df_zeros<-function(df){
  # create a dataframe of zeros from an input df
  cols<-as.numeric(length(colnames(df)))
  rws<-as.numeric(nrow(df))
  mat = data.frame(matrix(0, rws, cols))
  colnames(mat)<-as.character(colnames(df))
  return(mat)
}

# take results of com_df and update the values to 1 on the

res_compdf<-com_df(df1,df2)

s



subj_id<-rep(c("AA01","AA02","AA03"),2)
timepoint<-rep(c("SCREEN",'C1D1',"C1D1"),2)
timepointr<-rep(c("C1D1","C1D2","SCREEN"),2)

df1<-data.frame(SUBJID=subj_id,TIMEPOINT=timepoint)
df2<-data.frame(SUBJID=subj_id,TIMEPOINT=timepointr)

(res<-as.data.frame(summary(comparedf(x=df1,y=df2))["diffs.table"]))



find_col(df1,"SUBJIDz")


