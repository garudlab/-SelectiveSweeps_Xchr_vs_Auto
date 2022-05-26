library(tidyr)
library(dplyr)

args = commandArgs(TRUE)

file = args[1]

df<-read.csv(file,header = F,sep = "\t", fill = TRUE)
df<- drop_na(df)
write.table(df,file,sep="\t",row.names=FALSE,col.names = F)

