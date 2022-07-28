
rm(list=ls())

library(tidyr)
source("Indexes_formulas.r")
#########################Import data##########################
taxonomy<- read.csv("Taxonomy_18_01_2022.csv", header=T, sep=";", dec=".",na.strings = "NA")
samples<- read.csv("Bio_lit_Samples.csv", header=T, sep=";", dec=".")
macro_data<- read.csv("Bio_lit_UPDATED.csv", header=T, sep=";", dec=".")
names(taxonomy)
names(samples)
names(macro_data)

# merging the inputs
all_data<-merge(taxonomy,macro_data,by="GBIF")
all_data<-merge(all_data,samples,by="Sampling_Code")
all_data$Sampling_Code_Sub<-paste(all_data$Sampling_Code,all_data$Substrate, sep="-")

# select the relevant columns
sel_data<-all_data[,c(1,40,37,38,34,21,36)]
names(sel_data)

# widen the data frame
macro_habitat<-as.data.frame(pivot_wider(sel_data,id_cols=c(1:5),names_from="Name",values_from="N_individuals",values_fill =0))
head(macro_habitat)
# write.table(macro_habitat,"macro_habitat.csv", sep = ";", dec = ".", row.names = F,na = "", col.names = TRUE)


# CONSIDERING LAKE AND INLET/OUTLET
macro_lake_inout_sample<-aggregate(macro_habitat[,c(2:4)],by=list(macro_habitat$Sampling_Code),unique)
macro_lake_inout<-aggregate(macro_habitat[,7:ncol(macro_habitat)],by=list(macro_habitat$Sampling_Code),sum)
names(macro_lake_inout_sample)[1]<-"Sampling_Code"
names(macro_lake_inout)[1]<-"Sampling_Code"
head(macro_lake_inout)

bio<-t(macro_lake_inout[,-1])
colnames(bio)<-macro_lake_inout[,1]
sample_i<-data.frame(rownames(bio),bio[,1])
names(sample_i)<-c("taxa",colnames(bio)[1])
data<-sample_i
res<-ACI_R(sample_i)
for(i in 2:ncol(bio)){
sample_i<-data.frame(rownames(bio),bio[,i])
names(sample_i)<-c("taxa",colnames(bio)[i])
res<-rbind(res,ACI_R(sample_i))
}

res_lake_inout<-merge(res,macro_lake_inout_sample,by.x = "Sample", by.y = "Sampling_Code")
res_lake_inout<-merge(res_lake_inout,macro_lake_inout,by.x = "Sample", by.y = "Sampling_Code")
head(res_lake_inout)
nrow(res_lake_inout)

write.table(res_lake_inout,"res_lake_inout.csv", sep = ";", dec = ".", row.names = F, col.names = TRUE)



# CONSIDERING LAKE AS A WHOLE
macro_lake_sample<-aggregate(macro_habitat[,c(2:4)],by=list(factor(macro_habitat$Lake_Code):factor(macro_habitat$Date)),unique)
macro_lake<-aggregate(macro_habitat[,7:ncol(macro_habitat)],by=list(factor(macro_habitat$Lake_Code):factor(macro_habitat$Date)),sum)
names(macro_lake_sample)[1]<-"Sampling_Code_Lake"
names(macro_lake)[1]<-"Sampling_Code_Lake"
head(macro_lake)

bio<-t(macro_lake[,-1])
colnames(bio)<-macro_lake[,1]
sample_i<-data.frame(rownames(bio),bio[,1])
names(sample_i)<-c("taxa",colnames(bio)[1])
data<-sample_i
res<-ACI_R(sample_i)
for(i in 2:ncol(bio)){
sample_i<-data.frame(rownames(bio),bio[,i])
names(sample_i)<-c("taxa",colnames(bio)[i])
res<-rbind(res,ACI_R(sample_i))
}

res_lake<-merge(res,macro_lake_sample,by.x = "Sample", by.y = "Sampling_Code_Lake")
res_lake<-merge(res_lake,macro_lake,by.x = "Sample", by.y = "Sampling_Code_Lake")
head(res_lake)
nrow(res_lake)
res_lake <- apply(res_lake,2,as.character)
write.table(res_lake,"res_lake.csv", sep = ";", dec = ".", row.names = F, col.names = TRUE)
