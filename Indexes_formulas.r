
Shannon_ind <- function(data){
	pin <- data[,2]/sum(data[,2])
	lnpi <- log(data[,2]/sum(data[,2]))
	return(round(-sum(pin*lnpi),3))
	}
	
# AWIC index at family level
AWIC_fam_ASPT_ind <- function(data){
	AWIC_fam_table <- unique(data[!is.na(data$AWICfam_taxa),c(6,7)])
	ASPT <- NA
	if (length(AWIC_fam_table[,2]) > 0) {
	ASPT <- (round(sum(AWIC_fam_table$AWICfam)/length(AWIC_fam_table[,2]),3))
	} else {
	ASPT <- NA
	}	
	return(ASPT)
	}
	
# AWIC index at specie level
AWIC_sp_ASPT_ind <- function(data){
	AWIC_sp_table <- unique(data[!is.na(data$AWICsp_taxa),c(8,9)])
	ASPT <- NA
	if (length(AWIC_sp_table[,2]) > 0) {
	ASPT <- (round(sum(AWIC_sp_table$AWICsp)/length(AWIC_sp_table[,2]),3))
	} else {
	ASPT <- NA
	}	
	return(ASPT)
	}
	
# LAMM index
LAMM_ind <- function(data){
	LAMM_table <- data[!is.na(data$LAMM_Sk),c(2,12,13)]
	LAMM_table$perc_scoring<-LAMM_table[,1]/sum(LAMM_table[,1])
	LAMM_table$Hk<-ifelse(LAMM_table$perc_scoring<0.05,1,ifelse(LAMM_table$perc_scoring<0.2,3,5))
	LAMM <- NA
	if (length(LAMM_table[,2]) > 0) {
	LAMM <- round((sum(LAMM_table$LAMM_Sk*LAMM_table$LAMM_Wk*LAMM_table$Hk)/sum(LAMM_table$LAMM_Wk*LAMM_table$Hk)),3)
	} else {
	LAMM <- NA
	}	
	return(LAMM)
	}

#Braukmann index
Braukmann_ind <- function(data){
	# calcolo abbondanze relative
	data$rel_ab<-data[,2]/sum(data[,2])
	ab_1<-sum(data$rel_ab[data$Braukmann==1], na.rm=T)
	ab_2<-sum(data$rel_ab[data$Braukmann==1|data$Braukmann==2], na.rm=T)
	ab_3<-sum(data$rel_ab[data$Braukmann==1|data$Braukmann==2|data$Braukmann==3], na.rm=T)
	ab_4<-sum(data$rel_ab[data$Braukmann==1|data$Braukmann==2|data$Braukmann==3|data$Braukmann==4], na.rm=T)
	ab_5<-sum(data$rel_ab[data$Braukmann==1|data$Braukmann==2|data$Braukmann==3|data$Braukmann==4|data$Braukmann==5], na.rm=T)
	Brauk<-NA
	if (ab_1 >= 0.1) {
	Brauk<-1
	} else if (ab_2 >= 0.1){
	Brauk<-2
	} else if (ab_3 >= 0.1){
	Brauk<-3
	} else if (ab_4 >= 0.1){
	Brauk<-4
	} else if (ab_5 >= 0.1){
	Brauk<-5
	}
	
return(Brauk)
}

# Acidification indexes function
# "data" should be a data frame with two columns, "taxa" with the name 
# of the taxa to be matched in the taxonomy database, "sample code" with
# the number of individuals in the sample
ACI_R <- function(data){
	# load reference data, file "Taxonomy_18_01_2022.csv"
	base <- read.csv("Taxonomy_18_01_2022.csv", header=T, sep=";", dec=".",na.strings = "NA");
		
	# format data
	data <- data[which(data[,2]>0),]
	data[,3:19] <- base[match(data[,1],base$Name),c(22:32,9:11,14,16,17)]
	data<-droplevels(data)

	# indexes calculation
	N_Taxa <- length(data[,1])
	Ind_Tot <- sum(data[,2])
	N_Pleco_Taxa <- length(data[!is.na(data$Order) & data$Order=="Plecoptera",1])
	N_Ephem_Taxa <- length(data[!is.na(data$Order) & data$Order=="Ephemeroptera",1])
	N_Trico_Taxa <- length(data[!is.na(data$Order) & data$Order=="Trichoptera",1])
	N_EPT_Taxa <- length(data[!is.na(data$Order) & (data$Order=="Plecoptera"|data$Order=="Ephemeroptera"|data$Order=="Trichoptera"),1])
	N_Chiro_Taxa <- length(data[!is.na(data$Family) & data$Family=="Chironomidae",1])
	N_Oligo_Taxa <- length(data[!is.na(data$Class) & data$Class=="Oligochaeta",1])
	Fam_Pleco_Taxa <- length(unique(data[!is.na(data$Order) & data$Order=="Plecoptera","Family"]))
	Fam_Ephem_Taxa <- length(unique(data[!is.na(data$Order) & data$Order=="Ephemeroptera","Family"]))
	Fam_Trico_Taxa <- length(unique(data[!is.na(data$Order) & data$Order=="Trichoptera","Family"]))
	Fam_EPT_Taxa <- length(unique(data[!is.na(data$Order) & (data$Order=="Plecoptera"|data$Order=="Ephemeroptera"|data$Order=="Trichoptera"),"Family"]))
	Abu_Pleco <- sum(data[!is.na(data$Order) & data$Order=="Plecoptera",2])/Ind_Tot
	Abu_Ephem <- sum(data[!is.na(data$Order) & data$Order=="Ephemeroptera",2])/Ind_Tot
	Abu_Trico <- sum(data[!is.na(data$Order) & data$Order=="Trichoptera",2])/Ind_Tot
	Abu_EPT <- sum(data[!is.na(data$Order) & (data$Order=="Plecoptera"|data$Order=="Ephemeroptera"|data$Order=="Trichoptera"),2])/Ind_Tot
	Abu_Chiro <- sum(data[!is.na(data$Family) & data$Family=="Chironomidae",2])/Ind_Tot
	Abu_Oligo <- sum(data[!is.na(data$Class) & data$Class=="Oligochaeta",2])/Ind_Tot
	Shannon <- Shannon_ind(data)
	Raddum_1988 <- ifelse(max(data$Raddum_1988,na.rm=T)>0,max(data$Raddum_1988,na.rm=T),NA)
	Raddum_1990 <- ifelse(max(data$Raddum_1990,na.rm=T)>0,max(data$Raddum_1990,na.rm=T),NA)
	NIVA <- ifelse(min(data$NIVA,na.rm=T)<5,min(data$NIVA,na.rm=T),NA)
	AWIC_fam_ASPT <- AWIC_fam_ASPT_ind(data)
	AWIC_sp_ASPT <- AWIC_sp_ASPT_ind(data)
	TL<-ifelse(max(data$TL,na.rm=T)>0,max(data$TL,na.rm=T),NA)
	Braukmann<-Braukmann_ind(data)
	LAMM <- LAMM_ind(data)
	Raddum_1988_taxa<-toString(data[!is.na(data$Raddum_1988),1])
	Raddum_1990_taxa<-toString(data[!is.na(data$Raddum_1990),1])
	NIVA_taxa<-toString(data[!is.na(data$NIVA),1])
	AWIC_fam_taxa<-toString(unique(data$AWICfam_taxa[!is.na(data$AWICfam)]))
	AWIC_sp_taxa<-toString(unique(data$AWICsp_taxa[!is.na(data$AWICsp)]))
	TL_taxa<-toString(data[!is.na(data$TL),1])
	Braukmann_taxa<-toString(data[!is.na(data$Braukmann),1])
	LAMM_taxa<-toString(data[!is.na(data$LAMM_Sk),1])
	
	# print output: table with all output variables
	out <- data.frame(colnames(data)[2],N_Taxa,Ind_Tot,N_Pleco_Taxa,
	N_Ephem_Taxa,N_Trico_Taxa,N_EPT_Taxa,N_Chiro_Taxa,N_Oligo_Taxa,
	Fam_Pleco_Taxa,Fam_Ephem_Taxa,Fam_Trico_Taxa,Fam_EPT_Taxa,
	Abu_Pleco,Abu_Ephem,Abu_Trico,Abu_EPT,Abu_Chiro,Abu_Oligo,Shannon,
	Raddum_1988,Raddum_1990,NIVA,AWIC_fam_ASPT,AWIC_sp_ASPT,TL,Braukmann,
	LAMM,Raddum_1988_taxa,Raddum_1990_taxa,NIVA_taxa,AWIC_fam_taxa,
	AWIC_sp_taxa,TL_taxa,Braukmann_taxa,LAMM_taxa)
	names(out)[1] <- "Sample"
	out
}
