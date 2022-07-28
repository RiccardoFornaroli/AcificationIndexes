# AcificationIndexes
R script to calculate multiple acidification indexes

The tool was specifically developed to allow the calculation of acidification specific indexes based on the occurrence or on the relative abundance of macroinvertebrate taxa sensitive to acidification effects, the indices for which the calculation is already implemented are reported below: 
1.	Raddum 1988 index (Raddum et al., 1988), developed to evaluate the effects of acidification on rivers and lakes in Norway; 
2.	Raddum 1990 index (Fjellheim & Raddum, 1990), developed to evaluate the effects of acidification on rivers and lakes in Norway;  
3.	NIVA index (Bækken & Kjellberg, 2004), developed to evaluate the effects of acidification on humus-rich streams in eastern Norway; 
4.	AWICfam index (Davy-Bowker et al. 2003, 2005) developed to evaluate the effects of acidification on streams and rivers in England and Wales;
5.	AWICsp index (Davy-Bowker et al. 2003) developed to evaluate the effects of acidification on streams and rivers in England and Wales; 
6.	Braukmann index (Braukmann & Biss, 2004) developed to evaluate the effects of acidification on streams and rivers in Germany; 
7.	LAMM index (McFarland et al. 2010) developed to evaluate the effects of acidification on clear and humic lakes in the UK; 
8.	TL index (Hämäläinen & Huttunen, 1990) developed to evaluate the effects of acidification on streams and rivers in Finland;

List of files:
"Indexes_formulas.r" is the R script with the function "ACI_R" that evaluate the indexes for a sample.
"Acidification_indexes.r" is the R script that apply the function "ACI_R" to the Italian data as a test case.
"Taxonomy_18_01_2022.csv" is the taxonomic database used to classify taxa and store sensitivity level for the different index when available.
"Bio_lit_Samples.csv" is the list of Italian samples to be used as a test case.
"Bio_lit_UPDATED.csv" is the file with the abundance of taxa in each sample in long format.
