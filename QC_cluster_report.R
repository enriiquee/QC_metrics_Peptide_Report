# PRIDE Cluster clusters all MS/MS spectra submitted to PRIDE Archive repository release: 2015-04
# http://www.ebi.ac.uk/pride/cluster/
 
# Description: The present script provides a reliable QC (Quality control) report about peptides in PRIDE Cluster. 

# INPUT: The input files must be in the same directory as the script, with the names: 
# pride_cluster_peptides_ALL.tsv AND pride_cluster_peptides_ALL2.tsv, being the file "pride_cluster_peptides_ALL2.tsv" the new release to comparee.
knitr::spin 
# Upload packages
packages <- c("data.table", "dplyr", "ggplot2", "stringr", "reshape2")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))  }

library("data.table"); library("dplyr"); library("ggplot2"); library("stringr"); 
library(knitr); library("reshape2")
#Upload files using getwd
pride_cluster_peptides_ALL_version1  <- read.delim(file.path(getwd(),"pride_cluster_peptides_Human.tsv"), comment.char="#")
pride_cluster_peptides_ALL_version2 <- read.delim(file.path(getwd(),"pride_cluster_peptides_Human2.tsv"), comment.char="#")

#Split of the column in PEP and CPE
new_data_frame_PEP_v1 = subset(pride_cluster_peptides_ALL_version1, PEH == 'PEP')
new_data_frame_CPE_v1 = subset(pride_cluster_peptides_ALL_version1, PEH == 'CPE')

new_data_frame_PEP_v2 = subset(pride_cluster_peptides_ALL_version2, PEH == 'PEP')
new_data_frame_CPE_v2 = subset(pride_cluster_peptides_ALL_version2, PEH == 'CPE')

# Checking repeated sequences.
# Here, we call duplicated twice: first from the start of the sequence column to the 
# end and then from the end back to the start using fromLast. duplicated returns a 
# logical which is or'ed (i.e., |) so that we retrieve all the row indices that have 
# duplicates. We then subset rows of new_data_frame_PEP with respect to those. To 
# accomplish the same with more than one column so that we find all rows that have 
# duplicated values in both the sequence and modification columns, we need to select 
# those columns to pass to duplicated. This can be done using subset:

cat("===============================================================================\n")
cat( "REPEATED SEQUENCES RESULTS:\n")
cat("===============================================================================\n")

res_v1 <- new_data_frame_PEP_v1[duplicated(new_data_frame_PEP_v1$sequence) | duplicated(new_data_frame_PEP_v1$sequence, fromLast=TRUE),]
seq.mod_v1 <- subset(new_data_frame_PEP_v1, select=c("sequence","modifications"))
data_duplicate_v1 <- new_data_frame_PEP_v1[duplicated(seq.mod_v1) | duplicated(seq.mod_v1, fromLast=TRUE),]


if (length(data_duplicate_v1$sequence)== 0) {
    cat("No repeated values in release 1. \n")
    
} else {
    cat("Repeated values in release 1.\n")
    print(data_duplicate_v1)
}

#New data: 
res_v2 <- new_data_frame_PEP_v2[duplicated(new_data_frame_PEP_v2$sequence) | duplicated(new_data_frame_PEP_v2$sequence, fromLast=TRUE),]
seq.mod_v2 <- subset(new_data_frame_PEP_v2, select=c("sequence","modifications"))
data_duplicate_v2 <- new_data_frame_PEP_v2[duplicated(seq.mod_v2) | duplicated(seq.mod_v2, fromLast=TRUE),]


if (length(data_duplicate_v2$sequence)== 0) {
    cat("No repeated values in release 2.\n")
    
} else {
    cat("Repeated values in release 1.\n \n")
    print(data_duplicate_v2)
}
cat("===============================================================================\n")
cat( "TOTAL NUMBER OF PEPTIDES:\n")
cat("===============================================================================\n")

#Total number: 
cat("Release1:", length(new_data_frame_PEP_v1$sequence),"\nRelease2:", length(new_data_frame_PEP_v2$sequence),"\n")
#Conditional
if (length(new_data_frame_PEP_v1$sequence)>length(new_data_frame_PEP_v2$sequence)) {
    cat("Peptides reduced:", (length(new_data_frame_PEP_v1$sequence)-length(new_data_frame_PEP_v2$sequence)),"\n")
} else if ( length(new_data_frame_PEP_v1$sequence)<length(new_data_frame_PEP_v2$sequence)) {
    cat("Peptides increased:", (length(new_data_frame_PEP_v2$sequence)-length(new_data_frame_PEP_v1$sequence)),"\n")
} else
    cat("The number of peptides are equals \n")


cat("===============================================================================\n")
cat( "PEPTIDES WITH NO-MODIFICATIONS\n")
cat("===============================================================================\n")
cat("Release1:",length(new_data_frame_PEP_v1$modifications[new_data_frame_PEP_v1$modifications == "NULL"]), "Release2: ", length(new_data_frame_PEP_v2$modifications[new_data_frame_PEP_v2$modifications == "NULL"]))
if (length(new_data_frame_PEP_v1$modifications[new_data_frame_PEP_v1$modifications == "NULL"])>length(new_data_frame_PEP_v2$modifications[new_data_frame_PEP_v2$modifications == "NULL"])) {
    cat("Peptides with no-modification reduced:\n", (length(new_data_frame_PEP_v1$modifications[new_data_frame_PEP_v1$modifications == "NULL"])-length(new_data_frame_PEP_v2$modifications[new_data_frame_PEP_v2$modifications == "NULL"])))
} else if ( length(new_data_frame_PEP_v1$modifications)<length(new_data_frame_PEP_v2$modifications)) {
    cat("Peptides with no-modification increased:\n", (length(new_data_frame_PEP_v2$modifications[new_data_frame_PEP_v2$modifications == "NULL"])-length(new_data_frame_PEP_v1$modifications[new_data_frame_PEP_v1$modifications == "NULL"])))
} else
    cat("The number of peptides with no-modifications are equals.\n")




cat("===============================================================================\n")
cat( "PEPTIDES WITH MODIFICATIONS\n")
cat("===============================================================================\n")
cat("v1: ",length(new_data_frame_PEP_v1$modifications[new_data_frame_PEP_v1$modifications != "NULL"]), "v2: ", length(new_data_frame_PEP_v2$modifications[new_data_frame_PEP_v2$modifications != "NULL"]))
if (length(new_data_frame_PEP_v1$modifications[new_data_frame_PEP_v1$modifications != "NULL"])>length(new_data_frame_PEP_v1$modifications[new_data_frame_PEP_v2$modifications != "NULL"])) {
    cat("Peptides with modification reduced:\n", (length(new_data_frame_PEP_v1$modifications[new_data_frame_PEP_v1$modifications != "NULL"])-length(new_data_frame_PEP_v2$modifications[new_data_frame_PEP_v2$modifications != "NULL"])))
} else if ( length(new_data_frame_PEP_v1$modifications)<length(new_data_frame_PEP_v2$modifications)) {
    cat("Peptides with modification increased:\n", (length(new_data_frame_PEP_v2$modifications[new_data_frame_PEP_v2$modifications != "NULL"])-length(new_data_frame_PEP_v1$modifications[new_data_frame_PEP_v1$modifications != "NULL"])))
} else
    cat("The number of peptides with modifications are equals.\n")

cat("===============================================================================\n")
cat( "NEW PEPTIDES:\n")
cat("==============================================================================\n")

df <- new_data_frame_PEP_v1
df2 <- new_data_frame_PEP_v2

#First convert your peptide counts to numeric (they're a factor with numeric character labels, that's a bit messed up):

df$peptideNumberSpectra=as.numeric(as.character(df$peptideNumberSpectra))
df2$peptideNumberSpectra=as.numeric(as.character(df2$peptideNumberSpectra))


df_final<- df %>% 
    full_join(df2, by = c("sequence", "modifications"), suffix = c(".1", ".2")) %>%
    # Fix data to convert to character and numeric
    mutate_each(funs(as.numeric(as.character(.))), starts_with("pept")) %>%
    # See difference
    mutate(change = peptideNumberSpectra.2 - peptideNumberSpectra.1)

#How many new peptides there are in the new release. 
df_peptides <- df_final[is.na(df_final$PEH.1),]
df_peptides2 <- data.frame(sequences=df_peptides$sequence, modifications=df_peptides$modifications)

df_peptides3 <- df_final[is.na(df_final$PEH.2),]
df_peptides4 <- data.frame(sequences=df_peptides3$sequence, modifications=df_peptides3$modifications)

#To know if the new release has obtained. 
if (length(df_final[is.na(df_final$PEH.1),])==0) {
    cat("The new release has got the same peptides than release before\n")
} else {
    cat("The new release has  obtained new peptides\n")
    print(df_peptides2)
}

#To know if the release has lost peptides.
if (length(df_final[is.na(df_final$PEH.2),])!=0) {
    cat("The new release has lost peptides\n")
    print(df_peptides4)
} else {
    cat("The new release has not lost peptides\n")
}

cat("===============================================================================\n")
cat( "NEW PEPTIDES SPECTRA:\n")
cat("===============================================================================\n")


df$peptideNumberSpectra=as.numeric(as.character(df$peptideNumberSpectra))
df2$peptideNumberSpectra=as.numeric(as.character(df2$peptideNumberSpectra))


df_final<- df %>% 
    full_join(df2, by = c("sequence", "modifications"), suffix = c(".1", ".2")) %>%
    # Fix data to convert to character and numeric
    mutate_each(funs(as.numeric(as.character(.))), starts_with("pept")) %>%
    # See difference
    mutate(change = peptideNumberSpectra.2 - peptideNumberSpectra.1)


#df_finalX = subset(final_df, change != 'NA') Remove NA from the column
df_final2 <- subset(df_final, change == 0) 
df_final3 <- data.frame(sequence= df_final2$sequence, modifications= df_final2$modifications, Change=df_final2$change)



if (sum(df_final2$change == 0)) {
    cat("The number of spectras are the same\n")
} else{
    cat("The number of Spectras are different\n")
    print(df_final3)
}

cat("===============================================================================\n")
cat( "NEW PROJECTS:\n")
cat("===============================================================================\n")

df$numberProjects=as.numeric(as.character(df$numberProjects))
df2$numberProjects=as.numeric(as.character(df2$numberProjects))


df_final_project<- df %>% 
    full_join(df2, by = c("sequence", "modifications"), suffix = c(".1", ".2")) %>%
    # Fix data to convert to character and numeric
    mutate_each(funs(as.numeric(as.character(.))), starts_with("proyects")) %>%
    # See difference
    mutate(change = numberProjects.2 - numberProjects.1)


#df_finalX = subset(final_df, change != 'NA') Remove NA from the column
df_final_project2 <- subset(df_final_project, change != 0) 
df_final_project3 <- data.frame(sequence= df_final_project2$sequence, modifications= df_final_project2$modifications, Change=df_final_project2$change)

if (sum(df_final_project2$change == 0)) {
    cat("The number of spectras are the same\n")
} else {
    cat("The number of Spectras are different\n")
    print(df_final_project3)
}

cat("===============================================================================\n")
cat( "NEW CLUSTERS:\n")
cat("===============================================================================\n")

df$numberClusters=as.numeric(as.character(df$numberClusters))
df2$numberClusters=as.numeric(as.character(df2$numberClusters))


df_final_clusters<- df %>% 
    full_join(df2, by = c("sequence", "modifications"), suffix = c(".1", ".2")) %>%
    # Fix data to convert to character and numeric
    mutate_each(funs(as.numeric(as.character(.))), starts_with("clusters")) %>%
    # See difference
    mutate(change = numberClusters.2 - numberClusters.1)


#df_finalX = subset(final_df, change != 'NA') Remove NA from the column
df_final_clusters2 <- subset(df_final_clusters, change != 0) 
df_final_clusters3 <- data.frame(sequence= df_final_clusters2$sequence, modifications= df_final_clusters2$modifications, Change=df_final_clusters2$change)

if (sum(df_final_clusters2$change == 0)) {
    cat("The number of spectras are the same\n")
} else {
    cat("The number of Spectras are different\n")
    print(df_final_clusters3)
}

#HISTOGRAM: 
#Para preparar el histograma tenemos que clasificar los tipos de modificaciones, POSICION-DATABA-ID. 
# Para ello cogemos los datos que no sean NULL y eliminamos los espacios en banco. 

histo1 = subset(new_data_frame_PEP_v1, modifications != 'NULL')
histo1[histo1==""] <- NA
histo1 = subset(histo1, modifications != 'NA')

histo1_2 = subset(new_data_frame_PEP_v2, modifications != 'NULL')
histo1_2[histo1_2==""] <- NA
histo1_2 = subset(histo1_2, modifications != 'NA')

#Split dataset. 
histo2 <- data.frame(str_split_fixed(histo1$modifications, ",", 20))
histo2_2 <- data.frame(str_split_fixed(histo1_2$modifications, ",", 20))

#If you want to check how many columns are empty, you can use the code below: 
#columns_emply <- histo2[!sapply(histo2, function(x) all(x == ""))]
#columns_emply <- histo2_2[!sapply(histo2_2, function(x) all(x == ""))]

#Merge the columns in one. 
histo3 <- melt(setDT(histo2),                              # set df to a data.table
               measure.vars = list(c(1:20)),    # set column groupings
               value.name = 'V')[                      # set output name scheme
                   , -1, with = F]

histo3_2 <- melt(setDT(histo2_2),                              # set df to a data.table
                 measure.vars = list(c(1:20)),    # set column groupings
                 value.name = 'V')[                      # set output name scheme
                     , -1, with = F]

#Remove white rows.  
histo3[histo3==""] <- NA
histo3 = subset(histo3, V1 != " ")

histo3_2[histo3_2==""] <- NA
histo3_2 = subset(histo3_2, V1 != " ")

#Remove first part of the string [num]-
histo4 <- data.frame(modifications=gsub(" [A-Za-z] ", "", gsub("[0-9]*-", "", histo3$V1)))
histo4_2 <- data.frame(modifications.2=gsub(" [A-Za-z] ", "", gsub("[0-9]*-", "", histo3_2$V1)))

#Histograma: 
histo5 <- data.frame(table(histo4))
histo5_2 <- data.frame(table(histo4_2))

ggplot(data=histo5, aes(x=histo4, y=Freq, fill=histo4)) +
    geom_bar(stat="identity") + guides(fill=FALSE)+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data=histo5_2, aes(x=histo4_2, y=Freq, fill=histo4_2)) +
    geom_bar(stat="identity") + guides(fill=FALSE)+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


#Prueba para juntar dos barplot

#df <- cbind.data.frame(histo1, histo1_2[match(histo1$modifications, histo1_2$modifications), ])
#colnames(df) <- c("PEH", "sequence", "modifications", "bestRank","bestScore","peptideNumberSpectra", "numberProjects","numberClusters","species","projects","PEH.1", "sequence.1", "modifications.1","bestRank.1", "bestScore.1", "peptideNumberSpectra.1", "numberProjects.1", "numberClusters.1" ,  "species.1" , "projects.1")

df <- cbind.data.frame(histo4, histo4_2[match(histo4$modifications, histo4_2$modifications.2), ])
colnames(dataprueba) <- c("modifications.1","modifications.2")

ggplot(melt(df,measure.vars = names(df)), aes(x = value, fill = variable)) + 
    geom_bar(stat = "count", position = "dodge") + 
    theme(axis.text.x = element_text(angle = 20, hjust = 0.5, vjust = -0.1)) + 
    guides(fill=FALSE)+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
    scale_fill_manual(values = c('red', 'blue'))



