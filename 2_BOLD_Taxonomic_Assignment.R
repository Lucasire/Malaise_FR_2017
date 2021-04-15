##### ------------------------------------------------------------------------
## R script to do the taxonomic assignation using BOLDsystem

##########################################################################################################
##########################################################################################################
# Taxonomic assignment
##########################################################################################################
##########################################################################################################



### ------------------------------------------------------------------------
## Libraries



# remotes::install_github("ropensci/bold@async")



library("tidyverse")
library("bold")
library("plyr") # always load plyr before dplyr
library("seqinr") # used to read fasta file
library("dplyr")
library("readr")
library("stringi")
library("data.table")
library("magicfor")
library("naniar")



### ------------------------------------------------------------------------
## General options



options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE) # if you want to pass arguments to R from the shell (bash) command line
print(args)


otutablefolder <- args[1]


path_name <- file.path(paste0("PATH"))
setwd(path_name)
print(path_name)
getwd()



### ------------------------------------------------------------------------
## Files preparation



# Create the loop of files to loop over



fileNames = list.files(path = "C:/Users/lucas/Desktop/bioinfo/BOLD/tabletest", pattern="table_ABC_97_lulu*")



# Create a list for BOLD outputs



Taxonomy_ID <- vector('list')



### ------------------------------------------------------------------------
## BOLD assignement



# Loop for taxonomic assignation



for (i in 1:length(fileNames)) {
  
    print(i)
  
    print(Sys.time())
       
    file <- read.fasta(fileNames[i], seqtype = "DNA", set.attributes = FALSE, as.string = TRUE, forceDNAtolower = FALSE)
       
    boldoutput_ID <- bold_identify(file, db = "COX1")
     
    Taxonomy_ID <- c(Taxonomy_ID, boldoutput_ID)

    print(Sys.time())
  
}


# Remove any 'NULL' elements



Taxonomy_ID_Curated <- Taxonomy_ID[!sapply(Taxonomy_ID, is.null)]



# Compare number of 'NULL' elements removed



print("Number of OTUs processed through BOLD")

a <- length(Taxonomy_ID)

print(a)



print("Number of OTUs with taxonomic assignation")

b <- length(Taxonomy_ID_Curated)

print(b)



print("Number of OTUs without taxonomic assignation")

print(a-b)



### ------------------------------------------------------------------------
## BOLD parent taxonomic assignement



# Create a list for BOLD outputs



Taxonomy_Parents <- vector('list')



# Loop for parent taxonomic assignation



for (i in 1:length(Taxonomy_ID_Curated)) {

    print(i)

    print(Sys.time())
  
    boldoutput_Parents <- bold_identify_parents(Taxonomy_ID_Curated[i], wide = TRUE)

    Taxonomy_Parents <- c(Taxonomy_Parents, boldoutput_Parents)

    print(Sys.time())
  
}



### ------------------------------------------------------------------------
## Data handling and Curation



# Convert OTU taxonomy list to dataframe



Taxonomy_Parents.df <- ldply(Taxonomy_Parents, data.frame)
Taxonomy_Parents.df <- read.table("BOLD_ID_DEF.xlsx", header = T, na.strings = "NA")
help("read.table")

# Filter for similarity threshold and rearrange data



Taxonomy_Parents_97_to_100.df <- dplyr::filter(Taxonomy_Parents.df, similarity >= 0.97)



Taxonomy_Parents_97_to_100.df <- dplyr::arrange(Taxonomy_Parents_97_to_100.df, .id, desc(similarity)) # .id is the name of the sequence that was identified



OTU_list_table <- Taxonomy_Parents_97_to_100.df %>% distinct(Taxonomy_Parents_97_to_100.df$.id)



OTU_list <- c(OTU_list_table[,])

rm(OTU_list)


magic_for(print, silent = TRUE)



for (i in 1:length(OTU_list)) {
  
    print(i)
  
    Concatenation <- Taxonomy_Parents_97_to_100.df %>% filter(.id == OTU_list[i])
    
    Concatenation <- Concatenation %>% distinct(Concatenation$taxonomicidentification, .keep_all = TRUE)

    print(Concatenation)
}



Taxonomy_Parents_97_to_100_Curated <- magic_result_as_dataframe() # store the magic results in a dataframe with lists



magic_free()



Taxonomy_Parents_97_to_100_Curated.df <- ldply(Taxonomy_Parents_97_to_100_Curated$Concatenation, data.frame) # change the lists into dataframes and concatenate them together



Taxonomy_Parents_97_to_100_Curated.df$Concatenation.taxonomicidentification <- NULL



Taxonomy_Parents_97_to_100_Curated.df$taxonomicidentification <- NULL



Concatenation_unique <- Concatenation_unique %>% filter(phylum == "Arthropoda")



# Filter for finer taxonomic rank



magic_for(print, silent = TRUE)



for (i in 1:length(OTU_list)) {
  
  print(i)
  
  Concatenation_unique <- Taxonomy_Parents_97_to_100_Curated.df %>% filter(.id == OTU_list[i])
    
  if (length(unique(na.omit(Concatenation_unique$class), incomparables = FALSE)) > 1) {
    
    phylum <- unique(na.omit(Concatenation_unique$phylum, incomparables = FALSE))
    
    Concatenation_unique <- Concatenation_unique %>% replace_with_na_at(., .var = c("class_id", "class", "order_id", "order", "family_id", "family", "subfamily_id", "subfamily", "genus_id", "genus", "species_id", "species"), condition = ~.x != 0)
    
    Concatenation_unique <- Concatenation_unique %>% filter(phylum == phylum)
    
    Concatenation_unique <- Concatenation_unique %>% distinct(Concatenation_unique$phylum, .keep_all = TRUE)
    
    colnames(Concatenation_unique)[29] <- "Taxonomy"
    
    print(Concatenation_unique)
    
        } else if (length(unique(na.omit(Concatenation_unique$order), incomparables = FALSE)) > 1) {
      
        class <- unique(na.omit(Concatenation_unique$class, incomparables = FALSE))
      
        Concatenation_unique <- Concatenation_unique %>% replace_with_na_at(., .var = c("order_id", "order", "family_id", "family", "subfamily_id", "subfamily", "genus_id", "genus", "species_id", "species"), condition = ~.x != 0)
      
        Concatenation_unique <- Concatenation_unique %>% filter(class == class)
      
        Concatenation_unique <- Concatenation_unique %>% distinct(Concatenation_unique$class, .keep_all = TRUE)
      
        colnames(Concatenation_unique)[29] <- "Taxonomy"
      
        print(Concatenation_unique)
      
          } else if (length(unique(na.omit(Concatenation_unique$family), incomparables = FALSE)) > 1) {
        
          order <- unique(na.omit(Concatenation_unique$order), incomparables = FALSE)
        
          Concatenation_unique <- Concatenation_unique %>% replace_with_na_at(., .var = c("family_id", "family", "subfamily_id", "subfamily", "genus_id", "genus", "species_id", "species"), condition = ~.x != 0)
        
          Concatenation_unique <- Concatenation_unique %>% filter(order == order)
        
          Concatenation_unique <- Concatenation_unique %>% distinct(Concatenation_unique$order, .keep_all = TRUE)
        
          colnames(Concatenation_unique)[29] <- "Taxonomy"
        
          print(Concatenation_unique)
        
            } else if (length(unique(na.omit(Concatenation_unique$subfamily), incomparables = FALSE)) > 1) {
          
            family <- unique(na.omit(Concatenation_unique$family), incomparables = FALSE)
          
            Concatenation_unique <- Concatenation_unique %>% replace_with_na_at(., .var = c("subfamily_id", "subfamily", "genus_id", "genus", "species_id", "species"), condition = ~.x != 0)
          
            Concatenation_unique <- Concatenation_unique %>% filter(family == family)
          
            Concatenation_unique <- Concatenation_unique %>% distinct(Concatenation_unique$family, .keep_all = TRUE)
          
            colnames(Concatenation_unique)[29] <- "Taxonomy"
          
            print(Concatenation_unique)
          
              } else if (length(unique(na.omit(Concatenation_unique$genus), incomparables = FALSE)) > 1) {
            
                subfamily <- unique(na.omit(Concatenation_unique$subfamily), incomparables = FALSE)
            
                Concatenation_unique <- Concatenation_unique %>% replace_with_na_at(., .var = c("genus_id", "genus", "species_id", "species"), condition = ~.x != 0)
            
                Concatenation_unique <- Concatenation_unique %>% filter(subfamily == subfamily)
            
                Concatenation_unique <- Concatenation_unique %>% distinct(Concatenation_unique$subfamily, .keep_all = TRUE)
            
                colnames(Concatenation_unique)[29] <- "Taxonomy"
            
                print(Concatenation_unique)
            
                  } else if (length(unique(na.omit(Concatenation_unique$species), incomparables = FALSE)) > 1) {

                  genus <- unique(na.omit(Concatenation_unique$genus), incomparables = FALSE)
              
                  Concatenation_unique <- Concatenation_unique %>% replace_with_na_at(., .var = c("species_id", "species"), condition = ~.x != 0)
              
                  Concatenation_unique <- Concatenation_unique %>% filter(genus == genus)
              
                  Concatenation_unique <- Concatenation_unique %>% distinct(Concatenation_unique$genus, .keep_all = TRUE)
              
                  colnames(Concatenation_unique)[29] <- "Taxonomy"
              
                  print(Concatenation_unique)
              
                    } else {
                
                    species <- unique(na.omit(Concatenation_unique$species), incomparables = FALSE)
                
                    Concatenation_unique <- Concatenation_unique %>% filter(species == species)
                    
                    Concatenation_unique <- Concatenation_unique %>% distinct(Concatenation_unique$species, .keep_all = TRUE)
                    
                    colnames(Concatenation_unique)[29] <- "Taxonomy"
              
                    print(Concatenation_unique)
                
                    }
  
}



Taxonomy_Parents_97_to_100_Unique <- magic_result_as_dataframe() # store the magic results in a dataframe with lists



magic_free()



Taxonomy_Parents_97_to_100_Unique.df <- ldply(Taxonomy_Parents_97_to_100_Unique$Concatenation_unique, data.frame) # change the lists into dataframes and concatenate them together



Filtered_OTU_list <- Taxonomy_Parents_97_to_100_Unique.df %>% select(., .id)



### ------------------------------------------------------------------------
## Data back-up



# List of OTUs with taxonomic assignation at Arthropoda level


write.table(Filtered_OTU_list, file = "Filtered_OTU_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


# Excel copy


write.table(Filtered_OTU_list, file = "Filtered_OTU_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


write.csv(Taxonomy_Parents_97_to_100_Unique.df, file = "BOLD_ID_Unique_ABC_97.xlsx")



write.csv(Taxonomy_Parents_97_to_100_Curated.df, file = "BOLD_ID_ABC_97.xlsx")



# Hardware Disk copy



write_tsv(Taxonomy_Parents_97_to_100_Unique.df, "BOLD_ID_Unique_ABC_97.tsv")



write_tsv(Taxonomy_Parents_97_to_100_Curated.df, "BOLD_ID_ABC_97.tsv")



write_tsv(Taxonomy_Parents.df, "BOLD_ID_ABC.tsv")