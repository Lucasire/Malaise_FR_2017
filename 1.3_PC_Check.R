##### ------------------------------------------------------------------------
## R script to create de file for internal Blast of the positive control and decide the DAMe filtration parameters

### ------------------------------------------------------------------------
## Libraries

library(tidyverse)

### ------------------------------------------------------------------------
## General options

experiment <- c("ABC")

options(echo=TRUE) # To see commands in output file
args <- commandArgs(trailingOnly = TRUE) # To pass arguments to R from the shell (bash) command line
print(args)

otutablefolder <- args[1]

path_name <- file.path(paste0("PATH", otutablefolder, "/OTU_tables/"))
setwd(path_name)
print(path_name)

### ------------------------------------------------------------------------
## Formatting table

otutable_ABC_97 <- read.table("table_ABC_97_lulu.txt", header = T, sep = "\t")
otutable_ABC_97 <- otutable_ABC_97 %>% dplyr::mutate(OTU_PC_tot= C._1 + C._2 + C._3 + C._4 + C._5 + C._6) # To merge all positive controls together in a new column
PC_ABC_97 <- otutable_ABC_97 %>% dplyr::select("OTU", "OTU_PC_tot", "C._1", "C._2", "C._3", "C._4", "C._5", "C._6") # To select only the informative columns

### ------------------------------------------------------------------------
## Removing OTUs that are not in the positive controls

PC_ABC_97 <- PC_ABC_97 %>% filter(OTU_PC_tot > 0)
PC_OTU_ABC_97 <- PC_ABC_97 %>% select(OTU)
write.table(PC_OTU_ABC_97, file="PC_OTU_ABC_97.txt",row.names=FALSE, col.names = FALSE, quote=FALSE)

##### ------------------------------------------------------------------------
## END
