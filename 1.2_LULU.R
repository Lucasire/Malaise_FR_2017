##### ------------------------------------------------------------------------
## R script to run Lulu pipeline finalizing chimeras removal

### ------------------------------------------------------------------------
## Libraries

library(tidyverse) # includes all data-formatting packages as one big package (e.g. dplyr, tidyr, ggplot2, readr, readxl, tibble, and others)
# library("devtools")
# install_github("tobiasgf/lulu", force = TRUE)
library(lulu)
# sessionInfo()

### ------------------------------------------------------------------------
## General options

experiment <- c("ABC")

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE) # if you want to pass arguments to R from the shell (bash) command line
print(args)

otutablefolder <- args[1]

path_name <- file.path(paste0("PATH", otutablefolder, "/OTU_tables/"))
setwd(path_name)
print(path_name)

### ------------------------------------------------------------------------
## Lulu

for(i in experiment) # input the matchlist files
{
	for(j in 97:97)
	{
		assign(paste0("match_list_", i, "_", j), read.table(paste0("match_list_", i, "_", j, ".txt"), header = FALSE, as.is = TRUE, stringsAsFactors = FALSE))
	}
}

for(i in experiment) # input the sumaclust otu tables
{
	for(j in 97:97)
	{
		assign(paste0("table_", i, "_", j), read.table(paste0("table_", i, "_", j,".txt"), header = TRUE, as.is = TRUE, stringsAsFactors = FALSE))
	}
}

# names(table_*_97)

for(i in experiment) # select only the columns with OTU read numbers
{
	for(j in 97:97)
	{
	    assign(paste0("table_", i, "_", j, "_lulu"), dplyr::select(get(paste0("table_", i, "_", j)), -Seq)) # create OTU table without the Seq
	}
}

# move OTU names to rownames
for(i in experiment) # select only the columns with OTU read numbers
{
	for(j in 97:97)
	{
	    assign(paste0("table_", i, "_", j, "_lulu"), get(paste0("table_", i, "_", j, "_lulu")) %>% column_to_rownames(var = "OTU"))
	}
}

for(i in experiment) # select only the columns with OTU read numbers
{
	for(j in 97:97)
	{
	    assign("otutable_lulu", get(paste0("table_", i, "_", j, "_lulu")))
	    assign("match_list_lulu", get(paste0("match_list_", i, "_", j)))

	    curated_result <- lulu(otutable_lulu, match_list_lulu)
	    curated_table <- rownames_to_column(curated_result$curated_table, var = "OTU")

	    write.table(curated_table, file = paste0("table_", i, "_", j, "_lulu.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

	    # write.table(curated_result$curated_otus, file = paste0("otus_", i, "_", j, "_lulu.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
	    # write.table(curated_result$original_table, file = paste0("table_", i, "_", j, ".txt"), sep = "\t")
	    rm(otutable_lulu)
	    rm(match_list_lulu)
	}
}

##### ------------------------------------------------------------------------
## END
