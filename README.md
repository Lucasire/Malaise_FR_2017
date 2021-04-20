# Malaise_FR_2017
Files necessary to reproduce analyses from Sire et al. 'Climate-induced forest dieback drives compositional change in insect communities that is concentrated amongst rare species.' Preprint deposited in BioRxv. (2021)

Files labeled (1) are used in the bioinformatic demultiplexing process. Both (1.1) and (1.2) are called within the DAMe pipeline and paths should be changed accordingly for running the whole script.

Second step (2) is the loop created to assign taxonomy to generated MOTUs by using R package 'bold'. Manual curation of data have been made on the .xlsx result file to merge MOTUs with identical taxonomic match, resulting in a more conservative approach.

Final step corresponding to statistical analyses can be reproduced using files (3.1.) and (3.2.) with their associated datasets.
