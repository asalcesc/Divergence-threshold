##########################################################################
# SCRIPT 1. Applying genetic divergence threshold to establish de novo LMDH or PBS
## Created by Carmelo Ándujar 2020
## Supporting material of Salces-Castellano et al. (2020) Climate drives community‐wide divergence within species over a limited spatial scale: evidence from an oceanic island
##########################################################################


#####################################
### Execution
#####################################

### Setting working directory, containing complete alignment PUs with a code "XX + number"
setwd ("working-directory") 

### Loading libraries
library(ape)
library(stringr)
library(phangorn)

### Loading function
source("Function_PU_Congruence_from_alignment.R")

### Reading input alignment
alignment <- ("Alignment.fasta")

### Executing function
PU_Congruence_from_alignment(alignment, "TF", "K80", 0.04, print.subtrees.fasta = "NO")


##FUNCTION PARAMETERS
#PU_Congruence_from_alignment <- function (nameOFfasta,codePU, distmodel,limite1, print.subtrees.fasta = "YES") 
#IN THIS ORDER:  
#1. nameOFfasta #complete alignment PUs with a code "XX + number"
#2. codePU #letters before the PU number, all PUs should have the same letter e.g. "LG"
#3. distmodel #e.g.  "F84". A character string specifying the evolutionary model to be used; must be one of "raw", "N", "TS", "TV", "JC69", "K80" (the default), "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", or "indelblock". For details see help dist.dna {ape}
#4. limite1 #number indicating divergent (ratio in one unit) per lineage. eg. 0.04 means 8% divergence
#5. print.subtrees.fasta = "YES"  # logic with info for printing or not .fasta files for each subtrees (could be yes or no)


