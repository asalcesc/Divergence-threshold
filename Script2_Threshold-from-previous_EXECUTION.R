##########################################################################
# SCRIPT 2. Applying divergence threshold to define LMDH and PBS from previous ones
## Created by Carmelo Andújar 2020
## Supporting material of (i) Salces-Castellano et al. (2020) doi 10.1111/ele.13433 and (ii) Salces-Castellano et al. (in rev)
##########################################################################


#####################################
### Execution
#####################################

### Setting working directory, containing complete alignment PUs with a code "XX + number"
setwd ("working-directory") 

### Loading librarieslibrary(ape)
library(stringr)
library(ape)
library(phangorn)

### Loading function
source("Function_PBS_PU_Congruence.R")

### Loading input alignment
alignment = ("Alignment.fasta")

### Executing function
PBS_PU_Congruence(alignment, "EH", "PBS", "F84", 0.04, print.subtrees.fasta = "NO")


##FUNCTION PARAMETERS
#PBS_PU_Congruence <- function (nameOFfasta,codePU,codePBS, distmodel,limite1, print.subtrees = "YES", print.subtrees.fasta = "YES", print.singletons.fasta = "YES")
#IN THIS ORDER:  
#1. nameOFfasta # complete alignment, PBSs with a code "PBS+number), new PUs with a code "XX + number"
#2. codePU #letters before the PU number, all PUs should have the same letter e.g. "LG"
#3. codePBS #letters before the PBS number, all PBS should have the same letter e.g., "PBS"
#4. distmodel #e.g.  "F84". A character string specifying the evolutionary model to be used; must be one of "raw", "N", "TS", "TV", "JC69", "K80" (the default), "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", or "indelblock". For details see help dist.dna {ape}
#5. limite1 #number indicating divergent (ratio in one unit) per lineage. eg. 0.04 means 8% divergence
#6. print.subtrees.fasta = "YES"  # logic with info for printing or not .fasta files for each subtrees (could be yes or no)

