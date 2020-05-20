##########################################################################
# SCRIPT 1. Applying genetic divergence threshold to establish de novo LMDH or PBS
## Created by Carmelo Ándujar 2020
## Supporting material of Salces-Castellano et al. (2020) Climate drives community‐wide divergence within species over a limited spatial scale: evidence from an oceanic island
##########################################################################


#####################################
### Function PU_Congruence_from_alignment
#####################################

PU_Congruence_from_alignment <- function(nameOFfasta,codePU, distmodel,limite1, print.subtrees.fasta = "YES") 
  #nameOFfasta #complete alignment PUs with a code "XX + number"
  #codePU #letters before the PU number, all PUs should have the same letter e.g. "LG"
  #distmodel #e.g.  "F84". A character string specifying the evolutionary model to be used; must be one of "raw", "N", "TS", "TV", "JC69", "K80" (the default), "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", or "indelblock". For details see help dist.dna {ape}
  #limite1 #number indicating divergent (ratio in one unit) per lineage. eg. 0.04 means 8% divergence
  #print.subtrees.fasta = "YES"  # logic with info for printing or not .fasta files for each subtrees (could be yes or no)

{  
  
  #################################################################################################################################
  ########################         PART 1. Reading complete alignment and generating a list of parataxonomic      #################
  ########################         units (PU) and probable biological species (PBS) from names of sequences       #################
  #################################################################################################################################
  
  #Note: In our example, name of sequences has the pattern "TF" + number. This can be adjusted, or a list of PUs can be directly provided
