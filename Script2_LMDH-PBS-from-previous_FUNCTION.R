##########################################################################
# SCRIPT 2. Applying genetic divergence threshold to establish LMDH and PBS from previous ones
## Created by Carmelo Ándujar 2020
## Supporting material of Salces-Castellano et al. (2020) Climate drives community‐wide divergence within species over a limited spatial scale: evidence from an oceanic island
##########################################################################


#####################################
### Function PBS_PU_Congruence
#####################################

PBS_PU_Congruence <- function (nameOFfasta, codePU, codePBS, distmodel, limite1, print.subtrees.fasta="YES")
  #nameOFfasta # complete alignment, PBSs with a code "PBS+number), new PUs with a code "XX + number"
  #codePU #letters before the PU number, all PUs should have the same letter e.g. "LG"
  #codePBS #letters before the PBS number, all PBS should have the same letter e.g., "PBS"
  #distmodel #e.g.  "F84". A character string specifying the evolutionary model to be used; must be one of "raw", "N", "TS", "TV", "JC69", "K80" (the default), "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", or "indelblock". For details see help dist.dna {ape}
  #limite1 #number indicating divergent (ratio in one unit) per lineage. eg. 0.04 means 8% divergence
  #print.subtrees.fasta = "YES"  # logic with info for printing or not .fasta files for each subtrees (could be yes or no)
  
{  
  
  print.subtrees = "YES"
  print.singletons.fasta = "YES"
  
  FASTA = read.FASTA(nameOFfasta)
  fasta = read.table (nameOFfasta, header=FALSE) #Re-reading fasta
  
  
  #################################################################################################################################
  ########################         PART 1. Reading complete alignment and generating a list of parataxonomic      #################
  ########################         units (PU) and probable biological species (PBS) from names of sequences       #################
  #################################################################################################################################
 
  PU_fasta_specimens = str_extract(labels(FASTA), str_c(codePU,"\\d{1,10}")) #To extract PU codes from name of sequences for each sequence
  PBS_fasta_specimens = str_extract(labels(FASTA), str_c(codePBS,"\\d{1,10}")) #To extract PU codes from name of sequences for each sequence
  
  PU_fasta_specimens = c(PU_fasta_specimens, PBS_fasta_specimens)
  Tabla_PU_fasta = table(PU_fasta_specimens)  #To get a list of PUs
  Table_PU_values = Tabla_PU_fasta
  dim(Table_PU_values) = c(length(Tabla_PU_fasta),1) #To get a list of number of sequences(i.e. specimens) within each PU
  
  
  #################################################################################################################################
  ########################         PART 2. Making tree, and getting subtrees and singletons       #################################
  #################################################################################################################################
  
  ### PART 2a. Distances matrix and UPGMA tree
  
  dist_fasta = dist.dna(FASTA, model = distmodel, as.matrix = T)
  tree = upgma(dist_fasta)
  
  
  ### PART 2b. Subtrees and fasta files
  
  subt.t = subtrees(tree, wait=TRUE)
  tabla.all.limites = data.frame()
  lista = data.frame()
  for (i in 1:Nnode(tree))
  {
    if ((max(branching.times(subt.t[[i]])))<limite1)
    {
      for (x in 1:Ntip(subt.t[[i]]))
      {  
        caso = c("species", i, max(branching.times(subt.t[[i]])),(subt.t[[i]])[[3]][x],Ntip(subt.t[[i]]))
        dim(caso) = c(1,5)
        colnames(caso) = c("sp","spn","depth","taxa", "ntip")
        lista = rbind(lista,caso)
      }
    }
  }
  write.table(lista, file="lista.txt", row.names = FALSE, col.names= TRUE)  
  lista<-read.table ("lista.txt", header=TRUE)
  names = tree$tip.label
  maxnode = data.frame()
  for (j in 1:length(names))
  {
    subset <- which(lista$taxa == names[j])
    maximum <- max(lista$depth[subset])
    maxtips <- max(lista$ntip[subset])
    asignacion <- which(lista$taxa == names[j] & lista$depth == maximum & lista$ntip == maxtips)
    maxnode <- rbind(maxnode, lista[asignacion,])
  }
  
  #To generate a single file with subtree and fasta for each subtree
  #lista_species<-read.table ("lista_species.txt", header=TRUE)
  nodestab = table(maxnode[["spn"]])
  nodes_list = row.names(nodestab)
  all.used = data.frame()
  
  for (h in 1:length(nodes_list))
  {
    if (print.subtrees == "YES")
    {
      write.tree (subt.t[[(as.numeric(nodes_list[h]))]], file = str_c("node",h,"_lim",limite1, ".newick"))
    }
    subtreecase = subt.t[[(as.numeric(nodes_list[h]))]]
    subtreetips = subtreecase$tip.label #A table with two columns (name and seq) without header
    
    #Preparing fasta
    fasta.v = as.vector(fasta[,1])
    dim(fasta.v) = c(2,(length(fasta[,1])/2))
    fasta.v.t = t(fasta.v)
    seq_final = fasta.v.t
    seq_final[,1] 
    
    #Reading the list of external groups (GE) and generating fasta file without GE
    GE = subtreetips  #Table with a list of GE
    GE
    length(GE)
    
    id_seqs = data.frame()
    for (l in 1:length(GE))
    {
      id = which (seq_final[,1] == as.vector(str_c(">",as.vector(GE[l]))))
      id_seqs = c(id_seqs,id)       
    }
    id_seqs = as.numeric(id_seqs)
    id_seqs
    seq_final_sinGE = seq_final[id_seqs,]
    
    #To add data on a txt
    used.vector = rep(str_c("clade",h,"_limite",limite1),length(seq_final_sinGE[,1]))
    dim(used.vector) = c(length(seq_final_sinGE[,1]),1)
    namelist.used = cbind(seq_final_sinGE,used.vector)
    all.used = rbind(all.used, namelist.used)
    
    #Converting to fasta
    seq_fasta = data.frame()
    for (g in 1:length(seq_final_sinGE[,1]))
    {
      seq_name = as.vector(seq_final_sinGE[g,1])
      seq_fasta = c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
    }
    dim(seq_fasta)
    dim(seq_fasta) = c((2*length(seq_final_sinGE[,1])),1)
    
    #Printing the table in fasta format
    if (print.subtrees.fasta == "YES")
    {
      write.table(seq_fasta, file = str_c("node",h,"_lim",limite1, ".fas"), row.names = FALSE, col.names= FALSE) 
    }
  }
  
  #To generate a fasta file with sequences not included in subtrees (unused)
  id_seqs_unused = data.frame()
  for (k in 1:length(nodes_list))
  {
    #Preparingfasta
    fasta.v = as.vector(fasta[,1])
    dim(fasta.v) = c(2,(length(fasta[,1])/2))
    fasta.v.t = t(fasta.v)
    seq_final = fasta.v.t
    seq_final[,1] 
    
    subtreecase = subt.t[[(as.numeric(nodes_list[k]))]]
    subtreetips = subtreecase$tip.label
    
    id_seqs = data.frame()
    for (y in 1:length(subtreetips))
    {
      id = which (seq_final[,1] == as.vector(str_c(">",as.vector(subtreetips[y]))))
      id_seqs = c(id_seqs,id)        
    }
    id_seqs_unused = c(id_seqs_unused,id_seqs)
  }
  
  id_seqs_unused = as.numeric(id_seqs_unused)
  
  seq_final_sinGE = seq_final[-id_seqs_unused,]
  
  #To add data on a txt
  unused.vector = rep(str_c("unused_limite", limite1),length(seq_final_sinGE[,1]))
  dim(unused.vector) = c(length(seq_final_sinGE[,1]),1)
  namelist.unused = cbind(seq_final_sinGE,unused.vector)
  all.table.info.clade = rbind(all.used, namelist.unused)
  
  #Converting to fasta
  seq_fasta = data.frame()
  for (g in 1:length(seq_final_sinGE[,1]))
  {
    seq_name = as.vector(seq_final_sinGE[g,1])
    seq_fasta = c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
  }
  dim(seq_fasta)
  dim(seq_fasta) = c((2*length(seq_final_sinGE[,1])),1)
  
  #Printing the table in fasta format
  if (print.singletons.fasta == "YES")
  {
    write.table(seq_fasta, file = str_c("lim",limite1,"unused", ".fas"), row.names = FALSE, col.names= FALSE) 
  }
  
  #Creating final tables named "all.table.info.clade" and "all.table.info.cladeAZ"
  write.table(all.table.info.clade, file=str_c("all.table.info.clade","lim",limite1,".txt"))
  all.table.info.cladeAZ = read.table (str_c("all.table.info.clade","lim",limite1,".txt"), header=TRUE)
  ordenacion = order(all.table.info.cladeAZ[,1])
  all.table.info.cladeAZ = all.table.info.cladeAZ[ordenacion,]
  write.table(all.table.info.cladeAZ, file=str_c("all.table.info.clade_ordered","lim",limite1,".txt"))
  tabla.all.limites = all.table.info.cladeAZ
  
  
  ###########################################################################################################################################
  #############         PART 3. Running loop for each subtree with two or more taxon, and cheking for how many PUs are inside,    ###########
  #############         and for each PU inside, how many sequences considering total number of sequences within the PU            ###########
  ###########################################################################################################################################
  
  LISTA_COMPLETE = data.frame() #To generate the element where we will put all information at the end
  ficheros = dir(pattern = ".newick") #To get a vector with name of each file in the folder with ext .newick, which correspond to contig-trees obtained for each cluster at 8% with 3 or more sequences
  
  for (f in ficheros) #Loop to go through each contig-tree (cluster at 8% with 2 or more sequences)
  {
    tree = read.tree(f)
    subt.t = subtrees(tree, wait=TRUE)  #To get all subtrees within the contig-tree
    PU_specimens = str_extract(tree$tip.label, str_c(codePU,"\\d{1,10}")) #To get all PU codes for each tip of the contig-tree
    PBS_specimens = str_extract(tree$tip.label, str_c(codePBS,"\\d{1,10}"))
    PU_specimens = c(PU_specimens,PBS_specimens)
    Tabla_PU = table(PU_specimens) 
    PU_list = row.names(Tabla_PU)  #To generate a list with the PUs present in the contig-tree
    
    wdir = getwd()   #To store the path for the directory with all inputs
    d = str_c(f,"_d")
    dir.create(d)   #To generate a folder with the name of the contig-tree "FileName_d"
    localdir = str_c(wdir,"/",d)  
    setwd(localdir)  #Set working directory as the new forlder with the name of the contig-tree
    
    lista_def = data.frame() #To generate the element where we will put info for each contig-tree
    
    for (i in 1:length (PU_list)) #Loop to go throught each ENTITY (PU OR PBS) (in "PU_LIST") represented on each contig-tree F and generate a list with specimens (full names) of each PU
    #This list is used to find out the smallest subtree within the contig-tree including all the specimens of that PU on that contig-tree, checking if this is monphyletic or not. 
    {
      pos = which(row.names(Tabla_PU_fasta)==PU_list[i]) #Positions in initial table obtained from the complete fasta with all specimens that match a particular PU
      treesTips = data.frame() #To generate the element
      species = (str_c  ("\\.*",PU_list[i],"\\.*")) #To generate an string with the pateern of complete name in fasta and tips of contig-tree
      positions_species = which (str_detect(tree$tip.label, species) == TRUE) #To identify which tips (positions) correspond to a PU
      Species_list = tree$tip.label [positions_species] #List of tips for a PU
      for (n in 1: length (subt.t)) #Loop to go through every subtree for every PU represented on each contig-tree
      { 
        if  (length (which (str_detect(subt.t[[n]]$tip.label, species) == TRUE)) == length (positions_species)) #If the number of specimens for a PU is the same in the subtree than in all the contig-tree 
        #(i.e., the subtree includes all specimens for that PU in that contig-tree)
        {
          caso = c(n,subt.t[[n]]$Ntip)
          dim(caso) = c(1,2)
          treesTips = rbind(treesTips, caso)
          t = which (treesTips[,2]== min (as.vector(treesTips[,2])))
          Min_tree_number = as.vector(treesTips [t,1])
          Min_tree_tips = as.vector(treesTips [t,2])  #Lines to identify among all subtrees with all representatives from a particular PU, which is the smaller
        }
      } 
      
      #### If representatives for a PU within the contig-tree are monophyletic    
      if (Min_tree_tips == length (Species_list))  #it means that the representatives for a PU within the contig-tree are monophyletic
      {
        caso = c(f,length(tree$tip.label),PU_list[i],length(Species_list),Table_PU_values[pos,],Min_tree_tips,"yes",max(branching.times(subt.t[[Min_tree_number]])), "NA", "NA", "NA") #Combine a series of elements, named below
        dim(caso) = c(1,11)
        colnames(caso) = c("cluster","tamano_cluster","PU","n.specimens_in_cluster","total_PU_specimens","Min_tree_tips","congruence","max.branching.times", "inc_PU", "inc_specimens","indiv")
        lista_def = rbind(lista_def,caso)
        write.tree (subt.t[[Min_tree_number]], file = str_c(f, "_", PU_list[i],"_OK",".nwk") ) #To export the smallest tree with all PU representatives in a contig-tree 
      }
      
      #### If representatives for a PU within the contig-tree are NOT monophyletic i.e., there is no subtree with all and only representatives for that PU
      if (Min_tree_tips > length (Species_list)) #It means that the representatives for a PU within the contig-tree are NOT monophyletic
      {
        #To generate a list with specimens from others PU that are inlcuded in the smallest tree including all representatives for a PU within the contig-tree  
        positions_in_subtree = which (str_detect(subt.t[[Min_tree_number]]$tip.label, species) == FALSE)
        nconsistent = subt.t[[Min_tree_number]]$tip.label [positions_in_subtree]
        inc_specimens = "Inc. Specimens:"
        for (x in 1:length (nconsistent)) {
          inc_specimens = str_c (inc_specimens,nconsistent[x], sep =" ")}  
        
        #To generate a list the ENTITIES (PUs OR PBS) that are inlcuded in the smallest tree including all representatives for a PU within the contig-tree  
        PU_specimens_lineage = str_extract(nconsistent, str_c(codePU,"\\d{1,10}"))
        PBS_specimens_lineage = str_extract(nconsistent, str_c(codePBS,"\\d{1,10}"))
        PU_specimens_lineage = c(PU_specimens_lineage,PBS_specimens_lineage)
        Tabla_PU_lineage = table (PU_specimens_lineage)
        PU_list_lineage = row.names(Tabla_PU_lineage)
        inc_PU = "Inc. PU:"
        for (x in 1:length (PU_list_lineage) ) {
          inc_PU = str_c (inc_PU,PU_list_lineage[x], sep =" ")} 
        
        #To write all relevant information into the final table
        caso = c(f,length(tree$tip.label),PU_list[i],length(Species_list),Table_PU_values[pos,],Min_tree_tips,"NO",max(branching.times(subt.t[[Min_tree_number]])), inc_PU, inc_specimens, "NA") 
        dim(caso) = c(1,11)
        colnames(caso) = c("cluster","tamano_cluster","PU","n.specimens_in_cluster","total_PU_specimens","Min_tree_tips","congruence","max.branching.times", "inc_PU", "inc_specimens", "indiv")
        lista_def = rbind(lista_def,caso)
        write.tree(subt.t[[Min_tree_number]], file = str_c(f,"_",PU_list[i],"_FAIL",".nwk") )   #To export the smallest tree with all PU representatives in a contig-tree 
      }
    }   
    
    #To export tables for each contig-tree inside contig-tree folders
    Table_sort = lista_def [,1:8]
    Table_Inc_PU = lista_def [,1:9]
    Table_Inc_PU_Inc_specimens = lista_def
    
    write.table(Table_sort, file = str_c(f,"_Table_sort.txt"), row.names = FALSE, col.names= TRUE)
    write.table(Table_Inc_PU, file = str_c(f,"_Table_Inc_PU.txt"), row.names = FALSE, col.names= TRUE)
    write.table(Table_Inc_PU_Inc_specimens, file = str_c(f,"_Table_Inc_PU_Inc_specimens.txt"), row.names = FALSE, col.names= TRUE)
    
    #To generate final table with info for all contig-trees 
    LISTA_COMPLETE = rbind(LISTA_COMPLETE,lista_def)
    setwd(wdir)
  }
  
  #To export final table with info for all contig-trees
  write.table(LISTA_COMPLETE, file= "TABLA_allSubtrees.txt", row.names = FALSE, col.names= TRUE)
  
  
  #################################################################################################################################
  ########################         PART 4. Working with singletons                   ##############################################
  #################################################################################################################################
  
  singletons = read.FASTA(str_c("lim",limite1,"unused", ".fas"))
  singletons_PUs = str_extract(labels(singletons), str_c(codePU,"[A-B]")) #To extract PU codes from name of sequences for each sequence
  singletons_PBS = str_extract(labels(singletons), str_c(codePBS,"\\d{1,10}")) #To extract PBS codes from name of sequences for each sequence
  singletons_PUs = c(singletons_PUs,singletons_PBS) #To combine PUs and PBS
  Tabla_singletons_PUs = table (singletons_PUs) #To get a list of PUs and PBS
  Tabla_singletons_PUs_values = Tabla_singletons_PUs
  dim(Tabla_singletons_PUs_values) = c(length(Tabla_singletons_PUs),1) #To get a list of number of sequences(i.e. specimens) within each PU
  
  for (nn in 1:length(row.names(Tabla_singletons_PUs)))
  {
    pos = which(row.names(Tabla_PU_fasta)==row.names(Tabla_singletons_PUs)[nn])
    caso = c("singleton",1,row.names(Tabla_singletons_PUs)[nn],Tabla_singletons_PUs_values[nn,],Table_PU_values[pos,],"NA","NA","NA", "NA", "NA", "NA") 
    dim(caso) = c(1,11)
    colnames(caso) = c("cluster","tamano_cluster","PU","n.specimens_in_cluster","total_PU_specimens","Min_tree_tips","congruence","max.branching.times", "inc_PU", "inc_specimens", "indiv")
    LISTA_COMPLETE = rbind(LISTA_COMPLETE,caso) 
  }
  
  #To export final table with info for all contig-trees
  write.table(LISTA_COMPLETE, file= "TABLA_Subtrees_singletons.txt", row.names = FALSE, col.names= TRUE)
  
}
