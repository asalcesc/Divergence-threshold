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

  read.FASTA(nameOFfasta)->FASTA
  #read.table(nameOFfasta, header=FALSE)->fasta #re-reading fasta
  str_extract(labels(FASTA), str_c(codePU,"\\d{1,10}"))->TF_fasta_specimens #To extract PU codes from name of sequences for each sequence
  table (TF_fasta_specimens)->Tabla_TF_fasta  #To get a list of PUs
  Tabla_TF_fasta->Table_TF_values
  dim(Table_TF_values)<-c(length(Tabla_TF_fasta),1) #To get a list of number of sequences (i.e. specimens) within each PU


  #################################################################################################################################
  ########################         PART 2. Making tree, and getting subtrees and singletons       #################################
  #################################################################################################################################
  
  library(ape)
  library(phangorn)
  library(seqinr)

  ### PART 2a. Distances matrix and UPGMA tree

  dist.dna(FASTA, model = distmodel, as.matrix = T)->dist_fasta
  tree <- upgma(dist_fasta)

  ### PART 2b. Subtrees and fasta files

  subtrees(tree, wait=TRUE)->subt.t
  tabla.all.limites <-data.frame()
  lista<-data.frame()

  for (i in 1:Nnode(tree))
  {
    if ((max(branching.times(subt.t[[i]])))<limite1)
    {
      for (x in 1:Ntip(subt.t[[i]]))
      {  
        caso<-c("species", i, max(branching.times(subt.t[[i]])),(subt.t[[i]])[[3]][x],Ntip(subt.t[[i]]))
        dim(caso)<-c(1,5)
        colnames(caso)<-c("sp","spn","depth","taxa", "ntip")
        lista<-rbind(lista,caso)
      }
    }
  }
  write.table(lista, file="lista.txt", row.names = FALSE, col.names= TRUE)  

  lista<-read.table ("lista.txt", header=TRUE)
  tree$tip.label->names
  data.frame()->maxnode
  for (j in 1:length(names))
  {
    subset<-which (lista$taxa == names[j])
    maximum<-max(lista$depth[subset])
    maxtips<-max(lista$ntip[subset])
    asignacion<-which(lista$taxa == names[j] & lista$depth == maximum & lista$ntip == maxtips)
    maxnode<-rbind(maxnode, lista[asignacion,])
  }
  write.table(maxnode, file="lista_species.txt", row.names = FALSE, col.names= TRUE)

  lista_species<-read.table ("lista_species.txt", header=TRUE)
  table(lista_species[["spn"]])->nodestab
  row.names(nodestab)->nodes_list

  #To generate a single file with subtree and fasta for each subtree
  all.used<-data.frame()
  for (h in 1:length(nodes_list))
  {
    if (print.subtrees == "YES")
    {
      write.tree(subt.t[[(as.numeric(nodes_list[h]))]], file = str_c("node",h,"_lim",limite1, ".newick") )
    }
    subt.t[[(as.numeric(nodes_list[h]))]]->subtreecase
    subtreecase$tip.label->subtreetips #A table with two columns (name and seq) without header
  
    #Preparing fasta
    as.vector(fasta[,1])->fasta.v
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
  
    #Reading the list of external groups (GE) and generating fasta file without GE
    GE<-subtreetips  #Table with a list of GE
    GE
    length(GE)
    id_seqs<-data.frame()
    for (l in 1:length(GE))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(GE[l]))))->id
      id_seqs <-c(id_seqs,id)       
    }
    as.numeric(id_seqs)->id_seqs
    id_seqs
    seq_final_sinGE<-seq_final[id_seqs,]
  
    #To add data on a txt
    rep(str_c("clade",h,"_limite",limite1),length(seq_final_sinGE[,1]))->used.vector
    dim(used.vector)<-c(length(seq_final_sinGE[,1]),1)
    cbind(seq_final_sinGE,used.vector)->namelist.used
    all.used<-rbind(all.used, namelist.used)
  
    #Converting to fasta  
    seq_fasta<-data.frame()
    for (g in 1:length(seq_final_sinGE[,1]))
    {
      seq_name<-as.vector(seq_final_sinGE[g,1])
      seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
    }
    dim(seq_fasta)
    dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
  
    #Printing the table in fasta format
    if (print.subtrees.fasta == "YES")
    {
      write.table(seq_fasta, file = str_c("node",h,"_lim",limite1, ".fas"), row.names = FALSE, col.names= FALSE) 
    }
    #write.table(namelist.used, file=str_c("clade",h,"_limite",limite1,".txt"))
  }

  #To generate a fasta file with sequences not included in subtrees (unused)
  id_seqs_unused<-data.frame()
  for (k in 1:length(nodes_list))
  {
    #Preparing fasta
    as.vector(fasta[,1])->fasta.v
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
  
    subt.t[[(as.numeric(nodes_list[k]))]]->subtreecase
    subtreecase$tip.label->subtreetips
  
    id_seqs<-data.frame()
    for (y in 1:length(subtreetips))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(subtreetips[y]))))->id
      id_seqs <-c(id_seqs,id)        
    }
    id_seqs_unused <-c(id_seqs_unused,id_seqs)
  }

  as.numeric(id_seqs_unused)->id_seqs_unused
  id_seqs
  seq_final_sinGE<-seq_final[-id_seqs_unused,]

  #To add data on a txt
  rep(str_c("unused_limite", limite1),length(seq_final_sinGE[,1]))->unused.vector
  dim(unused.vector)<-c(length(seq_final_sinGE[,1]),1)
  cbind(seq_final_sinGE,unused.vector)->namelist.unused
  rbind(all.used, namelist.unused)->all.table.info.clade

  #Converting to fasta
  seq_fasta<-data.frame()
  for (g in 1:length(seq_final_sinGE[,1]))
  {
    seq_name<-as.vector(seq_final_sinGE[g,1])
    seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
  }
  dim(seq_fasta)
  dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)

  #Printing the table in fasta format
  if (print.singletons.fasta == "YES")
  {
    write.table(seq_fasta, file = str_c("lim",limite1,"unused", ".fas"), row.names = FALSE, col.names= FALSE) 
  }

  #Creating final tables named "all.table.info.clade" and "all.table.info.cladeAZ"
  write.table(all.table.info.clade, file=str_c("all.table.info.clade","lim",limite1,".txt"))
  read.table(str_c("all.table.info.clade","lim",limite1,".txt"), header=TRUE)->all.table.info.cladeAZ
  ordenacion=order(all.table.info.cladeAZ[,1])#Ordering data
  all.table.info.cladeAZ<-all.table.info.cladeAZ[ordenacion,]
  write.table(all.table.info.cladeAZ, file=str_c("all.table.info.clade_ordered","lim",limite1,".txt"))
  tabla.all.limites<-all.table.info.cladeAZ


  ###########################################################################################################################################
  #############         PART 3. Running loop for each subtree with two or more taxon, and cheking for how many PUs are inside,    ###########
  #############         and for each PU inside, how many sequences considering total number of sequences within the PU            ###########
  ###########################################################################################################################################
  
  LISTA_COMPLETE <- data.frame() #To generate the element where we will put all information at the end
  ficheros <- dir(pattern = ".newick") #To get a vector with name of each file in the folder with ext .newick, which correspond to contig-trees obtained for each cluster at 8% with 3 or more sequences

  for (f in ficheros) #Loop to go through each contig-tree (cluster at 8% with 3 or more sequences)
  {
    read.tree(f)->tree
    subtrees(tree, wait=TRUE)->subt.t  #To get all subtrees within the contig-tree
    str_extract(tree$tip.label, str_c(codePU,"\\d{1,10}"))->TF_specimens #To get all PU codes for each tip of the contig-tree
    table (TF_specimens)->Tabla_TF  
    row.names(Tabla_TF)->TF_list  #To generate a list with the PUs present in the contig-tree
  
    getwd()->wdir   #To store the path for the directory with all inputs
    str_c (f,"_d")->d
    dir.create(d)   #To generate a folder with the name of the contig-tree "FileName_d"
    str_c(wdir,"/",d)->localdir  
    setwd(localdir)  #Set working directory as the new forlder with the name of the contig-tree
  
    lista_def <- data.frame() #To generate the element where we will put info for each contig-tree
  
    for (i in 1:length (TF_list)) #Loop to go thorugh each PU (in "TF_LIST") represented on each contig-tree F and generate a list with specimens (full names) of each PU
    #This list is used to find out the smallest subtree within the contig-tree including all the specimens of that PU on that contig-tree, checking if this is monphyletic or not. 
    {
      which(row.names(Tabla_TF_fasta)==TF_list[i])->pos #Positions in initial table obtained from the complete fasta with all specimens that match a particular PU
      treesTips <- data.frame() #To generate the element
      species<-(str_c  ("\\.*",TF_list[i],"\\.*")) #To generate an string with the pateern of complete name in fasta and tips of contig-tree
    
      which (str_detect(tree$tip.label, species) == TRUE)->positions_species #To identify which tips (positions) correspond to a PU
      Species_list<- tree$tip.label [positions_species] #List of tips for a PU
    
      for (n in 1: length (subt.t)) #Loop to go through every subtree for every PU represented on each contig-tree
      { 
        if  (length (which (str_detect(subt.t[[n]]$tip.label, species) == TRUE)) == length (positions_species)) #If the number of specimens for a PU is the same in the subtree than in all the contig-tree 
        #(i.e., the subtree includes all specimens for that PU in that contig-tree)
        {
          caso <-c(n,subt.t[[n]]$Ntip)
          dim(caso)<-c(1,2)
          treesTips<- rbind(treesTips, caso)
        
          which (treesTips[,2]  == min (as.vector(treesTips[,2])))->t
          as.vector(treesTips [t,1])->Min_tree_number
          as.vector(treesTips [t,2])->Min_tree_tips  #Lines to identify among all subtrees with all representatives from a particular PU, which is the smaller
        }
      } 

      #### If representatives for a PU within the contig-tree are monophyletic    
      if (Min_tree_tips == length (Species_list))  #It means that the representatives for a PU within the contig-tree are monophyletic
      {
        caso<-c(f,length(tree$tip.label),TF_list[i],length (Species_list),Table_TF_values[pos,],Min_tree_tips,"yes",max(branching.times(subt.t[[Min_tree_number]])), "NA", "NA") #Combine a series of elements, named below
        dim(caso)<-c(1,10)
        colnames(caso)<-c("cluster","tamano_cluster","PU","n_specimens_in_cluster","total_PU_specimens","Min_tree_tips","congruence","max_branching_times", "inc_PU", "inc_specimens")
        lista_def<-rbind(lista_def,caso)
        write.tree (subt.t[[Min_tree_number]], file = str_c(f, "_", TF_list[i],"_OK",".nwk") ) #To export the smallest tree with all PU representatives in a contig-tree 
      }

      #### If representatives for a PU within the contig-tree are NOT monophyletic i.e., there is no subtree with all and only representatives for that PU
      if (Min_tree_tips > length (Species_list)) ##it means that the representatives for a PU within the contig-tree are NOT monophyletic
      {
        #To generate a list with specimens from others PU that are included in the smallest tree including all representatives for a PU within the contig-tree  
        which (str_detect(subt.t[[Min_tree_number]]$tip.label, species) == FALSE)->positions_in_subtree
        nconsistent <-subt.t[[Min_tree_number]]$tip.label [positions_in_subtree]
        inc_specimens <- "Inc. Specimens:"
        for (x in 1:length (nconsistent) ) {
          inc_specimens <- str_c (inc_specimens,nconsistent[x], sep =" ")}  
      
        #To generate a list the PUs that are inlcuded in the smallest tree including all representatives for a PU within the contig-tree
        str_extract(nconsistent, str_c(codePU,"\\d{1,10}"))->TF_specimens_lineage
        table (TF_specimens_lineage)->Tabla_TF_lineage
        row.names(Tabla_TF_lineage)->TF_list_lineage
        inc_PU<- "Inc. PU:"
        for (x in 1:length (TF_list_lineage) ) {
          inc_PU <- str_c (inc_PU,TF_list_lineage[x], sep =" ")} 
      
        #To write all relevant information into the final table
        caso<-c(f,length(tree$tip.label),TF_list[i],length (Species_list),Table_TF_values[pos,],Min_tree_tips,"NO",max(branching.times(subt.t[[Min_tree_number]])), inc_PU, inc_specimens) 
        dim(caso)<-c(1,10)
        colnames(caso)<-c("cluster","tamano_cluster","PU","n_specimens_in_cluster","total_PU_specimens","Min_tree_tips","congruence","max_branching_times", "inc_PU", "inc_specimens")
        lista_def<-rbind(lista_def,caso)
        write.tree (subt.t[[Min_tree_number]], file = str_c(f,"_",TF_list[i],"_FAIL",".nwk") )   #To export the smallest tree with all PU representatives in a contig-tree 
      }
    }   
  
    #To export tables for each contig-tree inside contig-tree folders
    Table_sort<- lista_def [,1:8]
    Table_Inc_PU<- lista_def [,1:9]
    Table_Inc_PU_Inc_specimens<- lista_def
  
    write.table(Table_sort, file = str_c(f,"_Table_sort.txt"), row.names = FALSE, col.names= TRUE)
    write.table(Table_Inc_PU, file= str_c(f,"_Table_Inc_PU.txt"), row.names = FALSE, col.names= TRUE)
    write.table(Table_Inc_PU_Inc_specimens, file= str_c(f,"_Table_Inc_PU_Inc_specimens.txt"), row.names = FALSE, col.names= TRUE)
  
    #To generate final table with info for all contig-trees 
    LISTA_COMPLETE<-rbind(LISTA_COMPLETE,lista_def)
    setwd(wdir)
  }

  #To export final table with info for all contig-trees
  write.table(LISTA_COMPLETE, file= "TABLA_CONTIGSmayor2.txt", row.names = FALSE, col.names= TRUE)


  #################################################################################################################################
  ########################         PART 4. Working with singletons                                #################################
  #################################################################################################################################

  read.FASTA(str_c("lim",limite1,"unused", ".fas"))->singletons

  str_extract(labels(singletons), str_c(codePU,"\\d{1,10}"))->singletons_PUs #To extract PU codes from name of sequences for each sequence

  table (singletons_PUs)->Tabla_singletons_PUs  #To get a list of PUs
  Tabla_singletons_PUs->Tabla_singletons_PUs_values
  dim(Tabla_singletons_PUs_values)<-c(length(Tabla_singletons_PUs),1) #To get a list of number of sequences(i.e. specimens) within each PU

  for (nn in 1:length(row.names(Tabla_singletons_PUs)))
  {
    which(row.names(Tabla_TF_fasta)==row.names(Tabla_singletons_PUs)[nn])->pos
    caso<-c("singleton",1,row.names(Tabla_singletons_PUs)[nn],Tabla_singletons_PUs_values[nn,],Table_TF_values[pos,],"NA","NA","NA", "NA", "NA") 
    dim(caso)<-c(1,10)
    colnames(caso)<-c("cluster","tamano_cluster","PU","n_specimens_in_cluster","total_PU_specimens","Min_tree_tips","congruence","max_branching_times", "inc_PU", "inc_specimens")
    LISTA_COMPLETE<-rbind(LISTA_COMPLETE,caso) 
  }

  #To export final table with info for all contig-trees
  write.table(LISTA_COMPLETE, file= "TABLA_subtrees_and_singletons.txt", row.names = FALSE, col.names= TRUE)
  
}

