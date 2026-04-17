#######################################################
# Example of downloading some FliG and MgtE 
# amino acid (AA) sequences, from sequence IDs
#
# Then, saving to FASTA, and manipulating sequence labels
# 
#######################################################

install_cmds='
install.packages("ape")
install.packages("seqinr")
install.packages("rvest")
install.packages("rexpokit")
install.packages("cladoRcpp")
install.packages("rentrez")
install.packages("phytools")
install.packages("openxlsx")
install.packages("queryup")

library(devtools)
devtools::install_github(repo="nmatzke/BioGeoBEARS")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("msa")
BiocManager::install("ginmappeR")

' # END installs

library(ape)
library(msa)
library(seqinr)
library(BioGeoBEARS)
library(rvest) # for html_table()
library(phytools) # for midpoint.root
library(openxlsx) # for read.xlsx
library(ginmappeR)
library(queryup)
library(rentrez)  # for rentrez::entrez_fetch(db="protein", id=MotA_seq_ids[1:170], rettype="fasta")


# Copy to str2phy
#sourceall("/GitHub/bioinfRhints/Rsrc/")
#source("/GitHub/bioinfRhints/Rsrc/protein_bioinf_v1.R") # for firstword
sourceall("~/GitHub/AA_to_3Di/Rsrc/")

wd = "~/GitHub/AA_to_3Di/01_get_seqs/"
setwd(wd)


# Manual list of 34 protein sequence IDs:
# FliG = bacterial flagellum, flagellar switch protein (inner membrane)
# MgtE = magnesium transporter, inner membrane; thought to be a remote homolog of part of FliG
# Also listed in: 
# /ex/z_other/FliG_MgtE_seqids.txt
# /ex/z_other/FliG_MgtE_hmmer_hits_Seqnames.txt
seqids_34 = c("WP_160516493.1","ACD91087.1","BAD39524.1","ABF88562.1","CCB86228.1","AUX26777.1","QJR34130.1","AHZ85048.1","ABY36902.1","QUV78475.1","UOY14042.1","AAQ00836.1","EKT86237.1","WP_002722827.1","ADE84642.1","AJJ33030.1","QEO40368.1","AGW37977.1","CCW35060.1","CCG56443.1","ADB18917.1","UOY14499.1","BAM05993.1","ACI20682.1","ABC76852.1","AAC75006.1","ATU66743.1","BAH37372.1","AJJ31778.1","AIQ93516.1","UOY14169.1","UPH48283.1","AAP98818.1","WP_011385860.1")

# Retrieve the sequences to FASTA
downloaded_seqs = rentrez::entrez_fetch(db="protein", id=seqids_34, rettype="fasta")
outfn = "FliG_MgtEs_downloaded_seqs_to34.fasta"
write(x=downloaded_seqs, file=outfn, sep="")

# Look at sequences
tmpseqs = read_FASTA_safe(outfn, type="AA")
tmpseqs
numseqs = length(tmpseqs)
numseqs

# Convert from AAbin (binary storage of AAs) to character data
seqs = sapply(X=as.character(tmpseqs), FUN=paste0, collapse="")
seqs[1:2]
length(seqs)

# Sequence labels
seqnames = names(seqs)
seqnames

# Retrieve the sequence IDs
seqids_from_seqnames = get_seqids(strings=seqnames, split=" ")
seqids_from_seqnames

# Retrieve the species names
spnames = get_spnames(list_of_strings=seqnames, replace_spaces=TRUE)
spnames

# Retrive the rest (middle of the label)
protein_names = get_label_without_seqID_or_species(seqnames=seqnames, split=" ", recodes=genbank_prefixes(), removetxt=c(">"), replace_spaces=FALSE, prflag=FALSE)
protein_names

# See the counts of different labels
rev(sort(table(protein_names)))

# The user could devise new/simplified labels at this point, if desired.
# Usually, the best labels for use across multiple alignment and phylogeny programs are:
# seqid_StandardizedProteinName_Species
# ...with all spaces replaced by underscores ("_"), and all other special characters removed
