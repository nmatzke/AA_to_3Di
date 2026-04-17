
#######################################################
# Using the downloaded raw sequences, I will run ChimeraX via script to get the AlphaFolds
#######################################################

install_cmds='
install.packages("ape")
install.packages("seqinr")
install.packages("rvest")
install.packages("rexpokit")
install.packages("cladoRcpp")
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

# Just for "linecount", not very important
#source("~/GitHub/bioinfRhints/R/_genericR_v3.R")

wd = "~/GitHub/AA_to_3Di/02_get_alphafolds/"
setwd(wd)

fasta_fn = "FliG_MgtEs_downloaded_seqs_to34.fasta"

tmpseqs = read_FASTA_safe(fasta_fn, type="AA")
numseqs = length(tmpseqs)
numseqs

# Look at sequences, as strings
seqs = sapply(X=as.character(tmpseqs), FUN=paste0, collapse="")
seqs[1:2]
length(seqs)

#######################################################
# Code to run ChimeraX from R
# *if* you have everything installed
#######################################################
outdfs = NULL; i=1

# outdfs = outdfs[1:124,]
#exclude_is = c(125,152, 159, 160)
exclude_is = c()  # You can skip certain sequences if you determine they are problematic
for (i in 1:numseqs)
	{
	txt = paste0("Getting alphafold structure for sequence #", i, "/", numseqs)
	cat("\n")
	cat(txt)
	cat("\n")
	
	if ((i %in% exclude_is) == TRUE)
		{
		seqs = ape::read.FASTA(fasta_fn, type="AA")
		header = names(seqs)[[i]]
		gid = firstword(header)
		
		txt = paste0("Sequence #", i, "/", numseqs, ", aka '", gid, "', skipped as it is known to cause a problem.")
		cat("\n")
		cat(txt)
		cat("\n")
		outdf = rep(NA, times=13)
		names(outdf) = c("seqnum", "AlphaFold Model", "Query Sequence", "Identity %", "Coverage %", "cif_fn", "html_fn", "table_fn", "dump_fn", "chain_names_fn", "di3_fn", "aa_fn", "both_fn")
		outdfs = rbind(outdfs, outdf)
		next()
		} # END if (i %in% exclude_is)
	
	# Just get the CIF files
	outdf = fasta_to_alphafold_cifs(fasta_fn=fasta_fn, seqnum=i, chimeraX_install=2025)
	#tmpdf = fasta_to_alphafold_cifs(fasta_fn=fasta_fn, seqnum=i, chimeraX_install=2024)
	print(outdf[1:5])
	outdfs = rbind(outdfs, outdf)
	
	
	
	outdfs_fn = "seqs_alphafolds.txt"
	write.table(x=outdfs, file=outdfs_fn, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	} # END for (i in 1:numseqs)



# Extract the information from the HTML tables
source("~/GitHub/AA_to_3Di/Rsrc/str2phy_v1.R")

outdfs2_fn = "seqs_alphafolds.txt"
outdfs2 = read.table(outdfs2_fn, quote="", sep="\t", header=TRUE, skip=0)
outdfs2
names_outdfs2 = c("seqnum", "AlphaFold Model", "Query Sequence", "Identity %", "Coverage %", "cif_fn", "html_fn", "table_fn", "dump_fn", "chain_names_fn", "di3_fn", "aa_fn", "both_fn")
names(outdfs2) = names_outdfs2
outdfs2
alphafold_dfs = NULL

# Re-assemble the table of HTML results, giving "Identity %" and "Coverage %" of each hit
for (i in 1:numseqs)
	{
	txt = paste0("Reading ChimeraX --alphafold table returned for sequence #", i, "/", numseqs)
	cat("\n")
	cat(txt)
	
	html_fn = outdfs2$html_fn[i]
	table_fn = outdfs2$table_fn[i]
	
	alphafold_df = alphafold_html_to_df(html_fn, table_fn, chimeraX_install=2025)
	alphafold_df
	
	alphafold_dfs = rbind(alphafold_dfs, alphafold_df)
	}
cat("\n...done\n")
alphafold_dfs
dim(outdfs2)


# Write over, if needed
outdfs3 = outdfs2
outdfs3[,"AlphaFold Model"] = alphafold_dfs[,"AlphaFold Model"]
outdfs3[,"Query Sequence"] = alphafold_dfs[,"Query Sequence"]
outdfs3[,"Identity %"] = alphafold_dfs[,"Identity %"]
outdfs3[,"Coverage %"] = alphafold_dfs[,"Coverage %"]
outdfs3 = dfnums_to_numeric(outdfs3)

outdfs3_fn = "seqs_alphafolds3.txt"
write.table(x=outdfs3, file=outdfs3_fn, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# Look at the output table
outdfs3



# Weak or no hits: 1 sequences
TF = outdfs3[,"Identity %"] < 80
sum(TF)
outdfs3_with_weak_hits = outdfs3[TF,]
outdfs3_with_weak_hits[,1:5]
dim(outdfs3_with_weak_hits)
# The first sequence was a weaker hit; user must decide what to do here;
# these might we worth doing by hand in Google CoLabFold
outdfs3_with_weak_hits$seqnum

# Look at the weak-hit sequence
seqs[outdfs3_with_weak_hits$seqnum]

outdfs3_with_weak_hits_fn = "outdfs3_with_weak_hits.fasta"
ape::write.FASTA(tmpseqs[outdfs3_with_weak_hits$seqnum], file=outdfs3_with_weak_hits_fn)
moref(outdfs3_with_weak_hits_fn)



# The user can decide to replace any AlphaFold structure ".cif"
# files by hand, e.g. with .cif files derived from novel
# AlphaFold structure predictions from Google colabfold.


# Once the .cif structure files are all available in:
# /GitHub/AA_to_3Di/02_get_alphafolds
# ...generate the 3Di characters with the 
# foldseek structureto3didescriptor command,
# run via R functions
#######################################################
# Convert alphafold structures to 3di & associated files
#######################################################
str_fns = outdfs3$cif_fn
fns_dfs = NULL
for (i in 1:length(str_fns))
	{
	fns_df = str_to_3di(str_fns[i], suffix="\\.cif")
	fns_dfs = rbind(fns_dfs, fns_df)
	}

fns_dfs





#######################################################
# YOU CAN STOP HERE; ADDITIONAL COMMANDS ARE IN
# TERMINAL TO MAKE AND RUN ALIGNMENTS ETC, THERE 
# ARE SEVERAL OPTIONS
#######################################################


 




#######################################################
# Align the 3di sequences, using famsa3di
#######################################################
cmds='
cd /GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/

famsa3di seqs_3di.fasta seqs_3di_famsa3diAligned.fasta
' # END cmds


# Align the AAs to the famsa3di alignment
dumps_wd = "~/GitHub/AA_to_3Di/ex/FliG_MgtE/2025-04-29_FliG_MgtE/02_get_alphafolds/"
wd = "~/GitHub/AA_to_3Di/ex/FliG_MgtE/2025-04-29_FliG_MgtE/03_famsa3di/"
setwd(wd)
source("~/GitHub/AA_to_3Di/Rsrc/str2phy_v1.R")


# Get all the AA and 3di sequences from the dump files in the directory
fns = list.files(path=dumps_wd, pattern="*.dump", full.names=TRUE)
TF = base::endsWith(fns, suffix=".dump")
fns = fns[TF]
fns = slashslash(fns)
fns

# Get all the AA and 3di sequences from the dump files in the directory
outfn_AAs = "seqs_AAs_fromStructs.fasta"
get_AA_FASTA_from_foldseek_dumps(fns, outfn=outfn_AAs)

outfn_3dis = "seqs_3dis_fromStructs.fasta"
get_3di_FASTA_from_foldseek_dumps(fns, outfn=outfn_3dis)


#######################################################
# Align with famsa3di
#######################################################
cmds='
cd ~/GitHub/AA_to_3Di/ex/FliG_MgtE/2025-04-29_FliG_MgtE/03_famsa3di/
ls
famsa3di seqs_3dis_fromStructs.fasta seqs_3di_famsa3diAligned.fasta

'

# Align the AAs against the famsa3di alignment of 3dis
wd = "~/GitHub/AA_to_3Di/ex/FliG_MgtE/2025-04-29_FliG_MgtE/03_famsa3di/"
setwd(wd)

aligned_fn = "seqs_3di_famsa3diAligned.fasta"
unaligned_fn = outfn_AAs
outfn = "seqs_AA_famsa3diAligned.fasta"
aligned_AAs = align_3dis_to_AAs(aligned_fn, unaligned_fn, outfn)
moref(outfn)






#######################################################
# Trim the alignments to just columns with >10% complete data
#######################################################
cmds='
cd ~/GitHub/AA_to_3Di/ex/FliG_MgtE/2025-04-29_FliG_MgtE/03_famsa3di/

trimal -in seqs_AA_famsa3diAligned.fasta -out seqs_AA_famsa3diAligned_trim05.fasta -gt 0.05 -colnumbering | tee seqs_AA_famsa3diAligned_trim05_cols.txt

trimal -in seqs_3di_famsa3diAligned.fasta -out seqs_3di_famsa3diAligned_trim05.fasta -gt 0.05 -colnumbering | tee seqs_3di_famsa3diAligned_trim05_cols.txt
' # END cmds




# Foldmason alignment of structures
cmds='
cd ~/GitHub/AA_to_3Di/ex/FliG_MgtE/2025-04-29_FliG_MgtE/04_foldmason

foldmason easy-msa ~/GitHub/AA_to_3Di/ex/FliG_MgtE/2025-04-29_FliG_MgtE/02_get_alphafolds/*.cif FliG_MgtE_cifs_FMalign.fasta tmpFolder --report-mode 1
'



#######################################################
# Reorder the famsa3di alignments to match the foldmason alignment order
#######################################################
wd = "~/GitHub/AA_to_3Di/ex/FliG_MgtE/2025-04-29_FliG_MgtE/03_famsa3di/"
setwd(wd)

aln_with_good_order_fn = "~/GitHub/AA_to_3Di/ex/FliG_MgtE/2025-04-29_FliG_MgtE/04_foldmason/FliG_MgtE_cifs_FMalign.fasta_aa.fa"
aln_with_good_order = read_FASTA_safe(aln_with_good_order_fn, type="AA")
tip_txt = names(aln_with_good_order)
tip_txt
AA_fn = "~/GitHub/AA_to_3Di/ex/FliG_MgtE/2025-04-29_FliG_MgtE/03_famsa3di/seqs_AA_famsa3diAligned.fasta"
output = reorder_fasta(fasta_fn=AA_fn, tip_txt, outfn=NULL, type="AA")

di3_fn = "~/GitHub/AA_to_3Di/ex/FliG_MgtE/2025-04-29_FliG_MgtE/03_famsa3di/seqs_3di_famsa3diAligned.fasta"
output = reorder_fasta(fasta_fn=di3_fn, tip_txt, outfn=NULL, type="AA")

AA_fn = "~/GitHub/AA_to_3Di/ex/FliG_MgtE/2025-04-29_FliG_MgtE/03_famsa3di/seqs_AA_famsa3diAligned_trim05.fasta"
output = reorder_fasta(fasta_fn=AA_fn, tip_txt, outfn=NULL, type="AA")

di3_fn = "~/GitHub/AA_to_3Di/ex/FliG_MgtE/2025-04-29_FliG_MgtE/03_famsa3di/seqs_3di_famsa3diAligned_trim05.fasta"
output = reorder_fasta(fasta_fn=di3_fn, tip_txt, outfn=NULL, type="AA")



#######################################################
#######################################################
# NJM says:
# NICK STOPPED HERE - 2025-04-29
# Files at: /GitHub/AA_to_3Di/ex/FliG_MgtE/2025-04-29_FliG_MgtE/
#######################################################
#######################################################













#######################################################
# Terminal commands for using USalign to do a structural alignment 
# using just the 3D-structure files, *not 3di*.
#######################################################

terminal_cmds='
#outdfs_fn = "3281_BRDs_alphafolds_WORKED_BKUP.txt"
#write.table(x=outdfs, file=outdfs_fn, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


# -mm 4         -- 4: alignment of multiple monomeric chains into a consensus alignment
# -ter 2        -- 2: (default) only align the first chain
# -TMscore 0 	  -- 0: (default) sequence independent structure alignment

cd /GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/
mkdir chains

# Put all the matching filenames into a text file
cp *.cif chains
ls chains/*.cif > list_of_chains.txt

# Add :A to each line, as these are all A chains
#sed -i '' 's/$/:A/'  list_of_chains.txt
# Remove chains/chains/ to chains/
sed -i '' 's/chains\///'  list_of_chains.txt
cp list_of_chains.txt list_of_chains_orig.txt
wc -l list_of_chains_orig.txt
head -18 list_of_chains_orig.txt > list_of_chains.txt
wc -l list_of_chains.txt
head list_of_chains.txt

tail list_of_chains.txt

USalign -dir chains/ list_of_chains.txt -mol prot -mm 4 | tee usalign_test1.txt &

more usalign_test1.txt
wc -l usalign_test1.txt  # 51 lines in file

head -13 usalign_test1.txt # 13 lines in header
wc -l usalign_test1.txt

# Get file minus first 13 lines
tail -38 usalign_test1.txt > usalign_test2.fasta

head usalign_test2.fasta
tail usalign_test2.fasta
wc -l usalign_test2.fasta

# Keep just first 36 lines
head -36 usalign_test2.fasta > usalign_test3.fasta
tail usalign_test3.fasta
wc -l usalign_test3.fasta  # 36, which equals 18 * 2
' # END terminal commands



#######################################################
# Align some 3di sequences to AA sequences already aligned
# (eg vs USalign structural alignment)
#######################################################

# 1. Read in a alignment (eg from USalign)
# 2. Line up the 3di sequences to that

# Convert an aa fasta file to a 3di fasta file
library(ape)
library(seqinr)
library(BioGeoBEARS)
sourceall("/GitHub/bioinfRhints/Rsrc")
source("/GitHub/AA_to_3Di/Rsrc/str2phy_v1.R")


wd = "/GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/"
setwd(wd)
list.files()

# Unaligned fasta
fasta_fn = "seqs.fasta"

# Aligned fasta
usaln_fn = "usalign_test3.fasta"
di3s_dir = "/GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/"
di3_outfn = "seqs_3di_aln.fasta"

di3_alignment = align_3dis(usaln_fn, di3s_dir, di3_outfn, pattern_3di_fn="*_3di.fasta")
di3_alignment












#######################################################
# Trim the alignments to just columns with >35% complete data
#######################################################
cmds='
cd /GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/

trimal -in seqs_AA_famsa3diAligned.fasta -out seqs_AA_famsa3diAligned_trim35.fasta -gt 0.35 -colnumbering | tee seqs_AA_famsa3diAligned_trim35_cols.txt

trimal -in seqs_3di_famsa3diAligned.fasta -out seqs_3di_famsa3diAligned_trim35.fasta -gt 0.35 -colnumbering | tee seqs_3di_famsa3diAligned_trim35_cols.txt
' # END cmds


#######################################################
# Horizontally concatenate 2 alignment files
#######################################################
fn1 = "/GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/seqs_AA_famsa3diAligned_trim35.fasta"
fn2 = "/GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/seqs_3di_famsa3diAligned_trim35.fasta"
fn3 = "seqs_AA3di_famsa3diAligned.fasta"
merged_df = hcat_align_fns(fn1, fn2, fn3)
moref(fn3)
merged_df


#######################################################
# Single-dataset (AA or 3di) runs
# (run in Command-line terminal, not R
#######################################################

# Trim35 = only keep columns with 35% or higher sites 

# AA sequences
cd /GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/
mkdir iqtree_AAsTrim35
cp seqs_AA_famsa3diAligned_trim35.fasta iqtree_AAsTrim35/seqs_AA_famsa3diAligned_trim35.fasta
cp seqs_AA_famsa3diAligned_trim35.fasta iqtree_AAsTrim35/seqs.fasta
cd /GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/iqtree_AAsTrim35/ 
iqtree -s seqs.fasta -mset 3DI,Blosum62,Dayhoff,DCMut,JTT,JTTDCMut,LG,Poisson,Poisson+FQ,Poisson,PMB,WAG,EX2,EX3,EHO,EX_EHO -mfreq FU,F -mrate E,G,R --ufboot 1000 -alrt 1000 -bnni --redo | tee seqs_so1.txt &


cd /GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/
mkdir iqtree_AAsTrim35_allmodels
cp seqs_AA_famsa3diAligned_trim35.fasta iqtree_AAsTrim35_allmodels/seqs_AA_famsa3diAligned_trim35.fasta
cp seqs_AA_famsa3diAligned_trim35.fasta iqtree_AAsTrim35_allmodels/seqs.fasta
cd /GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/iqtree_AAsTrim35_allmodels/ 
iqtree -s seqs.fasta -madd 3DI -mdef 3DI.nexus  --ufboot 1000 -alrt 1000 -bnni --redo | tee seqs_so1.txt &


# 3di sequences
cd /GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/
mkdir iqtree_3disTrim35
cp seqs_3di_famsa3diAligned_trim35.fasta iqtree_3disTrim35/seqs_3di_famsa3diAligned_trim35.fasta
cp seqs_3di_famsa3diAligned_trim35.fasta iqtree_3disTrim35/seqs.fasta
cd /GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/iqtree_3disTrim35/ 
iqtree -s seqs.fasta -mset 3DI,Blosum62,Dayhoff,DCMut,JTT,JTTDCMut,LG,Poisson,Poisson+FQ,Poisson,PMB,WAG,EX2,EX3,EHO,EX_EHO -mfreq FU,F -mrate E,G,R --ufboot 1000 -alrt 1000 -bnni --redo | tee seqs_so1.txt &

cd /GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/
mkdir iqtree_3disTrim35_allmodels
cp seqs_3di_famsa3diAligned_trim35.fasta iqtree_3disTrim35_allmodels/seqs_3di_famsa3diAligned_trim35.fasta
cp seqs_3di_famsa3diAligned_trim35.fasta iqtree_3disTrim35_allmodels/seqs.fasta
cd /GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/iqtree_3disTrim35_allmodels/ 
iqtree -s seqs.fasta -madd 3DI -mdef 3DI.nexus  --ufboot 1000 -alrt 1000 -bnni --redo | tee seqs_so1.txt &


# BOTH sequences
cd /GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/
cp seqs_AA3di_famsa3diAligned.fasta seqs_BOTH_famsa3diAligned.fasta

mkdir iqtree_BOTHsTrim35
cp seqs_BOTH_famsa3diAligned.fasta iqtree_BOTHsTrim35/
cp seqs_BOTH_famsa3diAligned.fasta iqtree_BOTHsTrim35/seqs.fasta
cp BOTHp.raxml iqtree_BOTHsTrim35/
cp 3DI.nexus iqtree_BOTHsTrim35/
cp 3DI iqtree_BOTHsTrim35/
cd /GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/iqtree_BOTHsTrim35/ 
iqtree -s seqs.fasta -spp BOTHp.raxml -m MFP+MERGE -madd 3DI -mdef 3DI.nexus -mset 3DI,Blosum62,Dayhoff,DCMut,JTT,JTTDCMut,LG,Poisson,Poisson+FQ,Poisson,PMB,WAG,EX2,EX3,EHO,EX_EHO -mfreq FU,F -mrate E,G,R --ufboot 1000 -alrt 1000 -bnni --redo | tee seqs_so1.txt &


cd /GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/
mkdir iqtree_BOTHsTrim35_allmodels
cp seqs_BOTH_famsa3diAligned.fasta iqtree_BOTHsTrim35_allmodels/seqs_BOTH_famsa3diAligned.fasta
cp seqs_BOTH_famsa3diAligned.fasta iqtree_BOTHsTrim35_allmodels/seqs.fasta
cp BOTHp.raxml iqtree_BOTHsTrim35_allmodels/
cp 3DI.nexus iqtree_BOTHsTrim35_allmodels/
cp 3DI iqtree_BOTHsTrim35_allmodels/
cd /GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/iqtree_BOTHsTrim35_allmodels/ 
iqtree -s seqs.fasta -spp BOTHp.raxml -m MFP+MERGE -madd 3DI -mdef 3DI.nexus -mfreq FU,F -mrate E,G,R  --ufboot 1000 -alrt 1000 -bnni --redo | tee seqs_so1.txt &





#######################################################
# Open in FigTree and view
# We can also read the IQtree output to tables
#######################################################

source("/GitHub/AA_to_3Di/Rsrc/parsing_iqtree_file_v1.R")

start_txt = "List of models sorted by BIC scores: "
calc_extra = TRUE  # back-calculate k and n
partitioned_TF=FALSE

# AAs
iqtree_fn = "/GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/iqtree_AAsTrim35/seqs.fasta.iqtree"
AICs_AA_somemodels_df = read_iqtree_model_scores(iqtree_fn, start_txt="List of models sorted by BIC scores: ", calc_extra=TRUE, partitioned_TF=FALSE)
conditional_format_table(AICs_AA_somemodels_df)
tail(cft(AICs_AA_somemodels_df))

iqtree_fn = "/GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/iqtree_AAsTrim35_allmodels/seqs.fasta.iqtree"
AICs_AA_allmodels_df = read_iqtree_model_scores(iqtree_fn, start_txt="List of models sorted by BIC scores: ", calc_extra=TRUE, partitioned_TF=FALSE)
conditional_format_table(AICs_AA_allmodels_df)
tail(cft(AICs_AA_allmodels_df))
printall(cft(AICs_AA_allmodels_df))


# 3dis
iqtree_fn = "/GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/iqtree_3disTrim35/seqs.fasta.iqtree"
AICs_3di_somemodels_df = read_iqtree_model_scores(iqtree_fn, start_txt="List of models sorted by BIC scores: ", calc_extra=TRUE, partitioned_TF=FALSE)
conditional_format_table(AICs_3di_somemodels_df)
tail(cft(AICs_3di_somemodels_df))

iqtree_fn = "/GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/iqtree_3disTrim35_allmodels/seqs.fasta.iqtree"
AICs_3di_allmodels_df = read_iqtree_model_scores(iqtree_fn, start_txt="List of models sorted by BIC scores: ", calc_extra=TRUE, partitioned_TF=FALSE)
conditional_format_table(AICs_3di_allmodels_df)
tail(cft(AICs_3di_allmodels_df))
printall(cft(AICs_3di_allmodels_df), 30)


# BOTH
start_txt = "List of models sorted by BIC scores: "
calc_extra = TRUE  # back-calculate k and n
partitioned_TF=TRUE

iqtree_fn = "/GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/iqtree_BOTHsTrim35/BOTHp.raxml.iqtree"
AICs_BOTH_somemodels_df = read_iqtree_model_scores(iqtree_fn, start_txt="List of best-fit models per partition:", calc_extra=TRUE, partitioned_TF=TRUE)
conditional_format_table(AICs_BOTH_somemodels_df)
tail(cft(AICs_BOTH_somemodels_df))

iqtree_fn = "/GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/iqtree_BOTHsTrim35_allmodels/BOTHp.raxml.iqtree"
AICs_BOTH_allmodels_df = read_iqtree_model_scores(iqtree_fn, start_txt="List of best-fit models per partition:", calc_extra=TRUE, partitioned_TF=TRUE)
conditional_format_table(AICs_BOTH_allmodels_df)
tail(cft(AICs_BOTH_allmodels_df))
printall(cft(AICs_BOTH_allmodels_df))








#######################################################
# Plotting trees with R scripts
#######################################################
wd = "/GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/"
setwd(wd)

iqtree_fn = "/GitHub/AA_to_3Di/ex/simple/02_ChimeraX_seqs/03_famsa3di/iqtree_BOTHsTrim35/BOTHp.raxml.treefile"

# Read the tree with APE
tr = ape::read.tree(iqtree_fn)
plot(tr)
add.scale.bar()

# Midpoint root the tree with phytools, ladderize for better display
tr1 = phytools::midpoint.root(tr)
tr1 = ladderize(tr1)
plot(tr1)
add.scale.bar()

# It would be better to have informative labels. 
# Generally these are best kept in e.g. Excel, since alignment & phylogenetics programs 
# often cut the label information. 
xlsfn = "/GitHub/AA_to_3Di/ex/simple/18_ZorABs.xlsx"
xls = openxlsx::read.xlsx(xlsxFile=xlsfn, sheet=1, startRow=2, colNames=TRUE)
head(xls)
dim(xls)


# Functions for editing tipnames etc.
source("/GitHub/AA_to_3Di/Rsrc/seqnames_v1.R")

gidskey = xls$protein.accession
newmatch = zora_names = fixnames(xls$nameA)
tipnames = tr$tip.label
i=2
tr2 = fix_tipnames_gid(gids=gidskey, newmatch=zora_names, tr=tr)

tr
tr2

cbind(tr$tip.label, tr2$tip.label)

plot(tr2)
add.scale.bar()


source("/GitHub/AA_to_3Di/Rsrc/seqnames_v1.R")

pdffn = "ZorA.pdf"
pdf(file=pdffn, width=18, height=12)

tr3 = ape::root(phy=tr2, outgroup="WP_171537395.1_z1_HYP_Acinetobacter_terrestris")

plot(tr3)
add.scale.bar()
title("18 ZorAs, AA+3di data")

plot_nodelabels_on_edges(tr3, cex=0.8)
plot_nodelabels_on_nodes(tr3, cex=0.8)

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)

