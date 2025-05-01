#######################################################
# Get protein information:
# Convert seqids to information
#######################################################
require(ginmappeR)
require(queryup)

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


get_uniprot_data_on_seqids <- function(seqids, runslow=TRUE, base_fn="seqids", version="v1")
	{
	ex='
	runslow=TRUE, 
	base_fn = "379_AQBs" # includes directory
	version = "v1"
	
	bigdf_outdf1 = get_uniprot_data_on_seqids(seqids, runslow=runslow, base_fn=base_fn, version=version)

	runslow = FALSE
	base_fn = "379_AQBs" # includes directory
	version = "v1"
	bigdf_outdf2 = get_uniprot_data_on_seqids(seqids, runslow=runslow, base_fn=base_fn, version=version)
	'
	
	
	out_table_fn = paste0(base_fn, "_wUniProtdata_", version, ".Rdata")
	
	
	# Get identical protein seqids for each seqid, using getNCBIIdenticalProteins
	identical_protIDs_fn = paste0(base_fn, "_identical_proteins_", version, ".Rdata")
	identical_protein_seqids_list = get_IDs_identical_proteins(seqids, runslow=runslow, identical_protIDs_fn =identical_protIDs_fn)
	# Counts
	counts_identical_proteins = sapply(FUN=length, X=identical_protein_seqids_list)
	rev(sort(counts_identical_proteins))
	sort(counts_identical_proteins)
	
	# Take a list of identical proteins, return just the ones that match genbank_prefixes()
	
	# Extract IDs with genbank prefixes
	recodes = genbank_prefixes()
	identical_protein_seqids_list_wGenBank_codes = sapply(X=identical_protein_seqids_list, FUN=reduce_identical_proteins_to_matching_codes, codes=recodes)
	seqids_wGenBank_codes_counts = sapply(X=identical_protein_seqids_list_wGenBank_codes, FUN=length)
	
	seqids_wGenBank_codes_txt = sapply(X=identical_protein_seqids_list_wGenBank_codes, FUN=paste0, split="", collapse=",")
	
	# Use queryup to search for UniProt for each entry
	library(queryup)
	query_fields$field
	#uniprot_xref_list = get_uniprot_xrefs(seqids)
	
	# fields you can return:
	printall(return_fields)
	
	query <- list("xref" = seqids)
	
	runslow = runslow
	#seqids_to_uniprot_accession_entries_fn = "seqids_to_uniprot_accession_entries_df_v1.Rdata"
	seqids_to_uniprot_accession_entries_fn = paste0(base_fn, "_seqids_to_uniprot_accession_entries_df_", version, ".Rdata")
	
	if (runslow)
		{
		seqids_to_uniprot_accession_entries_df = get_uniprot_accession_by_seqid(seqids, collapse_multiple_hits_to_1row=TRUE)
		save(seqids_to_uniprot_accession_entries_df, file=seqids_to_uniprot_accession_entries_fn)
		} else {
		# Loads to: seqids_to_uniprot_accession_entries_df
		load(file=seqids_to_uniprot_accession_entries_fn)
		} # END if (runslow)
	
	# Get lots of fields, but not linked to seqid
	bigdf = queryup::query_uniprot(query, columns=c("accession", "id", "gene_names", "protein_name", "organism_name", "organism_id", "reviewed", "xref_refseq", "xref_cdd", "xref_pfam", "organelle", "length", "mass"), show_progress=TRUE)
	# Do these separately, they crash when linked to other fields
	bigdf2_pdbsum = queryup::query_uniprot(query, columns=c("xref_pdbsum"), show_progress=TRUE)
	bigdf2_pdb = queryup::query_uniprot(query, columns=c("xref_pdb"), show_progress=TRUE)
	bigdf2_alphafold = queryup::query_uniprot(query, columns=c("xref_alphafolddb"), show_progress=TRUE)
	
	dim(seqids_to_uniprot_accession_entries_df)
	bigdf_outdf = link_query_uniprot_table_to_accession_table(seqids_to_uniprot_accession_entries_df, bigdf)
	dim(bigdf_outdf)
	
	bigdf_outdf = link_query_uniprot_table_to_accession_table(bigdf_outdf, bigdf2_pdbsum)
	dim(bigdf_outdf)
	bigdf_outdf = link_query_uniprot_table_to_accession_table(bigdf_outdf, bigdf2_pdb)
	dim(bigdf_outdf)
	bigdf_outdf = link_query_uniprot_table_to_accession_table(bigdf_outdf, bigdf2_alphafold)
	dim(bigdf_outdf)
	
	
	bigdf_outdf = cbind(bigdf_outdf, counts_identical_proteins, seqids_wGenBank_codes_txt, seqids_wGenBank_codes_counts)
	bigdf_outdf
	
	return(bigdf_outdf)
	}

get_IDs_identical_proteins <- function(seqids, runslow=TRUE, identical_protIDs_fn="protIDs_identical_to_seqIDs_v1.Rdata")
	{
	ex='
	seqids = c("AAC07752",
	"AAC73831",
	"AAC74960",
	"AAC76042",
	"AAG03587")
	
	runslow = TRUE
	identical_protIDs_fn = "protIDs_identical_to_MotA_seqIDs_v1.Rdata"
	' # END ex
	
	if (runslow == TRUE)
		{
		identical_protein_seqids_list = list()
		cat("\n")
		cat("\nBecause runslow=TRUE, get_IDs_identical_proteins() is using ginmappeR::getNCBIIdenticalProteins() to find seqids of identical proteins for ", length(seqids), " seqids...", sep="")
		cat("\n")
		
		for (i in 1:length(seqids))
			{
			
			# Print every 10th
			if ((i %% 10) == 0)
				{
				cat(i, "...", sep="")
				} # END if ((i %% 10) == 0
			result = try(getNCBIIdenticalProteins(seqids[i]))
			if ("try-error" %in% class(result))
				{
				identical_protein_seqids_list[[i]] = NA
				} else {
				identical_protein_seqids_list[[i]] = unique(result[[1]])
				}
			} # END for (i in 1:length(seqids))
		# Saves to identical_protein_seqids_list
		save(identical_protein_seqids_list, file=identical_protIDs_fn)
		} else {
		# runslow is FALSE
		# Loads to identical_protein_seqids_list
		cat("\nBecause runslow=FALSE, get_IDs_identical_proteins() is attempting to load results from identical_protIDs_fn='", identical_protIDs_fn, "'...", sep="")
		# Should load to "identical_protein_seqids_list"
		loaded_object_name = load(file=identical_protIDs_fn)
		cmdtxt = paste0("identical_protein_seqids_list = ", loaded_object_name)
		eval(expr=parse(text=cmdtxt))
		cat("...done. Results in: identical_protein_seqids_list .\n")
		} # END if (runslow)
	
	return(identical_protein_seqids_list)
	} # END get_IDs_identical_proteins <- function(seqids, runslow=TRUE, identical_protIDs_fn="protIDs_identical_to_MotA_seqIDs_v1.Rdata")



reduce_identical_proteins_to_matching_codes <- function(identical_seqids, codes=genbank_prefixes())
	{
	ex='
	codes=genbank_prefixes()
	'
	
	matchnums = match_grepl(patterns=codes, x=identical_seqids, return_counts=FALSE)
	if (all(is.na(matchnums)) == TRUE)
		{
		identical_seqids_matched = NULL
		} else {
		matchnums = matchnums[!is.na(matchnums)]
		identical_seqids_matched = identical_seqids[matchnums]
		}
	
	return(identical_seqids_matched)
	} # END reduce_identical_proteins_to_matching_codes <- function(identical_seqids, codes=genbank_prefixes())


# Return a reasonable list, when you have a huge list of matching identical sequence IDs
reduce_identical_proteins_to_reasonable_list <- function(identical_seqids, codes=genbank_prefixes(), num_other_IDs=5)
	{
	# Get the IDs that match GenBank
	matchnums = match_grepl(patterns=codes, x=identical_seqids, return_counts=FALSE)
	if (all(is.na(matchnums)) == TRUE)
		{
		identical_seqids_matched = NULL
		} else {
		matchnums = matchnums[!is.na(matchnums)]
		identical_seqids_matched = identical_seqids[matchnums]
		}

	# Get the shortest IDs
	numchars_by_ID = nchar(identical_seqids)
	num_other_IDs_to_include = min(c(num_other_IDs,length(identical_seqids)))
	shortest_IDs = identical_seqids[order(numchars_by_ID)][1:num_other_IDs_to_include]
	
	reasonable_list = unique(c(shortest_IDs, identical_seqids_matched))
	reasonable_list = reasonable_list[order(nchar(reasonable_list))]
	return(reasonable_list)
	} # END reduce_identical_proteins_to_reasonable_list(identical_seqids, codes=genbank_prefixes(), num_other_IDs=5)



# Quite slow, but the only way to get the seqid -> uniprot accession match;
# later queries can get the rest to match
get_uniprot_accession_by_seqid <- function(seqids, collapse_multiple_hits_to_1row=TRUE)
	{
	ex='
	# 
	seqids = c("AAC07752","AAC73831","AAC74960","AAC76042","AAG03587")
	collapse_multiple_hits_to_1row = TRUE
	seqids_to_uniprot_accession_entries_df = get_uniprot_accession_by_seqid(seqids)
	'
	
	cat("\nget_uniprot_accession_by_seqid() is retrieving UniProt xref records for ", length(seqids), " seqids...", sep="")
	uniprot_xref_list = list()
	for (i in 1:length(seqids))
		{
		# Print every 10th
		if ((i %% 1) == 0)
			{
			cat(i, ",", sep="")
			} # END if ((i %% 10) == 0
		
		seqid = seqids[i]
		query <- list("xref" = seqid)
		tmpdf = queryup::query_uniprot(query, columns=c("accession"), show_progress = FALSE)
		tmpdf2 = as.data.frame(cbind(rep(i, times=nrow(tmpdf)), rep(seqid, times=nrow(tmpdf)), tmpdf$Entry), stringsAsFactors=FALSE)
		names(tmpdf2) = c("i", "seqid", "Entry")

		# Fix when there is more than 1 accession hit
		if (collapse_multiple_hits_to_1row == TRUE)
			{
			if (nrow(tmpdf2) > 1)
				{
				new_accessions = tmpdf2$Entry
				numchars = nchar(new_accessions)
				new_accession_txt = paste0(new_accessions[order(numchars)], collapse=",", sep="")
				new_accession_txt
				tmpdf2 = tmpdf2[1,]
				tmpdf2$Entry = new_accession_txt
				} # END if (nrow(tmpdf2) > 1)
			} # END if (collapse_multiple_hits_to_1row == TRUE)
		uniprot_xref_list[[i]] = tmpdf2
		} # END for (i in 1:length(seqids))
	cat("...done. Processing into a big table:")
	
	uniprot_xref_mat = NULL
	for (i in 1:length(seqids))
		{
		seqid = seqids[i]
		if (nrow(uniprot_xref_list[[i]]) == 0)
			{
			accession = NA
			tmpline = c(i, seqid, accession)
			uniprot_xref_mat = rbind(uniprot_xref_mat, tmpline)
			} else {
			i_col = rep(i, times=length(uniprot_xref_list[[i]]$Entry))
			seqid_col = rep(seqid, times=length(uniprot_xref_list[[i]]$Entry))
			tmpline = cbind(i_col, seqid_col, uniprot_xref_list[[i]]$Entry)
			uniprot_xref_mat = rbind(uniprot_xref_mat, tmpline)
			}
		}
	
	seqids_to_uniprot_accession_entries_df = as.data.frame(uniprot_xref_mat, stringsAsFactors=FALSE)
	names(seqids_to_uniprot_accession_entries_df) = c("i", "seqid", "accession")
	row.names(seqids_to_uniprot_accession_entries_df) = NULL
	head(seqids_to_uniprot_accession_entries_df)
	dim(seqids_to_uniprot_accession_entries_df)
	
	return(seqids_to_uniprot_accession_entries_df)
	} # END get_uniprot_accession_by_seqid(seqids, collapse_multiple_hits_to_1row=TRUE)


link_query_uniprot_table_to_accession_table <- function(seqids_to_uniprot_accession_entries_df, bigdf)
	{
	matchnums_bigdf_to_seqids = match(x=bigdf$Entry, table=seqids_to_uniprot_accession_entries_df$accession)
	counts_of_matches = table(matchnums_bigdf_to_seqids)
	if (max(counts_of_matches, na.rm=TRUE) > 1)
		{
		txt = "STOP ERROR in link_bigdf_via_accession_Entry(): multiple entries in bigdf$Entry matched the same seqids_to_uniprot_accession_entries_df$accession."
		stop(txt)
		}
	
	# Find commas in accession; take the shortest accession
	tmp_accessions = seqids_to_uniprot_accession_entries_df$accession
	tmp_TF = grepl(pattern=",", x=tmp_accessions)
	if (sum(tmp_TF) > 0)
		{
		for (i in 1:length(tmp_TF))
			{
			if (tmp_TF[i] == TRUE)
				{
				words = strsplit(tmp_accessions[i], split=",")[[1]]
				numchars = nchar(words)
				words = words[order(numchars)]
				word = words[1]
				tmp_accessions[i] = word
				} # END if (tmp_TF[i] == TRUE)
			} # END for (i in 1:length(tmp_TF))
		} # END if (sum(tmp_TF) > 0)
	
	matchnums_seqids_to_bigdf = match(x=tmp_accessions, table=bigdf$Entry)
	counts_of_matches = table(matchnums_seqids_to_bigdf)
	if (max(counts_of_matches, na.rm=TRUE) > 1)
		{
		txt = "Warning in link_bigdf_via_accession_Entry(): multiple entries in seqids_to_uniprot_accession_entries_df$accession matched the same bigdf$Entry. This suggests you have multiple seqids that ginmappeR::getNCBIIdenticalProteins() matched to the same UniProt Entry/accession. Not a problem, just be aware of it."
		warning(txt)
		}
		
	bigdf_outdf = cbind(seqids_to_uniprot_accession_entries_df, matchnums_seqids_to_bigdf, bigdf[matchnums_seqids_to_bigdf,])
	return(bigdf_outdf)
	} # END link_query_uniprot_table_to_accession_table <- function(seqids_to_uniprot_accession_entries_df, bigdf)
