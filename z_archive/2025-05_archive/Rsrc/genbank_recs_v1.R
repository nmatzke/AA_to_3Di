
#######################################################
# Get db_xref (database cross-references) IDs for each seqid
#######################################################

get_db_xref_for_seqids <- function(seqids, recs)
	{
	ex = '
	
	# Some MotAs
	seqids = c("AAC07752","AAC73831","AAC74960","AAC76042","AAG03587","AAG04082","AAG04358","AAG04849","AAG06371","AAG08339")
	seqids = c("AAC07752","AAC73831")
	
	gb_files <- reutils::efetch(uid=seqids, db="protein", rettype = "gp", retmode = "text")
	tmpfn = "~/gb_files.gp"
	write(content(gb_files, "text"), file=tmpfn)
	moref(tmpfn)
	
	# parse (locally saved)
	recs <- biofiles::gbRecord(tmpfn, progress=TRUE)

	# parse the efetch object into a gbRecord instance.
	# (from online)
	recs <- biofiles::gbRecord(gb_files, progress=TRUE)
	dbxrefs_table = get_db_xref_for_seqids(seqids, recs)
	'
	
	CDS_db_xref_list = NULL
	protein_db_xref_list = NULL
	region_db_xref_list = NULL
	source_db_xref_list = NULL
	sig_peptide_db_xref_list = NULL
	
	db_xref_list = NULL
	nums = 1:10
	cat("\nGetting db_xref names for ", length(recs), " gbRecords in gbRecordList: ", sep="")
	for (i in 1:length(recs))
		{
		# Print every 10th i
		if ((i %% 10) == 0)
			{
			cat(i, ",", sep="")
			}
		
		r = recs[[i]]
		
		# traps for missing features
		if ("CDS" %in% names(biofiles::featureTable(r)))
			{
			CDS_db_xref_list[[i]] = biofiles::dbxref(r["CDS"])
			} else {
			CDS_db_xref_list[[i]] = NA
			}
		if ("Protein" %in% names(biofiles::featureTable(r)))
			{
			protein_db_xref_list[[i]] = biofiles::dbxref(r["Protein"])
			} else {
			protein_db_xref_list[[i]] = NA
			}
		if ("Region" %in% names(biofiles::featureTable(r)))
			{
			region_db_xref_list[[i]] = biofiles::dbxref(r["Region"])
			} else {
			region_db_xref_list[[i]] = NA
			}
		if ("source" %in% names(biofiles::featureTable(r)))
			{
			source_db_xref_list[[i]] = biofiles::dbxref(r["source"])
			} else {
			source_db_xref_list[[i]] = NA
			}
		if ("sig_peptide" %in% names(biofiles::featureTable(r)))
			{
			sig_peptide_db_xref_list[[i]] = biofiles::dbxref(r["sig_peptide"])
			} else {
			sig_peptide_db_xref_list[[i]] = NA
			}
		}
	cat("...done!\n")
	
	CDS_db_xref_names = unique(unlist(sapply(X=CDS_db_xref_list, FUN=names)))
	protein_db_xref_names = unique(unlist(sapply(X=protein_db_xref_list, FUN=names)))
	region_db_xref_names = unique(unlist(sapply(X=region_db_xref_list, FUN=names)))
	source_db_xref_names = unique(unlist(sapply(X=source_db_xref_list, FUN=names)))
	sig_peptide_db_xref_names = unique(unlist(sapply(X=sig_peptide_db_xref_list, FUN=names)))
	
	CDS_db_xref_names
	protein_db_xref_names
	region_db_xref_names
	source_db_xref_names
	sig_peptide_db_xref_names
	
	dbxref_names = unique(c(CDS_db_xref_names, protein_db_xref_names, region_db_xref_names, source_db_xref_names, sig_peptide_db_xref_names))
	#dbxref_names = dbxref_names[dbxref_names != "X1"]
	#catc(dbxref_names)
	
	# Build up a table of dbxrefs for every record
	dbxrefs_table = matrix(NA, nrow=length(recs), ncol=length(dbxref_names)+2)
	dbxrefs_table = as.data.frame(dbxrefs_table, stringsAsFactors=FALSE)
	names(dbxrefs_table) = c("i", "seqid", dbxref_names)
	
	dbxrefs_table$i = 1:length(recs)
	dbxrefs_table$seqid = seqids
	
	# Go through the lists and accumulate db_xrefs
	cat("\nGetting db_xref entries for ", length(recs), " gbRecords in gbRecordList: ", sep="")
	for (i in 1:length(recs))
		{
		# Print every 10th i
		if ((i %% 10) == 0)
			{
			cat(i, ",", sep="")
			}
		
		if (((length(unlist(CDS_db_xref_list[[i]]))==1) && (is.na(CDS_db_xref_list[[i]]) == TRUE)) == FALSE)
			{
			tmpnames = names(CDS_db_xref_list[[i]])
			if (nrow(CDS_db_xref_list[[i]]) > 1) {
				CDS_db_xref_list[[i]] = apply(X=CDS_db_xref_list[[i]], FUN=paste0, MARGIN=2, collapse=",") }
			dbxrefs_table[i, tmpnames] = CDS_db_xref_list[[i]]
			}
		if (((length(unlist(protein_db_xref_list[[i]]))==1) && (is.na(protein_db_xref_list[[i]]) == TRUE)) == FALSE)
			{
			tmpnames = names(protein_db_xref_list[[i]])
			if (nrow(protein_db_xref_list[[i]]) > 1) {
				protein_db_xref_list[[i]] = apply(X=protein_db_xref_list[[i]], FUN=paste0, MARGIN=2, collapse=",") }
			dbxrefs_table[i, tmpnames] = protein_db_xref_list[[i]]
			}
		if (((length(unlist(region_db_xref_list[[i]]))==1) && (is.na(region_db_xref_list[[i]]) == TRUE)) == FALSE)
			{
			tmpnames = names(region_db_xref_list[[i]])
			if (nrow(region_db_xref_list[[i]]) > 1) {
				region_db_xref_list[[i]] = apply(X=region_db_xref_list[[i]], FUN=paste0, MARGIN=2, collapse=",") }
			dbxrefs_table[i, tmpnames] = region_db_xref_list[[i]]
			}
		if (((length(unlist(source_db_xref_list[[i]]))==1) && (is.na(source_db_xref_list[[i]]) == TRUE)) == FALSE)
			{
			tmpnames = names(source_db_xref_list[[i]])
			if (nrow(source_db_xref_list[[i]]) > 1) {
				source_db_xref_list[[i]] = apply(X=source_db_xref_list[[i]], FUN=paste0, MARGIN=2, collapse=",") }
			dbxrefs_table[i, tmpnames] = source_db_xref_list[[i]]
			}
		if (((length(unlist(sig_peptide_db_xref_list[[i]]))==1) && (is.na(sig_peptide_db_xref_list[[i]]) == TRUE)) == FALSE)
			{
			tmpnames = names(sig_peptide_db_xref_list[[i]])
			if (nrow(sig_peptide_db_xref_list[[i]]) > 1) {
				sig_peptide_db_xref_list[[i]] = apply(X=sig_peptide_db_xref_list[[i]], FUN=paste0, MARGIN=2, collapse=",") }
			dbxrefs_table[i, tmpnames] = sig_peptide_db_xref_list[[i]]
			}
		} # END for (i in 1:length(recs))
	cat("...done!\n")
	
	# Remove the blank X1 field
	num = match("X1", names(dbxrefs_table))
	dbxrefs_table = dbxrefs_table[,-num]
	
	#head(dbxrefs_table)
	return(dbxrefs_table)
	} # END get_db_xref_for_seqids <- function(seqids, recs)
	
