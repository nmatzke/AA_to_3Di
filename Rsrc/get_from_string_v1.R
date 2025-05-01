#######################################################
# get_from_string
#
# A bunch of functions/shortcuts to get information from sequence labels
#######################################################


# Source for accession prefixes:
# These prefixes can start sequence IDs
# Source:
# https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly
genbank_prefixes <- function()
	{
	recodes = c("AC_", "NC_",	"NG_", "NT_", "NW_", "NZ_", "NM_", "NR_", "XM_", "XR_", 	"AP_", "NP_", "YP_", "XP_", "WP_", "GCF_", "GCA_")
	return(recodes)
	}


accession_prefixes <- function()
	{
	accession_prefixes = c("AC_", "NC_",	"NG_", "NT_", "NW_", "NZ_", "NM_", "NR_", "XM_", "XR_", 	"AP_", "NP_", "YP_", "XP_", "WP_", "GCF_", "GCA_")
	
	molecule_type = c("Genomic",
"Genomic",
"Genomic",
"Genomic",
"Genomic",
"Genomic",
"mRNA",
"RNA",
"mRNA",
"RNA",
"Protein",
"Protein",
"Protein",
"",
"Protein",
"RefSeq (RS) genome ID",
"GenBank (GB) genome ID")

	comment_txt = c("Complete genomic molecule, usually alternate assembly",
"Complete genomic molecule, usually reference assembly",
"Incomplete genomic region",
"Contig or scaffold, clone-based or WGS a",
"Contig or scaffold, primarily WGS a",
"Complete genomes and unfinished WGS data",
"Protein-coding transcripts (usually curated)",
"Non-protein-coding transcripts",
"Predicted model protein-coding transcript",
"Predicted model non-protein-coding transcript",
"Annotated on AC_ alternate assembly",
"Associated with an NM_ or NC_ accession",
"Annotated on genomic molecules without an instantiated transcript record",
"Predicted model, associated with an XM_ accession",
"Non-redundant across multiple strains and species",
"RefSeq (RS) genome ID",
"GenBank (GB) genome ID")

	footnote = c(
"",
"",
"a Whole Genome Shotgun sequence data.",
"a Whole Genome Shotgun sequence data.",
"b Whole Genome Shotgun sequence data.",
"",
"",
"c Computed.",
"c Computed.",
"",
"",
"c Computed.",
"c Computed.",
"c Computed.",
"",
"")
	
	mat = cbind(accession_prefixes, molecule_type, comment_txt, footnote)
	refseq_prefix_df = as.data.frame(mat, stringsAsFactors=FALSE)
	refseq_prefix_df
	} # END accession_prefixes <- function()


get_info_after_leading_seqid <- function(seqnames, split="_", recodes=genbank_prefixes(), removetxt=c(">"))
	{
	get_info_after_leading_seqid_from_name(seqnames, split=split, recodes=recodes, removetxt=removetxt)
	}


get_info_after_leading_seqid_from_name <- function(seqnames, split="_", recodes=genbank_prefixes(), removetxt=c(">"))
	{
	ex = '
	split=" "
	recodes=genbank_prefixes()
	removetxt=c(">")
	'
	
	seqids = get_leading_seqids_from_name(strings=seqnames, split=split, recodes=recodes, removetxt=removetxt)
	trailing_txt = rep("", times=length(seqnames))
	for (i in 1:length(seqnames))
		{
		word = seqids[i]
		txt_to_remove = paste0(word, split)
		trailing_txt[i] = gsub(pattern=txt_to_remove, replacement="", x=seqnames[i])
		}
	
	return(trailing_txt)
	} # END get_info_after_leading_seqid_from_name


#######################################################
# Get seqids from name (if they are leading the sequence name string)
#######################################################
get_seqids <- function(strings, split="_", recodes=genbank_prefixes(), removetxt=c(">"))
	{
	get_leading_seqids_from_name(strings, split=split, recodes=recodes, removetxt=removetxt)
	}


get_leading_seqids <- function(strings, split="_", recodes=genbank_prefixes(), removetxt=c(">"))
	{
	get_leading_seqids_from_name(strings, split=split, recodes=recodes, removetxt=removetxt)
	}
	
# Defaults are e.g. to get seqids out of a FASTA header
get_leading_seqids_from_name <- function(strings, split="_", recodes=genbank_prefixes(), removetxt=c(">"))
	{
	example_code='
	strings = c("NP_008289.1_Drosophila_melanogaster","NP_050081.1_Dictyostelium_discoideum","NP_059390.1_Paramecium_aurelia","NP_059359.1_Cyanidioschyzon_merolae","NP_037618.1_Phytophthora_infestans","NP_066498.1_Naegleria_gruberi","NP_044754.1_Reclinomonas_americana","NP_066325.1_Malawimonas_jakobiformis","NP_042570.1_Chlamydomonas_reinhardtii","NP_149403.1_Tetrahymena_thermophila")

	
	split=" "
	recodes = genbank_prefixes()
	removetxt=c(">")
	
	seqids = get_leading_seqids_from_name(strings, split="_", recodes=genbank_prefixes(), removetxt=c(">"))
	seqids
	
	tmptxt = ">QQS07318.1 MAG: MotA/TolQ/ExbB proton channel family protein [Fibrobacteres bacterium]"
	firstword(tmptxt, split=" ")
	'
	
	# Remove ">" etc.
	for (i in 1:length(removetxt))
		{
		TF = grepl(pattern=removetxt[i], x=strings)
		strings[TF] = gsub(pattern=removetxt[i], replacement="", x=strings[TF])
		}
	
	# Fix any firstwords that are XP_, YP_, etc.
	# Convert "WP_" etc into parsable
	recodes = recodes
	newcodes = NULL
	IDs1 = NULL
	# Make the newcodes for the recodes
	for (r in 1:length(recodes))
		{
		newcode = gsub(pattern="_", replacement="-", x=recodes[r])
		newcodes = c(newcodes, newcode)
		TF = startsWith(x=strings, prefix=recodes[r])
		length_recode = nchar(recodes[r])
		name_lengths = sapply(X=strings[TF], FUN=nchar)
		strings[TF] = paste0(newcode, mapply(FUN=substr, x=strings[TF], stop=name_lengths, MoreArgs=list(start=length_recode+1)))
		strings[TF]
		}
	
	# NOW, you can split on the delimiter "_"
	seqids = firstwords(strings=strings, split=split)
	seqids
	
	# Return to the recodes
	for (r in 1:length(recodes))
		{
		recode = recodes[r]
		newcode = newcodes[r]
		TF = startsWith(x=seqids, prefix=newcode)
		length_newcode = nchar(newcode)
		name_lengths = sapply(X=seqids[TF], FUN=nchar)
		seqids[TF] = paste0(recode, mapply(FUN=substr, x=seqids[TF], stop=name_lengths, MoreArgs=list(start=length_newcode+1)))
		}	
	
	return(seqids)
	} # END get_leading_seqids_from_name





# Defaults are e.g. to get seqids out of a FASTA header
get_second_seqids_from_name <- function(strings, split="_", recodes=genbank_prefixes(), removetxt=c(">"))
	{
	example_code='
	strings = c("NP_008289.1_Drosophila_melanogaster","NP_050081.1_Dictyostelium_discoideum","NP_059390.1_Paramecium_aurelia","NP_059359.1_Cyanidioschyzon_merolae","NP_037618.1_Phytophthora_infestans","NP_066498.1_Naegleria_gruberi","NP_044754.1_Reclinomonas_americana","NP_066325.1_Malawimonas_jakobiformis","NP_042570.1_Chlamydomonas_reinhardtii","NP_149403.1_Tetrahymena_thermophila")

	
	split=" "
	split="_"
	recodes = genbank_prefixes()
	removetxt=c(">")
	
	seqids = get_leading_seqids_from_name(strings, split="_", recodes=genbank_prefixes(), removetxt=c(">"))
	seqids
	
	tmptxt = ">QQS07318.1 MAG: MotA/TolQ/ExbB proton channel family protein [Fibrobacteres bacterium]"
	firstword(tmptxt, split=" ")
	'
	
	# Remove ">" etc.
	for (i in 1:length(removetxt))
		{
		TF = grepl(pattern=removetxt[i], x=strings)
		strings[TF] = gsub(pattern=removetxt[i], replacement="", x=strings[TF])
		}
	
	# Fix any firstwords that are XP_, YP_, etc.
	# Convert "WP_" etc into parsable
	recodes = recodes
	newcodes = NULL
	IDs1 = NULL
	# Make the newcodes for the recodes
	for (r in 1:length(recodes))
		{
		newcode = gsub(pattern="_", replacement="-", x=recodes[r])
		newcodes = c(newcodes, newcode)
		TF = startsWith(x=strings, prefix=recodes[r])
		length_recode = nchar(recodes[r])
		name_lengths = sapply(X=strings[TF], FUN=nchar)
		strings[TF] = paste0(newcode, mapply(FUN=substr, x=strings[TF], stop=name_lengths, MoreArgs=list(start=length_recode+1)))
		strings[TF]
		}
	
	# NOW, you can split on the delimiter "_"
	seqids = firstwords(strings=strings, split=split)
	seqids
	
	# Remove the first words, repeat on the 2nd words
	seqnames_after_first_seqid = strings
	for (i in 1:length(strings))
		{
		seqnames_after_first_seqid[i] = gsub(pattern=paste0(seqids[i], "_"), replacement="", x=strings[i], ignore.case=TRUE)
		}

	recodes2 = recodes
	newcodes2 = NULL
	IDs2 = NULL
	# Make the newcodes for the recodes
	for (r in 1:length(recodes2))
		{
		newcode2 = gsub(pattern="_", replacement="-", x=recodes2[r])
		newcodes2 = c(newcodes2, newcode2)
		TF = startsWith(x=seqnames_after_first_seqid, prefix=recodes2[r])
		length_recode2 = nchar(recodes2[r])
		name_lengths2 = sapply(X=seqnames_after_first_seqid[TF], FUN=nchar)
		seqnames_after_first_seqid[TF] = paste0(newcode2, mapply(FUN=substr, x=seqnames_after_first_seqid[TF], stop=name_lengths2, MoreArgs=list(start=length_recode2+1)))
		seqnames_after_first_seqid[TF]
		}

	seqids = firstwords(strings=seqnames_after_first_seqid, split=split)
	seqids
	
	
	# Return to the recodes
	for (r in 1:length(recodes))
		{
		recode = recodes[r]
		newcode = newcodes[r]
		TF = startsWith(x=seqids, prefix=newcode)
		length_newcode = nchar(newcode)
		name_lengths = sapply(X=seqids[TF], FUN=nchar)
		seqids[TF] = paste0(recode, mapply(FUN=substr, x=seqids[TF], stop=name_lengths, MoreArgs=list(start=length_newcode+1)))
		}	
	
	return(seqids)
	} # END get_second_seqids_from_name




firstwords <- function(strings, split=" ")
	{
	words = unname(sapply(X=strings, FUN=firstword, split=split))
	return(words)
	}

firstword <- function(string, split=" ")
	{
	example_code='
	tmptxt = ">QQS07318.1 MAG: MotA/TolQ/ExbB proton channel family protein [Fibrobacteres bacterium]"
	secondword(tmptxt, split=" ")
	'

	words = strsplit(gdata::trim(string), split=split)[[1]]
	return(gdata::trim(words[1]))
	}

secondword <- function(string, split=" ")
	{
	example_code='
	tmptxt = ">QQS07318.1 MAG: MotA/TolQ/ExbB proton channel family protein [Fibrobacteres bacterium]"
	secondword(tmptxt, split=" ")
	'

	words = strsplit(gdata::trim(string), split=split)[[1]]
	return(gdata::trim(words[2]))
	}





lastword <- function(string, split="/")
	{
	example_code='
	tmptxt = ">QQS07318.1 MAG: MotA/TolQ/ExbB proton channel family protein [Fibrobacteres bacterium]"
	lastword(tmptxt, split=" ")
	
	# Run on multiple inputs
	tmptxts = c(tmptxt, tmptxt)
	results = sapply(X=tmptxts, FUN=lastword, split=" ")
	results
	unname(results)
	'
	words = strsplit(gdata::trim(string), split=split)[[1]]
	return(gdata::trim(words[length(words)]))
	}



get_firstword <- function(OTU, split="_")
	{
	if (is.na(OTU) == TRUE)
		{
		return(NA)
		}
	firstword = strsplit(x=OTU, split=split)[[1]][1]
	return(firstword)
	}

get_lastword <- function(OTU, split="_")
	{
	if (is.na(OTU) == TRUE)
		{
		return(NA)
		}
	words = strsplit(x=OTU, split=split)[[1]]
	lastword = words[length(words)]
	return(lastword)
	}

get_lastwords <- function(OTUs, split="_")
	{
	lastwords = sapply(X=OTUs, FUN=get_lastword, split=split)
	lastwords = unname(lastwords)
	return(lastwords)
	}


get_2nd_to_lastword <- function(OTU, split="_")
	{
	if (is.na(OTU) == TRUE)
		{
		return(NA)
		}
	words = strsplit(x=OTU, split=split)[[1]]
	
	if (length(words) < 2)
		{
		warning("WARNING: get_2nd_to_lastword() reports that input 'OTU' splits to less than 2 words. Returning NA.")
		return(NA)
		}
	
	lastword = words[length(words)-1]
	return(lastword)
	}



get_secondword <- function(OTU, split="_")
	{
	if (is.na(OTU) == TRUE)
		{
		return(NA)
		}
	secondword = strsplit(x=OTU, split=split)[[1]][2]
	return(secondword)
	}
get_thirdword <- function(OTU, split="_")
	{
	if (is.na(OTU) == TRUE)
		{
		return(NA)
		}
	thirdword = strsplit(x=OTU, split=split)[[1]][3]
	return(thirdword)
	}
get_fourthword <- function(OTU, split="_")
	{
	if (is.na(OTU) == TRUE)
		{
		return(NA)
		}
	fourthword = strsplit(x=OTU, split=split)[[1]][4]
	return(fourthword)
	}

get_firstwords <- function(OTUs, split="_")
	{
	firstwords = sapply(X=OTUs, FUN=get_firstword, split=split)
	firstwords = unname(firstwords)
	return(firstwords)
	}

get_secondwords <- function(OTUs, split="_")
	{
	secondwords = sapply(X=OTUs, FUN=get_secondword, split=split)
	secondwords = unname(secondwords)
	return(secondwords)
	}


get_fourthwords <- function(OTUs, split="_")
	{
	fourthwords = sapply(X=OTUs, FUN=get_fourthword, split=split)
	fourthwords = unname(fourthwords)
	return(fourthwords)
	}



get_all_but_suffix <- function(tmptxt, split="\\.")
	{
	all_but_suffix(tmptxt, split=split)
	}

get_all_but_suffixes <- function(tmptxts, split="\\.")
	{
	all_but_suffixes(tmptxts, split=split)
	}



all_but_suffix <- function(tmptxt, split="\\.")
	{
	ex='
	tmptxt = c("29447-30982 atp1 ATP synthase F1 subunit alpha","31153-31583 rpL11","31606-32112 rpL10","32134-32505 rpS12","lcl|ORF31 32515-32994 rpS7","33239-33397 rpL32","33408-34031 ccmA heme ABC exporter ATP-binding protein","34034-34726 ccmB heme exporter protein","34743-35501 ccmC ABC transporter subunit cytochrome c biogenesis C","35523-37409 ccmF heme lyase CcmF-NrfE family subunit","complement(37398-37532) rpL34","complement(37593-39170) cox1")
	subtxt = all_but_prefix(tmptxt=tmptxt, split=" ")
	cbind(subtxt, tmptxt)
	
	tmptxt = "complement(48035-48451) atp8 ATP synthase F0 subunit 8"
	split = " "
	all_but_prefix(tmptxt=tmptxt, split=" ")
	all_but_prefixes(tmptxt=tmptxt, split=" ")

	tmptxt = "atp8 ATP synthase F0 subunit 8 complement(48035-48451)"
	split = " "
	all_but_suffix(tmptxt=tmptxt, split=" ")
	all_but_suffixes(tmptxt=tmptxt, split=" ")
	' # END example

	if (length(tmptxt) != 1)
		{
		txt = paste0("STOP ERROR in all_but_suffix(). The length of tmptxt= must be 1. Instead, it was ", length(tmptxt), ". Fix and re-run, or perhaps use all_but_suffixes().")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		}

	words = strsplit(tmptxt, split=split)[[1]]
	
	# This fails with e.g. () in tiplabel
	#word_to_remove = paste0(split, words[length(words)])
	#subtxt = gsub(pattern=word_to_remove, replacement="", x=tmptxt)

	subtxt = paste(words[1:(length(words)-1)], sep="", collapse=split)

	return(subtxt)
	} # END all_but_suffix <- function(tmptxt, split="\\.")


all_but_suffixes <- function(tmptxts, split="\\.")
	{
	ex='
	tmptxt = c("29447-30982 atp1 ATP synthase F1 subunit alpha","31153-31583 rpL11","31606-32112 rpL10","32134-32505 rpS12","lcl|ORF31 32515-32994 rpS7","33239-33397 rpL32","33408-34031 ccmA heme ABC exporter ATP-binding protein","34034-34726 ccmB heme exporter protein","34743-35501 ccmC ABC transporter subunit cytochrome c biogenesis C","35523-37409 ccmF heme lyase CcmF-NrfE family subunit","complement(37398-37532) rpL34","complement(37593-39170) cox1")
	subtxt = all_but_prefix(tmptxt=tmptxt, split=" ")
	cbind(subtxt, tmptxt)
	
	tmptxt = "complement(48035-48451) atp8 ATP synthase F0 subunit 8"
	split = " "
	all_but_prefix(tmptxt=tmptxt, split=" ")
	all_but_prefixes(tmptxt=tmptxt, split=" ")

	tmptxt = "atp8 ATP synthase F0 subunit 8 complement(48035-48451)"
	split = " "
	all_but_suffix(tmptxt=tmptxt, split=" ")
	all_but_suffixes(tmptxt=tmptxt, split=" ")
	' # END example

	subtxts = unname(sapply(X=tmptxts, FUN=all_but_suffix, split=split))
	return(subtxts)
	}





get_all_but_prefix <- function(tmptxt, split=" ")
	{
	all_but_prefix(tmptxt, split=split)
	}


get_all_but_prefixes <- function(tmptxts, split=" ")
	{
	all_but_prefixes(tmptxts, split=split)
	}



# Remove everything before the first delimiter
all_but_prefix <- function(tmptxt, split=" ")
	{
	ex='
	tmptxt = c("29447-30982 atp1 ATP synthase F1 subunit alpha","31153-31583 rpL11","31606-32112 rpL10","32134-32505 rpS12","lcl|ORF31 32515-32994 rpS7","33239-33397 rpL32","33408-34031 ccmA heme ABC exporter ATP-binding protein","34034-34726 ccmB heme exporter protein","34743-35501 ccmC ABC transporter subunit cytochrome c biogenesis C","35523-37409 ccmF heme lyase CcmF-NrfE family subunit","complement(37398-37532) rpL34","complement(37593-39170) cox1")
	subtxt = all_but_prefix(tmptxt=tmptxt, split=" ")
	cbind(subtxt, tmptxt)
	
	tmptxt = "complement(48035-48451) atp8 ATP synthase F0 subunit 8"
	split = " "
	all_but_prefix(tmptxt=tmptxt, split=" ")
	all_but_prefixes(tmptxt=tmptxt, split=" ")

	tmptxt = "atp8 ATP synthase F0 subunit 8 complement(48035-48451)"
	split = " "
	all_but_suffix(tmptxt=tmptxt, split=" ")
	all_but_suffixes(tmptxt=tmptxt, split=" ")
	' # END example
	
	if (length(tmptxt) != 1)
		{
		txt = paste0("STOP ERROR in all_but_prefix(). The length of tmptxt= must be 1. Instead, it was ", length(tmptxt), ". Fix and re-run, or perhaps use all_but_prefixes().")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		}
	words = strsplit(tmptxt, split=split)[[1]]
	
	# This fails on "complement(48035-48451) atp8 ATP synthase F0 subunit 8"
	#word_to_remove = paste0(words[1], split) 
	#subtxt = gsub(pattern=word_to_remove, replacement="", x=tmptxt)
	
	subtxt = paste(words[2:length(words)], sep="", collapse=split)
	
	return(subtxt)
	} # END all_but_prefix <- function(tmptxt, split=" ")




all_but_prefixes <- function(tmptxts, split=" ")
	{
	ex='
	tmptxt = c("29447-30982 atp1 ATP synthase F1 subunit alpha","31153-31583 rpL11","31606-32112 rpL10","32134-32505 rpS12","lcl|ORF31 32515-32994 rpS7","33239-33397 rpL32","33408-34031 ccmA heme ABC exporter ATP-binding protein","34034-34726 ccmB heme exporter protein","34743-35501 ccmC ABC transporter subunit cytochrome c biogenesis C","35523-37409 ccmF heme lyase CcmF-NrfE family subunit","complement(37398-37532) rpL34","complement(37593-39170) cox1")
	subtxt = all_but_prefix(tmptxt=tmptxt, split=" ")
	cbind(subtxt, tmptxt)
	
	tmptxt = "complement(48035-48451) atp8 ATP synthase F0 subunit 8"
	split = " "
	all_but_prefix(tmptxt=tmptxt, split=" ")
	all_but_prefixes(tmptxt=tmptxt, split=" ")

	tmptxt = "atp8 ATP synthase F0 subunit 8 complement(48035-48451)"
	split = " "
	all_but_suffix(tmptxt=tmptxt, split=" ")
	all_but_suffixes(tmptxt=tmptxt, split=" ")
	' # END example

	subtxts = unname(sapply(X=tmptxts, FUN=all_but_prefix, split=split))
	return(subtxts)
	} # END all_but_prefixes <- function(tmptxts, split=" ")


# Split off seqID as well as [Genus species strain]
get_label_without_seqID_or_species <- function(seqnames, split=" ", recodes=genbank_prefixes(), removetxt=c(">"), replace_spaces=FALSE, prflag=FALSE)
	{
	ex='
	split=" "
	recodes=genbank_prefixes()
	removetxt=c(">")
	replace_spaces=FALSE
	prflag=FALSE
	'
	
	all_but_prefixes = get_info_after_leading_seqid_from_name(seqnames, split=split, recodes=recodes, removetxt=removetxt )
	# Remove the spnames after brackets
	all_but_prefixes2 = all_but_prefixes
	for (i in 1:length(all_but_prefixes2))
		{
		words = strsplit(all_but_prefixes2[i], split="\\[")[[1]]
		all_but_prefixes2[i] = gdata::trim(words[1])
		}
	all_but_prefixes2
	return(all_but_prefixes2)
	}


get_spnames <- function(list_of_strings, replace_spaces=TRUE, prflag=FALSE)
	{
	extract_last_brackets(list_of_strings=list_of_strings, replace_spaces=replace_spaces, prflag=prflag)
	}


get_last_brackets <- function(list_of_strings, replace_spaces=TRUE, prflag=FALSE)
	{
	extract_last_brackets(list_of_strings=list_of_strings, replace_spaces=replace_spaces, prflag=prflag)
	}

# Handy function
extract_last_brackets <- function(list_of_strings, replace_spaces=TRUE, prflag=FALSE)
	{
	example_code='
	replace_spaces=TRUE
	prflag=FALSE
	
	tmptxt = ">QQS07318.1 MAG: MotA/TolQ/ExbB proton channel family protein [Fibrobacteres bacterium]"
	list_of_strings = tmptxt; replace_spaces=TRUE
	species_names = extract_last_brackets(list_of_strings=list_of_strings, replace_spaces=replace_spaces)
	species_names
	
	# NA check
	source("/GitHub/str2phy/Rsrc/seqnames_v1.R")
	replace_spaces=TRUE
	prflag=FALSE
	extract_last_brackets(list_of_strings=NA, replace_spaces=replace_spaces)
	
	'
	
	# NA error trap
	if ((length(list_of_strings) == 1) && (is.na(list_of_strings) == TRUE) )
		{
		return(NA)
		}

	species_names = rep("", length(list_of_strings))
	
	if (prflag == TRUE)
		{
		txt = paste0("\nextract_last_brackets() is processing ", length(list_of_strings), " strings. String #")
		cat(txt)
		}
	
	for (i in 1:length(list_of_strings))
		{
		if (prflag == TRUE)
			{
			cat(i, ",", sep="")
			} # END if
		tmptxt = list_of_strings[i]
		tmpstrs = unlist(regmatches(tmptxt, gregexpr("\\[.+?\\]", tmptxt)))
		tmpstr = tmpstrs[length(tmpstrs)] # take the last bracketed text, if more than 1
		
		# Remove "[", "]"
		tmpstr = gsub(pattern="\\[", replacement="", x=tmpstr)
		tmpstr = gsub(pattern="\\]", replacement="", x=tmpstr)
		
		# Replace spaces with "_"
		if (replace_spaces == TRUE)
			{
			tmpstr = gsub(pattern=" ", replacement="_", x=tmpstr)
			}
		species_names[i] = tmpstr
		}
	return(species_names)
	}

