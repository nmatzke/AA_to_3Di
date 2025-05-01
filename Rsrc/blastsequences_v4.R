setup='

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")

'


# 
# source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")
# 

# ALL, All Fields, All terms from all searchable fields
# UID, UID, Unique number assigned to each sequence
# FILT, Filter, Limits the records
# WORD, Text Word, Free text associated with record
# TITL, Title, Words in definition line
# KYWD, Keyword, Nonstandardized terms provided by submitter
# AUTH, Author, Author(s) of publication
# JOUR, Journal, Journal abbreviation of publication
# VOL, Volume, Volume number of publication
# ISS, Issue, Issue number of publication
# PAGE, Page Number, Page number(s) of publication
# ORGN, Organism, Scientific and common names of organism, and all higher levels of taxonomy
# ACCN, Accession, Accession number of sequence
# PACC, Primary Accession, Does not include retired secondary accessions
# GENE, Gene Name, Name of gene associated with sequence
# PROT, Protein Name, Name of protein associated with sequence
# ECNO, EC/RN Number, EC number for enzyme or CAS registry number
# PDAT, Publication Date, Date sequence added to GenBank
# MDAT, Modification Date, Date of last update
# SUBS, Substance Name, CAS chemical name or MEDLINE Substance Name
# PROP, Properties, Classification by source qualifiers and molecule type
# SQID, SeqID String, String identifier for sequence
# GPRJ, Genome Project, Genome Project
# SLEN, Sequence Length, Length of sequence
# FKEY, Feature key, Feature annotated on sequence
# RTYP, Replicon type, Replicon type
# RNAM, Replicon name, Replicon name
# ORGL, Organelle, Organelle

# E.g. mitochondrion[ORGL]
# E.g. srcdb refseq[Properties]
#      &srcdb%20refseq%5BProperties%5D

# http://www.ncbi.nlm.nih.gov/blast/html/blastcgihelp.html
# https://web.archive.org/web/20151230120228/http://www.ncbi.nlm.nih.gov/blast/html/blastcgihelp.html
#
# https://blast.auckland.ac.nz/blast/fasta.shtml
# https://blast.auckland.ac.nz/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp
# https://blast.auckland.ac.nz/doc/blast-topics/
#

# BLAST® Command Line Applications User Manual
# https://www.ncbi.nlm.nih.gov/books/NBK279690/
# https://www.ncbi.nlm.nih.gov/books/NBK279690/pdf/Bookshelf_NBK279690.pdf

# LIMITS
# (txid9608 [ORGN] OR txid9688 [ORGN]) NOT (txid9615 [ORGN] OR txid9606 [ORGN] OR (srcdb refseq model[prop] AND biomol rna[prop]) OR environmental samples[organism] OR metagenomes[orgn] OR txid32644[orgn])

# CODES
# ASCII Encoding Reference
# http://www.w3schools.com/tags/ref_urlencode.asp

# SPECIAL CHARACTERS IN URLs
# space  = %20
# [ = %5B
# ] = %5D
# \ = %5C
# ( = %28
# ) = %29

# OTHERS
# ! = %21
# ! = %21
# " = %22
# # = %23
# $ = %24
# % = %25
# & = %26



#######################################################
# Example search to get an RID:
#
# https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?QUERY=MQLEKMLTEVNIERSKFNTYNVLLQYIKAVTQKYNLNKVSSKYGKVLFIKDGVVKVSGLSQIKIGEKVEFVGKNLYGMALNLEATSVGIVIFGEDTAIYEGVIVKRCEQNFAIKVDKTMLGRVVDVLGQPIDGLGELKDTKTTRVMSVERKAPGIVTRKSVHESMLTGVKIVDALLPIGRGQRELIIGDRQTGKSAIAVDAILNQQVNKDIVCIYVAVGQKKSTVRRLVEMLNTKGALEYTIVVVSTASDAAPLQFLAPYTGCTIGEYFRDEGKHALIVYDDLSKHAVAYRQMSLLLRRPPGREAYPGDVFYIHSRLLERAAKLNEKYGCGSLTAFPIVETQAGDVSAYIPTNIISITDGQIFLEKELFNKGIRPAVNVGLSVSRVGSAAQSAVMKKLAGALKLELAQYRELARFEQFSSNADAVTTQILKKGKLTIELLKQVNNNPMSIGLEALMIYAMGTSYFQNLDLSLVRSEETKLLNYINSIASFKLYAACVDAVKAFNPKDTVFAEMCQSYVK&DATABASE=nr&HITLIST_SIZE=25&FILTER=L&EXPECT=10&PROGRAM=blastp&ENTREZ_QUERY=(txid176299[orgn] OR txid335992[orgn] OR txid366602[orgn] OR txid156889[orgn] OR txid426355[orgn] OR txid349102[orgn] OR txid269796[orgn] OR txid1105096[orgn] OR txid392499[orgn] OR txid163164[orgn] )&CMD=Put
#
# 
#######################################################


# Retrieve sequences from Genbank based on ID
# Retrieve sequences from NCBI based on ID
# Download sequences from NCBI based on ID
# _retrieve_
# _download_
# downloaded_seqs = rentrez::entrez_fetch(db="protein", id=ids, rettype="fasta")
# breaks down at large numbers of sequences

# ONLY GOOD FOR SMALL NUMBERS OF SEQUENCES;
# Try: https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_EFetch_ 
# ...for manual download

download_seqs_from_genbank <- function(ids, db="protein", rettype="fasta", maxnum=50, type="AA", return_what="sequences", outfn="tmpseqs.txt")
	{
	ex = '
	ids = c("AAK87078.2", "ACB24041.1", "Q2RU27.1")
	db="protein"
	rettype="fasta"
	maxnum=50
	type="AA"
	return_what = "sequences"
	#return_what = "string"
	#return_what = "text_filename")
	outfn="tmpseqs.txt"
	'
	cat("\n")
	downloaded_seqs_txt = NULL
	
	if (length(ids) > maxnum)
		{
		numchunks = ceiling(length(ids)/maxnum)
		} else {
		numchunks = 1
		}
	
	for (i in 1:numchunks)
		{
		txt = paste0("Downloading ", length(ids), " sequences from GenBank '", db,
"' database, chunk ", i, "/", numchunks, "...")
		
		cat("\n")
		cat(txt)
		
		startnum = ((i-1) * maxnum) + 1
		endnum = (i * maxnum)
		tmpids = ids[startnum:endnum]
		
		tmpseqs = rentrez::entrez_fetch(db=db, id=tmpids, rettype=rettype)
		downloaded_seqs_txt = c(downloaded_seqs_txt, tmpseqs)
		}
	cat("\n...done7!\n")
	
	# Paste the downloaded_seqs together into 1 text string:
	downloaded_seqs_txt = paste0(downloaded_seqs_txt, collapse="\n\n")

	if ((return_what == "string") || (return_what == "text_filename"))
		{
		# Remove double-returns:
		downloaded_seqs_txt = gsub(pattern="\n\n", replacement="\n", x=downloaded_seqs_txt)
		downloaded_seqs_txt = gsub(pattern="\n\n", replacement="\n", x=downloaded_seqs_txt)
		downloaded_seqs_txt = gsub(pattern="\n\n", replacement="\n", x=downloaded_seqs_txt)
		downloaded_seqs_txt = gsub(pattern="\n\n", replacement="\n", x=downloaded_seqs_txt)
		return(downloaded_seqs_txt)
		}

	# Write the string to a text file
	#tmpfn = "tmpseqs.txt"
	write(x=downloaded_seqs_txt, file=outfn, append=FALSE, sep="\n")
	
	if (return_what == "text_filename")
		{
		txt = paste0("Results of download_seqs_from_genbank() saved to the text file '", outfn, "'.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		return(outfn)
		}


	if (return_what == "sequences")
		{
		downloaded_seqs = ape::read.FASTA(file=outfn, type=type)
		return(downloaded_seqs)
		}

	if (return_what == "fasta_filename")
		{
		downloaded_seqs = ape::read.FASTA(file=outfn, type=type)
		ape::write.FASTA(x=downloaded_seqs, file=outfn, header=NULL, append=FALSE)
		txt = paste0("Results of download_seqs_from_genbank() saved to the fasta file '", outfn, "'.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		return(outfn)
		}

	
	# Shouldn't get here; but the default
	txt = paste0("WARNING: download_seqs_from_genbank() could not match your 'return_what=", return_what, "' input to any programmed options for ; returning results as text string.")
	cat("\n\n")
	cat(txt)
	cat("\n\n")
	warning(txt)
	return(downloaded_seqs_txt)
	}



convert_txt_to_NCBI <- function(txt)
	{
	defaults='
	txt = "Canis[orgn] AND (srcdb RefSeqGene[Properties] OR srcdb RefSeq[Properties])"
	convert_txt_to_NCBI(txt)
	'
	
	# Spaces
	txt2 = txt
	txt2 = gsub(pattern=" ", replace="%20", x=txt2)
	# [
	txt2 = gsub(pattern="\\[", replace="%5B", x=txt2)
	# ]
	txt2 = gsub(pattern="\\]", replace="%5D", x=txt2)
	# [
	txt2 = gsub(pattern="\\(", replace="%28", x=txt2)
	# ]
	txt2 = gsub(pattern="\\)", replace="%29", x=txt2)
	
	txt2
	return(txt2)
	}


convert_seqs_to_seqCollapse <- function(x)
	{
	# collapse seq into string
	seqCollapse <- paste(toupper(as.character(x)), collapse = "")
	return(seqCollapse)
	}



current_time_txt <- function()
	{
	currtime = as.character(Sys.time())
	currtime
	
	words = strsplit(currtime, split="\\.")
	currtime = words[[1]][1]
	
	currtime = gsub(pattern=" ", replacement="_", x=currtime)
	currtime = gsub(pattern=":", replacement=".", x=currtime)
	currtime
	return(currtime)
	}


current_time_txt_wd <- function(prefix=NULL)
	{
	currtime = as.character(Sys.time())
	currtime
	
	words = strsplit(currtime, split="\\.")
	currtime = words[[1]][1]
	
	currtime = gsub(pattern=" ", replacement="_", x=currtime)
	currtime = gsub(pattern=":", replacement=".", x=currtime)
	currtime
	
	# Add the working directory
	tmpwd = getwd()
	
	if (is.null(prefix))
		{
		currtime = slashslash(paste0(addslash(tmpwd), currtime))
		} else {
		currtime = slashslash(paste0(addslash(tmpwd), prefix, "_", currtime))
		}
	currtime
	return(currtime)
	}


# sapply version of check_blast_txid
check_blast_txids <- function(taxid_txts)
	{
	TFs = sapply(X=taxid_txts, FUN=check_blast_txid)
	return(TFs)
	}

get_RID_from_url <- function(search_url)
	{
	ex='
	search_url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?RID=SP6VDV24013&FORMAT_TYPE=XML&CMD=Get"
	rid = get_RID_from_url(search_url)
	rid
	
	search_url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
	rid = get_RID_from_url(search_url)
	rid
	'
	words = strsplit(search_url, split="?RID=")[[1]]
	str2 = words[2]
	words = strsplit(str2, split="&")[[1]]
	words
	rid = words[1]
	return(rid)
	}


#######################################################
# Check taxid
#######################################################
# For BLAST via URL, the [ORGN] field has to have this format:
# txid1234
#
# NOT:
# 1234
# taxid1234
# NA
check_blast_txid <- function(taxid_txt)
	{
	ex='
	library(openxlsx)
	library(ape)
	library(phytools)
	library(BioGeoBEARS)
	library(varhandle) # for check.numeric
	library(rentrez) # Fetching full records: entrez_fetch()
	
	#sourceall("/GitHub/bioinfRhints/Rsrc/")
	sourceall("/GitHub/str2phy/Rsrc/")
	source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")
	
	taxid_txt = NA
	varhandle::check.numeric(taxid_txt, na.rm=FALSE)
	varhandle::check.numeric(taxid_txt, na.rm=TRUE)
	
	taxid_txt = NA
	check_blast_txid(taxid_txt)
	taxid_txt = "NA"
	check_blast_txid(taxid_txt)
	taxid_txt = "na"
	check_blast_txid(taxid_txt)
	taxid_txt = NULL
	check_blast_txid(taxid_txt)
	taxid_txt = "1234"
	check_blast_txid(taxid_txt)
	taxid_txt = "taxid1234"
	check_blast_txid(taxid_txt)
	taxid_txt = "txid1234"
	check_blast_txid(taxid_txt)
	
	organism = c("txid1234","txid1235","txid1236")
	check_blast_txid(organism)

	organism = c("txid1234","txid1235","human")
	check_blast_txid(organism)

	
	taxid_txts = c("txid1234","txid1235","txid1236")
	check_blast_txids_comma_delimited(taxid_txts)
	taxid_txts = "taxid1234,txid1235,txid1236"
	check_blast_txids_comma_delimited(taxid_txts)
	taxid_txts = "txid1234,txid1235,NA"
	check_blast_txids_comma_delimited(taxid_txts)
	taxid_txts = "txid1234,txid1235,1236"
	check_blast_txids_comma_delimited(taxid_txts)
	taxid_txts = "txid1234,txid1235,txid1236a"
	check_blast_txids_comma_delimited(taxid_txts)
	taxid_txts = "txid1234,txid1235,txid1236"
	check_blast_txids_comma_delimited(taxid_txts)	
	' # END ex
	
	# Allow NULL inputs through
	if (is.null(taxid_txt) == TRUE)
		{
		return(TRUE)
		}
	
	if (length(taxid_txt) > 1)
		{
		# run sapply verison:
		TFs = check_blast_txids(taxid_txt)
		return(TFs)
		}
	
	# Check for NA
	naTF1 = is.na(taxid_txt)
	naTF2 = toupper(taxid_txt) == "NA"
	if ((naTF1 == TRUE) || (naTF2 == TRUE))
		{
		errortxt = paste0("STOP ERROR #1 in check_blast_txid(). Taxonomy ID txt='", taxid_txt, "' converts to 'NA', which will break a BLAST search URL. Fix your input taxids and re-try.")
		cat("\n")
		cat(errortxt)
		cat("\n")
		stop(errortxt)
		} # END if ((naTF1 == TRUE) || (naTF2 == TRUE))
	
	# Check for "txid" (the only valid text)
	txidTF = grepl(pattern="txid", x=taxid_txt)
	if (txidTF == TRUE)
		{
		taxid_txt2 = gsub(pattern="txid", replacement="", x=taxid_txt)
		numericTF = varhandle::check.numeric(taxid_txt2)
		if (numericTF == TRUE)
			{
			return(TRUE)
			} else {
			errortxt = paste0("STOP ERROR #2 in check_blast_txid(). Taxonomy ID txt='", taxid_txt, "' correctly has 'txid', but the rest should convert to pure numbers, but does not. Fix your input taxids and retry.")
			cat("\n")
			cat(errortxt)
			cat("\n")
			stop(errortxt)
			} # END if (numericTF == TRUE)
		} else {
		errortxt = paste0("STOP ERROR #3 in check_blast_txid(). Taxonomy ID txt='", taxid_txt, "' must have 'txid' (not 'taxid') before the number, but does not. Fix your input taxids and retry.")
		cat("\n")
		cat(errortxt)
		cat("\n")
		stop(errortxt)
		} # END if (txidTF == TRUE)
	
	# Shouldn't get here
	stop("STOP ERROR #4 in check_blast_txid(). Shouldn't get here.")
	} # END check_blast_txid <- function(taxid_txt)

# Check a comma-delimited list in a string
check_blast_txids_comma_delimited <- function(taxid_txts)
	{
	ex='
	library(openxlsx)
	library(ape)
	library(phytools)
	library(BioGeoBEARS)
	library(varhandle) # for check.numeric
	library(rentrez) # Fetching full records: entrez_fetch()
	
	#sourceall("/GitHub/bioinfRhints/Rsrc/")
	sourceall("/GitHub/str2phy/Rsrc/")
	source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")
	
	taxid_txt = NA
	varhandle::check.numeric(taxid_txt, na.rm=FALSE)
	varhandle::check.numeric(taxid_txt, na.rm=TRUE)
	
	taxid_txt = NA
	check_blast_txid(taxid_txt)
	taxid_txt = "NA"
	check_blast_txid(taxid_txt)
	taxid_txt = "na"
	check_blast_txid(taxid_txt)
	taxid_txt = NULL
	check_blast_txid(taxid_txt)
	taxid_txt = "1234"
	check_blast_txid(taxid_txt)
	taxid_txt = "taxid1234"
	check_blast_txid(taxid_txt)
	taxid_txt = "txid1234"
	check_blast_txid(taxid_txt)
	
	taxid_txts = c("txid1234","txid1235","txid1236")
	check_blast_txids_comma_delimited(taxid_txts)
	taxid_txts = "taxid1234,txid1235,txid1236"
	check_blast_txids_comma_delimited(taxid_txts)
	taxid_txts = "txid1234,txid1235,NA"
	check_blast_txids_comma_delimited(taxid_txts)
	taxid_txts = "txid1234,txid1235,1236"
	check_blast_txids_comma_delimited(taxid_txts)
	taxid_txts = "txid1234,txid1235,txid1236a"
	check_blast_txids_comma_delimited(taxid_txts)
	taxid_txts = "txid1234,txid1235,txid1236"
	check_blast_txids_comma_delimited(taxid_txts)	
	' # END ex



	# Allow NULL inputs through
	if (is.null(taxid_txts) == TRUE)
		{
		return(TRUE)
		}

	if (length(taxid_txts) > 1)
		{
		txt = paste0("STOP ERROR #1 in check_blast_txids_comma_delimited(). Input taxid_txts should have a length of 1, i.e. be a single character variable. Revise inputs and re-run.")
		catn()
		cat(txt)
		catn()
		stop(txt)
		}

	if (is.character(taxid_txts) == FALSE)
		{
		txt = paste0("STOP ERROR #2 in check_blast_txids_comma_delimited(). Input taxid_txts should be a string of type 'character', but is.character(taxid_txts) returns FALSE. Revise inputs and re-run.")
		catn()
		cat(txt)
		catn()
		stop(txt)
		}
	
	# Otherwise, continue with parsing on comma
	words = strsplit(x=taxid_txts, split=",")[[1]]
	for (i in 1:length(words))
		{
		check_blast_txid(taxid_txt=words[i])
		}
	# If you get to the end, return TRUE
	return(TRUE)
	} # END check_blast_txids_comma_delimited <- function(taxid_txts)

#######################################################
# Construct the URLs to submit a job to online BLAST
#
# These could be pasted into a command-line, or 
# submitted via another R function.
# 
# Database options:
# UniProtKB/Swiss-Prot (swissprot)
# Non-redundant protein sequences (nr)
# Reference proteins (refseq_protein)
# Reference Select proteins (refseq_select)
# Protein Data Bank proteins (pdb)
# 
# See also: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp
#######################################################
construct_blasturl <- function (searchseq, database="nr", hitListSize="10", 
                        filter="L", expect="10", program="blastp",
                        organism=NULL, not_organism=NULL, addl_url=NULL, 
                        txt=NULL, baseUrl="http://blast.ncbi.nlm.nih.gov/blast/Blast.cgi", assume_txid=TRUE, printout=TRUE)
#baseUrl="http://www.ncbi.nlm.nih.gov/blast/Blast.cgi")
	{
	defaults='
	#wd = "/drives/Dropbox/_njm/"
	#setwd(wd)
	fasta_fn = "/drives/Dropbox/_njm/Wayne_etal_1997_SysBio_canid_mtDNA_appendix_wolf_cytB.fasta"
	mySeq = read.fasta(fasta_fn)

	searchseq = convert_seqs_to_seqCollapse(x=mySeq[[1]])
	database="nr"
	hitListSize="10"
	filter="L"
	expect="10"
	program="blastp"
	attempts=10
	organism=NULL
	not_organism = NULL
	addl_url=NULL
	txt=NULL
	justRID_url=FALSE
	attempts=10
	retry_wait=5
	baseUrl = "https://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
	#baseUrl = "uoa"
	assume_txid = TRUE
	
	
	results = blastSeqKK2(searchseq=searchseq, database="nr", hitListSize="10", 
                        filter="L", expect="10", program="blastn",
                        organism=NULL, not_organism=NULL, addl_url=NULL, 
                        txt=txt, justRID_url=FALSE, attempts=10, retry_wait=5, baseUrl="default")
	' # END defaults
	
	# Check on input
	if (length(searchseq) != 1)
		{
		txt = paste0('STOP ERROR in construct_blasturl(). Input "searchseq" must have length 1, but instead it has length=', length(searchseq), '. Probably you have input the sequence as a list of characters instead of a single string. Try e.g. \n searchseq = paste0(as.character(searchseq), collapse="")\n....')
		
		catn()
		cat(txt)
		catn()
		stop(txt)
		} # END if (length(searchseq) != 1)
	
	# Default
	# baseUrl="http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
	# Other, non-blast, baseUrls
	# http://eutils.ncbi.nlm.nih.gov/blast/Blast.cgi
	# If you are on a university campus that auto-replaces the URL, you will have to customize this:
	if (baseUrl == "uoa")
		{
		baseUrl = "https://blast-ncbi-nlm-nih-gov.ezproxy.auckland.ac.nz/Blast.cgi"
		}
	if (baseUrl == "default")
		{
		#baseUrl = "https://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
		baseUrl = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
		}
	
	
	if (is.null(filter) == TRUE)
		{
		query = paste("QUERY=", as.character(searchseq), "&DATABASE=", database, 
			 "&HITLIST_SIZE=", hitListSize, "&EXPECT=", 
			 expect, "&PROGRAM=", program, sep = "")
		} else {
		# Include filter
		query = paste("QUERY=", as.character(searchseq), "&DATABASE=", database, 
			 "&HITLIST_SIZE=", hitListSize, "&FILTER=", filter, "&EXPECT=", 
			 expect, "&PROGRAM=", program, sep = "")
		}
	

	#######################################################
	# Search by including or excluding one or more
	# organisms, identified by name or txid (NOT taxid)
	#######################################################
	# organism, and not_organism, are vectors of strings
	
	# Error the user-specificed vector of taxids, if assume_txid==TRUE
	# This checks that every string starts with txid, etc.
	# Throws a STOP ERROR if not
	if (assume_txid == TRUE)
		{
		TFs = check_blast_txid(organism)
		TFs = check_blast_txid(not_organism)
		}

	# Put in organism limitation
	entrez_query = FALSE
	q2 = ""
	if (!is.null(organism) || !is.null(not_organism))
		{
		entrez_query = TRUE
		
		# Add e.g. "human[orgn]" to query
		# http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node25.html
		if (length(c(organism, not_organism)) == 1)
			{
			if (length(organism) == 1)
				{
				q2 = paste0(q2, organism, "[orgn]&")
				}
			if (length(not_organism) == 1)
				{
				q2 = paste0(q2, "NOT ", not_organism, "[orgn]&")
				}
			} else {
			
			if ( (length(organism) > 1) || (length(not_organism) > 1) )
				{
				if (length(organism) > 0)
					{
					orgn_txt = rep("[orgn]", times=length(organism))
					organisms_txt1 = paste(unlist(organism), orgn_txt, sep="")
					organisms_txt1
					organisms_txt2 = paste(organisms_txt1, collapse=" OR ", sep="")
					organisms_txt2
					} else {
					organisms_txt2 = NULL
					}
				if (length(not_organism) > 0)
					{
					NOT_txt = rep("NOT ", times=length(not_organism))
					orgn_txt = rep("[orgn]", times=length(not_organism))
					not_organisms_txt1 = paste(NOT_txt, unlist(not_organism), orgn_txt, sep="")
					not_organisms_txt1
					not_organisms_txt2 = paste(not_organisms_txt1, collapse=" ", sep="")
					not_organisms_txt2
					} else {
					not_organisms_txt2 = NULL
					}
				}
			
			txt3 = paste0("(", organisms_txt2, " ", not_organisms_txt2, ")")
			q2 = paste0(q2, txt3, "&")
			}
		} # END if (!is.null(organism) || !is.null(not_organism))



	# Put in additional (manual) limitations into the URL string
	if (!is.null(addl_url))
		{
		entrez_query = TRUE
		
		# Add e.g. "human[orgn]" to query
		# http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node25.html
		q2 = paste0(q2, addl_url, "&")
		} # END if (!is.null(addl_url))

	# Put in additional (text, in ENTREZ format) limitations on the search
	if (!is.null(txt))
		{
		entrez_query = TRUE
		
		# Convert to URL format...
		#addl_txt = convert_txt_to_NCBI(txt)
		# Add e.g. "human[orgn]" to query
		# http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node25.html
		q2 = paste0(q2, txt, "&")
		} # END if (!is.null(addl_url))
	
	# If you've added anything anywhere, add it to "query"
	if (entrez_query == TRUE)
		{
		q2 = paste0("&ENTREZ_QUERY=", q2)
		query = paste0(query, q2)
		} # END if (entrez_query = TRUE)
	q2

	# Convert the URL to HTML format
	query = gsub(pattern="&&", replacement="&", x=query)
	query_orig = query
	query = convert_txt_to_NCBI(query_orig)		
	RID_url = sprintf("%s?%s&CMD=Put", baseUrl, query)
	RID_url_orig = sprintf("%s?%s&CMD=Put", baseUrl, query_orig)
	RID_url = gsub(pattern="&&", replacement="&", x=RID_url)
	RID_url_orig = gsub(pattern="&&", replacement="&", x=RID_url_orig)
	
	if (printout == TRUE)
		{
		cat("\n\nconstruct_blasturl() says: The search URL, 'results$RID_url' has been assembled:\n")
		cat(RID_url_orig)
		cat("\n")
		}

	results = NULL
	results$RID_url_orig = RID_url_orig
	results$RID_url = RID_url
	return(results)
	} # END construct_blasturl 
	

# Given a results object from construct_blasturl(), or a RID_url,
# query BLAST to get:
# * an RID (Request ID) for a particular search
# * an rtoe (expected waiting time for the BLAST result)
#
# Note: This command does *not* run the search, just 
# sets it up on the NCBI BLAST website
# 
# RIDs are temporary and are deleted after a few hours
# 
# (They may be saveable, if you have an NCBI BLAST account)
# 
get_blast_RID_waiting_time <- function(results=NULL, RID_url=NULL, baseUrl="default")
	{
	defaults='
	RID_url = NULL
	baseUrl = "default"
	'
	
	
	require(XML)
	
	# Default
	# baseUrl="http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
	# Other, non-blast, baseUrls
	# http://eutils.ncbi.nlm.nih.gov/blast/Blast.cgi
	# If you are on a university campus that auto-replaces the URL, you will have to customize this:
	
	if (is.null(baseUrl))
		{
		baseUrl = "default"
		}

	if (baseUrl == "default")
		{
		#baseUrl = "https://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
		baseUrl = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
		}
	if (baseUrl == "uoa")
		{
		baseUrl = "https://blast-ncbi-nlm-nih-gov.ezproxy.auckland.ac.nz/Blast.cgi"
		}
	
	
	if (is.null(results) && is.null(RID_url))
		{
		txt = paste0("STOP ERROR in get_blast_RID_waiting_time(): both results and RID_url cannot both be NULL.")
		stop(txt)
		}
	
	if (is.null(results))
		{
		results = NULL
		results$RID_url = RID_url
		} else {
		RID_url = results$RID_url
		}
	
	
	# Get the search RID
	cat("\nDownloading HTML from NCBI 'Format Request' page...")
	RIDrequest_fn = paste0(current_time_txt_wd(), "_request_RID_before_BLAST.html")
	tmp = download.file(RID_url, destfile=RIDrequest_fn)
	cat("done1.\n")
	#moref(destfile)

	cat("Extracting RID (Request ID)...rid=")
	post = htmlTreeParse(RIDrequest_fn, useInternalNodes = TRUE)
	x = post[["string(//comment()[contains(., \"QBlastInfoBegin\")])"]]
	rid = sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
	cat(rid)
	cat("\n")
	
	# Put the search RID into a url for searching/downloading
	# RTOE = Estimated waiting time before the BLAST results will be ready
	rtoe = as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1", x))
	
	cat("Estimated waiting time for search is rtoe=")
	cat(as.character(rtoe))
	cat(" seconds.\n")
	
	# Search commands
	search_url = sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl, rid)

	terminal_cmd = paste0("open -n '", search_url, "'")


	results$RID_url = RID_url
	results$search_url = search_url
	results$terminal_cmd = terminal_cmd
	results$RIDrequest_fn = RIDrequest_fn
	results$rid = rid
	results$rtoe = rtoe
	
	return(results)
	} # END get_blast_RID_waiting_time

start_BLAST_websearch <- function(results=NULL, terminal_cmd=NULL, search_url=NULL)
	{
	if (is.null(results) && is.null(terminal_cmd) && is.null(search_url))
		{
		txt = paste0("STOP ERROR in start_BLAST_websearch(): all three inputs (results and terminal_cmd and search_url) cannot be NULL at once.")
		stop(txt)
		}
	
	if (is.null(results) && is.null(search_url))
		{
		results = NULL
		results$terminal_cmd = terminal_cmd
		rtoe = "(unknown)"
		words = strsplit(terminal_cmd, split="RID=")[[1]]
		str2 = words[2]
		words = strsplit(str2, split="&")[[1]]
		words
		rid = words[1]
		}
		
	if (is.null(terminal_cmd) && is.null(search_url))
		{
		terminal_cmd = results$terminal_cmd
		rtoe = results$rtoe
		rid = results$rid
		}
	
	# Option for just searching with the search_url (no RID provided)
	if (is.null(search_url) == FALSE)
		{
		terminal_cmd = paste0("open -n '", search_url, "'")
		rtoe = NA
		rid = NA
		}
	
	cat("Starting NCBI BLAST web search with RID=", rid, " via a Terminal command-line command. This search is projected to take approximately ", rtoe, " seconds.\n\nThis is the terminal command (will open a window in the default browser):\n", sep="")
	cat(terminal_cmd)
	cat("\n...")
	
	# terminal_cmd = paste0("open -n '", search_url, "'")
	system(terminal_cmd)
	
	cat("\n\nOnce the search is completed -- in Chrome, you will see something like:\n'This XML file does not appear to have any style information associated with it. The document tree is shown below.'\nOnce finished, you can use download_BLAST_xml() to download the completed XML file for parsing.\n")
		
	return(NULL)
	} # END start_BLAST_websearch



download_BLAST_xml <- function(results=NULL, search_url=NULL, attempts=1, retry_wait=1, blastXML_fn="temp.xml", already_saved=FALSE)
	{
	ex='
	attempts=2
	retry_wait=2
	search_url = "https://blast-ncbi-nlm-nih-gov.ezproxy.auckland.ac.nz/blast/Blast.cgi?RID=21131NB6016&FORMAT_TYPE=XML&CMD=Get"
	search_url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?RID=211ZD50W013&FORMAT_TYPE=XML&CMD=Get"
	blastXML_fn="temp.xml"
	already_saved=FALSE
	'
	
	if (is.null(results) && is.null(search_url))
		{
		txt = paste0("STOP ERROR in start_BLAST_websearch(): both results and terminal_cmd cannot both be NULL.")
		stop(txt)
		}
	
	if (is.null(results))
		{
		results = NULL
		results$search_url = search_url
		
		# Extract the RID from the search_url
		#tmpRID = gsub(pattern="https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi\\?RID=", replacement="", x=search_url)
		#tmpRID
		#tmpRID = gsub(pattern="&FORMAT_TYPE=XML&CMD=Get", replacement="", x=tmpRID)
		#tmpRID
		
		results$blastXML_fn = blastXML_fn
		}
		
	if (is.null(search_url))
		{
		search_url = results$search_url
		}

	cat("\nAttempting to download the completed BLAST XML file with .tryParseResult()'s download.file()...\n")
	
	# Wait for the waiting time, plus a bit
	#if (justRID_url == FALSE)
	#	{
	#	Sys.sleep(rtoe + runif(n=1, min=retry_wait, max=2*retry_wait))
	#	}
	
	# Try to parse
	tryParse_results = .tryParseResult(search_url, attempts, retry_wait=retry_wait, blastXML_fn=results$blastXML_fn, already_saved=already_saved)
	
	cat("done2.\n")
	
	# tmp = download.file(results$search_url, destfile=results$blastXML_fn)
	
	#moref(results$blastXML_fn)
	
	return(tryParse_results)
	} # END download_BLAST_xml






#######################################################
# Table approach:
# 1. Build table from a FASTA file
# 2. Construct BLAST urls for each
# 3. Getting RIDs for each search
# 4. Launching each search
# 5. Checking if they are done & downloading
# 6. Parsing, once they are done.
#######################################################

# The "search df" is a search dataframe/table containing the results
# of the BLAST search attempt on each sequence in the original FASTA file
read_search_df <- function(fasta_blast_table_fn)
	{
	search_df = read.table(file=fasta_blast_table_fn, header=TRUE, sep="\t", quote="", row.names=NULL, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors = FALSE)
	return(search_df)
	}

write_search_df <- function(search_df, fasta_blast_table_fn)
	{
	write.table(x=search_df, file=fasta_blast_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	return(fasta_blast_table_fn)
	}

subset_search_df <- function(search_df)
	{
	ex='
	fasta_blast_table_fn = "/GitHub/str2phy/Rsrc/gradual_blast_ex/7mitocogs_blast_table_v1.txt"
	search_df = read_search_df(fasta_blast_table_fn)
	subcols_search_df = subset_search_df(search_df)
	subcols_search_df
	
	search_df$RIDrequest_fn
	
	matdf = check_num_still_running(fasta_blast_table_fn)
	matdf
	'
	
	headings = c("i","cTF","rTF","lTF","sTF","pTF","eTF","ncbiWarn","pTF_numhits","seqid","taxname","len","rid", "blastXML_fn", "blastres_fns")
	newnames = c("i","cTF","rTF","lTF","sTF","pTF","eTF","ncbiWarn","nhit","seqid","taxname","len","rid", "blastXML_fn", "blastres_fns")
	subcols_search_df = search_df[,headings]
	names(subcols_search_df) = newnames
	return(subcols_search_df)
	}


#######################################################
# 1. Build table from a FASTA file
#######################################################
build_blastsearch_table_from_fasta <- function(fasta_fn, type="AA", fasta_blast_table_fn=NULL, overwrite=FALSE)
	{
	cmds='
	wd = "/GitHub/str2phy/ex/MitoCOGs/_03_alphafolds/"
	setwd(wd)
	fasta_fn = "133aas_wo_UniProt.fasta"
	type = "AA"
	search_df = build_blastsearch_table_from_fasta(fasta_fn, type=type)
	head(search_df)
	class(search_df)
	dim(search_df)
	'
	
	# Come up with filename
	if (is.null(fasta_blast_table_fn) == TRUE)
		{
		prefix = all_but_suffix(fasta_fn)
		fasta_blast_table_fn = paste0(prefix, "_blast_table_v1.txt")
		}
	
	# Load sequences
	seqs = ape::read.FASTA(fasta_fn, type=type)
	seqs_chars = as.character(seqs)
	seq_length = sapply(X=seqs_chars, FUN=length)
	seqs_strings = sapply(X=seqs_chars, FUN=paste0, collapse="")
	length(seqs)
	seqs_strings

	# Extract the seq ids, even though it is delimited by "_" (special code to identify NP_, WP_, etc.)
	seqnames = names(seqs)
	seqids = get_leading_seqids_from_name(strings=seqnames, split="_", recodes=genbank_prefixes(), removetxt=c(">"))
	seqids

	# Process seqnames to genus, species
	taxnames = get_info_after_leading_seqid_from_name(seqnames, split="_", recodes=genbank_prefixes(), removetxt=c(">"))

	# Make a big table
	ivals = 1:length(seqs)
	blank = rep("", times=length(seqs))
	blank0 = rep(0, times=length(seqs))
	blankF = rep(FALSE, times=length(seqs))
	tmptable = cbind(ivals, blank, blank, blank, blank, blank, blank, blankF, blank0, seqids, taxnames, seq_length, seqnames, blank, seqs_strings)
	search_df = as.data.frame(tmptable, stringsAsFactors=FALSE)
	names(search_df) = c("i", "cTF", "rTF", "lTF", "sTF", "pTF", "eTF", "ncbiWarn", "pTF_numhits", "seqid", "taxname", "len", "seqnames", "blastres_fns", "seqstring")
	row.names(search_df) = NULL
	# cTF = constrTF  = has the search URL been constructed?
	# rTF = ridTF     = has the RID & search time been retrived?
	# lTF = launchTF  = has the search been launched?
	# sTF = saveTF    = have the results been saved?
	# pTF = parseTF   = have the results been successfully parsed?
	# eTF = errorTF   = the program thinks the whole run was done & saved, but there was an error nontheless
	# pTF_numhits     = the number of hits when a blast result is processed
	# rtime           = time the RID results came back and indicated search time (usually a wild underestimate)
	# ltime           = time the BLAST search was launched
	# stime           = time the BLAST search was saved
	# ptime           = time the BLAST search was parsed

	# BLAST Request ID (RID) and expected times
	blastres_names = c("RID_url_orig","RID_url","search_url","terminal_cmd","RIDrequest_fn","blastXML_fn","rid","rtoe")

	# BLAST results names
	blast_fields = c("rtime","ltime","stime","ptime","matchto","From","Entry","Reviewed","Entry.Name","Protein.names","Gene.Names","Organism","Length","3D","taxnames")

	rid_df_tmp = matrix(data="", ncol=length(blastres_names), nrow=nrow(search_df))
	rid_df_tmp = as.data.frame(rid_df_tmp, stringsAsFactors=FALSE)
	names(rid_df_tmp) = blastres_names

	blast_df_tmp = matrix(data="", ncol=length(blast_fields), nrow=nrow(search_df))
	blast_df_tmp = as.data.frame(blast_df_tmp, stringsAsFactors=FALSE)
	names(blast_df_tmp) = blast_fields

	# Merge subtables
	search_df = cbind(search_df, rid_df_tmp, blast_df_tmp)
	
	# Write out
	if (overwrite == TRUE)
		{
		write.table(x=search_df, file=fasta_blast_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
		} else {
		# Overwrite is FALSE, only write to file if it doesn't exist yet
		if (file.exists(fasta_blast_table_fn) == TRUE)
			{
			warning("WARNING in build_blastsearch_table_from_fasta(): The file fasta_blast_table_fn='", fasta_blast_table_fn, "', already exists, and 'overwrite' is set to FALSE, so we are loading the previous file. Set to overwrite=TRUE to overwrite, but keep in mind this might overwrite previous results.")
			
			# Load from previous file
			search_df = read.table(file=fasta_blast_table_fn, header=TRUE, sep="\t", quote="", row.names=NULL, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors = FALSE)
			} else {
			write.table(x=search_df, file=fasta_blast_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
			}# END if (file.exists(fasta_blast_table_fn) == TRUE)
		} # END if (overwrite == FALSE)
	
	res = NULL
	res$search_df = search_df
	res$fasta_blast_table_fn = fasta_blast_table_fn

	extract='
	search_df = res$search_df
	fasta_blast_table_fn = res$fasta_blast_table_fn
	search_df = read.table(file=fasta_blast_table_fn, header=TRUE, sep="\t", quote="", row.names=NULL, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors = FALSE)
	' # end extract
	
	
	return(res)
	} # END build_blastsearch_table_from_fasta <- function(fasta_fn, type="AA")



#######################################################
# 2. Construct BLAST urls for each
#######################################################
construct_blasturls_for_table <- function(fasta_blast_table_fn, database="nr", hitListSize="10", 
                        filter="L", expect="10", program="blastn",
                        organism=NULL, not_organism=NULL, addl_url=NULL, 
                        txt=NULL, baseUrl="http://blast.ncbi.nlm.nih.gov/blast/Blast.cgi", assume_txid=TRUE)
	{
	ex='
	'
	
	# Read in the table
	search_df = read.table(file=fasta_blast_table_fn, header=TRUE, sep="\t", quote="", row.names=NULL, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors = FALSE)
	dim(search_df)
	
	for (i in 1:nrow(search_df))
		{
		searchseq = search_df$seqstring[i]
		
		taxid = search_df$taxid[i]
		if (is.null(taxid) == FALSE)
			{
			stop("Error in construct_blasturls_for_table(): batch search by taxids not yet programmed.")
			}
		
		results_urls = construct_blasturl(searchseq=searchseq, database=database, hitListSize=hitListSize, 
									 filter=filter, expect=expect, program=program,
									 organism=organism, not_organism=not_organism, addl_url=addl_url, 
									 txt=txt, baseUrl=baseUrl, assume_txid=assume_txid, printout=FALSE)
		
		search_df$RID_url_orig[i] = results_urls$RID_url_orig
		search_df$RID_url[i] = results_urls$RID_url
		search_df$cTF[i] = TRUE
		} # END for (i in 1:nrow(search_df))
	
	# Return the updated table
	# Write out
	write.table(x=search_df, file=fasta_blast_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	
	res = NULL
	res$search_df = search_df
	res$fasta_blast_table_fn = fasta_blast_table_fn

	extract='
	search_df = res$search_df
	fasta_blast_table_fn = res$fasta_blast_table_fn
	search_df = read.table(file=fasta_blast_table_fn, header=TRUE, sep="\t", quote="", row.names=NULL, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors = FALSE)
	' # end extract

	return(res)
	} # END construct_blasturls_for_table




#######################################################
# 3. Getting RIDs for each search
#######################################################
get_some_RIDs_for_table <- function(fasta_blast_table_fn, num_to_get=5, wait_between=1.0, baseUrl="default")
	{
	ex='
	fasta_fn = "133aas_wo_UniProt.fasta"
	fasta_blast_table_fn = "133aas_wo_UniProt_blast_table_v1.txt"
	num_to_get=5
	wait_between=1.0
	'
	
	# Read in the table
	search_df = read.table(file=fasta_blast_table_fn, header=TRUE, sep="\t", quote="", row.names=NULL, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors = FALSE)
	
	# Get the next 5 which don't have RIDs
	
	# blank the NAs
	search_df$rTF[is.na(search_df$rTF)] = FALSE
	TF = search_df$rTF != TRUE
	
	if (sum(TF) == 0)
		{
		res$search_df = search_df
		res$fasta_blast_table_fn = fasta_blast_table_fn
		cat("\nNote: get_RIDs_for_table() says all entries in fasta_blast_table_fn='", fasta_blast_table_fn, "'  have RIDs. Returning res.\n", sep="")
		return(res)		
		}
	
	TFnums = 1:sum(TF)
	rownums = 1:nrow(search_df)
	searchnums = rep(0, times=nrow(search_df))
	searchnums[TF] = TFnums
	
	# Keep the first 5
	TF1 = searchnums > 0
	TF2 = searchnums <= num_to_get
	keepTF = TF1 + TF2 == 2
	rownums_to_add_RIDs = rownums[keepTF]
	
	if (sum(keepTF) > 0)
		{
		# get the RID and waiting time
		rownums_txt = paste0(rownums_to_add_RIDs, sep="", collapse=", ")
		
		cat("\nget_some_RIDs_for_table() is retrieving ", length(rownums_to_add_RIDs), " blast RID numbers, and expected waiting times. The rownums are: ", rownums_txt, " (/", nrow(search_df), " total).\n", sep="")
		# construct the "results" field
		for (rownum in rownums_to_add_RIDs)
			{
			results = NULL
		
			results$RID_url_orig = search_df$RID_url_orig[rownum]
			results$RID_url = search_df$RID_url[rownum]
			
			# Get the RID for this 
			results = get_blast_RID_waiting_time(results=results, RID_url=NULL, baseUrl=baseUrl)
			search_df$rtime[rownum] = Sys.time()
			search_df$RID_url_orig[rownum] = results$RID_url_orig
			search_df$RID_url[rownum] = results$RID_url
			search_df$search_url[rownum] = results$search_url
			search_df$terminal_cmd[rownum] = results$terminal_cmd
			search_df$RIDrequest_fn[rownum] = results$RIDrequest_fn
			#search_df$blastXML_fn[rownum] = results$blastXML_fn # get from download BLAST
			search_df$rid[rownum] = results$rid
			search_df$rtoe[rownum] = results$rtoe
			search_df$rTF[rownum] = TRUE
			
			# Wait for a bit
			Sys.sleep(rnorm(1, wait_between, wait_between/10))
			} # END for (rownum in rownums_to_add_RIDs)
		} # END if (sum(keepTF) > 0)

	# Return the updated table
	# Write out
	write.table(x=search_df, file=fasta_blast_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	
	res = NULL
	res$search_df = search_df
	res$fasta_blast_table_fn = fasta_blast_table_fn

	extract='
	search_df = res$search_df
	fasta_blast_table_fn = res$fasta_blast_table_fn
	search_df = read.table(file=fasta_blast_table_fn, header=TRUE, sep="\t", quote="", row.names=NULL, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors = FALSE)
	' # end extract

	return(res)
	} # END get_some_RIDs_for_table



#######################################################
# 4. Launching each search
#######################################################
launch_some_blast_searches <- function(fasta_blast_table_fn, num_to_launch=5, wait_between=2.0)
	{
	ex='
	fasta_fn = "133aas_wo_UniProt.fasta"
	fasta_blast_table_fn = "133aas_wo_UniProt_blast_table_v1.txt"
	num_to_launch=5
	wait_between=1.0
	'
	
	# Read in the table
	search_df = read.table(file=fasta_blast_table_fn, header=TRUE, sep="\t", quote="", row.names=NULL, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors = FALSE)
	
	# Get the next 5 which have RIDs
	
	# blank the NAs
	search_df$rTF[is.na(search_df$rTF)] = FALSE
	search_df$lTF[is.na(search_df$lTF)] = FALSE
	TF1 = search_df$rTF == TRUE
	TF2 = search_df$lTF == FALSE
	TF = (TF1 + TF2) == 2
	
	if (sum(TF) == 0)
		{
		res$search_df = search_df
		res$fasta_blast_table_fn = fasta_blast_table_fn
		cat("\nNote: launch_some_blast_searches() finds 0 entries in fasta_blast_table_fn='", fasta_blast_table_fn, "' with the results RIDs generated. Run get_some_RIDs_for_table() to generate some RIDs. Returning res.\n", sep="")
		return(res)		
		}
	
	# Otherwise
	TFnums = 1:sum(TF)
	rownums = 1:nrow(search_df)
	searchnums = rep(0, times=nrow(search_df))
	searchnums[TF] = TFnums
	
	# Keep the first 5
	TF1 = searchnums > 0
	TF2 = searchnums <= num_to_launch
	keepTF = TF1 + TF2 == 2
	rownums_to_launch_blast_searches = rownums[keepTF]
	
	if (sum(keepTF) > 0)
		{
		# get the RID and waiting time
		rownums_txt = paste0(rownums_to_launch_blast_searches, sep="", collapse=", ")
		
		cat("\nlaunch_some_blast_searches() is launching ", length(rownums_to_launch_blast_searches), " blast searches. The rownums are: ", rownums_txt, " (/", nrow(search_df), " total).\n", sep="")
		# construct the "results" field
		for (rownum in rownums_to_launch_blast_searches)
			{
			# Assemble the results object to run 
			results = NULL
			results$RID_url_orig = search_df$RID_url_orig[rownum]
			results$RID_url = search_df$RID_url[rownum]
			results$search_url = search_df$search_url[rownum]
			results$terminal_cmd = search_df$terminal_cmd[rownum]
			results$RIDrequest_fn = search_df$RIDrequest_fn[rownum]
			results$blastXML_fn = search_df$blastXML_fn[rownum]
			results$rid = search_df$rid[rownum]
			results$rtoe = search_df$rtoe[rownum]
			
			# Nothing returned; this just launches the web search
			start_BLAST_websearch(results=results, terminal_cmd=NULL)
			search_df$ltime[rownum] = Sys.time()
			search_df$lTF[rownum] = TRUE

			# Wait for a bit
			Sys.sleep(rnorm(1, wait_between, wait_between/5))
			} # END for (rownum in rownums_to_add_RIDs)
		} # END if (sum(keepTF) > 0)

	# Return the updated table
	# Write out
	write.table(x=search_df, file=fasta_blast_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	
	res = NULL
	res$search_df = search_df
	res$fasta_blast_table_fn = fasta_blast_table_fn

	extract='
	search_df = res$search_df
	fasta_blast_table_fn = res$fasta_blast_table_fn
	search_df = read.table(file=fasta_blast_table_fn, header=TRUE, sep="\t", quote="", row.names=NULL, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors = FALSE)
	' # end extract

	return(res)
	}




get_blastXML_fn_from_search_df <- function(search_df, rownum)
	{
	rownum_txt = sprintf(fmt="%05.0f", rownum)
	rid = get_RID_from_url(search_url=search_df$search_url[rownum])
	gid = search_df$seqid[rownum]

	if (is.null(gid) == TRUE)
		{
		gid_txt = ""
		} else {
		gid_txt = paste0("_gid_", gid)
		}
	
	if (is.na(rid) == FALSE)
		{
		blastXML_fn = paste0(current_time_txt_wd(rownum_txt), "_RID_", rid, gid_txt, "_blastres_download.xml")
		} else {
		blastXML_fn = paste0(current_time_txt_wd(rownum_txt), "_RID_unk", gid_txt, "_blastres_download.xml")
		}
	
	return(blastXML_fn)
	}

# rid = get_RID_from_url(search_url=search_df$search_url[rownum])
get_blastXML_fn_from_search_df2 <- function(rownum, gid, rid)
	{
	rownum_txt = sprintf(fmt="%05.0f", rownum)
	#rid = get_RID_from_url(search_url=search_df$search_url[rownum])
	#gid = search_df$seqid[rownum]

	if (is.null(gid) == TRUE)
		{
		gid_txt = ""
		} else {
		gid_txt = paste0("_gid_", gid)
		}
	
	if (is.na(rid) == FALSE)
		{
		blastXML_fn = paste0(current_time_txt_wd(rownum_txt), "_RID_", rid, gid_txt, "_blastres_download.xml")
		} else {
		blastXML_fn = paste0(current_time_txt_wd(rownum_txt), "_RID_unk", gid_txt, "_blastres_download.xml")
		}
	
	return(blastXML_fn)
	}



#######################################################
# 5. Checking if they are done & downloading, if so
#######################################################
download_finished_blast_searches <- function(fasta_blast_table_fn, num_to_get=5, wait_between=1.0, min_time_since_launch_before_download=1*60, blastXML_fn="default", type="AA")
	{
	ex='
	fasta_fn = "133aas_wo_UniProt.fasta"
	fasta_blast_table_fn = "133aas_wo_UniProt_blast_table_v1.txt"
	num_to_launch=5
	wait_between=1.0
	'
	
	# Read in the table
	search_df = read.table(file=fasta_blast_table_fn, header=TRUE, sep="\t", quote="", row.names=NULL, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors = FALSE)
	
	# Get the next 5 which have RIDs
	
	# blank the NAs, check for lines where launch is true AND save is false
	search_df$lTF[is.na(search_df$lTF)] = FALSE
	search_df$sTF[is.na(search_df$sTF)] = FALSE
	search_df$eTF[is.na(search_df$eTF)] = FALSE
	TF1 = search_df$lTF == TRUE
	TF2 = search_df$sTF != TRUE
	TF3 = search_df$eTF != TRUE
	TF = TF1 + TF2 + TF3 == 3
	sum(TF)

	# Check those for if the minimum time since launch has expired
	possible_rows = (1:length(TF))[TF]
	time_since_launch = as.numeric(Sys.time()) - as.numeric(search_df$ltime[possible_rows])
	time_since_launch_TF = rep(FALSE, times=length(TF))
	time_since_launch_TF[time_since_launch > min_time_since_launch_before_download] = TRUE
	TF[time_since_launch_TF == FALSE] = FALSE
	
	if (sum(TF) == 0)
		{
		res$search_df = search_df
		res$fasta_blast_table_fn = fasta_blast_table_fn
		cat("\nNote: download_finished_blast_searches() finds 0 entries in fasta_blast_table_fn='", fasta_blast_table_fn, "' with launched results not yet downloaded. Run get_some_RIDs_for_table() and launch_some_blast_searches() to generate some RIDs to download. Returning res.\n", sep="")
		return(res)		
		}
	
	# Otherwise
	TFnums = 1:sum(TF)
	rownums = 1:nrow(search_df)
	searchnums = rep(0, times=nrow(search_df))
	searchnums[TF] = TFnums
	
	# Keep the first 5
	TF1 = searchnums > 0
	TF2 = searchnums <= num_to_get
	keepTF = TF1 + TF2 == 2
	rownums_to_try_downloads = rownums[keepTF]
	rownums_to_try_downloads
	
	blastxmls_df = NULL
	if (sum(keepTF) > 0)
		{
		# get the RID and waiting time
		rownums_txt = paste0(rownums_to_try_downloads, sep="", collapse=", ")
		
		cat("\nlaunch_some_blast_searches() is launching ", length(rownums_to_try_downloads), " blast searches. The rownums are: ", rownums_txt, " (/", nrow(search_df), " total).\n", sep="")
		# construct the "results" field
		for (rownum in rownums_to_try_downloads)
			{
			cat("\nCheckpoint 1: Attempting download of rownum=", rownum, "...", sep="")
			
			

			if (blastXML_fn == "default")
				{
				blastXML_fn = get_blastXML_fn_from_search_df(search_df, rownum)
				tryParse_results_justXML_ptr_or_error = download_BLAST_xml(results=NULL, search_url=search_df$search_url[rownum], attempts=1, retry_wait=1, blastXML_fn=blastXML_fn)
				catn("done2.4")
				} else {
				tryParse_results_justXML_ptr_or_error = NULL
				
				
				tryParse_results_justXML_ptr_or_error$blastXML_fn = blastXML_fn
				catn("done2.5")
				}
			
			catn("done2.6")

			# This can cause this error:
			# grepl(pattern="ERROR in download_BLAST_xml", x=tryParse_results) == FALSE
			# 
			# Error in as.vector(x, "character") : 
      #    cannot coerce type 'externalptr' to vector of type 'character'
			
			print(tryParse_results_justXML_ptr_or_error)
			#print(tryParse_results$blastXML_fn) # tryParse_results is just the XML results, not subsettable
			# Error in tryParse_results$blastXML_fn : 
  		# 	object of type 'externalptr' is not subsettable

			print(length(tryParse_results_justXML_ptr_or_error))
			print(names(tryParse_results_justXML_ptr_or_error))
			print(dput(tryParse_results_justXML_ptr_or_error))
			#grepl(pattern="ERROR in download_BLAST_xml", x=tryParse_results) == FALSE
			
			catn("done2.7")
			
			print("blastXML_fn")
			print(blastXML_fn)
			
			catn("done2.8")
			
			# Some kind of error
			if ((length(tryParse_results_justXML_ptr_or_error) == 1) && (grepl(pattern="ERROR in download_BLAST_xml", x=tryParse_results_justXML_ptr_or_error) == FALSE))
				{
				#print(tryParse_results$result)
				search_df$sTF[rownum] = FALSE
				search_df$stime[rownum] = Sys.time()
				} else {
				search_df$sTF[rownum] = TRUE
				search_df$stime[rownum] = Sys.time()
				cat("...save successful.\n")
				}

			# If no initial error, process downloaded results
			if (search_df$sTF[rownum] == TRUE)
				{
				#blastXML_fn = tryParse_results$blastXML_fn
				
				# Save the filename
				catn("done2.9")
				search_df$blastXML_fn[rownum] = blastXML_fn
				checks = check_saved_BLAST_XMLs(blastXML_fn=blastXML_fn)
				checks
				catn("done3.0")
				
				# Check for downloaded, but not finished
				if (checks$still_updating == TRUE)
					{
					search_df$sTF[rownum] = FALSE
					search_df$stime[rownum] = Sys.time()
					search_df$pTF[rownum] = FALSE
					search_df$eTF[rownum] = FALSE
					search_df$ncbiWarn[rownum] = checks$NCBI_CPUtime_warning
					
					# Remove the old updating file, if it exists
					print("rownum")
					print("search_df$blastXML_fn[rownum]")
					print(search_df$blastXML_fn[rownum])
					print("file.exists(search_df$blastXML_fn[rownum]")
					print(file.exists(search_df$blastXML_fn[rownum]))
					print("")
					print("blastXML_fn")
					print(blastXML_fn)
					print("file.exists(blastXML_fn)")
					print(file.exists(blastXML_fn))
					print("")
					
					if (file.exists(search_df$blastXML_fn[rownum]) == TRUE)
						{
						file.remove(search_df$blastXML_fn[rownum])
						search_df$blastXML_fn[rownum] = blastXML_fn
						} # END if (file.exists(search_df$blastXML_fn[rownum]) == TRUE)

					next()
					}
				
				cat("\nCheckpoint 2: Attempting parsing of rownum=", rownum, "...", sep="")
				tryParse_results=tryParse_results_justXML_ptr_or_error
				try_processing_result = try(process_blastresult(tryParse_results=tryParse_results_justXML_ptr_or_error, genera=NULL, type=type))

				if (class(try_processing_result) == "try-error")
					{
					catn("check6")
					search_df$pTF[rownum] = FALSE
					search_df$eTF[rownum] = TRUE
		
					warning_txt = paste0("\nWARNING: process_blastresult() on rownum=", rownum, " gave a try-error on rownum=", rownum, ". Inserting a blank line of NAs in blastres$seqdf.\n")
		
					cat("\n\n")
					cat(warning_txt)
					cat("\n\n")
		
					headers = c("Hit_num", "downID", "taxon", "gi1", "gi2", "gi3", "gi4", "pctIdent", "num_identical", "num_disagree", "num_resolved", "aln_len", "genus", "species", "subspecies", "gnsp", "gnspssp", "pctIdent", "pctSim", "mt", "numt", "cp", "cgenome", "pgenome", "pCDS", "partial", "Hit_id", "Hit_def", "Hit_accession", "accession2", "Hit_len", "Hsp_num", "Hsp_bit_score", "Hsp_score", "Hsp_evalue", "Hsp_query_from", "Hsp_hit_to", "Hsp_query_frame", "Hsp_hit_frame", "Hsp_identity", "Hsp_positive", "Hsp_gaps", "Hsp_align_len", "qseq", "hseq")
		
					blank_line = rep(NA, times=length(headers))
					blank_line_df = as.data.frame(matrix(blank_line, nrow=1), stringsAsFactors=FALSE)
					names(blank_line_df) = headers
					blank_line_df
		
					blastres = NULL
					blastres$aln = list()
					blastres$seqdf = blank_line_df				
					} else {
					catn("check7")
					# If successful...
					blastres = try_processing_result
					search_df$pTF[rownum] = TRUE
					search_df$eTF[rownum] = FALSE

					# Number of hits once processed
					number_of_blast_hits = length(blastres$seqdf$Hit_num) - sum(is.na(blastres$seqdf$Hit_num))
					search_df$pTF_numhits[rownum] = number_of_blast_hits

					cat("...parsing of rownum=", rownum, " successful.\n", sep="")
					} # END if (class(try_result) == "try-error")
				} # END if (search_df$sTF[rownum] == TRUE)
			
			if (search_df$sTF[rownum] == TRUE)
				{
				cat("\nCheckpoint 3: Storing parsed results for rownum=", rownum, ".\n", sep="")
				
				# Construct output table
				RID = search_df$rid[rownum]
				seqid_to_search = search_df$seqid[rownum]
				seqname = search_df$seqnames[rownum]
				taxnames = search_df$taxname[rownum]
				searchseq = search_df$seqstring[rownum]
			
				iS = rep(rownum, times=nrow(blastres$seqdf))
	#				abbrS = rep(abbr[rownum], times=nrow(blastres$seqdf))
				RIDsS = rep(RID, times=nrow(blastres$seqdf))
				seqids_to_searchS = rep(seqid_to_search, times=nrow(blastres$seqdf))
				seqnames_to_searchS = rep(seqname, times=nrow(blastres$seqdf))
				taxnamesS = rep(taxnames, times=nrow(blastres$seqdf))
				searchseqS = rep(searchseq, times=nrow(blastres$seqdf))

				outdf1 = cbind(iS, RIDsS, seqids_to_searchS, seqnames_to_searchS, taxnamesS, searchseqS)
				outdf2 = as.data.frame(outdf1, stringsAsFactors=FALSE)
				names(outdf2) = c("rownum", "RID", "seqid", "seqname", "taxname", "searchseq")
				tmpdf = cbind(outdf2, blastres$seqdf)
				blast_results_lines_fn = paste0("blast_results_rownum", sprintf(fmt="%05.0f", rownum), "_", seqid_to_search, ".txt")
				write.table(x=tmpdf, file=blast_results_lines_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
				search_df$blastres_fns[rownum] = blast_results_lines_fn

				blastxmls_df = rbind(blastxmls_df, tmpdf)

				# Wait for a bit
				Sys.sleep(rnorm(1, wait_between, wait_between/5))
				} # END if (search_df$sTF[rownum] == TRUE)
			} # END for (rownum in rownums_to_add_RIDs)
		} # END if (sum(keepTF) > 0)

	# Return the updated table
	# Write out
	write.table(x=search_df, file=fasta_blast_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	
	res = NULL
	res$search_df = search_df
	res$fasta_blast_table_fn = fasta_blast_table_fn
	res$blastxmls_df = blastxmls_df

	extract='
	search_df = res$search_df
	fasta_blast_table_fn = res$fasta_blast_table_fn
	blastxmls_df = res$blastxmls_df
	search_df = read.table(file=fasta_blast_table_fn, header=TRUE, sep="\t", quote="", row.names=NULL, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors = FALSE)
	' # end extract

	return(res)
	} # END download_finished_blast_searches <- function(fasta_blast_table_fn, num_to_get=5, wait_between=1.0, min_time_since_launch_before_download=1*60, blastXML_fn="default", type="AA")


#######################################################
# Run steps 1-6
#######################################################
fill_out_table_while_loop <- function(fasta_blast_table_fn, wait_between_loops=30, max_num_searches=5)
	{
	ex='
	fasta_blast_table_fn = "/GitHub/str2phy/ex/MitoCOGs/_03_alphafolds/44-88aa_of_133_wo_UniProt_blast_table_v1.txt"
	wait_between_loops=30
	max_num_searches=1
	'
	
	start_time = Sys.time()
	wait_time = round(runif(1, min=wait_between_loops-(wait_between_loops/8), max=wait_between_loops+wait_between_loops/8), digits=1)

	# Run this loop until number_rows_complete equals total_rows
	#wait_between_loops = 30
	while(matdf["number_rows_complete",] < matdf["total_rows",])
		{
		matdf = check_num_still_running(fasta_blast_table_fn)
	
		if (matdf["number_launched_not_downloaded",] > 0)
			{
			res = download_finished_blast_searches(fasta_blast_table_fn, num_to_get=1, wait_between=1.0)
			search_df = res$search_df
			Sys.sleep(wait_time)
			}
	
		matdf = check_num_still_running(fasta_blast_table_fn)
		print(matdf)

		#max_num_searches = 5
		searches_to_add = max_num_searches - matdf["number_launched_not_downloaded",]
		if (searches_to_add > 0)
			{
			res = get_some_RIDs_for_table(fasta_blast_table_fn, num_to_get=searches_to_add, wait_between=runif(1,1.1,2.5), baseUrl="default")
			res = launch_some_blast_searches(fasta_blast_table_fn, num_to_launch=searches_to_add, wait_between=rnorm(n=1, mean=3.8, sd=0.9))
			search_df = res$search_df
			head(search_df[,1:7])
			}

		matdf = check_num_still_running(fasta_blast_table_fn)
		print(matdf)
	
		cat("\nWaiting for wait_between_loops=~", wait_time, " seconds.\n", sep="")	
		} # END while-loop

	matdf = check_num_still_running(fasta_blast_table_fn)
	
	stop_time = Sys.time()
	time_elapsed = as.character(as.numeric(stop_time-start_time))
	
	tmpmat = matrix(data=c(as.character(start_time), as.character(stop_time), time_elapsed), ncol=1)
	tmpmat = as.data.frame(tmpmat, stringsAsFactors=FALSE)
	names(tmpmat) = "count"
	row.names(tmpmat) = c("start_time", "stop_time", "time_elapsed")
	
	matdf2 = cft(rbind(matdf, tmpmat))
	return(matdf2)
	} # END fill_out_table_while_loop



#######################################################
# Re-download possible errors (parsing or formal error 
# codes)
#######################################################
refresh_all_downloads_errors_in_blast_table <- function(fasta_blast_table_fn)
	{
	search_df = read_search_df(fasta_blast_table_fn)
	# Just set all saves (and errors/processing flags) to FALSE, to re-do downloads
	search_df$sTF = FALSE
	search_df$pTF = FALSE
	search_df$eTF = FALSE
	write_search_df(search_df, fasta_blast_table_fn)
	return(search_df)
	}

refresh_possible_errors_in_blast_table <- function(fasta_blast_table_fn)
	{
	search_df = read_search_df(fasta_blast_table_fn)
	search_df$pTF_numhits[is.na(search_df$pTF_numhits)] = 0 # reset NAs to 0
	search_df$sTF[search_df$pTF_numhits < 1] = FALSE
	search_df$pTF[search_df$pTF_numhits < 1] = FALSE
	search_df$eTF[search_df$pTF_numhits < 1] = FALSE
	write_search_df(search_df, fasta_blast_table_fn)
	return(search_df)
	}




#######################################################
# Check the number of searches still outstanding
#######################################################
check_num_still_running <- function(fasta_blast_table_fn)
	{
	# Read in the table
	search_df = read.table(file=fasta_blast_table_fn, header=TRUE, sep="\t", quote="", row.names=NULL, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors = FALSE)
	
	# blank the NAs
	search_df$rTF[is.na(search_df$rTF)] = FALSE
	search_df$lTF[is.na(search_df$lTF)] = FALSE
	search_df$sTF[is.na(search_df$sTF)] = FALSE
	search_df$pTF[is.na(search_df$pTF)] = FALSE
	search_df$eTF[is.na(search_df$eTF)] = FALSE
	TF = search_df$lTF == TRUE
	total_rows = length(TF)
	number_launched = sum(TF)
	number_not_launched = total_rows - sum(TF)

	TF = search_df$sTF == TRUE
	number_downloaded = sum(TF)
	
	# Launched not saved
	TF1 = search_df$lTF == TRUE
	TF2 = search_df$sTF != TRUE
	TF = (TF1 + TF2) == 2
	number_launched_not_downloaded = sum(TF)
	
	TF = search_df$sTF == TRUE
	number_downloaded = sum(TF)

	TF = search_df$eTF == TRUE
	number_downloaded_but_error = sum(TF)
	
	# Rows complete
	TF1 = search_df$pTF == TRUE
	TF2 = search_df$eTF == TRUE
	TF = (TF1 + TF2) >= 1
	number_rows_complete = sum(TF)
	
	item = c("total_rows", "number_launched", "number_not_launched", "number_downloaded", "number_launched_not_downloaded", "number_downloaded_but_error", "number_rows_complete")
	mat = matrix(data=c(total_rows, number_launched, number_not_launched, number_downloaded, number_launched_not_downloaded, number_downloaded_but_error, number_rows_complete), ncol=1)
	matdf = as.data.frame(mat, stringsAsFactors=FALSE)
	names(matdf) = "count"
	row.names(matdf) = item
	
	return(matdf)
	} # END check_num_still_running <- function(fasta_blast_table_fn)





# Using the blastSequences function
# 
# The blastSequences function appears to work well in most instances. For slow web 
# connections or for large query sequences, it appears to fail. This seem to be due 
# to either a low number of attempts to retrieve hits from the ncbi server, or the 
# short time the routine waits before it throws an error in R. It seem that if we 
# increase the number of attempts to retrieve hits, generally, we can successfully 
# run our BLAST.
# The below code re-defines the blastSequences function, with an additional argument 
# (attempts), to allow users to control how many times the routine should try to 
# retrieve results from the server.
# 
# http://rstudio-pubs-static.s3.amazonaws.com/12097_1352791b169f423f910d93222a4c2d85.html
#
# See also: http://www.inside-r.org/packages/cran/BoSSA/docs/blast
# 
# function definition
blastSeqKK2 <- function (searchseq, database="nr", hitListSize="10", 
                        filter="L", expect="10", program="blastn",
                        organism=NULL, not_organism=NULL, addl_url=NULL, 
                        txt=NULL, justRID_url=FALSE, attempts=10, retry_wait=5, baseUrl="http://blast.ncbi.nlm.nih.gov/blast/Blast.cgi")
#baseUrl="http://www.ncbi.nlm.nih.gov/blast/Blast.cgi")
	{
	defaults='
	#wd = "/drives/Dropbox/_njm/"
	#setwd(wd)
	fasta_fn = "/drives/Dropbox/_njm/Wayne_etal_1997_SysBio_canid_mtDNA_appendix_wolf_cytB.fasta"
	mySeq = read.fasta(fasta_fn)

	searchseq = convert_seqs_to_seqCollapse(x=mySeq[[1]])
	database="nr"
	hitListSize="10"
	filter="L"
	expect="10"
	program="blastp"
	attempts=10
	organism=NULL
	not_organism = NULL
	addl_url=NULL
	txt=NULL
	justRID_url=FALSE
	attempts=10
	retry_wait=5
	baseUrl = "https://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
	#baseUrl = "uoa"
	
	
	results = blastSeqKK2(searchseq=searchseq, database="nr", hitListSize="10", 
                        filter="L", expect="10", program="blastn",
                        organism=NULL, not_organism=NULL, addl_url=NULL, 
                        txt=txt, justRID_url=FALSE, attempts=10, retry_wait=5, baseUrl="default")
	' # END defaults
	
	# Default
	# baseUrl="http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
	# Other, non-blast, baseUrls
	# http://eutils.ncbi.nlm.nih.gov/blast/Blast.cgi
	# If you are on a university campus that auto-replaces the URL, you will have to customize this:
	if (baseUrl == "uoa")
		{
		baseUrl = "https://blast-ncbi-nlm-nih-gov.ezproxy.auckland.ac.nz/Blast.cgi"
		}
	if (baseUrl == "default")
		{
		#baseUrl = "https://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
		baseUrl = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
		}
	
	if (is.null(filter) == TRUE)
		{
		query = paste("QUERY=", as.character(searchseq), "&DATABASE=", database, 
			 "&HITLIST_SIZE=", hitListSize, "&EXPECT=", 
			 expect, "&PROGRAM=", program, sep = "")
		} else {
		# Include filter
		query = paste("QUERY=", as.character(searchseq), "&DATABASE=", database, 
			 "&HITLIST_SIZE=", hitListSize, "&FILTER=", filter, "&EXPECT=", 
			 expect, "&PROGRAM=", program, sep = "")
		
		}
	
	# Put in organism limitation
	entrez_query = FALSE
	q2 = ""
	if (!is.null(organism) || !is.null(not_organism))
		{
		entrez_query = TRUE
		
		# Add e.g. "human[orgn]" to query
		# http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node25.html
		if (length(c(organism, not_organism)) == 1)
			{
			if (length(organism) == 1)
				{
				q2 = paste0(q2, organism, "[orgn]&")
				}
			if (length(not_organism) == 1)
				{
				q2 = paste0(q2, "NOT ", not_organism, "[orgn]&")
				}
			} else {
			
			if ( (length(organism) > 1) || (length(not_organism) > 1) )
				{
				if (length(organism) > 0)
					{
					orgn_txt = rep("[orgn]", times=length(organism))
					organisms_txt1 = paste(unlist(organism), orgn_txt, sep="")
					organisms_txt1
					organisms_txt2 = paste(organisms_txt1, collapse=" OR ", sep="")
					organisms_txt2
					} else {
					organisms_txt2 = NULL
					}
				if (length(not_organism) > 0)
					{
					NOT_txt = rep("NOT ", times=length(not_organism))
					orgn_txt = rep("[orgn]", times=length(not_organism))
					not_organisms_txt1 = paste(NOT_txt, unlist(not_organism), orgn_txt, sep="")
					not_organisms_txt1
					not_organisms_txt2 = paste(not_organisms_txt1, collapse=" ", sep="")
					not_organisms_txt2
					} else {
					not_organisms_txt2 = NULL
					}
				}
			
			txt3 = paste0("(", organisms_txt2, " ", not_organisms_txt2, ")")
			q2 = paste0(q2, txt3, "&")
			}
		} # END if (!is.null(organism) || !is.null(not_organism))


	# Put in additional (manual) limitations into the URL string
	if (!is.null(addl_url))
		{
		entrez_query = TRUE
		
		# Add e.g. "human[orgn]" to query
		# http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node25.html
		q2 = paste0(q2, addl_url, "&")
		} # END if (!is.null(addl_url))

	# Put in additional (text, in ENTREZ format) limitations on the search
	if (!is.null(txt))
		{
		entrez_query = TRUE
		
		# Convert to URL format...
		#addl_txt = convert_txt_to_NCBI(txt)
		# Add e.g. "human[orgn]" to query
		# http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node25.html
		q2 = paste0(q2, txt, "&")
		} # END if (!is.null(addl_url))
	
	# If you've added anything anywhere, add it to "query"
	if (entrez_query == TRUE)
		{
		q2 = paste0("&ENTREZ_QUERY=", q2)
		query = paste0(query, q2)
		} # END if (entrez_query = TRUE)
	q2

	# Convert the URL to HTML format
	query = gsub(pattern="&&", replacement="&", x=query)
	query_orig = query
	query = convert_txt_to_NCBI(query_orig)		
	RID_url = sprintf("%s?%s&CMD=Put", baseUrl, query)
	RID_url_orig = sprintf("%s?%s&CMD=Put", baseUrl, query_orig)
	RID_url = gsub(pattern="&&", replacement="&", x=RID_url)
	RID_url_orig = gsub(pattern="&&", replacement="&", x=RID_url_orig)
	
	cat("\n\nThe search URL, 'results$RID_url' has been assembled:\n\n")
	cat(RID_url_orig)
	cat("\n")
	
	# Just get the URL?
	if (justRID_url == TRUE)
		{
		results = NULL
		results$RID_url_orig = RID_url_orig
		results$RID_url = RID_url
		return(results)
		} else {
		results = NULL
		}
	
	# Continue to retrive resultsif desired
	#results = tempfile()
	
	# Wait 1-2 seconds (if needed)
	Sys.sleep(runif(n=1,min=0.0, max=2.0))
		
	require(XML)
	
	# Get the search RID
	cat("\nDownloading HTML from NCBI 'Format Request' page...")
	destfile = "temp_blast_search_request.html"
	tmp = download.file(RID_url, destfile=destfile)
	cat("done3.\n")
	#moref(destfile)

	cat("Extracting RID (Request ID)... RID=")
	post = htmlTreeParse(destfile, useInternalNodes = TRUE)
	x = post[["string(//comment()[contains(., \"QBlastInfoBegin\")])"]]
	rid = sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
	cat(rid)
	cat("\n")
	
	
	# Put the search RID into a url for searching/downloading
	# RTOE = Estimated waiting time before the BLAST results will be ready
	rtoe = as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1", x))
	
	cat("Estimated waiting time for search is rtoe=")
	cat(as.character(rtoe))
	cat(" seconds.\n")
	
	# Search 
	search_url = sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl, rid)
	
	cat("Retrieving RID=", rid, " with this terminal command (will open a browser window):\n", sep="")
	cat(search_url)
	cat("\n...")
	
	terminal_cmd = paste0("open -n '", search_url, "'")
	system(terminal_cmd)
	
	cat("counting down expected 'rtoe' waiting time...")
	for (i in rev(seq(rtoe)))
		{
    cat("Countdown:", i, "s\n")
    Sys.sleep(1)
  	}
	cat("\n...done4. Will now try to download the completed XML file for parsing.\n")
	
	# Wait for the waiting time, plus a bit
	#if (justRID_url == FALSE)
	#	{
	#	Sys.sleep(rtoe + runif(n=1, min=retry_wait, max=2*retry_wait))
	#	}
	
	# Try to parse
	results$RID_url_orig = RID_url_orig
	results$RID_url = RID_url
	results$search_url = search_url
	
	rid = get_RID_from_url(search_url=search_url)
	gid = NULL
	blastXML_fn = get_blastXML_fn_from_search_df2(rownum, rid, gid)
	result = .tryParseResult(search_url, attempts, retry_wait=retry_wait, blastXML_fn=blastXML_fn)
	
	cat("done5.\n")
	
	results = NULL
	if (justRID_url == FALSE)
		{
		results$result = result
		} # END if (justRID_url == FALSE)

	return(results)
	} # END blastSeqKK




#######################################################
# Check for the type of file saved - a "will be updated" HTML, or a results XML?
#######################################################
check_saved_BLAST_XMLs <- function(blastXML_fn)
	{
	ex='
	library(openxlsx)
	library(ape)
	library(phytools)
	library(BioGeoBEARS)
	library(varhandle) # for check.numeric
	library(rentrez) # Fetching full records: entrez_fetch()
	library(XML)
	library(Biostrings) # For parsing BLAST searches downloaded

	#sourceall("/GitHub/bioinfRhints/Rsrc/")
	sourceall("/GitHub/str2phy/Rsrc/")
	source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")

	wd = "/GitHub/str2phy/Rsrc/blast_parsing_fns/"
	setwd(wd)
	
	blastXML_fn = "will_be_automatically_updated_2025-01-18_17.01.22_SP9SBRUP016_RID_blastres_download.xml"
	checks = check_saved_BLAST_XMLs(blastXML_fn)
	checks

	blastXML_fn = "BLAST_result_with_no_hits_2025-01-18_17.23.33_SP9SBRUP016_RID_blastres_download.xml"
	checks = check_saved_BLAST_XMLs(blastXML_fn)
	checks

	blastXML_fn = "BLAST_result_with_hits_2025-01-18_10.31.52_SNJPB69S016_temp_blast_search_request.xml"
	checks = check_saved_BLAST_XMLs(blastXML_fn)
	checks

	blastXML_fn = "html_parsing_error_2025-01-18_14.37.17_temp_blast_search_request.html"
	checks = check_saved_BLAST_XMLs(blastXML_fn)
	checks
	'
	
	checks = NULL
	checks$this_is_HTML_request_ID = FALSE
	checks$NCBI_CPUtime_warning = FALSE
	checks$still_updating = FALSE
	checks$BlastOutput_iterations = FALSE
	checks$BlastOutput_hits = FALSE
	
	lines = readLines(blastXML_fn)
	for (line in lines)
		{
		# Check for the NCBI warning about the IP address using too many resources
		checks$NCBI_CPUtime_warning = grepl(pattern="Searches from this IP address have consumed a large amount", x=line, ignore.case=TRUE)

		# Automatically updating
		TF = grepl(pattern='<label for="rid">Request ID</label>', x=line, ignore.case=TRUE)
		if (TF == TRUE)
			{
			checks$this_is_HTML_request_ID = TRUE
			return(checks)
			}

		
		# Automatically updating
		TF = grepl(pattern="This page will be automatically updated in", x=line, ignore.case=TRUE)
		if (TF == TRUE)
			{
			checks$still_updating = TRUE
			return(checks)
			}

		# Blast output
		TF = grepl(pattern="BlastOutput_iterations", x=line, ignore.case=TRUE)
		if (TF == TRUE)
			{
			checks$BlastOutput_iterations = TRUE
			# don't return, since you still want to check for Hit_num
			}

		# Blast output - no hits
		TF = grepl(pattern="<Hit_num>", x=line, ignore.case=TRUE)
		if (TF == TRUE)
			{
			checks$BlastOutput_hits = TRUE
			return(checks)
			}

		} # END for (line in lines)
	return(checks)
	}




# Setup a try function to download the result (which may still be coming in)
# Make sure to add the sequence ID
.tryParseResult <- function(search_url, attempts=1, retry_wait=1, blastXML_fn="temp.xml", already_saved=FALSE)
	{
	ex = '
	#search_url = ????????
	attempts=1
	retry_wait=1
	blastXML_fn="temp.xml"
	
	'
	
	
	total_pause_time = 0
	i=1
	
	rid = get_RID_from_url(search_url)
	
	# Try a number of attempts
	for (i in 1:attempts)
		{
		cat("\nTrying to download BLAST search result with...\n")
		txt = paste0("tmp = download.file(url='", search_url, "', destfile='", blastXML_fn, "')")
		cat(txt)
		cat("\n")
		
		starttime = Sys.time()
		tmp_result = NULL
		tmp_result = tryCatch(
			{
			#destxmlfile = "temp_blast_search_request.xml"
			
			if (already_saved == FALSE)
				{
				# Download
				tmp = download.file(url=search_url, destfile=blastXML_fn)
				}
#			Sys.sleep(2.0)
			
			# Check if its a "This page will be automatically updated in" file
			checks = check_saved_BLAST_XMLs(blastXML_fn)
			still_running_automatic_updates_TF = checks$still_updating
			
			# validate=FALSE was helpful for reading in downloaded BLAST XML results
			if (still_running_automatic_updates_TF == FALSE)
				{		
				tmp_result2 = xmlTreeParse(file=blastXML_fn, validate=FALSE, useInternalNodes=TRUE, isSchema=FALSE, error = xmlErrorCumulator(immediate=FALSE), isHTML=FALSE)
				#return(tmp_result2) ## Is this the problem?
				} else {
				warning(paste0("RID=", rid, " is still automatically updating."))
				tmp_result2 = NULL
				#tmp_result2$blastXML_fn = blastXML_fn
				#return(tmp_result2) # returns to tmp_result # Is this the problem?
				} # END if (still_running_automatic_updates_TF == FALSE)
			#result = htmlTreeParse(blastXML_fn, useInternalNodes=TRUE, error = xmlErrorCumulator(immediate=FALSE))
			}, # END of tryCatch, except for errors
			error = function(e){ message("An error occurred:\n", e) },
			warning = function(w){ message("A warning occured:\n", w) },
			finally = { message("...XML parsing finished.") } 
			#, error=function(err) NULL  # This last bit returns NULL
		) # END tmp_result = tryCatch


		endtime = Sys.time()
		search_time = endtime - starttime
		cat("\n...attempt finished after ", round(search_time,4), " seconds.", sep="")
		total_pause_time = total_pause_time + search_time
		
		# Wait if needed
		if (!is.null(tmp_result))
			{
			blast_download_summary_txt = paste0("Online BLAST search returned after ", i, " tries, totalling ", round(total_pause_time,4), " seconds.\nResult XML has summary(result)$numNodes=", summary(tmp_result)$numNodes, " nodes, ", summary(tmp_result)$nameCounts["Hit"], " hits, and ", summary(tmp_result)$nameCounts["Hsp"], " Hsps.\n")
			cat(blast_download_summary_txt)
			
			# Return results
			tryParse_results = NULL
			tryParse_results$result = tmp_result
			tryParse_results$blast_download_summary_txt = blast_download_summary_txt
			tryParse_results$blastXML_fn = blastXML_fn
			tryParse_results$total_pause_time = total_pause_time
			
			return(tryParse_results)
			} # END if (!is.null(result))
			
		# Wait longer after each try, just in case it comes through...
		if (is.null(tmp_result) == TRUE)
			{
			retry_wait_time = i * runif(n=1, min=retry_wait, max=2*retry_wait)
			Sys.sleep(retry_wait_time)
			}
		} # END for (i in 1:(attempts+1))
	
	return(warning(paste("download_BLAST_xml(): no results after ", attempts, 
		   " attempts; please try again later", sep = "")))
	} # END .tryParseResult <- function(url, attempts)





get_longest_seqrecs <- function(Hit_def)
	{
	ex='
	Hit_def = blastres$seqdf$Hit_def
	'
	
	# NA error trap
	if ( (length(Hit_def) == 1) && (is.na(Hit_def) == TRUE) )
		{
		longest_name_seqrec = NA
		return(longest_name_seqrec)
		}
	
	# Split strings
	list_of_splits = strsplit(Hit_def, split=">")
	list_of_splits
	
	num_identical_seqrecs = sapply(X=list_of_splits, FUN=length)
	num_identical_seqrecs
	
	# For each entry, get the LONGEST name between brackets, and the associated record
	# (this should be much more specific than the "MULTISPECIES" entries that come 
	#  of the "Hit_id" WP_ - type records
	
	longest_name_seqrec = sapply(X=list_of_splits, FUN=choose_a_string_w_longest_brackets)
	longest_name_seqrec
	return(longest_name_seqrec)
	}
	

choose_a_string_w_longest_brackets <- function(splits)
	{
	ex = '
	Hit_def = blastres$seqdf$Hit_def
	list_of_splits = strsplit(Hit_def, split=">")

	splits = list_of_splits[[1]]
	splits2 = splits

	splits = list_of_splits[[16]]
	splits2 = splits

	'
	
	length_brackets = rep(0, times=length(splits))
	brackets_texts = gdata::trim(extract_last_brackets(list_of_strings=splits, replace_spaces=FALSE))
	brackets_lengths = sapply(X=brackets_texts, FUN=nchar)
	
	# Choose the longest that is NOT multispecies
	multispecies_TF = grepl(pattern="multispecies", x=brackets_texts, ignore.case=TRUE)
	sum(multispecies_TF)
	
	if (sum(multispecies_TF) > 0)
		{
		brackets_texts = brackets_texts[multispecies_TF]
		brackets_lengths = brackets_lengths[multispecies_TF]
		splits2 = splits[multispecies_TF]
		} else {
		splits2 = splits
		}
	
	longest_TF = brackets_lengths == max(brackets_lengths)
	brackets_texts = brackets_texts[longest_TF]
	brackets_lengths = brackets_lengths[longest_TF]
	splits3 = splits2[longest_TF]
	splits3
	
	longest_name = splits3[1]
	longest_name
	return(longest_name)
	} # END choose_a_string_w_longest_brackets









DNAmult_matches <- function(aln)
	{
	defaults='
	aln = blastres$aln[[1]]
	'
	# Get the number of IDENTICAL, NUCLEOTIDE positions
	matches_by_nucleotide = consensusMatrix(aln, as.prob=TRUE)
	matches_by_nucleotide[1:18, 1:10]
	# Just take where the frequency of A, C, G, or T (rows 1-4) equals 1
	nummatches = colSums(matches_by_nucleotide[1:4,] == 1)
	match_TF = nummatches == 1
	num_identical = sum(match_TF)
	
	# Get the number of DISAGREEING, NUCLEOTIDE positions
	nummatches = colSums(matches_by_nucleotide[1:4,] == 0.5)
	match_TF = nummatches == 2
	num_disagree = sum(match_TF)
	
	# Number of resolved positions
	num_resolved = num_identical + num_disagree
	
	# Alignment length
	aln_len = dim(aln)[2]
	
	# percent identical
	pctIdent = round( 100 * (num_identical / num_resolved), digits=1)
	
	match_df = c(pctIdent, num_identical, num_disagree, num_resolved, aln_len)
	match_df
	return(match_df)
	}


make_aln_from_qseq_hseq <- function(qseq, hseq, type="AA")
	{
	# Make the pairwise alignment objects
	if (type == "AA")
		{
		aln = Biostrings::AAMultipleAlignment(
			c(qseq, hseq), 
			rowmask = as(IRanges(), "NormalIRanges"), 
			colmask = as(IRanges(), "NormalIRanges")
			)
		}
	if (type == "DNA")
		{
		aln = Biostrings::DNAMultipleAlignment(
			c(qseq, hseq), 
			rowmask = as(IRanges(), "NormalIRanges"), 
			colmask = as(IRanges(), "NormalIRanges")
			)
		}
	if (type == "RNA")
		{
		aln = Biostrings::RNAMultipleAlignment(
			c(qseq, hseq), 
			rowmask = as(IRanges(), "NormalIRanges"), 
			colmask = as(IRanges(), "NormalIRanges")
			)
		}

	return(aln)
	} # END make_aln_from_hseq_qseq <- function(qseq, hseq)


# Pretty damn slow!
DNAmult_matches_from_seqdf <- function(seqdf, type="AA")
	{
	defaults='
	seqdf = blastres$seqdf
	'
	
	names_match_dfs = c("pctIdent", "num_identical", "num_disagree", "num_resolved", "aln_len")
	
	# NA trp
	if ( (length(seqdf$Hit_def) == 1) && (is.na(seqdf$Hit_def) == TRUE) )
		{
		tmpdata = rep(NA, times=length(names_match_dfs))
		match_dfs = as.data.frame(matrix(data=tmpdata, nrow=1, ncol=length(names_match_dfs)), stringsAsFactors=FALSE)
		names(match_dfs) = names_match_dfs
		return(match_dfs)
		} # END if (is.na(seqdf$Hit_def) == TRUE)
	
	
	alns = mapply(FUN=make_aln_from_qseq_hseq, qseq=seqdf$qseq, hseq=seqdf$hseq, type=type)
	tmp_dfs = sapply(X=alns, FUN=DNAmult_matches)
	names(tmp_dfs) = NULL
	match_dfs = t(tmp_dfs)
	match_dfs = as.data.frame(match_dfs, row.names=FALSE, stringsAsFactors=FALSE)
	names(match_dfs) = names_match_dfs
	head(match_dfs)
	
	return(match_dfs)
	} # END DNAmult_matches_from_seqdf <- function(seqdf)


# Get the gi and the first gb from blastres$seqdf$Hit_id
get_first4_from_split_Hit_id <- function(split_Hit_id)
	{
	gi_gbs = split_Hit_id[1:4]
	return(gi_gbs)
	} # END get_first4_from_split_Hit_id <- function(split_Hit_id)
	
get_gi_gb <- function(Hit_ids)
	{
	defaults='
	Hit_ids = blastres$seqdf$Hit_id
	'
	split_Hit_ids = sapply(X=Hit_ids, FUN=strsplit, "\\|")
	gi_gbs = t(sapply(X=split_Hit_ids, FUN=get_first4_from_split_Hit_id))
	gi_gbs = as.data.frame(gi_gbs, row.names=FALSE, stringsAsFactors=FALSE)
	#names(gi_gbs) = c("giTF", "gi", "gbTF", "gb")
	names(gi_gbs) = c("gi1", "gi2", "gi3", "gi4")
	#gi_gbs$gi = suppressWarnings(as.numeric(gi_gbs$gi)) # NAs can result if <4 fields in names
	return(gi_gbs)
	} # END get_gi_gb <- function(Hit_ids)

get_sp_ssp_from_hit_id <- function(hit_id, genus)
	{
	require(gdata) # for trim
	
	# Get the start/end of the first hit
	startpos = regexpr(pattern=genus, text=hit_id, ignore.case=FALSE)
	endpos = startpos + nchar(genus) - 1
	str_length = nchar(hit_id)
	print(str_length)
	subtxt = substr(x=hit_id, start=endpos+1, stop=str_length)
	print(subtxt)
	end_string = gdata::trim(subtxt)
	
	# Error check - if genus is at the end, or genus not found
	if (((endpos+1) > str_length) || (startpos == -1))
		{
		genus = genus
		species = ""
		subspecies = ""
		gsp = matrix(c(genus, species, subspecies), nrow=1)
		gsp = as.data.frame(gsp, row.names=NULL, stringsAsFactors=FALSE)
		names(gsp) = c("genus", "species", "subspecies")
		return(gsp)
		} # END if ((endpos+1) > str_length)
	
	# Extract the species (and perhaps subspecies) after the genus
	words = strsplit(end_string, split=" ")[[1]]
	
	# Another error check
	if (length(words[[1]]) == 0)
		{
		genus = genus
		species = ""
		subspecies = ""
		gsp = matrix(c(genus, species, subspecies), nrow=1)
		gsp = as.data.frame(gsp, stringsAsFactors=FALSE, row.names=NULL)
		names(gsp) = c("genus", "species", "subspecies")
		return(gsp)
		} # END if (length(words[[1]]) == 0)
	
	if (length(words) >= 1)
		{
		species = words[[1]]
		
		# Check for equals
		if (grepl(pattern="=", x=species) == TRUE)
			{
			words2 = strsplit(species, split="=")[[1]]
			species = words2[1]
			}
		
		if (length(words) >= 2)
			{
			subspecies = words[[2]]
			} else {
			subspecies = ""
			} # END if (length(words) >= 2)
		} else {
		species = ""
		subspecies = ""
		} # END if (length(words) >= 2)

	gsp = matrix(c(genus, species, subspecies), nrow=1)
	gsp = as.data.frame(gsp, stringsAsFactors=FALSE, row.names=NULL)
	names(gsp) = c("genus", "species", "subspecies")
	return(gsp)
	} # END get_sp_ssp_from_hit_id


seqdf_enhance <- function(seqdf, genera=NULL, calc_match_percentages=TRUE, type="AA")
	{
	defaults='
	seqdf=blastres$seqdf
	genera = c("Canis", "Lycaon")
	calc_match_percentages=TRUE
	
	genera=NULL
	calc_match_percentages=TRUE
	type="AA"
	' # END defaults

	Hit_ids = seqdf$Hit_id
	if ( (length(Hit_ids) == 1) && (is.na(Hit_ids) == TRUE) )
		{
		gi_gbs = NA
		} else {
		gi_gbs = get_gi_gb(Hit_ids)
		}
	
	# Find mentions of "mitochond"
	mt1 = grepl(pattern="mitochond", x=seqdf$Hit_def)
	mt2 = grepl(pattern="mtDNA", x=seqdf$Hit_def)
	mt = (mt1 + mt2) > 0

	# Find mentions of "numt" (mtDNA transfered to the nucleus)
	numt = grepl(pattern="numt", x=seqdf$Hit_def, ignore.case=TRUE)

	# Find mentions of "chloroplast"
	cp1 = grepl(pattern="chloroplast", x=seqdf$Hit_def, ignore.case=TRUE)
	cp2 = grepl(pattern="cpDNA", x=seqdf$Hit_def, ignore.case=TRUE)
	cp3 = grepl(pattern="plastid", x=seqdf$Hit_def, ignore.case=TRUE)
	cp = (cp1 + cp2 + cp3) > 0

	# Find mentions of "complete genome"
	cgenome = grepl(pattern="complete genome", x=seqdf$Hit_def, ignore.case=TRUE)
	# Find mentions of "partial genome"
	pgenome = grepl(pattern="partial genome", x=seqdf$Hit_def, ignore.case=TRUE)
	pCDS = grepl(pattern="partial CDS", x=seqdf$Hit_def, ignore.case=TRUE)
	partial = grepl(pattern="partial", x=seqdf$Hit_def, ignore.case=TRUE)
	partial[pgenome==TRUE] = FALSE
	
	# Percentage matches
	pctIdent = round(100*seqdf$Hsp_identity/seqdf$Hsp_align_len, 1)
	pctSim = round(100*seqdf$Hsp_positive/seqdf$Hsp_align_len, 1)
	
	# Make the new fields into a data.frame
	tmpmat = cbind(pctIdent, pctSim, mt, numt, cp, cgenome, pgenome, pCDS, partial)
	typedf = as.data.frame(tmpmat, row.names=NULL, stringsAsFactors=FALSE)
	names(typedf) = c("pctIdent", "pctSim", "mt", "numt", "cp", "cgenome", "pgenome", "pCDS", "partial")
	
	
	# Calculate genus, species from matching genus names
	gsp_df = get_genera_from_seqdf(seqdf, genera=genera)


	# Calculate match percentages
	matches_df = matrix(data=NA, nrow=nrow(seqdf), ncol=5)
	matches_df = as.data.frame(matches_df, row.names=NULL, stringsAsFactors=FALSE)
	names(matches_df) = c("pctIdent", "num_identical", "num_disagree", "num_resolved", "aln_len")
	
	# Calculate the match percentages etc.
	# Slooooooooow....
	# Pretty damn slow!
	if (calc_match_percentages == TRUE)
		{
		txt = "\n\nCalculating pairwise match percentages...\n\n"
		matches_df = DNAmult_matches_from_seqdf(seqdf, type=type)
		} # END if (calc_match_percentages == TRUE)
	
	# Add the best-guess ID to download
	Hit_def = seqdf$Hit_def
	longest_name_seqrecs = get_longest_seqrecs(Hit_def)
	IDs_for_longest_seqrecs = gdata::trim(get_secondwords(longest_name_seqrecs, split="\\|"))
	
	# 2nd best: accession2 number
	seq_ids = seqdf$accession2
	IDs_for_longest_seqrecs[is.na(IDs_for_longest_seqrecs)] = seq_ids[is.na(IDs_for_longest_seqrecs)]
	
	# 3rd best: Hit_accession
	Hit_accession = seqdf$Hit_accession
	IDs_for_longest_seqrecs[is.na(IDs_for_longest_seqrecs)] = Hit_accession[is.na(IDs_for_longest_seqrecs)]
	downID = IDs_for_longest_seqrecs
	
	taxon = extract_last_brackets(list_of_strings=longest_name_seqrecs, replace_spaces=TRUE, prflag=FALSE)
	taxon
	
	# Combine
	row.names(seqdf) = NULL
	ncols_seqdf = ncol(seqdf)
	seqdf2 = cbind(seqdf[, 1], downID, taxon, gi_gbs, matches_df, gsp_df, typedf, seqdf[,2:ncols_seqdf])
	
	names(seqdf2)[1] = names(seqdf)[1]
	#seqdf2[1:5, 1:38]
	#cls.df(seqdf2)
	
	return(seqdf2)
	} # END seqdf_enhance <- function(seqdf)





count_blastresult <- function(result)
	{
	defaults='
	wd = "/drives/Dropbox/_njm/"
	setwd(wd)
	fasta_fn = "Wayne_etal_1997_SysBio_canid_mtDNA_appendix_wolf_cytB.fasta"
	mySeq = read.fasta(fasta_fn)
	x = convert_seqs_to_seqCollapse(x=mySeq[[1]])
	result = blastSeqKK2(x=x, database="nr", hitListSize="10", 
				filter="L", expect="10", program="blastn",
				attempts=10, organism=NULL, addl_url=NULL)
	' # END defaults
	
	# Extract information from the XML in "result"
	Hit_num = xpathApply(result, "//Hit_num", xmlValue)
	Hsp_num = xpathApply(result, "//Hsp_num", xmlValue)
	
	num_Hsps = length(Hsp_num)
	
	# Summarize
	txt = paste0("\n\nThe XML returned by BLAST has ", length(Hit_num), " hits and ", num_Hsps, " Hsps.\n\n")
	cat(txt)
	
	return(num_Hsps)
	} # END count_blastresult


get_genera_from_seqdf <- function(seqdf, genera=NULL)
	{
	# Pull out genus/species names
	gsp_df = matrix(data=NA, nrow=nrow(seqdf), ncol=5)
	gsp_df = as.data.frame(gsp_df, row.names=NULL, stringsAsFactors=FALSE)
	names(gsp_df) = c("genus", "species", "subspecies", "gnsp", "gnspssp")
	
	if (!is.null(genera))
		{
		for (g in 1:length(genera))
			{
			genus = genera[g]
			TF = grepl(pattern=genus, x=seqdf$Hit_def, ignore.case=FALSE)
			
			tmpdf = t(sapply(X=seqdf$Hit_def[TF], FUN=get_sp_ssp_from_hit_id, genus=genus))
			tmpdf = as.data.frame(tmpdf, stringsAsFactors=FALSE, row.names=NULL)
			for (colnum in 1:3)
				{
				gsp_df[TF,colnum] = unlist(tmpdf[1:nrow(tmpdf),colnum])
				} # END for (colnum in 1:3)
			} # END for (g in 1:length(genera))
		} # END if (!is.null(genera))
	
	# Add genus_species and genus_species_ssp
	# (i.e., gnsp and gnspssp)
	gnsp = rep(NA, times=nrow(gsp_df))
	notNA = is.na(gsp_df$genus) == FALSE
	gnsp[notNA] = paste(gsp_df$genus[notNA], gsp_df$species[notNA], sep="_")

	gnspssp = rep(NA, times=nrow(gsp_df))
	gnspssp[notNA] = paste(gsp_df$genus[notNA], gsp_df$species[notNA], gsp_df$subspecies[notNA], sep="_")
	
	gsp_df[,4] = gnsp
	gsp_df[,5] = gnspssp
	
	gsp_df = as.data.frame(gsp_df, row.names=NULL, stringsAsFactors=FALSE)
	names(gsp_df) = c("genus", "species", "subspecies", "gnsp", "gnspssp")
	
	# Convert to character, not factor
	gsp_df$genus = as.character(gsp_df$genus)
	gsp_df$species = as.character(gsp_df$species)
	gsp_df$subspecies = as.character(gsp_df$subspecies)
	gsp_df$gnsp = as.character(gsp_df$gnsp)
	gsp_df$gnspssp = as.character(gsp_df$gnspssp)
	
	return(gsp_df)
	} # END get_genera_from_seqdf <- function(seqdf, genera=NULL


process_blastresult <- function(tryParse_results, genera=NULL, type="AA")
	{
	defaults='
	wd = "/drives/Dropbox/_njm/"
	setwd(wd)
	fasta_fn = "Wayne_etal_1997_SysBio_canid_mtDNA_appendix_wolf_cytB.fasta"
	mySeq = read.fasta(fasta_fn)
	x = convert_seqs_to_seqCollapse(x=mySeq[[1]])
	result = blastSeqKK2(x=x, database="nr", hitListSize="10", 
				filter="L", expect="10", program="blastn",
				attempts=10, organism=NULL, addl_url=NULL)
				

	genera=NULL
	type="AA"
	' # END defaults
	
	result = tryParse_results$result
	RID_url = tryParse_results$RID_url 				# sends in the sequence, returns a BLAST
	search_url  = tryParse_results$search_url # searches, using an RID
		
	# Extract information from the XML in "result"
	Hit_num = xpathApply(result, "//Hit_num", xmlValue)
	Hit_id = xpathApply(result, "//Hit_id", xmlValue)
	Hit_def = xpathApply(result, "//Hit_def", xmlValue)
	Hit_accession = xpathApply(result, "//Hit_accession", xmlValue)
	Hit_len = xpathApply(result, "//Hit_len", xmlValue)
	Hsp_num = xpathApply(result, "//Hsp_num", xmlValue)
	Hsp_bit_score = xpathApply(result, "//Hsp_bit-score", xmlValue)
	Hsp_score = xpathApply(result, "//Hsp_score", xmlValue)
	Hsp_evalue = xpathApply(result, "//Hsp_evalue", xmlValue)
	Hsp_query_from = xpathApply(result, "//Hsp_query-from", xmlValue)
	Hsp_hit_to = xpathApply(result, "//Hsp_hit-to", xmlValue)
	Hsp_query_frame = xpathApply(result, "//Hsp_query-frame", xmlValue)
	Hsp_hit_frame = xpathApply(result, "//Hsp_hit-frame", xmlValue)
	Hsp_identity = xpathApply(result, "//Hsp_identity", xmlValue)
	Hsp_positive = xpathApply(result, "//Hsp_positive", xmlValue)
	Hsp_gaps = xpathApply(result, "//Hsp_gaps", xmlValue)
	Hsp_align_len = xpathApply(result, "//Hsp_align-len", xmlValue)
	qseq = xpathApply(result, "//Hsp_qseq", xmlValue)
	hseq = xpathApply(result, "//Hsp_hseq", xmlValue)
	
	

	# Summarize
	txt = paste0("The XML returned by BLAST has ", length(Hit_num), " hits and ", length(hseq), " Hsps.\n")
	cat(txt)
	
	require(Biostrings)
	aln = list()
	seqdf = NULL
	
	seqdf_names = c("Hit_num", "Hit_id", "Hit_def", "Hit_accession", "accession2","Hit_len", "Hsp_num", "Hsp_bit_score", "Hsp_score", "Hsp_evalue", "Hsp_query_from", "Hsp_hit_to", "Hsp_query_frame", "Hsp_hit_frame", "Hsp_identity", "Hsp_positive", "Hsp_gaps", "Hsp_align_len", "qseq", "hseq")
	
	# Make alignment for each sequence
	current_Hit_num = 1
	previous_Hsp_num = 0
	previous_hit_done = FALSE
	count = 1
	cat("Processing ", length(hseq), " Hsps...", sep="")
	
	if (length(qseq) == 0)
		{
		NAdata = rep(NA, times=length(seqdf_names))
		seqrow = matrix(data=NAdata, ncol=length(seqdf_names), nrow=1)
		seqdf = as.data.frame(seqrow, stringsAsFactors=FALSE)
		seqdf$Hit_num = 0
		seqdf$Hit_len = 0
		seqdf$Hsp_num = 0
		} else {
		for (i in seq_len(length(qseq)))
			{
			if (count < 100)
				{
				count = count + 1
				} else {
				cat(i, " ")
				count = 0
				} # END if (count == 100)
		
			cat(i, sep="")
			cat(",", sep="")
		
			current_Hsp_num = Hsp_num[[i]]
			if (previous_Hsp_num >= current_Hsp_num)
				{
				previous_hit_done = TRUE
				}
			if (previous_hit_done == TRUE)
				{
				current_Hit_num = current_Hit_num + 1
				previous_hit_done = FALSE
				}
			previous_Hsp_num = Hsp_num[[i]]
		
			# Make the pairwise alignment objects
			if (type == "AA")
				{
				aln[[i]] = Biostrings::AAMultipleAlignment(
					c(hseq[[i]], qseq[[i]]), 
					rowmask = as(IRanges(), "NormalIRanges"), 
					colmask = as(IRanges(), "NormalIRanges")
					)
				}
			if (type == "DNA")
				{
				aln[[i]] = Biostrings::DNAMultipleAlignment(
					c(hseq[[i]], qseq[[i]]), 
					rowmask = as(IRanges(), "NormalIRanges"), 
					colmask = as(IRanges(), "NormalIRanges")
					)
				}
			if (type == "RNA")
				{
				aln[[i]] = Biostrings::RNAMultipleAlignment(
					c(hseq[[i]], qseq[[i]]), 
					rowmask = as(IRanges(), "NormalIRanges"), 
					colmask = as(IRanges(), "NormalIRanges")
					)
				}
		
			# Extract the hit ids
			accession2 = get_secondwords(unlist(Hit_id), split="\\|")

			n = current_Hit_num
			tmprow = c(Hit_num[[n]], Hit_id[[n]], Hit_def[[n]], Hit_accession[[n]], accession2[n], Hit_len[[n]], Hsp_num[[i]], Hsp_bit_score[[i]], Hsp_score[[i]], Hsp_evalue[[i]], Hsp_query_from[[i]], Hsp_hit_to[[i]], Hsp_query_frame[[i]], Hsp_hit_frame[[i]], Hsp_identity[[i]], Hsp_positive[[i]], Hsp_gaps[[i]], Hsp_align_len[[i]], qseq[[i]], hseq[[i]])
			seqdf = rbind(seqdf, tmprow)
			} # END for (i in seq_len(length(qseq)))
		} # END if (length(qseq) == 0)
	
	seqdf = as.data.frame(seqdf, row.names=NULL, stringsAsFactors=FALSE)
	names(seqdf) = seqdf_names
	
	cat("done6.\n")
	
	# Make appropriate columns numeric
	# This is hella-slow
	#seqdf = dfnums_to_numeric(seqdf)
	seqdf$Hit_num = as.numeric(seqdf$Hit_num)
	seqdf$Hit_len = as.numeric(seqdf$Hit_len)
	seqdf$Hsp_num = as.numeric(seqdf$Hsp_num)
	seqdf$Hsp_bit_score = as.numeric(seqdf$Hsp_bit_score)
	seqdf$Hsp_score = as.numeric(seqdf$Hsp_score)
	seqdf$Hsp_evalue = as.numeric(seqdf$Hsp_evalue)
	seqdf$Hsp_query_from = as.numeric(seqdf$Hsp_query_from)
	seqdf$Hsp_hit_to = as.numeric(seqdf$Hsp_hit_to)
	seqdf$Hsp_query_frame = as.numeric(seqdf$Hsp_query_frame)
	seqdf$Hsp_hit_frame = as.numeric(seqdf$Hsp_hit_frame)
	seqdf$Hsp_identity = as.numeric(seqdf$Hsp_identity)
	seqdf$Hsp_positive = as.numeric(seqdf$Hsp_positive)
	seqdf$Hsp_gaps = as.numeric(seqdf$Hsp_gaps)
	seqdf$Hsp_align_len = as.numeric(seqdf$Hsp_align_len)
	
	cls.df(seqdf)
	seqdf[1:5, 1:17]

	
	seqdf2 = seqdf_enhance(seqdf=seqdf, genera=genera, type=type)
	dim(seqdf2)
	#seqdf2[1:5, 1:25]

	
	# Return BLAST results
	blastres = NULL
	blastres$aln = aln
	blastres$seqdf = seqdf2
	blastres$RID_url = RID_url
	blastres$search_url = search_url
	
	# Extract
	extract = '
	seqdf = blastres$seqdf
	'
	
	blastres
	return(blastres)
	} # END process_blastresult <- function(result)




list_of_XML_fns_fill_in_search_df <- function(xml_fns, fasta_blast_table_fn, type="AA", baseUrl="default")
	{
	ex='
	library(openxlsx)
	library(ape)
	library(phytools)
	library(BioGeoBEARS)
	library(varhandle) # for check.numeric
	library(rentrez) # Fetching full records: entrez_fetch()
	library(XML)
	library(Biostrings) # For parsing BLAST searches downloaded
	
	sourceall("/GitHub/str2phy/Rsrc/")
	source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")
	
	# Search on the 133aas_wo_UniProt.fasta file, s1 = search attempt #1
	wd = "/GitHub/str2phy/ex/MitoCOGs/_03_alphafolds/133aas_wo_UniProt_s1/"
	setwd(wd)

	
	xml_fns = list.files(path=".", pattern="*.xml")
	xml_fns
	
	type="AA"
	baseUrl="default"
	
	'
	# 133aas_wo_UniProt_blast_table_v1.txt
	# 133aas_wo_UniProt_blast_table_v1.txt
	search_df = read_search_df(fasta_blast_table_fn)
	i = 1
	for (i in 1:length(xmls_fns))
		{
		blastXML_fn = xml_fns[i]
		
		tmp_result2 = xmlTreeParse_blastXML_wErrorCatch(blastXML_fn=blastXML_fn, validate=FALSE, useInternalNodes=TRUE, isSchema=FALSE, error = xmlErrorCumulator(immediate=FALSE), isHTML=FALSE)
		
		if (is.null(tmp_result2) == TRUE)
			{
			next()
			}
		
		# Otherwise, 
		tryParse_results = NULL
		tryParse_results$result = tmp_result2
		tryParse_results$blast_download_summary_txt = "xmlTreeParse_blastXML_wErrorCatch() parsing of blastXML_fn"
		tryParse_results$blastXML_fn = blastXML_fn
		tryParse_results$total_pause_time = 0

		if (baseUrl == "uoa")
			{
			baseUrl = "https://blast-ncbi-nlm-nih-gov.ezproxy.auckland.ac.nz/Blast.cgi"
			}
		if (baseUrl == "default")
			{
			#baseUrl = "https://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
			baseUrl = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
			}

		
		# This URL puts up a sequence and returns an RID
		tryParse_results$RID_url = NA
	
		gid = get_gid_from_fn(fn=blastXML_fn)
		rid = get_RID_from_fn(fn=blastXML_fn)
		if (rid != "")
			{
			search_url = sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl, rid)
			} else {
			search_url = ""
			}
		
		#blastres = process_blastresult(tryParse_results, genera=NULL, type="AA")
		try_processing_result = try(process_blastresult(tryParse_results=tryParse_results))

		if (class(try_processing_result) == "try-error")
			{
			# Parsing issue
			next()
			} else {
			# If successful...
			blastres = try_processing_result
			
			
			search_df$pTF[rownum] = TRUE
			search_df$eTF[rownum] = FALSE
			
			# Number of hits once processed
			number_of_blast_hits = length(blastres$seqdf$Hit_num) - sum(is.na(blastres$seqdf$Hit_num))
			search_df$pTF_numhits[rownum] = number_of_blast_hits
			
			cat("...parsing of rownum=", rownum, " successful.\n", sep="")
			} # END if (class(try_result) == "try-error")

		
		} # END for (i in 1:length(xmls_fns))
	} # END list_of_XML_fns_fill_in_search_df <- function(xml_fns, fasta_blast_table_fn, type="AA")






# Split cells of a vector / spreadsheet on separation (eg comma)
split_cells_on_sep <- function(vec, split=", ")
	{
	veclist = lapply(X=vec, FUN=split_cell_on_sep, split=split)
	newvec = unlist(veclist)
	return(newvec)
	} # END split_cells_on_sep <- function(vec, sep=",")

split_cell_on_sep <- function(cell, split=",")
	{
	words = strsplit(cell, split=split)[[1]]
	return(words)
	} # END split_cells_on_sep <- function(vec, sep=",")



# Convert Excel column to a vector of items
excel_column_to_vector <- function(vec, ifcommas="takeall")
	{
	defaults='
	vec=xls$CytB_col1
	'
	
	# Check the column for cells with commas; if so, take first or last
	if (ifcommas != "takeall")
		{
		for (v in 1:length(vec))
			{
			TF = grepl(pattern=",", x=vec[v])
			if (TF == TRUE)
				{
				print(vec[v])
				words = strsplit(x=vec[v], split=",")[[1]]
				
				if (ifcommas == "first")
					{
					vec[v] = words[1]
					} # END if (ifcommas == "first")

				if (ifcommas == "last")
					{
					vec[v] = words[length(words)]
					} # END if (ifcommas == "last")
				} # END if (TF == TRUE)
			} # END for (v in 1:length(vec))
		} # END if (ifcommas != "takeall")
	
	
	# Remove blanks
	TF = isblank_TF(vec)==FALSE
	newvec = vec[TF]
	TF = newvec != ".--"
	newvec = newvec[TF]
	
	# Split on ", " or ","
	newvec = split_cells_on_sep(newvec, split=", ")
	newvec = split_cells_on_sep(newvec, split=",")
	
	return(newvec)
	}






#######################################################
# prflag
#######################################################
#' Utility function to conditionally print intermediate results
#'
#' Just a handy shortcut function, allowing other functions to optionally 
#' print, depending on the value of \code{printflag}.
#' 
#' @param x What to print.
#' @param printflag If TRUE, do the printing
#' @return nothing
#' @export
#' @seealso \code{\link{get_daughters}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
prflag <- function(x, printflag=TRUE)
	{
	# A standard function to print (or not) certain variables,
	#   based on a master printflag
	# This avoids having to comment in/out various code chunks
	#   while debugging.
	if (printflag == TRUE)
		{
		# CAT instead of PRINT if it's a string or numeric
		if (is.character(x))
			{
			cat(x, "\n", sep="")
			}
		if (is.numeric(x))
			{
			cat(x, "\n", sep="")
			} else {
			print(x)
			}
		}
	else
		{
		pass="BLAH"
		}
	}



get_numlines <- function(fn)
	{
	# Get the number of lines in a file
	
	# Check if the file exists
	TF = file.exists(fn)
	if (TF == FALSE)
		{
		txt = paste0("\nWARNING in get_numlines(): file '", fn, "' does not exist in directory\n'", getwd(), "'. Returning linecount=0.")
		cat(txt)
		linecount = 0
		return(linecount)
		}
	
	# intern = TRUE reports the result to R
	cmdstr = paste("wc -l ", fn, sep="")
	wc_return = system(cmdstr, intern=TRUE) 
	
	# split on whitespace (spaces and tabs)
	wc_list = strsplit_whitespace(wc_return)
	
	linecount = as.numeric(wc_list[1])
	return(linecount)
	} # END get_numlines <- function(fn)


#######################################################
# strsplit_whitespace
#######################################################
#' Split strings on whitespace
#' 
#' This function splits strings on whitespace (spaces and tabs), so you don't have
#' to remember the \code{regexp}/\code{grep} format codes.
#' 
#' @param tmpline A string containing text.
#' @return \code{list_of_strs} 
#' @export
#' @seealso \code{\link[base]{strsplit}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' tmpline = "Hello world see	my	tabs."
#' strsplit_whitespace(tmpline)
#' 
strsplit_whitespace <- function(tmpline)
	{
	# split on 1 or more whitespaces
	temp = strsplit(tmpline, "[ \t]+")
	
	# get the list
	list_of_strs = temp[[1]]
	
	# remove any leading/trailing ""
	list_of_strs = list_of_strs[list_of_strs != ""]
	
	return(list_of_strs)
	}


#######################################################
# moref
#######################################################
#' print to screen the header of a file
#' 
#' This does the rough equivalent of the \code{UNIX} function \code{more}, but within R.
#' 
#' @param fn A filename.
#' @param printnotcat If \code{TRUE}, use \code{\link[base]{print}} instead of \code{\link[base]{cat}}. Default \code{FALSE}.
#' @return Nothing returned.
#' @export
#' @seealso \code{\link[base]{scan}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
moref <- function(fn, printnotcat = FALSE)
	{
	lines = scan(file=fn, what="character", sep="\n")
	
	if (printnotcat == TRUE)
		{
		for (i in 1:length(lines))
			{
			print(lines[i])
			}
		}
	else
		{
		for (i in 1:length(lines))
			{
			cat(paste(lines[i], "\n", sep=""))
			}
		}
	}





# Convert blanks etc. to 0
convert_blanks_NAs_to_0 <- function(tipdates, newval=0)
	{
	defaults='
	newval=0
	'
	
	TF = tipdates == ""
	tipdates[TF] = newval

	TF = tipdates == " "
	tipdates[TF] = newval

	TF = tipdates == "\t"
	tipdates[TF] = newval

	TF = is.na(tipdates)
	tipdates[TF] = newval

	TF = is.nan(tipdates)
	tipdates[TF] = newval
	
	return(tipdates)
	}


isblank_TF <- function(items)
	{
	TF1 = items == ""
	TF2 = items == " "
	TF3 = items == "\t"
	TF4 = is.na(items)
	TF5 = is.nan(items)
	
	# Keep only the items where none of the above occur
	TFall = TF1 + TF2 + TF3 + TF4 + TF5
	
	# Correct for NA, NaNs
	TFall[is.na(TFall)] = 1
	TFall[is.nan(TFall)] = 1

	blank_TF = TFall > 0
	return(blank_TF)
	}


# Remove blanks etc. from a list
remove_blanks_NAs_etc <- function(items)
	{
	TF1 = items == ""
	TF2 = items == " "
	TF3 = items == "\t"
	TF4 = is.na(items)
	TF5 = is.nan(items)
	
	# Keep only the items where none of the above occur
	TF = TF1 + TF2 + TF3 + TF4 + TF5
	items = items[TF == 0]
	
	return(items)
	}




# Get RID from "_"-delimited filename (fn)
get_RID_from_fn <- function(fn)
	{
	ex = '
	fn = "2025-01-20_12.22.10_SUURJEGV013_RID_blastres_download.xml"
	rid = get_ID_from_fn(fn, word_before="RID")
	'
	
	rid = get_ID_from_fn(fn, split="_", word_before="RID")
	return(rid)
	}

# Get gid (genID or whatever) from "_"-delimited filename (fn)
get_gid_from_fn <- function(fn)
	{
	ex = '
	fn = "2025-01-20_12.22.10_SUURJEGV013_RID_blastres_download.xml"
	rid = get_gid_from_fn(fn, word_before="gid")
	'
	
	gid = get_ID_from_fn(fn, split="_", word_before="gid")
	return(gid)
	}


# Get gene id (or protein ID) from a sequence, delimited by split and identified by word_before
get_ID_from_fn <- function(fn, split="_", word_before="gid")
	{
	ex = '
	fn = "2025-01-20_12.22.10_RID_SUURJEGV013_gid_NP12345.1_blastres_download.xml"
	split="_"
	word_before = "gid"	
	get_ID_from_fn(fn, split=split, word_before=word_before)

	word_before = "rid"	
	get_ID_from_fn(fn, split=split, word_before=word_before)

	word_before = "RID"	
	get_ID_from_fn(fn, split=split, word_before=word_before)
	'
	
	words = strsplit(fn, split=split)[[1]]
	wordnums = 1:length(words)
	wordnum = match(x=word_before, table=words)
	if (is.na(wordnum) == TRUE)
		{
		id = ""
		} else {
		idnum = wordnum+1
		id = words[idnum]
		}
	return(id)
	}








xmlTreeParse_blastXML_wErrorCatch <- function(blastXML_fn, validate=FALSE, useInternalNodes=TRUE, isSchema=FALSE, error = xmlErrorCumulator(immediate=FALSE), isHTML=FALSE)
	{
	ex='
	library(openxlsx)
	library(ape)
	library(phytools)
	library(BioGeoBEARS)
	library(varhandle) # for check.numeric
	library(rentrez) # Fetching full records: entrez_fetch()
	library(XML)
	library(Biostrings) # For parsing BLAST searches downloaded
	
	sourceall("/GitHub/str2phy/Rsrc/")
	source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")
	
	# Search on the 133aas_wo_UniProt.fasta file, s1 = search attempt #1
	wd = "/GitHub/str2phy/Rsrc/blast_parsing_fns/"
	setwd(wd)
	
	validate=FALSE
	useInternalNodes=TRUE
	isSchema=FALSE
	error = xmlErrorCumulator(immediate=FALSE)
	isHTML=FALSE
	
	blastXML_fn = "BLAST_result_with_no_hits_2025-01-18_17.23.33_SP9SBRUP016_RID_blastres_download.xml"
	tmp_result = xmlTreeParse_blastXML_wErrorCatch(blastXML_fn=blastXML_fn, validate=validate, useInternalNodes=useInternalNodes, isSchema=isSchema, error=error, isHTML=isHTML)
	tmp_result

	blastXML_fn = "will_be_automatically_updated_2025-01-18_17.01.22_SP9SBRUP016_RID_blastres_download.xml"
	tmp_result = xmlTreeParse_blastXML_wErrorCatch(blastXML_fn=blastXML_fn, validate=validate, useInternalNodes=useInternalNodes, isSchema=isSchema, error=error, isHTML=isHTML)
	tmp_result

	blastXML_fn = "BLAST_result_with_hits_2025-01-18_10.31.52_SNJPB69S016_temp_blast_search_request.xml"
	tmp_result = xmlTreeParse_blastXML_wErrorCatch(blastXML_fn=blastXML_fn, validate=validate, useInternalNodes=useInternalNodes, isSchema=isSchema, error=error, isHTML=isHTML)
	tmp_result

	blastXML_fn = "RID_request_gives_html_parsing_error_2025-01-18_14.37.17_temp_blast_search_request.html"
	tmp_result = xmlTreeParse_blastXML_wErrorCatch(blastXML_fn=blastXML_fn, validate=validate, useInternalNodes=useInternalNodes, isSchema=isSchema, error=error, isHTML=isHTML)
	tmp_result

	'

	tmp_result = NULL
	tmp_result = tryCatch(
		{
		# Check if its a "This page will be automatically updated in" file
		checks = check_saved_BLAST_XMLs(blastXML_fn)
		still_running_automatic_updates_TF = checks$still_updating
		this_is_HTML_request_ID_TF = checks$this_is_HTML_request_ID

		if (this_is_HTML_request_ID_TF == TRUE)
			{
			post = htmlTreeParse(blastXML_fn, useInternalNodes = TRUE)
			x = post[["string(//comment()[contains(., \"QBlastInfoBegin\")])"]]
			rid = sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
			} else {
			rid = NULL
			}


		if ( (still_running_automatic_updates_TF == FALSE) && (this_is_HTML_request_ID_TF == FALSE) )
			{		
			# validate=FALSE was helpful for reading in downloaded BLAST XML results
			tmp_result2 = xmlTreeParse(file=blastXML_fn, validate=validate, useInternalNodes=useInternalNodes, isSchema=isSchema, error=error, isHTML=isHTML)
			return(tmp_result2)
			} else {
			if (this_is_HTML_request_ID_TF == TRUE)
				{
				warning(paste0("xmlTreeParse_blastXML_wErrorCatch() says: This is a search request HTML file for RID=", rid, ". Function can't tell if this search was ever launched."))
				return(NULL)
				}

			rid = get_RID_from_fn(fn=blastXML_fn)
			if (!is.null(rid))
				{
				warning(paste0("xmlTreeParse_blastXML_wErrorCatch() says: RID=", rid, " is still automatically updating."))
				return(NULL)
				} else {
				warning(paste0("xmlTreeParse_blastXML_wErrorCatch() cannot identify file, and no RID found."))
				return(NULL)					
				} # END if (!is.null(rid))

			warning("xmlTreeParse_blastXML_wErrorCatch(): should not get here.")
			return(NULL) # returns to tmp_result
			} # END if (still_running_automatic_updates_TF == FALSE)
		#result = htmlTreeParse(blastXML_fn, useInternalNodes=TRUE, error = xmlErrorCumulator(immediate=FALSE))
		}, # END of tryCatch, except for errors
		error = function(e){ message("An error occurred:\n", e) },
		warning = function(w){ message("A warning occured:\n", w) },
		finally = { message("...XML parsing finished.") } 
		#, error=function(err) NULL  # This last bit returns NULL
		) # END tmp_result = tryCatch
	return(tmp_result)
	} # END xmlTreeParse_blastXML_wErrorCatch <- function(blastXML_fn, validate=FALSE, useInternalNodes=TRUE, isSchema=FALSE, error = xmlErrorCumulator(immediate=FALSE), isHTML=FALSE)


	