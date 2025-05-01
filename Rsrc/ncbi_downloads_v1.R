require(BioGeoBEARS)

download_genome_w_esearch_wget <- function(genomeID)
	{
	examples ='
	library(ape)
	library(BioGeoBEARS)
	sourceall("/GitHub/str2phy/Rsrc/")
	source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")

	wd = "/Users/nmat471/Library/CloudStorage/Dropbox-UniofAuckland/flag/2025-03-06_genomes/tmp"
	setwd(wd)
	
	genomeID = "GCA_000353875.1"
	'
	esearch_cmd = paste0("esearch -query ", genomeID, " -db assembly | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank > ftp_url.txt")
	
	writeLines(esearch_cmd, con="esearch_cmd.sh")
	moref("esearch_cmd.sh")
	
	system("chmod a+x esearch_cmd.sh")
	
	sh_cmd = "./esearch_cmd.sh"
	system(sh_cmd)
	
	# Result
	ftp_url = readLines(con="ftp_url.txt")
	ftp_url
	
	# Genome name
	genome_name = get_lastword(ftp_url, split="/")
	genome_name
	
	#genomeID = get_leading_seqids_from_name(strings=genome_name, split="_", recodes=genbank_prefixes(), removetxt=c(">"))
	
	wget_cmd = paste0("wget -r –-no-host-directories --no-parent --no-directories --directory-prefix=", genomeID, " ", ftp_url)
	writeLines(wget_cmd, con="wget_ftp_cmd.sh")
	moref("wget_ftp_cmd.sh")
	
	system("chmod a+x wget_ftp_cmd.sh")

	sh_cmd = "./wget_ftp_cmd.sh"
	system(sh_cmd)

	file.copy(from="esearch_cmd.sh", to=genomeID, overwrite=TRUE)
	file.copy(from="ftp_url.txt", to=genomeID, overwrite=TRUE)
	file.copy(from="wget_ftp_cmd.sh", to=genomeID, overwrite=TRUE)
	file.remove("esearch_cmd.sh")
	file.remove("ftp_url.txt")
	file.remove("wget_ftp_cmd.sh")
	
	return(genome_name)
	}
	

get_fasta_seqs_from_seqids <- function(seqids, tmpfn="tmpseqs.fasta", db="protein", type="AA", more_than_500seqs=FALSE)
	{
	ex='
	seqids = c("AAC07752","AAC73831","AAC74960","AAC76042","AAG03587")
	tmpfn = "tmpseqs.fasta"
	db = "protein"
	type = "AA"
	more_than_500seqs = FALSE
	tmpseqs = get_fasta_seqs_from_seqids(seqids, tmpfn, db, type, more_than_500seqs)
	'
	
	if (more_than_500seqs == FALSE)
		{
		seqs_from_genbank = reutils::efetch(uid=seqids, db="protein", rettype="fasta", retmode="text")
		# Write to a temporary file
		write(reutils::content(seqs_from_genbank, "text"), file=tmpfn)
		}
		
	if (more_than_500seqs == TRUE)
		{
		# Allegedly, this is needed with more than 500 sequences
		seqs_from_genbank = reutils::efetch(uid=seqids, db="protein", rettype="fasta", retmode="text", outfile=tmpfn)
		# Write to a temporary file
		write(reutils::content(seqs_from_genbank, "text"), file=tmpfn)
		}
	
	# Collapse to strings and return
	seqs = sapply(X=as.character(read_FASTA_safe(file=tmpfn, type=type)), FUN=paste0, sep="", collapse="")
	return(seqs)
	}

