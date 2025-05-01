#######################################################
# seqnames: functions for standardizing sequence names, removing unusual characters, etc.
# source("/GitHub/str2phy/Rsrc/seqnames_v1.R")
#######################################################

# Cat with \n
catn <- function(x="\n")
	{
	cat(x, sep="\n")
	}

# Cat to c("","")
catm <- function(x="")
	{
	txt = paste0(x, collapse='","')
	txt = paste0('c("', txt, '")')
	cat(txt)
	}

# Cat to c("","")
catc <- function(x="")
	{
	txt = paste0(x, collapse='","')
	txt = paste0('c("', txt, '")')
	cat(txt)
	}

# Remove blank lines, or any other double-whitespace, from e.g. a FASTA file
# https://stackoverflow.com/questions/41540972/delete-empty-lines-from-a-text-file-via-bash-including-empty-spaces-characters
remove_blanklines <- function(fn, newfn=NULL, prflag=TRUE)
	{
	ex='
	wd = "/GitHub/str2phy/ex/MitoCOGs/blast_downloads_v1/"
	setwd(wd)
	fn = "17_downloaded_seqs.fasta"
	newfn=NULL
	prflag=TRUE
	'
	
	write_over = FALSE
	if (is.null(newfn))
		{
		write_over = TRUE
		newfn = "tmp.txt"
		}
	cmd = paste0("perl -ne 'print if !/^[[:blank:]]*$/' ", fn, " > ", newfn)
	system(cmd)
	
	if (write_over == TRUE)
		{
		file.rename(from=newfn, to=fn)
		}
	
	txt = paste0("\nremove_blanklines(): '", fn, "' -> '", newfn, "'.")
	cat(txt)
	
	return(newfn)
	}



# Needs work
# Check if a pattern starts or ends with bar (|")
starts_or_ends_with_pipe_TF <- function(pattern)
	{
	# Check for "|" alone
	if (pattern == "|")
		{
		return(TRUE)
		}
	
	TF1 = startsWith(pattern, prefix="|")
	TF2 = endsWith(pattern, suffix="|")
	TF = (TF1 + TF2) == 1
	TF
	return(TF)
	}

#stop("STOP ERROR Finish grepl_safe()")

# Needs work
grepl_safe <- function(pattern, x, ignore.case=FALSE, perl=FALSE, fixed=FALSE, useBytes=FALSE)
	{
	ex='
	# Problem: grepl will return TRUE if the search PATTERN
	# argument starts, or ends, with the pipe "|" (which means "or")
	
	# Check pattern for "|"
	pattern = "lcl|"
	tmptxt = "abcdef"
	
	grepl(pattern="|", x=tmptxt)
	startsWith(pattern, prefix="|")
	endsWith(pattern, suffix="|")

	grepl(pattern="|", x=tmptxt)
	# TRUE (incorrectly)
	
	grepl(pattern="asdf|", x=tmptxt)
	# TRUE (incorrectly)

	grepl(pattern="|asdf", x=tmptxt)
	# TRUE (incorrectly)

	grepl(pattern="as|df", x=tmptxt)
	# FALSE (correctly)
	
	grepl_result = grepl_safe(pattern="|", x=tmptxt, ignore.case=FALSE, perl=FALSE, fixed=FALSE, useBytes=FALSE)
	grepl_result
	# FALSE (correctly)

	grepl_result = grepl_safe(pattern="asdf|", x=tmptxt, ignore.case=FALSE, perl=FALSE, fixed=FALSE, useBytes=FALSE)
	grepl_result
	# FALSE (correctly)

	grepl_result = grepl_safe(pattern="|asdf", x=tmptxt, ignore.case=FALSE, perl=FALSE, fixed=FALSE, useBytes=FALSE)
	grepl_result
	# FALSE (correctly)


	
	' # END example
	
	
	
	# Check for starting or ending with "|"
	if ( starts_or_ends_with_pipe_TF(pattern) == TRUE )
		{
		newpattern = gsub(pattern="\\|", replacement="_P_I_P_E_", x=pattern)
		newx = gsub(pattern="\\|", replacement="_P_I_P_E_", x=x)
		} else {
		newpattern = pattern
		newx = x
		}
	
	grepl_result = grepl(newpattern, newx, ignore.case=ignore.case, perl=perl, fixed=fixed, useBytes=useBytes)
	
	return(grepl_result)
	} # END grepl_safe


# Needs work
# Remove any word that is a grepl hit; words defined by split
remove_word <- function(tmptxt, pattern, split=" ", ignore.case=TRUE)
	{
	ex='
	tmptxt = "lcl|ORF31 32515-32994 rpS7"
	pattern = "lcl|"
	split=" "
	ignore.case=TRUE
	newtxt = remove_word(tmptxt, pattern, split=split, ignore.case=ignore.case)
	newtxt
	newtxts = remove_word_sapply(tmptxt, pattern, split=split, ignore.case=ignore.case)
	newtxts
	
	' # END example
	

	
	words = strsplit(x=tmptxt, split=split)[[1]]
	TF = grepl_safe(pattern=pattern, x=words, ignore.case=TRUE)
	newwords = words[TF==FALSE]
	newtxt = paste(newwords, sep="", collapse=split)
	newtxt
	return(newtxt)
	}


# Needs work
remove_word_sapply <- function(tmptxts, pattern, split=" ", ignore.case=TRUE)
	{
	newtxts = sapply(X=tmptxts, FUN=remove_word, pattern=pattern, split=split, ignore.case=ignore.case)
	return(newtxts)
	}





#######################################################
# Take a fasta file, reduce sequence names to just gids
# (assumes these are first, and take the usual form of
#  WP_, NC_, etc...)
#######################################################
fasta_w_justgids <- function(fn, split="_")
	{
	ex='
	fn = "/GitHub/str2phy/substmat/try_mafft_custom3DI/usalign_test3.fasta"
	seqs = ape::read.FASTA(fn)
	split="_"
	
	outfn = fasta_w_justgids(fn, split="_")
	moref(outfn)
	'
	
	seqs = ape::read.FASTA(fn)
	tipnames = names(seqs)
	tipnames3 = names_w_justgids(tipnames, split="_")
	names(seqs) = tipnames3
	
	# Split on lastword
	fnwords = strsplit(gdata::trim(fn), split="\\.")[[1]]
	final_word = fnwords[length(fnwords)]
	
	outfn = gsub(pattern=paste0(".", final_word), replacement="", x=fn)
	outfn = paste0(outfn, "_justgids.", final_word)
	
	ape::write.FASTA(x=seqs, file=outfn)
	return(outfn)
	}

names_w_justgids <- function(tipnames, split="_")
	{
	ex='
	fn = "/GitHub/str2phy/substmat/try_mafft_custom3DI/usalign_test3.fasta"
	seqs = ape::read.FASTA(fn)
	tipnames = names(seqs)
	split="_"
	names_w_justgids(tipnames, split="_")
	'
	
	# Substitute things that will fail
	oldstr = c("WP_", "NC_")
	newstr = c("WP-", "NC-")
	
	tipnames2 = tipnames
	for (i in 1:length(oldstr))
		{
		tipnames2 = gsub(pattern=oldstr[i], replacement=newstr[i], x=tipnames2)
		}
	
	# Remove leading ">", if any
	tipnames2 = gsub(pattern=">", replacement="", x=tipnames2)

	# Split names
	for (i in 1:length(tipnames2))
		{
		words = strsplit(tipnames2[i], split=split)[[1]]
		tipnames2[i] = words[1]
		}
	
	# Replace substitutions
	tipnames3 = tipnames2
	for (i in 1:length(oldstr))
		{
		tipnames3 = gsub(pattern=newstr[i], replacement=oldstr[i], x=tipnames3)
		}
	
	return(tipnames3)
	} # END names_w_justgids <- function(tipnames, split="_")


fixnames <- function(seqnames)
	{
	seqnames = gsub(pattern=" ", replacement="_", x=seqnames)
	seqnames = gsub(pattern="\\,", replacement="_", x=seqnames)
	seqnames = gsub(pattern="\\;", replacement="_", x=seqnames)
	seqnames = gsub(pattern="\\[", replacement="_", x=seqnames)
	seqnames = gsub(pattern="\\]", replacement="_", x=seqnames)
	seqnames = gsub(pattern="\\(", replacement="_", x=seqnames)
	seqnames = gsub(pattern="\\)", replacement="_", x=seqnames)
	seqnames = gsub(pattern="\\/", replacement="_", x=seqnames)
	seqnames = gsub(pattern="\\\\)", replacement="_", x=seqnames)

	seqnames = gsub(pattern="\\._", replacement="_", x=seqnames)
	seqnames = gsub(pattern="__", replacement="_", x=seqnames)
	seqnames = gsub(pattern="__", replacement="_", x=seqnames)
	seqnames = gsub(pattern="__", replacement="_", x=seqnames)
	seqnames = gsub(pattern="__", replacement="_", x=seqnames)
	seqnames = gsub(pattern="__", replacement="_", x=seqnames)
	seqnames = gsub(pattern="__", replacement="_", x=seqnames)
	seqnames = gsub(pattern="__", replacement="_", x=seqnames)
	seqnames = gsub(pattern="__", replacement="_", x=seqnames)
	seqnames = gsub(pattern="__", replacement="_", x=seqnames)
	
	seqnames = remove_trailing_underscores(seqnames)
	
	return(seqnames)
	} # END fixnames <- function(seqnames)


remove_trailing_underscores <- function(seqnames, endchar="_")
	{
	# Remove ending "_"
	TF = base::endsWith(x=seqnames, suffix=endchar)
	if (sum(TF) > 0)
		{
		nums = (1:length(TF))[TF]
		
		for (num in nums)
			{
			numchars = nchar(seqnames[num])
			seqnames[num] = base::substr(x=seqnames[num], start=1, stop=numchars-1)
			}
		} # END if (sum(TF) > 0)
	
	return(seqnames)
	} # END remove_trailing_underscores <- function(seqnames)


# Change the tipnames on a phylogenetic tree based on a unique ID
fix_tipnames_gid <- function(gidskey, newmatch, tr)
	{
	tipnames = tr$tip.label
	
	for (i in 1:length(gidskey))
		{
		TF = grepl(pattern=gidskey[i], x=tipnames)
		tipnum = ((1:length(TF))[TF])
		
		if (length(tipnum) > 0)
			{
			tipnames[tipnum[1]] = newmatch[i]
			} # END if (length(tipnum) > 0)
		} # END for (i in 1:length(gid))
	
	tr$tip.label = tipnames
	return(tr)
	} # END fix_tipnames_gid(gid, newnames, tr)



# Fix duplicates in a list of names
fix_duplicate_names <- function(sn)
	{
	duplicate_seqnames_table = rev(sort(table(sn)))
	if (max(duplicate_seqnames_table) < 2)
		{
		return(sn) # exit unchanged
		}
	duplicate_seqnames_table = duplicate_seqnames_table[duplicate_seqnames_table > 1]
	duplicate_seqnames_table
	duplicate_seqnames = names(duplicate_seqnames_table)	
	
	for (i in 1:length(duplicate_seqnames))
		{
		duplicate_seqname = duplicate_seqnames[i]
		TF = sn == duplicate_seqname
		nums = (1:length(TF))[TF]
		nums
		
		print(nums)
		
		# Add _dup2 etc on any repeats
		addnum = 0
		for (j in nums[1]:nums[length(nums)])
			{
			addnum = addnum+1
			sn[j] = paste0(sn[j], "_part0", addnum)
			}
		}
	return(sn)
	}



# Like match, but for grepl
# returns first grepl match in x
match_grepl <- function(patterns, x, return_counts=TRUE)
	{
	nums = 1:length(x)
	matchnums = rep(NA, times=length(patterns))
	counts = rep(0, times=length(patterns))
	for (i in 1:length(patterns))
		{
		try_result = try(grepl(pattern=patterns[i], x=x))
		if (class(try_result) == "try-error")
			{
			txt = paste0("\nSTOP ERROR in match_grepl(): in input 'patterns', item i=", i, " gave an error. patterns[i]=\n")
			cat(txt)
			print(patterns[i])
			txt2 = paste0("\nThe try_result was:\n")
			cat(txt2)
			print(try_result)
			cat("\n")
			stop(txt)
			}
		
		# Check passed, so try_result is TRUE/FALSE, i.e. "TF"
		TF = try_result
		matchnums_all = nums[TF]
		num_matches = sum(TF, na.rm=TRUE)
		matchnum = matchnums_all[1] # first match
		matchnums[i] = matchnum
		counts[i] = num_matches
		}
		
	if (return_counts == FALSE)
		{
		return(matchnums)
		} else {
		match_df = cbind(matchnums, counts)
		match_df = as.data.frame(match_df, stringsAsFactors=FALSE)
		return(match_df)
		}
	
	txt = "STOP ERROR in match_grepl(): Shouldn't get to end."
	stop(txt)
	}




# Assumptions:
# 1. The merge_names specify unique final sequences
# 2. If exons are on separate lines, they are in sequence order (exon 1 comes before exon 2)
#    - manually check and re-order in AliView, if this is not the case!
# 3. Merged sequences are non-overlapping (earlier ones will be overwritten, if untrue)
# 4. All rows must have the same length (inluding "-" characters)

merge_seqs_by_name_grepl <- function(aln, merge_names, type="DNA")
{
ex='
merge_names = c("Human",
"Neanderthal",
"Gorilla",
"Bonobo",
"Chimpanzee",
"Sumatran Orangutan",
"Gibbon",
"Macaque_sp",
"Rhesus macaque",
"Crab eating macaque",
"Drill",
"Olive baboon",
"Gelada",
"Black snub nosed monkey",
"Ugandan red colobus",
"Angola Colobus",
"Golden snub-nosed monkey",
"Mas Night Monkey",
"Green Monkey",
"Sooty Mangabey",
"Pig-Tailed Monkey",
"Tarsier",
"White-tufted ear marmoset",
"Bolivian Squirrel Monkey",
"Capuchin",
"Grey_Mouse_Lemur",
"Coquerels_Sifaka",
"Galago",
"Rabbit",
"Pika Ochotona_princeps",
"Rat",
"Mouse",
"Golden_Hamster",
"Mongolian_Gerbil",
"Striped_Gopher",
"Chinchilla",
"Brazilian_Guinea_Pig",
"Domestic_Guinea_Pig",
"Red_Fox",
"Dog",
"Californian_Sea_Lion",
"Cat",
"Canadian_Lynx",
"Great_Roundleaf_Bat",
"Greater_Horseshore_Bat",
"Leschenaults_Rousette",
"Horse",
"Donkey",
"Pig",
"Goat",
"Camel Camelus_dromedarius",
"Cow",
"Sheep",
"Narwhal",
"Beluga_Whale",
"Sperm_Whale",
"Vaquita Phocoena_sinus",
"Armadillo Dasypus_novemcinctus",
"Possum")

library(ape)
library(gdata)
source("/GitHub/str2phy/Rsrc/seqnames_v1.R")

wd = "/GitHub/gulo/data/"
setwd(wd)

# Load the FASTA alignment
#  - originally from Alexander Mansueto
#  - Nick manually added & aligned all other Mansueto sequences
#  - Here, we have to merge the exons (should have maybe done before, but
#      nice to have them on separate lines at first)
fasta_fn = "GULO_Total_merge_v7.fasta"

# Load the alignment
aln = ape::read.FASTA(fasta_fn, type="DNA")
class(aln)
aln = as.character(aln)
class(aln)

i = 1
type = "DNA"

new_aln2 = merge_seqs_by_name_grepl(aln, merge_names, type="DNA")
' # END ex

# Setup
seqnames = names(aln)
remaining_seqnames = names(aln)
length_aln = length(aln[[1]])
new_aln = NULL
new_aln_pos = 0
new_aln_names = NULL

# Maximum non-"-" position
lengths = sapply(X=aln, FUN=length)
minlength = min(lengths)
max_seq_pos = lengths
for (i in 1:length(aln))
	{
	nodashTF = aln[[i]] != "-"
	nodash_nums = (1:length(nodashTF))[nodashTF]
	max_seq_pos[i] = max(nodash_nums)
	}

# The length for all output sequences will be the 
# maximum of max_seq_pos
max_max_seq_pos = max(max_seq_pos)

for (i in 1:length(merge_names))
	{
	# Remaining names to search against

	nm = merge_names[i]
	TF = grepl(pattern=nm, x=seqnames, ignore.case=TRUE)
	matchnums = (1:length(TF))[TF]
	
	if (length(matchnums) == 0)
		{
		txt = paste0("STOP ERROR in merge_seqs_by_name_grepl(): merge_names[i=", i, "] = '", merge_names[i], "' was not found in seqnames. Edit merge_names to fix, and re-run.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (length(matchnums) == 0)

	if (length(matchnums) == 1)
		{
		new_aln_pos = new_aln_pos + 1
		new_aln[[new_aln_pos]] = aln[[matchnums]]
		
		# Cut to max length
		new_aln[[new_aln_pos]] = new_aln[[new_aln_pos]][1:max_max_seq_pos]
		
		new_aln_names[[new_aln_pos]] = nm
		
		paste0(new_aln[[1]], collapse="", sep="")
		
		# Edit remaining names
		keepTF = (remaining_seqnames %in% seqnames[matchnums]) == FALSE
		remaining_seqnames = remaining_seqnames[keepTF]
		} # END if (length(matchnums) == 1)

	if (length(matchnums) > 1)
		{
		new_aln_pos = new_aln_pos + 1

		newseq = rep("-", times=length_aln)
		for (j in 1:length(matchnums))
			{
			seq_to_stick_into_merge = aln[[matchnums[j]]]
			nodashTF = seq_to_stick_into_merge != "-"
			
			# Stick those characters in (newer seqs will overwrite previous insertions
			newseq[nodashTF] = seq_to_stick_into_merge[nodashTF]

			# Edit remaining names
			keepTF = (remaining_seqnames %in% seqnames[matchnums[j]]) == FALSE
			remaining_seqnames = remaining_seqnames[keepTF]
			} # END building newseq
		
		new_aln[[new_aln_pos]] = newseq
		# Cut to max length
		new_aln[[new_aln_pos]] = new_aln[[new_aln_pos]][1:max_max_seq_pos]

		new_aln_names[[new_aln_pos]] = nm
		} # END if (length(matchnums) > 1)
	} # END for (i in 1:length(merge_names))

# Finish building alignment object
names(new_aln) = unlist(new_aln_names)
new_aln
if (type == "DNA")
	{
	new_aln2 = ape::as.DNAbin(new_aln)
	}
if (type == "AA")
	{
	new_aln2 = ape::as.AAbin(new_aln)
	}

if (length(remaining_seqnames) > 0)
	{
	txt = "Returning merged alignment as new_aln2.\nNote: 'remaining_seqnames' were not merged, they are printed below:\n"
	cat("\n\n")
	cat(txt)
	catn(remaining_seqnames)
	cat("\n\n")
	} 
	
return(new_aln2)
} # END 




full_seqnames_to_just_label <- function(list_of_strings, split=" ")
	{
	ex= '
	list_of_strings = c("AAC07752.1 TolQ-like protein [Aquifex aeolicus VF5]","AAC73831.1 Tol-Pal system protein TolQ [Escherichia coli str. K-12 substr. MG1655]","AAC74960.1 motility protein A [Escherichia coli str. K-12 substr. MG1655]","AAC76042.1 Ton complex subunit ExbB [Escherichia coli str. K-12 substr. MG1655]","AAG03587.1 transport protein ExbB [Pseudomonas aeruginosa PAO1]")
	split = " "
	'
	
	seqnames_without_label = all_but_prefixes(tmptxts=list_of_strings, split=split)
	species_names = extract_last_brackets(list_of_strings=seqnames_without_label, replace_spaces=FALSE, prflag=FALSE)
	# Remove the species names
	species_names_w_brackets = paste0("\\[", species_names, "\\]")
	species_names_w_brackets
	seqnames_without_label_wo_species = mapply(FUN=gsub, pattern=species_names_w_brackets, x=seqnames_without_label, MoreArgs=list(replacement="")) 
	seqnames_without_label_wo_species = gdata::trim(seqnames_without_label_wo_species)
	seqnames_without_label_wo_species = unname(seqnames_without_label_wo_species)
	return(seqnames_without_label_wo_species)
	}



