# AA_to_3Di
Matzke, Nicholas (2025). "AA to 3Di pipeline, version 1.0." <i>GitHub</i>, https://github.com/nmatzke/AA_to_3Di. Accessed April 29, 2025.

<b>Abstract:</b> This is a script archive giving a pipeline for converting amino acid sequences to protein structure "3Di" characters for each residue position. These scripts run in R, which then runs other command-line accessible programs such as ChimeraX and Foldseek. The scripts: 

<ul>
  <li>take a FASTA file of unaligned amino acid (AA) sequences</li>
  <li>download the closest-match AlphaFold structural predictions as .cif files, recording information on the percent identity and percent coverage of the match</li>
  <li>The scripts then convert each .cif file to 3Di characters using Foldseek</li>
  <li>assemble the AA and 3Di characters into two unaligned FASTA files where each sequence has AA and 3Di character data of identical length</li>
  <li>From there, users can align the AAs and/or 3Di characters using (non-R) programs such as famsa3di</li>
  <li>R functions are also given that can transfer an alignment structure from a 3Di alignment to an unaligned AA dataset, or vice versa, and/or</li>
  <li>horizontally concatenate an AA alignment and a 3Di alignment with matching sequence labels or matching sequence IDs in the labels</li>
</ul>

This pipeline has been submitted for the following publication:

Matzke, Nicholas J.; Puente-Lelievre, Caroline; Baker, Matthew A. B. (2025). "A Pipeline for Generating Datasets of 3-Dimensional Tertiary Interaction Characters for Model-Based Structural Phylogenetics." Chapter in: <i>Evolutionary Genomics: Methods and Protocols</i>. Series Title: <i>Methods in Molecular Biology</i>. Series Editor: Gustavo Caetano-Anollés. Submitted May 1, 2025.


An early version of this pipeline was used in:

Puente-Lelievre, Caroline; Malik, Ashar J.; Douglas, Jordan; Ascher, David; Baker, Matthew; Allison, Jane; Poole, Anthony; Lundin, Daniel; Fullmer, Matthew; Bouckert, Remco; Kim, Hyunbin; Obara, Masafumi; Steinegger, Martin; Matzke, Nicholas J. (2024) Tertiary-interaction characters enable fast, model-based structural phylogenetics beyond the twilight zone. <i>bioRxiv</i>:571181. doi: <a href="https://www.biorxiv.org/content/10.1101/2023.12.12.571181v2">10.1101/2023.12.12.571181</a>

Please cite these publications if you use this pipeline and/or the R functions contained therein.

