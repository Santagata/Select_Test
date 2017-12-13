# Polar Workshop Workflow
How you choose to create a dataset of orthologous gene alignments for your work is up to you. If you have questions about this step, I suggest that you look into the following methods. 

A) Transdecoder: https://github.com/TransDecoder/TransDecoder/wiki

B) OrthoFinder: https://github.com/davidemms/OrthoFinder

C) PhyloTreePruner:http://sourceforge.net/projects/phylotreepruner/
 

Your final alignments should onlu have each of your species represented once in each alignment (no paralogs).

Once you have your alignment peptide files there will likely be a need to trim the names of species using some form of sed command and/or perl script that keeps everyting up to the first space on the Species ID line (e.g., see Trimtospace.pl). You do not want to overwrite your original files. One possible solution is listed below.

for file in *.fa; do sed 's///g' $file > "$(basename "$file" .tex)_1"; done

for file in *.fa_1; do perl Trimtospace.pl $file > "$(basename "$file" .tex)_2"; done

Keep in mind that you cannot trim the names too much, as you will need to Species_Name_contig_ID_gene_ID to gather the correct cds sequences for each alignment file.

There are numerous ways to accomplish this using various unix/linux or perl scripts, below I have included one example.

Your peptide alignment files will likely have species_gene_IDs similar to the following.

>Camptaug|Gene_17689__TR19742_0_g1_i1__g_17689__m_17689
peptide sequence
>Ccrass|Gene_62889__TR4866_0_g1_i1__g_62889__m_62889
peptide sequence

Make a header list for each alignment using grep.

for file in *.fa_pruned.fa_2; do grep '>' $file > "$(basename "$file" .tex)_headlist.txt"; done

Using this header list you can use the provided perl script (subset_fasta.pl) to gather the necessary cds sequences for each alignment file. You will need to combine all the individual cds files for earch species into one file using the cat commnand.

cat *.cds.fa > ALLspecies_cds.fa

for file in *.fa_pruned.fa_2_headlist.txt; do perl subset_fasta.pl -i $file Allspecies_cds.fasta > "$(basename "$file" .tex)_cds.fa"; done

Now you will have a directory of alignment.peptide.fa files with the corresponding cds.fa files with the same alignment basenames.

OG0009035.pep.fa_2 OG0009035._cds.fa   #Remember the headers for each species must be identical within each peptide alignment and cds file.

You can now make two lists 1) List for the peptide alignments and 2) List for the corresponding cds files.

ls *.pep.fa_2 > ALN_list.txt
ls *._cds.fa > CDS_list.txt

You will need to install PAL2NAL. It can be downloaded here http://www.bork.embl.de/pal2nal/. Make sure the prgram is accessible from your account path.

Use the perl script CDS_AA_Hash_for_PAL2NAL.pl with the command usage. Be sure to understand the usage of tthe screen command (https://thesystemadministrator.net/cpanel/how-to-install-and-use-screen-in-linux) so that long processes do not stop when you exit your server.

perl CDS_AA_Hash_for_PAL2NAL.pl CDS_list.txt ALN_list.txt

the perl script is set up to output the CDS alignment file for PAML = Phyllip format with no gaps nor inframe stop codons in the alignment. You should run this to produce files for PAML and also modifiy the command line for pal2nal below to output fasta files for HYPHY.

print "perl pal2nal.pl $a_files[$i] $c_files[$i] -output paml -nogap -codontable 1 > $c_files[$i]_align.phy\n";

system "perl pal2nal.pl $a_files[$i] $c_files[$i] -output paml -nogap -codontable 1 > $c_files[$i]_align.phy\n";

change paml to fasta

Now you have all the aligned cds files that you will need for PAML:CODEML and HYPHY.

1)  Running CODEML through PAML

PAML may be downloaded here: http://abacus.gene.ucl.ac.uk/software/paml.html

We will be using the Branch-Site Model within codeml that allows for 2 or more foreground branches to be tested for positive (episodic) selection.

You will need to create a control file *.ctl for running CODEML using the null hypothesis where omega=1 on the foreground branches and also the alternative hypothesis where omega is allowed to be > 1 and estimated from the data provdided using maximum likelihood.

You will need to create a list of the aligned cds files you plan to test as described above. You should also creat a list that repeats the corresponding species tree file for as many aligned cds files that you wish to test. The perl script below will allow you to either use the same species tree for all your tests are try out different species trees. 

You will also need to label the foreground branches (Polar Species) on your species tree (using a newick format) with #1, for example see the tree file below. Your tree should include branch lengths and be UNROOTED. Particular nodes may also be labeled in the foreground. If you have trouble creating an unrooted tree consider using the included Afterphylo.pl perl script.

(((Bugner:0.159527,(Camptaug#1:0.077304,Ausvul#1:0.071029):0.030435):0.055321,(Cellarlat#1:0.116033,((WatersipFL:0.127171,Ccrass#1:0.070258):0.042494,SchizoFL:0.143287):0.020109):0.055403):0.0498,(Amvid:0.259907,(Horn#1:0.976799,(Alflab#1:0.28738,Alcyohauf:0.24686);

Do not forget to run the screen command so that you may detach from the screen when running long repetetive processes that require that you logout and wait for results. Screen can be downloaded here (https://thesystemadministrator.net/cpanel/how-to-install-and-use-screen-in-linux).

The perl script CODEML_TREE_ALN_NULL.pl will create the control file for running the null model.

The perl script CODEML_TREE_ALN_ALT.pl will create the control file for running the alternative model.

Again make a list of the CTL files that you wish to run, then use the perl script below.

>cat *.ctl > MCL_ALT_list.ctl
>cat *.ctl > MCL_NULL_list.ctl

The perl script CODEML_CTL_RUN.pl will run the control files for the null and alternative models. I suggest that you run all the ALT models first and all the NULL models second. You will need to modify the line below to include the complete path to your installation of CODEML.

system("/YOUR_PATH/codeml/paml4.9d/bin/codeml /YOUR_PATH/MCL_CTL_Files/$ctl_files[$i]\n");

>perl CODEML_CTL_RUN.pl MCL_ALT_list.ctl

>screen -d    #if you want to logout without stopping any processes.

Once complete you will have two CODEML MCL Results files for each of your CDS gene alignments (ALT versus NULL). Using grep you can now gather all of the lnL values for each model.

Example (degrees of freedom depends on the number of branches in your tree)
>grep lnL mlcnull_OG0009178.cdsaln.phy
>lnL(ntime: 21  np: 25):  -2592.532466      +0.000000

There are numerous ways to just grep the lnL from the null and alt models if you have them all in the same directory, but if you are unfamiliar with those options, just put them in separate directories.

grep lnL mlcnalt_OG000****.cdsaln.phy > ALT_list
grep lnL mlcnull_OG000****.cdsaln.phy > Null_list

These lists can be imported into excel or another spreadsheet program. I will provide a template for Excel.The difference between the ALT and NULL lnL can be used in a chi square test for level of significance.

Excel FIle Template = CODEML_Results
2(lnL1(ALT)-lnL0(NULL))=CHISQ.DIST(I10,1, FALSE)

The resulting P-values may be corrected for multiple comparisons using Benjamini and Hochberg (1995;2010) using R or simply use the website listed below.

https://www.sdmproject.com/utilities/?show=FDR

You will need to identify the resulting significant candidate genes under positive selection. If you are familiar with batch BLAST searches, use what you like but I suggest making a batch query file from a representative peptide sequence from each of your peptide alignment files. Again, there are many ways to do this, but one example is listed below.

Extract the first sequence from each peptide alignment file.

for file in *.pep.fa; do awk '/^>/{if(N)exit;++N;} {print;}' $file > "$(basename "$file" .tex)_seq1.fasta"; done

Then you can cat all of those sequences with their alignment IDs into one fasta file.

>cat *._seq1.fasta > candidategenes.fasta

blastp -query candidategenes.fasta -num_threads 8 -outfmt '6 qseqid sseqid pident length evalue bitscore staxids sscinames scomnames sskingdoms stitle' -max_target_seqs 1 -seg yes -evalue 0.001 -db nr -out Candidategenes_blastp 2> Candidategenes.err

NOTE: You may want to make your own blastp reference database of a subset of relevant sequences or use the entire blastp database if you have everything already installed on a local machine.
Obviously you can adjust many options in the above search, partiularly the number of maximum target sequences/evalue for genes with ambiguous identitities. 

makeblastdb -dbtype prot -in Relevant.fasta -out RelevantDBseqs.fasta


2) Testing your dataset using HyPHy: BUSTED

For all of the tests using HyPHy you will need to edit your tree file as to how foreground branches are designated. Instead of foreground branches designated with a #1 as with CODEML, you will need to use {Foreground}.

(((Bugner:0.159527,(Camptaug{Foreground}:0.077304,Ausvul{Foreground}:0.071029):0.030435):0.055321,(Cellarlat{Foreground}:0.116033,((WatersipFL:0.127171,Ccrass{Foreground}:0.070258):0.042494,SchizoFL:0.143287):0.020109):0.055403):0.0498,(Amvid:0.259907,(Horn{Foreground}:0.976799)

If you have trouble editing the newick tree file format, there is a helpful widget for selecting foreground braches and exporting the newick tree for Hyphy.

http://phylotree.hyphy.org

Remember all of your CDS alignment files must be in a fasta format made earlier using PAL2NAL.

You must install the latest version of HyPHy from here as some earlier versions have bugs in some of the selection tests we will use here.Latest version as of 12/13/17 is 2.3.7.

https://github.com/veg/hyphy/releases

Installation instructions can be found here. https://veg.github.io/hyphy-site/installation/

HYPHY installation depends on using cmake https://cmake.org/download/
https://cmake.org/install/

There are several ways to batch run files using Hyphy and parse the resulting *.json files from BUSTED, aBSREL, and MEME analyses. One interesting additional option if you prefer Biopython can be found at https://github.com/sjspielman/phyphy.

However I will suggest starting with a more simplistic approach using a shell script.

gen_busted.sh

Use a terminal window editor such as nano (or other) to edit the following line of the shell script.You need to designate the full path the to BUSTED.bf file depending on where you installed the program (example below is for usr/local). There are other versions of the *.bf files in the /TemplateBatchFiles/ directory but be sure to specify those analyses in the /TemplateBatchFiles/SelectionAnalyses/ directory.

fileToExe = "/usr/local/lib/hyphy/TemplateBatchFiles/SelectionAnalyses/BUSTED.bf";

To run BUSTED on one cds alignment file do the following.

./gen_busted.sh OG0008203.cdsaln.fa Speciestree_HYPHY.tre Foreground > OG0008203_script.bf

Then run that script using Hyphy

HYPHYMP OG0008203_script.bf

Results will be found in OG0008203.BUSTED.json. You can quickly assess whether the test was siginificant or not by using grep.

grep p-value OG0008203.BUSTED.json
>"p-value":1.269001392856239e-08

The results may also be parsed using Hyphy-Vision (http://vision.hyphy.org)

Obviously, you will want to set up and run all of these tests sequentially. To do this modify the above commands to the those listed below.

for file in *.fa; do ./gen_busted.sh $file SpeciesTree_Hyphy.tre Foreground > "$(basename "$file" .tex).BUSTED"; done

Then make a file list as described previously. 

cat *.BUSTED > Busted_list.txt

Then use the provided perl script to run all of the BUSTED tests on your files. You will need to edit the following line of the perl script below to match your full path to your *.BUSTED files that you just created.

system("HYPHYMP /full_path_to_your_Busted_files/BUSTED/$ctl_files[$i]\n");

>perl HYPHYMP_RUN.pl Busted_list.txt

>screen -d    #if you want to logout without stopping any processes.

3-4) Repeat this process for running HyPHy tests aBSREL and MEME using the following shell scripts

gen_aBSREL.sh

for file in *.fa; do ./gen_aBSREL.sh $file SpeciesTree_HYPHY.tre Foreground > "$(basename "$file" .tex).aBSREL"; done

gen_MEME2.sh

for file in *.fa; do ./gen_MEME2.sh $file SpeciesTree_HYPHY.tre Foreground > "$(basename "$file" .tex).MEME"; done

Parsing *.json files produced by aBSREL and MEME is not a simple one or two command process (see above for Biopython option). You may need to view all of these results using the Hyphy-Vision Option. 

When viewing MEME results click on the "p-value" column to bring all of the statistically significant sites to the top of the list.




































