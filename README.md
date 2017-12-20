# Polar Workshop Workflow
How you choose to create a dataset of orthologous gene alignments for your work is up to you. If you have questions about this step, I suggest that you look into the following methods. 

A) Transdecoder: https://github.com/TransDecoder/TransDecoder/wiki

B) OrthoFinder: https://github.com/davidemms/OrthoFinder

C) PhyloTreePruner:http://sourceforge.net/projects/phylotreepruner/
 

Your final alignments should only have each of your species represented once in each alignment (no paralogs).

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

Remember when editing or using perl and shell scripts for the first time you may need to change the permissions on the file to allow you to run it.

chmod 755 perl script or shell script should do the trick.

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

Below are some of the important results you should look at in your MCL_ALT_files if the likelihood ratio test is significant.

MLEs of dN/dS (w) for site classes (K=4)
site         class 0   1        2a         2b
proportion   0.74271   0.24369  0.01024    0.00336
background w 0.06038   1.00000  0.06038    1.00000
foreground w 0.06038   1.00000  176.87071  176.87071

M0 : Proportion of sites that are under purifying selection (ω0 < 1) on both foreground and background branches.
M1 : Proportion of sites that are under neutral evolution (ω1 = 1) on both foreground and background branches.
M2a: Proportion of sites that are under positive selection (ω2 ≥ 1) on the foreground branch and under purifying selection (ω0 < 1) on background branches.
M2b: Proportion of sites that are under positive selection (ω2 ≥ 1) on the foreground branch and under neutral evolution (ω1 = 1) on background branches.

Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005. Mol.
Biol. Evol. 22:1107-1118)
Positive sites for foreground lineages Prob(w>1):
33 T 0.975*
108 E 0.999**

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


5) GO-Panther Classifications for the Significant Genes.

All of you likely have a preferred way of gathering the GO and Panther classifications for your genes of interest, below I describe a relaitively simple and fast way to classify the molecular functions of your gene candidates.

First make a query file from the peptide sequences (one per candidate gene) in a fasta format.

Then go to Agbase http://agbase.msstate.edu/cgi-bin/tools/GOanna.cgi and load this file into the GOanna Tool. I suggest using the Agbase-Uniprot database but how you set up the stringency of your blastp search is up to you to get the most meaningful results. Depending on the taxa of choice you may need to adjust the default values for evalue, percent identity, and percent coverage.

The output from GOanna should be then be converted into a Gene Association Format using GOanna2ga

http://agbase.msstate.edu/cgi-bin/tools/GOanna2ga.cgi

The output from GOanna2ga provides Go classifications and Uniprot identifications for your genes. The Uniprto identification column will likely have repeat IDs from the blastp results, so it is best to creat a new text file list from this column (designated as "With (or) From"). 

Sort this text file to contain only the unique IDs

sort -u UniprotID.txt > UniprotID_s.txt

Then replace the new line breaks (/n) with commas

perl -p -e 's/\n/,/' UniprotID_s.txt > UniprotID_s2.txt

Use the comma deliminated Uniprot IDs to get Panther classifications for your genes.

http://pantherdb.org

Select the organism database(s) sets most applicable to your query and review the results for the five categories of Panther classifications: 1) Molecular Function, 2) Biological Processes, 3)Cellular Component, 4) Protein Class, and 5)Pathway.


OPTIONAL TESTS

1) SOWHAT program for automating the SOWH test (Church et al., 2015) to evaluate likelihood differences between habitat-constrained topologies for individual gene trees and phylogenomic-based species tree. Essentially the likelihood difference between the Best-Habitat-Constrained-Tree and the Best-Species-Tree provide a metric by which one may quantify the degree to which particular genes have characters that agree with Habitiat versus Species Relationship groupings. This test is based on RAxML and although RAxML can use either phylip or fasta alignments, it may be better to use phylip formats.

The sowhat program automates the SOWH phylogenetic topology test, which uses parametric bootstrapping to evaluate alternative tree topologies. It can be downloaded here (https://github.com/josephryan/SOWHAT). The process can be computationally intensive so I suggest running it on a subset of genes found significant by the other selection tests. A sample command line is shown below.

sowhat --constraint=Habitat_Tree.tre --aln=ALignment --name=Gene_Name --dir=Directory_Name --raxml_model=PROTGAMMAGTR --treetwo=Unrooted_Best_Species_Tree.tre



Select References

Benjamini Y. 2010. Discovering the false discovery rate. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 72:405–416.

Bielawski JP, Baker JL, Mingrone J. 2016. Inference of Episodic Changes in Natural Selection Acting on Protein Coding Sequences via CODEML. Curr Protoc Bioinformatics 54:6.15.1–6.15.32.

Church SH, Ryan JF, Dunn CW. 2015. Automation and Evaluation of the SOWH Test with SOWHAT. Syst Biol 64:1048–1058.

Delport W, Poon AFY, Frost SDW, Kosakovsky Pond SL. 2010. Datamonkey 2010: a suite of phylogenetic analysis tools for evolutionary biology. Bioinformatics 26:2455–2457.

Goldman, N. and Yang, Z. 1994. A codon-based model of nucleotide substitution for protein-coding DNA sequences. Mol. Biol. Evol. 11:725-736. 

Murrell B, Weaver S, Smith MD, Wertheim JO, Murrell S, Aylward A, Eren K, Pollner T, Martin DP, Smith DM, Scheffler K, Kosakovsky Pond SL. 2015. Gene-Wide Identification of Episodic Selection. Molecular Biology and Evolution 32:1365–1371.

Murrell B, Wertheim JO, Moola S, Weighill T, Scheffler K, Kosakovsky Pond SL. 2012. Detecting individual sites subject to episodic diversifying selection. PLoS Genet 8:1-10.

Muse, S.V. and Gaut, B.S. 1994. A likelihood approach for comparing synonymous and nonsynonymous nucleotide substitution rates, with application to the chloroplast genome. Mol. Biol. Evol. 11:715-724. 

Pond SLK, Frost SDW, Muse SV. 2005. HyPhy: hypothesis testing using phylogenies. Bioinformatics 21:676–679.
Revell LJ. 2011. phytools: an R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution 3:217–223.

Smith MD, Wertheim JO, Weaver S, Murrell B, Scheffler K, Kosakovsky Pond SL. 2015. Less is more: an adaptive branch-site random effects model for efficient detection of episodic diversifying selection. Molecular Biology and Evolution 32:1342–1353.

Thomas GWC, Hahn MW, Hahn Y. 2017. The effects of increasing the number of taxa on inferences of molecular convergence. Genome Biology and Evolution 9(1):213–221. 

Yang Z. 2007. PAML 4: Phylogenetic Analysis by Maximum Likelihood. Molecular Biology and Evolution 24:1586–1591.

Yang Z, Bielawski JP. 2000. Statistical methods for detecting molecular adaptation. Trends in Ecology & Evolution 15:496–503.

Yang Z, Wong WSW, Nielsen R. 2005. Bayes empirical bayes inference of amino acid sites under positive selection. Molecular Biology and Evolution 22:1107–1118.

Zhang J. 2005. Evaluation of an Improved Branch-Site Likelihood Method for Detecting Positive Selection at the Molecular Level. Molecular Biology and Evolution 22:2472–2479.





































