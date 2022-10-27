### Download seleted genome from Mycocosm
Noted we will use the non-redundant proteins from each species, including protein isoforms created incorrect orthogroups. 

### Identify orthogoups
```bash
module load OrthoFinder/2.5.4-foss-2020b

orthofinder -f DecolorGenomes/ -t 16
```
This should finish within an hour for proteomes.

### CAFE
The lastest version of CAFE on HPC is CAFExp
```bash
module load CAFExp/5.0b2-intel-2019b
```
I followed the CAFE tutorial (https://github.com/hahnlab/cafe_tutorial) to prepare input files and run analysis. 

First, take the output file `Results_Oct25/Orthogroups/Orthogroups.GeneCount.tsv` and filter gene families. Gene families that have large gene copy number variance can cause parameter estimates to be non-informative.
```bash
# I reformat the gene family table into the way as cafe input and run filtering
python cafe_tutorial/python_scripts/cafetutorial_clade_and_size_filter.py -i Orthogroups.GeneCount.tsv -o filtered_
cafe_input.txt -s
# This creates two files, filtered_cafe_input.txt and large_filtered_cafe_input.txt
```

Next, we need to make ultrametric trees, which is phylogenetic trees scale to time. We need to calculate the number of alignment sites on the species in the analysis

I learned Exidia glandulosa belongs to order Auriculariales, and others are Polyporales. The divergence time between these two orders are 247.9497 Mya (refer to https://doi.org/10.1186/s12862-018-1229-7). We will use this number as a scale factor.

I used the script from Orthofinder tools, convert the rooted tree into ultrametric tree with the scale factor
```bash
cd CAFE4
module load OrthoFinder
wget https://raw.githubusercontent.com/davidemms/OrthoFinder/master/tools/make_ultrametric.py
python make_ultrametric.py -r 247.9497 ../Results_Oct25/Species_Tree/SpeciesTree_rooted.txt
# Ultrametric tree written to: ../Results_Oct25/Species_Tree/SpeciesTree_rooted.txt.ultrametric.tre
```
We will use this and the filtered input for CAFExp

### Run CAFExp
Rename the species ID in both tree file and input file, names should be consistant.
```bash
module load CAFExp/5.0b2-intel-2019b
cafexp -i filtered_cafe_input.txt -t SpeciesTree_rooted.txt.ultrametric.tre
#Best match is: 0.0013069055279629
#Final -lnL: 107734.95348471
#48 values were attempted (0% rejected)

#Inferring processes for Base model
#Score (-lnL): 107734.95348471
#Computing pvalues...done!
#Starting reconstruction processes for Base model
#Done!
```
A `results` directory is created!

### Inplement the results
Find out how many gene families is significant
```bash
# total gene families
grep -c 'OG' Base_family_results.txt
# 6956
# Count significant gene family
grep -c 'y' Base_family_results.txt
192
```
Let's focus on these 192 families!

### Run CAFE on large gene families
```bash
cafexp -i large_filtered_cafe_input.txt -t SpeciesTree_rooted.txt.ultrametric.tre -l 0.0013068948299531 -o large_results
```
There is only one orthogroup in this category, not enough for a result.

Run a gamma modeling
```bash
cafexp -i filtered_cafe_input.txt -t SpeciesTree_rooted.txt.ultrametric.tre -o gamma_results -k 2
# find 120 significant families
```
Try different k values, not sure which one gives the best result.

### Functional annotation for gene families
We will only analyze the significant gene families. Obtain Significant family list from gamma_results
```bash
grep 'y' Gamma_family_results.txt > Significant_families.txt
```
Take the list of OG IDs and save it in a file `Sig_OG.txt`. Grab the orthogroup sequence files from Orthofinder results
```bash
mkdir SignificantOG
while IFS= read -r line; do mv Orthogroup_Sequences/$line.fa SignificantOG/; done < Sig_OG.txt
# Check if the file number match significant OG number
ls SignificantOG/ | wc -l
```
Take one sequence from each OG and submit to InterproScan to search for the gene functions.
```bash
module load SeqKit/0.10.0-linux-x86_64
touch sigOG_sequences.fa
for file in SignificantOG/*fa
do
seqkit head -n 1 $file >> sigOG_sequences.fa
done
# chech the number of sequences extracted
grep -c '>' sigOG_sequences.fa
```
Use a python script to replace the headers to OG. Interproscan doesn't take '*" as a stop codon, need to remove it. Then submit a batch job for interproscan.

