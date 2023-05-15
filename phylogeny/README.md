### Get data and database
ITS1 and ITS4 sequences obtained from chromatogram files and trimed with the script `trim_seq.py`.
```bash
python trim_seq.py <MTBxxx.TXT> <MTBxxxtrim.fasta>
# make sure Biopython is installed
```

Combine all the sequences:
```bash
cat *trim.fasta > MTB_all.fasta
```

Download fungal ITS datbase in fasta format from [UNITE](https://unite.ut.ee/). Create BLAST db:
```bash
module load BLAST+
makeblastdb -in sh_general_release_dynamic_s_16.10.2022.fasta -dbtype nucl -out UNITEdb
# error message:
# Error: (803.7) [makeblastdb] Blast-def-line-set.E.title
# Bad char [0xC3] in string at byte 23
# Aspicilia_simoënsis|EU057927|SH0055795.09FU|reps_singleton|k__Fungi;p__Ascomycota;c__Lecanoromycetes;o__Pertusariales;f__Megasporaceae;g__Aspicilia;s__Aspicilia_simoënsis
```
Characters like `ë` generate an error. We convert it to English 'e'. 
```bash
iconv -t ASCII//TRANSLIT sh_general_release_dynamic_s_16.10.2022.fasta > unite.seq.fasta
# run makeblastdb again
makeblastdb -in unite.seq.fasta -dbtype nucl -out UNITEdb
# now blastdb is created without error
```
### Run BLAST
```bash
blastn -query MTB_all.fasta -db UNITEdb -out ./BLAST.result.tsv -evalue 1e-10 -outfmt 6 -max_target_seqs 1
blastn -query MTB_all.fasta -db /scratch/data/bio/blastdb-2022.06.09/ITS_RefSeq_Fungi -out ./BLAST.Refseq.tsv -evalue 1e-10 -outfmt 6 -max_target_seqs 1 -num_threads 4
```

### Extract ITS1-2 sequences
```bash
module load SeqKit/0.10.0-linux-x86_64
seqkit grep -r -p '^(0[1-6]).*' MTB_all.fasta > ITS1-2.fasta
# sequences with headers 01-06 are ITS1-2, 07-12 are ITS4
```
### Run multiple sequence alignment
On terra, create a bash script
```bash
module load Clustal-Omega/1.2.4-intel-2018b

clustalo -i ITS/ITS1_species.fasta -o species-MSA.fa -v --outfmt=fa --threads=4
```
### Build maximum likelyhood tree with fasta file
```bash
input_dna_phy='species-MSA.fa'

######## PARAMETERS ########
model='GTRGAMMA'                    # see manual for many available options
algorithm='a'                       # see manual for many available options
p_seed=12345                        # random seed for parsimony inferences
num_runs=100                        # number of alternative runs on distinct starting trees

mpi_threads=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))   # number of raxmlHPC commands across all nodes
perhost=$SLURM_NTASKS_PER_NODE      # number of raxmlHPC commands per node
pthreads=$SLURM_CPUS_PER_TASK       # number of threads per raxmlHPC command

########## OUTPUTS #########
output_suffix="${model}_${algorithm}-${p_seed}-${num_runs}-hybrid-avx"

################################### COMMANDS ###################################

mpiexec -perhost $perhost -np $mpi_threads raxmlHPC -T $pthreads -m $model -p $p_seed -s $input_dna_phy -N $num_runs -n $output_suffix
```
