"""
Author: Jiali
This program take the chromatogram sequence file and extract the trimmed sequences.
Usage: python trim_seq.py <input file> <output fasta>
"""

import sys
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

file_in = sys.argv[1]
file_out = sys.argv[2]

def make_seq_record(record):
    return SeqRecord(seq=record.seq[40:620], \
            id=record.id, \
            description="trimmed")


sequences = SeqIO.parse(open(file_in), "fasta")
with open(file_out, "w") as out:
    for seq in sequences:
        #print(seq.description)
        if "CHROMAT" in seq.description:
            trim_Seq = make_seq_record(seq)
            SeqIO.write(trim_Seq, out, "fasta")
