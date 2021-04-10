#!/usr/bin/env python

# Author: Jacopo Martelossi

# This script will convert all input fasta files in nexus format, then it will concatenate them. 
#Output are: 
    #Concat file in fasta, nexus and phylip.
    #Partition file in nexus format.

# Named variables. Every run needs the following defined:
# 1) --in_dir - The directory containing the .fasta alignments that need to be merged.
# 2) --out - The full filepath and name you want for the output file.

import os
from Bio import SeqIO
from Bio.Nexus import Nexus
import argparse

# Argument Parser
parser = argparse.ArgumentParser(description = 'This script will convert all input fasta files in nexus format, then it will concatenate them in a supermatrix. Files must have the extension ".fasta"')
parser.add_argument('--in_dir', required=True, help='The input directory containing alignment files.')
parser.add_argument('--out', required=True, help='The filepath and filename of the output file.')
args = parser.parse_args() 

file_list = []

IN_DIR = args.in_dir
OUT = args.out
for filename in os.listdir(IN_DIR):
    if filename.endswith(".fasta"):
        PathFile=os.path.join(IN_DIR, filename)
        NewName=PathFile.replace(".fasta",".nex")
        count = SeqIO.convert(PathFile, "fasta", NewName, "nexus" ,molecule_type="DNA")
        print("Converted %i records" % count)
        file_list.append(NewName)

nexi = [(fname, Nexus.Nexus(fname)) for fname in file_list]
combined = Nexus.combine(nexi)
combined.write_nexus_data(filename=open('%s.nex' % OUT, 'w'),append_sets=False, mrbayes=False, codons_block=False)
combined.export_fasta(filename='%s.fa' % OUT, width=70)
combined.export_phylip(filename='%s.phy' % OUT)
with open('%s.txt' % OUT, 'w') as f:
    print(combined.append_sets(), file = f)
