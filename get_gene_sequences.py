#!/usr/bin/env python
from Bio import SeqIO
from Bio import SeqFeature
import argparse
import pprint
import re

# Prepare input arguments
parser = argparse.ArgumentParser(description=
                                 '''This script takes a list of gene names as .txt with one gene name
                                  per line and retrieves their sequences from an embl or genbank file''')
parser.add_argument('reference', help="Reference file in .embl/.gb format")
parser.add_argument('gene_list', help="Gene list .txt file")
parser.add_argument('-o', '--output', help = 'output file prefix')
args = parser.parse_args()
if args.output == None:
    args.output = args.reference

# Read in reference file (print error message if reference file wrong format)
if 'embl' in args.reference:
	reference_seqrecord = SeqIO.read(args.reference, 'embl')
elif 'gb' in args.reference:
	reference_seqrecord = SeqIO.read(args.reference, 'gb')
elif 'gbk' in args.reference:
	reference_seqrecord = SeqIO.read(args.reference, 'gbk')
else:
	sys.exit("Reference file must be in embl or genbank format with '.embl', '.gb' or '.gbk' file extension")
print("Read in reference file: " + args.reference)

# Read in gene names and put in list
file_handle = open(args.gene_list, 'r')
genes = set()
for gene in file_handle.readlines():
	genes.add(gene.rstrip())
print(genes)

# open output fasta file
outputfilename = (args.output + "_sequences.fa")
file = open(outputfilename, 'w')

#define function to save sequence from a particular gene
def get_sequence(gene):
    count = 0
    for feature in reference_seqrecord.features:
        if gene[0].isalnum == False:
        	continue
        if gene in str(feature.qualifiers.values()):
            seq = feature.extract(reference_seqrecord)
            info = ("\n" + ">" + gene + "\n" + seq.seq)
            print(info)
            file.writelines(info)
            count +=1
    if count ==0:
        print(gene + " not found in " + args.reference)

#loop over gene list
for gene in sorted(genes):
	print(gene)
	get_sequence(gene)

file.close()



print("Gene info written to " + outputfilename)