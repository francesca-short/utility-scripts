#!/usr/bin/env python
from Bio import SeqIO
from Bio import SeqFeature
import argparse
import pprint
import re

# Prepare input arguments
parser = argparse.ArgumentParser(description=
                                 '''This script takes a list of gene names as .txt with one gene name
                                  per line and retrieves their upstream and downstream sequences from an embl or genbank file''')
parser.add_argument('reference', help="Reference file in .embl/.gb format")
parser.add_argument('gene_list', help="Gene list .txt file")
parser.add_argument('-o', '--output', help = 'output file prefix')
parser.add_argument('-n', '--number_bases', help = 'number of bases upstream and downstream to retrieve', type = int)
args = parser.parse_args()
if args.output == None:
    args.output = args.reference
if args.number_bases == None:
	args.number_bases = 600

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
def get_flanking_sequences(gene):
    count = 0
    for feature in reference_seqrecord.features:
        if gene[0].isalnum == False:
        	continue
        if gene in str(feature.qualifiers.values()):
            upstream = reference_seqrecord[int(feature.location.start - args.number_bases):int(feature.location.start)]
            upstream_sequence = str(upstream.seq)
            downstream = reference_seqrecord[int(feature.location.end):int(feature.location.end + args.number_bases)]
            downstream_sequence = str(downstream.seq)
            up_info = ("\n" + ">" + gene + "_" + str(args.number_bases) + "bp_5'" + "\n" + upstream_sequence)
            down_info = ("\n" + ">" + gene + "_" + str(args.number_bases) + "bp_3'" + "\n" + downstream_sequence)
            print(up_info, down_info)
            line = (up_info, down_info)
            file.writelines(line)
            count +=1
    if count ==0:
        print(gene + " not found in " + args.reference)

#loop over gene list
for gene in sorted(genes):
	get_flanking_sequences(gene)

file.close()



print("Gene info written to " + outputfilename)