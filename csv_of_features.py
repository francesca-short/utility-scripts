#!/usr/bin/env python
from Bio import SeqIO
from Bio import SeqFeature
import sys
import argparse
import pprint
import csv
import re

# Prepare input arguments
parser = argparse.ArgumentParser(description=
                            '''This script reads in an embl or genbank file and outputs a csv file
                            containing the type, strand, start, end and all other information
                            for each feature. The extra columns (function, locus tag etc) are
                            the qualifiers from the genome file and will differ between genomes.''')
parser.add_argument('reference', help="Reference file in .embl/.gb format")
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
	sys.exit('''Reference file must be in embl or genbank format with 
	'.embl', '.gb' or '.gbk' file extension''')
print("Read in reference file: " + args.reference)

# open output csv file
outputfilename = (args.output + "_CDS_annotations.csv")
csvfile = open(outputfilename, 'w', encoding='utf8', newline='')
writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
header = ["Type","Strand","Start","End", "Gene"]
for feature in reference_seqrecord.features[0:10]:
    if feature.type == 'CDS':
        qualifiers = (feature.qualifiers.keys())
        break
    else:
        continue
for entry in qualifiers:
    header.append(entry)
writer.writerow(header)

# open output csv file
#outputfilename = (args.output + "_gene_insert_sites.csv")
#csvfile = open(outputfilename, 'w', encoding='utf8', newline='')
#writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
#header = (
#"locus_tag", "gene_name", "type", "start", "end", "strand", "read_count", "ins_index", "gene_length", "ins_count","function")
#writer.writerow(header)

def write_info_to_csv(gene):
    qualifier_values = []
    for key in header[4:]:
        value = str(gene.qualifiers.get(key))
        print(value)
        value = re.sub('(\[\')' ,'', re.sub('(\'\])', '', value))
        qualifier_values.append(value)
    print(qualifier_values)
    row = [gene.type, gene.location.strand, int(gene.location.start),
           int(gene.location.end), gene.qualifiers.get('gene', "")]
    for value in qualifier_values:
        row.append(value)
    print(row)
    writer.writerow(row)

#call info-getting function for every feature except for excluded 
#categories (repeats etc). Add categories to exclude here
#if your genome includes something weird
for gene in reference_seqrecord.features:
     if gene.type == 'source':
         continue
     if gene.type == 'repeat_region':
         continue
     if gene.type == 'mobile_element':
         continue
     write_info_to_csv(gene)

print("Gene info written to " + outputfilename)