#!/usr/bin/env python
from Bio import SeqIO
import argparse

# Prepare input arguments
parser = argparse.ArgumentParser(description=
                                 '''This script takes a list of embl files and converts to gb''')
parser.add_argument('files', help ="Text file with names of gb or embl files to merge")
args = parser.parse_args()

#Read file names
file_handle = open(args.files, 'r')
files = [line.rstrip() for line in file_handle.readlines()]

#Extract info from files, write as gb
sequence_entries = [SeqIO.read(file, 'embl') for file in files]
for seq in sequence_entries:
	print(seq.id)
	filename = seq.id + ".gb"
	SeqIO.write(seq, filename, 'gb')

##Open merged sequence record with first gb or embl record
#combined_record = sequence_entries[0]

#Add subsequent ones
#for seq in sequence_entries[1:]:
#	combined_record = combined_record + ("N"*50) + seq
	
#combined_record.id = args.output

# write merged output file
#SeqIO.write(combined_record, args.output, args.out_format)
#print("Merged entries written to " + args.output)