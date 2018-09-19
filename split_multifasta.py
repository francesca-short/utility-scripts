#!/usr/bin/env python

import re
import argparse

parser = argparse.ArgumentParser(description=
                                 '''This script reads in a multifasta file and splits it into
                                 individual fasta files. No sanity checking is done for the final
                                 file count. Split fasta files are saved according to the
                                 header in the file with no suffix or unique identifier, so if you have
                                 >1 entry with the same name then all but one will be overwritten. No changes
                                 are made to the whitespace or line breaks in the individual sequence entries.
                                 
                                 Input:
                                 1) multifasta file (.fa or .fasta, not compressed)
                                 ''')
parser.add_argument('infile', help="Target DNA sequences in .fasta format")
args = parser.parse_args()

if '.fa' in args.infile:
    file_name = args.infile
elif '.fasta' in args.infile:
    file_name = args.infile
else:
    sys.exit("Input doesn't have .fa or .fasta extension")

file = open("init.txt", 'w')

with open(file_name) as f:
    for line in f:
        if '>' in line:
            file.close()
            filename = re.sub(r'>', '', line.rstrip()) + ".fa"
            file_text = line
            file = open(filename, 'w')
            print("Writing " + filename)
            file.writelines(file_text)
        else:
            file.writelines(line)