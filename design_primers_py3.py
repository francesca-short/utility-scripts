#!/usr/bin/env python
import csv
import re
import argparse
import sys
import math
from tabulate import tabulate

parser = argparse.ArgumentParser(description='''This script reads in a DNA sequence and returns primers to use for that sequence that fit supplied parameters. Multiple targets can be used if they are in the same fasta file.
                              Inputs:
                                1) Target DNA sequence in .fasta format.
                                ''')
parser.add_argument('target', help="Target DNA sequence in .fasta format")
parser.add_argument('-o', '--output', help="Output file prefix")
parser.add_argument('-f', '--forward', type = int, help="Optional primer design Fwd location (default = start)", default = 0)
parser.add_argument('-r', '--reverse', type = int, help="Optional primer design Rev location (default = end)")
parser.add_argument('-tm', '--melt_temp', type = float, help="Optional primer design target melt temp (default = 55)", default = 55)
parser.add_argument('-tol', '--tolerance', type = float, help="Tm tolerance either side of optimum (default = 1)", default = 1)
parser.add_argument('-s', '--strict_location', type = bool, help="Strict location matching, options: 1 (will only use this exact site, returns the closest tm to optimum even outside tolerance), 0 (will test nearby sites until it finds a primer with tm within range)")
args = parser.parse_args()
if args.output is None:
    args.output = (args.target + "_primers")

# load target sequence, save name and sequence (no whitespace, upper case) as variables
targets = {}
with open(args.target, 'r') as f:
	target = f.readlines()
	target_name = re.sub('>', '', target[0].rstrip())
	target_sequence_processed = "".join([re.sub(r"\s+", "", item.rstrip().upper()) for item in target[1:]])
if args.reverse is None:
	args.reverse = len(target_sequence_processed)
print("Read in targets: " + args.target)
target_name = re.sub(r'.fa', '', args.target)

# Reverse complement function: changes each base to its complement, then returns the
# list of complemented bases in reverse order
def reverse_complement(seq):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    letters = [basecomplement[base] for base in seq]
    return "".join(letters[::-1])

# Load common restriction sites
restriction_enzymes = dict(AatII="GACGTC", Acc65I="GGTACC", AclI="AACGTT", AfeI="AGCGCT", AflII="CTTAAG", AgeI="ACCGGT",
                           ApaI="GGGCCC", ApaLI="GTGCAC", AscI="GGCGCGCC", AseI="ATTAAT", AvrII="CCTAGG",
                           BamHI="GGATCC", BclI="TGATCA", BglII="AGATCT", BmtI="GCTAGC", BsiWI="CGTACG", BspDI="ATCGAT",
                           ClaI="ATCGAT", BspEI="TCCGGA", BspHI="TCATGA", BsrGI="TGTACA", BssHII="GCGCGC",
                           BstBI="TTCGAA", BstZ17I="GTATAC", DraI="TTTAAA", EagI="CGGCCG", Eco53kI="GAGCTC",
                           EcoRI="GAATTC", EcoRV="GATATC", FseI="GGCCGGCC", FspI="TGCGCA", HindIII="AAGCTT",
                           HpaI="GTTAAC", KasI="GGCGCC", KpnI="GGTACC", MfeI="CAATTG", MluI="ACGCGT", MscI="TGGCCA",
                           NaeI="GCCGGC", NarI="GGCGCC", BbvCI="CCTCAGC", BsmI="GAATGC", BsrDI="GCAATG", BssSI="CACGAG",
                           BtsI="GCAGTG", NcoI="CCATGG", NdeI="CATATG", NgoMIV="GCCGGC", NheI="GCTAGC", NotI="GCGGCCGC",
                           NruI="TCGCGA", NsiI="ATGCAT", PacI="TTAATTAA", PciI="ACATGT", PluTI="GGCGCC",
                           PmeI="GTTTAAAC", PmlI="CACGTG", PsiI="TTATAA", PspOMI="GGGCCC", PstI="CTGCAG", PvuI="CGATCG",
                           PvuII="CAGCTG", SacI="GAGCTC", SacII="CCGCGG", SalI="GTCGAC", SbfI="CCTGCAGG", ScaI="AGTACT",
                           SfoI="GGCGCC", SmaI="CCCGGG", SnaBI="TACGTA", SpeI="ACTAGT", SphI="GCATGC", SrfI="GCCCGGGC",
                           SspI="AATATT", StuI="AGGCCT", SwaI="ATTTAAAT", TspMI="CCCGGG", XmaI="CCCGGG", XbaI="TCTAGA",
                           XhoI="CTCGAG", PaeR7I="CTCGAG", ZraI="GACGTC")

# Find restriction sites, returns sorted list of hits
def find_restriction_sites(sequence):
    restriction_sites = []
    for key, value in restriction_enzymes.items():
        if str(value) in sequence:
            restriction_sites.append(key)
        else:
            continue
    return sorted(restriction_sites)

# Calculate Tm, salt-adjusted assumes concentration 0.05M Na
def calc_tm(primer):
    count_C = primer.count("C")
    count_G = primer.count("G")
    return 100.5 + (41*(count_G + count_C)/len(primer)) - (820/len(primer)) + 16.6*math.log10(0.05)

fwd_site = args.forward

#Make primer design function, given a specific location
def design_primer(target, location):
	n = 16
	primer = target[location:(location + n)]
	while (float(calc_tm(primer)) < (args.melt_temp - args.tolerance)) & (n < 22):
		primer = target[location:(location + n)]
		n += 1
		if (args.melt_temp - args.tolerance) <= calc_tm(primer) <= (args.melt_temp + args.tolerance):
			return(primer)
		else:
			continue

primers = []
if design_primer(target_sequence_processed, fwd_site) is not None:
	fwd_primer = design_primer(target_sequence_processed, fwd_site)
	location = fwd_site + 1
	print("Found forward primer at exact location: " + fwd_primer)
	primers.append([target_name,"Fwd", fwd_primer, location, calc_tm(fwd_primer)])
else:
	while design_primer(target_sequence_processed, fwd_site) is None:
		fwd_site +=1
		if design_primer(target_sequence_processed, fwd_site) is not None:
			fwd_primer = design_primer(target_sequence_processed, fwd_site)
			location = fwd_site + 1
			print("Found forward primer at nearby location: " + fwd_primer)
			primers.append([target_name,"Fwd", fwd_primer, location, calc_tm(fwd_primer)])
		else:
			continue

rev_site = args.reverse
if rev_site == len(target_sequence_processed):
	revcomp_rev_site = 0
else:
	revcomp_rev_site = (len(target_sequence_processed) - rev_site)
rev_target = reverse_complement(target_sequence_processed)
if design_primer(rev_target, revcomp_rev_site) is not None:
	rev_primer = design_primer(rev_target, revcomp_rev_site)
	location = rev_site
	print("Found reverse primer at exact location: " + rev_primer)
	primers.append([target_name,"Rev", rev_primer, location, calc_tm(rev_primer)])
else:
	while design_primer(rev_target, revcomp_rev_site) is None:
		revcomp_rev_site +=1
		if design_primer(rev_target, revcomp_rev_site) is None:
			continue
		else:
			rev_primer = design_primer(rev_target, revcomp_rev_site)
			print("Found reverse primer at nearby location: " + rev_primer)
			location = (len(target_sequence_processed) - (revcomp_rev_site))
			primers.append([target_name,"Rev", rev_primer, location, calc_tm(rev_primer)])

#Open output csv file
outputfilename = (args.output + ".csv")
csvfile = open(outputfilename, 'w', encoding='utf8', newline='')
writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
header = ["Target", "Direction", "Sequence", "Start", "Tm"]
writer.writerow(header)

#Start table to print to screen
table = []
table.append(header)
for primer in primers:
	table.append(primer)
	writer.writerow(primer)
	
print(tabulate(table))

