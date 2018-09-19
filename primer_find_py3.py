#!/usr/bin/env python
import csv
import re
import argparse
import sys
import math

parser = argparse.ArgumentParser(description='''This script reads in a list of primers and a DNA sequence, and reports
                              primers with matches to the target sequence.
                              Inputs:
                                1) Target DNA sequence in .fasta format.
                                2) Primer list in csv format with names in first column and sequences in second column.
                                ''')
parser.add_argument('target', help="Target DNA sequence in .fasta format")
parser.add_argument('primers', help="List of primers to be scanned in .csv format")
parser.add_argument('-o', '--output', help="Output file prefix")
args = parser.parse_args()
if args.output is None:
    args.output = (args.target + "_primer_matches")

# open dictionary for primer names and sequences
primers = {}

# open csv file and save relevant info into dictionary, convert sequences to uppercase no whitespace
with open(args.primers, 'r') as f:
    reader = csv.reader(f, delimiter=',', dialect=csv.excel_tab,
                        quoting=csv.QUOTE_NONE)
    # for each row, add the primer name (column 1) as the dict key, and the processed sequence (column 2) as the value
    for row in reader:
        primers.update({row[0]: re.sub(r"\s+", "", row[1].upper())})
print("{0} primer sequences loaded.".format(len(primers)))

# Take out any primer sequences with non ATGC characters
primers_filtered = dict((primer_names, primer_sequences) for primer_names, primer_sequences in primers.items() if
                        re.search('[^ATGC]', primer_sequences) is None)
print("Removed {0} primer sequences with non-ATGC nucleotides.".format(
    str(len(primers) - len(primers_filtered))))

# load target sequence, save name and sequence (no whitespace, upper case) as variables
with open(args.target, 'r') as f:
    target = f.readlines()
    target_name = re.sub('>', '', target[0].rstrip())
    target_sequence_processed = "".join([re.sub(r"\s+", "", item.rstrip().upper()) for item in target[1:]])


# Reverse complement function: changes each base to its complement, then returns the
# list of complemented bases in reverse order
def reverse_complement(seq):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    letters = [basecomplement[base] for base in seq]
    return "".join(letters[::-1])


# Function to find matches to primer sequence. Starts with the whole primer, then progressively chops from the 5' end
#  (ie. to remove restriction sites and other trailing nucleotides). Reports matches that are 14 nucleotides or longer.
def find_binding_sites_fwd(primer):
    n = 0
    primer_seed = primer[n:]
    while len(primer) >= len(primer_seed) > 13:
        p = re.compile(str(primer_seed))
        m = p.search(target_sequence_processed)
        if m is not None:
            return m, primer_seed
        else:
            n += 1
            primer_seed = primer[n:]


# Same function but checks for matches to reverse complement of primer
def find_binding_sites_rev(primer):
    n = 0
    primer_seed = primer[n:]
    while len(primer) >= len(primer_seed) > 13:
        p = re.compile(str(reverse_complement(primer_seed)))
        m = p.search(target_sequence_processed)
        if m is not None:
            return m, reverse_complement(primer_seed)
        else:
            n += 1
            primer_seed = primer[n:]


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

# Tests each primer for matches to the target
matches = []
for key, value in primers_filtered.items():
    match_fwd = find_binding_sites_fwd(value)
    match_rev = find_binding_sites_rev(value)
    if match_fwd is None and match_rev is None:
        continue
    elif match_fwd is not None:
        sites = find_restriction_sites(value)
        tm = calc_tm(match_fwd[0].group())
        whole_tm = calc_tm(primers_filtered[key])
        matches.append([key, value, len(match_fwd[1]), match_fwd[0].span()[0], match_fwd[0].span()[1],
                        "Sense", ", ".join(sites), '{0:.2f}'.format(tm), '{0:.2f}'.format(whole_tm)])
    elif match_rev is not None:
        sites = find_restriction_sites(value)
        tm = calc_tm(match_rev[0].group())
        whole_tm = calc_tm(primers_filtered[key])
        matches.append([key, value, len(match_rev[1]), match_rev[0].span()[0], match_rev[0].span()[1],
                        "Complement", ", ".join(sites), '{0:.2f}'.format(tm), '{0:.2f}'.format(whole_tm)])

if matches is None:
    sys.exit("No matches to " + args.target + " found in " + args.primers)
elif matches is not None:
    print(str(len(matches)) + " matching primers found, writing to " + args.output + ".csv")

# open output csv file
outputfilename = (args.output + ".csv")
csvfile = open(outputfilename, 'w', encoding='utf8', newline='')
writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
header = ["Name", "Sequence", "Match length", "Start", "End", "Strand", "Restriction sites", "Match Tm", "Total Tm"]
writer.writerow(header)
print(header)

# sort matches alphabetically based on primer name and write each to a new line of a csv
matches = sorted(matches)
for match in matches:
    row = (match[0], match[1], match[2], match[3], match[4], match[5], match[6], match[7], match[8])
    print(row)
    writer.writerow(row)
