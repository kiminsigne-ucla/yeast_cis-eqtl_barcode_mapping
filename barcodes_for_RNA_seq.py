"""
This script selects barcodes that will be used for downstream analysis
in RNA-seq


Usage: python barcode_selection_for_rna_seq.py <reads file> <reference file> 
				<library> <output file>
reads file :    A file in FASTQ format 
reference file: A file in CSV format that contains the reference sequences. 
				  The first field must be the name of the sequence and the 
				  second field must be the sequence.
library:		  The name of the library to analyze. Valid arguments include: 
				  'all', 'sharon', 'firstTileAUG' and 'firstTileUTR'. These are
				  just the names of the sub-library as listed in the excel file,
				  so they just need to match the names in the reference file.
output file:    The name of the output file which will contain a list of the
				  final barcode
"""

import os
import sys
import collections
import re
import Levenshtein
import numpy
import random

#-------------------------------------------------------------------------------
def reverse_complement(string):
	"""
	Return the reverse complement of a string
	"""
	rc = ''
	reverse = string[::-1]
	for letter in reverse:
		if letter == 'A':
			rc += 'T'
		elif letter == 'T':
			rc += 'A'
		elif letter == 'C':
			rc += 'G'
		elif letter == 'G':
			rc += 'C'
	return rc

#-------------------------------------------------------------------------------
def bootstrap_levenshtein(reference, n):
	"""
	This function calculates the reference Levenshtein distribution. It randomly
	picks two sequences from the reference sequences and calculates the distance.
	"""
	distances = []
	# bootstrap n times
	for i in range(0, n):
		# randomly grab two sequences with replacement
		string1 = random.choice(reference)
		string2 = random.choice(reference)

		distances.append(Levenshtein.distance(string1, string2))
	
	# take cutoff at 1% percentile
	cutoff = numpy.percentile(distances, 1)

	# If the distribution consists of mainly large distances, the 1% percentile
	# will be large, so readjust the cutoff lower in this case
	
	if cutoff >= 50: cutoff = 11 
	# 11 is 1% percentile in control library, but choice is somewhat arbitrary
		 
	return cutoff


# read in arguments
all_reads = sys.argv[1]
reference = open(sys.argv[2], 'r')
library = sys.argv[3]
barcode_output = open(sys.argv[4], 'w')

# these output files outputs the full map of variants to barcodes and barcodes to
# variants, for use with statistics and plots in R
filename_barcode_map = library + '_variant_to_barcode.txt'
filename_variant_map = library + '_barcode_to_variant.txt'
variant_to_barcode_file = open(filename_barcode_map, 'w')
barcode_to_variant_file = open(filename_variant_map, 'w')

print "Processing input..."
# extract only the sequences from the fastq file for easier manipulation, aka
# every fourth line
command = 'sed -n \'2~4p\' '+ all_reads + ' > master_reads_seq_only.txt'
os.system(command)

all_reads = open('master_reads_seq_only.txt', 'r')

print "Reading in reference..."
# translate reference sequences to reverse complement
ref_seqs = collections.defaultdict(bool) # make this a dict for fast lookup
ref_list = [] # used later for bootstrapping

# read in reference sequences from specific library
if library == 'all':
	for line in reference:
		fields = line.split(',')
		seq = fields[1]

		ref_seq = reverse_complement(seq)
		ref_seqs[ref_seq] = True
		ref_list.append(seq)
else: 
	for line in reference:
		fields = line.split(',')
		name = fields[0]
		seq = fields[1]
		# only grab reference sequences from library of interest
		match = re.match(library, name)
		if match:
			ref_seq = reverse_complement(seq)
			ref_seqs[ref_seq] = True
			ref_list.append(seq)

print "Number of unique reference sequences: ", len(set(ref_list))
# grab reads that match a reference sequence
reference_reads = []

for read in all_reads:
	read = read.strip()
	# trim barcode and SbfI site, 20 + 8 first bp
	variant = read[28:]

	if ref_seqs[variant] == True:
		reference_reads.append(read)

all_reads.close()

# grab barcodes that map to a reference sequence
reference_barcodes = [string[0:20] for string in reference_reads]

print "Number of barcodes before filtering: ", len(set(reference_barcodes))
print "Filter by barcode frequency..."
# Count the frequency of barcodes 
barcode_counts = collections.defaultdict(int)
for barcode in reference_barcodes:
	barcode_counts[barcode] += 1

# Throw out barcodes that appear 1 or 2 times
# Make this a dictionary for fast lookup
potential_barcodes = collections.defaultdict(bool)

for barcode in barcode_counts.keys():
	if barcode_counts[barcode] >= 3:
		potential_barcodes[barcode] = True

print "Number of potential barcodes: ", len(potential_barcodes.keys())
print "Mapping reads to potential barcodes..."

# Go through all the reads again and map sequences to potential barcodes
barcode_map = collections.defaultdict(list)

# keep track which barcodes map to which variants, useful for statistics
variant_map = collections.defaultdict(list)


all_reads = open(sys.argv[1], 'r')

for read in all_reads:
	read = read.strip()
	barcode = read[0:20]
	variant = read[28:]

	if potential_barcodes[barcode] == True:

		variants = barcode_map[barcode]
		variants.append(variant)
		barcode_map[barcode] = variants

# bootstrap reference sequences to get a reference Levenshtein distribution 
# to determine cutoff
print "Bootstrapping reference sequences to obtain cutoff...", 
cutoff = bootstrap_levenshtein(ref_list, 100000)
print "cutoff is Levenshtein distance ", cutoff

# write header to output file
fields = ['barcode', 'unique_sequences', 'num_reads', 'num_reads_most_common', 
		  'is_reference', 'distance\n']
barcode_output.write('\t'.join(fields))


# for each barcode, calculate the Levenshtein distance between the variants 
# that map to that barcode
final_barcodes = []

print "Calculating Levenshtein distance to select final barcodes..."
for barcode in barcode_map.keys():

	variants = barcode_map[barcode]
	total_reads = len(variants)

	# calculate distance between promoter sequence only, so 
	# trim primer sites, the first 15bp and last 15bp
	variants = [string[15:] for string in variants]
	variants = [string[:-15] for string in variants]

	counts = collections.Counter(variants)
	most_common = counts.most_common(1)
	most_common_read = most_common[0][0]
	most_common_count = most_common[0][1]

	# remove first occurrence of most common sequence, use as reference
	variants.remove(most_common_read)

	# calculate Levenshtein distance between the most common sequence and all
	# other variants that map to this barcode
	distances = [Levenshtein.distance(most_common_read, variant) for variant in variants]

	max_dist = max(distances)
	# throw out any barcodes with a max Levenshtein distance greater than cutoff
	if max_dist <= cutoff:
		final_barcodes.append(barcode)
	else:
		del barcode_map[barcode]


print "Outputting results..."
# write barcodes to file
for barcode in final_barcodes:
	barcode_output.write(barcode+'\n')

print "Number of final barcodes: ", len(final_barcodes)


print "Writing full output files..."

# write headers to output
variant_to_barcode_file.write('barcode\tunique_sequences\tnum_reads\n')
barcode_to_variant_file.write('variant\tunique_barcodes\tnum_barcodes\n')


for barcode in barcode_map.keys():

	seqs = barcode_map[barcode]
	unique_seqs = set(seqs) # output only the set of unique sequences

	info = [barcode, ','.join(unique_seqs), str(len(seqs))]
	variant_to_barcode_file.write('\t'.join(info) + '\n')

# # Given these final barcodes, how many reference sequences do we see?
for read in reference_reads:

	barcode = read[0:20]
	variant = read[28:]

	if barcode in final_barcodes:
		if ref_seqs[variant] == True:
			barcodes = variant_map[variant]
			barcodes.append(barcode)
			variant_map[variant] = barcodes

for variant in variant_map.keys():
	barcodes = variant_map[variant]
	num_barcodes = len(barcodes)
	unique_barcodes = set(barcodes)
	info = [variant, ','.join(unique_barcodes), str(num_barcodes)]
	barcode_to_variant_file.write('\t'.join(info)+'\n')
















