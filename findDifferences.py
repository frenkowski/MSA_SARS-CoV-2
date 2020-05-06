from collections import OrderedDict # Remembers the order entries were added
import re
import os

## Returns for each alignment all the positions where it differs from the reference
def findDifferences(reference, alignments):
	ref_id = reference.keys()[0]
	ref = reference[ref_id]
	differences = OrderedDict()
	for align_id in alignments.keys():
		differences[align_id] = []
	for i in range(0, len(ref)):
		for align_id in alignments.keys():
			curr_base = alignments[align_id][i]
			curr_diff = dict()
			if curr_base != ref[i]: # A difference has been found
				# Store the value of reference
				curr_diff['ref'] = '*'
				# Store the value of the alignment
				curr_diff['seq'] = '*' 
				# Store the position where it happens
				curr_diff['start'] = i

				# Determine the difference type
				if ref[i] == 'N':	
					curr_diff['type']="na1"	# ref not available
				elif curr_base == 'N':
					curr_diff['type']="na2" # seq not available
				elif ref[i] == 'N' and curr_base == 'N':
					curr_diff['type']="na3" # both ref and seq not available 
				elif ref[i] == '-' and curr_base != '-':
					curr_diff['type']='ins' # insertion
				elif ref[i] != '-' and curr_base == '-':
					curr_diff['type']='del' # deletion
				else:
					curr_diff['type']='rep' # replacement

				if curr_diff['type'] == 'na1':
					# only store seq since ref is not available
					curr_diff['seq'] = curr_base
				elif curr_diff['type'] == 'na2':
					# only store ref since seq is not available
					curr_diff['ref'] = ref[i]
				#elif curr_diff['type'] == 'na3': # do nothing
				elif curr_diff['type'] == 'ins':
					# store insert in seq
				 	curr_diff['seq'] = curr_base
				elif curr_diff['type'] == 'del':
				 	# store deleted in ref
				 	curr_diff['ref'] = ref[i]
				else: # base changed
				 	# store both
				 	curr_diff['seq'] = curr_base
				 	curr_diff['ref'] = ref[i]

				# Add the newfound difference to the other ones
				differences[align_id].append(curr_diff)
	return differences

## Compact for each alignment all the consecutive differences of the same type
def compactDifferences(differences):
	new_Differences = OrderedDict()
	for align_id in differences.keys():
		new_Differences[align_id] = []

	for align_id in new_Differences.keys():
		diff = differences[align_id]

		curr_diff = dict()
		if diff: # Not empty
			curr_diff['start'] = diff[0]['start']
			curr_diff['type'] = diff[0]['type']
			curr_diff['ref'] = diff[0]['ref']
			curr_diff['seq'] = diff[0]['seq']
			curr_diff['length'] = 1

			for i in range(1, len(diff)):
				# Try to extend difference
				if diff[i]['start'] == diff[i-1]['start'] + 1 and diff[i]['type'] == diff[i-1]['type']:
					# Increase current dif length
					curr_diff['length'] = curr_diff['length'] + 1

					if curr_diff['type'] == 'na1':
						# update seq
						curr_diff['seq'] = curr_diff['seq'] + diff[i]['seq']
						curr_diff['ref'] = '*'
					elif curr_diff['type'] == 'na2':
						# update ref
						curr_diff['ref'] = curr_diff['ref'] + diff[i]['ref']
						curr_diff['seq'] = '*'
					elif curr_diff['type'] == 'na3':
						curr_diff['ref'] = '*'
						curr_diff['seq'] = '*'
					elif curr_diff['type'] == 'ins':
						# Update seq
						curr_diff['seq'] = curr_diff['seq'] + diff[i]['seq']
						curr_diff['ref'] = '*'
					elif curr_diff['type'] == 'del':
						# Update ref
						curr_diff['ref'] = curr_diff['ref'] + diff[i]['ref']
						curr_diff['seq'] = '*'
					else:
						curr_diff['seq'] = curr_diff['seq'] + diff[i]['seq']
						curr_diff['ref'] = curr_diff['ref'] + diff[i]['ref']

				else: # New difference
					# Close current difference
					if curr_diff:
						new_Differences[align_id].append(curr_diff)
						curr_diff = {}
					# Store new difference
					curr_diff['start'] = diff[i]['start']
					curr_diff['type'] = diff[i]['type']
					curr_diff['ref'] = diff[i]['ref']
					curr_diff['seq'] = diff[i]['seq']
					curr_diff['length'] = 1

			if curr_diff:	# Last element still in curr_diff
				new_Differences[align_id].append(curr_diff)
				curr_diff = {}

	return new_Differences

## Print the differences to terminal
def printDifferences(reference, compactDifferences):
	ref_id = reference.keys()[0]
	ref = "##Ref={}, len={}".format(ref_id, len(reference[ref_id]))
	fields = "START\tLENGTH\tTYPE\tREF\tSEQ\t"
	value ="{}\t{}\t{}\t{}\t{}\t"
	print(ref)
	for align_id in compactDifferences.keys():
		seq = "##Seq={}".format(align_id)
		print(seq)
		print(fields)
		for diff in compactDifferences[align_id]:
			temp = value.format(diff['start']+1,diff['length'],diff['type'],diff['ref'],diff['seq'])
			print(temp)

def findStats(reference, alignments):
	# Initialize dict
	stats = dict()
	for align_id in alignments.keys():
		stats[align_id] = dict()
		stats[align_id]['matches'] = 0
		stats[align_id]['mismatches'] = 0
		stats[align_id]['na'] = 0

	# Extract reference from dict
	ref_id = reference.keys()[0]
	ref = reference[ref_id]

	for i in range(0, len(ref)):
		for align_id in alignments.keys():
			curr_base = alignments[align_id][i]
			if curr_base == 'N':
				stats[align_id]['na'] = stats[align_id]['na'] + 1
			elif curr_base == ref[i]:
				stats[align_id]['matches'] = stats[align_id]['matches'] + 1
			else: # curr_base != ref[i]
				stats[align_id]['mismatches'] = stats[align_id]['mismatches'] + 1

	return stats



## Write the differences to file
def writeToFile(fileName, reference, compactDifferences, stats):
	ref_id = reference.keys()[0]
	ref = "##Ref={}, len={}\n".format(ref_id, len(reference[ref_id]))
	fields = "START\tLENGTH\tTYPE\tREF\tSEQ\t\n"
	value ="{}\t{}\t{}\t{}\t{}\t\n"
	path = "{}.mad".format(fileName) # MAD = Multiple Alignment Difference
	with open(path, "w") as f: 
		f.write(ref)
		for align_id in compactDifferences.keys():
				seq = "##Seq={},Matches={},Mismatches={},NA={}\n".format(align_id,
																		stats[align_id]['matches'],
																		stats[align_id]['mismatches'],
																		stats[align_id]['na'])
				f.write(seq)
				f.write(fields)
				for diff in compactDifferences[align_id]:
					temp = value.format(diff['start']+1,diff['length'],diff['type'],diff['ref'],diff['seq'])
					f.write(temp)

def parseClustal(referenceId, fileName):
	alignments = dict()
	reference = dict()
	path = "{}.clw".format(fileName)
	with open(path, "r") as f:

		# Remove header
		lines = f.readlines()[1:]
		number_of_aligns = 0
		for line in lines:

			# remove (as many as possible) spaces, *, \n from beginning and end of string
			line = line.strip(" *\n")
			# replace tabs with spaces for split
			line = line.replace('\t',' ')
			# replace multiple spaces with single space for split
			line = re.sub(' +', ' ', line)
			
			if line != '':
				line_to_list = line.split()
				key = line_to_list[0]
				align = line_to_list[1]
				if len(line_to_list) == 3: # Note: only some aligners show base count
					num_of_bases = line_to_list[2]
				if key not in alignments.keys(): # Create dict entry
					alignments[key] = ''
				alignments[key] = alignments[key] + align

		ref = alignments.pop(referenceId)
		reference[referenceId] = ref
		
		return reference, alignments		

def main():
	path_alignments = "./Alignments/"
	path_output = "./Outputs/"
	for file in os.listdir(path_alignments):
	    if file.endswith(".clw"): # for each .clw file
	    	# Obtain filename, will be used for both input and output
	    	fileName = os.path.splitext(file)[0]
	    	# Parse reference and alignments from file
	        reference, alignments = parseClustal("NC_045512", path_alignments+fileName)
	        # Find for each alignment all single differences
	        diff = findDifferences(reference,alignments)
	        # Compact consecutive difference
	        compact = compactDifferences(diff)
	        # Find alignments stats
	        stats = findStats(reference, alignments)
	        # Write differences to file
	        writeToFile(path_output+fileName,reference,compact,stats)

if __name__ == "__main__":
    main()

# reference, alignments = parseClustal("NC_045512", "ClustalOmega.clustal_num")	
# #print(reference)
# print("\n")
# #print(alignments)	 

# diff = findDifferences(reference,alignments)
# compact = compactDifferences(diff)
# printDifferences(reference, compact)
# writeToFile("difficile",reference,compact)


# ref = dict()
# ref['ref'] = "acgt-a"
# aligns = OrderedDict()
# aligns["seq1"] = "ccgta-"
# aligns["seq2"] = "gggt-a"
# diff = findDifferences(ref,aligns)
# compact = compactDifferences(diff)
# printDifferences(ref, compact)

# ref = dict()
# ref['ref'] = "gataggat--ggg-catacgataggcccgtgatgcc"
# aligns = OrderedDict()
# aligns["seq1"] = "gatatgattagggccatacgataggcccgt---gcc"
# aligns["seq2"] = "gataggatt-gggccatacga-aggcccgtgatgcc"
# aligns["seq3"] = "gataggattagggccatacgataggcccgtgatgcc"
# diff = findDifferences(ref,aligns)
# compact = compactDifferences(diff)
# printDifferences(ref, compact)
# writeToFile("esempio",ref,compact)


# print("\n")

# ref = "a"
# aligns = {"seq1":"c","seq2":"a"}
# diff = findDifferences(ref,aligns)
# compact = compactDifferences(diff)
# printDifferences(ref, aligns, compact)

# print("\n")

# ref = "a"
# aligns = {"seq1":"c","seq2":"t"}
# diff = findDifferences(ref,aligns)
# compact = compactDifferences(diff)
# printDifferences(ref, aligns, compact)

# print("\n")

# ref = "aa"
# aligns = {"seq1":"gg","seq2":"tt"}
# diff = findDifferences(ref,aligns)
# compact = compactDifferences(diff)
# printDifferences(ref, aligns, compact)
# writeToFile("prova",ref, aligns, compact)
