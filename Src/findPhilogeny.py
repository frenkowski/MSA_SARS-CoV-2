import numpy as np 

def createMatrix(alignments, diff_by_gene_relative):
	
	### FIRST CREATE AN EMPTY MATRIX

	# Find total numeber of differences
	differences_no = 0
	for gene_id in diff_by_gene_relative.keys():
		for diff in diff_by_gene_relative[gene_id]:
			differences_no = differences_no + 1

	# Find total number of sequences
	sequences_no = len(alignments.keys())

	# Create matrix having as rows the sequences, and as columns the differences
	matrix = np.zeros((sequences_no,differences_no))

	
	### NOW FILL THE MATRIX

	# Convert sequences to list
	sequences_list = list(alignments.keys())

	# For each diff
	column = 0
	for gene_id in diff_by_gene_relative.keys():
		for diff in diff_by_gene_relative[gene_id]:
			# For all sequences where the diff happens
			for seq in diff['where']:
				# Find position in sequences_list (= row)
				row = sequences_list.index(seq)

				matrix[row][column] = 1

			column = column + 1

	print("Matrix",matrix.astype(int))
	np.savetxt('test.out', matrix, fmt="%d", delimiter=',')

	return (sequences_list, matrix) 