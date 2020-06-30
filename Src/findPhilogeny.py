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
	# print("Sequences_no", sequences_no)
	# print("Differences_no",differences_no)

	
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

	# print("Matrix",matrix.astype(int))
	# print("Shape",matrix.shape)
	np.savetxt('test.out', matrix, fmt="%d", delimiter=',')

	return (sequences_list, matrix)


def containsForbidden(original_matrix):
	n_one = np.array([])
	aux_matrix = np.zeros(original_matrix.shape)
	
	# Transpose the original matrix so that we can work row-by-row
	# instead of column-by-column
	tr_matrix = original_matrix.transpose()

	# print('Matrice originale: ')
	# print(original_matrix)
	# print('\n')
	
	#Count numbers of 1 in each column
	for i in tr_matrix:
		cont = 0
		for j in i:
			if j == 1:
				cont = cont +1	
		n_one = np.append(n_one,cont)		

	#Sort columns by descending order of 1s 
	n_one_sort = sorted(n_one, reverse= True)
	k1 = 0
	for i in n_one_sort:
		k = 0
		for j in n_one:
			if i == j:
				n_one[k] = -1
				l = 0
				while l != len(tr_matrix[0]):
					aux_matrix[l][k1] = original_matrix[l][k]
					l = l + 1
				#print('\n')
				#aux_matrix[[k1]] = tr_matrix[[k]]
				k1 = k1 + 1
			k = k + 1
	

	#Transpose again to obtain the original matrix
	original_matrix = aux_matrix.transpose()
	# print('Matrice originale ordinata: ')
	# print(original_matrix)
	# print('\n')

	aux_matrix = np.zeros(original_matrix.shape)

	#Create auxillary matrix
	rows = 0
	for i in original_matrix:
		pos = 0
		cont = -2 #contains the position where the last 1 was found in the current row   
		for j in i:

			if j == 1:
				aux_matrix[rows][pos] = cont + 1
				cont = pos
			else:
				aux_matrix[rows][pos] = 0	
			pos = pos + 1
		rows = rows + 1
	# print('Matrice ausiliaria')
	# print(aux_matrix)


	#Check forbidden matrix 
	aux_matrix = aux_matrix.transpose()

	cont = False
	row = 0
	for i in aux_matrix:
		column = 0
		for j in i:
			column = column+1
			if j != 0:
				for k in i:
					if k != 0:
						if j != k:
							print("Row", column)
							print("Column", row)
							#cont = True
							#break
		row = row+1
		#if cont:
			#break
							

	if cont:
		return True
	else:
		return False

def createTree(sequences_list, original_matrix):
	pass