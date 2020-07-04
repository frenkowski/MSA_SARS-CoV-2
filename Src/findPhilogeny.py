import numpy as np 
from anytree import Node, RenderTree, find_by_attr
from anytree.dotexport import RenderTreeGraph

def createMatrix(alignments, diff_by_gene_relative):
	'''
		Returns the features matrix, i.e. a matrix having reads as rows, 
		and diffs as columns. Matrix[i][j] = 1 if the i-eth read has 
		the j-eth diff, 0 otherwise. In order to remember to which
		row each reads corresponds to, a list of sequence_ids is also returned,
		such that sequences_list[i] corresponds to the i-eth row.

		:param alignments:
            A dict of the form {align_id_1 : value, align_id_2 : value, etc.}

        :param diff_by_gene_relative:
        	A dict of differences {gene_id:[differences]} where information
            about the different codons involved in the difference is provided

        :returns:
	        A tuple where:
	                - the first element is the matrix as described above
	                
	                - the second element is a list of sequences_id,
	                  such that the i-eth row in the matrix corresponds to
	                  the i-eth sequence_id
	'''

	
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
	'''
		Returns wheter or not a given matrix contains the forbidden matrix

		:param original_matrix:
            A matrix of features 0/1

        :returns:
	        A boolean whose value is True if original_matrix contains
	        the forbidden matrix, False otherwise
	'''

	# if original_matrix is None:
	# 	return False

	n_one = np.array([])

	### STEP 1: Sort original_matrix by column, according to the number
	###		    of 1 in each column 

	transpose_matrix = original_matrix.transpose()
	#Count number of 1s in each column
	for original_col in transpose_matrix:
		cont = 0
		for val in original_col:
			if val == 1:
				cont = cont +1	
		n_one = np.append(n_one,cont)		

	# Sort by descending number of 1s
	n_one_sort = sorted(n_one, reverse= True)

	# Holds the column number in the sorted matrix
	pos_new = 0
	aux_matrix = np.zeros(original_matrix.shape)

	for i in n_one_sort:
		
		# Holds the column number in the original matrix
		pos_old = 0
		
		# Find the column number in the original matrix
		for j in n_one:
			if i == j:
				# Column has been found
				# -1 is a sentinel value
				n_one[pos_old] = -1

				# Copy the old column in the new one
				row = 0
				while row != len(transpose_matrix[0]):
					aux_matrix[row][pos_new] = original_matrix[row][pos_old]
					row = row + 1

				# Move to the next column
				pos_new = pos_new + 1
			
			pos_old = pos_old + 1

	#print("Sorted matrix\n", aux_matrix)

	### STEP 2: Find, for each 1 in the sorted matrix, the position
	###		    of the previous 1 in the same row

	#Create auxillary matrix
	row = 0

	aux_matrix = np.zeros(original_matrix.shape)

	for curr_line in original_matrix:
		col = 0
		
		#contains the last column where 1 was found
		cont = -1   
		
		for val in curr_line:
			if val == 1:
				# set to last column with 1
				aux_matrix[row][col] = cont
				# store the current column index
				cont = col + 1
			else:
				aux_matrix[row][col] = 0	
			col = col + 1
		
		row = row + 1

	#print("Aux matrix\n", aux_matrix)


	
	### STEP 3: Check if the forbidden matrix is actually present
	###         in aux_matrix

	is_prohibited = False

	# Check column by column
	for j in range(0, aux_matrix.shape[1]):
		# For each row
		for i in range(0, aux_matrix.shape[0]):
			if aux_matrix[i][j] == 0:
				continue
			else:
				# For all other rows
				for k in range (0, aux_matrix.shape[0]):
					if aux_matrix[i][j] != aux_matrix[k][j] and aux_matrix[k][j] != 0:
						is_prohibited = True
						break

		if is_prohibited:
			break

	return is_prohibited


def createTree(sequences_list, original_matrix):
	'''
		Returns the philogenetic tree obtained from a 0/1 matrix

		:param original_matrix:
			A matrix of features 0/1

		:param sequences_list:
            a list of sequences_id, such that the i-eth row
            in the matrix corresponds to the i-eth sequence_id

        :returns:
	        A list of edges (Node1, Node2)
	'''

	# Initialize root
	root = dict()
	root['id'] = "Root"
	root['reads'] = list() # Should always be empty
	root['edges'] = list()

	### STEP 1: Create tree

	# For each row
	for i in range(0,original_matrix.shape[0]):

		# The first 1 in each row starts from Root
		curr_node = root

		# For each column
		for j in range(0,original_matrix.shape[1]):
			# Initialize dict for current column
			node_to_add = dict()

			if original_matrix[i][j] == 1:
				
				# Create node
				node_to_add['id'] = "C" + str(j+1)
				node_to_add['reads'] = []
				node_to_add['reads'].append(sequences_list[i])
				node_to_add['edges'] = []

				# Remove read from curr_node
				# it will "move down" the tree
				if sequences_list[i] in curr_node['reads']:
						curr_node['reads'].remove(sequences_list[i])

				found = False
				# Check if edge already exists
				for adj_node in curr_node['edges']:
					if node_to_add['id'] == adj_node['id']:
						found = True
						# just add the read to the node
						adj_node['reads'].append(sequences_list[i])
						# restart from said node
						curr_node = adj_node
						break

				# Create new node
				if not found:
					# add an edge to new node
					curr_node['edges'].append(node_to_add)	
					# restart from said node	
					curr_node = node_to_add
				
				

	### STEP 2: Convert to list of edges (start,end)

	result = visit(root)

	#print("Root is",root)
	return result

def visit(tree):
	'''
		Performs a visit on the resulting tree, in order to obtain
		a list of pairs (starting_node, ending_node) to feed in
		to the tree visualization algorithm.

		:param tree:
            a tree where nodes are represented as dicts, and each has a 
            list edges, where the adjacent nodes are present

        :returns:
	        A list of pairs (starting_node, ending_node)

	'''
	result = []

	for read in tree['reads']:
		result.append((tree['id'],read))

	for adj_node in tree['edges']:
		result.append((tree['id'],adj_node['id']))
		result.extend(visit(adj_node))

	return result

def findBiggestNotForbidden(matrix):
	new_matrix = matrix[:,0:1].astype(int)

	for i in range(1,matrix.shape[1]):
		curr_column = matrix[:,i:i+1].astype(int)
		
		new_matrix = np.concatenate((new_matrix, curr_column),axis=1)

		if containsForbidden(new_matrix):
			new_matrix = new_matrix[:,:-1]
			
	return new_matrix

def printTree(edges_list):
	# create root
	root = Node("Root")

	for (node1, node2) in edges_list:
		Node(node2, parent=find_by_attr(root, node1))
	for pre, _, node in RenderTree(root):
		print("%s%s" % (pre, node.name))

	# Print to file --- GRAPVIZ MUST BE INSTALLED!!!
	#RenderTreeGraph(root).to_picture("./tree.png")


def main():
	# m1 = np.array([
	# 	[0,1],
	# 	[1,0],
	# 	[1,1]
	# ])

	# print("Is forbidden?",containsForbidden(m1))

	m2 = np.array([
		[1,1],
		[1,0],
		[1,1]
	])

	print("Is forbidden?",containsForbidden(m2))
	#print("Tree:",createTree(["read0","read1","read2"],m2))
	nodes = createTree(["read0","read1","read2"],m2)
	printTree(nodes)

if __name__ == "__main__":
    main()