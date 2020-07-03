import numpy as np
from anytree import Node, RenderTree, find_by_attr

def orderMatrix(original_matrix):
	# Contains the number of 1s for each column
	n_one = np.array([])

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
					# print("Row", row)
					# print("Pos_old", pos_old)
					# print("Pos_new", pos_new)
					row = row + 1

				# Move to the next column
				pos_new = pos_new + 1
			
			pos_old = pos_old + 1

	#Ho la matrice ordinata 
	ordered_original_matrix = aux_matrix
	print('Original Matrix Ordered: ')
	print(ordered_original_matrix)
	print('\n')

	matrix_aux(ordered_original_matrix, original_matrix)

def matrix_aux(ordered_original_matrix, original_matrix):
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
	
	print('Matrix Aux')
	print(aux_matrix, '\n')

	prohibited(aux_matrix)

def prohibited(aux_matrix):
	# Verifica matrice proibita 
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

	print("Is prohibited", is_prohibited)
	return is_prohibited

def tree():

	# USA UN DIZIONARIO DEL TIPO
	#dizionario = {node_start:"", node_end:"", value:""}, dizionario['node-start'] = "Ciao"
	#lista_dizionari = [{...},{...},{...}] lista_dizionari[0], lista_dizionari[0]['node-start']


	#lista temporale che contiene [arco_entrante arco_uscente e valore]
	edge_list = []
	i = 0
	for curr_line in original_matrix: 
		j = 0
		
		first_edge = False

		#contiene la lista degli archi di una riga della matrice
		edges = list()
		#salvo in length la prima posizione di dove si trova 1
		while j + 1 < len(original_matrix[0]) and curr_line[j] != 1:
			if j + 1 < len(original_matrix[0]):
				j = j + 1

		if curr_line[j] == 1:
			tmp_length = j + 1

			if not edges:		
				edges.append(str(tmp_length))
	
		
			while tmp_length + 1 < len(original_matrix[0])+ 1:

				if curr_line[tmp_length] == 1:
					edges.append(str(tmp_length+ 1))
					#first_edge = True

				tmp_length = tmp_length + 1

			#radice - nome_arco -nodo
			if len(edges) == 1:
				edge_list.append("Root "+"C"+ str(edges[len(edges)-1])+" "+ str(i))

		#arco_partenza(per ricavare il valore) - arco_nuovo - nodo
		if len(edges) > 1:
			edge_list.append("C" + str(edges[len(edges)-2])+ " C"+str(edges[len(edges)-1])+ ' ' + str(i))

		i = i + 1

	print(edge_list)
	app_list = ['None']
	while app_list:
		app_list = []
		tmp_list = []
		for edge_value in edge_list:	
			head = edge_value.split(' ')[0]		
			if head == 'Root':
				tmp_list.append(edge_value.split(' ')[0] + ' ' + edge_value.split(' ')[2])
			else:
				for j in edge_list:
					change = False
					h = edge_value.split(' ')[0]
					if j.split(' ')[1] == h:
						tmp_list.append(j.split(' ')[2] + ' ' + edge_value.split(' ')[2])
						change = True	
				if not change:
					app_list.append('Root' +' ' + edge_value.split(' ')[0] + ' None' )
		
		if app_list:
			edge_list = edge_list + app_list

	#Ordino la lista
	final_list = []
	for i in tmp_list:
		if i.split(' ')[0] == 'Root':
			final_list.append(i)
			i = i.split(' ')[1]	
			for j in tmp_list:
					if i == j.split(' ')[0]:
						final_list.append(j)
						i = j.split(' ')[1]
						for h in tmp_list:
							if i == h.split(' ')[0]:
								final_list.append(h)

	print(final_list,'\n')

	root = Node('Root')
	for line in final_list:
		line = line.split(" ")
		Node("".join(line[1:]).strip(), parent=find_by_attr(root, line[0]))
	for pre, _, node in RenderTree(root):
		print("%s%s" % (pre, node.name))

order_matrix(original_matrix)