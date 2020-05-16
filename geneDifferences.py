def findDifferencesbyGene(compactDifferences, genes):
    """
        Store all differences by gene
        :param compactDifferences:
            A dict of differences {'seq_id':[differences]}
        :param genes:
            A dict of genes {gene_id1 : [start_pos, end_pos]}
        :returns:
            A dict of differences by gene {'gene_id1 : {'align_id1' : [differences]}}
    """
    differences_by_gene = dict()

    # Initize dict for each gene
    for gene_id in genes.keys():
        differences_by_gene[gene_id] = dict()
        # Initiliaze dict for each alignment
        for seq_id in compactDifferences.keys():
            differences_by_gene[gene_id][seq_id] = list()

    for seq_id in compactDifferences.keys():
        for i in range(0, len(compactDifferences[seq_id])):
            curr_diff = compactDifferences[seq_id][i]
            for gene_id in genes.keys():
                gene_start = genes[gene_id][0]
                gene_stop = genes[gene_id][1]
                # check if difference is in gene range
                if (curr_diff['start']+1 >= gene_start) and (curr_diff['start']+1+curr_diff['length'] < gene_stop):
                    # initialize list if not already initialized 
                    #if seq_id not in differences_by_gene[gene_id].keys():
                    #    differences_by_gene[gene_id][seq_id] = list()
                    # check if curr_diff has not already been found
                    if curr_diff not in differences_by_gene[gene_id][seq_id]:
                        differences_by_gene[gene_id][seq_id].append(curr_diff)

    return differences_by_gene

    # Find differences by gene
    # dif_by_gene = findDifferencesbyGene(compact, {'gene1':[5000,6000]})