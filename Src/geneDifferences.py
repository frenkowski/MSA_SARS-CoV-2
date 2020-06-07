### WILL BE USED IN PART 2
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

def findDifferencesbyGene2(differences, genes):
    """
        Store all differences by gene
        :param differences:
            A list of differences, with 'where' field [differences]
        :param genes:
            A dict of genes {gene_id1 : [start_pos, end_pos]}
        :returns:
            A dict of differences by gene {'gene_id1 : [differences]}
    """
    differences_by_gene = dict()

    # Initize dict for each gene
    for gene_id in genes.keys():
        differences_by_gene[gene_id] = list()
    
    for curr_diff in differences:
        for gene_id in genes.keys():
            gene_start = genes[gene_id][0]
            gene_stop = genes[gene_id][1]            
            if (curr_diff['start']+1 >= gene_start) and (curr_diff['start']+1+curr_diff['length'] < gene_stop):
                differences_by_gene[gene_id].append(curr_diff)
                break

    # NOTE: differences not in any gene are removed
    return differences_by_gene

def transcribeSequence(sequence):
    """
        Performs the transcription of a given sequence
        :param sequence:
            The sequence to be transcribed
        :returns:
            The transcribed sequence
    """

    result = ""
    
    for base in sequence:
        if base == 'T':
            result += 'U'
        else:
            result += base

    return result

def splitSequenceByCds(sequence, genes):
    """
        Splits a transcribed sequence according to the CDSs
        :param sequence:
            The sequence to be splitted (preferably transcribed)
        :param genes:
            A dict of genes {gene_id:[start,end]}
        :returns:
            A dict of CDSs {gene_id:[TranscribedSubSequence]}
    """

    seqs_by_cds = dict()
    # Initize dict for each gene
    for gene_id in genes.keys():
        seqs_by_cds[gene_id] = ""

    for gene_id in genes.keys():
            gene_start = genes[gene_id][0]
            gene_stop = genes[gene_id][1]            
           
            # Split the sequence according to the gene CDS
            # -1 should be necessary since 1-based on FASTA annotation
            seqs_by_cds[gene_id] = sequence[gene_start-1:gene_stop-1]
    
    return seqs_by_cds

def findDifferencesRelativePos(diff_by_gene, genes):
    """
        Replaces 'start' in each difference with relative position with respect to gene CDS
        :param diff_by_gene:
            A dict of differences {gene_id:[differences]}
        :param genes:
            A dict of genes {gene_id:[start,end]}
        :returns:
            A dict of CDSs {gene_id:[TranscribedSubSequence]}
    """   
    
    for gene_id in diff_by_gene.keys():
        gene_start = genes[gene_id][0]
        gene_stop = genes[gene_id][1]
        for curr_diff in diff_by_gene[gene_id]:
            curr_diff['start'] = curr_diff['start'] - gene_start
    
    return diff_by_gene


def findTranscriptDifferences(seqs_by_cds, diff_by_gene_relative, genes):

    result = dict()
    for gene_id in genes.keys():
        result[gene_id] = list()

    codons = dict()
    for gene_id in genes.keys():
        codons[gene_id] = list()

    ### STEP 1: Split into codons (NOTE: this is the reference, aka "correct" codons)
    for gene_id in genes.keys():
        curr_seq = seqs_by_cds[gene_id]
        # Split into groups of 3, aka codons
        codons[gene_id] = [curr_seq[i:i+3] for i in range(0, len(curr_seq), 3)]

    ### STEP 2: obtain codons of differences
    alternative_codons = dict()
    for gene_id in genes.keys():
        alternative_codons[gene_id] = list()

    for gene_id in diff_by_gene_relative.keys():
        for diff in diff_by_gene_relative[gene_id]:
            if diff['type'] == 'rep':
                # Create dict for storing the result
                codon_difference = dict()
                
                # First replace bases in sequence
                curr_seq = seqs_by_cds[gene_id]
                ## Remove old bases and insert new bases
                ## keep everything before, add difference, keep everything after
                ## TODO: check if it actually works as it should
                curr_seq = curr_seq[:diff['start']] + diff['seq'] + curr_seq[diff['start']+diff['len']:]

                # Now split into codons
                temp = [curr_seq[i:i+3] for i in range(0, len(curr_seq), 3)]
                codon_difference['codons'] = temp

                # add to alternative codons
                alternative_codons[gene_id].append(codon_difference)
                # TODO: other types of differences

    ### STEP 3 : compare ref and seq codon by codon
    for gene_id in alternative_codons.keys():
        for cod_diff in alternative_codons[gene_id]:
            for i in range(0, len(cod_diff))    #should also work with len(codons[gene_id])
                if codons[gene_id][i] != cod_diff[i]:
                    new_diff = dict()
                    # A different codon has been found
                    new_diff['ref-codon'] = codons[gene_id][i]
                    new_diff['seq-codon'] = cod_diff[i]
                    # new_diff['ref-protein'] = create_protein(new_diff['ref-codon'])
                    # new_diff['seq-protein'] = create_protein(new_diff['seq-codon'])
                    result[gene_id].append(new_diff)
    
    return result