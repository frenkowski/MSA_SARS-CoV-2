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

def split_reference_by_cds(transcribedSequence, genes):
    """
        Splits a transcribed sequence according to the CDSs
        :param transcribedSequence:
            The transcribed sequence to be slit
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
            seqs_by_cds[gene_id] = transcribedSequence[gene_start-1:gene_stop-1]
    
    return seqs_by_cds


def findTranscriptDifferences(seqs_by_cds, diff_by_gene, genes):

    codons = dict()
    for gene_id in genes.keys():
        codons[gene_id] = list()

    ### STEP 1: Split into codons (NOTE: this is the reference, aka "correct" codons)
    for gene_id in genes.keys():
        curr_seq = seqs_by_cds[gene_id]
        # Split into groups of 3, aka codons
        codons[gene_id] = [curr_seq[i:i+3] for i in range(0, len(curr_seq), 3)]

    ### STEP 2: obtain differences in term of codons
    for gene_id in diff_by_gene.keys():




