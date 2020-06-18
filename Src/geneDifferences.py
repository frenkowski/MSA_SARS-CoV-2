### WILL BE USED IN PART 2

amino_acids = {'UUU':'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG':'Leu',
        'UCU':'Ser', 'UUC':'Ser', 'UCA':'Ser', 'UCG':'Ser',
        'UAU':'Tyr', 'UAC':'Tyr', 'UAA':'stop', 'UAG':'stop',
        'UGU':'Cys', 'UGC':'Cys', 'UGA':'stop', 'UGG':'stop',
        'CUU':'Leu', 'CUC':'Leu', 'CUA':'Leu', 'CUG':'Leu',
        'CCU':'Pro', 'CCC':'Pro', 'CCA':'Pro', 'CCG':'Pro',
        'CAU':'His', 'CAC':'His', 'CAA':'Gin', 'CAG':'Gin',
        'CGU':'Arg', 'CGC':'Arg', 'CGA':'Arg', 'CGG':'Arg',
        'AUU':'Ile', 'AUC':'Ile', 'AUA':'Ile', 'AUG':'Met',
        'ACU':'Thr', 'ACC':'Thr', 'ACA':'Thr', 'ACG':'Thr',
        'AAU':'Asn', 'AAC':'Asn', 'AAA':'Lys', 'AAG':'Lys',
        'AGU':'Ser', 'AGC':'Ser', 'AGA':'Arg', 'AGG':'Arg',
        'GUU':'Val', 'GUC':'Val', 'GUA':'Val', 'GUG':'Val',
        'GCU':'Ala', 'GCC':'Ala', 'GCA':'Ala', 'GCG':'Ala',
        'GAU':'Asp', 'GAC':'Asp', 'GAA':'Glu', 'GAG':'Glu',
        'GGU':'Gly', 'GGC':'Gly', 'GGA':'Gly', 'GGG':'Gly',
        'UCC':'Ser'}

def findDifferencesbyGene(compactDifferences, genes):
    """
        Store all differences by gene
        :param compactDifferences:
            A dict of differences {'seq_id':[differences]}
        :param genes:
            A dict of genes {gene_id1 : [start_pos, end_pos]}
        :returns:
            A dict of differences by gene {'gene_id1 : [differences]}
    """
    differences_by_gene = dict()

    # Initize dict for each gene
    for gene_id in genes.keys():
        differences_by_gene[gene_id] = dict()
        # Initiliaze list of diffs for each alignment
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
        if curr_diff['type'] != 'na2':
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

    #print('T' in result)

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
            # +1 on the right since Python considers the right index as non inclusive
            seqs_by_cds[gene_id] = sequence[gene_start-1:gene_stop-1+1]

            #print(gene_id)
            #print(len(seqs_by_cds[gene_id])%3==0)
    
    return seqs_by_cds

def findDifferencesRelativePos(diff_by_gene, genes):
    """
        Adds new 'start_rel' field in each difference with relative position with respect to gene CDS
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
            curr_diff['start_rel'] = curr_diff['start'] - gene_start
    
    return diff_by_gene


def findTranscriptDifferences(seqs_by_cds, diff_by_gene_relative, genes, all_seqs_id):
    
    # Check if ref is %3

    #result = dict()
    #for gene_id in genes.keys():
    #    result[gene_id] = list()

    codons = dict()
    for gene_id in genes.keys():
        codons[gene_id] = list()

    ### STEP 1: Split into codons (NOTE: this is the reference, aka "correct" codons)
    for gene_id in genes.keys():
        curr_seq = seqs_by_cds[gene_id]
        # Split into groups of 3, aka codons
        codons[gene_id] = [curr_seq[i:i+3] for i in range(0, len(curr_seq), 3)]

    #print("Codons: ",codons)

    ### STEP 2: obtain codons of differences
    alternative_codons = dict()
    for gene_id in genes.keys():
        alternative_codons[gene_id] = list()

    for gene_id in diff_by_gene_relative.keys():
        for diff in diff_by_gene_relative[gene_id]:
            if diff['type'] == 'rep':
                # Create dict for storing the result
                #codon_difference = dict()
                
                # First replace bases in sequence
                curr_seq = seqs_by_cds[gene_id]

                #Replace eventual 'T' in diff['seq']
                #Note: differences are taken from starting DNA and not RNA, so they can contain T
                if 'T' in diff['seq']:
                    diff['seq'] = diff['seq'].replace('T','U')

                ## Remove old bases and insert new bases
                ## keep everything before, add difference, keep everything after
                curr_seq = curr_seq[:diff['start_rel']] + diff['seq'] + curr_seq[diff['start_rel']+diff['length']:]

                # Now split into codons
                temp = [curr_seq[i:i+3] for i in range(0, len(curr_seq), 3)]
                #codon_difference['codons'] = temp

                # print("Temp is: ", temp)
                diff['codons'] = temp
                #print("Diff['codons'] is: ", diff['codons'])

                # add to alternative codons
                #alternative_codons[gene_id].append(temp)
            
            elif diff['type'] == 'del':
                # Create dict for storing the result
                #codon_difference = dict()
                
                # First replace bases in sequence
                curr_seq = seqs_by_cds[gene_id]
                # Add indels in seq
                curr_seq = curr_seq[:diff['start_rel']] + '-'*diff['length'] + curr_seq[diff['start_rel']+diff['length']:]

                # Now split into codons
                temp = [curr_seq[i:i+3] for i in range(0, len(curr_seq), 3)]
                #codon_difference['codons'] = temp

                diff['codons'] = temp
                # add to alternative codons
                #alternative_codons[gene_id].append(temp)



    #print("Diff by gene relative: ", diff_by_gene_relative)
    #print("Alternative_codons:",alternative_codons)
    ### STEP 3 : compare ref and seq codon by codon
    for gene_id in diff_by_gene_relative.keys():
        # Store codons of current gene for an easier comparison
        codons_ref = codons[gene_id]
        # Compare codons of current gene with codons modified by diffs
        for diff in diff_by_gene_relative[gene_id]:
            #print("Diff: ",diff)
            cod_diff = dict()
            # Compare codon by codon
            i = 0
            diff['cod-info'] = list()
            for ref, dif in zip(codons_ref,diff['codons']):
                if ref != dif:
                    cod_diff['ref-cod'] = ref
                    cod_diff['ref-aa'] = amino_acids[ref]
                    cod_diff['dif-cod'] = dif
                    if '-' not in dif:
                        cod_diff['dif-aa'] = amino_acids[dif]
                    else:
                        cod_diff['dif-aa'] = None
                    cod_diff['pos'] = i

                    diff['cod-info'].append(cod_diff)
                    cod_diff = {}
                i = i + len(dif) - dif.count('-') # Do not consider indels

            if i%3 == 0:
                diff['still-translatable'] = True
            else:
                diff['still-translatable'] = False


    ### STEP 4 : Obtain, for each sequence and for each CDS, all the new codons
    new_cds_by_seq = dict()
    for seq_id in all_seqs_id:
        # Contains, for each sequence, all its CDSs
        new_cds_by_seq[seq_id] = dict()
        for gene_id in genes.keys():
            # Contains a list of codons
            new_cds_by_seq[seq_id][gene_id] = list()

    # Initialize all codons with the same values as ref
    for seq_id in new_cds_by_seq.keys():
        for gene_id in diff_by_gene_relative.keys():
            new_cds_by_seq[seq_id][gene_id] = codons[gene_id]

    # Replace codons involved in differences
    for gene_id in diff_by_gene_relative.keys():
        for diff in diff_by_gene_relative[gene_id]:
            for seq_id in diff['where']:
                for i in range(0, len(codons)):
                    if codons[gene_id][i] != diff['codons'][i]:
                        new_cds_by_seq[seq_id][gene_id][i] = diff['codons'][i]

    # Remove list of codons, useless after this
    for gene_id in diff_by_gene_relative.keys():
        for diff in diff_by_gene_relative[gene_id]:
            diff.pop('codons')
    
    return (diff_by_gene_relative, new_cds_by_seq)