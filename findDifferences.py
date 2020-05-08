from collections import OrderedDict  # Remembers the order entries were added
import re
import os


def findDifferences(reference, alignments):
    """
        Find the differences between a reference and one or more alignments. A difference is a dict having keys: 'start','length','type','seq','ref'

        :param reference:
            A dict of the form {ref_id : reference}

        :param alignments:
            A dict of the form {align_id_1 : value, align_id_2 : value, etc.}

        :returns:
            A dict containing, for each alignment, all the differences with the reference. Note that all differences obtained will be of length 1.
    """
    ref_id = list(reference.keys())[0]
    ref = reference[ref_id]
    differences = OrderedDict()
    for align_id in alignments.keys():
        differences[align_id] = []
    for i in range(0, len(ref)):
        for align_id in alignments.keys():
            curr_base = alignments[align_id][i]
            curr_diff = dict()
            if curr_base != ref[i]:  # A difference has been found
                # Store the value of reference
                curr_diff['ref'] = '*'
                # Store the value of the alignment
                curr_diff['seq'] = '*'
                # Store the position where it happens
                curr_diff['start'] = i

                # Determine the difference type
                if ref[i] == 'N':
                    curr_diff['type'] = "na1"  # ref not available
                elif curr_base == 'N':
                    curr_diff['type'] = "na2"  # seq not available
                elif ref[i] == 'N' and curr_base == 'N':
                    curr_diff['type'] = "na3"  # both ref and seq not available
                elif ref[i] == '-' and curr_base != '-':
                    curr_diff['type'] = 'ins'  # insertion
                elif ref[i] != '-' and curr_base == '-':
                    curr_diff['type'] = 'del'  # deletion
                else:
                    curr_diff['type'] = 'rep'  # replacement

                if curr_diff['type'] == 'na1':
                    # only store seq since ref is not available
                    curr_diff['seq'] = curr_base
                elif curr_diff['type'] == 'na2':
                    # only store ref since seq is not available
                    curr_diff['ref'] = ref[i]
                # elif curr_diff['type'] == 'na3': # do nothing
                elif curr_diff['type'] == 'ins':
                    # store insert in seq
                    curr_diff['seq'] = curr_base
                elif curr_diff['type'] == 'del':
                    # store deleted in ref
                    curr_diff['ref'] = ref[i]
                else:  # base changed
                    # store both
                    curr_diff['seq'] = curr_base
                    curr_diff['ref'] = ref[i]

                # Add the newfound difference to the other ones
                differences[align_id].append(curr_diff)
    return differences


def compactDifferences1(differences):
    """
        Compact for each alignment all the consecutive differences of the same type

        :param differences:
            The differences dict obtained from findDifferences()

        :returns:
            A dict containing the "compacted" differences, all having length >= 1
    """
    new_Differences = OrderedDict()
    for align_id in differences.keys():
        new_Differences[align_id] = []

    for align_id in new_Differences.keys():
        diff = differences[align_id]

        curr_diff = dict()
        if diff:  # Not empty
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

                else:  # New difference
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

            if curr_diff:  # Last element still in curr_diff
                new_Differences[align_id].append(curr_diff)
                curr_diff = {}

    return new_Differences

def compactDifferences2(compactDifferences):
    index = dict()
    values = dict()
    for seq_id in compactDifferences.keys():
        # Check if there is at least one difference
        if compactDifferences[seq_id]:
            index[seq_id] = 0
            values[seq_id] = ''

    result = []
    empty = False
    while not empty:
        # Extract values for comparison in this iteration
        for seq_id in values.keys():
            curr_index = index[seq_id]
            values[seq_id] = compactDifferences[seq_id][curr_index]
        
        # Find the difference with the lowest starting position
        lowest_starting = dict()
        for difference in values.values():
            if not lowest_starting:
                lowest_starting = difference
            elif difference['start'] < lowest_starting['start']:
                lowest_starting = difference
            # else do nothing
        
        # Find if difference is contained in multiple alignments
        # NOTE: differences are sorted by starting pos
        where_field = list()
        for seq_id in values.keys():
            if lowest_starting == values[seq_id]:
                # store the seq_id
                where_field.append(seq_id)
                # increase the index
                index[seq_id] = index[seq_id]+1

        lowest_starting['where'] = where_field
        where_field = []

        result.append(lowest_starting)
        lowest_starting = {}

        # Check if all lists have been empty'd
        empty=True
        for seq_id in list(index.keys()): 
            # List necessary due to RuntimeError: dictionary changed size during iteration
            number_of_differences = len(compactDifferences[seq_id])
            if index[seq_id] < number_of_differences:
                empty=False
            else:  # all differences have been analyzed, remove from dicts
                index.pop(seq_id)
                values.pop(seq_id)

    return result

    def optmizeDifferences(differences):
        result = list()
        
        differences = []
        if differences:
            for diff1, diff2 in zip(differences, differences[1:]):
                if diff1['start'] == diff2['start']:
                    curr_differences.append(diff1)
                    pass
        
        return result


def findStats(reference, alignments):
    """
        For each alignment, find its stats (no. of Matches, no. of Mismatches, no. of NAs)

        :param reference:
            A dict of the form {ref_id : reference}

        :param alignments:
            A dict of the form {align_id_1 : value, align_id_2 : value, etc.}

        :returns:
            A dict containing, for each alignment, all of its stats
    """
    # Initialize dict
    stats = dict()
    for align_id in alignments.keys():
        stats[align_id] = dict()
        stats[align_id]['matches'] = 0
        stats[align_id]['mismatches'] = 0
        stats[align_id]['na'] = 0

    # Extract reference from dict
    ref_id = list(reference.keys())[0]
    ref = reference[ref_id]

    # Add entry in stats for reference
    stats[ref_id] = dict()
    stats[ref_id]['na'] = 0  # Count number of Ns

    for i in range(0, len(ref)):
        # Check if ref is Not Available
        if ref[i] == 'N':
            stats[ref_id]['na'] = stats[ref_id]['na'] + 1
        # Check each alignment
        for align_id in alignments.keys():
            curr_base = alignments[align_id][i]
            if curr_base == 'N':
                stats[align_id]['na'] = stats[align_id]['na'] + 1
            elif curr_base == ref[i]:
                stats[align_id]['matches'] = stats[align_id]['matches'] + 1
            else:  # curr_base != ref[i]
                stats[align_id]['mismatches'] = stats[align_id]['mismatches'] + 1

    return stats


def writeToFile(fileName, reference, compactDifferences, stats):
    """
        Print all differences to terminal
        :param fileName:
            The name of the file that will be used to store the differences (WITHOUT the extension)
        :param reference:
            A dict of the form {ref_id : reference}
        :param compactDifferences:
            The compactDifferences obtained from compactDifferences()
        :param stats:
            The stats obtained from findStats()
    """
    ref_id = list(reference.keys())[0]
    ref = "##Ref={},len={},NA={}\n".format(
        ref_id, len(reference[ref_id]), stats[ref_id]['na'])
    fields = "START\tLENGTH\tTYPE\tREF\tSEQ\t\n"
    value = "{}\t{}\t{}\t{}\t{}\t\n"
    path = "{}.mad".format(fileName)  # MAD = Multiple Alignment Difference
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
                temp = value.format(
                    diff['start']+1, diff['length'], diff['type'], diff['ref'], diff['seq'])
                f.write(temp)

def writeToFile2(fileName, reference, differences, stats):
    """
        Print all differences to terminal
        :param fileName:
            The name of the file that will be used to store the differences (WITHOUT the extension)
        :param reference:
            A dict of the form {ref_id : reference}
        :param differences:
            A list of differences with the 'where' field from compactDifferences2()
        :param stats:
            The stats obtained from findStats()
    """
    ref_id = list(reference.keys())[0]
    ref = "##Ref={},len={},NA={}\t\n".format(
        ref_id, len(reference[ref_id]), stats[ref_id]['na'])
    fields = "START\tLENGTH\tTYPE\tREF\tSEQ\tWHERE\t\n"
    value = "{}\t{}\t{}\t{}\t{}\t{}\t\n"
    path = "{}.mad2".format(fileName)  # MAD2 = Multiple Alignment Difference 2
    with open(path, "w") as f:
        f.write(ref)
        stats.pop(ref_id) # remove stats relative to ref
        for align_id in stats.keys():
            seq = "##Seq={},Matches={},Mismatches={},NA={}\n".format(align_id,
                                                                    stats[align_id]['matches'],
                                                                    stats[align_id]['mismatches'],
                                                                    stats[align_id]['na'])
            f.write(seq)
        
        f.write(fields)
        for diff in differences:
            temp = value.format(
                diff['start']+1, diff['length'], diff['type'], diff['ref'], diff['seq'],diff['where'])
            f.write(temp)


def findCommonDifferences(differencesByAligner):
    """
        Store all common differences
        :param differencesByAligner:
            A dict of the form {aligner_1 : differences, aligner_2 : differences}
    """
    common_diff = list()

    for aligner in differencesByAligner:
        for seq_id in aligner.keys():
            if not common_diff[seq_id]:
                common_diff[seq_id] = set()
            common_diff[seq_id] = common_diff[seq_id] & set(aligner[seq_id])
    return common_diff
            

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


def parseClustal(referenceId, fileName):
    """
        Parse Clustal file

        :param referenceId:
            The id of the reference as of the .clw file
        :param fileName:
            The name of the file that contains the alignments (WITHOUT the extension)
        :returns:
            A dict containing the reference, and a dict containing the alignments
    """
    alignments = dict()
    reference = dict()
    path = "{}.clw".format(fileName)
    with open(path, "r") as f:

        # Remove header
        lines = f.readlines()[1:]
        for line in lines:

            # remove (as many as possible) spaces, *, \n from beginning and end of string
            line = line.strip(" *\n")
            # replace tabs with spaces for split
            line = line.replace('\t', ' ')
            # replace multiple spaces with single space for split
            line = re.sub(' +', ' ', line)

            if line != '':
                line_to_list = line.split()
                key = line_to_list[0]
                align = line_to_list[1]
                # if len(line_to_list) == 3:  # Note: only some aligners show base count
                #    num_of_bases = line_to_list[2]
                if key not in alignments.keys():  # Create dict entry
                    alignments[key] = ''
                alignments[key] = alignments[key] + align

        ref = alignments.pop(referenceId)
        reference[referenceId] = ref

        return reference, alignments


def main():
    path_alignments = "./Alignments/"
    path_output = "./Outputs/"
    all_differences_by_aligner = dict()
    for file in os.listdir(path_alignments):
        if file.endswith(".clw"):  # for each .clw file
            # Obtain filename, will be used for both input and output
            fileName = os.path.splitext(file)[0]
            # Parse reference and alignments from file
            reference, alignments = parseClustal(
                "NC_045512", path_alignments+fileName)
            # Find for each alignment all single differences
            diff = findDifferences(reference, alignments)
            # Compact consecutive difference
            compact = compactDifferences1(diff)
            # Find alignments stats
            stats = findStats(reference, alignments)
            # Write differences to file
            writeToFile(path_output+fileName, reference, compact, stats)

            # START of part 2
            # Compact differences by column
            differences = compactDifferences2(compact)
            # Print differences to file
            writeToFile2(path_output+fileName, reference, differences, stats)
            
            # Find differences by gene
            dif_by_gene = findDifferencesbyGene(compact, {'gene1':[5000,6000]})
            
            #Store differences by aligner
            all_differences_by_aligner[fileName] = dict()
            all_differences_by_aligner[fileName] = compact 
    
    #writeToFile2(path_output+'result',all_differences_by_aligner)

if __name__ == "__main__":
    main()
