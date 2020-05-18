import os
from differencesIO import parseClustal, writeToFile, writeToFile2, writeToFile3
from findDifferences import findStats

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

    # Get the reference
    ref_id = list(reference.keys())[0]
    ref = reference[ref_id]
    
    # Initialize variables
    differences = list()
    active_differences = list()
    temp_differences = list()
    #for align_id in alignments.keys():
    #    differences[align_id] = []

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
                # Store the length
                curr_diff['length'] = 1
                # Store the sequence_id where it happens
                curr_diff['where'] = set()
                curr_diff['where'].add(align_id)

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

                # Set where field
                found = False
                for diff in temp_differences:
                    if curr_diff['seq'] == diff['seq']:
                        found = True
                        diff['where'] = diff['where'].union(curr_diff['where'])
                        break 
                if not found:
                    temp_differences.append(curr_diff)
                
                curr_diff = {}

        # differences.append(temp_differences)
        # temp_differences = []
        #extend = False
        for diff_to_add in temp_differences:
            extend = False

            # Try to add diff_to_add to active_differences
            for diff in active_differences:
                # Consecutive diff
                if diff_to_add['start'] == diff['start'] + diff['length']:
                    # Extend one of active differences
                    if (diff_to_add['where'] == diff['where'] and 
                        diff_to_add['type'] == diff['type']):
                        # Extend diff length
                        diff['length'] = diff['length'] + 1
                        # Set ref and seq accordingly
                        if diff_to_add['type'] == 'na1':
                            diff['seq'] = diff['seq'] + diff_to_add['seq']
                        elif diff_to_add['type'] == 'na2':
                            diff['ref'] = diff['ref'] + diff_to_add['ref']
                        # elif diff_to_add['type'] == 'na3': # do nothing
                        elif diff_to_add['type'] == 'ins':
                            diff['seq'] = diff['seq'] + diff_to_add['seq']
                        elif diff_to_add['type'] == 'del':
                            diff['ref'] = diff['ref'] + diff_to_add['ref']
                        else:
                            diff['seq'] = diff['seq'] + diff_to_add['seq']
                            diff['ref'] = diff['ref'] + diff_to_add['ref']

                        extend = True
                        temp_differences.remove(diff_to_add)
                        break
                    # else: # New diff that start at next pos
                    #     differences.append(diff)
                    #     # Remove diff to active differences
                    #     active_differences.remove(diff)
                    #     # Add diff_to_add to active differences
                    #     active_differences.append(diff_to_add)
            
            if not extend:
                active_differences.append(diff_to_add)
                temp_differences.remove(diff_to_add)

        # Remove inactive differences
        for diff in active_differences:
            # Can't extend difference since there is at least 1 mismatch
            if (diff['start'] + diff['length']) < i:
                differences.append(diff)
                active_differences.remove(diff)

        #if temp_differences != []:
        #    print(temp_differences)
        # Reset temp_differences for next iteration
        #temp_differences = []

    if active_differences:
        differences = differences + active_differences

    return differences

def compactDifferencesByAligner(all_differences_by_aligner):
    result = list()
    aligns_field = list()

    for aligner in all_differences_by_aligner.keys():
        for diff in all_differences_by_aligner[aligner]:
            # Check if diff is present in other alignments'diff
            for other_aligner in (all_differences_by_aligner.keys()-aligner):
                # Obtain the list that contains all its diffs
                other_aligner_diffs = all_differences_by_aligner[other_aligner]
                if diff in other_aligner_diffs:
                    # remove diff (since it's duplicated)
                    other_aligner_diffs.remove(diff)
                    # add diff to aligns
                    aligns_field.append(other_aligner)
            
            diff['aligns'] = aligns_field
            aligns_field = []

            result.append(diff)
            
    ## TODO: sort by starting pos?
    # result = result.sort(reverse=False, key=getStartingPos())
    result.sort(key=lambda x: x['start'])
    #print(result)
    return result

def main():
    #differences = findDifferences({'3':'a'},{'1':'a','2':'a'})
    #differences = findDifferences({'3':'a'},{'1':'b','2':'b'})
    #differences = findDifferences({'3':'aa'},{'1':'bb','2':'bb'})
    #differences = findDifferences({'3':'aa'},{'1':'ba','2':'bb'})
    #differences = findDifferences({'3':'aa'},{'1':'ba','2':'ba','3':'ca'})
    #differences = findDifferences({'3':'aaaa'},{'1':'aaab','2':'aaab'})
    #differences = findDifferences({'3':'aaaa'},{'1':'bbba','2':'bbba'})
    #differences = findDifferences({'3':'aaaa'},{'1':'abbb','2':'abbb'})
    #print(differences)

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
            # Find the differences
            differences = findDifferences(reference, alignments)
            # Find alignments stats
            stats = findStats(reference, alignments)
            # Write differences to file
            writeToFile2(path_output+fileName, reference, differences, stats)
            #Store differences by aligner
            all_differences_by_aligner[fileName] = dict()
            all_differences_by_aligner[fileName] = differences
    
    result = compactDifferencesByAligner(all_differences_by_aligner)
    writeToFile3(path_output+'finalResult',reference,result)

if __name__ == "__main__":
    main()