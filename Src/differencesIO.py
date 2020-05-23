import re

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
    path = "{}.aln".format(fileName)
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
                if key.endswith((".1",".2")):
                    key = key[:-2]
                if key not in alignments.keys():  # Create dict entry
                    alignments[key] = ''
                alignments[key] = alignments[key] + align

        ref = alignments.pop(referenceId)
        reference[referenceId] = ref

        return reference, alignments

def writeToFilePairwise(fileName, reference, compactDifferences, stats):
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
    fields = "#START\tLENGTH\tTYPE\tREF\tSEQ\t\n"
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
                    diff['start']+1, 
                    diff['length'], 
                    diff['type'], 
                    diff['ref'], 
                    diff['seq'])
                f.write(temp)

def writeToFileMSA(fileName, reference, differences, stats):
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
    ref = "##Ref={},len={},NA={}\n".format(
        ref_id, len(reference[ref_id]), stats[ref_id]['na'])
    fields = "#START\tLENGTH\tTYPE\tREF\tSEQ\tWHERE\n"
    value = "{}\t{}\t{}\t{}\t{}\t{}\n"
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
            where_to_string = str(diff['where']).strip('[]').replace('\'','')
            temp = value.format(
                diff['start']+1, 
                diff['length'], 
                diff['type'], 
                diff['ref'], 
                diff['seq'],
                where_to_string)
            f.write(temp)

def writeToFileFinal(fileName, reference, differences, stats):
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
    ref = "##Ref={},len={}\n".format(
        ref_id, len(reference[ref_id]))
    fields = "#START\tLENGTH\tTYPE\tREF\tSEQ\tWHERE\tTOOLS\n"
    value = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
    path = "{}.mad3".format(fileName)  # MAD3 = Multiple Alignment Difference 3
    with open(path, "w") as f:
        f.write(ref)
        
        #for tool in stats.keys():
            #print(stats[tool][ref_id])
            #stats[tool].pop(ref_id) # remove stats about ref
            #for align_id in stats[tool].keys():
        #print(stats)
        for align_id in stats.keys():
            #stats.pop(ref_id)
            for stat in stats[align_id]:
                seq = "##Tool={},Seq={},Matches={},Mismatches={},NA={}\n".format(
                                                                        stat['aligner'],
                                                                        align_id,
                                                                        stat['matches'],
                                                                        stat['mismatches'],
                                                                        stat['na'])
                f.write(seq)
        
        f.write(fields)
        for diff in differences:
            where_to_string = str(diff['where']).strip('{}').replace('\'','').replace(' ','')
            aligns_to_string = str(diff['aligns']).strip('[]').replace('\'','').replace(' ','')
            temp = value.format(
                diff['start']+1, 
                diff['length'], 
                diff['type'], 
                diff['ref'], 
                diff['seq'],
                where_to_string,
                aligns_to_string
            )
            f.write(temp)