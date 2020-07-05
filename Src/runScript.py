import os
from findDifferencesMSA import findDifferences, compactDifferencesByAligner
from differencesIO import parseFasta, parseClustal, writeToFileMSA, writeToFileFinal, writeCdsDifferencesToFile
from findDifferencesPairwise import findStats
from geneDifferences import findDifferencesbyGene2, transcribeSequence, findDifferencesRelativePos
from geneDifferences import findTranscriptDifferences, splitSequenceByCds
from findPhylogeny import createMatrix, containsForbidden, createTree, printTree, findBiggestNotForbidden
import numpy as np

def main():
    path_alignments = "../Alignments/"
    path_output = "../Outputs/"
    all_differences_by_aligner = dict()
    #all_stats_by_aligner = dict()
    all_stats_by_id = dict()
    for file in os.listdir(path_alignments):
        if file.endswith(".aln"):  # for each .clw file
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
            writeToFileMSA(path_output+"Part1/"+fileName, reference, differences, stats)
            #Store differences by aligner
            all_differences_by_aligner[fileName] = dict()
            all_differences_by_aligner[fileName] = differences
            #Store stats by aligner
            #all_stats_by_aligner[fileName] = dict()
            #all_stats_by_aligner[fileName] = stats
            #Store stats by align_id
            for align_id in alignments.keys():
                if align_id not in all_stats_by_id.keys():
                    all_stats_by_id[align_id] = list()
                # Store aligner name
                stats[align_id]['aligner'] = fileName
                all_stats_by_id[align_id].append(stats[align_id])
    
    result = compactDifferencesByAligner(all_differences_by_aligner)
    writeToFileFinal(path_output+"Part1/"+'finalResult',reference,result, all_stats_by_id)

    #### START OF PART 2 ####
    genes_NC_045512 = {"ORF1ab":[266,21555], "ORF3a":[25393,26220], "E":[26245,26472], "M":[26523,27191], 
         "ORF6":[27202,27387], "ORF7a":[27394,27759], "ORF7b":[27756,27887], "ORF8":[27894,28259],
         "N":[28274,29533], "ORF10":[29558,29674]}
    
    #reference = parseFasta("../Reference/reference")
    #reference = reference["NC_045512"]

    # Replace T->U in reference
    ref_trs = transcribeSequence(reference["NC_045512"]) # trs = transcribed

    # Obtain CDSs as substrings of reference
    ref_by_cds = splitSequenceByCds(ref_trs, genes_NC_045512)

    # Split the differences by gene position
    diff_by_gene = findDifferencesbyGene2(result, genes_NC_045512)

    # Replace start with start relative (according to genes/cds)
    diff_by_gene_relative = findDifferencesRelativePos(diff_by_gene, genes_NC_045512)

    # Detect differences in CDS
    (cds_differences, new_cds_by_seq) = findTranscriptDifferences(ref_by_cds, diff_by_gene_relative, genes_NC_045512, list(alignments.keys()))

    #print("New cds by seq")
    #print(new_cds_by_seq['MT459899'])

    # Write to file differences 
    writeCdsDifferencesToFile(path_output+'Part2/'+'codon_differences',cds_differences,new_cds_by_seq,genes_NC_045512)


    #### START OF PART 3 ####

    # Create matrix from codon differences
    (sequences_list, matrix) = createMatrix(reference, alignments, diff_by_gene_relative)
    print("A matrix of size {} was found".format(matrix.shape))

    # Find biggest possible matrix not containing forbidden
    matrix = findBiggestNotForbidden(matrix)
    print("The biggest matrix not containing the forbidden matrix has size {}".format(matrix.shape))
    
    # Obtains and prints the tree
    if not containsForbidden(matrix):   # Should always be True, just in case
        tree = createTree(sequences_list,matrix)
        printTree(tree)

if __name__ == "__main__":
    main()