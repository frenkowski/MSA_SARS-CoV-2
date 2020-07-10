# Multiple Sequence Alignment (MSA) Analysis of SARS-Cov2 sequences

## Team

Ali Manan (_mohammadalimanan+github@gmail.com_)  
Francesco Porto (_francescoporto97+github@gmail.com_)  
Francesco Stranieri (_frenkowski+github@gmail.com_)

## Abstract

In this report we analyze Multiple Sequence Alignments (MSAs) obtained from 18 samples of European COVID-19 patients and we compare against 3 samples from China, in order to figure out how similar or different they are. We also present our tool for detecting differences in a MSA and building a phylogenetic tree from a subset of these differences. By analyzing the results of our tool on the provided samples, we make an attempt to either confirm or deny various claims about the origin of the outbreak in Europe.

There are 3 parts to this project:
1. We look for differences between a reference sequence and other sequences in a MSA; we obtain a file containing information on these differences (i.e. position, type, etc.)
2. We compare codons in annotated genes to figure out whether CDS are still translatable or not. We also obtain a list of codon-differences, along with their relative amino-acids.
3. We obtain a phylogenetic tree from a character matrix obtained from the output of the previous part.

## Instructions

- Place a MSA in the /Alignments folder.
- Set the reference name in /Src/findDifferencesMSA.py (_line 255_).
- Run the following command:
  `python3 runScript.py`
- The resulting files should be /Outputs/Part_X.

**NOTE:** currently only works under \*NIX.
