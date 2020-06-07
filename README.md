# Multiple Sequence Alignment (MSA) of SARS-Cov2 sequences

## Team

Ali Manan (_mohammadalimanan+github@gmail.com_)  
Francesco Porto (_francescoporto97+github@gmail.com_)  
Francesco Stranieri (_frenkowski+github@gmail.com_)

## Abstract

In this report we analyze Multiple Sequence Alignments (MSAs) obtained from 18 samples of European COVID-19 patients and we compare against 3 samples from China, in order to figure out how similar or different they are. We also present our tool for detecting differences in a MSA and three ad-hoc formats for storing them. By analyzing the results of our tool on the provided samples, we make an attempt to either confirm or deny various claims about the origin of the outbreak in Europe.

## Instruction

- Place a MSA in the /Alignments folder.
- Set the reference name in /Src/findDifferencesMSA.py (_line 255_).
- Run the following command:
  `python3 findDifferencesMSA.py`
- The resulting MAD3 file should be /Outputs/finalResult.mad3.

**NOTE:** currently only works under \*NIX.
