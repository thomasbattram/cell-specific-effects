# cell-specific effects

This repo is to store notes and potentially work on the effectiveness of methods that look to ascertain cell-specific DNAm effects in EWAS.

[`methods-list.xlsx`](methods-list.xlsx) contains a list of the different methods with their PMID, paper title, first author, corresponding author + email, and how the method can be implemented.

Notes on the different methods from their respective papers can be found in the [`notes`](notes) folder.

[`previous-work`](previous-work) contains work done by a summer student of Matt's that assessed how well TCA and HIRE tend to estimate actual cell-types. 

The false positive rate of the methods can be tested by running [`test-methods.R`](scripts/test-methods.R) - this requires submission to bc3. To plot the results run through [`test-methods-results.R`](scripts/test-methods-results.R).

To do:

0. Streamline processes using snakemake (maybe)
1. Test the methods on real phenotypes to see if they tend to output the same results
2. Generate a report