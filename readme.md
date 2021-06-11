# cell-specific effects

This repo is to store notes and potentially work on the effectiveness of methods that look to ascertain cell-specific DNAm effects in EWAS.

[`methods-list.xlsx`](methods-list.xlsx) contains a list of the different methods with their PMID, paper title, first author, corresponding author + email, and how the method can be implemented.

Notes on the different methods from their respective papers can be found in the [`notes`](notes) folder.

[`previous-work`](previous-work) contains work done by a summer student of Matt's that assessed how well TCA and HIRE tend to estimate actual cell-types. 

The false positive rate of TCA and CellDMC can be tested by running [`test-methods.R`](scripts/test-methods.R) - currently this is just when both models are assuming DNAm is the outcome. To do this a submission script will need to be made for bluecrystal.