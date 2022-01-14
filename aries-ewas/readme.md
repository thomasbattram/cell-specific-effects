# ARIES EWAS

This section of the project looks to run cell-specific EWAS for real phenotypes within ARIES and then compare the results for these EWAS with those that were conducted as part of The EWAS Catalog.

Broad steps:

1. Extract ALSPAC data
2. Run cell-specific EWAS
3. Compare cell-specific results with normal EWAS

## Extract ALSPAC data

1. Move phenotype meta data from The EWAS Catalog RDSF to `data/ec-aries-metadata.tsv`
2. Extract ALSPAC phenotype data using [`pheno-data-extraction.R`](scripts/pheno-data-extraction.R)
3. Move that data to the RDSF and to BluePebble

__NOTE:__ Make sure no ALSPAC data is saved locally!