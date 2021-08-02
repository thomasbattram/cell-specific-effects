# Notes on each method as reading through the papers

__NOTE: THESE NOTES WERE MADE AFTER A QUICK READ THROUGH THE PAPERS__.

The maths behind some of the methods would take me a while to understand properly so this is not really discussed here. Take what you will from that.

## [CellDMC](https://pubmed.ncbi.nlm.nih.gov/30504870/)

The input of the method is just the DNAm matrix, covariates, and a cell type reference of the users choosing and the output is an effect estimate for each CpG in each cell type and which cell type seems to be driving the associations between DNAm and the phenotype of interest. 

The software uses a simple linear model with an interaction term (interaction between cell type and phenotype) to estimate these effects.

Simulations suggested the cell type fraction range had a large impact on sensitivity. Each person has a different fraction of each cell type and the cell type fraction range is simply the range of cell type fractions across all individuals. Low cell type fraction ranges gave very low sensitivities (see fig 2 of the paper), so for things like eosinophils CellDMC is unlikely to be able to identify cell-specific effects UNLESS the number of eosinophils in cases and controls differ substantially.

The paper uses 3 datasets to validate the method, one of which is the shite Liu et al. 2013 rheumatoid arthritis EWAS dataset on GEO so results of that analysis don't really say too much about CellDMC.

## [TCA](https://pubmed.ncbi.nlm.nih.gov/31366909/)

The input is the same as for CellDMC - DNAm matrix, covaraites, a cell type reference. The output also looks fairly similar. Currently it's a little less clear and will need to run it to check.

TCA is more complex than CellDMC. It creates something called "tensors" which represent the 3D data structure of cell-types by DNAm by individuals. It allows the users to choose whether DNAm is the outcome or exposure and can fit models accordingly. 

Simulations suggested increased power compared to CellDMC and empirical analyses suggested increased specificity. However, these results aren't too reliable (see below). 

TCA also provides an option to estimate or re-estimate cell types. However, Matt's student tested these re-estimated cell-types and they don't seem to be capturing real cell-types at all ([report.docx](../previous-work/report.docx)). This was also suggested by another [paper](https://www.biorxiv.org/content/10.1101/2021.02.28.433245v1.full.pdf) (figure 1).

## TCA vs. CellDMC

There has been lots of back and forth between these groups and this is how the timeline of these papers has gone:

1. CellDMC was proposed as the first (?) method to ascertain cell-specific associations from bulk tissue data for EWAS: [paper link](https://www.nature.com/articles/s41592-018-0213-x)
2. [TCA paper](https://www.nature.com/articles/s41467-019-11052-9) suggested that TCA has an improved performance over CellDMC
3. Rebuttal by CellDMC group ([paper link](https://www.biorxiv.org/content/10.1101/822940v2)) saying the default model of TCA shown in their paper was different to model used by CellDMC so comparisons were meaningless. CellDMC fits a X|Y model (i.e. DNAm is outcome) and the default TCA fits a Y|X model (i.e. DNAm is exposure). Their simulations also suggested that the positive predictive value (PPV), which is equivalent to 1-FDR, is fairly low when using the TCA Y|X model.
4. Rebuttal by TCA group ([paper link](https://www.biorxiv.org/content/10.1101/2021.02.14.431168v1.full.pdf)) saying comparison was fine and gave new results comparing TCA and CellDMC using an X|Y model. HOWEVER their new results suggested that actually TCA and CellDMC perform very similarly when using this model (figure 1). They also compare methods again when the TRUE (simulated) model is a Y|X model (and so TCA is setup this way) and this suggests TCA performs slightly better in general, but this depends on the number of cell types for which there is an effect at a given CpG site (figure 2).
5. Last rebuttal by CellDMC group ([paper link](https://www.biorxiv.org/content/10.1101/2021.02.28.433245v1.full.pdf)) essentially saying that the TCA Y|X model is shite (as argued in 3), and the TCA X|Y model is better, but has no real improvement over CellDMC. They also levy other criticisms at the TCA paper, for example the exact model that was used to demonstrate the performance of TCA differed depending on the simulation scenario. This is ofc a problem because the real scenario is unknown.

## [HIRE](https://pubmed.ncbi.nlm.nih.gov/31308366/)

Input is DNAm matrix and covariates. The user also needs to specify the number of cell types to estimate. It doesn't seem as though the user can specify the cell types before hand AND the results from the work of Matt's student ([report.docx](../previous-work/report.docx)) suggests the cells estimated by HIRE don't tend to recapitulate real (FACS sorted) cell-types

The output is the cell types (I think?) and results of the EWAS.

In the paper they only compare their results with EWAS that also use reference-free approaches of estimating cell-comp (e.g. SVA and RefFreeEWAS). They don't compare their cell-specific effects with those from CellDMC or TCA. The method for which they obtain these effects isn't entirely clear to me either, but that doesn't really matter because the cell types generated can't be trusted

## [omicwas](https://pubmed.ncbi.nlm.nih.gov/33752591/)

The input is the same as for CellDMC - DNAm matrix, covaraites, a cell type reference. The output is very similar too - coefficients and P values for each cell type.

There are several different options to test the association between DNAm within each cell type and the trait of interest. Within the paper they suggest using nonlinear ridge regression for DNAm. 

In the paper they suggest that running linear regressions and adding in interaction terms to estimate cell-specific effects induces multicollinearity (figure 1 in paper). This is the model that CellDMC and TOAST uses, so it's suggesting these models are incorrect.

Simulations involved using a publicly available DNAm dataset and simulating phenotypes for which a small percentage of the CpG sites were associated with. When comparing sensitivity, specificity and positive predictive value (PPV), the methods performed similarly. Each had very high specificity, generally low sensitivity, and varied PPV across cell types. It's unclear how they fitted the models, including whether they treated DNAm as the exposure or outcome - which is especially important for TCA.

## [TOAST](https://pubmed.ncbi.nlm.nih.gov/31484546/)

The input is the same as CellDMC - DNAm matrix, covariates, a cell type reference. 

TOAST seems to be first and foremost an algorithm to generate cell types without a reference. The paper focuses entirely on the ability of TOAST to do this and does not detail how their model of estimating cell-specific effects performs. Uses fairly simple linear model to estimate cell-specific effects and asks whether there is evidence that the effect in one cell type is different to the mean effect across all other cell types




