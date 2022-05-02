# Method failure notes

Notes on when methods fail and how.

## TOAST

After original EWAS - many of the "failed" EWAS were when using the TOAST method. All of these had the same error:

`simpleError in solve.default(t(W) %*% W): system is computationally singular: reciprocal condition number = XXX`, where XXX = some very small number, e.g. 4.87335e-26. 

In order to run TOAST, you need to setup your phenotype, methylation, and cell type data and run the following commands:

```r
Design_out <- makeDesign(phenotype_data, cell_count_data)
fitted_model <- fitModel(Design_out, methylation_data)
res <- csTest(fitted_model, coef = trait, 
              cell_type = NULL, contrast_matrix = NULL)
```

The error occurs after running the `fitModel()` function. The error indicates the design matrix is not invertible and therefore can't be used to develop a regression model. This results from linearly dependent columns, i.e. strongly correlated variables. After checking the correlation of covariates when running the analysis for one "failed" trait, it wasn't too high: max = 0.21 (absolute value).

One potential option that might help here, is to use fewer covariates in the analysis. Currently, the covariates are: age, top 20 SVs, top 10 genomic PCs.

## Other methods

No failures for TCAREG. The other methods seem to be suffering from a similar problem to TOAST when they break:

CellDMC:

`simpleError in CellDMC(beta.m = temp_meth, pheno.v = temp_phen[[phen]], frac.m = cc,     cov.mod = cov_mat): The design matrix is not full ranked.
This means that you coundn't make inference for all cell-types in your fraction matrix.
This is usally casued by fractions of a cell-type of one pheno type are all 0 or some fractions in one pheno type are paralle to that of another cell-type.
You might use which(colSums(model.matrix(~ frac.m + pheno.v:frac.m)[, -1]) == 0) to find the cell type.>`


TCA:

`simpleError in quadprog::solve.QP(Dmat, dvec, t(Amat), bvec, meq = meq): matrix D in quadratic function is not positive definite!`
