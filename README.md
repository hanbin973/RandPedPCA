# Randomized Pedigree Principal Component Analysis

Randomized Pedigree Principal Component Analysis (RandPedPCA) performs rapid principal component analysis (PCA) of pedigree genetic relatedness matrix (GRM) using randomized linear algebra.
[Henderson (1975)](https://doi.org/10.3168/jds.S0022-0302(75)84776-X) developed an efficient way to compute the lower Cholesky factor of the inverse GRM and
[Colleau (2002)](https://doi.org/10.1186/1297-9686-34-4-409) explicitly showed how to multiply pedigree GRM with an arbitrary vector efficiently.
RandPedPCA uses these two algorithmic ingredients to rapidly compute the principal components of pedigree GRM without forming the pedigree GRM.
This is achieved via the randomized singular value decomposition (rSVD) described in [Halko et al. (2011)](http://dx.doi.org/10.1137/090771806).
The resulting principal components can reveal the underlying population structure of a pedigree.

# R package

## Setup
To install the package from GitHub (it's not yet on CRAN), run:
```
devtools::install_github("HighlanderLab/RandPedPCA", subdir = "randPedPCA",
                         ref = "v0.9.8", build_vignettes = TRUE)
```
## First steps
For a demonstration, check out
```
vignette("pedigree-pca")
```

# Python prototype

Initial prototype was developed in Python. An example can be found in `notebook/Example.ipynb`, which uses `rppca` module in this repository.
