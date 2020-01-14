# PhEMD - "Phenotypic Earth Mover's Distance"
PhEMD is a package for relating a network of single-cell samples. By first modeling the cell-state landscape of a model system and then computing pairwise distances between samples using EMD, PhEMD constructs a low-dimensional embedding that uncovers intrinsic manifold structure.


## Installing PhEMD
PhEMD can currently be installed directly through Bioconductor as follows:
```
BiocManager::install("phemd")
```
Alternatively, the package can be cloned and installed directly from Git as below (although installation via Bioconductor is preferred).
```
git clone --recursive git://github.com/wschen/phemd.git
cd phemd
R CMD INSTALL .
```

## Running PhEMD
See https://bioconductor.org/packages/release/bioc/vignettes/phemd/inst/doc/phemd.html for a step-by-step tutorial on applying PhEMD to a multi-sample single-cell dataset.
