# PhEMD - "Phenotypic Earth Mover's Distance"
PhEMD is a package for relating a network of single-cell samples. By first modeling the cell-state landscape of a model system and then computing pairwise distances between samples using EMD, PhEMD constructs a low-dimensional embedding that uncovers intrinsic manifold structure.


## Installing PhEMD
PhEMD can currently be installed by cloning and installing the package directly from Git. We plan to make PhEMD available through Bioconductor soon.
git clone --recursive git://github.com/wschen/phemd.git
cd phemd
R CMD INSTALL .
