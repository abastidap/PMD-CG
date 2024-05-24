# Project Title

Bash script to generate structures of the AHSSHLKSKKGQSTSRHKKL peptide using the probabilistic Molecular Dynamics chain growth (PMD-CG) method


### Dependencies

* OS: Ubuntu 20.04.6 LTS
* chimera version 1.16
* pymol version 2.3.0 

### Installing

* download the tripeptides folder

This folder includes the statistical information of the tripeptides

* download the programs folder

This folder includes three fortran codes that should be compiled. An example of compilation (compile.sh) is included using gfortran but any other fortran compiler may be used.

* download the pmd-cg-v1.sh script

### Executing program

* ./pmd-cg-v1.sh first last path1 path2

"first" and "last" are the first and last number of the structures to be created. 

"path1" is the full path of the local tripeptides directory

"path2" is the full path of the local programs directory

* Example

./pmd-cg-v1.sh 1 100 $HOME/tripeptides $HOME/programs

## Authors

Contributors names and contact info

* Adolfo Bastida (bastida@um.es)

## Version History

* v1 

Initial Release

