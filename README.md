# Automated Parameter Assignment for AMOEBA(+) Model

## Introduction
This program is designed to automatically assign the bonded and nonbonded parameters for AMOEBA and AMOEBA+ models. 

### Bonded terms

It can automatically assign valence parameters to molecules based on current atomic types and parameter set. In addition, ranking tree and reversed searching algorithm have been implemented.
* Format of the typing tree

  The whole structure of typing tree has been documented in file typing_tree.log.
  The format for this file:
  
  level   upper_index   general_index   type_index
  
### Polarizability

Polarizability parameters have been derived for neutral molecules and common charged molecules. To assign this set of parameters, run the following command:

```shell
python lAssignAMOEBAplusPRM.py -potent polar -xyz your.xyz -key your.prm
```

### Other terms

Coming soon.
