# Automated Parameter Assignment for AMOEBA+ Model

## Introduction
This program (`lAssignAMOEBAplusPRM.py`), together with others, is designed to automatically assign the bonded and nonbonded parameters for the AMOEBA and AMOEBA+ models, based on chemical pattern matching. Here we will use phenol as an example (files can be found in example folder)

### Bonded terms

It can automatically assign valence parameters to molecules based on current atomic types and parameter set. 

```shell
python lAssignAMOEBAplusPRM.py -potent bonded -xyz phenol.xyz -key phenol.key -konly NO
```
In addition, ranking tree and reversed searching algorithm can be used to assign the best bonded parameters if no exact match available.
* Format of the typing tree

  The whole structure of typing tree has been documented in file typing_tree.log.
  The format for this file:
  
  level   upper_index   general_index   type_index
  
### Polarizability

Polarizability parameters have been derived for neutral molecules and common charged molecules. To assign this set of parameters, run the following command:

```shell
python lAssignAMOEBAplusPRM.py -potent polar -xyz phenol.xyz -key phenol.key
```

### Charge flux

Charge flux parameters can be assigned if `bond` and `angle` keywords exist in the `key` file. 

```shell
python lAssignAMOEBAplusPRM.py -potent CF -xyz phenol.xyz -key phenol.key
```

### Reference
* Polarizability parameters. In submission (2021)
* Valence parameters. Submitted (2021)
* CF parameters. J. Chem. Phys. 153, 064103 (2020); [link](https://doi.org/10.1063/5.0016376)
