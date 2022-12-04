# Software for nonnegative block-term decomposition with the beta-divergence:
This repository contains the software for joint HSR and unmixing using LL1 decomposition.

Copyright (c) 2022 Clemence Prevost, Valentin Leplat <br>
Contact: ```clemence.prevost@univ-lille.fr``` or ```v.leplat@skoltech.ru```


This MATLAB software reproduces the results from the following paper:

```
@unpublished{prevost:hal-03831661,
  TITLE = {{NONNEGATIVE BLOCK-TERM DECOMPOSITION WITH THE BETA-DIVERGENCE: JOINT DATA FUSION AND BLIND SPECTRAL UNMIXING}},
  AUTHOR = {Pr{\'e}vost, Cl{\'e}mence and Leplat, Valentin},
  URL = {https://hal.archives-ouvertes.fr/hal-03831661},
  NOTE = {working paper or preprint},
  YEAR = {2022},
  MONTH = Oct,
  KEYWORDS = {Nonnegative tensor factorization ; block-term decomposition ; beta-divergence ; blind spectral unmixing ; hyperspectral super-resolution},
  PDF = {https://hal.archives-ouvertes.fr/hal-03831661/file/Paper_LongVersion.pdf},
  HAL_ID = {hal-03831661},
  HAL_VERSION = {v1},
}
```

## Acknowledgements

The baseline algorithms used in the manuscript are courtesy of their respective authors.


## Content
 
 - /data : contains synthetic datasets and reference endmembers for the Samson and Jasper Ridge dataset.
 
 - /demos : contains demo files that produce tables and figures (including ```demo.m```).

 - main.m : codes that allow to run desired demos.
 
 - /src : contains helpful files to run the demos.

## Minimal requirements

In order to run the demos, you will need to add to your MATLAB path:
- Tensorlab 3.0: https://www.tensorlab.net

Additionnally, you will need:
- The HSMS data fusion toolbox: 
   - either here:  http://bitly.ws/xt2b
   - or here    :  https://openremotesensing.net/knowledgebase/hyperspectral-and-multispectral-data-fusion/ 
- The tensor HSR toolbox: https://drive.google.com/file/d/1dwLwwNTseMGIGhxS432BC3PREkJSMlo7/view?usp=share_link

Tips:
- Create a subdirectory "libraries",
- place the libraries "Tensorlab 3.0", "HSMS data fusion toolbox" and "tensor HSR toolbox" within,
- launch script "add_path_to_libraries.m".


## Demo file
 
 A demo with minimal requirements is available. To proceed, please type "1" when running the ```main.m``` file.
 
 ## Tuning the parameters
 
 There are several parameters that you can choose:
 - The noise type: Gaussian, Poisson or Gamma
 - The ranks Lr of the decomposition
 - The number of trials
 - The number of iterations and the threshold for the proposed approach.
 
For benchmarked approaches, the parameters have been tuned according to the original works.
 
  
  ## Reproduce figures from the paper
  
  To do so, you need to run the ```main.m``` file. Here, a menu is available and allows you to choose which figure or table you want to generate. Each number in the table below corresponds to a set of figures.

| Number | Content                                                         |
|--------|-----------------------------------------------------------------|
| 1      | Demo with minimal requirements                                  |
| 2      | Benchmarking on synthetic dataset                               |
| 3      | Benchmarking on semi-real dataset                               |
| 4      | Benchmarking on semi-real dataset with unknown operators        |

