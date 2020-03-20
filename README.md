# Local influence for the analysis of agreement

[![R](https://img.shields.io/badge/Made%20with-R%20under%20development-success)](https://cran.r-project.org/)
[![heavy](https://img.shields.io/badge/Using%20heavy-0.38.19-important)](https://cran.r-project.org/package=heavy)
[![DOI](https://img.shields.io/badge/DOI-10.1002/bimj.201800124-blue)](http://doi.org/10.1002/bimj.201800124)

Supplementary material to **Assessment of local influence for the analysis of agreement** by Carla Leal, Manuel Galea and Felipe Osorio (Biometrical Journal 61 (4), 955-972. DOI: [10.1002/bimj.20180012](https://doi.org/10.1002/bimj.201800124))

Code written by: Carla Leal and Felipe Osorio

Correspondence author: Felipe Osorio, Email: felipe.osorios@usm.cl

Code tested on:
- R under development (2018-02-21 r74285), running Linux Mint 18.3 (64 bits)
- R version 3.3.0, running OS X 10.13.4 (64 bits)
- R version 3.4.3, running Windows 10 (64 bits)

Attached base packages: stats, graphics, grDevices, utils, datasets, methods, base

Other attached packages: heavy_0.38.19

CONTENTS:
- case_study/case_study.R: R commands for the analysis of transient sleep disorder
  dataset (described/analyzed at Sections 2.1 and 5 from manuscript).
- code/ccc.R: R functions for estimation and confidence intervals for the CCC.
- code/ccc.influence.R: R functions to compute influence measures for the CCC.
- code/poa.R: R functions for estimation and confidence intervals for the probability of agreement.
- code/poa.influence.R: R functions to compute influence measures for the probability of agreement.
- data/PSG.rda: clinical trial on transient sleep disorder dataset.
- simulation/ccc.simul.R: R function to compute the outlier detection percentage using
  different influence measures: CCC as objective function.
- simulation/poa.simul.R: R function to compute the outlier detection percentage using
  different influence measures: probability of agreement as objective function.
- simulation/simulation.R: R commands to perform the simulation study described at
  Section 4 from manuscript.
- supplement/ccc.2points.R: R functions to compute the outlier detection percentage
  using CCC as objective function (see Supplementary material).
- supplement/poa.2points.R: R functions to compute the outlier detection percentage
  using the probability of agreement as objective function (see Supplementary material).
- supplement/additional.simul.R: R commands to perform the additional simulation experiment
  described in the Supplementary material.
- ccc-manual.pdf: reference manual describing the functions included in files 'ccc.R',
  'poa.R', 'ccc.influence.R' and 'poa.influence.R'
- README.txt: this file.
