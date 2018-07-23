Supplementary material to "Local influence for the concordance correlation coefficient" by Carla Leal, Manuel Galea and Felipe Osorio

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
- code/ccc.influence.R: R functions to compute influence measures for the CCC.
- code/ccc.R: R functions for estimation and confidence intervals for the CCC.
- data/PSG.rda: clinical trial on transient sleep disorder dataset.
- simulation/ccc.simul.R: R funtion to compute the outlier detection percentage using
  different influence measures (Section 4 from manuscript).
- simulation/simulation.R: R commands to perform the simulation study described at
  Section 4 from manuscript.
- ccc-manual.pdf: reference manual describing the functions included in files 'ccc.R'
  and 'ccc.influence.R'.
- README.md: this file.
