# ABO-PREG

This repository contains code for reproducing the Mendelian randomization examined in ' The effect of levels of BGAT protein on pregnancy related outcomes：a Mendelian randomization study'. In this study, we estimated the effect of the protein level of Histo-Blood Group ABO System Transferase (BGAT) on ten pregnancy related outcomes using genetic instruments respectively from deCODE, INTERVAL and Fenland. Our study highlighted BGAT as a putative causal protein for venous complications and haemorrhoids in pregnancy, which requires further investigation to validate it as a potential drug target.


## To start using the code, you need to install TwoSampleMR package:

install.packages("remotes")

remotes::install_github("MRCIEU/TwoSampleMR")

devtools::install_github("mrcieu/ieugwasr")

## The following R scripts are included in this repository:

1.	MR analysis: MR analysis: estimate the effect of protein level of BGAT on pregnancy outcomes using genetic instruments from deCODE
2.	MR analysis: estimate the effect of protein level of BGAT on pregnancy outcomes using genetic instruments from INTERVAL and Fenland
3.	Colocalisation analysis

## The following data are included in this repository:

1.	Genetic instruments for BGAT from datasets deCODE, INTERVAL and Fenland
2.	Sample input for Colocalisation analysis



