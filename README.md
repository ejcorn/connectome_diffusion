# connectome_diffusion

Code to reproduce all analysis in Henderson et al. 2019 ("α-Synuclein pathology spread through the brain connectome is modulated by selective vulnerability and genetic risk factors").

## Requirements:
  - MATLAB R2017a or later
  - R 3.3.3 or later. Requisite packages are listed in code/packages.R

## Directory structure

Master branch contains 2 major folders:
  - code/ contains folders, each of which contain scripts specific to certain analyses, i.e. ‘code/diffmodel’ contains code that uses linear diffusion models to predict spread of protein through structural connectome, 'code/G20vsNTG' deals with comparisons between G20 mice and NTG mice.
  - Data83018/ contains csv and xlsx files with 1) experimentally obtained pathology data and 2) parcellated Snca expression data and connectome data from Allen Brain Institute

## Input specification

The file ‘pipeline.R’ is located in the main directory. This file will coordinate the sequential execution of all scripts within the code/ folder, generating all the figures in the paper and more. Custom specification of the following inputs at the top of ‘pipeline.R’ is required:
  - basedir:  path to the main directory containing the 'code' and 'Data83018' folders 
  - matlab.path: path to MATLAB binary
  - opdir: name of output directory that contains all results, which will be housed in basedir
  - grps: character vector containing the name of groups in data file to test. For our data set, these were 'NTG' and 'G20'.

## Questions, suggestions, comments?

Please contact Eli Cornblath (Eli.Cornblath@pennmedicine.upenn.edu) with any questions regarding network analysis and code, and contact Mike Henderson (hendm@pennmedicine.upenn.edu) with any questions regarding experiments and data.
