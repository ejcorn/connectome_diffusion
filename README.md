# connectome_diffusion

Code to reproduce all analysis in Henderson et al. 2019 ("Quantitative α-synuclein pathology mapping and network analysis to model pathological protein spread and selective vulnerability").

## Requirements:
  - MATLAB R2017a or later
  - R 3.3.3 or later. Requisite packages are listed in code/misc/packages.R
  - Brain Connectivity Toolbox, MATLAB version: https://sites.google.com/site/bctnet/

## Directory structure

Master branch contains 4 major folders:
  - code/ contains folders, each of which contain scripts specific to certain analyses, i.e. ‘code/sncamodel’ contains code that ***
  - Data83018/ contains csv and xlsx files with ***

## Input specification

The file ‘pipeline.R’ is located in the main directory. This file will coordinate the sequential execution of all scripts within the code/ folder, generating all the figures in the paper and more. At the top of ‘pipeline.R,’ one must select their own home directory (path to the main directory containing the above 4 folders) and path to MATLAB binary (matlab.path), in addition to the following variables:
  - ***

## Questions, suggestions, comments?

Please contact Eli Cornblath (Eli.Cornblath@pennmedicine.upenn.edu) with any questions regarding network analysis and code, and contact Mike Henderson (hendm@pennmedicine.upenn.edu) with any questions regarding experiments and input data.
