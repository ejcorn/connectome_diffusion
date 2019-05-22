rm(list=ls())
basedir <- '~/Dropbox/Neurodegeneration/MouseDiffusion/'
setwd(basedir)
opdir <- 'asyndiffusion1/'
dir.create(opdir,recursive = T)
params <- list(basedir=basedir,
               opdir=opdir,
               matlab.path='/Applications/MATLAB/R2018b/bin/matlab',
               grps = c('NTG','G20'))

####################
### Process data ###
####################

source('code/process/process_pathconn.R')
source('code/process/getLout.R')

################################
### Run base diffusion model ###
################################

for(grp in params$grps){source('code/diffmodel/analyzespread.R')}

grp <- 'NTG' # run these analyses only with NTG
  source('code/diffmodel/analyzespreadtraintest.R')
  source('code/diffmodel/traintestplot.R')
  source('code/diffmodel/seedspecificity.R')
  source('code/diffmodel/examineseedspecificity.R')
  source('code/diffmodel/forwardpredict.R')

##########################
### Analyze cell death ###
##########################

  source('celldeath/celldeath.R')

###########################
### Run synuclein model ###
###########################

for(grp in params$grps){
  source('code/sncamodel/sncamodel.R')
  #source('code/sncamodel/sncamodelbyROI.R')
}

###################################
### Run NTG vs. G20 comparisons ###
###################################

for(grp in params$grp){source('code/G20vsNTG/G20vsNTGsyntimeconst.R')}
source('code/G20vsNTG/plotG20vsNTGsyntc.R')

#######################
### Run null models ###
#######################

# run matlab scripts to generate distance matrix and rewired connectivity matrix
mat.savedir <- paste(basedir,params$opdir,'processed',sep='')
mat.cmd <- paste(params$matlab.path,' -nojvm -r \'cd(\'',basedir,'code/nullmodels/\'); homedir = \'',basedir,params$opdir,'\'; ',
                 'savedir = \'', mat.savedir,'\'; run(\'makenullconn.m\'); exit;\'',sep='')
system(mat.cmd)
grp <- 'NTG'
  source('code/nullmodels/analyzespread_anterograde.R')
  source('nullmodels/nullmodels.R')
