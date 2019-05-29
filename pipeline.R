rm(list=ls())
basedir <- '~/Dropbox/Neurodegeneration/MouseDiffusion/connectome_diffusion/'
setwd(basedir)
opdir <- 'asyndiffusion3/'
dir.create(opdir,recursive = T)
params <- list(basedir=basedir,
               opdir=opdir,
               matlab.path='/Applications/MATLAB_R2018b.app/bin/matlab',
               grps = c('NTG','G20'))

source('code/packages.R')

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

###########################
### Run synuclein model ###
###########################

for(grp in params$grps){
  source('code/sncamodel/sncamodel.R')
  source('code/sncamodel/sncavspath.R')
  source('code/sncamodel/sncamodelbyROI.R')
}

###################################
### Run NTG vs. G20 comparisons ###
###################################

for(grp in params$grp){source('code/G20vsNTG/G20vsNTGsyntimeconst.R')}
source('code/G20vsNTG/plotG20vsNTGsyntc.R')
source('code/G20vsNTG/G20vsNTGvuln.R')

#######################
### Run null models ###
#######################

# run matlab scripts to generate distance matrix and rewired connectivity matrix
mat.savedir <- paste(params$basedir,params$opdir,'processed',sep='')
mat.cmd <- paste(params$matlab.path,' -nojvm -r \"cd(\'',params$basedir,'code/nullmodels/\'); homedir = \'',params$basedir,'\'; ',
                 'savedir = \'', mat.savedir,'\'; run(\'makenullconn.m\'); exit;\"',sep='')
system(mat.cmd)
grp <- 'NTG'
source('code/nullmodels/analyzespread_anterograde.R')
source('code/nullmodels/nullmodels.R')
