library(ggplot2)
library(stringr)
library(R.matlab)
library(cowplot)
library(expm)
library(viridis)
library(RColorBrewer)

#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/traintest/',sep='')
dir.create(savedir,recursive=T)
source('code/fitfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

nROIs <- length(conn.names)*2
ROInames <- c(sapply(conn.names, function(x) paste('R',x)),
              sapply(conn.names, function(x) paste('L',x)))
ROInames <- factor(ROInames, ordered = TRUE, levels = ROInames)

# get mean of each month
tp <- c(1,3,6)
Mice <- lapply(tp, function(M) path.data[path.data$`Time post-injection (months)` == M & path.data$Condition == grp,-(1:2)])
Grp.mean <- lapply(Mice, function(x) colMeans(x,na.rm = T))

##################################################
### Generate dist. of out-of-sample prediction ###
##################################################

load(file = paste(params$opdir,'processed/Lout.RData',sep=''))
# fits from iCPU
c.rng <- seq(0.01,10,length.out = 40) # scaling parameter range
log.path <- lapply(Grp.mean, function(x) log(x,base=10))

tp <- c(1,3,6)
path.data <- lapply(tp,function(M) path.data[path.data$`Time post-injection (months)` == M & path.data$Condition == grp,-(1:2)])
ROI <- 'R CPu'
n.regions <- length(ROInames)
r.SC <- list()
Xo <- matrix(0,nrow=n.regions)
Xo[which(ROInames == 'R CPu')] <- 1 # seed ROI with path

n.reps <- 100
tf <- 0.5
for(REP in 1:n.reps){
  print(paste('Rep',REP))
  # mark training samples for each time point
  train.idx <- lapply(path.data, function(X) sample(1:nrow(X),size = tf*nrow(X)))
  # split into non-overlapping training and testing sets
  path.data.train <- mapply(function(X,t.idx) {as.matrix(X[t.idx,])}, X=path.data,t.idx=train.idx)
  path.data.test <- mapply(function(X,t.idx) {as.matrix(X[-t.idx,])}, X=path.data,t.idx=train.idx)
  # compute train and test means
  Grp.mean.train <- lapply(path.data.train, function(x) log(colMeans(x,na.rm = T),base=10))
  Grp.mean.test <- lapply(path.data.test, function(x) log(colMeans(x,na.rm = T),base=10))
  
  # fit time scale parameter on training data
  c.train <- c.fit(Grp.mean.train,L.out,tp,ROI,c.rng)
  # predict path
  Xt <- lapply(1:length(tp), function(t) log(expm(-L.out*tp[t]*c.train) %*% Xo,base=10))
  # compute fit with test data
  mask <- lapply(Grp.mean.test, function(x) x != -Inf)
  r.SC[[REP]] <- sapply(1:length(tp), function(t)
    cor(Grp.mean.test[[t]][mask[[t]]],Xt[[t]][mask[[t]]]))
}

save(r.SC,file = paste(savedir,grp,'testfits',tf,'.RData',sep=''))
