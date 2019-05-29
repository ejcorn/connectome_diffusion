#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/forward/',sep='')
dir.create(savedir,recursive=T)
source('code/fitfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

nROIs <- length(conn.names)*2
ROInames <- c(sapply(conn.names, function(x) paste('i',x,sep='')),
              sapply(conn.names, function(x) paste('c',x,sep='')))
ROInames <- factor(ROInames, ordered = TRUE, levels = ROInames)

# get mean of each month
tp <- c(1,3,6)
Mice <- lapply(tp, function(M) path.data[path.data$`Time post-injection (months)` == M & path.data$Condition == grp,-(1:2)])
Grp.mean <- lapply(Mice, function(x) colMeans(x,na.rm = T))

##################
### Base Model ### 
##################

load(paste(params$opdir,'processed/Lout.RData',sep=''))

c.rng <- seq(0.01,10,length.out = 100) # scaling parameter
log.path <- lapply(Grp.mean, function(x) log(x,base=10))
c.Grp <- c.fit(log.path,L.out,tp,'iCPu',c.rng,ROInames)

ROIs <- c('iSN','iPir','iHipp','iM2','iCPu')

for(ROI in ROIs){
  Xo <- matrix(0,nrow=n.regions)
  Xo[which(ROInames == ROI)] <- 1 # seed ROI with path
  Xt.Grp <- do.call('cbind',lapply(tp, function(t) log(predict.Lout(L.out,Xo,c.Grp,t),base=10)))
  colnames(Xt.Grp) <- c('Month 1','Month 3','Month 6')
  rownames(Xt.Grp) <- ROInames
  Xt.Grp <- Xt.Grp[orig.order,]
  write.csv(x = Xt.Grp,file = paste(savedir,grp,'Base',ROI,'.csv',sep=''),row.names = T)
}


#######################
### Synuclein Model ### 
#######################

load(paste(params$opdir,'processed/Lout_syn.RData',sep=''))

c.rng <- seq(0.01,10,length.out = 100) # scaling parameter
log.path <- lapply(Grp.mean, function(x) log(x,base=10))
c.Grp <- c.fit(log.path,L.out,tp,'iCPu',c.rng,ROInames)

ROIs <- c('iSN','iPir','iHipp','iM2','iCPu')

for(ROI in ROIs){
  Xo <- matrix(0,nrow=n.regions)
  Xo[which(ROInames == ROI)] <- 1 # seed ROI with path
  Xt.Grp <- do.call('cbind',lapply(tp, function(t) log(predict.Lout(L.out,Xo,c.Grp,t),base=10)))
  colnames(Xt.Grp) <- c('Month 1','Month 3','Month 6')
  rownames(Xt.Grp) <- ROInames
  Xt.Grp <- Xt.Grp[orig.order,]
  write.csv(x = Xt.Grp,file = paste(savedir,grp,'Syn',ROI,'.csv',sep=''),row.names = T)
}