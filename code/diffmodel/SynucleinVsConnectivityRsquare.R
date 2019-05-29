#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/',sep='')
dir.create(savedir,recursive=T)
dir.create(paste(savedir,'roilevel',sep=''),recursive = T)

source('code/fitfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
Synuclein <- as.matrix(read.csv('Data83018/Snca.csv'))

nROIs <- length(conn.names)*2
ROInames <- c(sapply(conn.names, function(x) paste('R',x)),
              sapply(conn.names, function(x) paste('L',x)))
ROInames <- factor(ROInames, ordered = TRUE, levels = ROInames)

# get mean of each month
tp <- c(1,3,6)
Mice <- lapply(tp, function(M) path.data[path.data$`Time post-injection (months)` == M & path.data$Condition == grp,-(1:2)])
Grp.mean <- lapply(Mice, function(x) colMeans(x,na.rm = T))

##################################################
### Train model to predict path from iCPu seed ### 
##################################################

load(paste(params$opdir,'processed/Lout.RData',sep=''))

# Fit time scaling parameter
c.rng <- seq(0.01,10,length.out = 100) # scaling parameter
log.path <- lapply(Grp.mean, function(x) log(x,base=10))
c.Grp <- c.fit(log.path,L.out,tp,'R CPu',c.rng,ROInames)

########################
### Test base model  ###
########################

Xo <- make.Xo('R CPu',ROInames) # seed pathology in R CPu (actual injection site)
Xt.Grp <- lapply(as.list(tp), function(t) predict.Lout(L.out,Xo,c.Grp,t))

# make data frame of log path, log path prediction, and synuclein
df.by.month <- lapply(1:length(tp), function(M) data.frame(path = log(Grp.mean[[M]],base=10), Xt = log(Xt.Grp[[M]],base = 10),Syn=as.numeric(Synuclein)))
# remove regions with 0 pathology
mask.by.month <- lapply(df.by.month, function(X) X$path != -Inf & X$Xt != -Inf & !is.na(X$Xt))
# apply mask
df.by.month.mask <- lapply(1:length(tp), function(M) df.by.month[[M]][mask.by.month[[M]],])
# variance explained in pathology by connectivity
path.r.conn <- lapply(df.by.month.mask, function(X) cor(X$path,X$Xt)^2)
lapply(1:length(tp), function(M) print(paste('Month ',tp[M],', variance explained in pathology by connectivity: R^2 = ',path.r.conn[[M]],sep='')))

# variance explained in pathology by Synuclein
path.r.syn <- lapply(df.by.month.mask, function(X) cor(X$path,X$Syn)^2)
lapply(1:length(tp), function(M) print(paste('Month ',tp[M],', variance explained in pathology by Synuclein: R^2 = ',path.r.syn[[M]],sep='')))
# variance explained in pathology by connectivity and Synuclein
path.r.conn.syn <- lapply(df.by.month.mask, function(X) summary(lm.beta(lm(path~Xt+Syn,data=X))))
lapply(1:length(tp), function(M) print(paste('Month ',tp[M],', variance explained in pathology by Connectivity + Synuclein: R^2 = ',path.r.conn.syn[[M]]$r.squared,sep='')))

# variance explained in pathology by connectivity*Synuclein interaction
path.r.connbysyn <- lapply(df.by.month.mask, function(X) summary(lm.beta(lm(path~Xt*Syn,data=X))))
lapply(1:length(tp), function(M) print(paste('Month ',tp[M],', variance explained in pathology by Connectivity*Synuclein: R^2 = ',path.r.connbysyn[[M]]$r.squared,sep='')))

#############################
### Test synuclein model  ###
#############################

load(file = paste(params$opdir,'processed/Lout_syn.RData',sep=''))
c.rng <- seq(0.01,10,length.out = 100) # scaling parameter
log.path <- lapply(Grp.mean, function(x) log(x,base=10))
c.Grp <- c.fit(log.path,L.out,tp,'R CPu',c.rng)

Xo <- make.Xo('R CPu',ROInames) # seed pathology in R CPu (actual injection site)
Xt.Grp <- lapply(as.list(tp), function(t) predict.Lout(L.out,Xo,c.Grp,t))

# make data frame of log path, log path prediction, and synuclein
df.by.month <- lapply(1:length(tp), function(M) data.frame(path = log(Grp.mean[[M]],base=10), Xt = log(Xt.Grp[[M]],base = 10),Syn=as.numeric(Synuclein)))
# remove regions with 0 pathology
mask.by.month <- lapply(df.by.month, function(X) X$path != -Inf & X$Xt != -Inf & !is.na(X$Xt))
# apply mask
df.by.month.mask <- lapply(1:length(tp), function(M) df.by.month[[M]][mask.by.month[[M]],])
# variance explained in pathology by connectivity
path.r.conn <- lapply(df.by.month.mask, function(X) cor(X$path,X$Xt)^2)
lapply(1:length(tp), function(M) print(paste('Month ',tp[M],', variance explained in pathology by syn*connectivity: R^2 = ',path.r.conn[[M]],sep='')))

# variance explained in pathology by Synuclein
path.r.syn <- lapply(df.by.month.mask, function(X) cor(X$path,X$Syn)^2)
lapply(1:length(tp), function(M) print(paste('Month ',tp[M],', variance explained in pathology by Synuclein: R^2 = ',path.r.syn[[M]],sep='')))
# variance explained in pathology by connectivity and Synuclein
path.r.conn.syn <- lapply(df.by.month.mask, function(X) summary(lm.beta(lm(path~Xt+Syn,data=X))))
lapply(1:length(tp), function(M) print(paste('Month ',tp[M],', variance explained in pathology by Syn*Connectivity + Synuclein: R^2 = ',path.r.conn.syn[[M]]$r.squared,sep='')))
