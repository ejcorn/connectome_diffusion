#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'G20vsNTG/',sep='')
dir.create(savedir,recursive=T)

source('code/fitfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
Synuclein <- as.matrix(read.csv('Data83018/Snca.csv'))

nROIs <- length(conn.names)*2
ROInames <- c(sapply(conn.names, function(x) paste('R',x)),
              sapply(conn.names, function(x) paste('L',x)))
ROInames <- factor(ROInames, ordered = TRUE, levels = ROInames)

#####################
### Prepare model ###
#####################

W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
n.regions <- nrow(W)
W <- W * !diag(n.regions) # get rid of diagonal
W <- diag(as.numeric(Synuclein)) %*% W
W <- W / (max(Re(eigen(W)$values))) # scale to max eigenvalue

L.out <- get.L.out(W)

# fits from iCPU
c.rng <- seq(0.01,2,length.out = 50) # scaling parameter range

tp <- c(1,3,6)
path.data.Grp <- lapply(tp,function(M) path.data[path.data$`Time post-injection (months)` == M & path.data$Condition == grp,-(1:2)])
ROI <- 'R CPu'
Xo <- make.Xo(ROI,ROInames)

n.reps <- 100
tf <- 0.5 # training fraction
c.train.Grp <- r.Grp <- list()
for(REP in 1:n.reps){
  print(paste('Rep',REP))
  # mark training samples for each time point
  train.idx.Grp <- lapply(path.data.Grp, function(X) sample(1:nrow(X),size = tf*nrow(X)))
  # split into non-overlapping training and testing sets
  path.data.train.Grp <- mapply(function(X,t.idx) {list(as.matrix(X[t.idx,]))}, X=path.data.Grp,t.idx=train.idx.Grp)
  # compute train and test means
  Grp.mean.train <- lapply(path.data.train.Grp, function(x) log(colMeans(x,na.rm = T),base=10))
  # fit time scale parameter on training data
  Grp.fit <- c.fit.r(Grp.mean.train,L.out,tp,ROI,c.rng,ROInames)
  c.train.Grp[[REP]] <- Grp.fit$c
  r.Grp[[REP]] <- Grp.fit$r
  
}

save(c.train.Grp,r.Grp,file = paste(savedir,grp,'SyntimeconstantsTF',tf,'.RData',sep=''))
