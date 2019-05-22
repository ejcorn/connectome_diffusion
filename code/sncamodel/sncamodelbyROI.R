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
savedir <- paste(params$opdir,'G20vsNTG/',sep='')
dir.create(savedir,recursive=T)

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
### Generate dist. of out-of-sample prediction ###
##################################################

W <- readMat('diffmodel/W.mat')$W
n.regions <- nrow(W)
W <- W * !diag(n.regions) # get rid of diagonal
W <- diag(as.numeric(Synuclein)) %*% W
W <- W / (max(Re(eigen(W)$values))) # scale to max eigenvalue

L.out <- get.L.out(W)

# fits from iCPU
c.rng <- seq(0.01,2,length.out = 50) # scaling parameter range
log.path <- lapply(Grp.mean, function(x) log(x,base=10))
c.Grp <- c.fit(log.path,L.out,tp,'R CPu',c.rng)

##################################
### Generate model time series ###
##################################
p.ts <- function(df.plt,tp,ttl){
  df.plt$ROI <- as.character(df.plt$ROI)
  p <- ggplot(df.plt,(aes(y=Xt,x=t,color=ROI))) + geom_line(size=0.25,alpha=0.5) +
    ggtitle(ttl) + ylab('Predicted Path.') +
    theme_classic() + theme(legend.position = 'none') + geom_vline(data=NULL,xintercept=tp,linetype=2) +
    theme(plot.title = element_text(size=8,hjust=0.5))
  return(p)
}

Xo <- matrix(0,nrow=n.regions)
Xo[which(ROInames == 'R CPu')] <- 1
n.t <- 100 # number of time points to plot
t.rng <- seq(0.1,36,length.out=n.t)
L.out[49,49] <- 0.8
Xt.Grp <- lapply(t.rng, function(t) expm(-L.out*t*c.Grp) %*% Xo)
Xt.Grp <- do.call('cbind',Xt.Grp)
t.plt <- t(do.call('cbind',lapply(1:n.regions, function(x) t.rng)))
r.plt <- do.call('cbind',lapply(1:n.t, function(x) 1:n.regions))
df.plt <- data.frame(Xt=as.vector(Xt.Grp),t=as.vector(t.plt),ROI=as.vector(r.plt))
df.plt <- df.plt[df.plt$ROI != which(ROInames == 'R CPu'),] # remove injected region
df.plt <- df.plt[df.plt$ROI <= n.regions/2,] # remove contralateral side
p <- p.ts(df.plt,tp,grp)
p
ggsave(p,filename = paste(savedir,grp,'predictedROItimeseries.pdf',sep=''),units = 'in',height = 2,width = 3)

### Understanding model dynamics ###

Xo <- matrix(0,nrow=n.regions)
Xo[which(ROInames == 'R CPu')] <- 1
n.t <- 100 # number of time points to plot
t.rng <- seq(0.1,500,length.out=n.t)
Xt.Grp <- lapply(t.rng, function(t) expm(-L.out*t*c.Grp) %*% Xo)
Xt.Grp <- do.call('cbind',Xt.Grp)
t.plt <- t(do.call('cbind',lapply(1:n.regions, function(x) t.rng)))
r.plt <- do.call('cbind',lapply(1:n.t, function(x) 1:n.regions))
df.plt <- data.frame(Xt=as.vector(Xt.Grp),t=as.vector(t.plt),ROI=as.vector(r.plt))

# do regions that peak earlier also have higher estimated pathology?
roi.t.to.peak <- sapply(1:n.regions,function(i) t.rng[which.max(df.plt$Xt[df.plt$ROI==i])])[-which(ROInames == 'R CPu')]
roi.peak <- sapply(1:n.regions,function(i) max(df.plt$Xt[df.plt$ROI==i]))[-which(ROInames == 'R CPu')]
hemi <- (1:n.regions)[-which(ROInames == 'R CPu')] >= n.regions/2
hemi <- ifelse(hemi,yes='Contra',no='Ipsi')
df.plt <- data.frame(x=roi.t.to.peak,y=roi.peak,r=hemi)
p <- ggplot(df.plt,aes(x=x,y=y,color=r)) + geom_point(size=0.5,alpha=0.5) + ggtitle('Regional Dynamics') + ylab('Peak Pathology') +
  xlab('Months to Peak') + scale_color_manual(values=c('red','blue'),name='') +
  theme_classic() + theme(legend.position = c(0.8,0.9)) + theme(plot.title = element_text(size=8,hjust=0.5))
p
ggsave(p,filename = paste(savedir,grp,'PeakPathvsTimeToPeak.pdf',sep=''),units = 'in',height = 2,width = 3)

