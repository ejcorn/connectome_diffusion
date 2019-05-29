#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'nullmodels/',sep='')
dir.create(savedir,recursive=T)

source('code/fitfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

nROIs <- length(conn.names)*2
ROInames <- c(sapply(conn.names, function(x) paste('R',x)),
              sapply(conn.names, function(x) paste('L',x)))
ROInames <- factor(ROInames, ordered = TRUE, levels = ROInames)
ROI = 'R CPu'

# get mean of each month
tp <- c(1,3,6)
Mice <- lapply(tp, function(M) path.data[path.data$`Time post-injection (months)` == M & path.data$Condition == grp,-(1:2)])
Grp.mean <- lapply(Mice, function(x) colMeans(x,na.rm = T))

#####################
### Rewiring Null ### 
#####################

D <- readMat(paste(params$opdir,'processed/Wrewire.mat',sep=''))$Wrw
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
L.out.d <- get.L.out(D)
# Fit time scaling parameter
c.rng <- seq(0.01,5,length.out = 100) # scaling parameter
log.path <- lapply(Grp.mean, function(x) log(x,base=10))
c.Grp <- c.fit(log.path,L.out.d,tp,ROI,c.rng,ROInames) # fit time scale
mask <- lapply(log.path, function(x) x != -Inf)
Xo <- make.Xo(ROI,ROInames)

Xt.Grp <- lapply(tp, function(t) predict.Lout(L.out.d,Xo,c.Grp,t))
x <- do.call('cbind',Xt.Grp)
p.SC <- p.vuln <- list()
r.SC <- matrix(nrow=length(tp))
for(M in 1:length(tp)) {
  df <- data.frame(path = log(Grp.mean[[M]],base=10), Xt = log(Xt.Grp[[M]],base = 10))
  mask[[M]] <- df$path != -Inf & df$Xt != -Inf & !is.na(df$Xt)
  df <- df[mask[[M]],]
  print(paste(sum(mask[[M]]),'regions'))
  print(paste(grp,'Month',tp[M],'p =',cor.test(df$path,df$Xt)$p.value))
  r.SC[M] <- cor(df$path,df$Xt)
  p.SC[[M]] <- ggplot(df,aes(x=Xt,y=path)) + geom_smooth(color = '#007257',method ='lm',size=1) + geom_point(color = '#007257',size = 1,alpha=0.6,stroke=0) +
    annotate(geom='text',x=max(df$Xt) - 0.2*diff(range(df$Xt)),y=min(df$path) + 0.1,label = paste('r =',signif(r.SC[M],2)),size=2.5) +
    theme_classic() + xlab('log(Predicted)') + ylab('log(Path.)') + ggtitle(paste('Month',tp[M])) +
    theme(text = element_text(size=8),plot.title = element_text(hjust=0.5,size=8))
}

plts.SC <- plot_grid(plotlist = p.SC,nrow=1)
ggsave(plts.SC,filename = paste(savedir,grp,'DegSeqPreservepredictedpathcontinuous.pdf',sep=''),units = 'in',height = 1.5,width = 4.5)

#######################
### Distance Matrix ###
#######################

D <- readMat(paste(params$opdir,'processed/D.mat',sep=''))$D
L.out.d <- get.L.out(D)
# Fit time scaling parameter
c.rng <- seq(0.01,3,length.out = 100) # scaling parameter
log.path <- lapply(Grp.mean, function(x) log(x,base=10))
c.Grp <- c.fit(log.path,L.out.d,tp,'R CPu',c.rng,ROInames) # fit time scale
mask <- lapply(log.path, function(x) x != -Inf)
n.regions <- length(ROInames)
Xo <- make.Xo('R CPu',ROInames) # seed ROI with path
Xt.Grp <- lapply(tp, function(t) predict.Lout(L.out.d,Xo,c.Grp,t))
x <- do.call('cbind',Xt.Grp)
p.SC <- p.vuln <- list()
r.SC <- matrix(nrow=length(tp))
for(M in 1:length(tp)) {
  df <- data.frame(path = log(Grp.mean[[M]],base=10), Xt = log(Xt.Grp[[M]],base = 10))
  mask[[M]] <- df$path != -Inf & df$Xt != -Inf & !is.na(df$Xt)
  df <- df[mask[[M]],]
  df <- df[-which(rownames(df) == 'iCPu'),] # get rid of outlier, which is the injected region
  print(paste(sum(mask[[M]]),'regions'))
  print(paste(grp,'Month',tp[M],'p =',cor.test(df$path,df$Xt)$p.value))
  r.SC[M] <- cor(df$path,df$Xt)
  p.SC[[M]] <- ggplot(df,aes(x=Xt,y=path)) + geom_smooth(color = '#007257',method ='lm',size=1) + geom_point(color = '#007257',size = 1,alpha=0.6,stroke=0) +
    annotate(geom='text',x=max(df$Xt) - 0.2*diff(range(df$Xt)),y=min(df$path) + 0.1,label = paste('r =',signif(r.SC[M],2)),size=2.5) +
    theme_classic() + xlab('log(Predicted)') + ylab('log(Path.)') + ggtitle(paste('Month',tp[M])) +
    scale_x_continuous(labels = function (x) round(x,3)) +
    theme(text = element_text(size=8),plot.title = element_text(hjust=0.5,size=8))
}

plts.SC <- plot_grid(plotlist = p.SC,nrow=1)
plts.SC

ggsave(plts.SC,filename = paste(savedir,grp,'DistanceMatpredictedpathcontinuous.pdf',sep=''),units = 'in',height = 1.5,width = 4.5)

############################
### Symmetric Connectome ###
############################

W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
W <- 0.5*(W + t(W))
#Synuclein <- as.matrix(read.csv('Data83018/Snca.csv'))
#W <- diag(as.numeric(Synuclein)) %*% W

L.out.d <- get.L.out(W)

# Fit time scaling parameter
c.rng <- seq(0.01,10,length.out = 100) # scaling parameter
log.path <- lapply(Grp.mean, function(x) log(x,base=10))
c.Grp <- c.fit(log.path,L.out.d,tp,'R CPu',c.rng,ROInames) # fit time scale
mask <- lapply(log.path, function(x) x != -Inf)
Xo <- make.Xo(ROI,ROInames)
Xt.Grp <- lapply(tp, function(t) predict.Lout(L.out.d,Xo,c.Grp,t))
x <- do.call('cbind',Xt.Grp)
p.SC <- p.vuln <- list()
r.SC <- matrix(nrow=length(tp))
for(M in 1:length(tp)) {
  df <- data.frame(path = log(Grp.mean[[M]],base=10), Xt = log(Xt.Grp[[M]],base = 10))
  mask[[M]] <- df$path != -Inf & df$Xt != -Inf & !is.na(df$Xt)
  df <- df[mask[[M]],]
  
  r.SC[M] <- cor(df$path,df$Xt)
  print(paste(sum(mask[[M]]),'regions'))
  print(paste(grp,'Month',tp[M],'p =',cor.test(df$path,df$Xt)$p.value))
  p.SC[[M]] <- ggplot(df,aes(x=Xt,y=path)) + geom_smooth(color = '#007257',method ='lm',size=1) + geom_point(color = '#007257',size = 1,alpha=0.6,stroke=0) +
    annotate(geom='text',x=max(df$Xt) - 0.2*diff(range(df$Xt)),y=min(df$path) + 0.1,label = paste('r =',signif(r.SC[M],2)),size=2.5) +
    theme_classic() + xlab('log(Predicted)') + ylab('log(Path.)') + ggtitle(paste('Month',tp[M])) +
    theme(text = element_text(size=8),plot.title = element_text(hjust=0.5,size=8))
}

plts.SC <- plot_grid(plotlist = p.SC,nrow=1)
ggsave(plts.SC,filename = paste(savedir,grp,'SymmetricWpredictedpathcontinuous.pdf',sep=''),units = 'in',height = 1.5,width = 4.5)
