library(ggplot2)
library(stringr)
library(R.matlab)
library(cowplot)
library(expm)
library(viridis)
library(RColorBrewer)
library(lm.beta)
library(mgcv)

#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/seedspec/',sep='')
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

######################################################
### Examine seed regions with better fits than CPu ###
######################################################

load(file = paste(savedir,grp,'betterfitseeds.RData',sep=''))
alt.seeds <- which(!(ROInames == 'R CPu')) # identify all the regions which were not injected
W <- readMat('diffmodel/W.mat')$W
n.regions <- nrow(W)
W <- W * !diag(n.regions) # get rid of diagonal
in.deg <- colSums(W)
out.deg <- rowSums(W)

# test if fit specificity is explained by in degree, out degree
lapply(as.data.frame(null.cors), function(x) cor.test(out.deg[alt.seeds],x))
lapply(as.data.frame(null.cors), function(x) cor.test(in.deg[alt.seeds],x))
# non-parametric test for where better fit regions lie in terms of in deg, out deg, part. coef
lapply(unique(unlist(better.regions)), function(x) mean(in.deg[alt.seeds][x] > in.deg[alt.seeds][-x]))
lapply(unique(unlist(better.regions)), function(x) mean(out.deg[alt.seeds][x] > out.deg[alt.seeds][-x]))

# compute similarity of connectivity to iCPU
iCPu <- which(ROInames == 'R CPu')
conn.sim.out <- sapply(alt.seeds, function(R) cor(W[R,],W[iCPu,]))
conn.sim.in <- sapply(alt.seeds, function(R) cor(W[,R],W[,iCPu]))
conn.sim.out <- fisher.r.to.z(conn.sim.out)
conn.sim.in <- fisher.r.to.z(conn.sim.in)
plot(conn.sim.in,null.cors[,3])
cor(conn.sim.in,null.cors)

summary(lm.beta(lm(null.cors[,1] ~ conn.sim.in + conn.sim.out + out.deg[alt.seeds])))

# plot

p.spec <- list()
r.c <- list()
for(M in 1:length(tp)){
  df <- data.frame(conn = conn.sim.in,nullcor = null.cors[,M])
  m <- gam(nullcor ~ s(conn,k=3),data = df,method = 'REML')
  p.val <- signif(summary(m)$s.table[,'p-value'],2)
  print(p.val)
  r.sq <- signif(summary(m)$r.sq,2)
  lb <- paste('R^2 ==',r.sq)
  # "REML and ML are less prone to local minima and may be preferred" 
  # https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gam.selection.html
  p.spec[[M]] <- ggplot(data=df,aes(x=conn,y=nullcor)) + geom_point(color ='#5F4B8B',alpha = 0.5,stroke = 0,size = 1) + 
    geom_smooth(method='gam',formula=y~s(x,k=3),color ='#5F4B8B') +
    annotate(geom='text',label=lb,x=max(df$conn),y=min(df$nullcor),hjust=1,vjust=0,parse=T,size=2) +
    xlab('In-Projection Similarity \n to iCPu') + ylab('Fit to iCPu Injection') + ggtitle(paste('Month',tp[M])) +
    theme_classic() + theme(text=element_text(size=6),plot.title = element_text(size=8,hjust=0.5)) #+
    #theme(axis.title.x = element_text(size=8))
  
}

plts.spec <- plot_grid(plotlist = p.spec,nrow=1)
ggsave(plts.spec,filename = paste(savedir,grp,'AlternateSeedFitByInProjectionSimilarity.pdf',sep=''),units = 'in',height = 1.5,width = 4.5)

p.spec <- list()
r.c <- list()
for(M in 1:length(tp)){
  df <- data.frame(conn = conn.sim.out,nullcor = null.cors[,M])
  m <- gam(nullcor ~ s(conn,k=3),data = df,method = 'REML')
  p.val <- signif(summary(m)$s.table[,'p-value'],2)
  print(p.val)
  r.sq <- signif(summary(m)$r.sq,2)
  lb <- paste('R^2 ==',r.sq)
  p.spec[[M]] <- ggplot(data=df,aes(x=conn,y=nullcor)) + geom_point(color ='#5F4B8B',alpha = 0.5,stroke = 0,size = 1) + 
    geom_smooth(method='gam',formula=y~s(x,k=3),color ='#5F4B8B') +
    annotate(geom='text',label=lb,x=max(df$conn),y=min(df$nullcor),hjust=1,vjust=0,parse=T,size=2) +
    xlab('Out-Projection Similarity \n to iCPu') + ylab('Fit to iCPu Injection') + ggtitle(paste('Month',tp[M])) +
    theme_classic() + theme(text=element_text(size=6),plot.title = element_text(size=8,hjust=0.5)) #+
  #theme(axis.title.x = element_text(size=8))
  
}

plts.spec <- plot_grid(plotlist = p.spec,nrow=1)
ggsave(plts.spec,filename = paste(savedir,grp,'AlternateSeedFitByOutProjectionSimilarity.pdf',sep=''),units = 'in',height = 1.5,width = 4.5)
