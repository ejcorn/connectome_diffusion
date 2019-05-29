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

##############################################
### Compare to prediction from other seeds ###
##############################################

load(file = paste(params$opdir,'processed/Lout.RData',sep=''))

# fits from iCPU
c.rng <- seq(0.01,10,length.out = 100) # scaling parameter range
log.path <- lapply(Grp.mean, function(x) log(x,base=10))
c.Grp <- c.fit(log.path,L.out,tp,'R CPu',c.rng,ROInames) # fit time scale
mask <- lapply(log.path, function(x) x != -Inf) # find regions with 0 path at each time point to exclude
n.regions <- length(ROInames)
Xo <- make.Xo('R CPu',ROInames) # seed ROI with path
Xt <- lapply(1:length(tp), function(t) log(predict.Lout(L.out,Xo,c.Grp,tp[t]),base=10))
r.SC <- lapply(1:length(tp), function(t)
  cor(log.path[[t]][mask[[t]]],Xt[[t]][mask[[t]]]))

alt.seeds <- which(!(ROInames == 'R CPu'))
null.cors <- matrix(nrow=length(alt.seeds),ncol = length(tp))
c.rng <- seq(0.01,10,length.out = 40) # coarser range when fitting so many ROIs

for(S in 1:length(alt.seeds)){
  print(paste('Seed',S))
  seed <- alt.seeds[S]
  Xo.null <- matrix(0,nrow=n.regions)
  Xo.null[seed] <- 1
  c.null <- c.fit(log.path,L.out,tp,ROInames[seed],c.rng,ROInames)
  Xt.null <- lapply(c(1,3,6), function(t) predict.Lout(L.out,Xo.null,c.null,t))
  
  pred <- lapply(Xt.null, function(x) log(x,base=10))
  null.cors[S,] <- sapply(1:length(tp), function(M)
    cor(log.path[[M]][mask[[M]]],pred[[M]][mask[[M]]]))
}

# compute p-value as probability that other seeds predict observed path better than iCPu seed
all.fits <- rbind(null.cors,unlist(r.SC)) # join all fits to calculate the percentile of the fit
seed.region.pct <- unlist(sapply(1:length(tp), function(M) print(paste('Month',tp[M],'pctile =',100*calc.pctile(r.SC[[M]],all.fits[,M])))))
better.regions <- lapply(1:length(tp), function(M) which(r.SC[M] < null.cors[,M]))
save(better.regions,null.cors,file = paste(savedir,grp,'betterfitseeds.RData',sep=''))
better.regions <- unlist(better.regions)
as.character(ROInames[alt.seeds][better.regions])
sapply(1:length(tp), function(i) null.cors[better.regions[i],i])
months <- do.call('c',lapply(1:length(tp), function(M) rep(as.character(tp[M]),length(alt.seeds))))

p.null.seeds <- ggplot() + 
  geom_jitter(aes(x=months,y = as.vector(null.cors)),color ='#5F4B8B',alpha = 0.5,stroke = 0,size = 1, position = position_dodge(width=1)) +
  geom_point(aes(x=as.character(tp),y=as.numeric(r.SC)),shape = 18,color = 'black',size=2) + 
  #geom_text(aes(x=as.character(tp),y=0.8,label = seed.region.pvals),size=2.5) +
  xlab('Month') + ylab('Fit') + ggtitle('Actual vs. Random Seed') +
  theme(text = element_text(size=8),plot.title = element_text(hjust=0.5,size=8)) +
  theme(axis.text.x = element_text(size=8)) + theme(axis.text.y = element_text(size=8))
p.null.seeds
ggsave(p.null.seeds,filename = paste(savedir,grp,'actualvsrandomseed.pdf',sep=''),units = 'in',height = 2,width = 3)
