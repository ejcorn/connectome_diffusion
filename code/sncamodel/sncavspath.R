#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'sncamodel/',sep='')
dir.create(savedir,recursive=T)

source('code/fitfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
Synuclein <- as.matrix(read.csv('Data83018/Snca.csv'))

nROIs <- length(conn.names)*2
ROInames <- c(sapply(conn.names, function(x) paste('R',x)),
              sapply(conn.names, function(x) paste('L',x)))
ROInames <- factor(ROInames, ordered = TRUE, levels = ROInames)

tp <- c(1,3,6)
Mice <- lapply(tp, function(M) path.data[path.data$`Time post-injection (months)` == M & path.data$Condition == grp,-(1:2)])
Grp.mean <- lapply(Mice, function(x) colMeans(x,na.rm = T))


p.SC <- mask <- list()
r.SC <- matrix(nrow=length(tp))
for(M in 1:length(tp)) {
  df <- data.frame(path = log(Grp.mean[[M]],base=10),Syn = as.numeric(Synuclein))
  # exclude regions with 0 pathology at each time point for purposes of computing fit
  mask[[M]] <- df$path != -Inf & df$Syn != -Inf & !is.na(df$Syn)
  print(paste(sum(mask[[M]]),'regions')) # number of regions left after exclusion
  df <- df[mask[[M]],] 
  print(paste(grp, 'Month',tp[M],'p =',cor.test(df$path,df$Syn)$p.value))
  r.SC[M] <- cor(df$path,df$Syn)
  p.SC[[M]] <- ggplot(df,aes(x=Syn,y=path)) + geom_smooth(color = '#007257',method ='lm',size=1) + geom_point(color = '#007257',size = 1,alpha=0.6,stroke=0) +
    annotate(geom='text',x=max(df$Syn) - 0.2*diff(range(df$Syn)),y=min(df$path) + 0.1,label = paste('r =',signif(r.SC[M],2)),size=2.5) +
    theme_classic() + xlab('Snca expression') + ylab('log(Path.)') + ggtitle(paste('Month',tp[M])) +
    theme(text = element_text(size=8),plot.title = element_text(hjust=0.5,size=8))
}

plts.SC <- plot_grid(plotlist = p.SC,nrow=1)
ggsave(plts.SC,filename = paste(savedir,grp,'SynvsPath.pdf',sep=''),units = 'in',height = 1.5,width = 4.5)
