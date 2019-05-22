library(ggplot2)
library(stringr)
library(R.matlab)
library(cowplot)
library(expm)
library(viridis)

#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'celldeath/',sep='')
dir.create(savedir,recursive=T)
source('code/fitfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

nROIs <- length(conn.names)*2
ROInames <- c(sapply(conn.names, function(x) paste('i',x,sep='')),
              sapply(conn.names, function(x) paste('c',x,sep='')))
ROInames <- factor(ROInames, ordered = TRUE, levels = ROInames)

cell.death <- read.csv('Data83018/CellDeath.csv',header = T, check.names = F)

# remove time and condition
# compute cell death as 6M mean cell count - 3M mean cell count
# for NTG and G2019S

cell.death.Grp <- data.frame(M13=with(cell.death, 
  colMeans(cell.death[Condition == grp & Time == '1 MPI',-(1:2)]) -
  colMeans(cell.death[Condition == grp & Time == '3 MPI',-(1:2)])),
  M36=with(cell.death,
  colMeans(cell.death[Condition == grp & Time == '3 MPI',-(1:2)]) -
  colMeans(cell.death[Condition == grp & Time == '6 MPI',-(1:2)]))
)

# remove extra regions
names(cell.death.Grp[!rownames(cell.death.Grp) %in% ROInames])
names(ROInames) = NULL
ROInames[!ROInames %in% rownames(cell.death.Grp)]
# reorder names

cell.death.Grp <- cell.death.Grp[rownames(cell.death.Grp) %in% ROInames,]
idx <- order(match(rownames(cell.death.Grp),ROInames))
cell.death.Grp <- cell.death.Grp[idx,]
identical(rownames(cell.death.Grp),as.character(ROInames))
save(file = paste(savedir,'celldeath',grp,'.RData',sep=''),cell.death.Grp)

##############################################################
### assess covariance between cell death and path measures ###
##############################################################

# get mean of each month
tp <- c(1,3,6)
Mice <- lapply(tp, function(M) path.data[path.data$`Time post-injection (months)` == M & path.data$Condition == grp,-(1:2)])
Grp.mean <- lapply(Mice, function(x) colMeans(x,na.rm = T))
log.path <- as.data.frame(lapply(Grp.mean,function(x) log(x,base=10)))
mask <- lapply(log.path, function(x) x != -Inf)

p <- list()
lab <- c('M1-M3 Counts','M3-M6 Counts')
for(M in 1:2){
  df <- data.frame(x=log.path[,M],y=cell.death.Grp[,M])
  df <- df[mask[[M]],]
  df <- df[df$y > 0,]
  df <- df[-which(df$y==max(df$y)),] # remove outlier
  p[[M]] <- ggplot(df,aes(x=x,y=y)) + geom_smooth(color = '#007257',method ='lm',size=1) + 
    geom_point(color = '#007257',size = 1,alpha=0.6,stroke=0) + theme_classic() +
    annotate(geom='text',x=max(df$x) - 0.2*diff(range(df$x)),y=min(df$y) + 0.1,label = paste('r =',signif(cor(df$x,df$y),2)),size=2.5) +
    theme_classic() + xlab('log(Path.)') + ylab(lab[M]) + ggtitle(paste('Month',tp[M])) +
    theme(text = element_text(size=8),plot.title = element_text(hjust=0.5,size=8))
}

p.P.CD <- plot_grid(plotlist = p,nrow=1)
ggsave(p.P.CD,filename = paste(savedir,grp,'CellDeathDeltaVsPathM1M3DeathRegionsOnly.pdf',sep=''),units = 'in',height = 1.5,width = 3)

vuln <- read.csv(paste(params$opdir,'diffmodel/',grp,'vulnerabilityEli.csv',sep=''))
ROInames.RL <- c(sapply(conn.names, function(x) paste('R ',x,sep='')),
              sapply(conn.names, function(x) paste('L ',x,sep='')))
vuln <- calc.vuln(vulnerability = vuln,ROInames = ROInames.RL)
p <- list()
lab <- c('Month 3 Death','Month 6 Death')
for(M in 1:2){
  df <- data.frame(x=vuln,y=cell.death.Grp[,M])
  df <- df[df$y > 0,]
  df <- df[-which(df$y==max(df$y)),] # remove outlier
  p[[M]] <- ggplot(df,aes(x=x,y=y)) + geom_smooth(color = '#007257',method ='lm',size=1) + 
    geom_point(color = '#007257',size = 1,alpha=0.6,stroke=0) + theme_classic() +
    annotate(geom='text',x=max(df$x) - 0.2*diff(range(df$x)),y=min(df$y) + 0.1,
             label = paste('r =',signif(cor(df$x,df$y),2)),size=2.5,vjust=2) +
    theme_classic() + xlab('Vulnerability') + ylab(lab[M]) + ggtitle('') +
    theme(text = element_text(size=8),plot.title = element_text(hjust=0.5,size=8))
  print(cor.test(df$x,df$y))
}

p.V.CD <- plot_grid(plotlist = p,nrow=1)
ggsave(p.V.CD,filename = paste(savedir,grp,'CellDeathDeltaVsVulnerabilityM1M3DeathRegionsOnly.pdf',sep=''),units = 'in',height = 1.5,width = 3)