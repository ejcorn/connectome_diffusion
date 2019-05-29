#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'nullmodels/anterograde',sep='')
dir.create(savedir,recursive=T)
dir.create(paste(savedir,'roilevel',sep=''),recursive = T)

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

##################################################
### Train model to predict path from iCPu seed ### 
##################################################

# stepwise prediction, i.e. from one month to the next
# Where i is row element and j is column element
# Wij is a connection from region i to region j
# convention is the opposite, so when I transpose W
# I am capturing "anterograde" connectivity

W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
n.regions <- nrow(W)
W <- W * !diag(n.regions) # get rid of diagonal
W <- W / (max(Re(eigen(W)$values))) # scale to max eigenvalue

in.deg <- colSums(W)
out.deg <- rowSums(W)
L.out <- diag(x = out.deg) - t(W) # outdegree laplacian, transpose
save(L.out,file = paste(savedir,'Loutanterograde.RData',sep=''))

# Fit time scaling parameter
c.rng <- seq(0.01,10,length.out = 100) # scaling parameter
log.path <- lapply(Grp.mean, function(x) log(x,base=10))
c.Grp <- c.fit(log.path,L.out,tp,ROI,c.rng,ROInames)

###############################################
### Test model at observed points for group ###
###############################################

Xo <- make.Xo(ROI,ROInames)
Xt.Grp <- lapply(as.list(tp), function(t) predict.Lout(L.out,Xo,c.Grp,t))
p.SC <- p.vuln <- list()
r.SC <- matrix(nrow=length(tp))
vulnerability <- mask <- list()
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
  vulnerability[[M]] <- matrix(0,nrow= length(mask[[M]]))
  vulnerability[[M]][mask[[M]]] <- residuals(lm(path~Xt,data=df))
  df2 <- data.frame(x=ROInames,y=vulnerability[[M]])
  p.vuln[[M]] <- ggplot(df2,aes(x= x, y= y)) + geom_col(fill='#5F4B8B',alpha = 0.8) + theme_classic() +
    ggtitle(paste('Month ',tp[M])) + xlab('') + ylab('Vulnerability') + 
    theme(text=element_text(size=8),axis.text.x = element_text(angle=90,size=4,vjust=0.5,hjust = 1),plot.title = element_text(size=8,hjust = 0.5))
}
save(r.SC,file = paste(savedir,grp,'cor.RData',sep=''))
plts.SC <- plot_grid(plotlist = p.SC,nrow=1)
ggsave(plts.SC,filename = paste(savedir,grp,'predictedpathcontinuous.pdf',sep=''),units = 'in',height = 1.5,width = 4.5)
plts.vuln <- plot_grid(plotlist = p.vuln,nrow=3)
ggsave(plts.vuln,filename = paste(savedir,grp,'vulnerabilitycontinuous.pdf',sep=''),units = 'in',height = 4.5,width = 5)

# compare predicted ts to observed for lowest vulnerability region
# get regions with path data for all 3 time points
mask <- do.call(what = 'cbind',args = mask)
mask <- rowSums(mask) == 3
df <- as.data.frame(cbind(do.call(what = 'cbind',args = Grp.mean),
                          do.call(what = 'cbind',args = Xt.Grp)))

df <- log(df,base = 10)
df <- df[mask,]
df[,4:6] <- (df[,4:6] - (df[,4] - df[,1]))
df.v <- do.call(what = 'cbind',args= vulnerability)
df.v <- df.v[mask,]
# find regions with low average vulnerability at all time points
low.vuln <- which(rowMeans(abs(df.v)) < 0.2)
for(ROI in low.vuln){
  # shift to observed
  p <- ggplot() + geom_point(aes(x=tp,y=unlist(df[ROI,1:3])),color = 'blue') +
    geom_line(aes(x=tp,y=unlist(df[ROI,4:6])),color='red') + 
    scale_x_continuous(breaks = tp,labels = tp) +
    scale_y_continuous(limits = c(min(df),max(df))) +
    xlab('Month') + ylab('log(Path)') + ggtitle(paste(grp,ROInames[ROI])) +
    theme(plot.title = element_text(size=8,hjust=0.5))
  p
  ggsave(plot = p,filename = paste(savedir,'roilevel/',grp,'_',ROInames[ROI],'.pdf',sep=''),width = 3,height = 3,units = 'in')
}

write.csv(do.call('cbind',vulnerability),paste(savedir,grp,'vulnerabilityEli.csv',sep=''),row.names = F)
write.csv(do.call('cbind',Xt.Grp),paste(savedir,grp,'predictedpathEli.csv',sep=''),row.names = F)
v.names <- c(sapply(conn.names, function(x) paste('i',x,sep='')),
             sapply(conn.names, function(x) paste('c',x,sep='')))

vulnerability.df <- as.data.frame(vulnerability)[orig.order,]
colnames(vulnerability.df) <- sapply(tp, function(x) paste('Month',x))
rownames(vulnerability.df) <- v.names[orig.order]
write.csv(vulnerability.df,paste(savedir,grp,'vulnerability.csv',sep=''))
predicted.df <- as.data.frame(Xt.Grp)[orig.order,]
colnames(predicted.df) <- sapply(tp, function(x) paste('Month',x))
rownames(predicted.df) <- v.names[orig.order]
write.csv(predicted.df,paste(savedir,grp,'predictedpath.csv',sep=''))

### pathology by months ##
sapply(1:length(tp), function(M1) sapply(1:length(tp), function(M2) 
  cor(Grp.mean[[M1]],Grp.mean[[M2]],use='pairwise.complete.obs',method = 'spearman')))

sapply(1:length(tp), function(M1) sapply(1:length(tp), function(M2) 
  cor(log.path[[M1]][log.path[[M1]] != -Inf & log.path[[M2]] != -Inf],log.path[[M2]][log.path[[M1]] != -Inf & log.path[[M2]] != -Inf],use='pairwise.complete.obs')))
