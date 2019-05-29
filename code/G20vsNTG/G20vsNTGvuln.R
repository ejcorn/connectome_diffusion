#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'G20vsNTG/',sep='')
dir.create(savedir,recursive=T)
source('code/fitfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
Synuclein <- as.matrix(read.csv('Data83018/Snca.csv'))


nROIs <- length(conn.names)*2
ROInames <- c(sapply(conn.names, function(x) paste('R ',x,sep='')),
              sapply(conn.names, function(x) paste('L ',x,sep='')))
ROInames <- factor(ROInames, ordered = TRUE, levels = ROInames)

# get mean of each month
tp <- c(1,3,6)
NTG <- lapply(tp, function(M) path.data[path.data$`Time post-injection (months)` == M & path.data$Condition == 'NTG',-(1:2)])
NTG.mean <- lapply(NTG, function(x) colMeans(x,na.rm = T))
G20 <- lapply(tp, function(M) path.data[path.data$`Time post-injection (months)` == M & path.data$Condition == 'G20',-(1:2)])
G20.mean <- lapply(G20, function(x) colMeans(x,na.rm = T))

##############################################################
### Understand how regional vulnerability changes with G20 ###
##############################################################

NTG.vuln <- read.csv(paste(params$opdir,'diffmodel/NTGvulnerabilityEli.csv',sep=''))

# how does vulnerability predict changes observed with G20?

G20.path <- lapply(G20.mean,function(X) log(X,base=10))
NTG.path <- lapply(NTG.mean,function(X) log(X,base=10))
mask <- mapply(function(X,Y) list(X != -Inf & Y != -Inf), X=NTG.path,Y=G20.path)

# correlate vulnerability from base model 
c.test <-mapply(function(G20,NTG,vuln,mask) list(cor.test(G20[mask] - NTG[mask],vuln[mask])),
       G20=G20.path,NTG=NTG.path,vuln=as.list(NTG.vuln),mask=mask)

p.vuln <- list()
for(M in 1:length(tp)){
  df <- data.frame(dif= G20.path[[M]][mask[[M]]] - NTG.path[[M]][mask[[M]]],NTG = NTG.vuln[,M][mask[[M]]])
  print(paste('Month',tp[M]))
  print(c.test[[M]])
  lb <- paste('r ==',signif(c.test[[M]]$estimate,2))
  p.vuln[[M]] <- ggplot(data=df,aes(y=dif,x=NTG)) + geom_point(color ='#5F4B8B',alpha = 0.5,stroke = 0,size = 1) + 
    geom_smooth(method='lm',color ='#5F4B8B') +
    annotate(geom='text',label=lb,y=min(df$dif),x=min(df$NTG),hjust=0,vjust=0,parse=T,size=2) +
    xlab('NTG Vulnerability') + ylab('log(G20/NTG) Pathology') + ggtitle(paste('Month',tp[M])) +
    theme_classic() + theme(text=element_text(size=6),plot.title = element_text(size=8,hjust=0.5)) #+
}

# print corrected p-values
p.vals <- p.adjust(sapply(c.test, function(X) X$p.value),method='bonferroni')
sapply(1:length(tp), function(M) print(paste('Month',tp[M],'p_corr =',p.vals[M])))

plts.vuln<- plot_grid(plotlist = p.vuln,nrow=1)
ggsave(plts.vuln,filename = paste(savedir,'NTGVulnvsG20-NTGPath.pdf',sep=''),units = 'in',height = 1.5,width = 4.5)