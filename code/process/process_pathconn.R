#################
### Load data ###
#################

rm(list=setdiff(ls(),'params'))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'processed/',sep='')
dir.create(savedir,recursive=T)

data <- read.csv('Data83018/data.csv',header = TRUE,check.names = F)
connectivity.ipsi <- read.csv('Data83018/connectivity_ipsi.csv',row.names = 1,header = T,check.names = F)
connectivity.contra <- read.csv('Data83018/connectivity_contra.csv',row.names = 1,header = T, check.names = F)

############################
### Process region names ###
############################

parInd <- sapply(colnames(connectivity.contra), function(X) str_locate_all(pattern ='\\(',X)[[1]][1] - 2)
conn.names <- sapply(1:length(parInd), function(i) substring(colnames(connectivity.contra)[i],1,parInd[i]))

path.names <- colnames(data)[-(1:2)]
cInd <- which(sapply(path.names, function(X) substr(X,1,1) == 'c'))
iInd <- which(sapply(path.names, function(X) substr(X,1,1) == 'i'))
path.names.contra <- sapply(cInd, function(X) substr(path.names[X],2,nchar(path.names[X])))
path.names.ipsi <- sapply(iInd, function(X) substr(path.names[X],2,nchar(path.names[X])))

# order path and connectivity by same names

path.data.ipsi <- data[,-(1:2)][,iInd]
path.data.contra <- data[,-(1:2)][,cInd]
path.data <- cbind(data[,(1:2)],path.data.ipsi[,order(match(path.names.ipsi,conn.names))],
                   path.data.contra[,order(match(path.names.contra,conn.names))])
colnames(path.data)[2] <- 'Condition'

# tile matrix such that sources are rows, columns are targets (see Oh et al. 2014 Fig 4)
W <- rbind(cbind(connectivity.ipsi,connectivity.contra),cbind(connectivity.contra,connectivity.ipsi))
n.regions <- nrow(W)

# retain indices to reorder like original data variable for plotting on mouse brains
ROInames <- c(sapply(conn.names, function(x) paste('i',x,sep='')),
              sapply(conn.names, function(x) paste('c',x,sep='')))
orig.order <- order(match(ROInames,names(data)[-(1:2)]))

# save data
writeMat(paste(savedir,'W.mat',sep=''),W=as.matrix(W))
save(path.data, conn.names, orig.order, n.regions, file = paste(savedir,'pathdata.RData',sep=''))

###############################
### Process ROI Coordinates ###
###############################

ROInames <- c(sapply(conn.names, function(x) paste('i',x,sep='')),
              sapply(conn.names, function(x) paste('c',x,sep='')))

coor <- read.csv('Data83018/ROIcoords.csv')
coor.names <- coor$X
idx <- order(match(coor.names,ROInames))
write.csv(coor[idx,],file = paste(savedir,'coor.csv',sep=''),row.names = F)

###############################
### Process Snca Expression ###
###############################

Synuclein <- read.csv('Data83018/SncaExpression.csv',header = F,stringsAsFactors = F, row.names = 1)
Synuclein <- Synuclein[order(match(rownames(Synuclein),ROInames)),]
write.csv(Synuclein,file = paste(savedir,'Snca.csv',sep=''),row.names = F)

