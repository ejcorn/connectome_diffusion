rm(list= setdiff(ls(),c('params','grp')))

#####################
### Test accuracy ###
#####################

tp <- c(1,3,6)
tf = 0.5
load(file = paste(params$opdir,'diffmodel/traintest/',grp,'testfits',tf,'.RData',sep=''))
n.reps <- length(r.SC)
test.accuracy <- t(do.call(what = 'cbind',args = r.SC))
load(file = paste(params$opdir,'diffmodel/',grp,'cor.RData',sep=''))

# test null hypothesis that fits on all mice are not greater than expected under distribution of out of sample fits
# using one-sided non-parametric test, ask whether full sample fit is above 5th percentile
p.vals <- mapply(function(X,Y) {mean(Y<X)},X=as.data.frame(test.accuracy),Y=r.SC)
lb <- sapply(p.vals, function(x) paste('p =',signif(x,2)))
months <- as.vector(sapply(1:length(tp), function(M) rep(as.character(paste('Month',tp[M])),n.reps)))
p <- ggplot() + geom_violin(aes(x=months,y=as.vector(test.accuracy)),fill = '#5F4B8B',alpha = 0.5) +
  geom_point(aes(x=unique(months),y=r.SC),shape = 18,color = 'black',size=2) +
  annotate(geom='text',x=unique(months),y=0,label=lb,size=2)+
  scale_y_continuous(limits = c(0,1)) + ylab('Pearson\'s r') + xlab('') + ggtitle('Out of Sample Prediction') +
  theme_classic() + theme(text = element_text(size=8),plot.title = element_text(size=8,hjust=0.5)) 
p
ggsave(p,filename = paste(params$opdir,'diffmodel/traintest/',grp,'TestAccuracy',n.reps,'Reps',tf,'.pdf',sep=''),units = 'in',height = 2,width = 3)
