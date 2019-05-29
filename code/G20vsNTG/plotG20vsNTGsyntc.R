rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'G20vsNTG/',sep='')

tp <- c(1,3,6)
tf = 0.5

load(paste(savedir,'NTGSyntimeconstantsTF',tf,'.RData',sep=''))
c.train.NTG <- c.train.Grp
r.NTG <- r.Grp
load(paste(savedir,'G20SyntimeconstantsTF',tf,'.RData',sep=''))
c.train.G20 <- c.train.Grp
r.G20 <- r.Grp
n.reps <- length(c.train.Grp)

df <- data.frame(y=c(unlist(c.train.NTG),unlist(c.train.G20)),
                 x=c(rep('NTG',n.reps),rep('G20',n.reps)))

w.t <- wilcox.test(x=unlist(c.train.NTG),y=unlist(c.train.G20),conf.int=T)
print(w.t)
p <- ggplot(df,aes(y=y,x=x)) + geom_jitter(alpha=0.5,size=1,stroke=0) + ylab('Time Constant') + theme_classic()+
  theme(text=element_text(size=8)) + xlab('') + geom_text(aes(x='G20',y= 0.35,label = paste('p =',signif(w.t$p.value,2))),size=2)
p
ggsave(plot = p,filename = paste(savedir,'SynTCjitterTF',tf,'.pdf',sep=''), height=2,width=2,units='in')
