c.fit.old <- function(log.path,L.out,tp,ROI,c.rng){
  n.regions <- length(log.path[[1]])
  # initialize Xo vector, path at t=0
  Xo <- matrix(0,nrow=n.regions)
  # seed ROI with max path at month 1
  Xo[which(ROInames == ROI)] <- 1 #max(exp(log.path[[1]]))
  # exclusion mask... don't count regions with 0 path
  mask <- lapply(log.path, function(x) x != -Inf)
  # compute fit at each time point for range of time
  Xt.sweep <- lapply(1:length(tp), function(t) # continuous time
    sapply(c.rng, function(c)
      cor(log.path[[t]][mask[[t]]],log(expm(-L.out*tp[t]*c) %*% Xo,base=10)[mask[[t]]])))
  c.opt <- sapply(Xt.sweep, function(x) c.rng[which.max(x)])
  print(c.opt)
  c <- mean(c.opt)
  return(c)
}

c.fit <- function(log.path,L.out,tp,ROI,c.rng){
  n.regions <- length(log.path[[1]])
  # initialize Xo vector, path at t=0
  Xo <- matrix(0,nrow=n.regions)
  # seed ROI with max path at month 1
  Xo[which(ROInames == ROI)] <- 1
  # exclusion mask... don't count regions with 0 path
  mask <- lapply(log.path, function(x) x != -Inf)
  # compute fit at each time point for range of time
  Xt.sweep <- lapply(1:length(tp), function(t) # continuous time
    sapply(c.rng, function(c)
      cor(log.path[[t]][mask[[t]]],log(expm(-L.out*tp[t]*c) %*% Xo,base=10)[mask[[t]]])))
  Xt.sweep <- Reduce('+',Xt.sweep) / length(tp) # mean fit for each tp
  c <- c.rng[which.max(Xt.sweep)]
  print(c)
  return(c)
}

c.fit.r <- function(log.path,L.out,tp,ROI,c.rng){
  n.regions <- length(log.path[[1]])
  # initialize Xo vector, path at t=0
  Xo <- matrix(0,nrow=n.regions)
  # seed ROI with max path at month 1
  Xo[which(ROInames == ROI)] <- 1
  # exclusion mask... don't count regions with 0 path
  mask <- lapply(log.path, function(x) x != -Inf)
  Xt.sweep <- lapply(1:length(tp), function(t) # continuous time
    sapply(c.rng, function(c)
      cor(log.path[[t]][mask[[t]]],log(expm(-L.out*tp[t]*c) %*% Xo,base=10)[mask[[t]]])))
  r <- sapply(Xt.sweep, function(x) max(x)) # best fit for each time point
  Xt.sweep <- Reduce('+',Xt.sweep) / length(tp) # mean fit for each time point
  c <- c.rng[which.max(Xt.sweep)]
  print(c)
  
  return(list(c=c,r=r))
}

get.L.out <- function(W){
  # stepwise prediction, i.e. from one month to the next
  # Where i is row element and j is column element
  # Wij is a connection from region i to region j
  # convention is the opposite, so without transposing W
  # I am capturing "retrograde" connectivity
  
  in.deg <- colSums(W)
  out.deg <- rowSums(W)
  L.out <- diag(x = out.deg) - W # outdegree laplacian
  return(L.out)
}

make.Xo <- function(ROI,ROInames){
  n.regions = length(ROInames)
  Xo <- matrix(0,nrow=n.regions)
  Xo[which(ROInames == ROI)] <- 1 # seed ROI with path
  return(Xo)
}

c.fit.r.old <- function(log.path,L.out,tp,ROI,c.rng){
  n.regions <- length(log.path[[1]])
  # initialize Xo vector, path at t=0
  Xo <- matrix(0,nrow=n.regions)
  # seed ROI with max path at month 1
  Xo[which(ROInames == ROI)] <- 1
  # exclusion mask... don't count regions with 0 path
  mask <- lapply(log.path, function(x) x != -Inf)
  Xt.sweep <- lapply(1:length(tp), function(t) # continuous time
    sapply(c.rng, function(c)
      cor(log.path[[t]][mask[[t]]],log(expm(-L.out*tp[t]*c) %*% Xo,base=10)[mask[[t]]])))
  c.opt <- sapply(Xt.sweep, function(x) c.rng[which.max(x)])
  print(c.opt)
  c <- mean(c.opt)
  r <- sapply(Xt.sweep, function(x) max(x))
  return(list(c=c,r=r))
}

c.fit.Syn <- function(log.path,L.out,Syn,tp,ROI,c.rng){
  n.regions <- length(log.path[[1]])
  # initialize Xo vector, path at t=0
  Xo <- matrix(0,nrow=n.regions)
  # seed ROI with max path at month 1
  Xo[which(ROInames == ROI)] <- max(exp(log.path[[1]]))
  # exclusion mask... don't count regions with 0 path
  mask <- lapply(log.path, function(x) x != -Inf)
  Xt.sweep <- lapply(1:length(tp), function(t) # continuous time
    sapply(c.rng, function(c)
      cor(log.path[[t]],log(expm(-L.out*tp[t]*c) %*% log(Xo,base=10)[mask[[t]]]))))
  c.opt <- sapply(Xt.sweep, function(x) c.rng[which.max(x)])
  print(c.opt)
  c <- mean(c.opt)
  return(c)
}

calc.vuln <- function(vulnerability,ROInames){
  # calculate vulnerability as the mean vulnerability value across time and hemispheres
  vulnerability <- 0.5*(vulnerability[grep('R ',as.character(ROInames)),] + vulnerability[grep('L ',as.character(ROInames)),])
  vulnerability <- rowMeans(vulnerability) # average across time
  vulnerability <- c(vulnerability,vulnerability)
  return(vulnerability)  
}

fisher.r.to.z <- function(r){
  r <- 0.5*(log(1+r) - log(1-r))
  return(r)
}