clustersplmm <- function(x,y,z,grp,time,lam1,lam2,nCluster,nonpen.b=1,nonpen.L=1,penalty.b=c("lasso","scad"),
                         penalty.L=c("lasso","scad"),standardize=TRUE,control=splmmControl()){
  
  
  penalty.b <- match.arg(penalty.b)
  penalty.L <- match.arg(penalty.L)
  
  # transform data
  if(is.data.frame(x)) x = as.matrix(x)
  if(is.data.frame(z)) z = as.matrix(z)
  
  # do some checks
  if (!is.matrix(x)) stop("x has to be a matrix or data frame")
  if (any(is.na(x))) stop("Missing values in x not allowed")
  
  if (any(is.na(y))) stop("Missing values in y not allowed")
  if (!is.numeric(y)) stop("y has to be of type 'numeric'")
  if (nrow(x)!=length(y)) stop("x and y have not correct dimensions")
  
  if (!is.matrix(z)) stop("z has to be a matrix")
  if (any(is.na(z))) stop("Missing values in z not allowed")
  
  if (any(x[,1]!=rep(1,dim(x)[[1]]))) stop("first column is not the intercept")
  
  if (length(levels(grp))==1) stop("Only one group. No covariance parameters!")
  
  if (!all(nonpen.b%in%1:dim(x)[[2]])) stop("Error with the argument nonpen for beta")
  if (!all(nonpen.L%in%1:dim(z)[[2]])) stop("Error with the argument nonpen for L")
  
  if (length(which(lam1<0))>0) stop("lam1 must be positive")
  if (length(which(lam2<0))>0) stop("lam2 must be positive")
  
  if (length(lam1)==1) lam1 <- rep(lam1, nCluster)
  if (length(lam2)==1) lam2 <- rep(lam2, nCluster)
  
  
  ##### Standardize covariates
  
  if (standardize)
  {
    xOr <- x
    meanx <- apply(x[,-1,drop = FALSE],2,mean)
    sdx <- apply(x[,-1,drop = FALSE],2,sd)
    x <- cbind(1,scale(x[,-1],center=meanx,scale=sdx))
    
    zOr <- z
    meanz <- apply(z[,-1,drop = FALSE],2,mean)
    sdz <- apply(z[,-1,drop = FALSE],2,sd)
    z <- cbind(1,scale(z[,-1,drop = FALSE],center=meanz,scale=sdz))
  }
  
  ##### allocate variables
  
  grp <- factor(grp)
  N <- length(levels(grp)) # N is the number of groups
  p <- dim(x)[[2]]         # p is the number of covariates
  q <- dim(z)[[2]]         # q is the number of random effects variables
  ntot <- length(y)        # ntot is the total number of observations
  Q <- q*(q+1)/2           # maximum number of variance components parameters
  lambda1 <- lam1*(N/nCluster)
  lambda2 <- lam2*(N/nCluster)
  
  ###### initiation
  ###### save the grouped information as components of a list
  long <- cbind.data.frame(grp,time,y)
  wide <- long %>% 
    spread(time, y)
  wide <- as.data.frame(wide)
  #cld.obj <- cld(wide,timeInData=2:ncol(wide),time=1:(ncol(wide)-1))
  #kml(cld.obj,nbClusters=nCluster,nbRedrawing=1)
  #myLetters <- letters[1:26]
  wide[is.na(wide)] <- 0
  ini.fit <- longclust::longclustEM(as.matrix(wide[2:ncol(wide)]),Gmin = nCluster, Gmax = nCluster)
  #cluster.label <- match(tolower(getClusters(cld.obj,nCluster)),myLetters)
  #cluster.label[which(is.na(cluster.label))]=1
  
  ##ini.fit <- kmeans(y,nCluster)
  #cluster.label <- sort(rep(1:nCluster, times=(ntot/nCluster)))
  yGrp.mixture <- xGrp.mixture <- zGrp.mixture <- ziGrp.mixture <- zIdGrp.mixture <- list()
  memb.prob <- vector()
  betaStart <- covStart <- sigmaStart <- parsStart <- LStart <- DStart <- VInvGrp <- list()
  membership <- matrix(ncol = N, nrow = nCluster)
  
  yGrp <- split(y,grp)
  xGrp <- split.data.frame(x,grp)
  zGrp <- split.data.frame(z,grp)
  zIdGrp <- mapply(ZIdentity,zGrp)
  
  ##ll1 <- 1/2*ntot*log(2*pi)
  ll1 <- list()
  ntot.mix <- rep(1,ntot)
  ntot.Grp <- split(ntot.mix,grp)
  ntot.mixture <-list()
  N.mixture <- list()
  
  fctStart <- vector()
  
  
  for (i in 1:nCluster) {
    #zi <- ifelse(ini.fit$cluster==i,1,0)
    #zi <- ifelse(cluster.label==i,1,0)
    zi <- round(ini.fit$zbest[,i])
    #memb.prob[i] <- sum(zi)/ntot
    memb.prob[i] <- sum(zi)/N
    #lambda1[i] <- lam1[i]*sum(zi)
    #lambda2[i] <- lam2[i]*sum(zi)
    
    
    ##ziGrp <- ziGrp.mixture[[i]] <- split(zi,grp)
    ##membership[i,] <- unlist(lapply(ziGrp, mean))
    ##y_i <- yGrp.mixture[[i]] <- mapply(multiplication,ziGrp,yGrp, SIMPLIFY = FALSE)
    ##x_i <- xGrp.mixture[[i]] <- mapply(multiplication,ziGrp,xGrp, SIMPLIFY = FALSE)
    ##z_i <- zGrp.mixture[[i]] <- mapply(multiplication,ziGrp,zGrp, SIMPLIFY = FALSE)
    ##zIdGrp.mixture[[i]] <- mapply(multiplication, ziGrp,zIdGrp, SIMPLIFY = FALSE)
    
    membership[i,] <- zi
    ntot.mixture[[i]] <- sum(unlist(mapply(multiplication,membership[i,],ntot.Grp, SIMPLIFY = FALSE)))
    N.mixture[[i]] <- sum(zi)
    ll1[[i]] <- 1/2*ntot.mixture[[i]]*log(2*pi)
    ntot.mixture[[i]] <- sum(zi)
    
    y_i <- yGrp.mixture[[i]] <- mapply(multiplication,membership[i,],yGrp, SIMPLIFY = FALSE)
    x_i <- xGrp.mixture[[i]] <- mapply(multiplication,membership[i,],xGrp, SIMPLIFY = FALSE)
    z_i <- zGrp.mixture[[i]] <- mapply(multiplication,membership[i,],zGrp, SIMPLIFY = FALSE)
    zIdGrp.mixture[[i]] <- mapply(multiplication, membership[i,],zIdGrp, SIMPLIFY = FALSE)
    
    ###### beta
    #init <- glmnet(x=zi*x,y=zi*y,lambda=lam1[i])
    #init <- glmnet(x=matrix(unlist(xGrp.mixture), ncol = p),unlist(yGrp.mixture),lambda = lam1[i])
    ###yi <- unlist(yGrp.mixture)
    ###xi <- matrix(unlist(xGrp.mixture), ncol = p)
    ###init <- lm(yi~xi[,-1])
    #betaStart[[i]] <- init$beta[,1]
    ###betaStart[[i]] <- init$coefficients
    ###betaStart[[i]][is.na(betaStart[[i]])] <- 0
    ##init <- optL1(yi,xi[,-1],model="linear",fold=10,trace=FALSE)
    ##betaStart[[i]] <- c(init$fullfit@unpenalized,init$fullfit@penalized)
    ###### covariance
    #covStart[[i]] <- covStartingValues(x_i,y_i,z_i,zIdGrp,betaStart[[i]],ntot.mixture[[i]],N.mixture[[i]])
    
    yi <- unlist(yGrp.mixture)
    xi <- matrix(unlist(xGrp.mixture), ncol = p)
    init <- optL1(yi,xi[,-1],model="linear",fold=10,trace=FALSE)
    betaStart[[i]] <- c(init$fullfit@unpenalized,init$fullfit@penalized)
    
    covStart[[i]] <- covStartingValues(x_i,y_i,z_i,zIdGrp,betaStart[[i]],ntot,N)
    sigmaStart[[i]] <- covStart[[i]]$sigma
    parsStart[[i]] <- vecli(covStart[[i]]$tau*diag(q))
    
    DStart[[i]] <- crossprod(triang(parsStart[[i]],q))
    LStart[[i]] <- chol(DStart[[i]])
    VInvGrp[[i]] <- mapply(VInv,Z=z_i,ZId=zIdGrp.mixture[[i]],MoreArgs=list(D=DStart[[i]],sigma=sigmaStart[[i]]))
    
  
    # --- Calculate objective function for the starting values ---
    # ------------------------------------------------------------
    
    fctStart[i] <- ObjFunction(xGroup=x_i,yGroup=y_i,LGroup=VInvGrp[[i]],b_nonpen=betaStart[[i]][-nonpen.b],L_nonpen=LStart[[i]][-nonpen.L,,drop=FALSE],
                               lambda1=lambda1[i],penalty_b=penalty.b,lambda2=lambda2[i],penalty_L=penalty.L,ll1=ll1[[i]],ntot = N.mixture[[i]]) 
    
    #fctStart[i] <- objFunction(xGroup=x_i,yGroup=y_i,zGroup=z_i,zIdGrp=zIdGrp.mixture[[i]], b=betaStart[[i]], L=LStart[[i]], sigma=sigmaStart[[i]],
    #                           lambda1=lambda1[i],penalty.b=penalty.b,lambda2=lambda2[i],nonpen.L=nonpen.L, nonpen.b=nonpen.b, penalty.L=penalty.L,ll1=ll1, ntot=1) 
  }
  
  
  
  # some necessary allocations:
  betaIter <- betaStart
  sigmaIter <- sigmaStart
  LIter <- LStart
  DIter <- DStart
  
  LvecIter <- lapply(LStart, function(x) x[lower.tri(x,diag = TRUE)])
  convPar <- max(unlist(lapply(betaIter, function(x) crossprod(x))))
  convCov <- max(mapply(function(x,y) crossprod(c(x,y)), sigmaStart, LvecIter))
  
  
  
  fctIter <- convFct <- fctStart
  hessian0 <- rep(0,p)
  mat0 <- matrix(0,ncol=p,nrow=N)
  covIter <- mapply(function(x,y) c(x,y), sigmaIter, LvecIter, SIMPLIFY = FALSE)
  
  ##### algorithm parameters
  
  stopped <- FALSE
  doAll <- FALSE
  converged <- 0
  counterIn <- 0
  counter <- 0       # counts the number of outer iterations
  
  convFct2 <- -10
  
  #while (counter<2) {
  
  
  while((counter<control$maxIter)&(convFct2<0|counter<1)&((convPar>control$tol|convFct>control$tol|convCov>control$tol|!doAll ))) {
    #while((counter<control$maxIter)&((convPar>control$tol|convFct>control$tol|convCov>control$tol|!doAll ))) {
    counter <- counter + 1 
    
    
    betaIterOld <- betaIter
    LIterOld <- LIter
    fctIterOld <- fctIter
    covIterOld <- covIter
    
    
    
    activeSet <- lapply(betaIter, function(x) which(x!=0))
    #activeSet <- lapply(activeSet, function(x) 1:p)
    #if ((length(activeSet)>min(p,ntot))&(lambda1>0)&(counter>2)) {stopped <- TRUE ; break}
    
    if (counterIn==0 | counterIn>control$number)
    {
      doAll <- TRUE
      activeSet <- lapply(activeSet, function(x) 1:p)
      counterIn <- 1    
    } else
    {
      doAll <- FALSE
      counterIn <- counterIn+1
    }
    
    
    for (i in 1:nCluster) {
      # --- optimization w.r.t the fixed effects vector beta ---
      # --------------------------------------------------------
      
      
      HessIter <- HessIterTrunc <- HessianMatrix(xGroup=xGrp.mixture[[i]],LGroup=VInvGrp[[i]],activeSet=activeSet[[i]],N=N,hessian=hessian0,mat=mat0[,activeSet[[i]],drop=FALSE])
      HessIter[activeSet[[i]]] <- pmin(pmax(HessIter[activeSet[[i]]],control$lower),control$upper)
      LxGrp <- as1(xGrp.mixture[[i]],VInvGrp[[i]],activeSet[[i]],N=N)
      ll2 <- nlogdet(LGroup=VInvGrp[[i]])
      
      for (j in activeSet[[i]])
      {
        cut1 <- as2(x=do.call("rbind",xGrp.mixture[[i]]),y=unlist(yGrp.mixture[[i]]),b=betaIter[[i]],j=j,activeSet=activeSet[[i]],group=grp,sGroup=LxGrp)
        JinNonpen <- j%in%nonpen.b
        
        # optimum can be calculated analytically
        if (HessIterTrunc[j]==HessIter[j])
        {
          if (JinNonpen) {betaIter[[i]][j] <- cut1/HessIter[j]} else {
            if(penalty.b=="scad"){
              scada = 3.7
              betaIter[[i]][j] <- SoftThreshold(cut1,lambda1[i])/(HessIter[j]*(1-1/scada))
              
            }else if(penalty.b=="lasso"){
              betaIter[[i]][j] <- SoftThreshold(cut1,lambda1[i])/HessIter[j]
            }
            
          }
        }else
          
          # optimimum is determined by the armijo rule
        {
          armijo <- ArmijoRule_b(xGroup=xGrp.mixture[[i]],yGroup=yGrp.mixture[[i]],LGroup=VInvGrp[[i]],b=betaIter[[i]],j=j,cut=cut1,HkOldJ=HessIterTrunc[j],
                                 HkJ=HessIter[j],JinNonpen=JinNonpen,lambda=lambda1[i],nonpen=nonpen.b,penalty=penalty.b,
                                 ll1=ll1[[i]],ll2=ll2,converged=converged,control=control)
          
          
          betaIter[[i]] <- armijo$b
          converged <- armijo$converged
          fctIter[i] <- armijo$fct
        }
        
      }
      betaIter[[i]][abs(betaIter[[i]])<0.05]=0
      
      # --- optimization w.r.t the variance components parameters ---
      # -------------------------------------------------------------
      
      # calculations before the covariance optimization
      activeSet[[i]] <- which(betaIter[[i]]!=0)
      resGrp <- ResAsSplit(x=do.call("rbind",xGrp.mixture[[i]]),y=unlist(yGrp.mixture[[i]]),b=betaIter[[i]],f=grp,activeset=activeSet[[i]])
      ll4 <- lambda1[i]*sum(abs(betaIter[[i]][-nonpen.b]))
      
      # optimization of L
      
      activeSet.L = which(rowSums(abs(LIter[[i]]))!=0)
      
      
      # calculate the hessian matrices for k in the activeSet
      
      D.grad = D_Gradient(xGroup=xGrp.mixture[[i]],zGroup=zGrp.mixture[[i]],LGroup=VInvGrp[[i]],yGroup=yGrp.mixture[[i]],b=betaIter[[i]],N=N,q=q)
      L.grad = t(LIter[[i]]%*%(D.grad+t(D.grad)))
      
      D.hessian = D_HessianMatrix(xGroup=xGrp.mixture[[i]],zGroup=zGrp.mixture[[i]],LGroup=VInvGrp[[i]],yGroup=yGrp.mixture[[i]],b=betaIter[[i]],N=N,q=q)
      D.hessian.submatrix = matsplitter(D.hessian,q)
      
      for (k in 1:q) {
        
        L.hessian.sub = sapply(D.hessian.submatrix[((k-1)*q+1):(k*q)],function(x) x[,k])
        L.hessian = diag(2*D.grad[k,k],q,q)+2*LIter[[i]]%*%(matrix(unlist(D.hessian.submatrix[[(k-1)*q+k]]),ncol = q,byrow = TRUE)+L.hessian.sub)%*%t(LIter[[i]])
        
        for (l in intersect(k:q,activeSet.L)) {
          
          L.lk.grad <- L.grad[l,k]
          L.lk.Hess <- L.hessian[l,l]
          L.lk.Hess <- min(max(L.lk.Hess,control$lower),control$upper)
          
          linNonpen <- l%in%nonpen.L
          
          ##armijo <- ArmijoRule_L(xGroup=xGrp.mixture[[i]][which(unlist(lapply(zGrp.mixture[[i]], function(x) sum(abs(x))!=0)))], yGroup=yGrp.mixture[[i]][which(unlist(lapply(zGrp.mixture[[i]], function(x) sum(abs(x))!=0)))], zGroup=zGrp.mixture[[i]][which(unlist(lapply(zGrp.mixture[[i]], function(x) sum(abs(x))!=0)))], L=LIter[[i]], l=l-1,k=k-1,grad=L.lk.grad,hessian=L.lk.Hess, 
          ##                       b=betaIter[[i]],sigma=sigmaIter[[i]],zIdGrp=zIdGrp.mixture[[i]], linNonpen=linNonpen, nonpen=nonpen.L-1, lambda=lambda2[i], penalty=penalty.L, 
          ##                       ll1=ll1[[i]], gamma=control$gamma, maxArmijo=control$maxArmijo, a_init=control$a_init, delta=control$delta, rho=control$rho, converged=converged, ntot=1)
          armijo <- armijoRule_L(xGroup=xGrp.mixture[[i]][which(unlist(lapply(zGrp.mixture[[i]], function(x) sum(abs(x))!=0)))], yGroup=yGrp.mixture[[i]][which(unlist(lapply(zGrp.mixture[[i]], function(x) sum(abs(x))!=0)))], zGroup=zGrp.mixture[[i]][which(unlist(lapply(zGrp.mixture[[i]], function(x) sum(abs(x))!=0)))], L=LIter[[i]], l=l,k=k,grad=L.lk.grad,hessian=L.lk.Hess, 
                                 b=betaIter[[i]],sigma=sigmaIter[[i]],zIdGrp=zIdGrp[which(unlist(lapply(zGrp.mixture[[i]], function(x) sum(abs(x))!=0)))], linNonpen=linNonpen,nonpen=nonpen.L, lambda=lambda2[i], penalty=penalty.L, 
                                 ll1=ll1[[i]],converged=converged,control=control, ntot=N.mixture[[i]])
          
          LIter[[i]] <- armijo$L
          converged <- armijo$converged
          fctIter[i] <- armijo$fct
          
        }
        
        
      }
      
      LIter[[i]][abs(LIter[[i]]) < 1e-2] <- 0
      
      DIter[[i]] = LIter[[i]]%*%t(LIter[[i]])
      LvecIter[[i]] = LIter[[i]][lower.tri(LIter[[i]],diag = TRUE)]
      
      
      # optimization of the error variance \sigma^2
      covParOpt <- MLsigma(zGroup=zGrp.mixture[[i]][which(unlist(lapply(zGrp.mixture[[i]], function(x) sum(abs(x))!=0)))],zIdGroup=zIdGrp.mixture[[i]][which(unlist(lapply(zGrp.mixture[[i]], function(x) sum(abs(x))!=0)))],resGroup=resGrp[which(unlist(lapply(zGrp.mixture[[i]], function(x) sum(abs(x))!=0)))],q=q,ll1=ll1[[i]],ll4=ll4,true.sigma=sigmaIter[[i]],D=DIter[[i]],
                           trace=control$trace,CovOpt="optimize",VarInt=control$VarInt)
      
      
      sigmaIter[[i]] <- covParOpt$sigma
      fctIter[i] <- covParOpt$fct
      
      VInvGrp[[i]] <- mapply(VInv,Z=zGrp.mixture[[i]],ZId=zIdGrp.mixture[[i]],MoreArgs=list(D=DIter[[i]],sigma=sigmaIter[[i]]))
      covIter[[i]] <- c(sigmaIter[[i]],LvecIter[[i]])
      
      
    }
    #betaIter <- lapply(betaIter,function(x) x[x<1e-4]=0)
    
    ##### clustering
    
    ##### update mixture membership
    
    ind.prob = matrix(nrow = nCluster, ncol = N)
    ll1.tmp <- 1/2*log(2*pi)
    
    for(i in 1:nCluster){
      VInvGrp.tmp <- mapply(VnotInv,Z=zGrp,ZId=zIdGrp,MoreArgs=list(D=DIter[[i]],sigma=sigmaIter[[i]]))
      cluster.den <- mapply(function(y,x,beta,V){
        dmvnorm(x=y,mean = x%*%beta, sigma = V)
      },y=yGrp, x=xGrp, V=VInvGrp.tmp, MoreArgs = list(beta=betaIter[[i]]))
      #cluster.den <- mapply(DensityFunction, x=xGrp, y=yGrp, V=VInvGrp.tmp, MoreArgs = list(b=betaIter[[i]], ll1=ll1.tmp))
      ind.prob[i,] <- cluster.den*memb.prob[i]
      #ind.prob[i,] <- exp(cluster.den)*memb.prob[i]
    }
    
    
    ind.prob[,which(colSums(ind.prob)==0)] <- 1/nCluster
    membership <- apply(ind.prob, 2, function(x) x/sum(x))
    memb.prob <- rowSums(membership)/sum(membership)
    
    ##### update x,y,z
    ntot.mixture[[i]] <- sum(unlist(mapply(multiplication,membership[i,],ntot.Grp, SIMPLIFY = FALSE)))
    N.mixture[[i]] <- sum(zi)
    ll1[[i]] <- 1/2*ntot.mixture[[i]]*log(2*pi)
    for (i in 1:nCluster) {
      yGrp.mixture[[i]] <- mapply(multiplication,membership[i,],yGrp, SIMPLIFY = FALSE)
      xGrp.mixture[[i]] <- mapply(multiplication,membership[i,],xGrp, SIMPLIFY = FALSE)
      zGrp.mixture[[i]] <- mapply(multiplication,membership[i,],zGrp, SIMPLIFY = FALSE)
      zIdGrp.mixture[[i]] <- mapply(multiplication, membership[i,],zIdGrp, SIMPLIFY = FALSE)
    }
    
    #lambda1 <- lam1*(rowSums(membership))
    #lambda2 <- lam2*(rowSums(membership))
    
    ##### reorder cluster
    
    
    
    
    # --- check convergence ---
    
    convPar <- max(mapply(function(x,y) sqrt(crossprod(x-y))/(1+sqrt(crossprod(x))), x=betaIter, y=betaIterOld))
    #convPar <- max(mapply(function(x,y) sqrt(crossprod(x-y)),x=betaIter, y=betaIterOld))
    convFct <- abs((sum(fctIterOld)-sum(fctIter))/(1+abs(sum(fctIter))))
    #convFct <- abs((sum(fctIterOld)-sum(fctIter)))
    convFct2 <- sum(fctIter) - sum(fctIterOld)
    #convCov <- max(mapply(function(x,y) sqrt(crossprod(x-y))/(1+sqrt(crossprod(x))), x=covIter, y=covIterOld))
    convCov <- max(mapply(function(x,y) sqrt(crossprod(x-y)), x=covIter, y=covIterOld))
    
    if ((convPar <= control$tol) & (convFct <= control$tol) & (convCov <= control$tol)) counterIn <- 0
    
  }
  
  
  if (standardize)
  {
    betaIter <- lapply(betaIter, function(x) {
      x[-1] <- x[-1]/sdx
      x[1] <- x[1] - sum(meanx*x[-1])
      return(x)
    })
    
    x <- xOr
    xGrp <- split.data.frame(xOr,grp)
    
    z <- zOr
    zGrp <- split.data.frame(zOr,grp)
  }
  
  
  # --- summary information ---
  # ---------------------------
  ntot.mix <- rep(1,ntot)
  ntot.Grp <- split(ntot.mix,grp)
  ntot.mixture <-list()
  
  for (i in 1:nCluster) {
    yGrp.mixture[[i]] <- mapply(multiplication,membership[i,],yGrp, SIMPLIFY = FALSE)
    xGrp.mixture[[i]] <- mapply(multiplication,membership[i,],xGrp, SIMPLIFY = FALSE)
    zGrp.mixture[[i]] <- mapply(multiplication,membership[i,],zGrp, SIMPLIFY = FALSE)
    ntot.mixture[[i]] <- mapply(multiplication,membership[i,],ntot.Grp, SIMPLIFY = FALSE)
    zIdGrp.mixture[[i]] <- mapply(multiplication, membership[i,],zIdGrp, SIMPLIFY = FALSE)
    VInvGrp[[i]] <- mapply(VInv,Z=zGrp.mixture[[i]],ZId=zIdGrp.mixture[[i]],MoreArgs=list(D=DIter[[i]],sigma=sigmaIter[[i]]))
  }
  
  ntot.mixture <- matrix(unlist(ntot.mixture), nrow=nCluster, byrow=T)
  ntot.mixture <- apply(ntot.mixture, 1, sum)
  
  N.mixture <- apply(membership,1,sum)
  ntot.mixture <- lapply(ntot.mixture, function(x) sum(x))
  
  ntot.mixture <- lapply(ntot.mixture, function(x){
    x[x<1] <- 1
    return(x)
  })
  
  N.mixture[N.mixture<1] <- 1
  D <- DIter
  #npar <- sum(unlist(lapply(betaIter, function(x) sum(x!=0)))) + sum(unlist(lapply(D, function(x) sum(diag(x)!=0)))) + 1
  npar<- sum(unlist(lapply(betaIter, function(x) sum(x!=0)))) + length(c(unlist(sigmaIter),unlist(LvecIter)))
  logLik <- sum(mapply(MLloglik, xGroup=xGrp.mixture,yGroup=yGrp.mixture,LGroup=VInvGrp,b=betaIter,activeSet=activeSet, ntot=ntot.mixture,MoreArgs = list(N=N)))
  deviance <- -2*logLik
  aic <- -2* logLik + 2*npar
  bic <- -2* logLik + log(ntot)*npar
  
  p <- sum(unlist(lapply(betaIter, function(x) sum(x!=0))))
  q <- sum(unlist(lapply(D, function(x) sum(diag(x)!=0))))
  bbic <- -2*logLik + max(1,log(log(p+q)))*log(ntot)*npar
  ebic <- -2*logLik + (log(ntot)+2*log(p+q))*npar
  
  if (converged>0) cat("Algorithm does not properly converge.","\n")
  if (stopped) {cat("|activeSet|>=min(p,ntot): Increase lambda or set stopSat=FALSE.","\n")
    ; sigmaIter <- LvecIter <- nlogLik <- aic <- bic <- NA ; betaIter <- rep(NA,p) ; bi <- fitted <- residuals <- NULL}
  
  out <- list(data=list(x=x,y=y,z=z,grp=grp),membership=membership, coefInit=list(betaStart=betaStart,parsStart=parsStart,sigmaStart=sigmaStart),penalty.b=penalty.b,penalty.L=penalty.L,
              nonpen.b=nonpen.b,nonpen.L=nonpen.L,lambda1=lambda1,lambda2=lambda2,sigma=sigmaIter,Lvec=LvecIter,coefficients=betaIter,D=D,converged=converged,logLik=logLik,npar=npar,deviance=deviance,
              aic=aic,bic=bic,bbic=bbic,ebic=ebic,counter=counter,control=control,call=match.call(),
              stopped=stopped,objective=fctIter)
  
  out
  structure(out,class="splmm")
}