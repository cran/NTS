# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#' Estimation of Autoregressive Conditional Mean Models
#'
#' Estimation of autoregressive conditional mean models with exogenous variables.
#' @param y time series of counts.
#' @param order the order of ACM model.
#' @param X matrix of exogenous variables.
#' @param cond.dist conditional distributions. "po" for Poisson, "nb" for negative binomial, "dp" for double Poisson.
#' @param ini initial parameter estimates designed for use in "nb" and "dp".
#' @return ACMx returns a list with components:
#' \item{data}{time series.}
#' \item{X}{matrix of exogenous variables.}
#' \item{estimates}{estimated values.}
#' \item{residuals}{residuals.}
#' \item{sresi}{standardized residuals.}
#' @examples
#' x=rnorm(1000)*0.1
#' y=matrix(0,1000,1)
#' y[1]=2
#' lambda=matrix(0,1000,1)
#' for (i in 2:1000){
#' 	lambda[i]=2+0.2*y[i-1]/exp(x[i-1])+0.5*lambda[i-1]
#' 	y[i]=rpois(1,exp(x[i])*lambda[i])
#' }
#' ACMx(y,order=c(1,1),x,"po")
#' @import stats
#' @export
"ACMx" <- function(y,order=c(1,1),X=NULL,cond.dist="po",ini=NULL){
  beta=NULL; k=0; withX=FALSE; nT=length(y)
  if(!is.null(X)){
    withX=TRUE
    if(!is.matrix(X))X=as.matrix(X)
    T1=dim(X)[1]
    if(nT > T1)nT=T1
    if(nT < T1){T1=nT; X=X[1:T1,]}
    ##
    k=dim(X)[2]
    ### Preliminary estimate of regression coefficients
    m1=glm(y~X,family=poisson)
    m11=summary(m1)
    beta=m11$coefficients[-1,1]; se.beta=m11$coefficients[-1,2]
    glmresi=m1$residuals
  }
  ### obtain initial estimate of the constant
  mm=glm(y[2:nT]~y[1:(nT-1)],family=poisson)
  m22=summary(mm)
  ome=m22$coefficients[1,1]; se.ome=m22$coefficients[1,2]
  ###
  p=order[1]; q=order[2]; S=1.0*10^(-6); params=NULL; loB=NULL; upB=NULL
  ###
  ###### Poisson model
  if(cond.dist=="po"){
    if(withX){
      params=c(params,beta=beta); loB=c(loB,beta=beta-2*abs(beta)); upB=c(upB,beta=beta+2*abs(beta))
    }
    params=c(params,omega=ome); loB=c(loB,omega=S); upB=c(upB,omega=ome+10*ome)
    if(p > 0){
      a1=rep(0.05,p); params=c(params,alpha=a1); loB=c(loB,alpha=rep(S,p)); upB=c(upB,alpha=rep(0.5,p))
    }
    if(q > 0){
      b1=rep(0.5,q)
      params=c(params,gamma=b1); loB=c(loB,gamma=rep(S,q)); upB=c(upB,gamma=rep(1-S,q))
    }
    ##mm=optim(params,poX,method="L-BFGS-B",hessian=T,lower=loB,upper=upB)
    #
    cat("Initial estimates: ",params,"\n")
    cat("loB: ",loB,"\n")
    cat("upB: ",upB,"\n")
    fit=nlminb(start = params, objective= poX, lower=loB, upper=upB, PCAxY=y, PCAxX=X, PCAxOrd=order,
               control=list(rel.tol=1e-6))
    #control=list(trace=3,rel.tol=1e-6))

    epsilon = 0.0001 * fit$par
    npar=length(params)
    Hessian = matrix(0, ncol = npar, nrow = npar)
    for (i in 1:npar) {
      for (j in 1:npar) {
        x1 = x2 = x3 = x4  = fit$par
        x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
        x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
        x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
        x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
        Hessian[i, j] = (poX(x1,PCAxY=y,PCAxX=X,PCAxOrd=order)-poX(x2,PCAxY=y,PCAxX=X,PCAxOrd=order)
                         -poX(x3,PCAxY=y,PCAxX=X,PCAxOrd=order)+poX(x4,PCAxY=y,PCAxX=X,PCAxOrd=order))/
          (4*epsilon[i]*epsilon[j])
      }
    }
    cat("Maximized log-likehood: ",-poX(fit$par,PCAxY=y,PCAxX=X,PCAxOrd=order),"\n")
    # Step 6: Create and Print Summary Report:
    se.coef = sqrt(diag(solve(Hessian)))
    tval = fit$par/se.coef
    matcoef = cbind(fit$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
    dimnames(matcoef) = list(names(tval), c(" Estimate",
                                            " Std. Error", " t value", "Pr(>|t|)"))
    cat("\nCoefficient(s):\n")
    printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
    ###### Compute the residuals
    est=fit$par; ist=1
    bb=rep(1,length(y)); y1=y
    if(withX){beta=est[1:k]; bb=exp(X%*%matrix(beta,k,1))}#; y1=y/bb}
    icnt=k+1
    ome=est[icnt]
    p=order[1]; q=order[2]; nT=length(y)
    if(p > 0){a1=est[(icnt+1):(icnt+p)]; icnt=icnt+p
    nobe=nT-p; rate=rep(ome,nobe); ist=p+1
    for (i in 1:p){
      rate=rate+a1[i]*y1[(ist-i):(nT-i)]
    }
    }
    r1=ome
    if(q > 0){g1=est[(icnt+1):(icnt+q)]
    r1=filter(rate,g1,"r",init=rep(mean(y/bb),q))
    }
    rate=bb[ist:nT]*r1
    resi=y[ist:nT]-rate
    sresi=resi/sqrt(rate)
  }
  #### Negative binomial ######################################################
  if(cond.dist=="nb"){
    ##
    if(is.null(ini)){
      if(withX){
        params=c(params,beta=rep(0.4,k)); loB=c(loB,beta=rep(-3,k)); upB=c(upB,beta=rep(2,k))
      }
      params=c(params,omega=ome); loB=c(loB,omega=S); upB=c(upB,omega=ome+6*se.ome)
      if(p > 0){
        a1=rep(0.05,p); params=c(params,alpha=a1); loB=c(loB,alpha=rep(S,p)); upB=c(upB,alpha=rep(0.5,p))
      }
      if(q > 0){
        b1=rep(0.7,q)
        params=c(params,gamma=b1); loB=c(loB,gamma=rep(S,q)); upB=c(upB,gamma=rep(1-S,q))
      }
      ## given initial estimates
    }
    else{
      se.ini=abs(ini/10); jst= k
      if(withX){
        params=c(params,beta=ini[1:k]); loB=c(loB,beta=ini[1:k]-4*se.ini[1:k]); upB=c(upB,beta=ini[1:k]+4*se.ini[1:k])
      }
      jst=k+1
      params=c(params,omega=ini[jst]); loB=c(loB,omega=ini[jst]-3*se.ini[jst]); upB=c(upB,omega=ini[jst]+3*se.ini[jst])
      if(p > 0){a1=ini[(jst+1):(jst+p)]; params=c(params,alpha=a1); loB=c(loB,alpha=rep(S,p))
      upB=c(upB,alpha=a1+3*se.ini[(jst+1):(jst+p)]); jst=jst+p
      }
      if(q > 0){
        b1=ini[(jst+1):(jst+q)]
        params=c(params,gamma=b1); loB=c(loB,gamma=rep(S,q)); upB=c(upB,gamma=rep(1-S,q))
      }
    }
    ### size of the negative binomial parameter
    meanY=mean(y); varY=var(y); th1=meanY^2/(varY-meanY)*2; if(th1 < 0)th1=meanY*0.1
    params=c(params,theta=th1); loB=c(loB,theta=0.5); upB=c(upB,theta=meanY*1.5)
    #
    cat("initial estimates: ",params,"\n")
    cat("loB: ",loB,"\n")
    cat("upB: ",upB,"\n")
    fit=nlminb(start = params, objective= nbiX, lower=loB, upper=upB, PCAxY=y,PCAxX=X,PCAxOrd=order,
               control=list(rel.tol=1e-6))
    #control=list(trace=3,rel.tol=1e-6))

    epsilon = 0.0001 * fit$par
    npar=length(params)
    Hessian = matrix(0, ncol = npar, nrow = npar)
    for (i in 1:npar) {
      for (j in 1:npar) {
        x1 = x2 = x3 = x4  = fit$par
        x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
        x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
        x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
        x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
        Hessian[i, j] = (nbiX(x1,PCAxY=y,PCAxX=X,PCAxOrd=order)-nbiX(x2,PCAxY=y,PCAxX=X,PCAxOrd=order)
                         -nbiX(x3,PCAxY=y,PCAxX=X,PCAxOrd=order)+nbiX(x4,PCAxY=y,PCAxX=X,PCAxOrd=order))/
          (4*epsilon[i]*epsilon[j])
      }
    }
    cat("Maximized log-likehood: ",-nbiX(fit$par,PCAxY=y,PCAxX=X,PCAxOrd=order),"\n")
    # Step 6: Create and Print Summary Report:
    se.coef = sqrt(diag(solve(Hessian)))
    tval = fit$par/se.coef
    matcoef = cbind(fit$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
    dimnames(matcoef) = list(names(tval), c(" Estimate",
                                            " Std. Error", " t value", "Pr(>|t|)"))
    cat("\nCoefficient(s):\n")
    printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
    #### Compute the residuals
    est=fit$par
    bb=rep(1,length(y)); y1=y
    if(withX){beta=est[1:k]; bb=exp(X%*%matrix(beta,k,1))}#; y1=y/bb}
    icnt=k+1
    ome=est[icnt]
    p=order[1]; q=order[2]; nT=length(y)
    if(p > 0){a1=est[(icnt+1):(icnt+p)]; icnt=icnt+p
    nobe=nT-p; rate=rep(ome,nobe); ist=p+1
    for (i in 1:p){
      rate=rate+a1[i]*y1[(ist-i):(nT-i)]
    }
    }
    if(q > 0){g1=est[(icnt+1):(icnt+q)]; icnt=icnt+q
    r1=filter(rate,g1,"r",init=rep(mean(y/bb),q))
    }
    rate=bb[ist:nT]*r1
    resi=y[ist:nT]-rate
    theta=est[icnt+1]
    pb=theta/(theta+rate)
    v1=theta*(1-pb)/(pb^2)
    sresi=resi/sqrt(v1)

  }
  ### Double Poisson model ###################################
  if(cond.dist=="dp"){
    ##
    if(withX){
      params=c(params,beta=beta); loB=c(loB,beta=beta-2*abs(beta)); upB=c(upB,beta=beta+2*abs(beta))
    }
    params=c(params,omega=ome*0.3); loB=c(loB,omega=S); upB=c(upB,omega=ome+4*se.ome)
    if(p > 0){
      a1=rep(0.05,p); params=c(params,alpha=a1); loB=c(loB,alpha=rep(S,p)); upB=c(upB,alpha=rep(1-S,p))
    }
    if(q > 0){
      b1=rep(0.6,q)
      params=c(params,gamma=b1); loB=c(loB,gamma=rep(S,q)); upB=c(upB,gamma=rep(1-S,q))
    }
    ### size of the negative binomial parameter
    meanY=mean(y); varY = var(y); t1=meanY/varY
    params=c(params,theta=t1); loB=c(loB,theta=S); upB=c(upB,theta=2)
    #
    cat("initial estimates: ",params,"\n")
    cat("loB: ",loB,"\n")
    cat("upB: ",upB,"\n")
    fit=nlminb(start = params, objective= dpX, lower=loB, upper=upB,PCAxY=y,PCAxX=X,PCAxOrd=order,
               control=list(rel.tol=1e-6))
    #control=list(trace=3,rel.tol=1e-6))

    epsilon = 0.0001 * fit$par
    npar=length(params)
    Hessian = matrix(0, ncol = npar, nrow = npar)
    for (i in 1:npar) {
      for (j in 1:npar) {
        x1 = x2 = x3 = x4  = fit$par
        x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
        x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
        x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
        x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
        Hessian[i, j] = (dpX(x1,PCAxY=y,PCAxX=X,PCAxOrd=order)-dpX(x2,PCAxY=y,PCAxX=X,PCAxOrd=order)
                         -dpX(x3,PCAxY=y,PCAxX=X,PCAxOrd=order)+dpX(x4,PCAxY=y,PCAxX=X,PCAxOrd=order))/
          (4*epsilon[i]*epsilon[j])
      }
    }
    cat("Maximized log-likehood: ",-dpX(fit$par,PCAxY=y,PCAxX=X,PCAxOrd=order),"\n")
    # Step 6: Create and Print Summary Report:
    se.coef = sqrt(diag(solve(Hessian)))
    tval = fit$par/se.coef
    matcoef = cbind(fit$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
    dimnames(matcoef) = list(names(tval), c(" Estimate",
                                            " Std. Error", " t value", "Pr(>|t|)"))
    cat("\nCoefficient(s):\n")
    printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
    #### Compute the residuals
    est=fit$par
    bb=rep(1,length(y)); y1=y
    if(withX){beta=est[1:k]; bb=exp(X%*%matrix(beta,k,1))}#; y1=y/bb}
    icnt=k+1
    ome=est[icnt]
    p=order[1]; q=order[2]; nT=length(y)
    if(p > 0){a1=est[(icnt+1):(icnt+p)]; icnt=icnt+p
    nobe=nT-p; rate=rep(ome,nobe); ist=p+1
    for (i in 1:p){
      rate=rate+a1[i]*y1[(ist-i):(nT-i)]
    }
    }
    if(q > 0){g1=est[(icnt+1):(icnt+q)]; icnt=icnt+q
    r1=filter(rate,g1,"r",init=rep(mean(y/bb),q))
    }
    theta=est[icnt+1]
    rate=bb[ist:nT]*r1
    resi=y[ist:nT]-rate
    v1=rate/theta
    sresi=resi/sqrt(rate)
  }

  #### end of the program
  ACMx <- list(data=y,X=X,estimates=est,residuals=resi,sresi=sresi)
}



"nbiX" <- function(par,PCAxY,PCAxX,PCAxOrd){
  ## compute the log-likelihood function of negative binomial distribution
  y <- PCAxY; X <- PCAxX; order <- PCAxOrd
  withX = F; k = 0
  #
  if(!is.null(X)){
    withX=T; k=dim(X)[2]
  }
  bb=rep(1,length(y)); y1=y
  if(withX){beta=par[1:k]; bb=exp(X%*%matrix(beta,k,1))}#; y1=y/bb}
  icnt=k+1
  ome=par[icnt]
  p=order[1]; q=order[2]; nT=length(y)
  if(p > 0){a1=par[(icnt+1):(icnt+p)]; icnt=icnt+p
  nobe=nT-p; rate=rep(ome,nobe); ist=p+1
  for (i in 1:p){
    rate=rate+a1[i]*y1[(ist-i):(nT-i)]
  }
  }
  if(q > 0){g1=par[(icnt+1):(icnt+q)]; icnt=icnt+q
  r1=filter(rate,g1,"r",init=rep(mean(y/bb),q))
  }
  theta=par[icnt+1]
  rate=bb[ist:nT]*r1
  pb=theta/(theta+rate)
  lnNBi=dnbinom(y[ist:nT],size=theta,prob=pb,log=T)
  nbiX <- -sum(lnNBi)
}



"poX" <- function(par,PCAxY,PCAxX,PCAxOrd){
  y <- PCAxY; X <- PCAxX; order <- PCAxOrd
  withX = F; k=0
  ##
  if(!is.null(X)){
    withX=T; k=dim(X)[2]
  }
  bb=rep(1,length(y)); y1=y
  if(withX){beta=par[1:k]; bb=exp(X%*%matrix(beta,k,1))}#; y1=y/bb}
  icnt=k+1; ist=1; nT=length(y); nobe=nT
  ome=par[icnt]; rate=rep(ome,nobe)
  p=order[1]; q=order[2]
  if(p > 0){a1=par[(icnt+1):(icnt+p)]; icnt=icnt+p
  nobe=nT-p; rate=rep(ome,nobe); ist=p+1
  for (i in 1:p){
    rate=rate+a1[i]*y1[(ist-i):(nT-i)]
  }
  }
  r1=rate
  if(q > 0){g1=par[(icnt+1):(icnt+q)]
  r1=filter(rate,g1,"r",init=rep(mean(y/bb),q))
  }

  rate=bb[ist:nT]*r1
  lnPoi=dpois(y[ist:nT],rate,log=T)
  poX <- -sum(lnPoi)
}




"dpX" <- function(par,PCAxY,PCAxX,PCAxOrd){
  # compute the log-likelihood function of a double Poisson distribution
  y <- PCAxY; X <- PCAxX; order <- PCAxOrd
  withX=F; k = 0
  #
  if(!is.null(X)){
    withX=TRUE; k=dim(X)[2]
  }
  bb=rep(1,length(y)); y1=y
  if(withX){beta=par[1:k]; bb=exp(X%*%matrix(beta,k,1))}#; y1=y/bb}
  icnt=k+1
  ome=par[icnt]
  p=order[1]; q=order[2]; nT=length(y)
  if(p > 0){a1=par[(icnt+1):(icnt+p)]; icnt=icnt+p
  nobe=nT-p; rate=rep(ome,nobe); ist=p+1
  for (i in 1:p){
    rate=rate+a1[i]*y1[(ist-i):(nT-i)]
  }
  }
  #plot(rate,type='l')

  if(q > 0){g1=par[(icnt+1):(icnt+q)]; icnt=icnt+q
  r1=filter(rate,g1,"r",init=rep(mean(y/bb),q))
  }
  theta=par[icnt+1]
  rate=bb[ist:nT]*r1
  yy=y[ist:nT]
  rate1 = rate*theta
  cinv=1+(1-theta)/(12*rate1)*(1+1/rate1)
  lcnt=-log(cinv)
  lpd=lcnt+0.5*log(theta)-rate1
  idx=c(1:nobe)[yy > 0]
  tp=length(idx)
  d1=apply(matrix(yy[idx],tp,1),1,factorialOwn)
  lpd[idx]=lpd[idx]-d1-yy[idx]+yy[idx]*log(yy[idx])+theta*yy[idx]*(1+log(rate[idx])-log(yy[idx]))
  #plot(lpd,type='l')
  #cat("neg-like: ",-sum(lpd),"\n")

  dpX <- -sum(lpd)
}



"factorialOwn" <- function(n,log=T){
  x=c(1:n)
  if(log){
    x=log(x)
    y=cumsum(x)
  }
  else{
    y=cumprod(x)
  }
  y[n]
}



#' Generate a CFAR(1) Process
#'
#' Generate a convolutional functional autoregressive process with order 1.
#' @param tmax length of time.
#' @param rho parameter for O-U process (noise process).
#' @param phi_func convolutional function. Default is density function of normal distribution with mean 0 and standard deviation 0.1.
#' @param grid the number of grid points used to construct the functional time series. Default is 1000.
#' @param sigma the standard deviation of O-U process. Default is 1.
#' @param ini the burn-in period.
#' @references
#' Liu, X., Xiao, H., and Chen, R. (2016) Convolutional autoregressive models for functional time series. \emph{Journal of Econometrics}, 194, 263-282.
#' @return The function returns a list with components:
#' \item{cfar1}{a tmax-by-(grid+1) matrix following a CFAR(1) process.}
#' \item{epsilon}{the innovation at time tmax.}
#' @examples
#' phi_func= function(x)
#' {
#'  	return(dnorm(x,mean=0,sd=0.1))
#' }
#' y=g_cfar1(1000,5,phi_func,grid=1000,sigma=1,ini=100)
#' @export
g_cfar1 <- function(tmax=1001,rho=5,phi_func=NULL,grid=1000,sigma=1,ini=100){
  #################################
  ###Simulate CFAR(1) processes
  #############################
  ###parameter setting
  #rho	### parameter for error process epsilon_t, O-U process
  #sigma is the standard deviation of OU-process
  #tmax	### length of time
  #iter	### number of replications in the simulation, iter=1
  #grid	### the number of grid points used to construct the functional time series X_t and epsilon_t
  # phi_func:   User supplied convolution function phi(.). The following default is used if not specified.
  if(is.null(phi_func)){
    phi_func= function(x){
      return(dnorm(x,mean=0,sd=0.1))
    }
  }
  if(is.null(grid)){
    grid=1000
  }
  if(is.null(sigma)){
    sigma=1
  }
  if(is.null(ini)){
    ini=100
  }
  #########################################################################################################
  #### OUTPUTS
  #######################
  #### A matrix f_grid is generated, with (tmax) rows and (grid+1) columns.
  ########################################################################################################

  #################################################################################################
  x_grid<- seq(0, 1, by= 1/grid); ### the grid points for input variables of functional time series in [0,1]
  b_grid<- seq(-1, 1, by=1/grid);	### the grid points for input variables of b-spline functions in [-1,1]
  eps_grid <- matrix(0, 1, grid+1);	### the grid points for input variables of error process in [0,1]
  fi <- exp(-rho/grid);		### the correlation of two adjacent points for an O-U process
  phi_b <- phi_func(b_grid)		### phi(s) when s in [-1,1]
  ### the convolutional function values phi(s-u), for different s in [0,1]
  phi_matrix <- matrix(0, (grid+1), (grid+1));
  for (i in 1:(grid+1)){
    phi_matrix[i,]= phi_b[seq((i+grid),i,by=-1)]/(grid+1)
  }

  ##############################
  ### Generate data
  ### functional time series f_grid is X_t in the paper
  #### A matrix f_grid is generated, with (tmax*iter) rows and (grid+1) columns.
  #### It contains iter CFAR(1) processes. The i-th process is in rows (i-1)*tmax+1:tmax.

  f_grid=matrix(0,tmax+ini,grid+1)
  i=1
  eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt((1-fi^2))*rnorm(n,0,1))  ###error process
  f_grid[1,]= eps_grid;
  for (i in 2:(tmax+ini)){
    eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt(1-fi^2)*rnorm(n,0,1))
    f_grid[i,]= phi_matrix%*%f_grid[i-1,]+ eps_grid;
  }
  g_cfar1 <- list(cfar1=f_grid[(ini+1):(tmax+ini),]*sigma,epsilon=eps_grid*sigma)
  return(g_cfar1)
}


#' Generate a CFAR(2) Process
#'
#' Generate a convolutional functional autoregressive process with order 2.
#' @param tmax length of time.
#' @param rho parameter for O-U process (noise process).
#' @param phi_func1 the first convolutional function. Default is 0.5*x^2+0.5*x+0.13.
#' @param phi_func2 the second convolutional function. Default is 0.7*x^4-0.1*x^3-0.15*x.
#' @param grid the number of grid points used to construct the functional time series. Default is 1000.
#' @param sigma the standard deviation of O-U process. Default is 1.
#' @param ini the burn-in period.
#' @references
#' Liu, X., Xiao, H., and Chen, R. (2016) Convolutional autoregressive models for functional time series. \emph{Journal of Econometrics}, 194, 263-282.
#' @return The function returns a list with components:
#' \item{cfar2}{a tmax-by-(grid+1) matrix following a CFAR(1) process.}
#' \item{epsilon}{the innovation at time tmax.}
#' @examples
#' phi_func1= function(x){
#'	return(0.5*x^2+0.5*x+0.13)
#' }
#' phi_func2= function(x){
#'	return(0.7*x^4-0.1*x^3-0.15*x)
#' }
#' y=g_cfar2(1000,5,phi_func1,phi_func2,grid=1000,sigma=1,ini=100)
#' @export
g_cfar2 <- function(tmax=1001,rho=5,phi_func1=NULL, phi_func2=NULL,grid=1000,sigma=1,ini=100){
  #################################
  ###Simulate CFAR(2) processes

  #########################################################################################
  #### OUTPUTS
  #######################
  #### A matrix f_grid is generated, with (tmax) rows and (grid+1) columns.
  ###########################################################################################

  #############################
  ###parameter setting
  ## rho=5	### parameter for error process epsilon_t, O-U process
  ##tmax= 1001;	### length of time
  ##iter= 100;	### number of replications in the simulation
  ##grid=1000;	### the number of grid points used to construct the functional time series X_t and epsilon_t
  ## phi_func1 and phi_func2: User specified convolution functions. The following defaults are used if not given
  if(is.null(phi_func1)){
    phi_func1= function(x){
      return(0.5*x^2+0.5*x+0.13)
    }
  }
  if(is.null(phi_func2)){
    phi_func2= function(x){
      return(0.7*x^4-0.1*x^3-0.15*x)
    }
  }
  if(is.null(grid)){
    grid=1000
  }
  if(is.null(sigma)){
    sigma=1
  }
  if(is.null(ini)){
    ini=100
  }
  #######################################################################################################
  x_grid= seq(0, 1, by= 1/grid);	### the grid points for input variables of functional time series in [0,1]
  b_grid= seq(-1, 1, by=1/grid);	### the grid points for input variables of b-spline functions in [-1,1]
  eps_grid= matrix(0, 1, grid+1);	### the grid points for input variables of error process in [0,1]
  fi = exp(-rho/grid);		### the correlation of two adjacent points for an O-U process
  phi_b1= phi_func1(b_grid)	### phi_1(s) when s in [-1,1]
  phi_b2= phi_func2(b_grid)	### phi_2(s) when s in [-1,1]

  ### the convolutional function values phi_1(s-u) and phi_2(s-u), for different s in [0,1]
  phi_matrix1= matrix(0, (grid+1), (grid+1));
  for (i in 1:(grid+1)){
    phi_matrix1[i,]= phi_b1[seq((i+grid),i,by=-1)]/(grid+1)
  }
  phi_matrix2= matrix(0, (grid+1), (grid+1));
  for (i in 1:(grid+1)){
    phi_matrix2[i,]= phi_b2[seq((i+grid),i,by=-1)]/(grid+1)
  }
  ##############################
  ### Generate data
  ### functional time series f_grid is X_t in the paper
  #### A matrix f_grid is generated, with (tmax*iter) rows and (grid+1) columns.
  ##  It contains iter CFAR(2) processes. The i-th process is in rows (i-1)*tmax+1:tmax.

  f_grid=matrix(0,tmax+ini,grid+1)
  i=1
  eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt(1-fi^2)*rnorm(n,0,1))
  f_grid[1,]= eps_grid;
  i=2
  eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt(1-fi^2)*rnorm(n,0,1))
  f_grid[i,]= phi_matrix1%*%as.matrix(f_grid[i-1,])+ eps_grid;
  for (i in 3:tmax){
    eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt(1-fi^2)*rnorm(n,0,1))
    f_grid[i,]= phi_matrix1%*%as.matrix(f_grid[i-1,])+ phi_matrix2%*%as.matrix(f_grid[i-2,])+ eps_grid;
  }
  ### The last iteration of epsilon_t is returned.
  g_cfar2 <- list(cfar2=f_grid*sigma,epsilon=eps_grid*sigma)
}


#' Generate a CFAR Process
#'
#' Generate a convolutional functional autoregressive process.
#' @param tmax length of time.
#' @param rho parameter for O-U process (noise process).
#' @param phi_list the convolutional function(s). Default is the density function of normal distribution with mean 0 and standard deviation 0.1.
#' @param grid the number of grid points used to construct the functional time series. Default is 1000.
#' @param sigma the standard deviation of O-U process. Default is 1.
#' @param ini the burn-in period.
#' @references
#' Liu, X., Xiao, H., and Chen, R. (2016) Convolutional autoregressive models for functional time series. \emph{Journal of Econometrics}, 194, 263-282.
#' @return The function returns a list with components:
#' \item{cfar}{a tmax-by-(grid+1) matrix following a CFAR(p) process.}
#' \item{epsilon}{the innovation at time tmax.}
#' @export
g_cfar <- function(tmax=1001,rho=5,phi_list=NULL, grid=1000,sigma=1,ini=100){
  #################################
  ###Simulate CFAR processes

  #########################################################################################
  #### OUTPUTS
  #######################
  #### A matrix f_grid is generated, with (tmax) rows and (grid+1) columns.
  ###########################################################################################

  #############################
  ###parameter setting
  ## rho=5	### parameter for error process epsilon_t, O-U process
  ##tmax= 1001;	### length of time
  ##iter= 100;	### number of replications in the simulation
  ##grid=1000;	### the number of grid points used to construct the functional time series X_t and epsilon_t
  ## phi_func User specified convolution functions. The following defaults are used if not given
  if(is.null(phi_list)){
    phi_func1= function(x){
      return(dnorm(x,mean=0,sd=0.1))
    }
    phi_list=phi_func1
  }
  if(is.null(grid)){
    grid=1000
  }
  if(is.null(sigma)){
    sigma=1
  }
  if(is.null(ini)){
    ini=100
  }
  #######################################################################################################
  x_grid= seq(0, 1, by= 1/grid);	### the grid points for input variables of functional time series in [0,1]
  b_grid= seq(-1, 1, by=1/grid);	### the grid points for input variables of b-spline functions in [-1,1]
  eps_grid= matrix(0, 1, grid+1);	### the grid points for input variables of error process in [0,1]
  fi = exp(-rho/grid);		### the correlation of two adjacent points for an O-U process
  p=length(phi_list)
  if(is.list(phi_list)){
    phi_b=sapply(phi_list,mapply,b_grid)
  }else{
    phi_b=matrix(phi_list(b_grid),2*grid+1,1)
  }
  ### the convolutional function values phi_1(s-u) and phi_2(s-u), for different s in [0,1]
  phi_matrix= matrix(0, (grid+1), p*(grid+1));
  for(j in 1:p){
    for (i in 1:(grid+1)){
      phi_matrix[i,((j-1)*(grid+1)+1):(j*(grid+1))]= phi_b[seq((i+grid),i,by=-1),j]/(grid+1)
    }
  }

  ##############################
  ### Generate data
  ### functional time series f_grid is X_t in the paper
  #### A matrix f_grid is generated, with (tmax*iter) rows and (grid+1) columns.
  ##  It contains iter CFAR(2) processes. The i-th process is in rows (i-1)*tmax+1:tmax.

  f_grid=matrix(0,tmax+ini,grid+1)
  i=1
  eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt(1-fi^2)*rnorm(n,0,1))
  f_grid[1,]= eps_grid;
  if(p>1){
    for(i in 2:p){
      eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt(1-fi^2)*rnorm(n,0,1))
      f_grid[i,]= phi_matrix[,1:((i-1)*(grid+1))]%*%matrix(t(f_grid[(i-1):1,]),(i-1)*(grid+1),1)+ eps_grid;
    }
  }
  for (i in (p+1):(tmax+ini)){
    eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt(1-fi^2)*rnorm(n,0,1))
    f_grid[i,]= phi_matrix%*%matrix(t(f_grid[(i-p):(i-1),]),p*(grid+1),1)+ eps_grid;
  }
  ### The last iteration of epsilon_t is returned.
  g_cfar <- list(cfar=f_grid[(ini+1):(tmax+ini),]*sigma,epsilon=eps_grid*sigma)
}



#' Generate a CFAR(2) Process with Heteroscedasticity and Irregular Observation Locations
#'
#' Generate a convolutional functional autoregressive process of order 2 with heteroscedasticity, irregular observation locations.
#' @param tmax length of time.
#' @param grid the number of grid points used to construct the functional time series.
#' @param rho parameter for O-U process (noise process).
#' @param min_obs the minimum number of observations at each time.
#' @param pois the mean for Poisson distribution. The number of observations at each follows a Poisson distribution plus min_obs.
#' @param phi_func1 the first convolutional function. Default is 0.5*x^2+0.5*x+0.13.
#' @param phi_func2 the second convolutional function. Default is 0.7*x^4-0.1*x^3-0.15*x.
#' @param weight the weight function to determine the standard deviation of O-U process (noise process). Default is 1.
#' @param ini the burn-in period.
#' @references
#' Liu, X., Xiao, H., and Chen, R. (2016) Convolutional autoregressive models for functional time series. \emph{Journal of Econometrics}, 194, 263-282.
#' @return The function returns a list with components:
#' \item{cfar2}{a tmax-by-(grid+1) matrix following a CFAR(1) process.}
#' \item{epsilon}{the innovation at time tmax.}
#' @examples
#' phi_func1= function(x){
#'	return(0.5*x^2+0.5*x+0.13)
#' }
#' phi_func2= function(x){
#'	return(0.7*x^4-0.1*x^3-0.15*x)
#' }
#' y=g_cfar2h(200,1000,1,40,5,phi_func1=phi_func1,phi_func2=phi_func2)
#' @export
g_cfar2h <- function(tmax=1001,grid=1000,rho=1,min_obs=40, pois=5,phi_func1=NULL, phi_func2=NULL, weight=NULL, ini=100){
  #################################
  ###Simulate CFAR(2) processes with heteroscedasticity, irregular observation locations
  ###################################################################################

  ######################################
  ### We need to have f_grid to record the functional time series
  ### We need to have x_pos_full to record the observation positions
  ######################################

  ###parameter setting
  #rho=1	### parameter for error process epsilon_t, O-U process
  #tmax= 1001;	### length of time
  #iter= 100;	### number of replications in the simulation
  #grid=1000;	### the number of grid points used to construct the functional time series X_t and epsilon_t\
  #min_obs=40	### the minimal number of observations at each time
  # phi_func1, phi_func2, weight: User specified convolution functions and weight function. The following defaults
  #  are used if not given
  if(is.null(phi_func1)){
    phi_func1= function(x){
      return(0.5*x^2+0.5*x+0.13)
    }
  }
  if(is.null(phi_func2)){
    phi_func2= function(x){
      return(0.7*x^4-0.1*x^3-0.15*x)
    }
  }

  if(is.null(weight)){
    ###heteroscedasticity weight function
    x_grid <- seq(0,1,by=1/grid)
    weight0= function(w){
      return((w<=0.6)*exp(-10*w)+(w>0.6)*(exp(-6)+0.2*(w-0.6))+0.1);
    }
    const=sum(weight0(x_grid)/(grid+1))
    weight= function(w){
      return(((w<=0.6)*exp(-10*w)+(w>0.6)*(exp(-6)+0.2*(w-0.6))+0.1)/const)
    }
  }
  if(is.null(ini)){
    ini=100
  }
  num_full=rpois(tmax,pois)+min_obs	###number of observations at time t follows a Poisson distribution
  ###plus a number, minimal observations is required
  #########################################################################
  x_grid= seq(0, 1, by= 1/grid);	### the grid points for input variables of functional time series in [0,1]
  b_grid= seq(-1, 1, by=1/grid);	### the grid points for input variables of b-spline functions in [-1,1]
  eps_grid= matrix(0, 1, grid+1);	### the grid points for input variables of error process in [0,1]
  fi = exp(-rho/grid);		### the correlation of two adjacent points for an O-U process
  phi_b1= phi_func1(b_grid)	### phi_1(s) when s in [-1,1]
  phi_b2= phi_func2(b_grid)	### phi_2(s) when s in [-1,1]

  ### the convolutional function values phi_1(s-u) and phi_2(s-u), for different s in [0,1]
  phi_matrix1= matrix(0, (grid+1), (grid+1));
  for (i in 1:(grid+1)){
    phi_matrix1[i,]= phi_b1[seq((i+grid),i,by=-1)]/(grid+1)
  }
  phi_matrix2= matrix(0, (grid+1), (grid+1));
  for (i in 1:(grid+1)){
    phi_matrix2[i,]= phi_b2[seq((i+grid),i,by=-1)]/(grid+1)
  }

  n=max(num_full)
  x_pos_full=matrix(0,tmax,n)	###observation positions
  for(k in 1:(tmax)){
    x_pos_full[k,1:num_full[k]]=sort(runif(num_full[k]))
  }
  f_grid=matrix(0,tmax+ini,grid+1)
  i=1
  eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt(1-fi^2)*rnorm(n,0,1))
  f_grid[1,]= eps_grid*weight(x_grid);
  i=2
  f_grid[i,]= phi_matrix1 %*% as.matrix(f_grid[i-1,])+ eps_grid*weight(x_grid);
  for (i in 3:(tmax+ini)){
    eps_grid= arima.sim(n= grid+1, model=list(ar=fi),rand.gen= function(n)sqrt(1-fi^2)*rnorm(n,0,1))
    f_grid[i,]= phi_matrix1%*%as.matrix(f_grid[i-1,]) + phi_matrix2%*%as.matrix(f_grid[i-2,])+ eps_grid*weight(x_grid);
  }

  ### Functional time series, number of observations at different time of periods,
  ### observation points at different time of periods, the last iteration of epsilon_t is returned.
  g_cfar2h <- list(cfar2h=f_grid[(ini+1):(tmax+ini),], num_obs=num_full, x_pos=x_pos_full,epsilon=eps_grid)
}

#' Estimation of a CFAR Process
#'
#' Estimation of a CFAR process.
#' @param f the functional time series.
#' @param p CFAR order.
#' @param df_b the degrees of freedom for natural cubic splines. Default is 10.
#' @param grid the number of gird points used to construct the functional time series and noise process. Default is 1000.
#' @references
#' Liu, X., Xiao, H., and Chen, R. (2016) Convolutional autoregressive models for functional time series. \emph{Journal of Econometrics}, 194, 263-282.
#' @return The function returns a list with components:
#' \item{phi_coef}{estimated spline coefficients for convolutional function values, a (2*grid+1)-by-p matrix.}
#' \item{phi_func}{estimated convolutional function(s), a (df_b+1)-by-p matrix.}
#' \item{rho}{estimated rho for O-U process (noise process).}
#' \item{sigma}{estimated sigma for O-U process (noise process).}
#' @import splines
#' @export
est_cfar <- function(f,p=3,df_b=10,grid=1000){
  if(!is.matrix(f))f <- as.matrix(f)
  t <- nrow(f)
  n=dim(f)[2]-1

  if(is.null(grid)){
    grid=1000
  }

  if(is.null(df_b)){
    df_b=10
  }

  ######################################################################
  ###INPUTS
  ####################
  ###parameter setting
  #rho=5	### parameter for error process epsilon_t, O-U process
  #t= 1000;	### length of time
  #iter= 100;	### number of replications in the simulation
  #grid=1000;	### the number of grid points used to construct the functional time series X_t and epsilon_t
  #df_b=10		### number of degrees for splines
  #n=100		### number of observations for each X_t, in the paper it is 'N'

  ############################################################################
  x_grid= seq(0, 1, by= 1/grid);	### the grid points for input variables of functional time series in [0,1]
  b_grid= seq(-1, 1, by=1/grid);	### the grid points for input variables of b-spline functions in [-1,1]
  eps_grid= matrix(0, 1, grid+1);	### the grid points for input variables of error process in [0,1]
  #############################
  ### It shows the best approximation of phi using b-spline with df=df_b

  coef_b=df_b+1;
  b_grid_sp= ns(b_grid, df=df_b);

  ###############################
  ###Estimation
  x= seq(0, 1, by=1/n);
  index= 1:(n+1)*(grid/n)-(grid/n)+1
  x_grid_sp= bs(x_grid, df=n, degree=1);	###basis function values if we interpolate discrete observations of x_t
  x_sp= x_grid_sp[index,];			###basis function values at each observation point
  df=n;
  coef= df+1

  ###convolutions of b-spline functions (df=df) for phi_1 and phi_2 and basis functions (df=n) for x
  x_grid_full= cbind(rep(1,grid+1), x_grid_sp);
  b_grid_full= cbind(rep(1,2*grid+1), b_grid_sp);
  b_matrix=matrix(0,grid+1,grid+1)
  bstar_grid2= matrix(0, grid+1, coef*coef_b)
  for (j in 1:(coef_b)){
    for (i in 1:(grid+1)){
      b_matrix[i,]= b_grid_full[seq((i+grid),i, by=-1),j]/(grid+1);
    }
    bstar_grid2[, (j-1)*coef+(1:coef)]= b_matrix %*% x_grid_full
  }

  bstar2= bstar_grid2[index,]
  bstar3= matrix(0, (n+1)*coef_b, coef)
  for (i in 1:coef_b){
    bstar3[(i-1)*(n+1)+(1:(n+1)),]=bstar2[, (i-1)*coef+(1:coef)]
  }


  indep=matrix(0,(n+1)*(t-1),coef_b)	### design matrix for b-spline coefficient of phi_1
  dsg=cbind(rep(1,n+1),x_sp);
  dsg_mat= solve(t(dsg)%*%dsg)%*%t(dsg)
  ahat= f%*%t(dsg_mat)	### b-spline coefficient when we interpolate x_t
  for (i in 2:t){
    M_t= matrix(ahat[i-1,]%*% t(bstar3), nrow=n+1, ncol=df_b+1)
    indep[(i-2)*(n+1)+(1:(n+1)),]=M_t	### design matrix for b-spline coefficient of phi_1
  }

  indep.tmp=NULL
  for(i in 1:p){
    tmp=indep[((i-1)*(n+1)+1):((n+1)*(t-p-1+i)),]
    indep.tmp=cbind(tmp,indep.tmp)
  }
  indep.p=indep.tmp

  pdf4=function(para4)
  {
    psi= para4;	### correlation between two adjacent observation point exp(-rho/n)
    psi_matinv= matrix(0, n+1, n+1)	### correlation matrix of error process at observation points
    psi_matinv[1,1]=1
    psi_matinv[n+1,n+1]=1;
    psi_matinv[1,2]=-psi;
    psi_matinv[n+1,n]=-psi;
    for (i in 2:n){
      psi_matinv[i,i]= 1+psi^2;
      psi_matinv[i,i+1]= -psi;
      psi_matinv[i,i-1]= -psi;
    }
    mat1= matrix(0, coef_b*p, coef_b*p);
    mat2= matrix(0, coef_b*p, 1)
    for(i in (p+1):t){
      tmp.n=n+1
      M_t=indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
      tmp_mat= t(M_t) %*% psi_matinv;
      mat1= mat1+ tmp_mat %*% M_t;
      mat2= mat2+ tmp_mat %*% f[i,];
    }
    phihat= solve(mat1)%*%mat2
    ehat=0
    log.mat=0
    for (i in (p+1):t){
      tmp.n=n+1
      M_t= indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
      eps= f[i,1:tmp.n]- M_t%*% phihat
      epart= t(eps)%*%psi_matinv %*%eps
      ehat= epart+ehat
      log.mat=log.mat+ log(det(psi_matinv))
    }
    l=-(t-p)*n/2*log(ehat)+1/2*log.mat
    return(-l);
  }
  para4= exp(-5)
  result4=optim(para4,pdf4, lower=0.001, upper=0.9999, method='L-BFGS-B')
  psi=result4$par
  psi_matinv= matrix(0, n+1, n+1)
  psi_matinv[1,1]=1
  psi_matinv[n+1,n+1]=1;
  psi_matinv[1,2]=-psi;
  psi_matinv[n+1,n]=-psi;
  for (i in 2:n){
    psi_matinv[i,i]= 1+psi^2;
    psi_matinv[i,i+1]= -psi;
    psi_matinv[i,i-1]= -psi;
  }
  mat1= matrix(0, p*(df_b+1), p*(df_b+1));	### matrix under the first brackets in formula (15)
  mat2= matrix(0, p*(df_b+1), 1);	### matrix under the second brackets in formula (15)
  for (i in (p+1):t){
    M_t= indep.p[(i-p-1)*(n+1)+(1:(n+1)),]
    tmp_mat= t(M_t) %*% psi_matinv;
    mat1= mat1+ tmp_mat %*% M_t;
    mat2= mat2+ tmp_mat %*% f[i,];
  }
  phihat= solve(mat1)%*%mat2
  phihat_func=matrix(0,p,(1+2*grid))
  for(i in 1:p){
    phihat_func[i,]= cbind(rep(1,2*grid+1),b_grid_sp)%*%phihat[(1+(coef_b*(i-1))):(i*coef_b)];		###b-spline coefficient of estimated phi_1
  }
  ehat=0
  for (i in (p+1):t){
    tmp.n= n+1
    M_t= indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
    eps= f[i,1:tmp.n]- M_t%*% phihat
    epart= t(eps)%*%psi_matinv %*%eps
    ehat= epart+ehat
  }
  rho_hat= -log(psi)*n
  sigma_hat=sqrt(ehat/((t-p)*(n+1))/(1-psi^2))
  phi_coef=t(matrix(phihat,coef_b,p))
  est_cfar <- list(phi_coef=phi_coef,phi_func=phihat_func,rho=rho_hat,sigma=sigma_hat)
  return(est_cfar)
}


#' F Test for a CFAR Process
#'
#' F test for a CFAR process to specify CFAR order.
#' @param f the functional time series.
#' @param p.max maximum CFAR order. Default is 6.
#' @param df_b the degrees of freedom for natural cubic splines. Default is 10.
#' @param grid the number of gird points used to construct the functional time series and noise process. Default is 1000.
#' @references
#' Liu, X., Xiao, H., and Chen, R. (2016) Convolutional autoregressive models for functional time series. \emph{Journal of Econometrics}, 194, 263-282.
#' @return The function outputs F test statistics and their p-values.
#' @export
F_test_cfar <- function(f,p.max=6,df_b=10,grid=1000){
  #	library(MASS)
  #	library(splines)
  ##########################
  ### F test: CFAR(0)vs CFAR(1), CFAR(1)vs CFAR(2), CFAR(2)vs CFAR(3) when rho is known.
  #
  if(!is.matrix(f))f <- as.matrix(f)
  t <- nrow(f)
  n=dim(f)[2]-1

  if(is.null(p.max)){
    p.max=6
  }
  if(is.null(df_b)){
    df_b=10
  }
  if(is.null(grid)){
    grid=1000
  }
  ##################################################################################
  ###parameter setting
  #rho=5	### parameter for error process epsilon_t, O-U process
  #t= 1000;	### length of time
  ##iter= 100;	### number of replications in the simulation
  #grid=1000;	### the number of grid points used to construct the functional time series X_t and epsilon_t
  #df_b=10		### number of degrees for splines
  #n=100		### number of observations for each X_t, in the paper it is 'N'
  ###############################################################################
  x_grid= seq(0, 1, by= 1/grid);	### the grid points for input variables of functional time series in [0,1]
  b_grid= seq(-1, 1, by=1/grid);	### the grid points for input variables of b-spline functions in [-1,1]
  eps_grid= matrix(0, 1, grid+1);	### the grid points for input variables of error process in [0,1]
  coef_b=df_b+1;				### number of df for b-spline
  b_grid_sp= ns(b_grid, df=df_b);	###b-spline functions
  x= seq(0, 1, by=1/n);			###observations points for X_t process
  index= 1:(n+1)*(grid/n)-(grid/n)+1
  x_grid_sp= bs(x_grid, df=n, degree=1);	###basis function values if we interpolate discrete observations of x_t
  x_sp= x_grid_sp[index,];			###basis function values at each observation point
  df=n;			### number of interpolation of X_t
  coef= df+1;
  ###convolutions of b-spline functions (df=df_b) for phi and basis functions (df=n) for x
  x_grid_full= cbind(rep(1,grid+1), x_grid_sp);
  b_grid_full= cbind(rep(1,2*grid+1), b_grid_sp);
  b_matrix=matrix(0,grid+1,grid+1)
  bstar_grid2= matrix(0, grid+1, coef*coef_b)
  for (j in 1:(coef_b)){
    for (i in 1:(grid+1)){
      b_matrix[i,]= b_grid_full[seq((i+grid),i, by=-1),j]/(grid+1);
    }
    bstar_grid2[, (j-1)*coef+(1:coef)]= b_matrix %*% x_grid_full
    ###	print(j)
  }
  bstar2= bstar_grid2[index,]
  bstar3= matrix(0, (n+1)*coef_b, coef)
  for (i in 1:coef_b){
    bstar3[(i-1)*(n+1)+(1:(n+1)),]=bstar2[, (i-1)*coef+(1:coef)]
  }
  dsg=cbind(rep(1,n+1),x_sp);
  dsg_mat= solve(t(dsg)%*%dsg)%*%t(dsg)
  ahat= f%*%t(dsg_mat)
  indep=matrix(0,(n+1)*(t-1),coef_b)
  for (i in 2:t){
    M_t= matrix(ahat[i-1,]%*% t(bstar3), nrow=n+1, ncol=df_b+1)
    indep[(i-2)*(n+1)+(1:(n+1)),]=M_t
  }
  #########################################
  ######### AR(0) and AR(1)
  #########################################
  #######################
  ####cat("Data set: ",q,"\n")
  ### Fit cfar(1) model
  ###	f= f_grid[(q-1)*tmax+(1:t)+(tmax-t), index];
  p=1
  dsg=cbind(rep(1,n+1),x_sp);
  dsg_mat= solve(t(dsg)%*%dsg)%*%t(dsg)
  ahat= f%*%t(dsg_mat)
  indep=matrix(0,(n+1)*(t-1),coef_b)
  for (i in 2:t){
    M_t= matrix(ahat[i-1,]%*% t(bstar3), nrow=n+1, ncol=df_b+1)
    indep[(i-2)*(n+1)+(1:(n+1)),]=M_t
  }
  est=est_cfar(f,1,df_b=df_b,grid)
  phihat= est$phi_coef
  psi=exp(-est$rho/n)
  mat1= matrix(0, coef_b*p, coef_b*p);
  mat2= matrix(0, coef_b*p, 1);
  phi_matinv=matrix(0, n+1, n+1)
  phi_matinv[1,1]=1;
  phi_matinv[n+1,n+1]=1;
  phi_matinv[1,2]= -psi;
  phi_matinv[n+1,n]= -psi
  for (i in 2:n){
    phi_matinv[i,i]= 1+psi^2;
    phi_matinv[i,i+1]= -psi;
    phi_matinv[i,i-1]= -psi
  }

  predict= t(matrix(indep%*%matrix(phihat,df_b+1,1),ncol=t-1,nrow=n+1))
  ssr1= sum(diag((predict)%*%phi_matinv%*% t(predict)))
  sse1= sum(diag((f[2:t,]-predict[1:(t-1),])%*%phi_matinv%*% t(f[2:t,]-predict[1:(t-1),])))
  sse0=sum(diag(f[2:t,]%*%phi_matinv%*% t(f[2:t,])))
  statistic= (sse0-sse1)/coef_b/sse1*((n+1)*(t-2)-coef_b)
  pval=1-pf(statistic,coef_b,(t-2)*(n+1)-coef_b)
  cat("Test and p-value of Order 0 vs Order 1: ","\n")
  print(c(statistic,pval))

  p=1
  sse.pre= sum(diag((f[(p+2):t,]-predict[2:(t-p),])%*%phi_matinv%*% t(f[(p+2):t,]-predict[2:(t-p),])))

  for(p in 2:p.max){
    indep.tmp=NULL
    for(i in 1:p){
      tmp=indep[((i-1)*(n+1)+1):((n+1)*(t-p-1+i)),]
      indep.tmp=cbind(tmp,indep.tmp)
    }
    indep.p=indep.tmp
    pdf4=function(para4)
    {
      mat1= matrix(0, coef_b*p, coef_b*p);
      mat2= matrix(0, coef_b*p, 1);
      psi= para4;
      phi_matinv=matrix(0, n+1, n+1)
      phi_matinv[1,1]=1;
      phi_matinv[n+1,n+1]=1;
      phi_matinv[1,2]= -psi;
      phi_matinv[n+1,n]= -psi
      for (i in 2:n){
        phi_matinv[i,i]= 1+psi^2;
        phi_matinv[i,i+1]= -psi;
        phi_matinv[i,i-1]= -psi
      }
      psi_matinv=phi_matinv
      for(i in (p+1):t){
        tmp.n=n+1
        tmp.index=which(!is.na(f[i,]))
        M_t=indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
        tmp_mat= t(M_t) %*% psi_matinv;
        mat1= mat1+ tmp_mat %*% M_t;
        mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];
      }
      phihat= solve(mat1)%*%mat2
      ehat=0
      log.mat=0
      for (i in (p+1):t){
        tmp.n=n+1
        M_t= indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
        eps= f[i,1:tmp.n]- M_t%*% phihat
        epart= t(eps)%*%psi_matinv %*%eps
        ehat= epart+ehat
        log.mat=log.mat+ log(det(psi_matinv))
      }
      l=-((t-p)*(n+1))/2*log(ehat)+1/2*log.mat
      return(-l);
    }
    para4= exp(-5)
    result4=optim(para4,pdf4, lower=0.001, upper=0.999, method='L-BFGS-B')

    mat1= matrix(0, p*coef_b, p*coef_b);
    mat2= matrix(0, coef_b*p, 1);
    psi_invall=matrix(0,(n+1)*(t-1),(n+1))
    psi=result4$par
    phi_matinv=matrix(0, n+1, n+1)
    phi_matinv[1,1]=1;
    phi_matinv[n+1,n+1]=1;
    phi_matinv[1,2]= -psi;
    phi_matinv[n+1,n]= -psi
    psi_matinv=phi_matinv
    for (i in 2:n){
      phi_matinv[i,i]= 1+psi^2;
      phi_matinv[i,i+1]= -psi;
      phi_matinv[i,i-1]= -psi
    }
    psi_matinv=phi_matinv
    for(i in (p+1):t){
      tmp.n=n+1
      M_t= indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
      tmp_mat= t(M_t) %*% psi_matinv;
      mat1= mat1+ tmp_mat %*% M_t;
      mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];
    }
    phihat= solve(mat1)%*%mat2
    predict= t(matrix(indep.p%*%phihat,ncol=t-p,nrow=n+1))
    ssr.p= sum(diag((predict)%*%phi_matinv%*% t(predict)))
    sse.p= sum(diag((f[(p+1):t,]-predict)%*%phi_matinv%*% t(f[(p+1):t,]-predict)))
    statistic= (sse.pre-sse.p)/coef_b/sse.p*((n+1)*(t-p-1)-p*coef_b)
    pval=1-pf(statistic,coef_b,(t-p)*(n+1)-p*coef_b)
    cat("Test and  p-value of Order", p-1, "vs Order",p,": ","\n")
    print(c(statistic,pval))
    sse.pre= sum(diag((f[(p+2):t,]-predict[2:(t-p),])%*%phi_matinv%*% t(f[(p+2):t,]-predict[2:(t-p),])))
  }
}



#' F Test for a CFAR Process with Heteroscedasticity and Irregular Observation Locations
#'
#' F test for a CFAR process with heteroscedasticity and irregular observation locations to specify the CFAR order.
#' @param f the functional time series.
#' @param weight the covariance functions for noise process.
#' @param p.max the maximum CFAR order. Default is 3.
#' @param grid the number of gird points used to construct the functional time series and noise process. Default is 1000.
#' @param df_b the degrees of freedom for natural cubic splines. Default is 10.
#' @param num_obs the numbers of observations. It is a t-by-1 vector, where t is the length of time.
#' @param x_pos the observation location matrix. If the locations are regular, it is a t-by-(n+1) matrix with all entries 1/n.
#' @references
#' Liu, X., Xiao, H., and Chen, R. (2016) Convolutional autoregressive models for functional time series. \emph{Journal of Econometrics}, 194, 263-282.
#' @return The function outputs F test statistics and their p-values.
#' @export
F_test_cfarh <- function(f,weight,p.max=3,grid=1000,df_b=10,num_obs=NULL,x_pos=NULL){
  #	library(MASS)
  #	library(splines)
  ########################

  if(!is.matrix(f))f <- as.matrix(f)
  t <- nrow(f)

  if(is.null(p.max)){
    p.max=3
  }

  if(is.null(num_obs)){
    num_obs=dim(f)[2]
    n=dim(f)[2]
  }else{
    n=max(num_obs)
  }

  if(length(num_obs)!=t){
    num_obs=rep(n,t)
  }
  if(is.null(df_b)){
    df_b=10
  }
  if(is.null(grid)){
    grid=1000
  }
  if(is.null(x_pos)){
    x_pos=matrix(rep(seq(0,1,by=1/(num_obs[1]-1)),each=t),t,num_obs[1])
  }

  ##################
  ###parameter setting
  #grid=1000;	### the number of grid points used to construct the functional time series X_t and epsilon_t\
  #coef_b=6	### k=6
  # num_obs is a t by 1 vector which records N_t for each time

  #######################################################K
  x_grid=seq(0,1,by=1/grid)
  coef_b=df_b+1
  b_grid=seq(-1,1,by=1/grid)
  b_grid_sp = ns(b_grid,df=df_b);
  b_grid_full= cbind(rep(1,2*grid+1), b_grid_sp);
  indep= matrix(0, (n+1)*(t-1),coef_b)
  bstar_grid= matrix(0, grid+1, coef_b)
  b_matrix= matrix(0, (grid+1), (grid+1))
  bstar_grid_full= matrix(0, (grid+1)*t,coef_b)
  index=matrix(0,t,n)
  for(i in 1:t){
    for(j in 1:num_obs[i]){
      index[i,j]=which(abs(x_pos[i,j]-x_grid)==min(abs(x_pos[i,j]-x_grid)))
    }
  }
  for (i in 1:(t-1)){
    rec_x= approx(x_pos[i,1:num_obs[i]],f[i,1:num_obs[i]],xout=x_grid,rule=2,method='linear')
    ### rec_x is the interpolation of x_t
    for(j in 1:coef_b){
      for (k in 1:(grid+1)){
        b_matrix[k,]= b_grid_full[seq((k+grid),k, by=-1),j]/(grid+1);
      }
      bstar_grid[, j]= b_matrix %*% matrix(rec_x$y,ncol=1,nrow=grid+1)
    }
    bstar_grid_full[(i-1)*(grid+1)+1:(grid+1),]=bstar_grid	###convolution of basis spline function and x_t
    tmp=bstar_grid[index[i+1,1:num_obs[i+1]],]
    indep[(i-1)*(n+1)+1:num_obs[i+1],]=tmp
  }




  #################AR(1)
  pdf4=function(para4)	###Q function in formula (12)
  {
    mat1= matrix(0, df_b+1, df_b+1);
    mat2= matrix(0, df_b+1, 1);
    ### correlation matrix of error process at observation points
    psi_invall=matrix(0,(n+1)*(t-1),(n+1))
    for(i in 2:t){
      psi= para4;
      tmp.n=num_obs[i]
      psi_mat= matrix(1, tmp.n, tmp.n)
      for(k in 1:(tmp.n)){
        for(j in 1:(tmp.n)){
          psi_mat[k,j]= psi^(abs(x_pos[i,j]-x_pos[i,k]))
        }
      }
      psi_matinv2=solve(diag(weight(x_pos[i,1:tmp.n]))%*%psi_mat%*%diag(weight(x_pos[i,1:tmp.n])))
      psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]=psi_matinv2
      M_t=indep[(i-2)*(n+1)+(1:tmp.n),]
      tmp_mat= t(M_t) %*% psi_matinv2;
      mat1= mat1+ tmp_mat %*% M_t;			### matrix under the first brackets in formula (15)
      mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];	### matrix under the second brackets in formula (15)
    }
    phihat= solve(mat1)%*%mat2
    ehat=0	###e(beta,phi) in formula (12)
    log.mat=0
    for (i in 2:t){
      tmp.n=num_obs[i]
      M_t= indep[(i-2)*(n+1)+(1:tmp.n),]
      eps= f[i,1:tmp.n]- M_t%*% phihat
      psi_matinv= psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
      epart= t(eps)%*%psi_matinv %*%eps
      ehat= epart+ehat
      log.mat=log.mat+ log(det(psi_matinv))
    }
    l=-sum(num_obs[2:t])/2*log(ehat)+1/2*log.mat	### the number defined in formula (14)
    return(-l);
  }
  para4= exp(-1)
  result4=optim(para4,pdf4, lower=0.001, upper=0.999, method='L-BFGS-B')
  mat1= matrix(0, df_b+1, df_b+1);
  mat2= matrix(0, df_b+1, 1);
  psi_invall=matrix(0,(n+1)*(t-1),(n+1))
  psi=result4$par
  for(i in 2:t){
    tmp.n=num_obs[i]
    psi_mat= matrix(1, tmp.n, tmp.n)
    for(k in 1:tmp.n){
      for(j in 1:tmp.n){
        psi_mat[k,j]=psi^(abs(x_pos[i,j]-x_pos[i,k]))
      }
    }
    psi_matinv2=solve(diag(weight(x_pos[i,1:tmp.n]))%*%psi_mat%*%diag(weight(x_pos[i,1:tmp.n])))
    psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]= psi_matinv2
    M_t= indep[(i-2)*(n+1)+(1:tmp.n),]
    tmp_mat= t(M_t) %*% psi_matinv2;
    mat1= mat1+ tmp_mat %*% M_t;
    mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];
  }
  phihat= solve(mat1)%*%mat2
  ehat=0		###e(beta,phi) in formula (12)
  for (i in 2:t){
    tmp.n=num_obs[i]
    M_t= indep[(i-2)*(n+1)+(1:tmp.n),]
    eps= f[i,1:tmp.n]- M_t%*% phihat
    psi_matinv=psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
    epart= t(eps)%*%psi_matinv %*%eps
    ehat= epart+ehat
  }
  rho_hat= -log(psi)
  sigma_hat=ehat/sum(num_obs[2:t])*2*rho_hat

  sigma_hat1=sigma_hat
  rho_hat1=rho_hat
  predict1= t(matrix(indep%*%phihat,ncol=t-1,nrow=n+1))
  ssr1=0
  sse1=0
  for (i in 2:t){
    tmp.n=num_obs[i]
    psi_matinv= psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
    ssr1= ssr1+ sum((predict1[i-1,1:tmp.n] %*% psi_matinv) * (predict1[i-1,1:tmp.n]))
    sse1= sse1+ sum(((f[i,1:tmp.n]- predict1[i-1,1:tmp.n]) %*% psi_matinv) *(f[i,1:tmp.n]-predict1[i-1,1:tmp.n]))
  }

  ###############################
  ######   No AR term
  psi=exp(-rho_hat1)
  mat1= matrix(0, df_b+1, df_b+1);
  mat2= matrix(0, df_b+1, 1);
  phihat=matrix(0,coef_b,1)
  ehat=0
  for (i in 2:t){
    tmp.n=num_obs[i]
    M_t= indep[(i-2)*(n+1)+(1:tmp.n),]
    eps= f[i,1:tmp.n]- M_t%*% phihat
    psi_matinv=psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
    epart= t(eps)%*%psi_matinv %*%eps
    ehat= epart+ehat
  }
  rho_hat= -log(psi)
  sigma_hat=ehat/sum(num_obs[2:t])*2*rho_hat
  rho_hat0=rho_hat
  sigma_hat0=sigma_hat
  sse0=0
  test=0
  for (i in 2:t){
    tmp.n=num_obs[i]
    psi_matinv= psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
    sse0= sse0+ sum((f[i,1:tmp.n] %*% psi_matinv) *(f[i,1:tmp.n]))
  }

  statistic= (sse0-sse1)/coef_b/sse1*(sum(num_obs[2:t])-coef_b)
  pval=1-pf(statistic,coef_b,sum(num_obs[2:t])-coef_b)
  cat("Test and  p-value of Order 0 vs Order 1: ","\n")
  print(c(statistic,pval))


  indep.p=indep
  if (p.max>1){
    for(p in 2:p.max){
      indep.pre=indep.p[(n+2):((n+1)*(t-p+1)),]
      indep.tmp=indep
      for (i in 1:(t-p)){
        for(j in 1:num_obs[i+p]){
          index[i+p,j]= which(abs(x_pos[i+p,j]-x_grid)==min(abs(x_pos[i+p,j]-x_grid)))
        }
        tmp=bstar_grid_full[(i-1)*(grid+1)+index[i+p,1:num_obs[i+p]],]
        indep.tmp[(i-1)*(n+1)+1:num_obs[i+p],]=tmp
      }
      indep.test=cbind(indep.p[(n+2):((n+1)*(t-p+1)),], indep.tmp[1:((n+1)*(t-p)),])
      indep.p=indep.test

      pdf4=function(para4)	###Q function in formula (12)
      {
        mat1= matrix(0, p*(df_b+1), p*(df_b+1));
        mat2= matrix(0, p*(df_b+1), 1);
        ### correlation matrix of error process at observation points
        psi_invall=matrix(0,(n+1)*(t-1),(n+1))
        for(i in (p+1):t){
          psi= para4;
          tmp.n=num_obs[i]
          psi_mat= matrix(1, tmp.n, tmp.n)
          for(k in 1:(tmp.n)){
            for(j in 1:(tmp.n)){
              psi_mat[k,j]= psi^(abs(x_pos[i,j]-x_pos[i,k]))
            }
          }
          psi_matinv2=solve(diag(weight(x_pos[i,1:tmp.n]))%*%psi_mat%*%diag(weight(x_pos[i,1:tmp.n])))
          psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]=psi_matinv2
          M_t=indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
          tmp_mat= t(M_t) %*% psi_matinv2;
          mat1= mat1+ tmp_mat %*% M_t;			### matrix under the first brackets in formula (15)
          mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];	### matrix under the second brackets in formula (15)
        }
        phihat= solve(mat1)%*%mat2
        ehat=0	###e(beta,phi) in formula (12)
        log.mat=0
        for (i in (p+1):t){
          tmp.n=num_obs[i]
          M_t= indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
          eps= f[i,1:tmp.n]- M_t%*% phihat
          psi_matinv= psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
          epart= t(eps)%*%psi_matinv %*%eps
          ehat= epart+ehat
          log.mat=log.mat+ log(det(psi_matinv))
        }
        l=-sum(num_obs[(p+1):t])/2*log(ehat)+1/2*log.mat	### the number defined in formula (14)
        return(-l);
      }
      para4= exp(-1)
      result4=optim(para4,pdf4, lower=0.001, upper=0.999, method='L-BFGS-B')
      mat1= matrix(0, p*(df_b+1), p*(df_b+1));
      mat2= matrix(0, p*(df_b+1), 1);
      psi_invall=matrix(0,(n+1)*(t-1),(n+1))
      psi=result4$par
      for(i in (p+1):t){
        tmp.n=num_obs[i]
        psi_mat= matrix(1, tmp.n, tmp.n)
        for(k in 1:tmp.n){
          for(j in 1:tmp.n){
            psi_mat[k,j]=psi^(abs(x_pos[i,j]-x_pos[i,k]))
          }
        }
        psi_matinv2=solve(diag(weight(x_pos[i,1:tmp.n]))%*%psi_mat%*%diag(weight(x_pos[i,1:tmp.n])))
        psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]= psi_matinv2
        M_t= indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
        tmp_mat= t(M_t) %*% psi_matinv2;
        mat1= mat1+ tmp_mat %*% M_t;
        mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];
      }
      phihat= solve(mat1)%*%mat2
      ehat=0		###e(beta,phi) in formula (12)
      for (i in (p+1):t){
        tmp.n=num_obs[i]
        M_t= indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
        eps= f[i,1:tmp.n]- M_t%*% phihat
        psi_matinv=psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
        epart= t(eps)%*%psi_matinv %*%eps
        ehat= epart+ehat
      }
      rho_hat= -log(psi)
      sigma_hat=ehat/sum(num_obs[(p+1):t])*2*rho_hat

      sigma_hat1=sigma_hat
      rho_hat1=rho_hat
      predict= t(matrix(indep.p%*%phihat,ncol=t-p,nrow=n+1))
      ssr.p= 0
      sse.p=0
      for (i in (p+1):t){
        tmp.n=num_obs[i]
        psi_matinv= psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
        ssr.p= ssr.p+ sum((predict[i-p,1:tmp.n] %*% psi_matinv) * (predict[i-p,1:tmp.n]))
        sse.p= sse.p+ sum(((f[i,1:tmp.n]- predict[i-p,1:tmp.n]) %*% psi_matinv) *(f[i,1:tmp.n]-predict[i-p,1:tmp.n]))
      }

      mat1= matrix(0, coef_b*(p-1), coef_b*(p-1));
      mat2= matrix(0, coef_b*(p-1), 1);
      psi_invall=matrix(0,(n+1)*(t-1),(n+1))
      for(i in (p+1):t){
        tmp.n=num_obs[i]
        tmp.index=which(!is.na(f[i,]))
        psi_mat= matrix(1, tmp.n, tmp.n)
        for(k in 1:tmp.n){
          for(j in 1:tmp.n){
            psi_mat[k,j]=psi^(abs(x_pos[i,j]-x_pos[i,k]))
          }
        }
        psi_matinv2=solve(diag(weight(x_pos[i,1:tmp.n]))%*%psi_mat%*%diag(weight(x_pos[i,1:tmp.n])))
        psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]= psi_matinv2
        M_t= indep.pre[(i-p-1)*(n+1)+(1:tmp.n),]
        tmp_mat= t(M_t) %*% psi_matinv2;
        mat1= mat1+ tmp_mat %*% M_t;
        mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];
      }
      phihat= solve(mat1)%*%mat2
      predict= t(matrix(indep.pre%*%phihat,ncol=t-p,nrow=n+1))
      ssr.pre=0
      sse.pre=0.0000
      for (i in (p+1):t){
        tmp.n=num_obs[i]
        tmp.index=which(!is.na(f[i,]))
        psi_matinv= psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
        ssr.pre= ssr.pre+ sum((predict[i-p,1:tmp.n] %*% psi_matinv) * (predict[i-p,1:tmp.n]))
        sse.pre= sse.pre+ sum(((f[i,1:tmp.n]- predict[i-p,1:tmp.n]) %*% psi_matinv) *(f[i,1:tmp.n]-predict[i-p,1:tmp.n]))
      }
      statistic= (sse.pre-sse.p)/coef_b/sse.p*(sum(num_obs[(p+1):t])-coef_b)
      pval=1-pf(statistic,coef_b,sum(num_obs[(p+1):t])-coef_b)
      cat("Test and  p-value of Order",p-1," vs Order", p,": ","\n")
      print(c(statistic,pval))
    }
  }
}




#' Estimation of a CFAR Process with Heteroscedasticity and Irregualar Observation Locations
#'
#' Estimation of a CFAR process with heteroscedasticity and irregualar observation locations.
#' @param f the functional time series.
#' @param weight the covariance functions of noise process.
#' @param p CFAR order.
#' @param grid the number of gird points used to construct the functional time series and noise process. Default is 1000.
#' @param df_b the degrees of freedom for natural cubic splines. Default is 10.
#' @param num_obs the numbers of observations. It is a t-by-1 vector, where t is the length of time.
#' @param x_pos the observation location matrix. If the locations are regular, it is a t-by-(n+1) matrix with all entries 1/n.
#' @references
#' Liu, X., Xiao, H., and Chen, R. (2016) Convolutional autoregressive models for functional time series. \emph{Journal of Econometrics}, 194, 263-282.
#' @return The function returns a list with components:
#' \item{phi_coef}{estimated spline coefficients for convolutional function(s).}
#' \item{phi_func}{estimated convolutional function(s).}
#' \item{rho}{estimated rho for O-U process (noise process).}
#' \item{sigma}{estimated sigma for O-U process (noise process).}
#' @export
est_cfarh <- function(f,weight,p=2,grid=1000,df_b=5, num_obs=NULL,x_pos=NULL){
  #	library(MASS)
  #	library(splines)
  ########################
  ### CFAR(2) with processes with heteroscedasticity, irregular observation locations
  ### Estimation of phi(), rho and sigma, and F test

  if(!is.matrix(f))f <- as.matrix(f)
  t <- nrow(f)

  if(is.null(num_obs)){
    num_obs=dim(f)[2]
    n=dim(f)[2]
  }else{
    n=max(num_obs)
  }

  if(length(num_obs)!=t){
    num_obs=rep(n,t)
  }

  if(is.null(x_pos)){
    x_pos=matrix(rep(seq(0,1,by=1/(num_obs[1]-1)),each=t),t,num_obs[1])
  }
  if(is.null(df_b)){
    df_b=10
  }
  if(is.null(grid)){
    grid=1000
  }
  ##################
  ###parameter setting
  #rho=1	### parameter for error process epsilon_t, O-U process
  #tmax= 1001;	### maximum length of time to generated data
  #t=1000	### length of time
  #iter= 100;	### number of replications in the simulation
  #grid=1000;	### the number of grid points used to construct the functional time series X_t and epsilon_t\
  #coef_b=6	### k=6
  #min_obs=40	### the minimal number of observations at each time
  #######################################################K
  x_grid=seq(0,1,by=1/grid)
  coef_b=df_b+1
  b_grid=seq(-1,1,by=1/grid)
  b_grid_sp = ns(b_grid,df=df_b);
  b_grid_full= cbind(rep(1,2*grid+1), b_grid_sp);
  indep= matrix(0, (n+1)*(t-1),coef_b)
  bstar_grid= matrix(0, grid+1, coef_b)
  b_matrix= matrix(0, (grid+1), (grid+1))
  bstar_grid_full= matrix(0, (grid+1)*t,coef_b)
  index=matrix(0,t,n)
  for(i in 1:t){
    for(j in 1:num_obs[i]){
      index[i,j]=which(abs(x_pos[i,j]-x_grid)==min(abs(x_pos[i,j]-x_grid)))
    }
  }
  for (i in 1:(t-1)){
    rec_x= approx(x_pos[i,1:num_obs[i]],f[i,1:num_obs[i]],xout=x_grid,rule=2,method='linear')
    ### rec_x is the interpolation of x_t
    for(j in 1:coef_b){
      for (k in 1:(grid+1)){
        b_matrix[k,]= b_grid_full[seq((k+grid),k, by=-1),j]/(grid+1);
      }
      bstar_grid[, j]= b_matrix %*% matrix(rec_x$y,ncol=1,nrow=grid+1)
    }
    bstar_grid_full[(i-1)*(grid+1)+1:(grid+1),]=bstar_grid	###convolution of basis spline function and x_t
    tmp=bstar_grid[index[i+1,1:num_obs[i+1]],]
    indep[(i-1)*(n+1)+1:num_obs[i+1],]=tmp
  }


  if(p==1){
    #################AR(1)
    pdf4=function(para4)	###Q function in formula (12)
    {
      mat1= matrix(0, df_b+1, df_b+1);
      mat2= matrix(0, df_b+1, 1);
      ### correlation matrix of error process at observation points
      psi_invall=matrix(0,(n+1)*(t-1),(n+1))
      for(i in 2:t){
        psi= para4;
        tmp.n=num_obs[i]
        psi_mat= matrix(1, tmp.n, tmp.n)
        for(k in 1:(tmp.n)){
          for(j in 1:(tmp.n)){
            psi_mat[k,j]= psi^(abs(x_pos[i,j]-x_pos[i,k]))
          }
        }
        psi_matinv2=solve(diag(weight(x_pos[i,1:tmp.n]))%*%psi_mat%*%diag(weight(x_pos[i,1:tmp.n])))
        psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]=psi_matinv2
        M_t=indep[(i-2)*(n+1)+(1:tmp.n),]
        tmp_mat= t(M_t) %*% psi_matinv2;
        mat1= mat1+ tmp_mat %*% M_t;			### matrix under the first brackets in formula (15)
        mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];	### matrix under the second brackets in formula (15)
      }
      phihat= solve(mat1)%*%mat2
      ehat=0	###e(beta,phi) in formula (12)
      log.mat=0
      for (i in 2:t){
        tmp.n=num_obs[i]
        M_t= indep[(i-2)*(n+1)+(1:tmp.n),]
        eps= f[i,1:tmp.n]- M_t%*% phihat
        psi_matinv= psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
        epart= t(eps)%*%psi_matinv %*%eps
        ehat= epart+ehat
        log.mat=log.mat+ log(det(psi_matinv))
      }
      l=-sum(num_obs[2:t])/2*log(ehat)+1/2*log.mat	### the number defined in formula (14)
      return(-l);
    }
    para4= exp(-1)
    result4=optim(para4,pdf4, lower=0.001, upper=0.999, method='L-BFGS-B')
    mat1= matrix(0, df_b+1, df_b+1);
    mat2= matrix(0, df_b+1, 1);
    psi_invall=matrix(0,(n+1)*(t-1),(n+1))
    psi=result4$par
    for(i in 2:t){
      tmp.n=num_obs[i]
      psi_mat= matrix(1, tmp.n, tmp.n)
      for(k in 1:tmp.n){
        for(j in 1:tmp.n){
          psi_mat[k,j]=psi^(abs(x_pos[i,j]-x_pos[i,k]))
        }
      }
      psi_matinv2=solve(diag(weight(x_pos[i,1:tmp.n]))%*%psi_mat%*%diag(weight(x_pos[i,1:tmp.n])))
      psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]= psi_matinv2
      M_t= indep[(i-2)*(n+1)+(1:tmp.n),]
      tmp_mat= t(M_t) %*% psi_matinv2;
      mat1= mat1+ tmp_mat %*% M_t;
      mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];
    }
    phihat= solve(mat1)%*%mat2
    phihat_func=matrix(0,p,(1+2*grid))
    for(i in 1:p){
      phihat_func[i,]= cbind(rep(1,2*grid+1),b_grid_sp)%*%phihat[(1+(coef_b*(i-1))):(i*coef_b)];		###b-spline coefficient of estimated phi_1
    }
    ehat=0		###e(beta,phi) in formula (12)
    for (i in 2:t){
      tmp.n=num_obs[i]
      M_t= indep[(i-2)*(n+1)+(1:tmp.n),]
      eps= f[i,1:tmp.n]- M_t%*% phihat
      psi_matinv=psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
      epart= t(eps)%*%psi_matinv %*%eps
      ehat= epart+ehat
    }
    rho_hat= -log(psi)
    sigma_hat=ehat/sum(num_obs[2:t])*2*rho_hat

  }
  indep.p=indep
  if(p>1){
    for(q in 2:p){
      indep.tmp=indep
      for (i in 1:(t-q)){
        for(j in 1:num_obs[i+q]){
          index[i+q,j]= which(abs(x_pos[i+q,j]-x_grid)==min(abs(x_pos[i+q,j]-x_grid)))
        }
        tmp=bstar_grid_full[(i-1)*(grid+1)+index[i+q,1:num_obs[i+q]],]
        indep.tmp[(i-1)*(n+1)+1:num_obs[i+q],]=tmp
      }
      indep.test=cbind(indep.p[(n+2):((n+1)*(t-q+1)),], indep.tmp[1:((n+1)*(t-q)),])
      indep.p=indep.test

    }

    pdf4=function(para4)	###Q function in formula (12)
    {
      mat1= matrix(0, p*(df_b+1), p*(df_b+1));
      mat2= matrix(0, p*(df_b+1), 1);
      ### correlation matrix of error process at observation points
      psi_invall=matrix(0,(n+1)*(t-1),(n+1))
      for(i in (p+1):t){
        psi= para4;
        tmp.n=num_obs[i]
        psi_mat= matrix(1, tmp.n, tmp.n)
        for(k in 1:(tmp.n)){
          for(j in 1:(tmp.n)){
            psi_mat[k,j]= psi^(abs(x_pos[i,j]-x_pos[i,k]))
          }
        }
        psi_matinv2=solve(diag(weight(x_pos[i,1:tmp.n]))%*%psi_mat%*%diag(weight(x_pos[i,1:tmp.n])))
        psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]=psi_matinv2
        M_t=indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
        tmp_mat= t(M_t) %*% psi_matinv2;
        mat1= mat1+ tmp_mat %*% M_t;			### matrix under the first brackets in formula (15)
        mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];	### matrix under the second brackets in formula (15)
      }
      phihat= solve(mat1)%*%mat2
      ehat=0	###e(beta,phi) in formula (12)
      log.mat=0
      for (i in (p+1):t){
        tmp.n=num_obs[i]
        M_t= indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
        eps= f[i,1:tmp.n]- M_t%*% phihat
        psi_matinv= psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
        epart= t(eps)%*%psi_matinv %*%eps
        ehat= epart+ehat
        log.mat=log.mat+ log(det(psi_matinv))
      }
      l=-sum(num_obs[(p+1):t])/2*log(ehat)+1/2*log.mat	### the number defined in formula (14)
      return(-l);
    }
    para4= exp(-1)
    result4=optim(para4,pdf4, lower=0.001, upper=0.999, method='L-BFGS-B')
    mat1= matrix(0, p*(df_b+1), p*(df_b+1));
    mat2= matrix(0, p*(df_b+1), 1);
    psi_invall=matrix(0,(n+1)*(t-1),(n+1))
    psi=result4$par
    for(i in (p+1):t){
      tmp.n=num_obs[i]
      psi_mat= matrix(1, tmp.n, tmp.n)
      for(k in 1:tmp.n){
        for(j in 1:tmp.n){
          psi_mat[k,j]=psi^(abs(x_pos[i,j]-x_pos[i,k]))
        }
      }
      psi_matinv2=solve(diag(weight(x_pos[i,1:tmp.n]))%*%psi_mat%*%diag(weight(x_pos[i,1:tmp.n])))
      psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]= psi_matinv2
      M_t= indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
      tmp_mat= t(M_t) %*% psi_matinv2;
      mat1= mat1+ tmp_mat %*% M_t;
      mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];
    }
    phihat= solve(mat1)%*%mat2
    phihat_func=matrix(0,p,(1+2*grid))
    for(i in 1:p){
      phihat_func[i,]= cbind(rep(1,2*grid+1),b_grid_sp)%*%phihat[(1+(coef_b*(i-1))):(i*coef_b)];		###b-spline coefficient of estimated phi_1
    }
    ehat=0		###e(beta,phi) in formula (12)
    for (i in (p+1):t){
      tmp.n=num_obs[i]
      M_t= indep.p[(i-p-1)*(n+1)+(1:tmp.n),]
      eps= f[i,1:tmp.n]- M_t%*% phihat
      psi_matinv=psi_invall[(i-2)*(n+1)+(1:tmp.n),1:tmp.n]
      epart= t(eps)%*%psi_matinv %*%eps
      ehat= epart+ehat
    }
    rho_hat= -log(psi)
    sigma_hat=sqrt(ehat/sum(num_obs[(p+1):t]))
  }

  est_cfarh <- list(phi_coef=t(matrix(phihat,df_b+1,p)),phi_func=phihat_func, rho=rho_hat, sigma=sigma_hat)
  return(est_cfarh)
}


#' Prediction of CFAR Processes
#'
#' Prediction of CFAR processes.
#' @param model CFAR model.
#' @param f the functional time series data.
#' @param m the forecast horizon.
#' @references
#' Liu, X., Xiao, H., and Chen, R. (2016) Convolutional autoregressive models for functional time series. \emph{Journal of Econometrics}, 194, 263-282.
#' @return The function returns a prediction of the CFAR process.
#' @examples
#' phi_func= function(x)
#' {
#'  	return(dnorm(x,mean=0,sd=0.1))
#' }
#' y=g_cfar1(100,5,phi_func)
#' f_grid=y$cfar
#' index=seq(1,1001,by=10)
#' f=f_grid[,index]
#' est=est_cfar(f,1)
#' pred=p_cfar(est,f,1)
#' @export
p_cfar <- function(model, f, m=3){
  ###p_cfar: this function gives us the prediction of functional time series. There are 3 input variables:
  ###1. model is the estimated model obtained from est_cfar function
  ###2. f is the matrix which contains the data
  ###3. m is the forecast horizon

  ####################
  ###parameter setting
  ##t= 100;	### length of time
  ##grid=1000;	### the number of grid points used to construct the functional time series X_t and epsilon_t
  #p is the cfar order
  ##n=100		### number of observations for each X_t, in the paper it is 'N'
  ###############################################################################
  if(!is.matrix(f))f <- as.matrix(f)
  ###################
  ###Estimation
  n=dim(f)[2]-1
  t=dim(f)[1]
  p=dim(model$phi_coef)[1]
  grid=(dim(model$phi_func)[2]-1)/2
  pred=matrix(0,t+m,grid+1)
  x_grid=seq(0,1,by=1/grid)
  x=seq(0,1,by=1/n)
  for(i in (t-p):t){
    pred[i,]= approx(x,f[i,],xout=x_grid,rule=2,method='linear')$y
  }
  phihat_func=model$phi_func
  for (j in 1:m){
    for (i in 1:(grid+1)){
      pred[j+t,i]= sum(phihat_func[1:p,seq((i+grid),i, by=-1)]*pred[(j+t-1):(j+t-p),])/grid
    }
  }
  pred_cfar=pred[(t+1):(t+m),]
  return(pred_cfar)
}


#' Partial Curve Prediction of CFAR Processes
#'
#' Partial prediction for CFAR processes. t curves are given and we want to predit the curve at time t+1, but we know the first n observations in the curve, to predict the n+1 observation.
#' @param model CFAR model.
#' @param f the functional time series data.
#' @param new.obs the given first \code{n} observations.
#' @references
#' Liu, X., Xiao, H., and Chen, R. (2016) Convolutional autoregressive models for functional time series. \emph{Journal of Econometrics}, 194, 263-282.
#' @return The function returns a prediction of the CFAR process.
#' @export
p_cfar_part  <- function(model, f, new.obs){
  ##p-cfar_part: this function gives us the prediction of x_t(s), s>s_0, give x_1,\ldots, x_{t-1}, and x(s), s <s_0. There are three input variables.
  ##1. model is the estimated model obtained from est_cfar function
  ##2. f is the matrix which contains the data
  ##3. new.obs is x_t(s), s>s_0.
  t=dim(f)[1]
  p=dim(model$phi_coef)[1]
  grid=(dim(model$phi_func)[2]-1)/2
  f_grid=matrix(0,t,grid+1)
  x_grid=seq(0,1,by=1/grid)
  n=dim(f)[2]-1
  index=seq(1,grid+1,by=grid/n)
  x=seq(0,1,by=1/n)
  for(i in (t-p):t){
    f_grid[i,]=approx(x,f[i,],xout=x_grid,rule=2,method='linear')$y
  }
  pred=matrix(0,grid+1,1)
  phihat_func=model$phi_func
  for (i in 1:(grid+1)){
    pred[i]= sum(phihat_func[1:p,seq((i+grid),i, by=-1)]*f_grid[t:(t-p+1),])/grid
  }
  pred_new=matrix(0,n+1,0)
  pred_new[1]=pred[1]
  rho=model$rho
  pred_new[2:(n+1)]=pred[index[1:n]]+(new.obs[1:n]-pred[index[1:n]])*exp(-rho/n)
  return(pred_new)
}


#' Generate Univariate SETAR Models
#'
#' Generate univariate SETAR model for up to 3 regimes.
#' @param nob number of observations.
#' @param arorder AR-order for each regime. The length of arorder controls the number of regimes.
#' @param phi a 3-by-p matrix. Each row contains the AR coefficients for a regime.
#' @param d delay for threshold variable.
#' @param thr threshold values.
#' @param sigma standard error for each regime.
#' @param cnst constant terms.
#' @param ini burn-in period.
#' @return uTAR.sim returns a list with components:
#' \item{series}{a time series following SETAR model.}
#' \item{at}{innovation of the time seres.}
#' \item{arorder}{AR-order for each regime.}
#' \item{thr}{threshold value.}
#' \item{phi}{a 3-by-p matrix. Each row contains the AR coefficients for a regime.}
#' \item{cnst}{constant terms}
#' \item{sigma}{standard error for each regime.}
#' @examples
#' arorder=rep(1,2)
#' ar.coef=matrix(c(0.7,-0.8),2,1)
#' y=uTAR.sim(100,arorder,ar.coef,1,0)
#' @export
"uTAR.sim" <- function(nob,arorder,phi,d=1,thr=c(0,0),sigma=c(1,1,1),cnst=rep(0,3),ini=500){
  p <- max(arorder)
  nT <- nob+ini
  ist <- max(p,d)+1
  et <- rnorm(nT)
  nregime <- length(arorder)
  if(nregime > 3)nregime <- 3
  ## For two-regime models, regime 3 is set to the same as regime 2.
  if(nregime == 2)arorder=c(arorder[1],arorder[2],arorder[2])
  k1 <- length(cnst)
  if(k1 < 3)cnst <- c(cnst,rep(cnst[k1],(3-k1)))
  k1 <- length(sigma)
  if(k1 < 3) sigma <- c(sigma,rep(sigma[k1],(3-k1)))
  if(!is.matrix(phi))phi <- as.matrix(phi)
  k1 <- nrow(phi)
  if(k1 < 3){
    for (j in (k1+1):3){
      phi <- rbind(phi,phi[k1,])
    }
  }
  k1 <- length(thr)
  if(k1 < 2) thr=c(thr,10^7)
  if(nregime == 2)thr[2] <- 10^7
  zt <- et[1:ist]*sigma[1]
  at <- zt
  ##
  for (i in ist:nT){
    if(zt[i-d] <= thr[1]){
      at[i] = et[i]*sigma[1]
      tmp = cnst[1]
      if(arorder[1] > 0){
        for (j in 1:arorder[1]){
          tmp <- tmp + zt[i-j]*phi[1,j]
        }
      }
      zt[i] <- tmp + at[i]
    }
    else{ if(zt[i-d] <= thr[2]){
      at[i] <- et[i]*sigma[2]
      tmp <- cnst[2]
      if(arorder[2] > 0){
        for (j in 1:arorder[2]){
          tmp <- tmp + phi[2,j]*zt[i-j]
        }
      }
      zt[i] <- tmp + at[i]
    }
      else{ at[i] <- et[i]*sigma[3]
      tmp <- cnst[3]
      if(arorder[3] > 0){
        for (j in 1:arorder[3]){
          tmp <- tmp + phi[3,j]*zt[i-j]
        }
      }
      zt[i] <- tmp + at[i]
      }
    }

  }
  uTAR.sim <- list(series = zt[(ini+1):nT], at=at[(ini+1):nT], arorder=arorder,thr=thr,phi=phi,cnst=cnst,sigma=sigma)
}



#' Estimation of a Univariate Two-Regime SETAR Model
#'
#' Estimation of a univariate two-regime SETAR model, including threshold value.
#' The procedure of Li and Tong (2016) is used to search for the threshold.
#' @param y a vector of time series.
#' @param p1,p2 AR-orders of regime 1 and regime 2.
#' @param d delay for threshold variable, default is 1.
#' @param thrV threshold variable. If thrV is not null, it must have the same length as that of y.
#' @param Trim lower and upper quantiles for possible threshold values.
#' @param k0 the maximum number of threshold values to be evaluated. If the sample size is large (> 3000), then k0 = floor(nT*0.5). The default is k0=300. But k0 = floor(nT*0.8) if nT < 300.
#' @param include.mean a logical value indicating whether constant terms are included.
#' @param thrQ lower and upper quantiles to search for threshold value.
#' @references
#' Li, D., and Tong. H. (2016) Nested sub-sample search algorithm for estimation of threshold models. \emph{Statisitca Sinica}, 1543-1554.
#' @return uTAR returns a list with components:
#' \item{data}{the data matrix, y.}
#' \item{arorder}{AR orders of regimes 1 and 2.}
#' \item{delay}{the delay for threshold variable.}
#' \item{residuals}{estimated innovations.}
#' \item{coef}{a 2-by-(p+1) matrices. The first row shows the estimation results in regime 1, and the second row shows these in regime 2.}
#' \item{sigma}{estimated innovational covariance matrices of regimes 1 and 2.}
#' \item{nobs}{numbers of observations in regimes 1 and 2.}
#' \item{model1,model2}{estimated models of regimes 1 and 2.}
#' \item{thr}{threshold value.}
#' \item{D}{a set of threshold values.}
#' \item{RSS}{RSS}
#' \item{AIC}{AIC value}
#' \item{cnst}{logical values indicating whether the constant terms are included in regimes 1 and 2.}
#' \item{sresi}{standardized residuals.}
#' @examples
#' phi=t(matrix(c(-0.3, 0.5,0.6,-0.3),2,2))
#' y=uTAR.sim(nob=2000, arorder=c(2,2), phi=phi, d=2, thr=0.2, cnst=c(1,-1),sigma=c(1, 1))
#' est=uTAR(y=y$series,p1=2,p2=2,d=2,k0=50,thrQ=c(0,1),Trim=c(0.1,0.9),include.mean=TRUE)
#' @export
"uTAR" <- function(y,p1,p2,d=1,thrV=NULL,Trim=c(0.1,0.9),k0=300,include.mean=TRUE,thrQ=c(0,1)){

  if(is.matrix(y))y=y[,1]
  p = max(p1,p2)
  if(p < 1)p=1
  ist=max(p,d)+1
  nT <- length(y)
  ### regression framework
  if(k0 > nT)k0 <- floor(nT*0.8)
  if(nT > 3000)k0 <- floor(nT*0.5)
  ### built in checking
  nT1 <- nT
  if(!is.null(thrV))nT1 <- length(thrV)
  if(nT != nT1){
    cat("Input error at thrV. Reset to a SETAR model","\n")
    thrV <- y
  }
  yt <- y[ist:nT]
  nobe <- nT-ist+1
  x1 <- NULL
  for (i in 1:p){
    x1=cbind(x1,y[(ist-i):(nT-i)])
  }
  if(length(thrV) < nT){
    tV <- y[(ist-d):(nT-d)]
  }else{
    tV <- thrV[(ist-d):(nT-d)]
  }
  D <- NeSS(yt,x1,x2=NULL,thrV=tV,Trim=Trim,k0=k0,include.mean=include.mean,thrQ=thrQ)
  #cat("Threshold candidates: ",D,"\n")
  ll = length(D)
  ### house keeping
  RSS=NULL
  resi <- NULL
  sresi <- NULL
  m1a <- NULL
  m1b <- NULL
  Size <- NULL
  Phi <- matrix(0,2,p+1)
  Sigma <- NULL
  infc <- NULL
  thr <- NULL
  ##
  if(ll > 0){
    for (i in 1:ll){
      thr = D[i]
      idx = c(1:nobe)[tV <= thr]
      X = x1[,c(1:p1)]
      if(include.mean){X=cbind(rep(1,nobe),X)}
      m1a <- lm(yt~-1+.,data=data.frame(X),subset=idx)
      X <- x1[,c(1:p2)]
      if(include.mean){X=cbind(rep(1,nobe),X)}
      m1b <- lm(yt~-1+.,data=data.frame(X),subset=-idx)
      RSS=c(RSS,sum(m1a$residuals^2,m1b$residuals^2))
    }
    ##cat("RSS: ",RSS,"\n")
    j=which.min(RSS)
    thr <- D[j]
    cat("Estimated Threshold: ",thr,"\n")
    ### Final estimation results
    resi = rep(0,nobe)
    sresi <- rep(0,nobe)
    idx = c(1:nobe)[tV <= thr]
    n1 <- length(idx)
    n2 <- length(yt)-n1
    Size <- c(n1,n2)
    X = x1[,c(1:p1)]
    if(include.mean){X=cbind(rep(1,(nT-ist+1)),X)}
    m1a <- lm(yt~-1+.,data=data.frame(X),subset=idx)
    resi[idx] <- m1a$residuals
    X <- as.matrix(X)
    cat("Regime 1: ","\n")
    Phi[1,1:ncol(X)] <- m1a$coefficients
    coef1 <- summary(m1a)
    sigma1 <- coef1$sigma
    sresi[idx] <- m1a$residuals/sigma1
    print(coef1$coefficients)
    sig1 <- sigma1^2*(n1-ncol(X))/n1
    AIC1 <- n1*log(sig1) + 2*ncol(X)
    cat("nob1, sigma1 and AIC: ",c(n1, sigma1, AIC1), "\n")
    X <- x1[,c(1:p2)]
    if(include.mean){X=cbind(rep(1,nobe),X)}
    m1b <- lm(yt~-1+.,data=data.frame(X),subset=-idx)
    resi[-idx] <- m1b$residuals
    X <- as.matrix(X)
    cat("Regime 2: ","\n")
    Phi[2,1:ncol(X)] <- m1b$coefficients
    coef2 <- summary(m1b)
    sigma2 <- coef2$sigma
    sresi[-idx] <- m1b$residuals/sigma2
    print(coef2$coefficients)
    sig2 <- sigma2^2*(n2-ncol(X))/n2
    AIC2 <- n2*log(sig2) + 2*ncol(X)
    cat("nob2, sigma2 and AIC: ",c(n2,sigma2, AIC2),"\n")
    fitted.values <- yt-resi
    Sigma <- c(sigma1,sigma2)
    ### pooled estimate
    sigmasq = (n1*sigma1^2+n2*sigma2^2)/(n1+n2)
    AIC <- AIC1+AIC2
    cat(" ","\n")
    cat("Overal MLE of sigma: ",sqrt(sigmasq),"\n")
    cat("Overall AIC: ",AIC,"\n")
  }else{cat("No threshold found: Try again with a larger k0.","\n")}

  uTAR <- list(data=y,arorder=c(p1,p2),delay=d,residuals = resi, coefs = Phi, sigma=Sigma,
               nobs = Size, model1 = m1a, model2 = m1b,thr=thr, D=D, RSS=RSS,AIC = AIC,
               cnst = rep(include.mean,2), sresi=sresi)
}


#' Search for Threshold Value of A Two-Regime SETAR Model
#'
#' Search for the threshold of a SETAR model for a given range of candidates for threshold values,
#' and perform recursive LS estimation.
#' The program uses a grid to search for threshold value.
#' It is a conservative approach, but might be more reliable than the Li and Tong (2016) procedure.
#' @param y a vector of time series.
#' @param p1,p2 AR-orders of regime 1 and regime 2.
#' @param d delay for threshold variable, default is 1. (Apply to the external threshold variable too.)
#' @param thrV threshold variable. if it is not null, thrV must have the same length as that of y.
#' @param thrQ lower and upper limits for the possible threshold values.
#' @param Trim lower and upper trimming to control the sample size in each regime.
#' @param include.mean a logical value for including constant term.
#' @return uTAR.grid returns a list with components:
#' \item{data}{the data matrix, y.}
#' \item{arorder}{AR orders of regimes 1 and 2.}
#' \item{residuals}{estimated innovations.}
#' \item{coefs}{a 2-by-(p+1) matrices. The first row shows the estimation results in regime 1, and the second row shows these in regime 2.}
#' \item{sigma}{estimated innovational covariance matrices of regimes 1 and 2.}
#' \item{nobs}{numbers of observations in regimes 1 and 2.}
#' \item{delay}{the delay for threshold variable.}
#' \item{model1,model2}{estimated models of regimes 1 and 2.}
#' \item{cnst}{logical values indicating whether the constant terms are included in regimes 1 and 2.}
#' \item{thr}{threshold value.}
#' \item{D}{a set of possible threshold values.}
#' \item{RSS}{residual sum of squares.}
#' \item{information}{information criterion.}
#' \item{sresi}{standardized residuals.}
#' @examples
#' phi=t(matrix(c(-0.3, 0.5,0.6,-0.3),2,2))
#' y=uTAR.sim(nob=2000, arorder=c(2,2), phi=phi, d=2, thr=0.2, cnst=c(1,-1),sigma=c(1, 1))
#' est1=uTAR.grid(y=y$series, p1=2, p2=2, d=2, thrQ=c(0,1),Trim=c(0.1,0.9),include.mean=TRUE)
#' @export
"uTAR.grid" <- function(y,p1,p2,d=1,thrV=NULL,thrQ=c(0,1),Trim=c(0.1,0.9),include.mean=T){
  if(is.matrix(y))y=y[,1]
  p = max(p1,p2)
  if(p < 1)p=1
  ist=max(p,d)+1
  nT <- length(y)
  yt <- y[ist:nT]
  nobe <- nT-ist+1
  Phi <- matrix(0,2,p+1)
  x1 <- NULL
  for (i in 1:p1){
    x1=cbind(x1,y[(ist-i):(nT-i)])
  }
  x1 <- as.matrix(x1)
  if(include.mean){x1 <- cbind(rep(1,nobe),x1)}
  k1 <- ncol(x1)
  x2 <- NULL
  for (i in 1:p2){
    x2 <- cbind(x2,y[(ist-i):(nT-i)])
  }
  x2 <- as.matrix(x2)
  if(include.mean){x2 <- cbind(rep(1,nobe),x2)}
  k2 <- ncol(x2)
  #
  if(length(thrV) < nT){
    tV <- y[(ist-d):(nT-d)]
  }else{
    tV <- thrV[ist:nT]
  }
  ### narrow down possible threshold range
  q1=c(min(tV),max(tV))
  if(thrQ[1] > 0)q1[1] <- quantile(tV,thrQ[1])
  if(thrQ[2] < 1)q1[2] <- quantile(tV,thrQ[2])
  idx <- c(1:length(tV))[tV < q1[1]]
  jdx <- c(1:length(tV))[tV > q1[2]]
  trm <- c(idx,jdx)
  if(length(trm)>0){
    yt <- yt[-trm]
    if(k1 == 1){x1 <- x1[-trm]
    }else{x1 <- x1[-trm,]}
    if(k2 == 1){x2 <- x2[-trm]
    }else{
      x2 <- x2[-trm,]}
    tV <- tV[-trm]
  }
  ### sorting
  DD <- sort(tV,index.return=TRUE)
  D <- DD$x
  Dm <- DD$ix
  yt <- yt[Dm]
  if(k1 == 1){x1 <- x1[Dm]
  }else{x1 <- x1[Dm,]}
  if(k2 == 1){x2 <- x2[Dm]
  }else{x2 <- x2[Dm,]}
  nobe <- length(yt)
  q2 = quantile(tV,Trim)
  jst <- max(c(1:nobe)[D < q2[1]])
  jend <- min(c(1:nobe)[D > q2[2]])
  cat("jst, jend: ",c(jst,jend),"\n")
  nthr <- jend-jst+1
  RSS=NULL
  ### initial estimate with thrshold D[jst]
  if(k1 == 1){
    X <- as.matrix(x1[1:jst])
  }else{ X <- x1[1:jst,]}
  Y <- as.matrix(yt[1:jst])
  XpX = t(X)%*%X
  XpY = t(X)%*%Y
  Pk1 <- solve(XpX)
  beta1 <- Pk1%*%XpY
  resi1 <- Y - X%*%beta1
  if(k2 == 1){
    X <- as.matrix(x2[(jst+1):nobe])
  }else{ X <- x2[(jst+1):nobe,]}
  Y <- as.matrix(yt[(jst+1):nobe])
  XpX <- t(X)%*%X
  XpY <- t(X)%*%Y
  Pk2 <- solve(XpX)
  beta2 <- Pk2%*%XpY
  resi2 <- Y - X%*%beta2
  RSS <- sum(c(resi1^2,resi2^2))
  #### recursive
  #### Moving one observations from regime 2 to regime 1
  for (i in (jst+1):jend){
    if(k1 == 1){xk1 <- matrix(x1[i],1,1)
    }else{
      xk1 <- matrix(x1[i,],k1,1)
    }
    if(k2 == 1){xk2 <- matrix(x2[i],1,1)
    }else{
      xk2 <- matrix(x2[i,],k2,1)
    }
    tmp <- Pk1%*%xk1
    deno <- 1 + t(tmp)%*%xk1
    er <- t(xk1)%*%beta1 - yt[i]
    beta1 <- beta1 - tmp*c(er)/c(deno[1])
    Pk1 <- Pk1 - tmp%*%t(tmp)/c(deno[1])
    if(k1 == 1){X <- as.matrix(x1[1:i])
    }else{X <- x1[1:i,]}
    resi1 <- yt[1:i] - X%*%beta1
    tmp <- Pk2 %*% xk2
    deno <- t(tmp)%*%xk2 - 1
    er2 <- t(xk2)%*%beta2 - yt[i]
    beta2 <- beta2 - tmp*c(er2)/c(deno[1])
    Pk2 <- Pk2 - tmp%*%t(tmp)/c(deno[1])
    if(k2 == 1){X <- as.matrix(x2[(i+1):nobe])
    }else{X <- x2[(i+1):nobe,]}
    resi2 <- yt[(i+1):nobe] - X%*%beta2
    RSS <- c(RSS, sum(c(resi1^2,resi2^2)))
  }
  ##cat("RSS: ",RSS,"\n")
  j=which.min(RSS)+(jst-1)
  thr <- D[j]
  cat("Estimated Threshold: ",thr,"\n")
  ### Final estimation results
  resi = rep(0,nobe)
  sresi <- rep(0,nobe)
  idx <- c(1:j)
  n1 <- j
  n2 <- nobe-j
  X <- x1
  m1a <- lm(yt~-1+.,data=data.frame(x1),subset=idx)
  resi[idx] <- m1a$residuals
  cat("Regime 1: ","\n")
  Phi[1,1:k1] <- m1a$coefficients
  coef1 <- summary(m1a)
  sresi[idx] <- m1a$residuals/coef1$sigma
  print(coef1$coefficients)
  cat("nob1 & sigma1: ",c(n1, coef1$sigma), "\n")
  m1b <- lm(yt~-1+.,data=data.frame(x2),subset=-idx)
  resi[-idx] <- m1b$residuals
  Phi[2,1:k2] <- m1b$coefficients
  cat("Regime 2: ","\n")
  coef2 <- summary(m1b)
  sresi[-idx] <- m1b$residuals/coef2$sigma
  print(coef2$coefficients)
  cat("nob2 & sigma2: ",c(n2,coef2$sigma),"\n")
  fitted.values <- yt-resi
  ### pooled estimate
  sigma2 = (n1*coef1$sigma+n2*coef2$sigma)/(n1+n2)
  d1 = log(sigma2)
  aic = d1 + 2*(p1+p2)/nT
  bic = d1 + log(nT)*(p1+p2)/nT
  hq = d1 + 2*log(log(nT))*(p1+p2)/nT
  infc = c(aic,bic,hq)
  cat(" ","\n")
  cat("The fitted TAR model ONLY uses data in the specified threshold range!!!","\n")
  cat("Overal MLE of sigma: ",sqrt(sigma2),"\n")
  cat("Overall information criteria(aic,bic,hq): ",infc,"\n")
  uTAR.grid <- list(data=y,arorder=c(p1,p2), residuals = resi, coefs = Phi,
                    sigma=c(coef1$sigma,coef2$sigma), nobs=c(n1,n2),delay=d,model1 = m1a,cnst = rep(include.mean,2),
                    model2 = m1b,thr=thr, D=D, RSS=RSS,information=infc,sresi=sresi)
}

#' General Estimation of TAR Models
#'
#' General estimation of TAR models with known threshold values.
#' It perform LS estimation of a univariate TAR model, and can handle multiple regimes.
#' @param y time series.
#' @param arorder AR order of each regime. The number of regime is the length of arorder.
#' @param thr given threshold(s). There are k-1 threshold for a k-regime model.
#' @param d delay for threshold variable, default is 1.
#' @param thrV external threshold variable if any. If it is not NULL, thrV must have the same length as that of y.
#' @param include.mean a logical value indicating whether constant terms are included. Default is TRUE.
#' @param output a logical value for output. Default is TRUE.
#' @return uTAR.est returns a list with components:
#' \item{data}{the data matrix, y.}
#' \item{k}{the dimension of y.}
#' \item{arorder}{AR orders of regimes 1 and 2.}
#' \item{coefs}{a (p*k+1)-by-(2k) matrices. The first row shows the estimation results in regime 1, and the second row shows these in regime 2.}
#' \item{sigma}{estimated innovational covariance matrices of regimes 1 and 2.}
#' \item{thr}{threshold value.}
#' \item{residuals}{estimated innovations.}
#' \item{sresi}{standardized residuals.}
#' \item{nobs}{numbers of observations in different regimes.}
#' \item{delay}{delay for threshold variable.}
#' \item{cnst}{logical values indicating whether the constant terms are included in different regimes.}
#' \item{AIC}{AIC value.}
#' @examples
#' phi=t(matrix(c(-0.3, 0.5,0.6,-0.3),2,2))
#' y=uTAR.sim(nob=2000, arorder=c(2,2), phi=phi, d=2, thr=0.2, cnst=c(1,-1),sigma=c(1, 1))
#' thr.est=uTAR.grid(y=y$series, p1=2, p2=2, d=2, thrQ=c(0,1),Trim=c(0.1,0.9))
#' est=uTAR.est(y=y$series, arorder=c(2,2), thr=thr.est$thr, d=2)
#' @export
"uTAR.est" <- function(y,arorder=c(1,1),thr=c(0),d=1,thrV=NULL,include.mean=c(TRUE,TRUE),output=TRUE){
  if(is.matrix(y))y <- y[,1]
  p <- max(arorder)
  if(p < 1)p <-1
  k <- length(arorder)
  ist <- max(p,d)+1
  nT <- length(y)
  Y <- y[ist:nT]
  X <- NULL
  for (i in 1:p){
    X <- cbind(X,y[(ist-i):(nT-i)])
  }
  X <- as.matrix(X)
  nobe <- nT-ist+1
  if(length(thrV) < nT){
    thrV <- y[(ist-d):(nT-d)]
    if(output) cat("Estimation of a SETAR model: ","\n")
  }else{thrV <- thrV[ist:nT]}
  k1 <- length(thr)
  ### house keeping
  coefs <- NULL
  sigma <- NULL
  nobs <- NULL
  resi <- rep(0,nobe)
  sresi <- rep(0,nobe)
  npar <- NULL
  AIC <- 99999
  #
  if(k1 > (k-1)){
    cat("Only the first: ",k1," thresholds are used.","\n")
    THr <- min(thrV[1:nobe])-1
    THr <- c(THr,thr[1:(k-1)],(max(thrV[1:nobe])+1))
  }
  if(k1 < (k-1)){cat("Too few threshold values are given!!!","\n")}
  if(length(include.mean) != k)include.mean=rep("TRUE",k)
  if(k1 == (k-1)){
    THr <- c((min(thrV[1:nobe])-1),thr,c(max(thrV[1:nobe])+1))
    if(output) cat("Thresholds used: ",thr,"\n")
    for (j in 1:k){
      idx <- c(1:nobe)[thrV <= THr[j]]
      jdx <- c(1:nobe)[thrV > THr[j+1]]
      jrm <- c(idx,jdx)
      n1 <- nobe-length(jrm)
      nobs <- c(nobs,n1)
      kdx <- c(1:nobe)[-jrm]
      x1 <- as.matrix(X)
      if(p > 1)x1 <- x1[,c(1:arorder[j])]
      cnst <- rep(1,nobe)
      if(include.mean[j]){x1 <- cbind(cnst,x1)}
      np <- ncol(x1)
      npar <- c(npar,np)
      beta <- rep(0,p+1)
      sig <- NULL
      if(n1 <= ncol(x1)){
        cat("Regime ",j," does not have sufficient observations!","\n")
      }
      else{
        mm <- lm(Y~-1+.,data=data.frame(x1),subset=kdx)
        c1 <- summary(mm)
        if(output){
          cat("Estimation of regime: ",j," with sample size: ",n1,"\n")
          print(c1$coefficients)
          cat("Sigma estimate: ",c1$sigma,"\n")
        }
        resi[kdx] <- c1$residuals
        sresi[kdx] <- c1$residuals/c1$sigma
        sig <- c1$sigma
        beta[1:np] <- c1$coefficients[,1]
        coefs <- rbind(coefs,beta)
        sigma <- c(sigma,sig)
      }
    }
    AIC <- 0
    for (j in 1:k){
      sig <- sigma[j]
      n1 <- nobs[j]
      np <- npar[j]
      if(!is.null(sig)){s2 <- sig^2*(n1-np)/n1
      AIC <- AIC + log(s2)*n1 + 2*np
      }
    }
    if(output) cat("Overal AIC: ",AIC,"\n")
  }
  uTAR.est <- list(data=y,k = k, arorder=arorder,coefs=coefs,sigma=sigma,thr=thr,residuals=resi,
                   sresi=sresi,nobs=nobs,delay=d,cnst=include.mean, AIC=AIC)
}



#' Prediction of A Fitted Univariate TAR Model
#'
#' Prediction of a fitted univariate TAR model.
#' @param model univariate TAR model.
#' @param orig forecast origin.
#' @param h forecast horizon.
#' @param iterations number of iterations.
#' @param ci confidence level.
#' @param output a logical value for output, default is TRUE.
#' @return uTAR.pred returns a list with components:
#' \item{model}{univariate TAR model.}
#' \item{pred}{prediction.}
#' \item{Ysim}{fitted y.}
#' @examples
#' phi=t(matrix(c(-0.3, 0.5,0.6,-0.3),2,2))
#' y=uTAR.sim(nob=2000, arorder=c(2,2), phi=phi, d=2, thr=0.2, cnst=c(1,-1),sigma=c(1, 1))
#' thr.est=uTAR.grid(y=y$series, p1=2, p2=2, d=2, thrQ=c(0,1),Trim=c(0.1,0.9))
#' est=uTAR.est(y=y$series, arorder=c(2,2), thr=thr.est$thr, d=2)
#' uTAR.pred(mode=est, orig=2000,h=1,iteration=100,ci=0.95,output=TRUE)
#' @export
"uTAR.pred" <- function(model,orig,h=1,iterations=3000,ci=0.95,output=TRUE){
  ## ci: probability of pointwise confidence interval coefficient
  ### only works for self-exciting uTAR models.
  ###
  y <- model$data
  arorder <- model$arorder
  coefs <- model$coefs
  sigma <- model$sigma
  thr <- c(model$thr)
  include.mean <- model$cnst
  d <- model$delay
  p <- max(arorder)
  k <- length(arorder)
  nT <- length(y)
  if(orig < 1)orig <- nT
  if(orig > nT) orig <- nT
  if(h < 1) h <- 1
  ### compute the predictions. Use simulation if the threshold is predicted
  Ysim <- matrix(0,h,iterations)
  et <- rnorm(h*iterations)
  for (it in 1:iterations){
    etcnt <- (it-1)*h
    yp <- y[1:orig]
    for (ii in 1:h){
      t <- orig+ii
      thd <- yp[(t-d)]
      JJ <- 1
      for (j in 1:(k-1)){
        if(thd > thr[j]){JJ <- j+1}
      }
      ino <- etcnt+ii
      at <- et[ino]*sigma[JJ]
      x <- NULL
      pJ <- arorder[JJ]
      npar <- pJ
      if(include.mean[JJ]){ x <- 1
      npar <- npar+1
      }
      for (i in 1:pJ){
        x <- c(x,yp[(t-i)])
      }
      yhat <- sum(x*coefs[JJ,1:npar])+at
      yp <- rbind(yp,yhat)
      Ysim[ii,it] <- yhat
    }
  }
  ### summary
  pred <- NULL
  CI <- NULL
  pr <- (1-ci)/2
  prob <- c(pr, 1-pr)
  for (ii in 1:h){
    pred <- rbind(pred,mean(Ysim[ii,]))
    int <- quantile(Ysim[ii,],prob=prob)
    CI <- rbind(CI,int)
  }
  if(output){
    cat("Forecast origin: ",orig,"\n")
    cat("Predictions: 1-step to ",h,"-step","\n")
    Pred <- cbind(1:h,pred)
    colnames(Pred) <- c("step","forecast")
    print(Pred)
    bd <- cbind(c(1:h),CI)
    colnames(bd) <- c("step", "Lowb","Uppb")
    cat("Pointwise ",ci*100," % confident intervals","\n")
    print(bd)
  }
  uTAR.pred <- list(data=y,pred = pred,Ysim=Ysim)
}

"NeSS" <- function(y,x1,x2=NULL,thrV=NULL,Trim=c(0.1,0.9),k0=50,include.mean=TRUE,thrQ=c(0,1),score="AIC"){
  ### Based on method proposed by Li and Tong (2016) to narrow down the candidates for
  ### threshold value.
  ### y: dependent variable (can be multivariate)
  ### x1: the regressors that allow for coefficients to change
  ### x2: the regressors for explanatory variables
  ### thrV: threshold variable. The default is simply the time index
  ### Trim: quantiles for lower and upper bound of threshold
  ### k0: the maximum number of candidates for threshold at the output
  ### include.mean: switch for including the constant term.
  ### thrQ: lower and upper quantiles to search for threshold
  ####
  #### return: a set of possible threshold values
  ####
  if(!is.matrix(y))y <- as.matrix(y)
  ky <- ncol(y)
  if(!is.matrix(x1))x1 <- as.matrix(x1)
  if(length(x2) > 0){
    if(!is.matrix(x2))x2 <- as.matrix(x2)
  }
  n <- nrow(y)
  n1 <- nrow(x1)
  n2 <- n1
  if(!is.null(x2))n2 <- nrow(x2)
  n <- min(n,n1,n2)
  ##
  k1 <- ncol(x1)
  k2 <- 0
  if(!is.null(x2)) k2 <- ncol(x2)
  k <- k1+k2
  ##
  if(length(thrV) < 1){thrV <- c(1:n)
  }else{thrV <- thrV[1:n]}
  ### set threshold range
  idx <- NULL
  jdx <- NULL
  if(thrQ[1] > 0){low <- quantile(thrV,thrQ[1])
  idx <- c(1:length(thrV))[thrV < low]
  }
  if(thrQ[2] < 1){upp <- quantile(thrV,thrQ[2])
  jdx <- c(1:length(thrV))[thrV > upp]
  }
  jrm <- c(idx,jdx)
  if(!is.null(jrm)){
    if(ky == 1){y <- as.matrix(y[-jrm])
    }else{y <- y[-jrm,]}
    if(k1 == 1){x1 <- as.matrix(x1[-jrm])
    }else{x1 <- x1[-jrm,]}
    if(k2 > 0){
      if(k2 == 1){x2 <- as.matrix(x2[-jrm])
      }else{x2 <- x2[-jrm,]}
    }
    thrV <- thrV[-jrm]
  }
  thr <- sort(thrV)
  n <- length(thrV)
  ### use k+2 because of including the constant term.
  n1 <- floor(n*Trim[1])
  if(n1 < (k+2)) n1 = k+2
  n2 <- floor(n*Trim[2])
  if((n-n2) > (k+2)){n2 = n -(k+2)}
  D <- thr[n1:n2]
  ####
  X=cbind(x2,x1)
  if(include.mean){k = k+1
  X=cbind(rep(1,n),X)}
  qprob=c(0.25,0.5,0.75)
  #### Scalar case
  if(ky == 1){
    X=data.frame(X)
    Y <- y[,1]
    m1 <- lm(Y~-1+.,data=X)
    R1 <- sum(m1$residuals^2)
    Jnr <- rep(R1,length(qprob))
    #cat("Jnr: ",Jnr,"\n")
    while(length(D) > k0){
      qr <- quantile(D,qprob)
      ##  cat("qr: ",qr,"\n")
      for (i in 1:length(qprob)){
        idx <- c(1:n)[thrV <= qr[i]]
        m1a <- lm(Y~-1+.,data=X,subset=idx)
        m1b <- lm(Y~-1+.,data=X,subset=-idx)
        Jnr[i] = - sum(c(m1a$residuals^2,m1b$residuals^2))
      }
      ##     cat("in Jnr: ",Jnr,"\n")
      if(Jnr[1] >= max(Jnr[-1])){
        D <- D[D <= qr[2]]
      }
      else{ if(Jnr[2] >= max(Jnr[-2])){
        D <- D[((qr[1] <= D) & (D <= qr[3]))]
      }
        else{ D <- D[D >= qr[2]]}
      }
      ## cat("n(D): ",length(D),"\n")
    }
  }
  #### Multivariate case
  if(ky > 1){
    ic = 2
    if(score=="AIC")ic=1
    m1 <- MlmNeSS(y,X,SD=FALSE,include.mean=include.mean)
    R1 <- m1$score[ic]
    Jnr <- rep(R1,length(qprob))
    while(length(D) > k0){
      qr <- quantile(D,qprob)
      for (i in 1:length(qprob)){
        idx <- c(1:n)[thrV <= qr[i]]
        m1a <- MlmNeSS(y,X,subset=idx,SD=FALSE,include.mean=include.mean)
        m1b <- MlmNeSS(y,X,subset=-idx,SD=FALSE,include.mean=include.mean)
        Jnr[i] <- - sum(m1a$score[ic],m1b$score[ic])
      }
      if(Jnr[1] >= max(Jnr[-1])){
        D <- D[D <= qr[2]]
      }
      else{ if(Jnr[2] >= max(Jnr[-2])){
        D <- D[D <= qr[3]]
        D <- D[D >= qr[1]]
      }
        else{ D <- D[D >= qr[2]]}
      }
    }
    ### end of Multivariate case
  }
  D
}

#' Generate Two-Regime (VAR) Models
#'
#' Generates two-regime multivariate vector auto-regressive models.
#' @param nob number of observations.
#' @param thr threshold value.
#' @param phi1 VAR coefficient matrix of regime 1.
#' @param phi2 VAR coefficient matrix of regime 2.
#' @param sigma1 innovational covariance matrix of regime 1.
#' @param sigma2 innovational covariance matrix of regime 2.
#' @param c1 constant vector of regime 1.
#' @param c2 constant vector of regime 2.
#' @param delay two elements (i,d) with "i" being the component index and "d" the delay for threshold variable.
#' @param ini burn-in period.
#' @return mTAR.sim returns a list with following components:
#' \item{series}{a time series following the two-regime multivariate VAR model.}
#' \item{at}{innovation of the time series.}
#' \item{threshold}{threshold value.}
#' \item{delay}{two elements (i,d) with "i" being the component index and "d" the delay for threshold variable.}
#' \item{n1}{number of observations in regime 1.}
#' \item{n2}{number of observations in regime 2.}
#' @examples
#' phi1=matrix(c(0.5,0.7,0.3,0.2),2,2)
#' phi2=matrix(c(0.4,0.6,0.5,-0.5),2,2)
#' sigma1=matrix(c(1,0,0,1),2,2)
#' sigma2=matrix(c(1,0,0,1),2,2)
#' c1=c(0,0)
#' c2=c(0,0)
#' delay=c(1,1)
#' y=mTAR.sim(100,0,phi1,phi2,sigma1,sigma2,c1,c2,delay,ini=500)
#' @export
"mTAR.sim" <- function(nob,thr,phi1,phi2,sigma1,sigma2=NULL,c1=NULL,c2=NULL,delay=c(1,1),ini=500){
  if(!is.matrix(phi1))phi1=as.matrix(phi1)
  if(!is.matrix(phi2))phi2=as.matrix(phi2)
  if(!is.matrix(sigma1))sigma1=as.matrix(sigma1)
  if(is.null(sigma2))sigma2=sigma1
  if(!is.matrix(sigma2))sigma2=as.matrix(sigma2)
  k1 <- nrow(phi1)
  k2 <- nrow(phi2)
  k3 <- nrow(sigma1)
  k4 <- nrow(sigma2)
  k <- min(k1,k2,k3,k4)
  if(is.null(c1))c1=rep(0,k)
  if(is.null(c2))c2=rep(0,k)
  c1 <- c1[1:k]
  c2 <- c2[1:k]
  p1 <- ncol(phi1)%/%k
  p2 <- ncol(phi2)%/%k
  ###
  "mtxroot" <- function(sigma){
    if(!is.matrix(sigma))sigma=as.matrix(sigma)
    sigma <- (t(sigma)+sigma)/2
    m1 <- eigen(sigma)
    P <- m1$vectors
    L <- diag(sqrt(m1$values))
    sigmah <- P%*%L%*%t(P)
    sigmah
  }
  ###
  s1 <- mtxroot(sigma1[1:k,1:k])
  s2 <- mtxroot(sigma2[1:k,1:k])
  nT <- nob+ini
  et <- matrix(rnorm(k*nT),nT,k)
  a1 <- et%*%s1
  a2 <- et%*%s2
  ###
  d <- delay[2]
  if(d <= 0)d <- 1
  p <- max(p1,p2,d)
  ### set starting values
  zt <- matrix(1,p,1)%*%matrix(c1,1,k)+a1[1:p,]
  resi <- matrix(0,nT,k)
  ist = p+1
  for (i in ist:ini){
    wk = rep(0,k)
    if(zt[(i-d),delay[1]] <= thr){
      resi[i,] <- a1[i,]
      wk <- wk + a1[i,] + c1
      if(p1 > 0){
        for (j in 1:p1){
          idx <- (j-1)*k
          phi <- phi1[,(idx+1):(idx+k)]
          wk <- wk + phi%*%matrix(zt[i-j,],k,1)
        }
      }
    }else{ resi[i,] <- a2[i,]
    wk <- wk+a2[i,] + c2
    if(p2 > 0){
      for (j in 1:p2){
        idx <- (j-1)*k
        phi <- phi2[,(idx+1):(idx+k)]
        wk <- wk + phi%*%matrix(zt[i-j,],k,1)
      }
    }
    }
    zt <- rbind(zt,c(wk))
  }
  #### generate data used
  n1 = 0
  for (i in (ini+1):nT){
    wk = rep(0,k)
    if(zt[(i-d),delay[1]] <= thr){
      n1 = n1+1
      resi[i,] <- a1[i,]
      wk <- wk + a1[i,] + c1
      if(p1 > 0){
        for (j in 1:p1){
          idx <- (j-1)*k
          phi <- phi1[,(idx+1):(idx+k)]
          wk <- wk + phi%*%matrix(zt[i-j,],k,1)
        }
      }
    }else{ resi[i,] <- a2[i,]
    wk <- wk+a2[i,] + c2
    if(p2 > 0){
      for (j in 1:p2){
        idx <- (j-1)*k
        phi <- phi2[,(idx+1):(idx+k)]
        wk <- wk + phi%*%matrix(zt[i-j,],k,1)
      }
    }
    }
    zt <- rbind(zt,c(wk))
  }
  mTAR.sim <- list(series = zt[(ini+1):nT,],at = resi[(ini+1):nT,],threshold=thr,
                   delay=delay, n1=n1, n2 = (nob-n1))
}


#' Estimation of a Multivariate Two-Regime SETAR Model
#'
#' Estimation of a multivariate two-regime SETAR model, including threshold.
#' The procedure of Li and Tong (2016) is used to search for the threshold.
#' @param y a (\code{nT}-by-\code{k}) data matrix of multivariate time series, where \code{nT} is the sample size and \code{k} is the dimension.
#' @param p1 AR-order of regime 1.
#' @param p2 AR-order of regime 2.
#' @param thr threshold variable. Estimation is needed if \code{thr} = NULL.
#' @param thrV vector of threshold variable. If it is not null, thrV must have the same sample size of that of y.
#' @param delay two elements (i,d) with "i" being the component and "d" the delay for threshold variable.
#' @param Trim lower and upper quantiles for possible threshold value.
#' @param k0 the maximum number of threshold values to be evaluated.
#' @param include.mean logical values indicating whether constant terms are included.
#' @param score the choice of criterion used in selection threshold, namely (AIC, det(RSS)).
#' @references
#' Li, D., and Tong. H. (2016) Nested sub-sample search algorithm for estimation of threshold models. \emph{Statisitca Sinica}, 1543-1554.
#' @return mTAR returns a list with the following components:
#' \item{data}{the data matrix, y.}
#' \item{beta}{a (\code{p*k+1})-by-(\code{2k}) matrices. The first \code{k} columns show the estimation results in regime 1, and the second \code{k} columns show these in regime 2.}
#' \item{arorder}{AR orders of regimes 1 and 2.}
#' \item{sigma}{estimated innovational covariance matrices of regimes 1 and 2.}
#' \item{residuals}{estimated innovations.}
#' \item{nobs}{numbers of observations in regimes 1 and 2.}
#' \item{model1, model2}{estimated models of regimes 1 and 2.}
#' \item{thr}{threshold value.}
#' \item{delay}{two elements (\code{i},\code{d}) with "\code{i}" being the component and "\code{d}" the delay for threshold variable.}
#' \item{thrV}{vector of threshold variable.}
#' \item{D}{a set of positive threshold values.}
#' \item{RSS}{residual sum of squares.}
#' \item{information}{overall information criteria.}
#' \item{cnst}{logical values indicating whether the constant terms are included in regimes 1 and 2.}
#' \item{sresi}{standardized residuals.}
#' @examples
#' phi1=matrix(c(0.5,0.7,0.3,0.2),2,2)
#' phi2=matrix(c(0.4,0.6,0.5,-0.5),2,2)
#' sigma1=matrix(c(1,0,0,1),2,2)
#' sigma2=matrix(c(1,0,0,1),2,2)
#' c1=c(0,0)
#' c2=c(0,0)
#' delay=c(1,1)
#' Trim=c(0.2,0.8)
#' include.mean=TRUE
#' y=mTAR.sim(1000,0,phi1,phi2,sigma1,sigma2,c1,c2,delay,ini=500)
#' est=mTAR(y$series,1,1,0,y$series,delay,Trim,300,include.mean,"AIC")
#' est2=mTAR(y$series,1,1,NULL,y$series,delay,Trim,300,include.mean,"AIC")
#' @export
"mTAR" <- function(y,p1,p2,thr=NULL,thrV=NULL,delay=c(1,1),Trim=c(0.1,0.9),k0=300,include.mean=TRUE,score="AIC"){
  if(!is.matrix(y))y <- as.matrix(y)
  p = max(p1,p2)
  if(p < 1)p=1
  d = delay[2]
  ist = max(p,d)+1
  nT <- nrow(y)
  ky <- ncol(y)
  nobe <- nT-ist+1
  ### Built in checking
  nT1 <- nT
  if(!is.null(thrV))nT1 <- length(thrV)
  if(nT != nT1){
    cat("Input error in thrV. Reset to a SETAR model","\n")
    thrV <- y[,delay[1]]
  }
  ##
  if(length(thrV) < nT){
    thrV <- y[(ist-d):(nT-d),delay[1]]
  }else{
    thrV <- thrV[(ist-d):(nT-d)]
  }
  beta <- matrix(0,(p*ky+1),ky*2)
  ### set up regression framework
  yt <- y[ist:nT,]
  x1 <- NULL
  for (i in 1:p){
    x1=cbind(x1,y[(ist-i):(nT-i),])
  }
  X <- x1
  if(include.mean){X=cbind(rep(1,nobe),X)}
  ### search for threshold if needed
  if(is.null(thr)){
    D <- NeSS(yt,x1,thrV=thrV,Trim=Trim,k0=k0,include.mean=include.mean,score=score)
    ##cat("Threshold candidates: ",D,"\n")
    ll = length(D)
    RSS=NULL
    ic = 2
    #### If not AIC, generalized variance is used.
    if(score=="AIC")ic <- 1
    for (i in 1:ll){
      thr = D[i]
      idx = c(1:nobe)[thrV <= thr]
      m1a <- MlmNeSS(yt,X,subset=idx,SD=FALSE,include.mean=include.mean)
      m1b <- MlmNeSS(yt,X,subset=-idx,SD=FALSE,include.mean=include.mean)
      RSS=c(RSS,sum(m1a$score[ic],m1b$score[ic]))
    }
    #cat("RSS: ",RSS,"\n")
    j=which.min(RSS)
    thr <- D[j]
    cat("Estimated Threshold: ",thr,"\n")
  }else{
    cat("Threshold: ",thr,"\n")
    RSS <- NULL
  }
  #### The above ends the search of threshold ###

  ### Final estimation results
  ###cat("Coefficients are in vec(beta), where t(beta)=[phi0,phi1,...]","\n")
  resi = matrix(0,nobe,ncol(yt))
  sresi <- matrix(0,nobe,ncol(yt))
  idx = c(1:nobe)[thrV <= thr]
  n1 <- length(idx)
  n2 <- nrow(yt)-n1
  k <- ncol(yt)
  X = x1
  if(include.mean){X=cbind(rep(1,nobe),X)}
  m1a <- MlmNeSS(yt,X,subset=idx,SD=TRUE,include.mean=include.mean)
  np <- ncol(X)
  resi[idx,] <- m1a$residuals
  infc = m1a$information
  beta[1:nrow(m1a$beta),1:ky] <- m1a$beta
  cat("Regime 1 with sample size: ",n1,"\n")
  coef1 <- m1a$coef
  coef1 <- cbind(coef1,coef1[,1]/coef1[,2])
  colnames(coef1) <- c("est","s.e.","t-ratio")
  for (j in 1:k){
    jdx=(j-1)*np
    cat("Model for the ",j,"-th component (including constant, if any): ","\n")
    print(round(coef1[(jdx+1):(jdx+np),],4))
  }
  cat("sigma: ", "\n")
  print(m1a$sigma)
  cat("Information(aic,bix,hq): ",infc,"\n")
  m1 <- eigen(m1a$sigma)
  P <- m1$vectors
  Di <- diag(1/sqrt(m1$values))
  siv <- P%*%Di%*%t(P)
  sresi[idx,] <- m1a$residuals%*%siv
  ###
  m1b <- MlmNeSS(yt,X,subset=-idx,SD=TRUE,include.mean=include.mean)
  resi[-idx,] <- m1b$residuals
  beta[1:nrow(m1b$beta),(ky+1):(2*ky)] <- m1b$beta
  cat(" ","\n")
  cat("Regime 2 with sample size: ",n2,"\n")
  coef2 <- m1b$coef
  coef2 <- cbind(coef2,coef2[,1]/coef2[,2])
  colnames(coef2) <- c("est","s.e.","t-ratio")
  for (j in 1:k){
    jdx=(j-1)*np
    cat("Model for the ",j,"-th component (including constant, if any): ","\n")
    print(round(coef2[(jdx+1):(jdx+np),],4))
  }
  cat("sigma: ","\n")
  print(m1b$sigma)
  cat("Information(aic,bic,hq): ",m1b$information,"\n")
  m2 <- eigen(m1b$sigma)
  P <- m2$vectors
  Di <- diag(1/sqrt(m2$values))
  siv <- P%*%Di%*%t(P)
  sresi[-idx,] <- m1b$residuals%*%siv
  #
  fitted.values <- yt-resi
  cat(" ","\n")
  cat("Overall pooled estimate of sigma: ","\n")
  sigma = (n1*m1a$sigma+n2*m1b$sigma)/(n1+n2)
  print(sigma)
  infc = m1a$information+m1b$information
  cat("Overall information criteria(aic,bic,hq): ",infc,"\n")
  Sigma <- cbind(m1a$sigma,m1b$sigma)

  mTAR <- list(data=y,beta=beta,arorder=c(p1,p2),sigma=Sigma,residuals = resi,nobs=c(n1,n2),
               model1 = m1a, model2 = m1b,thr=thr, delay=delay, thrV=thrV, D=D, RSS=RSS, information=infc,
               cnst=rep(include.mean,2),sresi=sresi)
}



#' Estimation of Multivariate TAR Models
#'
#' Estimation of multivariate TAR models with given thresholds. It can handle multiple regimes.
#' @param y vector time series.
#' @param arorder AR order of each regime. The number of regime is length of arorder.
#' @param thr threshold value(s). There are k-1 threshold for a k-regime model.
#' @param delay two elements (i,d) with "i" being the component and "d" the delay for threshold variable.
#' @param thrV external threshold variable if any. If thrV is not null, it must have the same number of observations as y-series.
#' @param include.mean logical values indicating whether constant terms are included. Default is TRUE for all.
#' @param output a logical value indicating four output. Default is TRUE.
#' @return mTAR.est returns a list with the following components:
#' \item{data}{the data matrix, \code{y}.}
#' \item{k}{the dimension of \code{y}.}
#' \item{arorder}{AR orders of regimes 1 and 2.}
#' \item{beta}{a (\code{p*k+1})-by-(\code{2k}) matrices. The first \code{k} columns show the estimation results in regime 1, and the second \code{k} columns show these in regime 2.}
#' \item{sigma}{estimated innovational covariance matrices of regimes 1 and 2.}
#' \item{thr}{threshold value.}
#' \item{residuals}{estimated innovations.}
#' \item{sresi}{standardized residuals.}
#' \item{nobs}{numbers of observations in different regimes.}
#' \item{cnst}{logical values indicating whether the constant terms are included in different regimes.}
#' \item{AIC}{AIC value.}
#' \item{delay}{two elements (\code{i,d}) with "\code{i}" being the component and "\code{d}" the delay for threshold variable.}
#' \item{thrV}{values of threshold variable.}
#' @examples
#' phi1=matrix(c(0.5,0.7,0.3,0.2),2,2)
#' phi2=matrix(c(0.4,0.6,0.5,-0.5),2,2)
#' sigma1=matrix(c(1,0,0,1),2,2)
#' sigma2=matrix(c(1,0,0,1),2,2)
#' c1=c(0,0)
#' c2=c(0,0)
#' delay=c(1,1)
#' y=mTAR.sim(100,0,phi1,phi2,sigma1,sigma2,c1,c2,delay,ini=500)
#' est=mTAR.est(y$series,c(1,1),0,delay)
#' @export
"mTAR.est" <- function(y,arorder=c(1,1),thr=c(0),delay=c(1,1),thrV=NULL,include.mean=c(TRUE,TRUE),output=TRUE){
  if(!is.matrix(y)) y <- as.matrix(y)
  p <- max(arorder)
  if(p < 1)p <-1
  k <- length(arorder)
  d <- delay[2]
  ist <- max(p,d)+1
  nT <- nrow(y)
  ky <- ncol(y)
  Y <- y[ist:nT,]
  X <- NULL
  for (i in 1:p){
    X <- cbind(X,y[(ist-i):(nT-i),])
  }
  X <- as.matrix(X)
  nobe <- nT-ist+1
  if(length(thrV) < nT){
    thrV <- y[(ist-d):(nT-d),delay[1]]
    if(output) cat("Estimation of a multivariate SETAR model: ","\n")
  }else{thrV <- thrV[ist:nT]}
  ### make sure that the thresholds are in increasing order!!!
  thr <- sort(thr)
  k1 <- length(thr)
  ### house keeping
  beta <- NULL
  sigma <- NULL
  nobs <- NULL
  resi <- matrix(0,nobe,ky)
  sresi <- matrix(0,nobe,ky)
  npar <- NULL
  AIC <- 99999
  #
  if(k1 > (k-1)){
    cat("Only the first: ",k1," thresholds are used.","\n")
    THr <- min(thrV)-1
    thr <- thr[1:(k-1)]
    THr <- c(THr,thr,(max(thrV)+1))
  }
  if(k1 < (k-1)){cat("Too few threshold values are given!!!","\n")}
  if(length(include.mean) != k)include.mean=rep("TRUE",k)
  if(k1 == (k-1)){
    THr <- c((min(thrV[1:nobe])-1),thr,c(max(thrV[1:nobe])+1))
    if(output) cat("Thresholds used: ",thr,"\n")
    for (j in 1:k){
      idx <- c(1:nobe)[thrV <= THr[j]]
      jdx <- c(1:nobe)[thrV > THr[j+1]]
      jrm <- c(idx,jdx)
      n1 <- nobe-length(jrm)
      nobs <- c(nobs,n1)
      kdx <- c(1:nobe)[-jrm]
      x1 <- X
      if(include.mean[j]){x1 <- cbind(rep(1,nobe),x1)}
      np <- ncol(x1)*ky
      npar <- c(npar,np)
      beta <- matrix(0,(p*ky+1),k*ky)
      sig <- diag(rep(1,ky))
      if(n1 <= ncol(x1)){
        cat("Regime ",j," does not have sufficient observations!","\n")
        sigma <- cbind(sigma,sig)
      }
      else{
        mm <- MlmNeSS(Y,x1,subset=kdx,SD=TRUE,include.mean=include.mean[j])
        jst = (j-1)*k
        beta1 <- mm$beta
        beta[1:nrow(beta1),(jst+1):(jst+k)] <- beta1
        Sigma <- mm$sigma
        if(output){
          cat("Estimation of regime: ",j," with sample size: ",n1,"\n")
          cat("Coefficient estimates:t(phi0,phi1,...) ","\n")
          colnames(mm$coef) <- c("est","se")
          print(mm$coef)
          cat("Sigma: ","\n")
          print(Sigma)
        }
        sigma <- cbind(sigma,Sigma)
        resi[kdx,] <- mm$residuals
        m1 <- eigen(Sigma)
        P <- m1$vectors
        Di <- diag(1/sqrt(m1$values))
        Sighinv <- P%*%Di%*%t(P)
        sresi[kdx,] <- mm$residuals %*% Sighinv
      }
    }
    AIC <- 0
    for (j in 1:k){
      sig <- sigma[,((j-1)*k+1):(j*k)]
      d1 <- det(sig)
      n1 <- nobs[j]
      np <- npar[j]
      if(!is.null(sig)){s2 <- sig^2*(n1-np)/n1
      AIC <- AIC + log(d1)*n1 + 2*np
      }
    }
    if(output) cat("Overal AIC: ",AIC,"\n")
  }
  mTAR.est <- list(data=y,k = k, arorder=arorder,beta=beta,sigma=sigma,thr=thr,residuals=resi,
                   sresi=sresi,nobs=nobs,cnst=include.mean, AIC=AIC, delay=delay,thrV=thrV)
}



#' Refine A Fitted 2-Regime Multivariate TAR Model
#'
#' Refine a fitted 2-regime multivariate TAR model using "thres" as threshold for t-ratios.
#' @param m1 a fitted mTAR object.
#' @param thres threshold value.
#' @return ref.mTAR returns a list with following components:
#' \item{data}{data matrix, \code{y}.}
#' \item{arorder}{AR orders of regimes 1 and 2.}
#' \item{sigma}{estimated innovational covariance matrices of regimes 1 and 2.}
#' \item{beta}{a (\code{p*k+1})-by-(\code{2k}) matrices. The first \code{k} columns show the estimation results in regime 1, and the second \code{k} columns shows these in regime 2.}
#' \item{residuals}{estimated innovations.}
#' \item{sresi}{standard residuals.}
#' \item{criteria}{overall information criteria.}
#' @examples
#' phi1=matrix(c(0.5,0.7,0.3,0.2),2,2)
#' phi2=matrix(c(0.4,0.6,0.5,-0.5),2,2)
#' sigma1=matrix(c(1,0,0,1),2,2)
#' sigma2=matrix(c(1,0,0,1),2,2)
#' c1=c(0,0)
#' c2=c(0,0)
#' delay=c(1,1)
#' y=mTAR.sim(100,0,phi1,phi2,sigma1,sigma2,c1,c2,delay,ini=500)
#' est=mTAR.est(y$series,c(1,1),0,delay)
#' ref.mTAR(est,0)
#' @export
"ref.mTAR" <- function(m1,thres=1.0){
  y <- m1$data
  arorder <- m1$arorder
  sigma <- m1$sigma
  nobs <- m1$nobs
  thr <- m1$thr
  delay <- m1$delay
  thrV <- m1$thrV
  include.mean <- m1$cnst
  mA <- m1$model1
  mB <- m1$model2
  vName <- colnames(y)
  if(is.null(vName)){
    vName <- paste("V",1:ncol(y))
    colnames(y) <- vName
  }
  ###
  p = max(arorder)
  d = delay[2]
  ist = max(p,d)+1
  nT <- nrow(y)
  ky <- ncol(y)
  nobe <- nT-ist+1
  ### set up regression framework
  yt <- y[ist:nT,]
  x1 <- NULL
  xname <- NULL
  vNameL <- paste(vName,"L",sep="")
  for (i in 1:p){
    x1=cbind(x1,y[(ist-i):(nT-i),])
    xname <- c(xname,paste(vNameL,i,sep=""))
  }
  colnames(x1) <- xname
  X <- x1
  if(include.mean[1]){X <- cbind(rep(1,nobe),X)
  xname <- c("cnst",xname)
  }
  colnames(X) <- xname
  ####
  ### Final estimation results
  ###cat("Coefficients are in vec(beta), where t(beta)=[phi0,phi1,...]","\n")
  npar <- NULL
  np <- ncol(X)
  beta <- matrix(0,np,ky*2)
  resi = matrix(0,nobe,ky)
  sresi <- matrix(0,nobe,ky)
  idx = c(1:nobe)[thrV <= thr]
  n1 <- length(idx)
  n2 <- nrow(yt)-n1
  ### Regime 1, component-by-component
  cat("Regime 1 with sample size: ",n1,"\n")
  coef <- mA$coef
  RES <- NULL
  for (ij in 1:ky){
    tra <- rep(0,np)
    i1 <- (ij-1)*np
    jj <- c(1:np)[coef[(i1+1):(i1+np),2] > 0]
    tra[jj] <- coef[(i1+jj),1]/coef[(i1+jj),2]
    ###   cat("tra: ",tra,"\n")
    kdx <- c(1:np)[abs(tra) >= thres]
    npar <- c(npar,length(kdx))
    ##
    cat("Significant variables: ",kdx,"\n")
    if(length(kdx) > 0){
      xreg <- as.matrix(X[idx,kdx])
      colnames(xreg) <- colnames(X)[kdx]
      m1a <- lm(yt[idx,ij]~-1+xreg)
      m1aS <- summary(m1a)
      beta[kdx,ij] <- m1aS$coefficients[,1]
      resi[idx,ij] <- m1a$residuals
      cat("Component: ",ij,"\n")
      print(m1aS$coefficients)
    }else{
      resi[idx,ij] <- yt[idx,ij]
    }
    RES <- cbind(RES,c(resi[idx,ij]))
  }
  ###
  sigma1 <- crossprod(RES,RES)/n1
  cat("sigma: ", "\n")
  print(sigma1)
  # cat("Information(aic,bix,hq): ",infc,"\n")
  m1 <- eigen(sigma1)
  P <- m1$vectors
  Di <- diag(1/sqrt(m1$values))
  siv <- P%*%Di%*%t(P)
  sresi[idx,] <- RES%*%siv
  ###  Regime 2 ####
  cat("Regime 2 with sample size: ",n2,"\n")
  coef <- mB$coef
  RES <- NULL
  for (ij in 1:ky){
    tra <- rep(0,np)
    i1 <- (ij-1)*np
    jj <- c(1:np)[coef[(i1+1):(i1+np),2] > 0]
    tra[jj] <- coef[(i1+jj),1]/coef[(i1+jj),2]
    ###   cat("tra: ",tra,"\n")
    kdx <- c(1:np)[abs(tra) >= thres]
    ##
    cat("Significant variables: ",kdx,"\n")
    npar <- c(npar,length(kdx))
    if(length(kdx) > 0){
      xreg <- as.matrix(X[-idx,kdx])
      colnames(xreg) <- colnames(X)[kdx]
      m1a <- lm(yt[-idx,ij]~-1+xreg)
      m1aS <- summary(m1a)
      beta[kdx,ij+ky] <- m1aS$coefficients[,1]
      resi[-idx,ij] <- m1a$residuals
      cat("Component: ",ij,"\n")
      print(m1aS$coefficients)
    }else{
      resi[-idx,ij] <- yt[-idx,ij]
    }
    RES <- cbind(RES,c(resi[-idx,ij]))
  }
  ###
  sigma2 <- crossprod(RES,RES)/n2
  cat("sigma: ", "\n")
  print(sigma2)
  # cat("Information(aic,bix,hq): ",infc,"\n")
  m1 <- eigen(sigma2)
  P <- m1$vectors
  Di <- diag(1/sqrt(m1$values))
  siv <- P%*%Di%*%t(P)
  sresi[-idx,] <- RES%*%siv
  #
  fitted.values <- yt-resi
  cat(" ","\n")
  cat("Overall pooled estimate of sigma: ","\n")
  sigma = (n1*sigma1+n2*sigma2)/(n1+n2)
  print(sigma)
  d1 <- log(det(sigma))
  nnp <- sum(npar)
  TT <- (n1+n2)
  aic <- TT*d1+2*nnp
  bic <- TT*d1+log(TT)*nnp
  hq <- TT*d1+2*log(log(TT))*nnp
  cri <- c(aic,bic,hq)
  ### infc = m1a$information+m1b$information
  cat("Overall information criteria(aic,bic,hq): ",cri,"\n")
  Sigma <- cbind(sigma1,sigma2)
  ref.mTAR <- list(data=y,arorder=arorder,sigma=Sigma,beta=beta,residuals=resi,sresi=sresi,criteria=cri)
}


#' Prediction of A Fitted Multivariate TAR Model
#'
#' Prediction of a fitted multivariate TAR model.
#' @param model multivariate TAR model.
#' @param orig forecast origin.
#' @param h forecast horizon.
#' @param iterations number of iterations.
#' @param ci confidence level.
#' @param output a logical value for output.
#' @return mTAR.pred returns a list with components:
#' \item{model}{the multivariate TAR model.}
#' \item{pred}{prediction.}
#' \item{Ysim}{fitted \code{y}.}
#' @examples
#' phi1=matrix(c(0.5,0.7,0.3,0.2),2,2)
#' phi2=matrix(c(0.4,0.6,0.5,-0.5),2,2)
#' sigma1=matrix(c(1,0,0,1),2,2)
#' sigma2=matrix(c(1,0,0,1),2,2)
#' c1=c(0,0)
#' c2=c(0,0)
#' delay=c(1,1)
#' y=mTAR.sim(100,0,phi1,phi2,sigma1,sigma2,c1,c2,delay,ini=500)
#' est=mTAR.est(y$series,c(1,1),0,delay)
#' pred=mTAR.pred(est,100,1,300,0.90,TRUE)
#' @export
"mTAR.pred" <- function(model,orig,h=1,iterations=3000,ci=0.95,output=TRUE){
  ## ci: probability of pointwise confidence interval coefficient
  ### only works for self-exciting mTAR models.
  ###
  y <- model$data
  arorder <- model$arorder
  beta <- model$beta
  sigma <- model$sigma
  thr <- model$thr
  include.mean <- model$cnst
  delay <- model$delay
  p <- max(arorder)
  k <- length(arorder)
  d <- delay[2]
  nT <- nrow(y)
  ky <- ncol(y)
  if(orig < 1)orig <- nT
  if(orig > nT) orig <- nT
  if(h < 1) h <- 1
  ### compute the predictions. Use simulation if the threshold is predicted
  Sigh <- NULL
  for (j in 1:k){
    sig <- sigma[,((j-1)*ky+1):(j*ky)]
    m1 <- eigen(sig)
    P <- m1$vectors
    Di <- diag(sqrt(m1$values))
    sigh <- P%*%Di%*%t(P)
    Sigh <- cbind(Sigh,sigh)
  }
  Ysim <- array(0,dim=c(h,ky,iterations))
  for (it in 1:iterations){
    yp <- y[1:orig,]
    et <- matrix(rnorm(h*ky),h,ky)
    for (ii in 1:h){
      t <- orig+ii
      thd <- yp[(t-d),delay[1]]
      JJ <- 1
      for (j in 1:(k-1)){
        if(thd > thr[j]){JJ <- j+1}
      }
      Jst <- (JJ-1)*ky
      at <- matrix(et[ii,],1,ky)%*%Sigh[,(Jst+1):(Jst+ky)]
      x <- NULL
      if(include.mean[JJ]) x <- 1
      pJ <- arorder[JJ]
      phi <- beta[,(Jst+1):(Jst+ky)]
      for (i in 1:pJ){
        x <- c(x,yp[(t-i),])
      }
      yhat <- matrix(x,1,length(x))%*%phi[1:length(x),]
      yhat <- yhat+at
      yp <- rbind(yp,yhat)
      Ysim[ii,,it] <- yhat
    }
  }
  ### summary
  pred <- NULL
  upp <- NULL
  low <- NULL
  pr <- (1-ci)/2
  pro <- c(pr, 1-pr)
  for (ii in 1:h){
    fst <- NULL
    lowb <- NULL
    uppb <- NULL
    for (j in 1:ky){
      ave <- mean(Ysim[ii,j,])
      quti <- quantile(Ysim[ii,j,],prob=pro)
      fst <- c(fst,ave)
      lowb <- c(lowb,quti[1])
      uppb <- c(uppb,quti[2])
    }
    pred <- rbind(pred,fst)
    low <- rbind(low,lowb)
    upp <- rbind(upp,uppb)
  }
  if(output){
    colnames(pred) <- colnames(y)
    cat("Forecast origin: ",orig,"\n")
    cat("Predictions: 1-step to ",h,"-step","\n")
    print(pred)
    cat("Lower bounds of ",ci*100," % confident intervals","\n")
    print(low)
    cat("Upper bounds: ","\n")
    print(upp)
  }
  mTAR.pred <- list(data=y,pred = pred,Ysim=Ysim)
}


"MlmNeSS" <- function(y,z,subset=NULL,SD=FALSE,include.mean=TRUE){
  # z: design matrix, including 1 as its first column if constant is needed.
  # y: dependent variables
  # subset: locators for the rows used in estimation
  ## Model is y = z%*%beta+error
  #### include.mean is ONLY used in counting parameters. z-matrix should include column of 1 if constant
  #### is needed.
  #### SD: switch for computing standard errors of the estimates
  ####
  zc = as.matrix(z)
  yc = as.matrix(y)
  if(!is.null(subset)){
    zc <- zc[subset,]
    yc <- yc[subset,]
  }
  n=nrow(zc)
  p=ncol(zc)
  k <- ncol(y)
  coef=NULL
  infc = rep(9999,3)
  if(n <= p){
    cat("too few data points in a regime: ","\n")
    beta = NULL
    res = yc
    sig = NULL
  }
  else{
    ztz=t(zc)%*%zc/n
    zty=t(zc)%*%yc/n
    ztzinv=solve(ztz)
    beta=ztzinv%*%zty
    res=yc-zc%*%beta
    ### Sum of squares of residuals
    sig=t(res)%*%res/n
    npar=ncol(zc)
    if(include.mean)npar=npar-1
    score1=n*log(det(sig))+2*npar
    score2=det(sig*n)
    score=c(score1,score2)
    ###
    if(SD){
      sigls <- sig*(n/(n-p))
      s <- kronecker(sigls,ztzinv)
      se <- sqrt(diag(s)/n)
      coef <- cbind(c(beta),se)
      #### MLE estimate of sigma
      d1 = log(det(sig))
      npar=nrow(coef)
      if(include.mean)npar=npar-1
      aic <- n*d1+2*npar
      bic <- n*d1 + log(n)*npar
      hq <- n*d1 + 2*log(log(n))*npar
      infc=c(aic,bic,hq)
    }
  }
  #
  MlmNeSS <- list(beta=beta,residuals=res,sigma=sig,coef=coef,information=infc,score=score)
}

#' Generate Univariate 2-regime Markov Switching Models
#'
#' Generate univariate 2-regime Markov switching models.
#' @param nob number of observations.
#' @param order AR order for each regime.
#' @param phi1,phi2 AR coefficients.
#' @param epsilon transition probabilities (switching out of regime 1 and 2).
#' @param sigma standard errors for each regime.
#' @param cnst constant term for each regime.
#' @param ini burn-in period.
#' @return MSM.sim returns a list with components:
#' \item{series}{a time series following SETAR model.}
#' \item{at}{innovation of the time series.}
#' \item{state}{states for the time series.}
#' \item{epsilon}{transition probabilities (switching out of regime 1 and 2).}
#' \item{sigma}{standard error for each regime.}
#' \item{cnst}{constant terms.}
#' \item{order}{AR-order for each regime.}
#' \item{phi1, phi2}{the AR coefficients for two regimes.}
#' @examples
#' y=MSM.sim(100,c(1,1),0.7,-0.5,c(0.5,0.6),c(1,1),c(0,0),500)
#' @export
"MSM.sim" <- function(nob,order=c(1,1),phi1=NULL,phi2=NULL,epsilon=c(0.1,0.1),sigma=c(1,1),cnst=c(0,0),ini=500){
  nT <- nob+ini
  et <- rnorm(nT)
  p <- max(order)
  if(p < 1) p <- 1
  if(order[1] < 0)order[1] <- 1
  if(order[2] < 0)order[2] <- 1
  if(sigma[1] < 0)sigma[1] <- 1
  if(sigma[2] < 0)sigma[2] <- 1
  ist <- p+1
  y <- et[1:p]
  at <- et
  state <- rbinom(p,1,0.5)+1
  for (i in ist:nT){
    p1 <- epsilon[1]
    sti <- state[i-1]
    if(sti == 2)p1 <- epsilon[2]
    sw <- rbinom(1,1,p1)
    st <- sti
    if(sti == 1){
      if(sw == 1)st <- 2
    }
    else{
      if(sw == 1)st <- 1
    }
    if(st == 1){
      tmp <- cnst[1]
      at[i] <- et[i]*sigma[1]
      if(order[1] > 0){
        for (j in 1:order[1]){
          tmp = tmp + phi1[j]*y[i-j]
        }
      }
      y[i] <- tmp+at[i]
    }
    else{ tmp <- cnst[2]
    at[i] <- et[i]*sigma[2]
    if(order[2] > 0){
      for (j in 1:order[2]){
        tmp <- tmp + phi2[j]*y[i-j]
      }
    }
    y[i] <- tmp + at[i]
    }
    state <- c(state,st)
  }
  MSM.sim <- list(series=y[(ini+1):nT],at = at[(ini+1):nT], state=state[(ini+1):nT], epsilon=epsilon,
                  sigma=sigma,cnst=cnst,order=order,phi1=phi1,phi2=phi2)
}






#' Threshold Nonlinearity Test
#'
#' Threshold nonlinearity test.
#' @references
#' Tsay, R. (1989) Testing and Modeling Threshold Autoregressive Processes. \emph{Journal of the American Statistical Associations} \strong{84}(405), 231-240.
#'
#' @param y a time series.
#' @param p AR order.
#' @param d delay for the threshold variable.
#' @param thrV threshold variable.
#' @param ini initial number of data to start RLS estimation.
#' @param include.mean a logical value for including constant terms.
#' @return \code{thr.test} returns a list with components:
#' \item{F-ratio}{F statistic.}
#' \item{df}{the numerator and denominator degrees of freedom.}
#' \item{ini}{initial number of data to start RLS estimation.}
#' @examples
#' phi=t(matrix(c(-0.3, 0.5,0.6,-0.3),2,2))
#' y=uTAR.sim(nob=2000, arorder=c(2,2), phi=phi, d=2, thr=0.2, cnst=c(1,-1),sigma=c(1, 1))
#' thr.test(y$series,p=2,d=2,ini=40,include.mean=TRUE)
#' @export
"thr.test" <- function(y,p=1,d=1,thrV=NULL,ini=40,include.mean=T){
  if(is.matrix(y))y <- y[,1]
  nT <- length(y)
  if(p < 1)p <-1
  if(d < 1) d <-1
  ist <- max(p,d)+1
  nobe <- nT-ist+1
  if(length(thrV) < nobe){
    cat("SETAR model is entertained","\n")
    thrV <- y[(ist-d):(nT-d)]
  }else{thrV <- thrV[1:nobe]}
  #
  Y <- y[ist:nT]
  X <- NULL
  for (i in 1:p){
    X <- cbind(X,y[(ist-i):(nT-i)])
  }
  if(include.mean){X <- cbind(rep(1,nobe),X)}
  X <- as.matrix(X)
  kx <- ncol(X)
  m1 <- sort(thrV,index.return=TRUE)
  idx <- m1$ix
  Y <- Y[idx]
  if(kx == 1){X <- X[idx]
  }else{ X <- X[idx,]}
  while(kx > ini){ini <- ini+10}
  if(kx == 1){
    XpX = sum(X[1:ini]^2)
    XpY <- sum(X[1:ini]*Y[1:ini])
    betak <- XpY/XpX
    Pk <- 1/XpX
    resi <- Y[1:ini] - X[1:ini]*betak
  }
  else{ XpX <- t(X[1:ini,])%*%X[1:ini,]
  XpY <- t(X[1:ini,])%*%as.matrix(Y[1:ini])
  Pk <- solve(XpX)
  betak <- Pk%*%XpY
  resi <- Y[1:ini]-X[1:ini,]%*%betak
  }
  ### RLS to obtain standardized residuals
  sresi <- resi
  if(ini < nobe){
    for (i in (ini+1):nobe){
      if(kx == 1){xk <- X[i]
      er <- xk*betak-Y[i]
      deno <- 1+xk^2/Pk
      betak <- betak -Pk*er/deno
      tmp <- Pk*xk
      Pk <- tmp*tmp/deno
      sresi <- c(sresi,er/sqrt(deno))
      }
      else{
        xk <- matrix(X[i,],kx,1)
        tmp <- Pk%*%xk
        er <- t(xk)%*%betak - Y[i]
        deno <- 1  + t(tmp)%*%xk
        betak <- betak - tmp*c(er)/c(deno)
        Pk <- Pk -tmp%*%t(tmp)/c(deno)
        sresi <- c(sresi,c(er)/sqrt(c(deno)))
      }
    }
    ###cat("final esimates: ",betak,"\n")
    ### testing
    Ys <- sresi[(ini+1):nobe]
    if(kx == 1){x <- data.frame(X[(ini+1):nobe])
    }else{ x <- data.frame(X[(ini+1):nobe,])}
    m2 <- lm(Ys~.,data=x)
    Num <- sum(Ys^2)-sum(m2$residuals^2)
    Deno <- sum(m2$residuals^2)
    h <- max(1,p+1-d)
    df1 <- p+1
    df2 <- nT-d-ini-p-h
    Fratio <- (Num/df1)/(Deno/df2)

    cat("Threshold nonlinearity test for (p,d): ",c(p,d),"\n")
    pv <- pf(Fratio,df1,df2,lower.tail=FALSE)
    cat("F-ratio and p-value: ",c(Fratio,pv),"\n")
  }
  else{cat("Insufficient sample size to perform the test!","\n")}

  thr.test <- list(F.ratio = Fratio,df=c(df1,df2),ini=ini)
}



#' Tsay Test for Nonlinearity
#'
#' Perform Tsay (1986) nonlinearity test.
#' @references
#' Tsay, R. (1986) Nonlinearity tests for time series. \emph{Biometrika} \strong{73}(2), 461-466.
#' @param y time series.
#' @param p AR order.
#' @return The function outputs the F statistic, p value, and the degrees of freedom. The null hypothesis is there is no nonlinearity.
#' @examples
#' phi=t(matrix(c(-0.3, 0.5,0.6,-0.3),2,2))
#' y=uTAR.sim(nob=2000, arorder=c(2,2), phi=phi, d=2, thr=0.2, cnst=c(1,-1),sigma=c(1, 1))
#' Tsay(y$series,2)
#' @export
"Tsay" <- function(y,p=1){
  if(is.matrix(y))y=y[,1]
  nT <- length(y)
  if(p < 1)p <- 1
  ist <- p+1
  Y <- y[ist:nT]
  ym <- scale(y,center=TRUE,scale=FALSE)
  X <- rep(1,nT-p)
  for (i in 1:p){
    X <- cbind(X,ym[(ist-i):(nT-i)])
  }
  colnames(X) <- c("cnst",paste("lag",1:p))
  m1 <- lm(Y~-1+.,data=data.frame(X))
  resi <- m1$residuals
  XX <- NULL
  for (i in 1:p){
    for (j in 1:i){
      XX <- cbind(XX,ym[(ist-i):(nT-i)]*ym[(ist-j):(nT-j)])
    }
  }
  XR <- NULL
  for (i in 1:ncol(XX)){
    mm <- lm(XX[,i]~.,data=data.frame(X))
    XR <- cbind(XR,mm$residuals)
  }
  colnames(XR) <- paste("crossp",1:ncol(XX))
  m2 <- lm(resi~-1+.,data=data.frame(XR))
  resi2 <- m2$residuals
  SS0 <- sum(resi^2)
  SS1 <- sum(resi2^2)
  df1 <- ncol(XX)
  df2 <- nT-p-df1-1
  Num <- (SS0-SS1)/df1
  Deno <- SS1/df2
  F <- Num/Deno
  pv <- 1-pf(F,df1,df2)
  cat("Non-linearity test & its p-value: ",c(F,pv),"\n")
  cat("Degrees of freedom of the test: ",c(df1,df2),"\n")
}


#' Backtest
#'
#' Backtest for an ARIMA time series model.
#' @param m1 an ARIMA time series model object.
#' @param rt the time series.
#' @param orig forecast origin.
#' @param h forecast horizon.
#' @param xre the independent variables.
#' @param fixed parameter constraint.
#' @param include.mean a logical value for constant term of the model. Default is TRUE.
#' @return The function returns a list with following components:
#' \item{orig}{the starting forecast origin.}
#' \item{err}{observed value minus fitted value.}
#' \item{rmse}{RMSE of out-of-sample forecasts.}
#' \item{mabso}{mean absolute error of out-of-sample forecasts.}
#' \item{bias}{bias of out-of-sample forecasts.}
#' @examples
#' data=arima.sim(n=100,list(ar=c(0.5,0.3)))
#' model=arima(data,order=c(2,0,0))
#' backtest(model,data,orig=70,h=1)
#' @export
"backtest" <- function(m1,rt,orig,h,xre=NULL,fixed=NULL,include.mean=TRUE){
  # m1: is a time-series model object
  # orig: is the starting forecast origin
  # rt: the time series
  # xre: the independent variables
  # h: forecast horizon
  # fixed: parameter constraint
  # inc.mean: flag for constant term of the model.
  #
  regor=c(m1$arma[1],m1$arma[6],m1$arma[2])
  seaor=list(order=c(m1$arma[3],m1$arma[7],m1$arma[4]),period=m1$arma[5])
  nT=length(rt)
  if(orig > nT)orig=nT
  if(h < 1) h=1
  rmse=rep(0,h)
  mabso=rep(0,h)
  bias <- rep(0,h)
  nori=nT-orig
  err=matrix(0,nori,h)
  jlast=nT-1
  for (n in orig:jlast){
    jcnt=n-orig+1
    x=rt[1:n]
    if (is.null(xre))
      pretor=NULL else pretor=xre[1:n,]
    mm=arima(x,order=regor,seasonal=seaor,xreg=pretor,fixed=fixed,include.mean=include.mean)
    if (is.null(xre)){nx=NULL}
    else {nx=matrix(xre[(n+1):(n+h),],h,ncol(xre))}
    fore=predict(mm,h,newxreg=nx)
    kk=min(nT,(n+h))
    # nof is the effective number of forecasts at the forecast origin n.
    nof=kk-n
    pred=fore$pred[1:nof]
    obsd=rt[(n+1):kk]
    err[jcnt,1:nof]=obsd-pred
  }
  #
  for (i in 1:h){
    iend=nori-i+1
    tmp=err[1:iend,i]
    mabso[i]=sum(abs(tmp))/iend
    rmse[i]=sqrt(sum(tmp^2)/iend)
    bias[i]=mean(tmp)
  }
  print("RMSE of out-of-sample forecasts")
  print(rmse)
  print("Mean absolute error of out-of-sample forecasts")
  print(mabso)
  print("Bias of out-of-sample forecasts")
  print(bias)
  backtest <- list(origin=orig,error=err,rmse=rmse,mabso=mabso,bias=bias)
}



#' Backtest for Univariate TAR Models
#'
#' Perform back-test of a univariate SETAR model.
#' @param model SETAR model.
#' @param orig forecast origin.
#' @param h forecast horizon.
#' @param iter number of iterations.
#' @return \code{backTAR} returns a list of components:
#' \item{model}{SETAR model.}
#' \item{error}{prediction errors.}
#' \item{State}{predicted states.}
#' @export
"backTAR" <- function(model,orig,h=1,iter=3000){
  y <- model$data
  order <- model$arorder
  thr1 <- c(model$thr)
  thr2 <- c(min(y)-1,thr1,max(y)+1)
  cnst <- model$cnst
  d1 <- model$delay
  nT <- length(y)
  if(orig < 1)orig <- nT-1
  if(orig >= nT) orig <- nT-1
  Err <- NULL
  Y <- c(y,rep(mean(y),h))
  nregime <- length(thr1)+1
  State <- rep(nregime,(nT-orig))
  for (n in orig:(nT-1)){
    y1 <- y[1:n]
    v1 <- y[n-d1]
    jj = nregime
    for (j in 1:nregime){
      if((v1 > thr2[j]) && (v1 <= thr2[j+1]))jj=j
    }
    State[n-orig+1] = jj
    m1 <- uTAR.est(y1,arorder=order,thr=thr1,d=d1,include.mean=cnst,output=FALSE)
    m2 <- uTAR.pred(m1,n,h=h,iterations=iter,output=FALSE)
    er <- Y[(n+1):(n+h)]-m2$pred
    if(h == 1) {Err=c(Err,er)
    }
    else{Err <- rbind(Err,matrix(er,1,h))}
  }
  ## summary statistics
  MSE <- NULL
  MAE <- NULL
  Bias <- NULL
  if(h == 1){
    MSE = mean(Err^2)
    MAE = mean(abs(Err))
    Bias = mean(Err)
    nf <- length(Err)
  }else{
    nf <- nrow(Err)
    for (j in 1:h){
      err <- Err[1:(nf-j+1),j]
      MSE <- c(MSE,mean(err^2))
      MAE <- c(MAE,mean(abs(err)))
      Bias <- c(Bias,mean(err))
    }
  }
  cat("Starting forecast origin: ",orig,"\n")
  cat("1-step to ",h,"-step out-sample forecasts","\n")
  cat("RMSE: ",sqrt(MSE),"\n")
  cat(" MAE: ",MAE,"\n")
  cat("Bias: ",Bias,"\n")
  cat("Performance based on the regime of forecast origins: ","\n")
  for (j in 1:nregime){
    idx=c(1:nf)[State==j]
    if(h == 1){
      Errj=Err[idx]
      MSEj = mean(Errj^2)
      MAEj = mean(abs(Errj))
      Biasj = mean(Errj)
    }else{
      MSEj <- NULL
      MAEj <- NULL
      Biasj <- NULL
      for (i in 1:h){
        idx=c(1:(nf-i+1))[State[1:(nf-i+1)]==j]
        err = Err[idx,i]
        MSEj <- c(MSEj,mean(err^2))
        MAEj <- c(MAEj,mean(abs(err)))
        Biasj <- c(Biasj,mean(err))
      }
    }
    cat("Summary Statistics when forecast origins are in State: ",j,"\n")
    cat("Number of forecasts used: ",length(idx),"\n")
    cat("RMSEj: ",sqrt(MSEj),"\n")
    cat(" MAEj: ",MAEj,"\n")
    cat("Biasj: ",Biasj,"\n")
  }
  backTAR <- list(model=model,error=Err,State=State)
}


#' Rank-Based Portmanteau Tests
#'
#' Performs rank-based portmanteau statistics.
#' @param zt time series.
#' @param lag the maximum lag to calculate the test statistic.
#' @param output a logical value for output. Default is TRUE.
#' @return \code{rankQ} function outputs the test statistics and p-values for Portmanteau tests, and returns a list with components:
#' \item{Qstat}{test statistics.}
#' \item{pv}{p-values.}
#' @examples
#' phi=t(matrix(c(-0.3, 0.5,0.6,-0.3),2,2))
#' y=uTAR.sim(nob=2000, arorder=c(2,2), phi=phi, d=2, thr=0.2, cnst=c(1,-1),sigma=c(1, 1))
#' rankQ(y$series,10,output=TRUE)
#' @export
"rankQ" <- function(zt,lag=10,output=TRUE){
  nT <- length(zt)
  Rt <- rank(zt)
  erho <- vrho <- rho <- NULL
  rbar <- (nT+1)/2
  sse <- nT*(nT^2-1)/12
  Qstat <- NULL
  pv <- NULL
  rmRt <- Rt - rbar
  deno <- 5*(nT-1)^2*nT^2*(nT+1)
  n1 <- 5*nT^4
  Q <- 0
  for (i in 1:lag){
    tmp <- crossprod(rmRt[1:(nT-i)],rmRt[(i+1):nT])
    tmp <- tmp/sse
    rho <- c(rho,tmp)
    er <- -(nT-i)/(nT*(nT-1))
    erho <- c(erho,er)
    vr <- (n1-(5*i+9)*nT^3+9*(i-2)*nT^2+2*i*(5*i+8)*nT+16*i^2)/deno
    vrho <- c(vrho,vr)
    Q <- Q+(tmp-er)^2/vr
    Qstat <- c(Qstat,Q)
    pv <- c(pv,1-pchisq(Q,i))
  }
  if(output){
    cat("Results of rank-based Q(m) statistics: ","\n")
    Out <- cbind(c(1:lag),rho,Qstat,pv)
    colnames(Out) <- c("Lag","ACF","Qstat","p-value")
    print(round(Out,3))
  }

  rankQ <- list(rho = rho, erho = erho, vrho=vrho,Qstat=Qstat,pv=pv)
}


#' F Test for Nonlinearity
#'
#' Compute the F-test statistic for nonlinearity
#' @param x time series.
#' @param order AR order.
#' @param thres threshold value.
#' @return The function outputs the test statistic and its p-value, and return a list with components:
#' \item{test.stat}{test statistic.}
#' \item{p.value}{p-value.}
#' \item{order}{AR order.}
#' @examples
#' y=rnorm(100)
#' F.test(y,2,0)
#' @export
"F.test" <- function(x,order,thres=0.0){
  if (missing(order)) order=ar(x)$order
  m <- order
  ist <- m+1
  nT <- length(x)
  y <- x[ist:nT]
  X <- NULL
  for (i in 1:m) X <- cbind(X,x[(ist-i):(nT-i)])
  lm1 <- lm(y~X)
  a1 <- summary(lm1)
  coef <- a1$coefficient[-1,]
  idx <- c(1:m)[abs(coef[,3]) > thres]
  jj <- length(idx)
  if(jj==0){
    idx <- c(1:m)
    jj <- m
  }
  for (j in 1:jj){
    for (k in 1:j){
      X <- cbind(X,x[(ist-idx[j]):(nT-idx[j])]*x[(ist-idx[k]):(nT-idx[k])])
    }
  }
  lm2 <- lm(y~X)
  a2 <- anova(lm1,lm2)
  list(test.stat = signif(a2[[5]][2],4),p.value=signif(a2[[6]][2],4),order=order)
}

#' ND Test
#'
#' Compute the ND test statistic of Pena and Rodriguez (2006, JSPI).
#' @param x time series.
#' @param m the maximum number of lag of correlation to test.
#' @param p AR order.
#' @param q MA order.
#' @references
#' Pena, D., and Rodriguez, J. (2006) A powerful Portmanteau test of lack of fit for time series. series. \emph{Journal of American Statistical Association}, 97, 601-610.
#' @return \code{PRnd} function outputs the ND test statistic and its p-value.
#' @examples
#' y=arima.sim(n=500,list(ar=c(0.8,-0.6,0.7)))
#' PRnd(y,10,3,0)
#' @export
"PRnd" <- function(x,m=10,p=0,q=0){
  pq <- p+q
  nu <- 3*(m+1)*(m-2*pq)^2
  de <- (2*m*(2*m+1)-12*(m+1)*pq)*2
  alpha <- nu/de
  #cat("alpha: ",alpha,"\n")
  beta <- 3*(m+1)*(m-2*pq)
  beta <- beta/(2*m*(m+m+1)-12*(m+1)*pq)
  #
  n1 <- 2*(m/2-pq)*(m*m/(4*(m+1))-pq)
  d1 <- 3*(m*(2*m+1)/(6*(m+1))-pq)^2
  r1 <- 1-(n1/d1)
  lambda <- 1/r1
  cat("lambda: ",lambda,"\n")
  if(is.matrix(x))x <- c(x[,1])
  nT <- length(x)
  adjacf <- c(1)
  xwm <- scale(x,center=T,scale=F)
  deno <- crossprod(xwm,xwm)
  if(nT > m){
    for (i in 1:m){
      nn <- crossprod(xwm[1:(nT-i)],xwm[(i+1):nT])
      adjacf <- c(adjacf,(nT+2)*nn/((nT-i)*deno))
    }
    ## cat("adjacf: ",adjacf,"\n")
  }
  else{
    adjacf <- c(adjacf,rep(0,m))
  }
  Rm <- adjacf
  tmp <- adjacf
  for (i in 1:m){
    tmp <- c(adjacf[i+1],tmp[-(m+1)])
    Rm <- rbind(Rm,tmp)
  }
  #cat("Rm: ","\n")
  #print(Rm)

  tst <- -(nT/(m+1))*log(det(Rm))
  a1 <- alpha/beta
  nd <- a1^(-r1)*(lambda/sqrt(alpha))*(tst^r1-a1^r1*(1-(lambda-1)/(2*alpha*lambda^2)))
  pv <- 2*(1-pnorm(abs(nd)))
  cat("ND-stat & p-value ",c(nd,pv),"\n")
}


#' Create Dummy Variables for High-Frequency Intraday Seasonality
#'
#' Create dummy variables for high-frequency intraday seasonality.
#' @param int length of time interval in minutes.
#' @param Fopen number of dummies/intervals from the market open.
#' @param Tend number of dummies/intervals to the market close.
#' @param days number of trading days in the data.
#' @param pooled a logical value indicating whether the data are pooled.
#' @param skipmin the number of minites omitted from the opening.
#' @examples
#' x=hfDummy(5,Fopen=4,Tend=4,days=2,skipmin=15)
#' @export
"hfDummy" <- function(int=1,Fopen=10,Tend=10,days=1,pooled=1,skipmin=0){
  nintval=(6.5*60-skipmin)/int
  Days=matrix(1,days,1)
  X=NULL
  ini=rep(0,nintval)
  for (i in 1:Fopen){
    x1=ini
    x1[i]=1
    Xi=kronecker(Days,matrix(x1,nintval,1))
    X=cbind(X,Xi)
  }
  for (j in 1:Tend){
    x2=ini
    x2[nintval-j+1]=1
    Xi=kronecker(Days,matrix(x2,nintval,1))
    X=cbind(X,Xi)
  }
  X1=NULL
  if(pooled > 1){
    nh=floor(Fopen/pooled)
    rem=Fopen-nh*pooled
    if(nh > 0){
      for (ii in 1:nh){
        ist=(ii-1)*pooled
        y=apply(X[,(ist+1):(ist+pooled)],1,sum)
        X1=cbind(X1,y)
      }
    }
    if(rem > 0){
      X1=cbind(X1,X[,(Fopen-rem):Fopen])
    }
    nh=floor(Tend/pooled)
    rem=Tend-nh*pooled
    if(nh > 0){
      for (ii in 1:nh){
        ist=(ii-1)*pooled
        y=apply(X[,(Fopen+ist+1):(Fopen+ist+pooled)],1,sum)
        X1=cbind(X1,y)
      }
    }
    if(rem > 0){
      X1=cbind(X1,X[,(Fopen+Tend-rem):(Fopen+Tend)])
    }

    X=X1
  }

  hfDummy <- X
}


"factorialOwn" <- function(n,log=T){
  x=c(1:n)
  if(log){
    x=log(x)
    y=cumsum(x)
  }
  else{
    y=cumprod(x)
  }
  y[n]
}


#' Estimate Time-Varying Coefficient AR Models
#'
#' Estimate time-varying coefficient AR models.
#' @param x a time series of data.
#' @param lags the lagged variables used, e.g. lags=c(1,3) means lag-1 and lag-3 are used as regressors.
#' It is more flexible than specifying an order.
#' @param include.mean a logical value indicating whether the constant terms are included.
#' @return \code{trAR} function returns the value from function \code{dlmMLE}.
#' @examples
#' t=50
#' x=rnorm(t)
#' phi1=matrix(0.4,t,1)
#' for (i in 2:t){
#'    phi1[i]=0.7*phi1[i-1]+rnorm(1,0,0.1)
#' 	x[i]=phi1[i]*x[i-1]+rnorm(1)
#' }
#' est=tvAR(x,1)
#' @import dlm
#' @export
"tvAR" <- function(x,lags=c(1),include.mean=TRUE){
  if(is.matrix(x))x <- c(x[,1])
  nlag <-  length(lags)
  if(nlag > 0){p <- max(lags)
  }else{
    p=0; inlcude.mean=TRUE}
  ist <- p+1
  nT <- length(x)
  nobe <- nT-p
  X <- NULL
  if(include.mean)X <- rep(1,nobe)
  if(nlag > 0){
    for (i in 1:nlag){
      ii = lags[i]
      X <- cbind(X,x[(ist-ii):(nT-ii)])
    }
  }
  X <- as.matrix(X)

  if(p > 0){c1 <- paste("lag",lags,sep="")
  }else{c1 <- NULL}
  if(include.mean) c1 <- c("cnt",c1)
  colnames(X) <- c1
  k <- ncol(X)
  y <- x[ist:nT]
  m1 <- lm(y~-1+X)
  coef <- m1$coefficients
  #### Estimation
  build <- function(parm){
    if(k == 1){
      dlm(FF=matrix(rep(1,k),nrow=1),V=exp(parm[1]), W=exp(parm[2]),
          GG=diag(k),m0=coef, C0 = diag(k),JFF=matrix(c(1:k),nrow=1),X=X)
    }else{
      dlm(FF=matrix(rep(1,k),nrow=1),V=exp(parm[1]), W=diag(exp(parm[-1])),
          GG=diag(k),m0=coef, C0 = diag(k),JFF=matrix(c(1:k),nrow=1),X=X)
    }
  }
  mm <- dlmMLE(y,rep(-0.5,k+1),build)
  par <- mm$par; value <- mm$value; counts <- mm$counts

  cat("par:", par,"\n")
  tvAR <-list(par=par,value=value,counts=counts,convergence=mm$convergence,message=mm$message)
}

#' Filtering and Smoothing for Time-Varying AR Models
#'
#' This function performs forward filtering and backward smoothing for a fitted time-varying AR model with parameters in 'par'.
#' @param x a time series of data.
#' @param lags the lag of AR order.
#' @param par the fitted time-varying AR models. It can be an object returned by function. \code{tvAR}.
#' @param include.mean a logical value indicating whether the constant terms are included.
#' @examples
#' t=50
#' x=rnorm(t)
#' phi1=matrix(0.4,t,1)
#' for (i in 2:t){
#'    phi1[i]=0.7*phi1[i-1]+rnorm(1,0,0.1)
#' 	x[i]=phi1[i]*x[i-1]+rnorm(1)
#' }
#' est=tvAR(x,1)
#' tvARFiSm(x,1,FALSE,est$par)
#' @return \code{trARFiSm} function return values returned by function \code{dlmFilter} and \code{dlmSmooth}.
#' @import dlm
#' @export
"tvARFiSm" <- function(x,lags=c(1),include.mean=TRUE,par){
  if(is.matrix(x))x <- c(x[,1])
  nlag <-  length(lags)
  if(nlag > 0){p <- max(lags)
  }else{
    p=0; inlcude.mean=TRUE}
  ist <- p+1
  nT <- length(x)
  nobe <- nT-p
  X <- NULL
  if(include.mean)X <- rep(1,nobe)
  if(nlag > 0){
    for (i in 1:nlag){
      ii = lags[i]
      X <- cbind(X,x[(ist-ii):(nT-ii)])
    }
  }
  X <- as.matrix(X)
  k <- ncol(X)
  y <- x[ist:nT]
  m1 <- lm(y~-1+X)
  coef <- m1$coefficients
  ### Model specification
  if(k == 1){
    tvAR <- dlm(FF=matrix(rep(1,k),nrow=1),V=exp(par[1]), W=exp(par[2]),
                GG=diag(k),m0=coef, C0 = diag(k),JFF=matrix(c(1:k),nrow=1),X=X)
  }else{
    tvAR <- dlm(FF=matrix(rep(1,k),nrow=1),V=exp(par[1]), W=diag(exp(par[-1])),
                GG=diag(k),m0=coef, C0 = diag(k),JFF=matrix(c(1:k),nrow=1),X=X)
  }
  mF <- dlmFilter(y,tvAR)
  mS <- dlmSmooth(mF)

  tvARFiSm <- list(filter=mF,smooth=mS)
}

#' Estimating of Random-Coefficient AR Models
#'
#' Estimate random-coefficient AR models.
#' @param x a time series of data.
#' @param lags the lag of AR models. This is more flexible than using order. It can skip unnecessary lags.
#' @param include.mean a logical value indicating whether the constant terms are included.
#' @return \code{rcAR} function returns a list with following components:
#' \item{par}{estimated parameters.}
#' \item{se.est}{standard errors.}
#' \item{residuals}{residuals.}
#' \item{sresiduals}{standardized residuals.}
#' @examples
#' t=50
#' x=rnorm(t)
#' phi1=matrix(0.4,t,1)
#' for (i in 2:t){
#'    phi1[i]=0.7*phi1[i-1]+rnorm(1,0,0.1)
#' 	x[i]=phi1[i]*x[i-1]+rnorm(1)
#' }
#' est=rcAR(x,1,FALSE)
#' @import dlm
#' @export
"rcAR" <- function(x,lags=c(1),include.mean=TRUE){
  if(is.matrix(x))x <- c(x[,1])
  if(include.mean){
    mu <- mean(x)
    x <- x-mu
    cat("Sample mean: ",mu,"\n")
  }
  nlag <-  length(lags)
  if(nlag  < 1){lags <- c(1); nlag <- 1}
  p <- max(lags)
  ist <- p+1
  nT <- length(x)
  nobe <- nT-p
  X <- NULL
  for (i in 1:nlag){
    ii = lags[i]
    X <- cbind(X,x[(ist-ii):(nT-ii)])
  }
  X <- as.matrix(X)
  k <- ncol(X)
  y <- x[ist:nT]
  m1 <- lm(y~-1+X)
  par <- m1$coefficients
  par <- c(par,rep(-3.0,k),-0.3)
  ##cat("initial estimates: ",par,"\n")
  ##
  rcARlike <- function(par,y=y,X=X){
    k <- ncol(X)
    nobe <- nrow(X)
    sigma2 <- exp(par[length(par)])
    if(k > 1){
      beta <- matrix(par[1:k],k,1)
      yhat <- X%*%beta
      gamma <- matrix(exp(par[(k+1):(2*k)]),k,1)
      sig <- X%*%gamma
    }else{
      beta <- par[1]; gamma <- exp(par[2])
      yhat <- X*beta
      sig <- X*gamma
    }
    at <- y-yhat
    v1 <- sig+sigma2
    ll <- sum(dnorm(at,mean=rep(0,nobe),sd=sqrt(v1),log=TRUE))
    rcARlike <- -ll
  }
  #
  mm <- optim(par,rcARlike,y=y,X=X,hessian=TRUE)
  est <- mm$par
  H <- mm$hessian
  Hi <- solve(H)
  se.est <- sqrt(diag(Hi))
  tratio <- est/se.est
  tmp <- cbind(est,se.est,tratio)
  cat("Estimates:","\n")
  print(round(tmp,4))
  ### Compute residuals and standardized residuals
  sigma2 <- exp(est[length(est)])
  if(k > 1){
    beta <- matrix(est[1:k],k,1)
    yhat <- X%*%beta
    gamma <- matrix(exp(est[(k+1):(2*k)]),k,1)
    sig <- X%*%gamma
  }else{
    beta <- est[1]; gamma <- exp(est[2])
    yhat <- X*beta
    sig <- X*gamma
  }
  at <- y-yhat
  v1 <- sig+sigma2
  sre <- at/sqrt(v1)
  #
  rcAR <- list(par=est,se.est=se.est,residuals=at,sresiduals=sre)
}



#' Generic Sequential Monte Carlo Method
#'
#' Function of generic sequential Monte Carlo method with delay weighting not using full information proposal distribution.
#' @param Sstep a function that performs one step propagation using a proposal distribution.
#' Its input includes \code{(mm,xx,logww,yyy,par,xdim,ydim)}, where
#' \code{xx} and \code{logww} are the last iteration samples and log weight. \code{yyy} is the
#' observation at current time step. It should return \code{xx} (the samples xt) and
#' \code{logww} (their corresponding log weight).
#' @param nobs the number of observations \code{T}.
#' @param yy the observations with \code{T} columns and \code{ydim} rows.
#' @param mm the Monte Carlo sample size.
#' @param par a list of parameter values to pass to \code{Sstep}.
#' @param xx.init the initial samples of \code{x_0}.
#' @param xdim the dimension of the state variable \code{x_t}.
#' @param ydim the dimension of the observation \code{y_t}.
#' @param resample.sch a binary vector of length \code{nobs}, reflecting the resampling schedule. resample.sch[i]= 1 indicating resample should be carried out at step \code{i}.
#' @param delay the maximum delay lag for delayed weighting estimation. Default is zero.
#' @param funH a user supplied function \code{h()} for estimation \code{E(h(x_t) | y_t+d}). Default
#' is identity for estimating the mean. The function should be able to take vector or matrix as input and operates on each element of the input.
#' @references
#' Tsay, R. and Chen, R. (2018). Nonlinear Time Series Analysis. John Wiley & Sons, New Jersey.
#' @return The function returns \code{xhat}, an array with dimensions \code{(xdim; nobs; delay+1)},
#' and the scaled log-likelihood value \code{loglike}. If \code{loglike} is needed, the log weight
#' calculation in the \code{Sstep} function should retain all constants that are related to
#' the parameters involved. Otherwise, \code{Sstep} function may remove all constants
#' that are common to all the Monte Carlo samples. It needs a utility function
#' \code{circular2ordinal}, also included in the \code{NTS} package, for efficient memory management.
#' @examples
#' nobs= 100; pd= 0.95; ssw= 0.1; ssv= 0.5;
#' xx0= 0; ss0= 0.1; nyy= 50;
#' yrange= c(-80,80); xdim= 2; ydim= nyy;
#' mm= 10000
#' yr=yrange[2]-yrange[1]
#' par=list(ssw=ssw,ssv=ssv,nyy=nyy,pd=pd,yr=yr)
#' simu=simuTargetClutter(nobs,pd,ssw,ssv,xx0,ss0,nyy,yrange)
#' xx.init=matrix(nrow=2,ncol=mm)
#' xx.init[1,]=yrange[1]+runif(mm)*yr
#' xx.init[2,]=rep(0.1,mm)
#' resample.sch=rep.int(1,nobs)
#' out= SMC(Sstep.Clutter,nobs,simu$yy,mm,par,xx.init,xdim,ydim,resample.sch)
#' @import tensor
#' @export
SMC=function(Sstep,nobs,yy,mm,par,xx.init,xdim,ydim,
             resample.sch,delay=0,funH=identity){
  #---------------------------------------------------------------
  xxd <- array(dim=c(xdim,mm,delay+1))
  delay.front = 1
  loglike=0
  xx <- xx.init
  if(xdim==1) xx=matrix(xx,nrow=1,ncol=mm)
  xxd[,,1] <- funH(xx)
  xhat <- array(dim=c(xdim,nobs,delay+1))
  logww <- rep(0,mm)
  for(i in 1:nobs){
    if(ydim==1){
      yyy <- yy[i]
    }else{
      yyy <- yy[,i]
    }
    step <- Sstep(mm,xx,logww,yyy,par,xdim,ydim)
    xx <- step$xx
    logww <- step$logww-max(step$logww)
    loglike <- loglike+max(step$logww)
    ww <- exp(logww)
    sumww <- sum(ww)
    ww <- ww/sumww
    xxd[,,delay.front] <- funH(xx)
    delay.front <- delay.front%%(delay+1) + 1
    order <- circular2ordinal(delay.front, delay+1)
    if(xdim==1) ww=as.vector(ww)  ## correction: insert
    xhat[,i,] <- tensor(xxd,ww,2,1)[,order]
    if(resample.sch[i]==1){
      r.index <- sample.int(mm,size=mm,replace=T,prob=ww)
      xx[,] <- xx[,r.index]
      xxd[,,] <- xxd[,r.index,]
      logww <- logww*0
      loglike <- loglike+log(sumww)
    }
  }
  if(resample.sch[nobs]==0){
    ww <- exp(logww)
    sumww=sum(ww)
    loglike=loglike+log(sumww)
  }
  if(delay > 0){
    for(dd in 1:delay){
      xhat[,1:(nobs-dd),dd+1] <- xhat[,(dd+1):nobs,dd+1]
      xhat[,(nobs-dd+1):nobs,dd+1] <- NA
    }
  }
  return(list(xhat=xhat,loglike=loglike))
}


ar2natural=function(par.ar){
  par.natural=par.ar
  par.natural[1]=par.ar[1]/(1-par.ar[2])
  par.natural[3]=par.ar[3]/sqrt(1-par.ar[2]**2)
  return(par.natural)
}
natural2ar=function(par.natural){
  par.ar=par.natural
  par.ar[1]=par.natural[1]*(1-par.natural[2])
  par.ar[3]=par.natural[3]*sqrt(1-par.natural[2]**2)
  return(par.ar)
}


circular2ordinal=function(circular.front, circule.length){
  output=c()
  if(circular.front>1){
    output = c(output, seq(circular.front-1,1,-1))
  }
  output = c(output, seq(circule.length, circular.front,-1))
  return(output)
}

Sstep.SV=function(mm,xx,logww,yyy,par,xdim,ydim){
  xx[1,] <- par[1]+par[2]*xx[1,]+rnorm(mm)*par[3]
  logww=logww-0.5*(yyy**2)*exp(-2*xx[1,])-xx[1,]
  return(list(xx=xx,logww=logww))
}


#' Generic Sequential Monte Carlo Smoothing with Marginal Weights
#'
#' Generic sequential Monte Carlo smoothing with marginal weights.
#' @param SISstep a function that performs one propagation step using a proposal distribution.
#' Its input includes \code{(mm,xx,logww,yyy,par,xdim,ydim)}, where
#' \code{xx} and \code{logww} are the last iteration samples and log weight. \code{yyy} is the
#' observation at current time step. It should return {xx} (the samples xt) and
#' {logww} (their corresponding log weight).
#' @param SISstep.Smooth the function for backward smoothing step.
#' @param nobs the number of observations \code{T}.
#' @param yy the observations with \code{T} columns and \code{ydim} rows.
#' @param mm the Monte Carlo sample size \code{m}.
#' @param par a list of parameter values.
#' @param xx.init the initial samples of \code{x_0}.
#' @param xdim the dimension of the state variable \code{x_t}.
#' @param ydim the dimension of the observation \code{y_t}.
#' @param resample.sch a binary vector of length \code{nobs}, reflecting the resampling schedule. resample.sch[i]= 1 indicating resample should be carried out at step \code{i}.
#' @param funH a user supplied function \code{h()} for estimation \code{E(h(x_t) | y_1,...,y_T}). Default
#' is identity for estimating the mean. The function should be able to take vector or matrix as input and operates on each element of the input.
#' @references
#' Tsay, R. and Chen, R. (2018). Nonlinear Time Series Analysis. John Wiley & Sons, New Jersey.
#' @return The function returns the smoothed values.
#' @export

SMC.Smooth=function(SISstep,SISstep.Smooth,nobs,yy,mm,par,
                    xx.init,xdim,ydim,resample.sch,funH=identity){
  #---------------------------------------------------------------
  xxall <- array(dim=c(xdim,mm,nobs))
  wwall <- matrix(nrow=mm,ncol=nobs)
  xhat <- matrix(nrow=xdim,ncol=nobs)
  xx <- xx.init
  if(xdim==1) xx=matrix(xx,nrow=1,ncol=mm)  ## correction: insert
  logww <- rep(0,mm)
  for(i in 1:nobs){
    # print(c('forward',i));flush.console();
    if(ydim==1){
      yyy <- yy[i]
    }else{
      yyy <- yy[,i]
    }
    step <- SISstep(mm,xx,logww,yyy,par,xdim,ydim)
    xx <- step$xx
    logww <- step$logww-max(step$logww)
    ww <- exp(logww)
    ww <- ww/sum(ww)
    xxall[,,i] <- xx
    wwall[,i] <- ww
    if(resample.sch[i]==1){
      r.index <- sample.int(mm,size=mm,replace=T,prob=ww)
      xx[,] <- xx[,r.index]
      logww <- logww*0
    }
  }
  vv <- wwall[,nobs]
  xhat[,nobs] <- funH(xxall[,,nobs])%*%vv
  for(i in (nobs-1):1){
    # print(c('backward',i));flush.console();
    xxt <- xxall[,,i]
    xxt1 <- xxall[,,i+1]
    ww <- wwall[,i]
    step2 <- SISstep.Smooth(mm,xxt,xxt1,ww,vv,par)
    vv <- step2$vv
    vv <- vv/sum(vv)
    xhat[,i] <- funH(xxall[,,i])%*%vv
  }
  return(list(xhat=xhat))
}




#' Simulate A Moving Target in Clutter
#'
#' The function simulates a target signal under clutter environment.
#' @param nobs the number observations.
#' @param pd the probability to observe the true signal.
#' @param ssw the standard deviation in the state equation.
#' @param ssv the standard deviation for the observation noise.
#' @param xx0 the initial location.
#' @param ss0 the initial speed.
#' @param nyy the dimension of the data.
#' @param yrange the range of data.
#' @references
#' Tsay, R. and Chen, R. (2018). Nonlinear Time Series Analysis. John Wiley & Sons, New Jersey.
#' @return The function returns a list with components:
#' \item{xx}{the location.}
#' \item{ss}{the speed.}
#' \item{ii}{the indicators for whether the observation is the true signal.}
#' \item{yy}{the data.}
#' @examples
#' data=simuTargetClutter(30,0.5,0.5,0.5,0,0.3,3,c(-30,30))
#' @export
simuTargetClutter=function(nobs,pd,ssw,ssv,xx0,ss0,nyy,yrange){
  xx <- 1:nobs
  ss <- 1:nobs
  xx[1] <- xx0
  ss[1] <- ss0
  ww <- rnorm(nobs)*ssw
  vv <- rnorm(nobs)*ssv
  for(i in 2:nobs){
    xx[i] <- xx[i-1]+ss[i-1]+0.5*ww[i]
    ss[i] <- ss[i-1]+ww[i]
  }
  ii <- floor(runif(nobs)+pd)
  temp <- sum(ii)
  ii[ii==1] <- floor(runif(temp)*nyy+1)
  yy <- matrix(runif(nobs*nyy),ncol=nobs,nrow=nyy)
  yy <- yrange[1]+yy*(yrange[2]-yrange[1])
  for(i in 1:nobs){
    if(ii[i]>0) yy[ii[i],i] <- xx[i]+vv[i]
  }
  return(list(xx=xx,ss=ss,ii=ii,yy=yy))
}


#' Sequential Monte Carlo for A Moving Target under Clutter Environment
#'
#' The function performs one step propagation using the sequential Monte Carlo method with partial state proposal for tracking in clutter problem.
#' @param mm the Monte Carlo sample size \code{m}.
#' @param xx the sample in the last iteration.
#' @param logww the log weight in the last iteration.
#' @param yyy the observations.
#' @param par a list of parameter values \code{(ssw,ssv,pd,nyy,yr)}, where \code{ssw} is the standard deviation in the state equation,
#' \code{ssv} is the standard deviation for the observation noise, \code{pd} is the probability to observe the true signal, \code{nyy} the dimension of the data,
#' and \code{yr} is the range of the data.
#' @param xdim the dimension of the state varible.
#' @param ydim the dimension of the observation.
#' @references
#' Tsay, R. and Chen, R. (2018). Nonlinear Time Series Analysis. John Wiley & Sons, New Jersey.
#' @return The function returns a list with the following components:
#' \item{xx}{the new sample.}
#' \item{logww}{the log weights.}
#' @examples
#' nobs <- 100; pd <- 0.95; ssw <- 0.1; ssv <- 0.5;
#' xx0 <- 0; ss0 <- 0.1; nyy <- 50;
#' yrange <- c(-80,80); xdim <- 2; ydim <- nyy;
#' simu <- simuTargetClutter(nobs,pd,ssw,ssv,xx0,ss0,nyy,yrange)
#' resample.sch <- rep(1,nobs)
#' mm <- 10000
#' yr <- yrange[2]-yrange[1]
#' par <- list(ssw=ssw,ssv=ssv,nyy=nyy,pd=pd,yr=yr)
#' yr<- yrange[2]-yrange[1]
#' xx.init <- matrix(nrow=2,ncol=mm)
#' xx.init[1,] <- yrange[1]+runif(mm)*yr
#' xx.init[2,] <- rep(0.1,mm)
#' out <- SMC(Sstep.Clutter,nobs,simu$yy,mm,par,xx.init,xdim,ydim,resample.sch)
#' @export
Sstep.Clutter=function(mm,xx,logww,yyy,par,xdim,ydim){
  ee <- rnorm(mm)*par$ssw
  xx[1,] <- xx[1,]+xx[2,]+0.5*ee
  xx[2,] <- xx[2,]+ee
  uu <- 1:mm
  for(j in 1:mm){
    temp <- sum(exp(-(yyy-xx[1,j])**2/2/par$ssv**2))
    uu[j] <- temp/par$ssv/sqrt(2*pi)*par$pd/par$nyy+(1-par$pd)/par$yr
  }
  logww <- logww+log(uu)
  return(list(xx=xx,logww=logww))
}

#' Kalman Filter for Tracking in Clutter
#'
#' This function implements Kalman filter to track a moving target under clutter environment with known indicators.
#' @param nobs the number of observations.
#' @param ssw the standard deviation in the state equation.
#' @param ssv the standard deviation for the observation noise.
#' @param yy the data.
#' @param ii the indicators.
#' @references
#' Tsay, R. and Chen, R. (2018). Nonlinear Time Series Analysis. John Wiley & Sons, New Jersey.
#' @return The function returns a list with the following components:
#' \item{xhat}{the fitted location.}
#' \item{shat}{the fitted speed.}
#' @examples
#' nobs <- 100; pd <- 0.95; ssw <- 0.1; ssv <- 0.5;
#' xx0 <- 0; ss0 <- 0.1; nyy <- 50;
#' yrange <- c(-80,80); xdim <- 2; ydim <- nyy;
#' simu <- simuTargetClutter(nobs,pd,ssw,ssv,xx0,ss0,nyy,yrange)
#' outKF <- clutterKF(nobs,ssw,ssv,simu$yy,simu$ii)
#' @export
clutterKF=function(nobs,ssw,ssv,yy,ii){
  #---------------------------------------------------------------
  xhat <- 1:nobs;shat <- 1:nobs
  bb <- c(0,0);cc <- 0
  HH <- matrix(c(1,1,0,1),ncol=2,nrow=2,byrow=TRUE)
  GG <- matrix(c(1,0),ncol=2,nrow=1,byrow=TRUE)
  WW <- matrix(c(0.5*ssw,ssw),ncol=1,nrow=2,byrow=TRUE)
  VV <- ssv
  mu <- rep(0,2)
  SS <- diag(c(1,1))*1000
  for(i in 1:nobs){
    outP <- KF1pred(mu,SS,HH,bb,WW)
    mu <- outP$mu
    SS <- outP$SS
    if(ii[i] > 0){
      outU <- KF1update(mu,SS,yy[ii[i],i],GG,cc,VV)
      mu <- outU$mu
      SS <- outU$SS
    }
    xhat[i] <- mu[1]
    shat[i] <- mu[2]
  }
  return(list(xhat=xhat,shat=shat))
}

KF1pred=function(mu,SS,HH,bb,WW){
  mu=HH%*%mu+bb
  SS=HH%*%SS%*%t(HH)+WW%*%t(WW)
  return(list(mu=mu,SS=SS))
}




#' Simulate A Sample Trajectory
#'
#' The function generates a sample trajectory of the target and the corresponding observations with sensor locations at (0,0)
#' and (20,0).
#' @param nn sample size.
#' @param q contains the information about the covariance of the noise.
#' @param r contains the information about \code{V}, where \code{V*t(V)} is the covariance matrix of the observation noise.
#' @param start the initial value.
#' @param seed the seed of random number generator.
#' @return The function returns a list with components:
#' \item{xx}{the state data.}
#' \item{yy}{the observed data.}
#' \item{H}{the state coefficient matrix.}
#' \item{W}{ \code{W*t(W)} is the state innovation covariance matrix.}
#' \item{V}{\code{V*t(V)} is the observation noise covariance matrix.}
#' @examples
#' s2 <- 20 #second sonar location at (s2,0)
#' q <- c(0.03,0.03)
#' r <- c(0.02,0.02)
#' nobs <- 200
#' start <- c(10,10,0.01,0.01)
#' H <- c(1,0,1,0,0,1,0,1,0,0,1,0,0,0,0,1)
#' H <- matrix(H,ncol=4,nrow=4,byrow=TRUE)
#' W <- c(0.5*q[1], 0,0, 0.5*q[2],q[1],0,0,q[2])
#' W <- matrix(W,ncol=2,nrow=4,byrow=TRUE)
#' V <- diag(r)
#' mu0 <- start
#' SS0 <- diag(c(1,1,1,1))*0.01
#' simu_out <- simPassiveSonar(nobs,q,r,start,seed=20)
#' yy<- simu_out$yy
#' tt<- 100:200
#' plot(simu_out$xx[1,tt],simu_out$xx[2,tt],xlab='x',ylab='y')
#' @export
simPassiveSonar=function(nn=200,q,r,start,seed){
  set.seed(seed)
  s2 <- 20  #secone sonar location at (s2,0)
  H <- c(1,0,1,0,
         0,1,0,1,
         0,0,1,0,
         0,0,0,1)
  H <- matrix(H,ncol=4,nrow=4,byrow=TRUE)
  W <- c(0.5*q[1], 0,
         0, 0.5*q[2],
         q[1],0,
         0,q[2])
  W <- matrix(W,ncol=2,nrow=4,byrow=TRUE)
  V <- diag(r)
  x <- matrix(nrow=4,ncol=nn)
  y <- matrix(nrow=2,ncol=nn)
  for(ii in 1:nn){
    if(ii == 1) x[,ii] <- start
    if(ii > 1) x[,ii] <- H%*%x[,ii-1]+W%*%rnorm(2)
    y[1,ii] <- atan(x[2,ii]/x[1,ii])
    y[2,ii] <- atan(x[2,ii]/(x[1,ii]-s2))
    y[,ii] <- y[,ii]+V%*%rnorm(2)
    y[,ii] <- (y[,ii]+0.5*pi)%%pi-0.5*pi
  }
  return(list(xx=x,yy=y,H=H,W=W,V=V))
}


#' Full Information Propagation Step under Mixture Kalman Filter
#'
#' This function implements the full information propagation step under mixture Kalman filter with full information proposal distribution and Rao-Blackwellization, no delay.
#' @param MKFstep.Full.RB a function that performs one step propagation under mixture Kalman filter, with full information proposal distribution.
#' Its input includes \code{(mm,II,mu,SS,logww,yyy,par,xdim,ydim)}, where
#' \code{II}, \code{mu}, and \code{SS} are the indicators and its corresponding mean and variance matrix of the Kalman filter components in the last iterations.
#' \code{logww} is the log weight of the last iteration. \code{yyy} is the
#' observation at current time step. It should return the Rao-Blackwellization estimation of the mean and variance.
#' @param nobs the number of observations \code{T}.
#' @param yy the observations with \code{T} columns and \code{ydim} rows.
#' @param mm the Monte Carlo sample size \code{m}.
#' @param par a list of parameter values to pass to \code{Sstep}.
#' @param II.init the initial indicators.
#' @param mu.init the initial mean.
#' @param SS.init the initial variance.
#' @param xdim the dimension of the state varible \code{x_t}.
#' @param ydim the dimension of the observation \code{y_t}.
#' @param resample.sch a binary vector of length \code{nobs}, reflecting the resampling schedule. resample.sch[i]= 1 indicating resample should be carried out at step \code{i}.
#' @references
#' Tsay, R. and Chen, R. (2018). Nonlinear Time Series Analysis. John Wiley & Sons, New Jersey.
#' @return The function returns a list with components:
#' \item{xhat}{the fitted value.}
#' \item{xhatRB}{the fitted value using Rao-Blackwellization.}
#' \item{Iphat}{the estimated indicators.}
#' \item{IphatRB}{the estimated indicators using Rao-Blackwellization.}
#' @export
MKF.Full.RB=function(MKFstep.Full.RB,nobs,yy,mm,par,II.init,
                     mu.init,SS.init,xdim,ydim,resample.sch){
  #---------------------------------------------------------------
  mu <- mu.init
  SS= SS.init
  II=II.init
  xhat <- matrix(nrow=xdim,ncol=nobs)
  xhatRB <- matrix(nrow=xdim,ncol=nobs)
  Iphat <- 1:nobs
  IphatRB <- 1:nobs

  logww <- rep(0,mm)
  for(i in 1:nobs){
    if(ydim==1){
      yyy <- yy[i]
    } else{
      yyy <- yy[,i]
    }
    step <- MKFstep.Full.RB(mm,II,mu,SS,logww,yyy,par,xdim,ydim,
                            resample.sch[i])
    mu <- step$mu
    SS <- step$SS
    II <- step$II
    xhatRB[,i] <- step$xhatRB
    xhat[,i] <- step$xhat
    IphatRB[i] <- step$IphatRB
    Iphat[i] <- step$Iphat
    logww <- step$logww-max(step$logww)
  }
  return(list(xhat=xhat,xhatRB=xhatRB,Iphat=Iphat,IphatRB=IphatRB))
}


KFoneLike=function(mu,SS,yy,HH,GG,WW,VV){
  mu=HH%*%mu
  SS=HH%*%SS%*%t(HH)+WW%*%t(WW)
  SSy=GG%*%SS%*%t(GG)+VV%*%t(VV)
  KK=t(solve(SSy,t(SS%*%t(GG))))
  res=yy-GG%*%mu
  mu=mu+KK%*%res
  SS=SS-KK%*%GG%*%SS
  like=0.5*log(det(SSy))+t(res)%*%solve(SSy,res)
  return(list(mu=mu,SS=SS,res=res,like=like))
}



#' One Propagation Step under Mixture Kalman Filter for Fading Channels
#'
#' This function implements the one propagation step under mixture Kalman filter for fading channels.
#' @param mm the Monte Carlo sample size.
#' @param II the indicators.
#' @param mu the mean in the last iteration.
#' @param SS the covariance matrix of the Kalman filter components in the last iteration.
#' @param logww is the log weight of the last iteration.
#' @param yyy the observations with \code{T} columns and \code{ydim} rows.
#' @param par a list of parameter values. \code{HH} is the state coefficient matrix, \code{WW*t(WW)} is the state innovation covariance matrix,
#' \code{VV*t(VV)} is the covariance matrix of the observation noise, \code{GG1} and \code{GG2} are the observation coefficient matrix.
#' @param xdim the dimension of the state variable \code{x_t}.
#' @param ydim the dimension of the observation \code{y_t}.
#' @param resample a binary vector of length \code{obs}, reflecting the resampling schedule. resample.sch[i]= 1 indicating resample should be carried out at step \code{i}.
#' @references
#' Tsay, R. and Chen, R. (2018). Nonlinear Time Series Analysis. John Wiley & Sons, New Jersey.
#' @return The function returns a list with components:
#' \item{xhat}{the fitted value.}
#' \item{xhatRB}{the fitted value using Rao-Blackwellization.}
#' \item{Iphat}{the estimated indicators.}
#' \item{IphatRB}{the estimated indicators using Rao-Blackwellization.}
#' @export
MKFstep.fading=function(mm,II,mu,SS,logww,yyy,par,xdim,ydim,resample){
  HH <- par$HH; WW <- par$WW;
  GG1 <- par$GG1; GG2 <- par$GG2; VV <- par$VV;
  prob1 <- 1:mm; prob2 <- 1:mm; CC <- 1:mm
  mu.i <- array(dim=c(xdim,2,mm)); SS.i <- array(dim=c(xdim,xdim,2,mm));
  xxRB <- matrix(nrow=xdim,ncol=mm)
  IpRB <- 1:mm
  for(jj in 1:mm){
    out1 <- KFoneLike(mu[,jj],SS[,,jj],yyy,HH,GG1,WW,VV)
    out2 <- KFoneLike(mu[,jj],SS[,,jj],yyy,HH,GG2,WW,VV)
    #print(c(jj,II[jj]));flush.console()
    prob1[jj] <- exp(-out1$like)
    prob2[jj] <- exp(-out2$like)
    CC[jj] <- prob1[jj]+prob2[jj]
    mu.i[,1,jj] <- out1$mu
    mu.i[,2,jj] <- out2$mu
    SS.i[,,1,jj] <- out1$SS
    SS.i[,,2,jj] <- out2$SS
    xxRB[,jj] <- (prob1[jj]*mu.i[,1,jj]+prob2[jj]*mu.i[,2,jj])/CC[jj]
    IpRB[jj] <- prob1[jj]/CC[jj]*(2-II[jj])+(1-prob1[jj]/CC[jj])*(II[jj]-1)
    logww[jj] <- logww[jj]+log(CC[jj])
  }
  ww <- exp(logww-max(logww))
  ww <- ww/sum(ww)
  xhatRB <- sum(xxRB*ww)
  IphatRB <- sum(IpRB*ww)
  if(resample==1){
    r.index <- sample.int(mm,size=mm,replace=T,prob=ww)
    mu.i[,,] <- mu.i[,,r.index]
    SS.i[,,,] <- SS.i[,,,r.index]
    prob1 <- prob1[r.index]
    CC <- CC[r.index]
    II <- II[r.index]
    logww <- logww*0
  }
  II.new <- (runif(mm) > prob1/CC)+1
  for(jj in 1:mm){
    mu[,jj] <- mu.i[,II.new[jj],jj]
    SS[,,jj] <- SS.i[,,II.new[jj],jj]
  }
  xhat <- sum(mu*ww)
  Iphat <- sum(((2-II)*(2-II.new)+(II-1)*(II.new-1))*ww)
  return(list(mu=mu,SS=SS,II=II.new,logww=logww,xhat=xhat,
              xhatRB=xhatRB,Iphat=Iphat,IphatRB=IphatRB))
}


KF1update=function(mu,SS,yy,GG,cc,VV){
  KK=t(solve(GG%*%SS%*%t(GG)+VV%*%t(VV),t(SS%*%t(GG))))
  mu=mu+KK%*%(yy-cc-GG%*%mu)
  SS=SS-KK%*%GG%*%SS
  return(list(mu=mu,SS=SS))
}

#' Sequential Importance Sampling Step for A Target with Passive Sonar
#'
#' This function implements one step of the sequential importance sampling method for a target with passive sonar.
#' @param mm the Monte Carlo sample size \code{m}.
#' @param xx the sample in the last iteration.
#' @param logww the log weight in the last iteration.
#' @param yy the observations with \code{T} columns and \code{ydim} rows.
#' @param par a list of parameter values. \code{H} is the state coefficient matrix, \code{W*t(W)} is the state innovation covariance matrix,
#' \code{V*t(V)} is the covariance matrix of the observation noise, \code{s2} is the second sonar location.
#' @param xdim the dimension of the state variable \code{x_t}.
#' @param ydim the dimension of the observation \code{y_t}.
#' @references
#' Tsay, R. and Chen, R. (2018). Nonlinear Time Series Analysis. John Wiley & Sons, New Jersey.
#' @return The function returns a list with the following components:
#' \item{xx}{the new sample.}
#' \item{logww}{the log weights.}
#' @export
Sstep.Sonar=function(mm,xx,logww,yy,par,xdim=1,ydim=1){
  H <- par$H; W <- par$W; V <- par$V; s2 <- par$s2;
  xx <- H%*%xx+W%*%matrix(rnorm(2*mm),nrow=2,ncol=mm)
  y1 <- atan(xx[2,]/xx[1,])
  y2 <- atan(xx[2,]/(xx[1,]-s2))
  res1 <- (yy[1]-y1+0.5*pi)%%pi-0.5*pi
  res2 <- (yy[2]-y2+0.5*pi)%%pi-0.5*pi
  uu <- -res1**2/2/V[1,1]**2-res2**2/2/V[2,2]**2
  logww <- logww+uu
  return(list(xx=xx,logww=logww))
}

#' Sequential Importance Sampling for A Target with Passive Sonar
#'
#' This function uses the sequential importance sampling method to deal with a target with passive sonar for smoothing.
#' @param mm the Monte Carlo sample size \code{m}.
#' @param xxt the sample in the last iteration.
#' @param xxt1 the sample in the next iteration.
#' @param ww  the forward filtering weight.
#' @param vv the backward smoothing weight.
#' @param par a list of parameter values. \code{H} is the state coefficient matrix, and \code{W*t(W)} is the state innovation covariance matrix.
#' @references
#' Tsay, R. and Chen, R. (2018). Nonlinear Time Series Analysis. John Wiley & Sons, New Jersey.
#' @return The function returns a list with the following components:
#' \item{xx}{the new sample.}
#' \item{logww}{the log weights.}
#' @export
Sstep.Smooth.Sonar=function(mm,xxt,xxt1,ww,vv,par){
  H <- par$H; W <- par$W;
  uu <- 1:mm
  aa <- 1:mm
  xxt1p <- H%*%xxt
  for(i in 1:mm){
    res1 <- (xxt1[1,i]-xxt1p[1,])/W[1,1]
    res2 <- (xxt1[2,i]-xxt1p[2,])/W[2,2]
    aa[i] <- sum(exp(-0.5*res1**2-0.5*res2**2)*ww)
  }
  for(j in 1:mm){
    res1 <- (xxt1[1,]-xxt1p[1,j])/W[1,1]
    res2 <- (xxt1[2,]-xxt1p[2,j])/W[2,2]
    uu[j] <- sum(exp(-0.5*res1**2-0.5*res2**2)*vv/aa)
  }
  vv <- ww*uu
  return(list(vv=vv))
}




#' Sequential Importance Sampling Step for Fading Channels
#'
#' This function implements one step of the sequential importance sampling method for fading channels.
#' @param mm the Monte Carlo sample size \code{m}.
#' @param xx the sample in the last iteration.
#' @param logww the log weight in the last iteration.
#' @param yyy the observations with \code{T} columns and \code{ydim} rows.
#' @param par a list of parameter values.  \code{HH} is the state coefficient model, \code{WW*t(WW)} is the state innovation covariance matrix,
#' \code{VV*t(VV)} is the covariance of the observation noise, \code{GG} is the observation model.
#' @param xdim2 the dimension of the state variable \code{x_t}.
#' @param ydim the dimension of the observation \code{y_t}.
#' @references
#' Tsay, R. and Chen, R. (2018). Nonlinear Time Series Analysis. John Wiley & Sons, New Jersey.
#' @return The function returns a list with the following components:
#' \item{xx}{the new sample.}
#' \item{logww}{the log weights.}
#' @export
SISstep.fading=function(mm,xx,logww,yyy,par,xdim2,ydim){
  HH <- par$HH; WW <- par$WW; VV <- par$VV; GG <- par$GG
  xxx <- xx[1:4,]
  xxx <- HH%*%xxx+WW%*%matrix(rnorm(mm),nrow=1,ncol=mm)
  alpha <- GG%*%xxx
  alpha <- as.vector(alpha)
  #print(c(alpha,yyy));flush.console()
  pp1 <- 1/(1+exp(-2*yyy*alpha/VV**2))
  temp1 <- -(yyy-alpha)**2/2/VV**2
  temp2 <- -(yyy+alpha)**2/2/VV**2
  temp.max <- apply(cbind(temp1,temp2),1,max)
  #print(c('temp max',temp.max));flush.console()
  temp1 <- temp1-temp.max
  temp2 <- temp2-temp.max
  loguu <- temp.max+log(exp(temp1)+exp(temp2))
  #print(c('loguu',loguu));flush.console()
  logww <- logww+loguu
  logww <- as.vector(logww)-max(logww)
  #print(c('logww',logww));flush.console()
  II.new <- (runif(mm) > pp1)+1
  xx[6,] <- (2-xx[5,])*(2-II.new)+(xx[5,]-1)*(II.new-1)
  xx[5,] <- II.new
  xx[1:4,] <- xxx
  return(list(xx=xx,logww=logww))
}

#' Simulate Signals from A System with Rayleigh Flat-Fading Channels
#'
#' The function generates a sample from a system with Rayleigh flat-fading channels.
#' @param nobs sample size.
#' @param par a list with following components: \code{HH} is the state coefficient matrix; \code{WW}, \code{WW*t(WW)} is the state innovation covariance matrix;
#' \code{VV}, \code{VV*t(VV)} is the observation noise covariance matrix; \code{GG} is the observation model.
#' @examples
#' HH <- matrix(c(2.37409, -1.92936, 0.53028,0,1,0,0,0,0,1,0,0,0,0,1,0),ncol=4,byrow=TRUE)
#' WW <- matrix(c(1,0,0,0),nrow=4)
#' GG <- matrix(0.01*c(0.89409,2.68227,2.68227,0.89409),nrow=1)
#' VV <- 1.3**15*0.001
#' par <- list(HH=HH,WW=WW,GG=GG,VV=VV)
#' set.seed(1)
#' simu <- simu_fading(200,par)
#' @export
simu_fading=function(nobs,par){
  HH <- par$HH
  WW <- par$WW
  GG <- par$GG
  VV <- par$VV
  x <- rnorm(4)
  xx <- matrix(nrow=4,ncol=nobs)
  yy <- 1:nobs
  alpha <- 1:nobs
  for(i in 1:100) x <- HH%*%x+WW*rnorm(1)
  ss <- 2*floor(runif(nobs)+0.5)-1
  xx[,1] <- HH%*%x+WW*rnorm(1)
  alpha[1] <- GG%*%xx[,1]
  yy[1] <- GG%*%xx[,1]*ss[1]+VV*rnorm(1)
  for(i in 2:nobs){
    xx[,i] <- HH%*%xx[,i-1]+WW*rnorm(1)
    alpha[i] <- GG%*%xx[,i]
    yy[i] <- GG%*%xx[,i]*ss[i]+VV*rnorm(1)
  }
  return(list(xx=xx,yy=yy,ss=ss,alpha=alpha))
}




#' Fitting Univariate Autoregressive Markov Switching Models
#'
#' Fit autoregressive Markov switching models to a univariate time series using the package MSwM.
#' @param y a time series.
#' @param p AR order.
#' @param nregime the number of regimes.
#' @param include.mean a logical value for including constant terms.
#' @param sw logical values for whether coefficients are switching. The length of \code{sw} has to be equal to the number of coefficients in the model plus include.mean.
#' @return \code{MSM.fit} returns an object of class code{MSM.lm} or \code{MSM.glm}, depending on the input model.
#' @importFrom MSwM msmFit
#' @export
"MSM.fit" <- function(y,p,nregime=2,include.mean=T,sw=NULL){
  #require(MSwM)
  if(is.matrix(y))y <- y[,1]
  if(p < 0)p <- 1
  ist <- p+1
  nT <- length(y)
  X <- y[ist:nT]
  if(include.mean) X <- cbind(X,rep(1,(nT-p)))
  for (i in 1:p){
    X <- cbind(X,y[(ist-i):(nT-i)])
  }
  if(include.mean){
    colnames(X) <- c("y","cnst",paste("lag",1:p))
  }else{
    colnames(X) <- c("y",paste("lag",1:p))
  }
  npar <- ncol(X)
  X <- data.frame(X)
  mo <- lm(y~-1+.,data=X)
  ### recall that dependent variable "y" is a column in X.
  if(is.null(sw))sw = rep(TRUE,npar)
  mm <- msmFit(mo,k=nregime,sw=sw)
  mm
}



#' Sequential Importance Sampling under Clutter Environment
#'
#' This function performs one step propagation using the sequential importance sampling with full information proposal distribution under clutter environment.
#' @param mm the Monte Carlo sample size \code{m}.
#' @param xx the samples in the last iteration.
#' @param logww the log weight in the last iteration.
#' @param yyy the observations.
#' @param par a list of parameter values \code{(ssw,ssv,pd,nyy,yr)}, where \code{ssw} is the standard deviation in the state equation,
#' \code{ssv} is the standard deviation for the observation noise, \code{pd} is the probability to observe the true signal, \code{nyy} the dimension of the data,
#' and \code{yr} is the range of the data.
#' @param xdim the dimension of the state variable \code{x_t}.
#' @param ydim the dimension of the observation \code{y_t}.
#' @param resample.sch a binary vector of length \code{obs}, reflecting the resampling schedule. resample.sch[i]= 1 indicating resample should be carried out at step \code{i}.
#' @references
#' Tsay, R. and Chen, R. (2018). Nonlinear Time Series Analysis. John Wiley & Sons, New Jersey.
#' @return The function returns a list with the following components:
#' \item{xx}{the new sample.}
#' \item{logww}{the log weights.}
#' \item{r.index}{resample index, if \code{resample.sch=1}.}
#' @export
Sstep.Clutter.Full=function(mm,xx,logww,yyy,par,xdim,ydim,resample.sch){
  ssw2 <- par$ssw**2; ssv2 <- par$ssv**2
  nyy <- par$nyy; pd <- par$pd; yr <- par$yr
  mu.i <- matrix(nrow=nyy,ncol=mm)
  cc.i <- matrix(nrow=nyy,ncol=mm)
  CC <- 1:mm
  tau2 <- 1/(4/ssw2+1/ssv2); tau <- sqrt(tau2)
  aa <- 0.5*(-log(pi)+log(tau2)-log(ssw2)-log(ssv2))
  #print("where");flush.console()
  for(jj in 1:mm){
    mu.i[,jj] <- (4*(xx[1,jj]+xx[2,jj])/ssw2+yyy/ssv2)*tau2
    temp <- -2*(xx[1,jj]+xx[2,jj])**2/ssw2-yyy**2/2/ssv2
    cc.i[,jj] <- temp +mu.i[,jj]**2/2/tau2+aa
    cc.i[,jj] <- exp(cc.i[,jj])*pd/nyy
    CC[jj] <- sum(cc.i[,jj])+(1-pd)/yr
    logww[jj] <- logww[jj]+log(CC[jj])
  }
  if(resample.sch==1){
    logpp <- logww-max(logww)
    pp <- exp(logpp)
    pp <- pp/sum(pp)
    r.index <- sample.int(mm,size=mm,replace=T,prob=pp)
    xx[,] <- xx[,r.index]
    mu.i[,] <- mu.i[,r.index]
    cc.i[,] <- cc.i[,r.index]
    CC <- CC[r.index]
  }
  for(jj in 1:mm){
    #print(c(length(cc.i[jj]),nyy+1));flush.console()
    ind <- sample.int(nyy+1,1,prob=c((1-pd)/yr,cc.i[,jj])/CC[jj])
    if(ind == 1){
      ee <- rnorm(1)
      xx[1,jj] <- xx[1,jj]+xx[2,jj]+0.5*ee*par$ssw
      xx[2,jj] <- xx[2,jj]+ee*par$ssw
    }
    if(ind >1){
      xx[2,jj] <- -2*xx[1,jj]-xx[2,jj]
      xx[1,jj] <- rnorm(1)*tau+mu.i[ind-1,jj]
      xx[2,jj] <- xx[2,jj]+2*xx[1,jj]
    }
  }
  return(list(xx=xx,logww=logww,r.index=r.index))
}



#' Generic Sequential Monte Carlo Using Full Information Proposal Distribution
#'
#' Generic sequential Monte Carlo using full information proposal distribution.
#' @param SISstep.Full a function that performs one step propagation using a proposal distribution.
#' Its input includes \code{(mm,xx,logww,yyy,par,xdim,ydim,resample)}, where
#' \code{xx} and \code{logww} are the last iteration samples and log weight. \code{yyy} is the
#' observation at current time step. It should return \code{xx} (the samples xt) and
#' \code{logww} (their corresponding log weight), \code{resample} is a binary value for resampling.
#' @param nobs the number of observations \code{T}.
#' @param yy the observations with \code{T} columns and \code{ydim} rows.
#' @param mm the Monte Carlo sample size \code{m}.
#' @param par a list of parameter values to pass to \code{Sstep}.
#' @param xx.init the initial samples of \code{x_0}.
#' @param xdim the dimension of the state varible \code{x_t}.
#' @param ydim the dimension of the observation \code{y_t}.
#' @param resample.sch a binary vector of length \code{nobs}, reflecting the resampling schedule. resample.sch[i]= 1 indicating resample should be carried out at step \code{i}.
#' @param delay the maximum delay lag for delayed weighting estimation. Default is zero.
#' @param funH a user supplied function \code{h()} for estimation \code{E(h(x_t) | y_t+d}). Default
#' is identity for estimating the mean. The function should be able to take vector or matrix as input and operates on each element of the input.
#' @references
#' Tsay, R. and Chen, R. (2018). Nonlinear Time Series Analysis. John Wiley & Sons, New Jersey.
#' @return The function returns a list with the following components:
#' \item{xhat}{the fitted values.}
#' \item{loglike}{the log-likelihood.}
#' @import tensor
#' @export
SMC.Full=function(SISstep.Full,nobs,yy,mm,par,xx.init,xdim,ydim,
                  resample.sch,delay=0,funH=identity){
  #---------------------------------------------------------------
  xxd <- array(dim=c(xdim,mm,delay+1))
  delay.front = 1
  xx <- xx.init
  if(xdim==1) xx=matrix(xx,nrow=1,ncol=mm)
  xhat <- array(dim=c(xdim,nobs,delay+1))
  logww <- rep(0,mm)
  loglike <- 0
  for(i in 1:nobs){
    if(ydim==1){
      yyy <- yy[i]
    }else{
      yyy <- yy[,i]
    }
    step <- SISstep.Full(mm,xx,logww,yyy,par,xdim,ydim,resample.sch[i])
    xx <- step$xx
    logww <- step$logww-max(step$logww)
    loglike <- loglike+max(step$logww)
    if(resample.sch[i]==1){
      xxd[,,] <- xxd[,step$r.index,]
      loglike=loglike+log(sum(exp(logww)))
      logww=logww*0
    }
    xxd[,,delay.front] <- funH(xx)
    delay.front <- delay.front%%(delay+1) + 1
    order <- circular2ordinal(delay.front, delay+1)
    ww <- exp(logww)
    ww <- ww/sum(ww)
    if(xdim==1) ww=as.vector(ww)  # correction: insert
    xhat[,i,] <- tensor(xxd,ww,2,1)[,order]
  }
  if(delay > 0){
    for(dd in 1:delay){
      xhat[,1:(nobs-dd),dd+1] <- xhat[,(dd+1):nobs,dd+1]
      xhat[,(nobs-dd+1):nobs,dd+1] <- NA
    }
  }
  return(list(xhat=xhat,loglike=loglike))
}

#' Generic Sequential Monte Carlo Using Full Information Proposal Distribution and Rao-Blackwellization
#'
#' Generic sequential Monte Carlo using full information proposal distribution with Rao-Blackwellization estimate, and delay is 0.
#' @param SISstep.Full.RB a function that performs one step propagation using a proposal distribution.
#' Its input includes \code{(mm,xx,logww,yyy,par,xdim,ydim,resample)}, where
#' \code{xx} and \code{logww} are the last iteration samples and log weight. \code{yyy} is the
#' observation at current time step. It should return \code{xx} (the samples xt) and
#' \code{logww} (their corresponding log weight), \code{resample} is a binary value for resampling.
#' @param nobs the number of observations \code{T}.
#' @param yy the observations with \code{T} columns and \code{ydim} rows.
#' @param mm the Monte Carlo sample size \code{m}.
#' @param par a list of parameter values to pass to \code{Sstep}.
#' @param xx.init the initial samples of \code{x_0}.
#' @param xdim the dimension of the state varible \code{x_t}.
#' @param ydim the dimension of the observation \code{y_t}.
#' @param resample.sch a binary vector of length \code{nobs}, reflecting the resampling schedule. resample.sch[i]= 1 indicating resample should be carried out at step \code{i}.
#' @references
#' Tsay, R. and Chen, R. (2018). Nonlinear Time Series Analysis. John Wiley & Sons, New Jersey.
#' @return The function returns a list with the following components:
#' \item{xhat}{the fitted values.}
#' \item{xhatRB}{the fitted values using Rao-Blackwellization.}
#' @export

SMC.Full.RB=function(SISstep.Full.RB,nobs,yy,mm,par,xx.init,xdim,ydim,
                     resample.sch){
  #---------------------------------------------------------------
  xx <- xx.init
  if(xdim==1) xx=matrix(xx,nrow=1,ncol=mm)  ## correction: insert
  xhat <- matrix(nrow=xdim,ncol=nobs)
  xhatRB <- matrix(nrow=xdim,ncol=nobs)
  logww <- rep(0,mm)
  for(i in 1:nobs){
    if(ydim==1){
      yyy <- yy[i]
    } else{
      yyy <- yy[,i]
    }
    step <- SISstep.Full.RB(mm,xx,logww,yyy,par,xdim,ydim,resample.sch[i])
    xx <- step$xx
    xhatRB[,i] <- step$xhatRB
    xhat[,i] <- step$xhat
    logww <- step$logww-max(step$logww)
  }
  return(list(xhat=xhat,xhatRB=xhatRB))
}



#' Sequential Importance Sampling under Clutter Environment
#'
#' This function performs one step propagation using the sequential importance sampling with full information proposal distribution and returns Rao-Blackwellization estimate of mean under clutter environment.
#' @param mm the Monte Carlo sample size \code{m}.
#' @param xx the samples in the last iteration.
#' @param logww the log weight in the last iteration.
#' @param yyy the observations.
#' @param par a list of parameter values \code{(ssw,ssv,pd,nyy,yr)}, where \code{ssw} is the standard deviation in the state equation,
#' \code{ssv} is the standard deviation for the observation noise, \code{pd} is the probability to observe the true signal, \code{nyy} the dimension of the data,
#' and \code{yr} is the range of the data.
#' @param xdim the dimension of the state variable \code{x_t}.
#' @param ydim the dimension of the observation \code{y_t}.
#' @param resample.sch a binary vector of length \code{obs}, reflecting the resampling schedule. resample.sch[i]= 1 indicating resample should be carried out at step \code{i}.
#' @references
#' Tsay, R. and Chen, R. (2018). Nonlinear Time Series Analysis. John Wiley & Sons, New Jersey.
#' @return The function returns a list with the following components:
#' \item{xx}{the new sample.}
#' \item{logww}{the log weights.}
#' \item{xhat}{the fitted vlaues.}
#' \item{xhatRB}{the fitted values using Rao-Blackwellization.}
#' @export
Sstep.Clutter.Full.RB=function(mm,xx,logww,yyy,par,xdim,ydim,
                               resample.sch){
  #---------------------------------------------------------------
  ssw2 <- par$ssw**2; ssv2 <- par$ssv**2
  nyy <- par$nyy; pd <- par$pd; yr <- par$yr
  mu.i <- matrix(nrow=nyy,ncol=mm)
  cc.i <- matrix(nrow=nyy,ncol=mm)
  CC <- 1:mm
  tau2 <- 1/(4/ssw2+1/ssv2); tau <- sqrt(tau2)
  aa <- 0.5*(-log(pi)+log(tau2)-log(ssw2)-log(ssv2))
  xxRB <- 1:mm
  #print("where");flush.console()
  for(jj in 1:mm){
    mu.i[,jj] <- (4*(xx[1,jj]+xx[2,jj])/ssw2+yyy/ssv2)*tau2
    cc.i[,jj] <- -2*(xx[1,jj]+xx[2,jj])**2/ssw2-yyy**2/2/ssv2+mu.i[,jj]**2/2/tau2+aa
    cc.i[,jj] <- exp(cc.i[,jj])*pd/nyy
    CC[jj] <- sum(cc.i[,jj])+(1-pd)/yr
    logww[jj] <- logww[jj]+log(CC[jj])
    xxRB[jj] <- (sum(mu.i[,jj]*cc.i[,jj])+(xx[1,jj]+xx[2,jj])*(1-pd)/yr)/CC[jj]
  }
  logww <- logww-max(logww)
  ww <- exp(logww)
  ww <- ww/sum(ww)
  xhatRB <- sum(xxRB*ww)
  if(resample.sch==1){
    r.index <- sample.int(mm,size=mm,replace=T,prob=ww)
    xx[,] <- xx[,r.index]
    mu.i[,] <- mu.i[,r.index]
    cc.i[,] <- cc.i[,r.index]
    CC <- CC[r.index]
    logww <- logww*0
  }
  for(jj in 1:mm){
    #print(c(length(cc.i[jj]),nyy+1));flush.console()
    ind <- sample.int(nyy+1,1,prob=c((1-pd)/yr,cc.i[,jj])/CC[jj])
    if(ind == 1){
      ee <- rnorm(1)
      xx[1,jj] <- xx[1,jj]+xx[2,jj]+0.5*ee*par$ssw
      xx[2,jj] <- xx[2,jj]+ee*par$ssw
    }
    if(ind >1){
      xx[2,jj] <- -2*xx[1,jj]-xx[2,jj]
      xx[1,jj] <- rnorm(1)*tau+mu.i[ind-1,jj]
      xx[2,jj] <- xx[2,jj]+2*xx[1,jj]
    }
  }
  xhat <- sum(xx[1,]*ww)
  return(list(xx=xx,logww=logww,xhat=xhat,xhatRB=xhatRB))
}




#' Sequential Monte Carlo Using Sequential Importance Sampling for Stochastic Volatility Models
#'
#' The function implements the sequential Monte Carlo method using sequential importance sampling for stochastic volatility models.
#' @param par.natural contains three parameters in AR(1) model. The first one is the stationary mean, the second is the AR coefficient, and the third is stationary variance.
#' @param yy the data.
#' @param mm the Monte Carlo sample size.
#' @param setseed the seed number.
#' @param resample the logical value indicating for resampling.
#' @references
#' Tsay, R. and Chen, R. (2018). Nonlinear Time Series Analysis. John Wiley & Sons, New Jersey.
#' @return The function returns the log-likelihood of the data.
#' @export
wrap.SMC=function(par.natural,yy, mm, setseed=T,resample=T){
  if(setseed)	set.seed(1)
  if(resample){
    resample.sch <- rep(c(0,0,0,0,1), length=nobs)
  }
  else {
    resample.sch <- rep(0,nobs)
  }
  par=natural2ar(par.natural)
  xx.init <- par.natural[1]+par.natural[3]*rnorm(mm)
  out <- SMC(Sstep.SV, nobs, yy, mm,
             par, xx.init, 1, 1, resample.sch)
  return(out$loglike)
}




#' Setting Up The Predictor Matrix in A Neural Network for Time Series Data
#'
#' The function sets up the predictor matrix in a neural network for time series data.
#' @param zt data matrix, including the dependent variable \code{Y(t)}.
#' @param locY location of the dependent variable (column number).
#' @param nfore number of out-of-sample prediction (1-step ahead).
#' @param lags a vector containing the lagged variables used to form the x-matrix.
#' @param include.lagY indicator for including lagged \code{Y(t)} in the predictor matrix.
#' @references
#' Tsay, R. and Chen, R. (2018). Nonlinear Time Series Analysis. John Wiley & Sons, New Jersey.
#' @return The function returns a list with following components.
#' \item{X}{\code{x}-matrix for training a neural network.}
#' \item{y}{\code{y}-output for training a neural network.}
#' \item{predX}{\code{x}-matrix for the prediction subsample.}
#' \item{predY}{\code{y}-output for the prediction subsample.}
#' @export
"NNsetting" <- function(zt,locY=1,nfore=0,lags=c(1:5),include.lagY=TRUE){
  if(!is.matrix(zt))zt <- as.matrix(zt)
  if(length(lags) < 1)lags=c(1)
  ist <- max(lags)+1
  nT <- nrow(zt)
  if(nfore <= 0) nfore <- 0
  if(nfore > nT) nfore <- 0
  predX <- NULL
  predY <- NULL
  X <- NULL
  train <- nT-nfore
  y <- zt[ist:train,locY]
  if(include.lagY){z1 <- zt
  }else{
    z1 <- zt[,-locY]
  }
  ##
  for (i in 1:length(lags)){
    jj <- lags[i]
    X <- cbind(X,z1[(ist-jj):(train-jj),])
  }
  ###
  if(train < nT){
    predY <- zt[(train+1):nT,locY]
    for (i in 1:length(lags)){
      jj <- lags[i]
      predX <- cbind(predX,z1[(train+1-jj):(nT-jj),])
    }
  }


  NNsetting <- list(X=X,y=y,predX=predX,predY=predY)
}



#' Check linear models with cross validation
#'
#' The function checks linear models with cross-validation (out-of-sample prediction).
#' @param y dependent variable.
#' @param x design matrix (should include constant if it is needed).
#' @param subsize sample size of subsampling.
#' @param iter number of iterations.
#' @references
#' Tsay, R. and Chen, R. (2018). Nonlinear Time Series Analysis. John Wiley & Sons, New Jersey.
#' @return The function returns a list with following components.
#' \item{rmse}{root mean squares of forecast errors for all iterations.}
#' \item{mae}{mean absolute forecast errors for all iterations.}
#' @export
"cvlm" <- function(y,x,subsize,iter=100){
  ### Use cross-validation (out-of-sample prediction) to check linear models
  # y: dependent variable
  # x: design matrix (should inclde constant if it is needed)
  # subsize: sample size of subsampling
  # iter: number of iterations
  #
  if(!is.matrix(x))x <- as.matrix(x)
  if(is.matrix(y)) y <- c(y[,1])
  rmse <- NULL; mae <- NULL
  nT <- length(y)
  if(subsize >= nT)subsize <- nT-10
  if(subsize < ncol(x))subsize <- floor(nT/2)+1
  nfore <- nT-subsize
  ##
  for (i in 1:iter){
    idx <- sample(c(1:nT),subsize,replace=FALSE)
    y1 <- y[idx]
    x1 <- x[idx,]
    m1 <- lm(y1~-1+.,data=data.frame(x1))
    yobs <- y[-idx]
    x2 <- data.frame(x[-idx,])
    yp <- predict(m1,x2)
    err <- yobs-yp
    rmse <- c(rmse,sum(err^2)/nfore)
    mae <- c(mae,sum(abs(err))/nfore)
  }
  RMSE <- mean(rmse)
  MAE <- mean(mae)
  cat("Average RMSE: ",RMSE,"\n")
  cat("Average MAE: ",MAE,"\n")
  cvlm <- list(rmse= rmse,mae=mae)
}

