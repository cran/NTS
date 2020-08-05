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
