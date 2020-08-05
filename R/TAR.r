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
#' Estimation of a univariate two-regime SETAR model, including threshold value, performing recursive least squares method or nested sub-sample search algorithm.
#' The procedure of Li and Tong (2016) is used to search for the threshold.
#' @param y a vector of time series.
#' @param p1,p2 AR-orders of regime 1 and regime 2.
#' @param d delay for threshold variable, default is 1.
#' @param thrV threshold variable. If thrV is not null, it must have the same length as that of y.
#' @param thrQ lower and upper quantiles to search for threshold value.
#' @param Trim lower and upper quantiles for possible threshold values.
#' @param include.mean a logical value indicating whether constant terms are included.
#' @param method "RLS": estimate the model by conditional least squares method implemented by recursive least squares; "NeSS": estimate the model by conditional least squares method implemented by Nested sub-sample search (NeSS) algorithm.
#' @param k0 the maximum number of threshold values to be evaluated, when the nested sub-sample search (NeSS) method is used. If the sample size is large (> 3000), then k0 = floor(nT*0.5). The default is k0=300. But k0 = floor(nT*0.8) if nT < 300.
#' @references
#' Li, D., and Tong. H. (2016) Nested sub-sample search algorithm for estimation of threshold models. \emph{Statisitca Sinica}, 1543-1554.
#' @return uTAR returns a list with components:
#' \item{data}{the data matrix, y.}
#' \item{arorder}{AR orders of regimes 1 and 2.}
#' \item{delay}{the delay for threshold variable.}
#' \item{residuals}{estimated innovations.}
#' \item{sresi}{standardized residuals.}
#' \item{coef}{a 2-by-(p+1) matrices. The first row shows the estimation results in regime 1, and the second row shows these in regime 2.}
#' \item{sigma}{estimated innovational covariance matrices of regimes 1 and 2.}
#' \item{nobs}{numbers of observations in regimes 1 and 2.}
#' \item{model1,model2}{estimated models of regimes 1 and 2.}
#' \item{thr}{threshold value.}
#' \item{D}{a set of threshold values.}
#' \item{RSS}{RSS}
#' \item{AIC}{AIC value}
#' \item{cnst}{logical values indicating whether the constant terms are included in regimes 1 and 2.}
#' @examples
#' phi=t(matrix(c(-0.3, 0.5,0.6,-0.3),2,2))
#' y=uTAR.sim(nob=2000, arorder=c(2,2), phi=phi, d=2, thr=0.2, cnst=c(1,-1),sigma=c(1, 1))$series
#' est=uTAR(y=y,p1=2,p2=2,d=2,thrQ=c(0,1),Trim=c(0.1,0.9),include.mean=TRUE,method="NeSS",k0=50)
#' @export
"uTAR" <- function(y, p1, p2, d=1, thrV=NULL, thrQ=c(0,1), Trim=c(0.1,0.9), include.mean=TRUE, method="RLS", k0=300){

  if(is.matrix(y))y=y[,1]
  p = max(p1,p2)
  if(p < 1)p=1
  ist=max(p,d)+1
  nT <- length(y)
  ### built in checking
  nT1 <- nT
  if(!is.null(thrV))nT1 <- length(thrV)
  if(nT != nT1){
    cat("Input error at thrV. Reset to a SETAR model","\n")
    thrV <- y
  }
  yt <- y[ist:nT]
  nobe <- nT-ist+1
  Phi <- matrix(0,2,p+1)
  x1 <- NULL
  for (i in 1:p){
    x1=cbind(x1,y[(ist-i):(nT-i)])
  }
  x1 <- as.matrix(x1)
  if(length(thrV) < nT){
    tV <- y[(ist-d):(nT-d)]
  }else{
    tV <- thrV[(ist-d):(nT-d)]
  }
  if(method=="NeSS"){
    if(nT > 3000)k0 <- floor(nT*0.5)
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
      cat("nob1 & sigma1:",c(n1, sigma1), "\n")
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
      cat("nob2 & sigma2: ",c(n2,sigma2),"\n")
      fitted.values <- yt-resi
      Sigma <- c(sigma1,sigma2)
      ### pooled estimate
      sigmasq = (n1*sigma1^2+n2*sigma2^2)/(n1+n2)
      AIC <- AIC1+AIC2
      cat(" ","\n")
      cat("Overall MLE of sigma: ",sqrt(sigmasq),"\n")
      cat("Overall AIC: ",AIC,"\n")
    }else{cat("No threshold found: Try again with a larger k0.","\n")}
  }
  if(method=="RLS"){
    if(include.mean){x1 <- cbind(rep(1,nobe),x1)}
    k1 <- ncol(x1)
    x2 <- NULL
    for (i in 1:p2){
      x2 <- cbind(x2,y[(ist-i):(nT-i)])
    }
    x2 <- as.matrix(x2)
    if(include.mean){x2 <- cbind(rep(1,nobe),x2)}
    k2 <- ncol(x2)
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

    sig1 = coef1$sigma^2*(n1-ncol(X))/n1
    AIC1 = n1*log(sig1)+2*ncol(X)
    sig2 = coef2$sigma^2*(n2-ncol(X))/n2
    AIC2 = n2*log(sig2)+2*ncol(X)
    AIC=AIC1+AIC2
    sigmasq=(n1*coef1$sigma^2+n2*coef2$sigma^2)/(n1+n2)
    cat(" ","\n")
    cat("Overall MLE of sigma: ",sqrt(sigmasq),"\n")
    cat("Overall AIC: ",AIC,"\n")
    Size=c(n1,n2)
    Sigma=c(coef1$sigma,coef2$sigma)
  }
  uTAR <- list(data=y,arorder=c(p1,p2),delay=d,residuals = resi, sresi=sresi, coefs = Phi, sigma=Sigma,
                 nobs = Size, model1 = m1a, model2 = m1b,thr=thr, D=D, RSS=RSS,AIC = AIC,
                 cnst = rep(include.mean,2))
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
#' \item{k}{the number of regimes.}
#' \item{arorder}{AR orders of regimes 1 and 2.}
#' \item{coefs}{a k-by-(p+1) matrices, where \code{k} is the number of regimes. The i-th row shows the estimation results in regime i.}
#' \item{sigma}{estimated innovational covariances for all the regimes.}
#' \item{thr}{threshold value.}
#' \item{residuals}{estimated innovations.}
#' \item{sresi}{standardized residuals.}
#' \item{nobs}{numbers of observations in different regimes.}
#' \item{delay}{delay for threshold variable.}
#' \item{cnst}{logical values indicating whether the constant terms are included in different regimes.}
#' \item{AIC}{AIC value.}
#' @examples
#' phi=t(matrix(c(-0.3, 0.5,0.6,-0.3),2,2))
#' y=uTAR.sim(nob=200, arorder=c(2,2), phi=phi, d=2, thr=0.2, cnst=c(1,-1),sigma=c(1, 1))
#' thr.est=uTAR(y=y$series, p1=2, p2=2, d=2, thrQ=c(0,1),Trim=c(0.1,0.9), method="RLS")
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
    if(output) cat("Overall AIC: ",AIC,"\n")
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
#' y=uTAR.sim(nob=2000, arorder=c(2,2), phi=phi, d=2, thr=0.2, cnst=c(1,-1), sigma=c(1, 1))
#' thr.est=uTAR(y=y$series, p1=2, p2=2, d=2, thrQ=c(0,1), Trim=c(0.1,0.9), method="RLS")
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

#' Generate Two-Regime (TAR) Models
#'
#' Generates multivariate two-regime threshold autoregressive models.
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
#' \item{series}{a time series following the multivariate two-regime  VAR model.}
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

