


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
#' VV <- 1.3**15*0.0001
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

