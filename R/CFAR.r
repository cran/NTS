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
#' y=g_cfar1(100,5,phi_func,grid=1000,sigma=1,ini=100)
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
#' y=g_cfar2(100,5,phi_func1,phi_func2,grid=1000,sigma=1,ini=100)
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
#' @param p the CFAR order.
#' @param df_b the degrees of freedom for natural cubic splines. Default is 10.
#' @param grid the number of gird points used to construct the functional time series and noise process. Default is 1000.
#' @references
#' Liu, X., Xiao, H., and Chen, R. (2016) Convolutional autoregressive models for functional time series. \emph{Journal of Econometrics}, 194, 263-282.
#' @return The function returns a list with components:
#' \item{phi_coef}{the estimated spline coefficients for convolutional function values, a (2*grid+1)-by-p matrix.}
#' \item{phi_func}{the estimated convolutional function(s), a (df_b+1)-by-p matrix.}
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

  pdf=function(para)
  {
    psi= para^(1/n);	### correlation between two adjacent observation point exp(-para/n)
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
    l=-(t-p)*(n+1)/2*log(ehat)+1/2*log.mat
    return(-l);
  }
  para= exp(-5)
  result=optim(para,pdf, lower=0.0001, upper=0.9999, method='L-BFGS-B')
  para=result$par		####exp(-rho)
  psi=para^(1/n)		####exp(-rho/n)
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
  rho_hat= -log(para)
  sigma_hat=sqrt(ehat/((t-p)*(n+1))/(1-psi^2))
  phi_coef=t(matrix(phihat,coef_b,p))
  est_cfar <- list(phi_coef=phi_coef,phi_func=phihat_func,rho=rho_hat,sigma=sigma_hat)
  return(est_cfar)
}


#' F Test for a CFAR Process
#'
#' F test for a CFAR process to specify the CFAR order.
#' @param f the functional time series.
#' @param p.max the maximum CFAR order. Default is 6.
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
  statistic= (sse0-sse1)/coef_b/sse1*((n+1)*(t-1)-coef_b)
  pval=1-pf(statistic,coef_b,(t-1)*(n+1)-coef_b)
  cat("Test and p-value of Order 0 vs Order 1: ","\n")
  print(c(statistic,pval))

  indep.p=indep
  for(p in 2:p.max){
    indep.pre=indep.p[(n+2):((n+1)*(t-p+1)),]
    indep.tmp=NULL
    for(i in 1:p){
      tmp=indep[((i-1)*(n+1)+1):((n+1)*(t-p-1+i)),]
      indep.tmp=cbind(tmp,indep.tmp)
    }
    indep.p=indep.tmp
    pdf=function(para)
    {
      mat1= matrix(0, coef_b*p, coef_b*p);
      mat2= matrix(0, coef_b*p, 1);
      psi= para^(1/n);
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
    para= exp(-5)
    result=optim(para,pdf, lower=0.0001, upper=0.9999, method='L-BFGS-B')

    mat1= matrix(0, p*coef_b, p*coef_b);
    mat2= matrix(0, coef_b*p, 1);
    psi_invall=matrix(0,(n+1)*(t-1),(n+1))
    para=result$par
    psi=para^(1/n)
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



    mat1= matrix(0, (p-1)*coef_b, (p-1)*coef_b);
    mat2= matrix(0, coef_b*(p-1), 1);
    psi_invall=matrix(0,(n+1)*(t-1),(n+1))
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
      M_t= indep.pre[(i-p-1)*(n+1)+(1:tmp.n),]
      tmp_mat= t(M_t) %*% psi_matinv;
      mat1= mat1+ tmp_mat %*% M_t;
      mat2= mat2+ tmp_mat %*% f[i,1:tmp.n];
    }
    phihat= solve(mat1)%*%mat2
    predict= t(matrix(indep.pre%*%phihat,ncol=t-p,nrow=n+1))
    sse.pre= sum(diag((f[(p+1):t,]-predict)%*%phi_matinv%*% t(f[(p+1):t,]-predict)))


    statistic= (sse.pre-sse.p)/coef_b/sse.p*((n+1)*(t-p)-p*coef_b)
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
    n=dim(f)[2]-1
  }else{
    n=max(num_obs)-1
  }

  if(length(num_obs)!=t){
    num_obs=rep(n+1,t)
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
  index=matrix(0,t,1+n)
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
  para4= exp(-5)
  result4=optim(para4,pdf4, lower=0.0001, upper=0.9999, method='L-BFGS-B')
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
  sigma_hat=sqrt(ehat/sum(num_obs[2:t]))

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
  sigma_hat=sqrt(ehat/sum(num_obs[2:t]))
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
      index=matrix(0,t,1+n)
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
      para4= exp(-5)
      result4=optim(para4,pdf4, lower=0.0001, upper=0.9999, method='L-BFGS-B')
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
      sigma_hat=sqrt(ehat/sum(num_obs[(p+1):t]))

      sigma_hat1=sigma_hat
      rho_hat1=rho_hat
      predict= t(matrix(indep.p%*%phihat,ncol=t-p,nrow=n+1))
      ssr.p=0
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
      statistic= (sse.pre-sse.p)/coef_b/sse.p*(sum(num_obs[(p+1):t])-p*coef_b)
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
#' @param p the CFAR order.
#' @param grid the number of gird points used to construct the functional time series and noise process. Default is 1000.
#' @param df_b the degrees of freedom for natural cubic splines. Default is 10.
#' @param num_obs the numbers of observations. It is a t-by-1 vector, where t is the length of time.
#' @param x_pos the observation location matrix. If the locations are regular, it is a t-by-(n+1) matrix with all entries 1/n.
#' @references
#' Liu, X., Xiao, H., and Chen, R. (2016) Convolutional autoregressive models for functional time series. \emph{Journal of Econometrics}, 194, 263-282.
#' @return The function returns a list with components:
#' \item{phi_coef}{the estimated spline coefficients for convolutional function(s).}
#' \item{phi_func}{the estimated convolutional function(s).}
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
    n=dim(f)[2]-1
  }else{
    n=max(num_obs)-1
  }

  if(length(num_obs)!=t){
    num_obs=rep(n+1,t)
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
  index=matrix(0,t,n+1)
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
    para4= exp(-5)
    result4=optim(para4,pdf4, lower=0.0001, upper=0.9999, method='L-BFGS-B')
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
    sigma_hat=sqrt(ehat/sum(num_obs[2:t]))
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
    para4= exp(-5)
    result4=optim(para4,pdf4, lower=0.0001, upper=0.9999, method='L-BFGS-B')
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
#' index=seq(1,1001,by=50)
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
