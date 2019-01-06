library(survival)
library(Rcpp)
library(parallel)
library(Brobdingnag)


rm(list = ls())
setwd("D:/Amir/Study/SurvivalAnalysis/")
# setwd("/a/home/cc/students/math/nirkeret/frailty_LTRC/")
sourceCpp("frailty_noLT_age_weights_brob.cpp", verbose=TRUE)

n <- 1000


r <- 0.61 #type1 ceonsoring - administrative censoring
timestocheck <- seq(0,r,0.01) #time grid to examine the estimated cumulative hazard functions

Hres12 <- Hres13 <- Hres23 <- par_res <- NULL  #arrays to save the results

g12 <- c(2,0.2,0.05)  #regression effects
g13 <- c(0.05,0.04)
g23 <- c(1,0.15)
 
recL <- 0.05   #piecewise constant hazards change points
recU <- 0.15

theta <- 1  #frailty distribution parameter
c12_1 <- 0.005 #constant hazard12 below recL
c12_2 <- 1 #constant hazard12 between recL and recU
c12_3 <- 1 #constant hazard12 above recU
c13_1 <- 0.5 #constant hazard13 below recL
c13_2 <- 1 #constant hazard13 between recL and recU
c13_3 <- 5 #constant hazard13 above recU
c23_1 <- 0  #constant hazard23 below recL
c23_2 <- 1 #constant hazard23 between recL and recU
c23_3 <- 1 #constant hazard23 above recU


reltol <- 1e-04  #relative tolerance stopping criterion


H012_real <- function(t) 
{
  ifelse(t < recL,c12_1*t,ifelse(t < recU,c12_1*recL + c12_2*(t - recL),c12_1*recL + c12_2*(recU - recL) + c12_3*(t - recU)))
}
H013_real <- function(t)
{
  ifelse(t < recL,c13_1*t,ifelse(t < recU,c13_1*recL + c13_2*(t - recL),c13_1*recL + c13_2*(recU - recL) + c13_3*(t - recU)))
}

H023_real <- function(t)
{
  ifelse(t < recL,c23_1*t,ifelse(t < recU,c23_1*recL + c23_2*(t - recL),c23_1*recL + c23_2*(recU - recL) + c23_3*(t - recU)))
}

#laplace transform and its derivatives for the gamma distribution

phi <- function(theta,s) (1+theta*s)^(-1/theta)
phi_tag <- function(theta,s) -(1+theta*s)^(-(theta + 1) / theta)
phi_tag2 <- function(theta,s) (theta+1) * (1+theta*s)^(-(2*theta + 1)/theta)


#getting A() at min(V,tau). tau = maximal observed failure time
get_A_col <- function(g,z,theta_hat,H,times,sortedTimes) 
{
	if(length(g) == 1) gz_mat <- z %*% t(rep(g,1))
	else gz_mat <- z %*% matrix(rep(g,1),nrow=length(g))
	a <- theta_hat*(exp(gz_mat))
	b <- H(sortedTimes)
	ret <- sum_e_sqr(a,b,times,sortedTimes)
	return(exp(gz_mat)*ret)
}

#getting ln of A at min(V,tau), for handling large values
get_ln_A_col <- function(g,z,theta_hat,H,times,sortedTimes) 
{
  if(length(g) == 1) gz_mat <- z %*% t(rep(g,1))
  else gz_mat <- z %*% matrix(rep(g,1),nrow=length(g))
  a <- theta_hat*(exp(gz_mat))
  b <- H(sortedTimes)
  ret <- ln_sum_e_sqr(a,b,times,sortedTimes)
  return(gz_mat+ret)
}

# getting ln of alpha star 
get_lnastar <- function(g,z,theta_hat,H,times,index_lst) {
	Htimes = H(times)[index_lst[,2]]
  	if(length(g) == 1) gz_mat = z[index_lst[,1]] %*% t(as.vector(g))
	else gz_mat = z[index_lst[,1],] %*% as.vector(g)
	return(gz_mat + theta_hat*exp(gz_mat)*Htimes)
}

#derivative of ln(alpha_star) w.r.t g
deri_lnastar_by_g <- function(g,z,theta_hat,H,times,index_lst) {
  Htimes = H(times)[index_lst[,2]]
  if(length(g) == 1) tmp_z = z[index_lst[,1]]
  else tmp_z = z[index_lst[,1],]
  if(length(g) == 1) gz_mat = tmp_z %*% t(as.vector(g))
  else gz_mat = tmp_z %*% as.vector(g)
  return(tmp_z + tmp_z*matrix(rep(theta_hat*exp(gz_mat)*Htimes,length(g)),ncol=length(g)))
}


deri_lnastar_by_theta <- function(g,z,theta_hat,H,times,index_lst) {
  Htimes = H(times)[index_lst[,2]]
  if(length(g) == 1) gz_mat = z[index_lst[,1]] %*% t(as.vector(g))
  else gz_mat = z[index_lst[,1],] %*% as.vector(g)
  return(exp(gz_mat)*Htimes)	
}


deri_A_by_g <- function(dvA,Ai,theta_hat,z_mat,N_t=N_tau,delta=NULL) {
  if(is.null(delta)) {delta <- rep(1,length(Ai))}
  dA = as.matrix(dvA[,-c(1,2)])
  ret_mat <- NULL
  ret <-  (N_t == 0)* -1*(1/(1 + theta_hat*Ai)) + (N_t == 1) * -1*((theta_hat + 1)/(1 + theta_hat*Ai)) +(N_t == 2) * -1*((2*theta_hat + 1)/(1 + theta_hat*Ai))
  for(j in 1:ncol(dA))
  {
    ret_mat <- cbind(ret_mat,as.numeric(sign(z_mat[,j])*delta*ret*as.brob(exp(1))^dA[,j]))
  }
  return(ret_mat)
}



deri_A_by_theta <- function(ds,s,theta_hat,N_t=N_tau) {
  ret = (N_t == 0)* (
    (log(1 + theta_hat*s)/(theta_hat**2)) - (theta_hat*ds + s)/(theta_hat*(1 + theta_hat*s))
  )
  ret = ret + (N_t == 1)* (
    (-(theta_hat+1)*theta_hat^2*ds-theta_hat*s*(-log(1 + theta_hat*s) + theta_hat + 1) + log(1 + theta_hat*s))/
      ((theta_hat**2)*(1 + theta_hat*s))
  )  
  ret = ret + (N_t == 2)* (
    -(2*theta_hat + 1)*(theta_hat*ds + s)/(theta_hat*(1 + theta_hat*s)) -
      (-(2*theta_hat +1)/(theta_hat**2) + 2/theta_hat)*log(1 + theta_hat*s) + 
      1/(theta_hat+1)
  )               
  return(ret)
}


Likelihood <- function(par)
{
  theta_hat <- par[1]
  g12_gscore_hat <- par[2]
  g12_escore_hat <- par[3]
  g12_family_hat <- par[4]
  
  g12_hat <- c(g12_gscore_hat,g12_escore_hat,g12_family_hat)
  
  g13_gscore_hat <- par[5]
  g13_escore_hat <- par[6]
  
  g13_hat <- c(g13_gscore_hat,g13_escore_hat)
  
  g23_gscore_hat <- par[7]
  g23_V_hat <- par[8]
  
  g23_hat <- c(g23_gscore_hat,g23_V_hat)
  
  z12 <- z[,1:3] ; z13 <- z[,1:2] ; z23 <- z[,c(1,4)]
  
  lnastar12 <- get_lnastar(g12_hat,z12,theta_hat,H012_estimated,sorted_ev_12, cbind(which(delta1),rank(V[delta1]))) 
  lnastar13 <- get_lnastar(g13_hat,z13,theta_hat,H013_estimated,sorted_ev_13, cbind(which(delta2),rank(V[delta2]))) 
  lnastar23 <- get_lnastar(g23_hat,z23,theta_hat,H023_estimated,sorted_ev_23, cbind(which(delta3),rank(W[delta3]))) 
  
  A12_V <- as.brob(exp(1)) ^ as.vector(get_ln_A_col(g12_hat,z12,theta_hat,H012_estimated,V,sorted_ev_12))
  A13_V <- as.brob(exp(1)) ^ as.vector(get_ln_A_col(g13_hat,z13,theta_hat,H013_estimated,V,sorted_ev_13))
  A23_V <- as.brob(exp(1)) ^ as.vector(get_ln_A_col(g23_hat,z23,theta_hat,H023_estimated,V,sorted_ev_23))
  A23_W <- as.brob(exp(1)) ^ as.vector(get_ln_A_col(g23_hat,z23,theta_hat,H023_estimated,W,sorted_ev_23))
  
  A_all_num <- (A12_V+A13_V+A23_W-delta1*A23_V) 
  
  phi_num_0 <- log(phi(theta_hat,A_all_num[N_tau == 0])) * weightsN0
  phi_num_1 <- log(-(phi_tag(theta_hat,A_all_num[N_tau == 1]))) * weightsN1
  phi_num_2 <- log(phi_tag2(theta_hat,A_all_num[N_tau == 2])) * weightsN2
  
  ret <- sum(lnastar12,lnastar13,lnastar23,phi_num_0,phi_num_1,phi_num_2) 
  
  return(-(ret))  #return the minus of the log-likelihood so that the optim will actually find the maximum
}


#function to get the derivatives of the log-likelihood w.r.t theta, gamma
Likelihood_grad <- function(par)
{
  theta_hat <- par[1]
  g12_gscore_hat <- par[2]
  g12_escore_hat <- par[3]
  g12_family_hat <- par[4]
  
  g12_hat <- c(g12_gscore_hat,g12_escore_hat,g12_family_hat)
  
  g13_gscore_hat <- par[5]
  g13_escore_hat <- par[6]
  
  g13_hat <- c(g13_gscore_hat,g13_escore_hat)
  
  g23_gscore_hat <- par[7]
  g23_V_hat <- par[8]
  
  g23_hat <- c(g23_gscore_hat,g23_V_hat)
  
  z12 <- z[,1:3] ; z13 <- z[,1:2] ; z23 <- z[,c(1,4)]
  
  by_g12 =                  colSums(weights12 * deri_lnastar_by_g(g12_hat,z12,theta_hat,H012_estimated,sorted_ev_12, cbind(which(delta1),rank(V[delta1]))))
  by_theta =            sum(weights12 * deri_lnastar_by_theta(g12_hat,z12,theta_hat,H012_estimated,sorted_ev_12, cbind(which(delta1),rank(V[delta1]))))
  by_g13 =                  colSums(weights13 * deri_lnastar_by_g(g13_hat,z13,theta_hat,H013_estimated,sorted_ev_13, cbind(which(delta2),rank(V[delta2]))))
  by_theta = by_theta + sum(weights13 * deri_lnastar_by_theta(g13_hat,z13,theta_hat,H013_estimated,sorted_ev_13, cbind(which(delta2),rank(V[delta2]))))
  by_g23 =                  colSums(weights23 * deri_lnastar_by_g(g23_hat,z23,theta_hat,H023_estimated,sorted_ev_23, cbind(which(delta3),rank(W[delta3]))))
  by_theta = by_theta + sum(weights23 * deri_lnastar_by_theta(g23_hat,z23,theta_hat,H023_estimated,sorted_ev_23, cbind(which(delta3),rank(W[delta3]))))

  #the function derive A combined returns the ln of the derivative of A w.r.t theta and g. the first column return is A itself
  A12_derive = derive_A_combined(theta_hat, as.matrix(z12), g12_hat, H012_estimated(sorted_ev_12), V, sorted_ev_12)
  A13_derive = derive_A_combined(theta_hat, as.matrix(z13), g13_hat, H013_estimated(sorted_ev_13), V, sorted_ev_13) 
  A23_derive_W = derive_A_combined(theta_hat, as.matrix(z23), g23_hat, H023_estimated(sorted_ev_23), W, sorted_ev_23) 
  A23_derive_V = derive_A_combined(theta_hat, as.matrix(z23), g23_hat, H023_estimated(sorted_ev_23), V, sorted_ev_23) 

  A_all <- as.brob(exp(1))^A12_derive[,1] + as.brob(exp(1))^A13_derive[,1] + as.brob(exp(1))^A23_derive_W[,1]	- delta1*as.brob(exp(1))^(A23_derive_V[,1])
  
  by_g12 = by_g12 + colSums(weights * deri_A_by_g(A12_derive,A_all,theta_hat,z12))
  by_g13 = by_g13 + colSums(weights * deri_A_by_g(A13_derive,A_all,theta_hat,z13))
  by_g23 = by_g23 + colSums(weights * deri_A_by_g(A23_derive_W,A_all,theta_hat,z23))
  by_g23 = by_g23 - colSums(weights * deri_A_by_g(A23_derive_V,A_all,theta_hat,z23,delta = delta1))
  by_theta = by_theta + sum(weights * deri_A_by_theta(as.brob(exp(1))^A12_derive[,2] + as.brob(exp(1))^A13_derive[,2] + as.brob(exp(1))^A23_derive_W[,2] - delta1*as.brob(exp(1))^A23_derive_V[,2],A_all,theta_hat))
  
  ret = -1*c(as.numeric(by_theta), by_g12, by_g13, by_g23)

  return(ret)
}


########################
seed <- 200734754
sim=1
##############################################
for (sim in 1:100) {

  start <- Sys.time()
  set.seed(sim+seed)
  
  ####sampling begins####
  
  z1 <- runif(n) ; z2 <- runif(n) ; z3 <- runif(n) 
  z <- cbind(z1,z2,z3)
  
  
  omega <- rgamma(n,shape = 1/theta,scale = theta)  #sampling the frailty variates
  u12 <- runif(n)
  u13 <- runif(n)
  u23 <- runif(n)
  #sampling the failure times:
  q12 <- log(1 - theta*log(u12)/omega) / (theta * as.vector(exp(z[,1:3] %*% g12)))
  T12 <- ifelse(q12 < c12_1*recL,q12/c12_1,ifelse(q12 < c12_1*recL + c12_2*(recU-recL),(q12 + (c12_2-c12_1)*recL)/c12_2,(q12 + (c12_2 - c12_1)*recL + (c12_3 - c12_2)*recU)/c12_3))
  q13 <- log(1 - theta*log(u13)/omega) / (theta * as.vector(exp(z[,1:2] %*% g13)))
  T13 <- ifelse(q13 < c13_1*recL,q13/c13_1,ifelse(q13 < c13_1*recL + c13_2*(recU-recL),(q13 + (c13_2-c13_1)*recL)/c13_2,(q13 + (c13_2 - c13_1)*recL + (c13_3 - c13_2)*recU)/c13_3))
  
  #sampling the failure time for 23 from the truncated distribution
  z <- cbind(z,T12)
  u23_trunc <- u23 * exp(-omega * (exp(theta * exp(z[,c(1,4)] %*% g23)*H023_real(T12)) - 1) / theta)
  q23 <- log(1 - theta*log(u23_trunc)/omega) / (theta * as.vector(exp(z[,c(1,4)] %*% g23)))
  T23 <- ifelse(q23 < c23_1*recL,q23/c23_1,ifelse(q23 < c23_1*recL + c23_2*(recU-recL),(q23 + (c23_2-c23_1)*recL)/c23_2,(q23 + (c23_2 - c23_1)*recL + (c23_3 - c23_2)*recU)/c23_3))

  C <- rexp(n,2)  #independent censoring times
  V <- pmin(T12,T13,C,r)  #observed times
  delta1 <- V == T12 
  delta2 <- V == T13
  W <- ifelse(!delta1,0,pmin(T23,C,r))  #death-after-disease times
  delta3 <- as.vector((W == T23) & delta1)
  
  N_1_tau <- delta1 + delta2   
  N_tau <- delta1 + delta2 + delta3
  
  
  ####sampling ends####

  n_ev_12 <- sum(delta1)
  n_ev_13 <- sum(delta2)
  n_ev_23 <- sum(delta3)
  
  sorted_ev_12 <- sort(V[delta1])
  sorted_ev_13 <- sort(V[delta2])
  sorted_ev_23 <- sort(W[delta3])
  
  #weights:
  weights <- rep(1,n)  #if no other weights should be used
  weights12 <- weights[delta1] ; weights13 <- weights[delta2] ; weights23 <- weights[delta3]
  weightsN0 <- weights[N_tau == 0] ; weightsN1 <- weights[N_tau == 1] ; weightsN2 <- weights[N_tau == 2] 
  weights_events_12 <- weights[match(sorted_ev_12,V)]
  weights_events_13 <- weights[match(sorted_ev_13,V)]
  weights_events_23 <- weights[match(sorted_ev_23,W)]

  #initialization
  
  #truth:
g12_hat <- g12 ; g13_hat <- g13; g23_hat <- g23; theta_hat <- theta


  H012_estimated <- H012_real
  H013_estimated <- H013_real
  H023_estimated <- H023_real
  
  #naive cox models:
  # model12_0 <- coxph(Surv(V,event = delta1) ~ z[,1:3])
  # model13_0 <- coxph(Surv(V,event = delta2) ~ z[,1:2])
  # model23_0 <- coxph(Surv(W[delta1],event = delta3[delta1]) ~ z[delta1,1])
  # g12_hat <- coef(model12_0)
  # g13_hat <- coef(model13_0)
  # g23_hat <- coef(model23_0)
  # theta_hat <- 0.01

  #naive Breslow estimators
  
  # baseHazCox12 <- basehaz(model12_0,centered = F)
  # H012_estimated <- stepfun(baseHazCox12$time,c(0,baseHazCox12$hazard))
  # baseHazCox13 <- basehaz(model13_0,centered = F)
  # H013_estimated <- stepfun(baseHazCox13$time,c(0,baseHazCox13$hazard))
  # baseHazCox23 <- basehaz(model23_0,centered = F)
  # H023_estimated <- stepfun(baseHazCox23$time,c(0,baseHazCox23$hazard))
  iter <- 0
  conv_cond <- T

  while(conv_cond) {
    iter <- iter+1

    #estimating the hazard functions jumps
    dH12 = estimate_H_JK(z[,1:3] %*% g12_hat, z[,1:2] %*% g13_hat, theta_hat, sorted_ev_12, c(0,sorted_ev_13), H012_estimated(c(0,sorted_ev_12)), H013_estimated(c(0,sorted_ev_13)), V , N_1_tau,weights_subjects = weights,weights_events = weights_events_12)
    dH13 = estimate_H_JK(z[,1:2] %*% g13_hat, z[,1:3] %*% g12_hat, theta_hat, sorted_ev_13, c(0,sorted_ev_12), H013_estimated(c(0,sorted_ev_13)), H012_estimated(c(0,sorted_ev_12)) , V , N_1_tau,weights_subjects = weights,weights_events = weights_events_13)
    A_1_tau = get_A_col(g12_hat,z[,1:3],theta_hat,H012_estimated,V,sorted_ev_12) + get_A_col(g13_hat,z[,1:2],theta_hat,H013_estimated,V,sorted_ev_13) - get_A_col(g23_hat,z[,c(1,4)],theta_hat,H023_estimated,V,sorted_ev_23)
    dH23 = estimate_H_23(z[,c(1,4)] %*% g23_hat, theta_hat, sorted_ev_23, H023_estimated(c(0,sorted_ev_23)), W,V, delta3, N_1_tau, A_1_tau,weights_subjects = weights,weights_events = weights_events_23)


    last_loglike <- Likelihood(c(theta_hat,g12_hat,g13_hat,g23_hat))


    H012_estimated <- stepfun(sorted_ev_12,c(0,cumsum(dH12)))
    H013_estimated <- stepfun(sorted_ev_13,c(0,cumsum(dH13)))
    H023_estimated <- stepfun(sorted_ev_23,c(0,cumsum(dH23)))

    error_switch <- T
    while(error_switch)
    {
    res <- try(optim(par=c(theta_hat,g12_hat,g13_hat,g23_hat),Likelihood,method ="L-BFGS-B",lower = (c(0.00001,rep(-Inf,7))), gr = Likelihood_grad)$par)
      if(class(res) != "try-error")
      {
        error_switch <- F
      } else   #if an error was returned, probably due to high values, we decrease theta and g23[2]. It should not happen often because of Brobdingnag.
      {
        theta_hat <- 0.9*theta_hat ; g23_hat[2] <- 0.9*g23_hat[2]
      }
    }

    print(c(sim,res))
    
    #updating the estimates
    theta_hat <- res[1]
    g12_hat <- res[2:4]
    g13_hat <- res[5:6]
    g23_hat <- res[7:8]
    cat("iteration",iter)

    #calculating the convergence criterion
    new_loglike <- Likelihood(c(theta_hat,g12_hat,g13_hat,g23_hat))
    loglikelihood_rel_change <- abs(new_loglike - last_loglike)/abs(last_loglike)
    print("relative change in loglikelihood is:")
    print(loglikelihood_rel_change)
    print(paste0('Sim is: ', sim))
    conv_cond <- loglikelihood_rel_change > reltol

  }

  Hres12 <- rbind(Hres12,H012_estimated(timestocheck))
  Hres13 <- rbind(Hres13,H013_estimated(timestocheck))
  Hres23 <- rbind(Hres23,H023_estimated(timestocheck))

  par_res <- rbind(par_res,c(theta_hat,g12_hat,g13_hat,g23_hat))
}

###############


