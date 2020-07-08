################### This file gives examples of functions for computing  
######## rank-based central sequences based on the vdW, sign and Spearman scores
######## for VAR(1) models

################### It also gives simulation examples showing 
######## how to obtain the QMLE and R-estimators 
######## under various error distributions, how to compute bias and MSE 
######## and generate their boxplots.

install.packages("clue") #Hugarian algorithm
library("clue")

install.packages("expm") #to get matrix power
library("expm")


############## Two-dimensional case
k = 2  #dimension
A = matrix(c(0.2, -0.6, 0.3, 1.1), k, k) #true parameter
n = 1000 #sample size


###coordiates of the grid points
nR = 25 #number of circles
nS = n/nR #number of lines

Y = matrix(NA, 2, n) #two dims (all points)
YY = array(NA, dim = c(2, nS, nR))
for(i in 1:nR){
  for(j in 1:nS){
    YY[1, j, i] = i/(nR + 1)*cos(2*pi*(j-1)/nS)
    YY[2, j, i] = i/(nR + 1)*sin(2*pi*(j-1)/nS)
  }
  Y[, ((i-1)*nS+1):(i*nS)] = YY[, , i]
}

######function to get the optimal assigment for the grid points (2 dims)
opt_assign = function(data, Yn){
  size = length(data[1, ])
  
  d = matrix(NA, size, size)
  for(i in 1:size){
    for(j in 1:size){
      d[i, j] = (data[1, i] - Yn[1, j])^2 + (data[2, i] - Yn[2, j])^2
    }
  }
  solve_LSAP(d) #optimal assignment of the Hungarian algorithm
}


############################## vdW R-estimation
################### in computation, we will treat \theta, \tau as matrices
######### Delta fun for the vdW
Delta.vdW = function(theta, data){
  
  size = length(data[1, ])
  
  #######residuals
  res_theta = matrix(0, k, size) 
  res_theta[, 1] = data[, 1]
  for(i in 2:size){
    res_theta[, i] = data[ ,i] - theta%*%data[ ,(i - 1)]
  }
  res_hat = res_theta
  
  ####emiprical dist., ranks and signs
  Yn = Y
  emp = matrix(0, k, size)
  emp = Yn[ , opt_assign(data = res_hat, Yn)]
  emp = as.matrix(emp)
  
  R_nR = rep(0, size)
  
  R_nR = sqrt(emp[1,]^2 + emp[2, ]^2)
  
  Sign = matrix(0, k, size)
  
  Sign[1, ]= emp[1, ]/R_nR
  Sign[2, ]= emp[2, ]/R_nR
  
  
  GAM_f = array(0, dim = c(k, k, size - 1))
  a = matrix(NA, nrow = k*k*(size - 1), ncol = 1) #vector Gamma
  
  #####compute cross-covariance matrix Gamma for the vdW score
  for(i in 1:(size - 1)){
    for(j in (i+1):size){
      GAM_f[, , i] = GAM_f[, , i] + sqrt(qchisq(R_nR[j], df = k))*sqrt(qchisq(R_nR[j-i], df = k))*Sign[,j]%*%t(Sign[,(j-i)]) 
    }
    GAM_f[, , i] = 1/(size - i)*GAM_f[, , i]
    a[((i-1)*k*k+1):(i*k*k), 1] = sqrt(size - i)*as.matrix(as.vector(GAM_f[, , i])) #vector S
  }
  
  #####compute M'P' and Qn'
  qn = matrix(NA, k, (size - 1)*k)
  MP_tran = kronecker(solve(theta), diag(k))
  
  for(i in 1: (size-1)){
    qn[, ((i - 1)*k + 1):(i*k)] = theta%^%i
  }
  
  Qn_tran = kronecker(qn, diag(k))
  
  ###central sequence
  Delta = MP_tran%*%Qn_tran%*%a
  
  return(Delta)
}


############################## sign R-estimation
################### in computation, we will treat \theta, \tau as matrices
######### Delta fun for the sign
Delta.sign = function(theta, data){
  
  size = length(data[1, ])
  
  #######residuals
  res_theta = matrix(0, k, size) 
  res_theta[, 1] = data[, 1]
  for(i in 2:size){
    res_theta[, i] = data[ ,i] - theta%*%data[ ,(i - 1)]
  }
  res_hat = res_theta
  
  ####emiprical dist., ranks and signs
  Yn = Y
  emp = matrix(0, k, size)
  emp = Yn[ , opt_assign(data = res_hat, Yn)]
  emp = as.matrix(emp)
  
  R_nR = rep(0, size)
  
  R_nR = sqrt(emp[1,]^2 + emp[2, ]^2)
  
  Sign = matrix(0, k, size)
  
  Sign[1, ]= emp[1, ]/R_nR
  Sign[2, ]= emp[2, ]/R_nR
  
  
  GAM_f = array(0, dim = c(k, k, size - 1))
  a = matrix(NA, nrow = k*k*(size - 1), ncol = 1) #vector Gamma
  
  #####compute cross-covariance matrix Gamma for the sign score
  for(i in 1:(size - 1)){
    for(j in (i+1):size){
      GAM_f[, , i] = GAM_f[, , i] + Sign[,j]%*%t(Sign[,(j-i)]) 
    }
    GAM_f[, , i] = 1/(size - i)*GAM_f[, , i]
    a[((i-1)*k*k+1):(i*k*k), 1] = sqrt(size - i)*as.matrix(as.vector(GAM_f[, , i])) #vector S
  }
  
  #####compute M'P' and Qn'
  qn = matrix(NA, k, (size - 1)*k)
  MP_tran = kronecker(solve(theta), diag(k))
  
  for(i in 1: (size-1)){
    qn[, ((i - 1)*k + 1):(i*k)] = theta%^%i
  }
  
  Qn_tran = kronecker(qn, diag(k))
  
  ###central sequence
  Delta = MP_tran%*%Qn_tran%*%a
  
  return(Delta)
}


############################## Spearman R-estimation
################### in computation, we will treat \theta, \tau as matrices
######### Delta fun for the Spearman
Delta.Sp = function(theta, data){
  
  size = length(data[1, ])
  
  #######residuals
  res_theta = matrix(0, k, size) 
  res_theta[, 1] = data[, 1]
  for(i in 2:size){
    res_theta[, i] = data[ ,i] - theta%*%data[ ,(i - 1)]
  }
  res_hat = res_theta
  
  ####emiprical dist., ranks and signs
  Yn = Y
  emp = matrix(0, k, size)
  emp = Yn[ , opt_assign(data = res_hat, Yn)]
  emp = as.matrix(emp)
  
  R_nR = rep(0, size)
  
  R_nR = sqrt(emp[1,]^2 + emp[2, ]^2)
  
  Sign = matrix(0, k, size)
  
  Sign[1, ]= emp[1, ]/R_nR
  Sign[2, ]= emp[2, ]/R_nR
  
  
  GAM_f = array(0, dim = c(k, k, size - 1))
  a = matrix(NA, nrow = k*k*(size - 1), ncol = 1) #vector Gamma
  
  #####compute cross-covariance matrix Gamma for the Spearman R-score
  for(i in 1:(size - 1)){
    for(j in (i+1):size){
      GAM_f[, , i] = GAM_f[, , i] + emp[, j]%*%t(emp[, (j - i)]) 
    }
    GAM_f[, , i] = 1/(size - i)*GAM_f[, , i]
    a[((i-1)*k*k+1):(i*k*k), 1] = sqrt(size - i)*as.matrix(as.vector(GAM_f[, , i])) #vector S
  }
  
  #####compute M'P' and Qn'
  qn = matrix(NA, k, (size - 1)*k)
  MP_tran = kronecker(solve(theta), diag(k))
  
  for(i in 1: (size-1)){
    qn[, ((i - 1)*k + 1):(i*k)] = theta%^%i
  }
  
  Qn_tran = kronecker(qn, diag(k))
  
  ###central sequence
  Delta = MP_tran%*%Qn_tran%*%a
  
  return(Delta)
}



################## Simulation under various error distributions

install.packages("mvtnorm") #for generating multivariate normal dist.
library("mvtnorm")

install.packages("MTS") #to get the QMLE for the VARMA 
library("MTS")

install.packages("clue") #Hungarian algorithm
library("clue")

install.packages("expm") #to get matrix power
library("expm")

install.packages("sn") # generate multivariate skew normal and skew-t distributions
library(sn)


k = 2  #dimension
A = matrix(c(0.2, -0.6, 0.3, 1.1), k, k) #true parameter
n = 1000 #sample size


########### generate VAR(1) model
nn = 500 #burn-in sample size
N = n + nn

R = 300 # number of replications
XXt = array(0, dim = c(k, N, R))
Xt = array(0, dim = c(k, n, R))


##################generate a large sample to compute the mean of skew-normal dist.
#a = rmst(100000, xi=rep(0, k), Omega = matrix(c(7,4,4,5), k, k), alpha = c(5, 2), nu=Inf)
#mean_skewnorm = apply(a, 2, mean)
#mat_mean_skewnorm = matrix(NA, N, k)
#for(i in 1:N){
#  mat_mean_skewnorm[i, ] = mean_skewnorm
#}
#rm(a)


##################generate a large sample to compute the mean of skew-t3 dist.
#a = rmst(100000, xi=rep(0, k), Omega = matrix(c(7,4,4,5), k, k), alpha = c(5, 2), nu=3)
#mean_skewt3 = apply(a, 2, mean)
#mat_mean_skewt3 = matrix(NA, N, k)
#for(i in 1:N){
#  mat_mean_skewt3[i, ] = mean_skewt3
#}
#rm(a)



for(r in 1:R){
  ###normal dist
  err = rmvnorm(N, mean = rep(0, k), sigma = diag(k))
  
  ###t(3) dist
  #err = rmvt(N, sigma = diag(k), df = 3)
  
  ###skew normal dist
  #err = rmst(N, xi=rep(0, k), Omega = matrix(c(7,4,4,5), k, k), alpha = c(5, 2), nu=Inf)
  #err = err - mat_mean_skewnorm
  
  ###skew t(3) dist
  #err = rmst(N, xi=rep(0, k), Omega = matrix(c(7,4,4,5), k, k), alpha = c(5, 2), nu=3)
  #err = err - mat_mean_skewt3
  
  ###mixture of 3 normal dists
  #err1 = rmvnorm(N, mean = c(-5, 0), sigma = matrix(c(7, 5, 5, 5), 2, 2))
  #err2 = rmvnorm(N, mean = c(5, 0), sigma = matrix(c(7, -6, -6, 6), 2, 2))
  #err3 = rmvnorm(N, mean = c(0, 0), sigma = matrix(c(4, 0, 0, 3), 2, 2))
  #u = runif(N)
  #err = (u<3/8)*err1 + ((u>=3/8)&(u<3/4))*err2 + (u>=3/4)*err3
  
  for(i in 2:N){
    XXt[, i, r] = A%*%XXt[, (i - 1), r] + err[i, ]
  }
}
Xt = XXt[, (nn + 1): N, ] 
rm(XXt)


###########################additive outliers
#prop = 0.05  #rate of contamination
#num.addit = prop*n  #number of contamination
#interval.len = n/num.addit #interval length between two contaminations 
#size = 4 #size of contamination
#addit.size = rep(c(rep(0, interval.len - 1), 1), num.addit)*size
#for(r in 1:R){
#  for(j in 1:k){
#    Xt[j, 1:n, r] = Xt[j, 1:n, r] + addit.size
#    Xt[j, 1:n, r] = Xt[j, 1:n, r] - size*prop #substract from mean
#  }
#}


################################ QMLE
#res = array(NA, dim = c(k, n, R))
A_QMLE = array(NA, dim = c(k, k, R))
for(r in 1:R){
  fit_QMLE = VARMA(t(Xt[, ,r]), p = 1, q=0, include.mean = F)
  A_QMLE[, , r] = t(fit_QMLE$coef)
  #res[, , r] = t(fit_QMLE$residuals)
}



###################### compute vdW R-estimator
###################### we need to estimate the (negative) cross information matrix $-\Upsilon$
######function to create canonical  basis
make_basis = function(place, dimen = k^2) replace(numeric(dimen), place, 1)

Upsilon = array(NA, dim = c(k^2, k^2, R)) 
Delta1 = matrix(NA, k^2, R)
Delta2 = matrix(NA, k^2, R)
change_vdW = matrix(NA, k^2, R) #change in each iteration of the one-step procedure
theta_vdW = array(NA, dim = c(k, k, R)) #R-estimate
###set initial values
theta_vdW = A_QMLE #initial value
}

################for each replications, compute vdW R-estimator
for(i in 1:R){
  Delta1[, i] = Delta.vdW(A_QMLE[, , i], Xt[, , i])
  
  for(l in 1:k^2){
    btau = make_basis(l)
    Delta2[ ,i] = Delta.vdW((A_QMLE[, , i]+ 1/sqrt(n)*matrix(btau, 2, 2)), Xt[, , i])
    Upsilon[, l, i] = (Delta2[ ,i] - Delta1[ ,i])
  }
  
  solve_Upsilon = solve(Upsilon[, , i])
  
  #################iterations for the one-step procedure
  for(j in 1:5){
    change_vdW[ , i] = -1/sqrt(n)*solve_Upsilon%*%Delta.vdW(theta_vdW[ , , i], Xt[, , i])
    theta_vdW[ , , i] = theta_vdW[ , , i] + matrix(change_vdW[ , i], 2, 2)
  }
}


########################compute sign R-estimator
################# we need to estimate the (negative) cross information matrix $-\Upsilon$
###function to create canonical basis
make_basis = function(place, dimen = k^2) replace(numeric(dimen), place, 1)

Upsilon = array(NA, dim = c(k^2, k^2, R)) 
Delta1 = matrix(NA, k^2, R)
Delta2 = matrix(NA, k^2, R)
change_sign = matrix(NA, k^2, R)

theta_sign = array(NA, dim = c(k, k, R))

theta_sign = A_QMLE #initial estimator


###############compute sign R-estimator for each replication
for(i in 1:R){
  Delta1[, i] = Delta.sign(A_QMLE[, , i], Xt[, , i])
  
  for(l in 1:k^2){
    btau = make_basis(l)
    Delta2[ ,i] = Delta.sign((A_QMLE[, , i]+ 1/sqrt(n)*matrix(btau, 2, 2)), Xt[, , i])
    Upsilon[, l, i] = (Delta2[ ,i] - Delta1[ ,i])
  }
  
  solve_Upsilon = solve(Upsilon[, , i])
  
  #################iterations for the one-step procedure
  for(j in 1:5){
    change_sign[ , i] = -1/sqrt(n)*solve_Upsilon%*%Delta.sign(theta_sign[ , , i], Xt[, , i])
    theta_sign[ , , i] = theta_sign[ , , i] + matrix(change_sign[ , i], 2, 2)
  }
}


########################compute Spearman R-estimator
#################  we need to estimate the (negative) cross information matrix $-\Upsilon$
###function to create canonical basis
make_basis = function(place, dimen = k^2) replace(numeric(dimen), place, 1)

Upsilon = array(NA, dim = c(k^2, k^2, R)) 
Delta1 = matrix(NA, k^2, R)
Delta2 = matrix(NA, k^2, R)
change_Sp = matrix(NA, k^2, R)

theta_Sp = array(NA, dim = c(k, k, R))

theta_Sp = A_QMLE #initial estimator
}

#######compute Spearman R-estimator for each replication
for(i in 1:R){
  Delta1[, i] = Delta.Sp(A_QMLE[, , i], Xt[, , i])
  
  for(l in 1:k^2){
    btau = make_basis(l)
    Delta2[ ,i] = Delta.Sp((A_QMLE[, , i]+ 1/sqrt(n)*matrix(btau, 2, 2)), Xt[, , i])
    Upsilon[, l, i] = (Delta2[ ,i] - Delta1[ ,i])
  }
  
  solve_Upsilon = solve(Upsilon[, , i])
  
  for(j in 1:5){
    change_Sp[ , i] = -1/sqrt(n)*solve_Upsilon%*%Delta.Sp(theta_Sp[ , , i], Xt[, , i])
    theta_Sp[ , , i] = theta_Sp[ , , i] + matrix(change_Sp[ , i], 2, 2)
  }
}


################ Once we have obtained the QMLE and R-estimators,
################ we do boxplot and compute the bias and MSE
QMLE = A_QMLE
vdW = theta_vdW
Sign = theta_sign
Spear = theta_Sp

#### compute bias and MSE
bias_QMLE = MSE_QMLE = bias_vdW =MSE_vdW = bias_sign = MSE_sign = bias_Sp = MSE_Sp = rep(NA, 4)
for(i in 1:2){
  bias_QMLE[i] = mean(QMLE[i,1,]) - A[i,1]
  MSE_QMLE[i] = mean((QMLE[i,1,] - A[i,1])^2)
  
  bias_vdW[i] = mean(vdW[i,1,]) - A[i,1]
  MSE_vdW[i] = mean((vdW[i,1,] - A[i,1])^2)
  
  bias_sign[i] = mean(Sign[i,1,]) - A[i,1]
  MSE_sign[i] = mean((Sign[i,1,] - A[i,1])^2)
  
  bias_Sp[i] = mean(Spear[i,1,]) - A[i,1]
  MSE_Sp[i] = mean((Spear[i,1,] - A[i,1])^2)
}


for(i in 1:2){
  bias_QMLE[i+2] = mean(QMLE[i,2,]) - A[i,2]
  MSE_QMLE[i+2] = mean((QMLE[i,2,] - A[i,2])^2)
  
  bias_vdW[i+2] = mean(vdW[i,2,]) - A[i,2]
  MSE_vdW[i+2] = mean((vdW[i,2,] - A[i,2])^2)
  
  bias_sign[i+2] = mean(Sign[i,2,]) - A[i,2]
  MSE_sign[i+2] = mean((Sign[i,2,] - A[i,2])^2)
  
  bias_Sp[i+2] = mean(Spear[i,2,]) - A[i,2]
  MSE_Sp[i+2] = mean((Spear[i,2,] - A[i,2])^2)
}


##################boxplots

#####get values for boxplot
names=c(rep("QMLE", R) , rep("vdW", R), rep("Sign", R), rep("Spearman", R))
value=c(QMLE[1, 1,],  vdW[1, 1, ], Sign[1 ,1, ], Spear[1, 1, ])
data=data.frame(names,value)


####Draw the boxplot, with the number of individuals per group
#boxplot for $a_11$
par(mfcol = c(2,2))
a=boxplot(data$value~data$names, main = "a11", lwd = 0.8)
abline(h = 0.2, col = "red")

#boxplot for $a_21$
value=c(QMLE[2, 1,], vdW[2, 1,], Sign[2, 1,], Spear[2, 1,])
data=data.frame(names,value)

a=boxplot(data$value~data$names, main = "a21", lwd = 0.8)
abline(h = -0.6, col = "red")

#boxplot for $a_12$
value=c(QMLE[1, 2,], vdW[1, 2, ], Sign[1, 2, ], Spear[1, 2, ])
data=data.frame(names,value)

a=boxplot(data$value~data$names, main = "a12", lwd = 0.8)
abline(h = 0.3, col = "red")

#boxplot for $a_22$
value=c(QMLE[2, 2,], vdW[2, 2, ], Sign[2, 2, ], Spear[2, 2, ])
data=data.frame(names,value)

a=boxplot(data$value~data$names, main = "a22", lwd = 0.8)
abline(h = 1.1, col = "red")




############################ Examples of R-estimation functions for 
################# three-dimensional case are given as follows

k = 3  # dimension
n = 1000

############################# R-estimation
########################First, using mvmesh package to create grid
install.packages("mvmesh")
library(mvmesh)

nR = 15
nS = 66
n0 = n - nR*nS
grid_unitSphere = t(UnitSphere(n = k, k = 2)$V)

Y = matrix(NA, k, n) #k dims (all points)
for(j in 1:nR){
  Y[, ((j-1)*nS+1):(j*nS)] = j/(nR + 1)*grid_unitSphere
}
Y[, (nR*nS+1):n] = matrix(0, k, n0)


###function to get the optimal grid points (3 dims)
opt_assign = function(data, Yn){
  size = length(data[1, ])
  
  d = matrix(NA, size, size)
  for(i in 1:size){
    for(j in 1:size){
      d[i, j] = sum((data[, i] - Yn[, j])^2)
    }
  }
  solve_LSAP(d) #optimal assignment of the Hungarian algorithm
}


############################## vdW R-estimation
################### in computation, we will treat \theta, \tau as matrices
######### Delta fun for the vdW
Delta.vdW = function(theta, data){
  #theta = matrix(theta, ncol = k)
  
  size = length(data[1, ])
  
  #######residuals
  res_theta = matrix(0, k, size) 
  res_theta[, 1] = data[, 1]
  for(i in 2:size){
    res_theta[, i] = data[ ,i] - theta%*%data[ ,(i - 1)]
  }
  res_hat = res_theta
  
  ####emiprical dist., ranks and signs
  Yn = Y
  emp = matrix(0, k, size)
  emp = Yn[ , opt_assign(data = res_hat, Yn)]
  emp = as.matrix(emp)
  
  R_nR = rep(0, size)  ##Rt/(nR + 1)
  
  for(i in 1:size){
    R_nR[i] = sqrt(sum((emp[ ,i])^2))
  }
  
  Sign = matrix(0, k, size)
  
  i = which(R_nR == 0)
  index = 1:size
  index = index[-i]
  
  for(i in 1:k){
    Sign[i, index]= emp[i, index]/R_nR[index]
  }
  
  GAM_f = array(0, dim = c(k, k, size - 1))
  a = matrix(NA, nrow = k*k*(size - 1), ncol = 1) #vector Gamma
  
  
  #####compute cross-covariance matrix Gamma and vector S for three types of R-scores
  for(i in 1:(size - 1)){
    for(j in (i+1):size){
      GAM_f[, , i] = GAM_f[, , i] + sqrt(qchisq(R_nR[j], df = k))*sqrt(qchisq(R_nR[j-i], df = k))*Sign[,j]%*%t(Sign[,(j-i)]) 
    }
    GAM_f[, , i] = 1/(size - i)*GAM_f[, , i]
    a[((i-1)*k*k+1):(i*k*k), 1] = sqrt(size - i)*as.matrix(as.vector(GAM_f[, , i])) #vector S
  }
  
  #####compute M'P' and Qn'
  qn = matrix(NA, k, (size - 1)*k)
  MP_tran = kronecker(solve(theta), diag(k))
  
  for(i in 1: (size-1)){
    qn[, ((i - 1)*k + 1):(i*k)] = theta%^%i
  }
  
  Qn_tran = kronecker(qn, diag(k))
  
  ###central sequence
  Delta = MP_tran%*%Qn_tran%*%a
  
  return(Delta)
}


######### Delta fun for the sign
Delta.sign = function(theta, data){
  #theta = matrix(theta, ncol = k)
  
  size = length(data[1, ])
  
  #######residuals
  res_theta = matrix(0, k, size) 
  res_theta[, 1] = data[, 1]
  for(i in 2:size){
    res_theta[, i] = data[ ,i] - theta%*%data[ ,(i - 1)]
  }
  res_hat = res_theta
  
  ####emiprical dist., ranks and signs
  Yn = Y
  emp = matrix(0, k, size)
  emp = Yn[ , opt_assign(data = res_hat, Yn)]
  emp = as.matrix(emp)
  
  R_nR = rep(0, size)  ##Rt/(nR + 1)
  
  for(i in 1:size){
    R_nR[i] = sqrt(sum((emp[ ,i])^2))
  }
  
  Sign = matrix(0, k, size)
  
  i = which(R_nR == 0)
  index = 1:size
  index = index[-i]
  
  for(i in 1:k){
    Sign[i, index]= emp[i, index]/R_nR[index]
  }
  
  
  GAM_f = array(0, dim = c(k, k, size - 1))
  a = matrix(NA, nrow = k*k*(size - 1), ncol = 1) #vector Gamma
  
  #####compute cross-covariance matrix Gamma and vector S for three types of R-scores
  for(i in 1:(size - 1)){
    for(j in (i+1):size){
      GAM_f[, , i] = GAM_f[, , i] + Sign[,j]%*%t(Sign[,(j-i)]) 
    }
    GAM_f[, , i] = 1/(size - i)*GAM_f[, , i]
    a[((i-1)*k*k+1):(i*k*k), 1] = sqrt(size - i)*as.matrix(as.vector(GAM_f[, , i])) 
  }
  
  #####compute M'P' and Qn'
  qn = matrix(NA, k, (size - 1)*k)
  MP_tran = kronecker(solve(theta), diag(k))
  
  for(i in 1: (size-1)){
    qn[, ((i - 1)*k + 1):(i*k)] = theta%^%i
  }
  
  Qn_tran = kronecker(qn, diag(k))
  
  ###central sequence
  Delta = MP_tran%*%Qn_tran%*%a
  
  return(Delta)
}


######### Delta fun for the Spearman
Delta.Sp = function(theta, data){
  #theta = matrix(theta, ncol = k)
  
  size = length(data[1, ])
  
  #######residuals
  res_theta = matrix(0, k, size) 
  res_theta[, 1] = data[, 1]
  for(i in 2:size){
    res_theta[, i] = data[ ,i] - theta%*%data[ ,(i - 1)]
  }
  res_hat = res_theta
  
  ####emiprical dist., ranks and signs
  Yn = Y
  emp = matrix(0, k, size)
  emp = Yn[ , opt_assign(data = res_hat, Yn)]
  emp = as.matrix(emp)
  
  R_nR = rep(0, size)  ##Rt/(nR + 1)
  
  for(i in 1:size){
    R_nR[i] = sqrt(sum((emp[ ,i])^2))
  }
  
  Sign = matrix(0, k, size)
  
  i = which(R_nR == 0)
  index = 1:size
  index = index[-i]
  
  for(i in 1:k){
    Sign[i, index]= emp[i, index]/R_nR[index]
  }
  
  GAM_f = array(0, dim = c(k, k, size - 1))
  a = matrix(NA, nrow = k*k*(size - 1), ncol = 1) #vector Gamma
  
  #####compute cross-covariance matrix Gamma and vector S for three types of R-scores
  for(i in 1:(size - 1)){
    for(j in (i+1):size){
      GAM_f[, , i] = GAM_f[, , i] + emp[, j]%*%t(emp[, (j - i)]) 
    }
    GAM_f[, , i] = 1/(size - i)*GAM_f[, , i]
    a[((i-1)*k*k+1):(i*k*k), 1] = sqrt(size - i)*as.matrix(as.vector(GAM_f[, , i])) #vector S
  }
  
  #####compute M'P' and Qn'
  qn = matrix(NA, k, (size - 1)*k)
  MP_tran = kronecker(solve(theta), diag(k))
  
  for(i in 1: (size-1)){
    qn[, ((i - 1)*k + 1):(i*k)] = theta%^%i
  }
  
  Qn_tran = kronecker(qn, diag(k))
  
  ###central sequence
  Delta = MP_tran%*%Qn_tran%*%a
  
  return(Delta)
}

