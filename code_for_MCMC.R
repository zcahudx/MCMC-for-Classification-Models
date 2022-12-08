#######################################################################
#######################################################################
#######################################################################
######                                                           ######
######            STAT0017 ICA1: Script for report output        ######
######                                                           ######
###### This script is designed to generate any numerical and     ######
###### graphical output which is gonna be used for our report    ######
###### on the Monte Carlo coursework.                            ######
######                                                           ######
#######################################################################
#######################################################################
#######################################################################

#########################################################################
# 1. Use the below code to generate the covariate data and true ?? values.
#########################################################################

library(mvtnorm)
library(coda)
library(matrixStats)

set.seed(7466) # change this to the last 4 digits of your student number

n <- 150 # number of observations
d <- 10 # number of beta parameters

# create matrix to populate with covariates
X <- matrix(nrow = n, ncol = d)
X[,1] <- rep(1, n) # first column is an intercept

# create base uncorrelated random numbers to turn into x_i's
z <- matrix(rnorm(n*(d-1)), nrow = n, ncol = d-1)

# create x_i's (ith column of matrix corresponds to variable x_i)
X[,2] <- z[,1]
X[,3] <- z[,1] + 0.2*z[,2]
X[,4] <- 0.5*z[,3]
X[,5] <- z[,4]
X[,6] <- 2*z[,5] + 20*z[,6]
X[,7] <- z[,6]
X[,8] <- 0.5 * (z[,7] + z[,4] + z[,8] + z[,1])
X[,9] <- z[,8] + 10*z[,4]
X[,10] <- z[,5] + 0.5*z[,9]

# create true beta values
beta <- seq(-2,2, length = 10)

##########################################################################
#2. Using the latent variable representation given above and the in-built 
#rt sampling function in R, draw n Yi's from the Cauchit regression model.
##########################################################################

# Note here that standard cauchy distribution is exactly the same as the 
# student t distribution with 1 degree of freedom as they have the same
# pdf. We just need to use rt for sampling from standard cauchy distribution.

# Calculate the linear term Y_hat_cauchy (a 150*1 vector)
Y_hat_cauchy<-X%*%beta

#now sample n values from the student t distribution with 1 degree of freedom
epsilon_cauchy<-rt(n, 1)

#Calculate the latent Y_star_cauchy
Y_star_cauchy<-Y_hat_cauchy+epsilon_cauchy

# Now return the n samples Y_i from Cauchit regression model
Y_cauchy<-as.numeric(Y_star_cauchy>0)

############################################################################
#3. Use rejection sampling techniques to sample from standard logistic
#distribution using standard cauchy distribution as a candidate distribution.
############################################################################

# Write functions that returns the density and CDF of standard logistic and
# standard cauchy distribution.

f_logistic<-function(x){
  return(exp(-x)/(1+exp(-x))^2)
}

f_cauchy<-function(x){
  return(1/((1+x^2)*pi))
}

F_cauchy<-function(x){
  return((1/pi)*atan(x)+(1/2))
}

F_logistic<-function(x){
  return(1/(1+exp(-x)))
}


# plot the two density functions between x in (-6, 6).
par(mfrow=c(1,1))
par(mar=c(4,4,1.5,1.5))

grid<-seq(-6,6,length.out=1000)

plot(grid, f_cauchy(grid),type="l", col="green", lwd=2, xlab="x", ylab="density",
     cex.lab=1.6)

lines(grid, f_logistic(grid),type="l",col="red",lwd=2)

legend("topleft", legend=c("Cauchy", "Logistic"),col=c("green", "red"), lwd = 2,
       cex = 1.2,bty = 'n')

dev.copy(pdf,"Cauchy_vs_logistic.pdf")
dev.off()

# Indeed, Cauchy distribution has heavy tails than logistic distributions, so the 
# value of f_logistic(x)/f_cauchy(x) will gradually decrease to 0 as x goes to 
# positive and negative infinity.

# Now we need to determine the value of M used in rejection sampling, the general
# principle is that M is larger than or equal to 1 and satisfies 
# (f_logistic(x)/f_cauchy(x))<=M. Also, Given that the probability that a sample 
# from the candidaite distribution will be accepted is 1/M, We want this M to be
# as small as possible and hence the ideal value of M should be something like
# the supremum of f_logistic(x)/f_cauchy(x) over the real line. However, here
# we are not able to find analytically the position where the function 
# f_logistic(x)/f_cauchy(x) takes its largest value via the first order derivative
# methods as f_logistic(x)/f_cauchy(x) is quite complex. What I will do is to do
# a plot of f_logistic(x)/f_cauchy(x) over (-6, 6) to approximately have an idea of
# what its maximum (or supremum) value is.

plot(grid, f_logistic(grid)/f_cauchy(grid),type="l", lwd=2, xlab="x", ylab="density",
    ylim=c(0,1.8),yaxt="n",cex.lab=1.6)

ytick<-seq(0, 1.8, by=0.1)

axis(side=2, at=ytick, labels = FALSE)

text(par("usr")[1], ytick,  
     labels = ytick, srt = 0, pos = 2, xpd = TRUE)

abline(h = 1.65, lty = 1, lwd = 2, col = "red")

dev.copy(pdf,"logistic_cauchy_ratio_plot.pdf")
dev.off()

# Indeed, the maximum value lies between 1.6 and 1.7 and hence we will set M to
# be 1.7 first and do the rejection sampling.

# Now write rejection sampling function for logistic distribution
# Note here we use dlogis and dcauchy in order to avoid numerical issues.
rejection_sampler_logistic<-function(n, M){
  
  # n is the number of samples intended for rejection sampling
  # M is the efficiency parameter of the rejection sampling algorithm
  
  samples <- rep(NA, n)
  i <- 1
  count <- 0
  
  while (i <= n){
    
    candidate<-rt(1,1)
    ratio<- dlogis(candidate)/(M*dcauchy(candidate))
    
    u <- runif(1)
    if (u < ratio) {
      samples[i] <- candidate
      i <- i + 1
    }
    count <- count + 1
  }
  
  return(list(samples = samples, count = count))
  
}



# Sample n values from logistic distribution use rejection_sampler_logistic()
output_logistic<-rejection_sampler_logistic(n, 1.7)
epsilon_logistic<-output_logistic$samples
count_M_1.7<-output_logistic$count

# divide the number of samples n by count to calculate efficiency
cat("The proportion of samples accepted when M is 1.7:", n/count_M_1.7 )
cat("The theoretical proportion of samples accepted when M is 1.7:",1/1.7)

# about 58,8% of the proposed samples are accepted, which is quite reasonable.

# Now try some different values for M, we choose 20 different values of M from
# 1.8 to 3.7, with a step of 0.1. We will then Sample n values from logistic 
# distribution use rejection_sampler_logistic() with each of the value M and 
# calculate the proportion of samples accepted in each M and create a table.
# For different M, call rejection_sampler_logistic() and calculte the efficiency
# of the rejection sampling via dividing n by the total counts in the rejection
# sampling procedure.
M_vec<-numeric(0)
for (m in seq(1.8,3.7,by=0.1)){
  Count_M_single<-rejection_sampler_logistic(n, m)$count
  M_vec<-c(M_vec,n/Count_M_single)
}

# Creat a table for efficieny of different M
tab_for_M<-data.frame(M=seq(1.8,3.7,by=0.1), efficiency=M_vec)

tab_for_M

# From the table above, we can see that the trend for the proportion of samples
# accepted is going down when M gets larger, which indeed suggests that the smaller
# the M is, the more efficient the rejection sampling is (as each proposed point
# is more likely to be accepted). Thus we will stick with M=1.7 and use the variable
# epsilon_logistic shown above as the 150 epsilon values sampled from the logistic
# distribution.

# Now we will use epsilon_logistic to creat n Y values from logistic regression
# model with the same idea as above.

# Calculate the linear term Y_hat_logistic (a 150*1 vector)
Y_hat_logistic<-X%*%beta

#Calculate the latent Y_star_logistic
Y_star_logistic<-Y_hat_cauchy+epsilon_logistic

# Now return the n samples Y_i from logistic regression model
Y_logistic<-as.numeric(Y_star_logistic>0)

##############################################################################
# 4. Use MCMC methods for generating estimates of beta value for both logistic
# and cauchit regression models. We need to first assume i.i.d standard normal 
# prior for each beta and then assume a Unit Information Prior for the beta
# vector.
##############################################################################

# Firstly, we need to write a function for the unnormalised posterior log-likelihood 
# of both the logistic and cauchit regression model under the i.i.d standard normal 
# prior.

logpi_std_normal_posterior<-function(beta, model="logistic"){
  
  # beta should be a vector of length 10
  
  # model should be either "logistic" or "cauchit" depending on the regression
  # model fitted. The default input is "logistic".
  d<-length(beta)
  log_prior<-dmvnorm(beta,mean=rep(0,d),sigma=diag(d),log = T)
  
  if (model=="logistic"){
    
    p_mat_logistic<-cbind(rep(0,150),as.vector(X%*%beta))
    log_likelihood<-sum(-apply(-p_mat_logistic[Y_logistic==1,],1,logSumExp))+
                    sum(-apply(p_mat_logistic[Y_logistic==0,],1,logSumExp))
    
    
  } else if (model=="cauchit"){
    p_vec_cauchy<-as.vector(pcauchy(X%*%beta))
    log_likelihood<-sum(dbinom(Y_cauchy,1,p_vec_cauchy,log=T))
  }
  
  return(log_likelihood+log_prior)
}

# Now write the posterior log-likelihood of both the logistic and cauchit 
# regression model under the Unit Information Prior for beta. Note here we
# have applied the logSumExp function to overcome the numerical underflow
# issue in logistic likelihood.

logpi_uni_info_posterior<-function(beta, model="logistic"){
  
  # beta should be a vector of length 10
  
  # model should be either "logistic" or "cauchit" depending on the regression
  # model fitted. The default input is "logistic".
  d<-length(beta)
  Cov_mat<-150*solve(t(X)%*%X)
  log_prior<-dmvnorm(beta,mean=rep(0,d),sigma=Cov_mat,log = T)
  
  if (model=="logistic"){
    
    p_mat_logistic<-cbind(rep(0,150),as.vector(X%*%beta))
    log_likelihood<-sum(-apply(-p_mat_logistic[Y_logistic==1,],1,logSumExp))+
                    sum(-apply(p_mat_logistic[Y_logistic==0,],1,logSumExp))
    
  } else if (model=="cauchit"){
    p_vec_cauchy<-as.vector(pcauchy(X%*%beta))
    log_likelihood<-sum(dbinom(Y_cauchy,1,p_vec_cauchy,log=T))
  }
  
  return(log_likelihood+log_prior)

}

# Having written up functions for unnormalised posterior density of beta in both
# prior, we will start to implement Metropolis-Hastings algorithms for both the 
# posterior distributions. In general, we will use Random Walk Metropolis-Hastings
# methods.

# We will first apply the vanilla Random Walk Metropolis-Hastings method to i.i.d 
# beta prior.

RWM<-function(logpi, model="logistic" ,nits, h, beta_curr, 
              V = diag(rep(1,length(beta_curr)))){
  
  # logpi is the posterior log-likelihood function of beta.
  
  # model is the regression model to be fitted, it should be either "logistic"
  # or "cauchit" and the default value is "logistic".
  
  # nits is the number of iterations.
  # h is the step size.
  # beta_curr is the starting point of beta, which is a vector of length 10.
  # V is the covariance function of the candidate multivariate normal distribution.
  
  logpi_curr <- logpi(beta_curr,model=model)
  d <- length(beta_curr)
  accepted <- 0
  beta_store <- matrix(nrow = nits, ncol = d)
  
  for (i in 1:nits) {
    # propose a candidate move
    beta_prop <- beta_curr + h*as.vector(rmvnorm(1, sigma = V))
    logpi_prop <- logpi(beta_prop,model=model)
    
    # accept-reject
    loga <- logpi_prop - logpi_curr
    u <- runif(1)
    if (log(u) < loga) {
      beta_curr <- beta_prop
      logpi_curr <- logpi_prop
      accepted <- accepted + 1
    }
    beta_store[i,] <- beta_curr
  }
  
  return(list(beta_store = beta_store, a_rate = accepted/nits))
}

# We will first try with i.i.d beta prior with Random Walk Metropolis-Hastings
# methods with cauchit regression model. First try Vanilla Random Walk
# Metropolis-Hastings with V being a 10*10 identity matrix.

# Firstly we need to carefully adjust for h so that the acceptance rate is around
# 23%. After trying a few values of h with nits being 10000, we found that h=0.12 
# seems to be a good choice, which has an acceptance rate of 24.04%, which is good
# according to the Goldilocks principle.
first_beta_std_normal_cauchit<-RWM(logpi_std_normal_posterior,
                           "cauchit",30000,0.12,seq(-2,2,length=10))

cat("The acceptance rate of the Vanilla RWM for cauchit regression with i.i.d 
    standard normal prior is:", first_beta_std_normal_cauchit$a_rate)

# We will then use coda packages for further assessment of the output
par(mfrow=c(2,5))
par(mar=c(2.5,2.5,1.5,1.5))
first_beta_std_normal_cauchit_x<-mcmc(first_beta_std_normal_cauchit$beta_store)

# We have a look at the traceplot and autocorrelation graph of the chain, the plots
# has shown that the Vanilla  RWM algorithm does not mix well at all and hence we
# need to do something else.
summary(first_beta_std_normal_cauchit_x)
traceplot(first_beta_std_normal_cauchit_x)
autocorr.plot(first_beta_std_normal_cauchit_x,auto.layout = F)

# We will try the strategy of pre-conditioning, that is, using the sampled beta
# values in the Vanilla to estimate the covariance matrix of the posterior 
# distribution of beta and use this as the input argument V of our RWM function.

Cov_std_normal_cauchit<-cov(first_beta_std_normal_cauchit$beta_store)

# Now generate the second RWM output using the pre-conditioned covariance matrix
# Cov_std_normal_logistic. We also need to tune the value of h so that the 
# acceptance rate is about 23%. h=0.1 seems to be really small here. We found that
# h=0.77 would be optimal that it generates an acceptance rate of 24.11%, which is
# great! Also, I decide increase the number of samples to be 30000.
second_beta_std_normal_cauchit<-RWM(logpi_std_normal_posterior,
                                    "cauchit",30000,0.77,seq(-2,2,length=10),
                                    Cov_std_normal_cauchit)

cat("The acceptance rate of the second pre-conditioned RWM for cauchit regression 
    with i.i.d standard normal prior is:", second_beta_std_normal_cauchit$a_rate)

# Again, use the coda package to produce the trace plots, autocorrelation graphs
# of the second sets of beta for logistic regression with i.i.d standard normal 
# prior.
second_beta_std_normal_cauchit_x<-mcmc(second_beta_std_normal_cauchit$beta_store)

par(mfrow=c(2,5))
par(mar=c(2.5,2.5,1.5,1.5))
summary(second_beta_std_normal_cauchit_x)
traceplot(second_beta_std_normal_cauchit_x)
autocorr.plot(second_beta_std_normal_cauchit_x,lag.max = 50,auto.layout = F)

# There are already huge improvements compared to the first set of beta for 
# cauchit regression with i.i.d standard normal prior. The traceplots are mixing
# much better and the autocorrelation graph also seems to be much lower than the
# first one. But still, there still exists places for improvement as the correlation 
# is still a bit high.
# So we will use the estimated covariance matrix for the second set as the input
# V for the RWM function and do a third MCMC sampling. We will use h=0.75 in the
# third MCMC sampling of beta for cauchit regression with i.i.d standard normal 
# prior. The acceptance rate is around 22.77%, which is good.

Cov_std_normal_cauchit<-cov(second_beta_std_normal_cauchit$beta_store)

third_beta_std_normal_cauchit<-RWM(logpi_std_normal_posterior,
                                     "cauchit",30000,0.75,seq(-2,2,length=10),
                                     Cov_std_normal_cauchit)

cat("The acceptance rate of the third pre-conditioned RWM for cauchit regression 
    with i.i.d standard normal prior is:", third_beta_std_normal_cauchit$a_rate)

# Still, we want to access the performance using the coda package. This time, other
# than traceplots and autocorrelation graphs. We will also use the Gelman-Rubin 
# diagnostic by sampling the MCMC from three very different starting points and
# compute R_hat.
third_beta_std_normal_cauchit_1st<-RWM(logpi_std_normal_posterior,
                                           "cauchit",30000,0.75,rep(35,10),
                                           Cov_std_normal_cauchit)

third_beta_std_normal_cauchit_2nd<-RWM(logpi_std_normal_posterior,
                                       "cauchit",30000,0.75,rep(-35,10),
                                       Cov_std_normal_cauchit)


third_beta_std_normal_cauchit_3rd<-RWM(logpi_std_normal_posterior,
                                       "cauchit",30000,0.75,rep(15,10),
                                       Cov_std_normal_cauchit)


third_beta_std_normal_cauchit_x<-mcmc(third_beta_std_normal_cauchit$beta_store)
third_beta_std_normal_cauchit_x_1st<-mcmc(third_beta_std_normal_cauchit_1st$beta_store)
third_beta_std_normal_cauchit_x_2nd<-mcmc(third_beta_std_normal_cauchit_2nd$beta_store)
third_beta_std_normal_cauchit_x_3rd<-mcmc(third_beta_std_normal_cauchit_3rd$beta_store)

# tracplots and autocorrelation plots. 
par(mfrow=c(2,5))
par(mar=c(2.5,2.5,1.5,1.5))
summary(third_beta_std_normal_cauchit_x)
traceplot(third_beta_std_normal_cauchit_x)
autocorr.plot(third_beta_std_normal_cauchit_x,lag.max = 50,auto.layout = F)

# Save the trceplot of the 1st component of the sampled beta.
par(mfrow=c(1,1))
par(mar=c(4,4,1.5,1.5))
plot(third_beta_std_normal_cauchit$beta_store[,1],type="l",xlab="Iteration",
     ylab="Beta_t",cex.lab=1.6)
dev.copy(pdf,"traceplot_norm_cauchit_1st_beta.pdf")
dev.off()

# Save the autocorrelation plot for the 1st component of sampled
# beta value.
par(mfrow=c(1,1))
par(mar=c(4,4,1.5,1.5))
plot(autocorr(third_beta_std_normal_cauchit_x[,1], lags = 1:50), type = 'b', 
     col = "red", pch = 16, ylim = c(-0.1, 1),xlab="lag",ylab="correlation coefficient",
     cex.lab=1.6)
abline(h = 0.1, lty = 3, lwd = 2, col = "black")
dev.copy(pdf, "autocorr_norm_cauchit_1st_beta.pdf")
dev.off()

# the Gelman-Rubin diagnostic for our third set of beta for cauchit regression 
# with i.i.d standard normal prior.
chains_third_beta_std_normal_cauchit<- mcmc.list(third_beta_std_normal_cauchit_x, 
                                                 third_beta_std_normal_cauchit_x_1st,
                                                 third_beta_std_normal_cauchit_x_2nd,
                                                 third_beta_std_normal_cauchit_x_3rd
                                                  )

par(mfrow=c(2,5))
par(mar=c(2.5,2.5,1.5,1.5))
gelman.diag(chains_third_beta_std_normal_cauchit)

gelman.plot(chains_third_beta_std_normal_cauchit, ylim = c(0,100),auto.layout = F)

# Saving the Gelman-Rubin diagnostic plots for 6th component of the
# beta vector.
par(mfrow=c(1,1))
par(mar=c(4,4,1.5,1.5))
gelman.plot(chains_third_beta_std_normal_cauchit[,6], ylim = c(0,100),
            auto.layout = F,cex.lab=1.6)
dev.copy(pdf,"Gelman_Rubin_norm_cauchit_6th_beta.pdf")
dev.off()

# According to Gelman-Rubin diagnostics, the chain is mixing well and the mutivariate
# R_hat is 1.01, which is less than 1.01 and the Gelman-Rubin diagnostics plots shows
# that the R_hat for all components goes to 0 at the end of 30000 samples. The R_hat
# for the 6th component goes to 0 fairly slower than others, approximately after the
# 25000th iteration, so we will discard the first 25000 samples and only calculate
# the mean value of the last 5000 samples as the estimator for beta_cauchit with standard 
# normal prior.

beta_cauchit_std_normal<-colMeans(third_beta_std_normal_cauchit$beta_store[25001:30000,])



# Having done the first set of sampling. We can apply a similar strategy for cauchit 
# regression with i.i.d standard normal prior. This time nits will just be 30000.
# Firstly, we will still apply the vanilla RWM method so that we can estimate the
# covariance matrix of the posterior distribution of beta for cauchit regression 
# with i.i.d standard normal prior. We set h=0.1 where we have an acceptance rate of
# 23.46%

first_beta_std_normal_logistic<-RWM(logpi_std_normal_posterior,
                                    "logistic",30000,0.1,seq(-2,2,length=10))

cat("The acceptance rate of the Vanilla RWM for logistic regression with i.i.d 
    standard normal prior is:", first_beta_std_normal_logistic$a_rate)

# Now estimate the posterior covariance matrix of beta using the mcmc samples
# from the first round and use this as the input V of the RWM function for the
# second round pre-conditioned RWM. h=0.78 and the acceptance rate is 25.53%.
Cov_std_normal_logistic<-cov(first_beta_std_normal_logistic$beta_store)

second_beta_std_normal_logistic<-RWM(logpi_std_normal_posterior,
                                   "logistic",30000,0.78,seq(-2,2,length=10),
                                   Cov_std_normal_logistic)



cat("The acceptance rate of the second pre-conditioned RWM for logistic regression 
    with i.i.d standard normal prior is:", second_beta_std_normal_logistic$a_rate)

# Estimate the covariance matrix using the second sampled data
Cov_std_normal_logistic<-cov(second_beta_std_normal_logistic$beta_store)

# Now do the third pre-conditioned RWM for logistic regression with i.i.d
# standard normal prior. h=0.78 and acceptance rate is 24.11%
third_beta_std_normal_logistic<-RWM(logpi_std_normal_posterior,
                                    "logistic",30000,0.78,seq(-2,2,length=10),
                                    Cov_std_normal_logistic)

cat("The acceptance rate of the third pre-conditioned RWM for logistic regression 
    with i.i.d standard normal prior is:", third_beta_std_normal_logistic$a_rate)

# Sample 3 other set of data using different start points of the RWM method for
# logistic regression with i.i.d standard normal prior and do Gelman-Rubin
# diagnostics.
third_beta_std_normal_logistic_1st<-RWM(logpi_std_normal_posterior,
                                     "logistic",30000,0.78, 
                                      rep(-45,10),
                                      Cov_std_normal_logistic)


third_beta_std_normal_logistic_2nd<-RWM(logpi_std_normal_posterior,
                                         "logistic",30000,0.78,rep(45,10),
                                         Cov_std_normal_logistic)


third_beta_std_normal_logistic_3rd<-RWM(logpi_std_normal_posterior,
                                        "logistic",30000,0.78,rep(20,10),
                                        Cov_std_normal_logistic)

third_beta_std_normal_logistic_x<-mcmc(third_beta_std_normal_logistic$beta_store)
third_beta_std_normal_logistic_x_1st<-mcmc(third_beta_std_normal_logistic_1st$beta_store)
third_beta_std_normal_logistic_x_2nd<-mcmc(third_beta_std_normal_logistic_2nd$beta_store)
third_beta_std_normal_logistic_x_3rd<-mcmc(third_beta_std_normal_logistic_3rd$beta_store)
chains_third_beta_std_normal_logistic<- mcmc.list(third_beta_std_normal_logistic_x, 
                                                 third_beta_std_normal_logistic_x_1st,
                                                 third_beta_std_normal_logistic_x_2nd,
                                                 third_beta_std_normal_logistic_x_3rd
                                                 )

# traceplots, autocorrelation plots and summary
par(mfrow=c(2,5))
par(mar=c(2.5,2.5,1.5,1.5))
summary(third_beta_std_normal_logistic_x)
traceplot(third_beta_std_normal_logistic_x)
autocorr.plot(third_beta_std_normal_logistic_x,lag.max = 50,auto.layout = F)

# Save the traceplot of the 1st component of the sampled beta.
par(mfrow=c(1,1))
par(mar=c(4,4,1.5,1.5))
plot(third_beta_std_normal_logistic$beta_store[,1],type="l",xlab="Iteration",
     ylab="Beta_t",cex.lab=1.6)
dev.copy(pdf,"traceplot_norm_logistic_1st_beta.pdf")
dev.off()

# Save the autocorrelation plot for the 1st component of sampled
# beta value.
par(mfrow=c(1,1))
par(mar=c(4,4,1.5,1.5))
plot(autocorr(third_beta_std_normal_logistic_x[,1], lags = 1:50), type = 'b', 
     col = "blue", pch = 16, ylim = c(-0.1, 1),xlab="lag",ylab="correlation coefficient",
     cex.lab=1.6)
abline(h = 0.1, lty = 3, lwd = 2, col = "black")
dev.copy(pdf, "autocorr_norm_logistic_1st_beta.pdf")
dev.off()

# Gelman-Rubin diagnostics for sampled beta
par(mfrow=c(2,5))
par(mar=c(2.5,2.5,1.5,1.5))
gelman.diag(chains_third_beta_std_normal_logistic)
gelman.plot(chains_third_beta_std_normal_logistic, ylim = c(0,100),auto.layout = F)

# save the Gelman-Rubin diagnostics plot for the 9th component of the sampled
# beta.
par(mfrow=c(1,1))
par(mar=c(4,4,1.5,1.5))
gelman.plot(chains_third_beta_std_normal_logistic[,9], ylim = c(0,100),auto.layout = F,
            ask=F,cex.lab=1.6)
dev.copy(pdf,"Gelman_Rubin_norm_logistic_9th_beta.pdf")
dev.off()

# Again, we now have pretty good traceplots and autocorrelation plot and also the 
# Gelman-Rubin diagnostics suggest that multivariate R_hat is less than 1.1, which
# means that the chain is mixing well and from Gelman-Rubin diagnostics plots.  
# From the Gelman-Rubin diagnostic plots I decided to discard the first 15000 points 
# and use the rest to calculated estimated value for beta_logistic under i.i.d standard 
# normal prior.
beta_logistic_std_normal<-colMeans(third_beta_std_normal_logistic$beta_store[15001:30000,])


# Now we will do logistic regression with unit information prior, the strategy is similar
# to what was done above.

# First Vanilla RWM with h=0.056 and acceptance rate being 24.06%
first_beta_uni_info_logistic<-RWM(logpi_uni_info_posterior,
                                    "logistic",30000,0.056,seq(-2,2,length=10))

cat("The acceptance rate of the Vanilla RWM for logistic regression with unit
    information prior is:", first_beta_uni_info_logistic$a_rate)

# estimate the posterior covariance matrix and do the second round pre-conditioned
# RWM sampling.

Cov_uni_info_logistic<-cov(first_beta_uni_info_logistic$beta_store)

# Second pre-conditioned RWM for logistic regression with unit information prior.
# h=0.89 and acceptance rate is 23.11%
second_beta_uni_info_logistic<-RWM(logpi_uni_info_posterior,
                                  "logistic",30000,0.89,seq(-2,2,length=10),
                                  Cov_uni_info_logistic)

cat("The acceptance rate of the second pre-conditioned RWM for logistic regression 
    with unit information prior is:", second_beta_uni_info_logistic$a_rate)

# Update the estimated covariance matrix
Cov_uni_info_logistic<-cov(second_beta_uni_info_logistic$beta_store)

# Now do the third pre-conditioned RWM for logistic regression with unit 
# information prior. h=0.77 and acceptance rate is 24.91%
third_beta_uni_info_logistic<-RWM(logpi_uni_info_posterior,
                                   "logistic",30000,0.77,seq(-2,2,length=10),
                                   Cov_uni_info_logistic)

cat("The acceptance rate of the third pre-conditioned RWM for logistic regression.
    with unit information prior is:", third_beta_uni_info_logistic$a_rate)

# Sample 3 RWM samples from different starting points. Again, still due to numerical
# underflow associated with logistic likelihood, I still do not choose starting
# points that are too extreme.
third_beta_uni_info_logistic_1st<-RWM(logpi_uni_info_posterior,
                                  "logistic",30000,0.77,rep(45,10),
                                  Cov_uni_info_logistic)

third_beta_uni_info_logistic_2nd<-RWM(logpi_uni_info_posterior,
                                  "logistic",30000,0.77,rep(-45,10),
                                  Cov_uni_info_logistic)

third_beta_uni_info_logistic_3rd<-RWM(logpi_uni_info_posterior,
                                      "logistic",30000,0.77,rep(20,10),
                                      Cov_uni_info_logistic)

third_beta_uni_info_logistic_x<-mcmc(third_beta_uni_info_logistic$beta_store)
third_beta_uni_info_logistic_x_1st<-mcmc(third_beta_uni_info_logistic_1st$beta_store)
third_beta_uni_info_logistic_x_2nd<-mcmc(third_beta_uni_info_logistic_2nd$beta_store)
third_beta_uni_info_logistic_x_3rd<-mcmc(third_beta_uni_info_logistic_3rd$beta_store)
chains_third_beta_uni_info_logistic<- mcmc.list(third_beta_uni_info_logistic_x, 
                                                third_beta_uni_info_logistic_x_1st,
                                                third_beta_uni_info_logistic_x_2nd,
                                                third_beta_uni_info_logistic_x_3rd
)

# traceplots, autocorrelation plots and summary
par(mfrow=c(2,5))
par(mar=c(2.5,2.5,1.5,1.5))
summary(third_beta_uni_info_logistic_x)
traceplot(third_beta_uni_info_logistic_x)
autocorr.plot(third_beta_uni_info_logistic_x,lag.max = 50,auto.layout = F)

# Save traceplot of the 1st component of the sampled beta.
par(mfrow=c(1,1))
par(mar=c(4,4,1.5,1.5))
plot(third_beta_uni_info_logistic$beta_store[,1],type="l",xlab="Iteration",
     ylab="Beta_t",cex.lab=1.6)
dev.copy(pdf,"traceplot_info_logistic_1st_beta.pdf")
dev.off()

# Save the autocorrelation plot for the 1st component of sampled
# beta value.
par(mfrow=c(1,1))
par(mar=c(4,4,1.5,1.5))
plot(autocorr(third_beta_uni_info_logistic_x[,1], lags = 1:50), type = 'b', 
     col = "brown", pch = 16, ylim = c(-0.1, 1),xlab="lag",ylab="correlation coefficient",
     cex.lab=1.6)
abline(h = 0.1, lty = 3, lwd = 2, col = "black")
dev.copy(pdf, "autocorr_info_logistic_1st_beta.pdf")
dev.off()


# Gelman-Rubin diagnostics for sampled beta.
par(mfrow=c(2,5))
par(mar=c(2.5,2.5,1.5,1.5))
gelman.diag(chains_third_beta_uni_info_logistic)
gelman.plot(chains_third_beta_uni_info_logistic, ylim = c(0,100),auto.layout = F)

# Save the Gelman-Rubin diagnostics plot for the 9th component of the 
# mcmc outcome of beta.
par(mfrow=c(1,1))
par(mar=c(4,4,1.5,1.5))
gelman.plot(chains_third_beta_uni_info_logistic[,9], ylim = c(0,100),
            auto.layout = F,cex.lab=1.6)
dev.copy(pdf,"Gelman_Rubin_info_logistic_9th_beta.pdf")
dev.off()

# Now we have a reasonable set of autocorrelation plots, traceplots and also the 
# Gelman-Rubin diagnostics R_hat that suggest really good mixing of the chain. 
# From the Gelman-Rubin diagnostics plots, it seems that the R_hat for most  
# components get below 1.1 after 20000 iterations, so we will only use the last 
# 10000 data of the third RWM output as the estimated beta_logistic with unit 
# information prior.
beta_logistic_uni_info<-colMeans(third_beta_uni_info_logistic$beta_store[20001:30000,])


# Now the last bit is the cauchit regression with unit information prior.
# Still, we will start with the vanilla RWM.

# First Vanilla RWM with h=0.055 and acceptance rate being 23.54%
first_beta_uni_info_cauchit<-RWM(logpi_uni_info_posterior,
                                  "cauchit",30000,0.055,seq(-2,2,length=10))

cat("The acceptance rate of the Vanilla RWM for cauchit regression with unit
    information prior is:", first_beta_uni_info_cauchit$a_rate)

# Use the samples from the first round RWM for estimating the covariance matrix
# of the posterior distribution of cauchit regression parameters beta with unit
# information prior.

Cov_uni_info_cauchit<-cov(first_beta_uni_info_cauchit$beta_store)

# Now do the second pre-conditioned RWM with h=0.87 and acceptance rate being
# 23.88%
second_beta_uni_info_cauchit<-RWM(logpi_uni_info_posterior,
                                 "cauchit",30000,0.87,seq(-2,2,length=10),
                                 Cov_uni_info_cauchit)

cat("The acceptance rate of the second pre-conditioned RWM for cauchit regression.
    with unit information prior is:", second_beta_uni_info_cauchit$a_rate)

# Use the sampled points in the second RWM to estimate another covariance matrix

Cov_uni_info_cauchit<-cov(second_beta_uni_info_cauchit$beta_store)

# Now do the third pre-conditioned RWM sampling for cauchit regression with unit
# information prior. h=0.76 and acceptance rate is 24.63%
third_beta_uni_info_cauchit<-RWM(logpi_uni_info_posterior,
                                  "cauchit",30000,0.76,seq(-2,2,length=10),
                                  Cov_uni_info_cauchit)

cat("The acceptance rate of the third pre-conditioned RWM for cauchit regression.
    with unit information prior is:", third_beta_uni_info_cauchit$a_rate)

# Sample 3 RWM samples from different starting points. Again, still due to numerical
# underflow associated with logistic likelihood, I still do not choose starting
# points that are too extreme.
third_beta_uni_info_cauchit_1st<-RWM(logpi_uni_info_posterior,
                                      "cauchit",30000,0.76,rep(45,10),
                                      Cov_uni_info_cauchit)

third_beta_uni_info_cauchit_2nd<-RWM(logpi_uni_info_posterior,
                                      "cauchit",30000,0.76,rep(-45,10),
                                      Cov_uni_info_cauchit)

third_beta_uni_info_cauchit_3rd<-RWM(logpi_uni_info_posterior,
                                      "cauchit",30000,0.76,rep(20,10),
                                      Cov_uni_info_cauchit)

third_beta_uni_info_cauchit_x<-mcmc(third_beta_uni_info_cauchit$beta_store)
third_beta_uni_info_cauchit_x_1st<-mcmc(third_beta_uni_info_cauchit_1st$beta_store)
third_beta_uni_info_cauchit_x_2nd<-mcmc(third_beta_uni_info_cauchit_2nd$beta_store)
third_beta_uni_info_cauchit_x_3rd<-mcmc(third_beta_uni_info_cauchit_3rd$beta_store)
chains_third_beta_uni_info_cauchit<- mcmc.list(third_beta_uni_info_cauchit_x, 
                                                third_beta_uni_info_cauchit_x_1st,
                                                third_beta_uni_info_cauchit_x_2nd,
                                                third_beta_uni_info_cauchit_x_3rd
)

# traceplots, autocorrelation plots and summary
par(mfrow=c(2,5))
par(mar=c(2.5,2.5,1.5,1.5))
traceplot(third_beta_uni_info_cauchit_x)
summary(third_beta_uni_info_cauchit_x)
autocorr.plot(third_beta_uni_info_cauchit_x,lag.max = 50,auto.layout = F)

# Save the traceplot of the 1st component of the sampled beta.
par(mfrow=c(1,1))
par(mar=c(4,4,1.5,1.5))
plot(third_beta_uni_info_cauchit$beta_store[,1],type="l",xlab="Iteration",
     ylab="Beta_t",cex.lab=1.6)
dev.copy(pdf,"traceplot_info_cauchit_1st_beta.pdf")
dev.off()

# Save the autocorrelation plot for 1st component of sampled
# beta value.
par(mfrow=c(1,1))
par(mar=c(4,4,1.5,1.5))
plot(autocorr(third_beta_uni_info_cauchit_x[,1], lags = 1:50), type = 'b', 
     col = "green", pch = 16, ylim = c(-0.1, 1),xlab="lag",ylab="correlation coefficient",
     cex.lab=1.6)
abline(h = 0.1, lty = 3, lwd = 2, col = "black")
dev.copy(pdf, "autocorr_info_cauchit_1st_beta.pdf")
dev.off()

# Gelman-Rubin diagnostics for sampled beta.
par(mfrow=c(2,5))
par(mar=c(2.5,2.5,1.5,1.5))
gelman.diag(chains_third_beta_uni_info_cauchit)
gelman.plot(chains_third_beta_uni_info_cauchit, ylim = c(0,100),auto.layout = F)

# Save the Gelman-Rubin diagnostic plot for the 9th components of the 
# mcmc outcome of beta.
par(mfrow=c(1,1))
par(mar=c(4,4,1.5,1.5))
gelman.plot(chains_third_beta_uni_info_cauchit[,9], ylim = c(0,100),
            auto.layout = F,cex.lab=1.6)
dev.copy(pdf,"Gelman_Rubin_info_cauchit_9th_beta.pdf")
dev.off()

# Again, all autocorrelation plots, traceplots and the Gelman-Rubin diagnostics 
# suggest a very good mixing of the chain. From the Gelman-Rubin diagnostics plots,  
# the R_hat of all components goes below 1.1 after 20000 iterations. So we will
# only use the last 10000 iterations for estimating the beta_cauchit with unit
# information prior.
beta_cauchit_uni_info<-colMeans(third_beta_uni_info_cauchit$beta_store[20001:30000,])


################################################################################
# 5. Fitting the logistic and cauchit regression model using the four sets of 
# beta sampled from the posterior distribution of beta under different prior.
# Note that here we simply use the mean of the mcmc samples of beta as the 
# point estimates. I use the Brier score for assessing the fit of each model.
################################################################################

# Calculate the Brier Score.
Brier_std_normal_logistic<-sum((Y_logistic-plogis(X%*%beta_logistic_std_normal))^2)/n

Brier_uni_info_logistic<-sum((Y_logistic-plogis(X%*%beta_logistic_uni_info))^2)/n

Brier_std_normal_cauchit<-sum((Y_cauchy-pcauchy(X%*%beta_cauchit_std_normal))^2)/n

Brier_uni_info_cauchit<-sum((Y_cauchy-pcauchy(X%*%beta_cauchit_uni_info))^2)/n

# Calculate log density value of true beta under both i.i.d standard 
# normal prior and unit information prior

# i.i.d standard normal
dmvnorm(beta,mean=rep(0,10),sigma=diag(10),log = T)

# unit information 
dmvnorm(beta,mean=rep(0,10),sigma=150*solve(t(X)%*%X),log = T)



