#Establishing

rm(list=ls())
# Load required libraries
library(MASS)  # for generating correlated data
library(lmtest) # for regression analysis
library(dplyr)  # for data manipulation
library(car)  #For test of multicolinearity
library(magrittr)#for data pipeline
library(glmnet) #for ridge regression
library(ggplot2)#for Visualization
library(coda)#For Markov Chain Monte Carlo

#Number of Independent Variables
n_predictors<-6

# Introduce correlations to create multicolinearity
n_observation <- c(50, 100, 500, 1000)

#rho values
rho_values<- c(0.3, 0.5, 0.6, 0.8, 0.95)

# Set the number of simulations
n_simulations <- 1000

#Setting Seed for reproducibility
set.seed(123)

#Initializing a list to store
correlated_data_list<-list() #List of Correlated Data
results <- list() #Coefficient
vifs<-list()  #Variance inflation factor

true_coeffs <- c(2.6, 3.0, 1.9, 4.2, 2.1, 1.0)  #True Coeffiecient

#Generation Multicolinear Predictor
for (size in n_observation){
  cat("sample size:",size, "\n")
for (rho in rho_values){
  # Set the number of observations and predictors
  correlation_matrix<-matrix(rho,nrow=n_predictors, ncol=n_predictors)
  diag(correlation_matrix)<-1
  # Create a data frame with the correlated variables
  correlated_data <- mvrnorm(n = size, mu = rep(0, n_predictors), Sigma = correlation_matrix)
  
  #Monte Carlo Simulation with OLS
  # Initialize a matrix to store coefficients for each simulation
  coef_mat <- matrix(0, nrow = n_predictors, ncol = n_simulations)
  for(i in n_simulations){
    error <- rnorm(n = size, mean = 0, sd = 1)
    Y <- correlated_data%*% true_coeffs + error
    
    # Combine the variables into a data frame
    data<- data.frame(Y, correlated_data)
    model<- lm(Y ~ ., data = data)
    coefficients <- coef(model)
    coef_mat[, i] <- coefficients[1:n_predictors]
  }
  
  # Store regression results in a list
  results[[as.character(rho)]]<- list(summary=summary(model))
  vifs[[as.character(rho)]]<-list(VIF=vif(model))
}

for (rho in rho_values){
  cat("average coefficient for rho=",rho,"\n")
  print(results[[as.character(rho)]])
  cat("\n")
}
for (rho in rho_values){
  cat("Variance_Inflation_Factor for rho=",rho,"\n")
  print(vifs[[as.character(rho)]])
  cat("\n")
}

}

#Monte Carlo simulation with Ridge Regression

# Vector of lambda values (Ridge regularization parameters)
lambda<-0.1

# Initialize data frame to store results
ridge_results=list()
loss_df1<-data.frame(
  Sample=numeric(),
  Rho = numeric(),
  loss=numeric()
)
results_df1 <- data.frame(
  Sample=numeric(),
  Rho = numeric(),
  Lambda = numeric(),
  Bias = numeric(),
  MSE = numeric(),
  RMSE = numeric()
)
for (size in n_observation){
for (rho in rho_values){
  # Set the number of observations and predictors
  correlation_matrix<-matrix(rho,nrow=n_predictors, ncol=n_predictors)
  diag(correlation_matrix)<-1
  # Create a data frame with the correlated variables
  X <- mvrnorm(n = size, mu = rep(0, n_predictors), Sigma = correlation_matrix)
  
  # Set the number of observations and predictors
  # Fit multiple linear regression
    # Fit a Ridge regression mode# Generate correlated data based on the correlation matrix
    result=matrix (0,nrow=n_simulations,ncol=n_predictors)
    for(i in 1: n_simulations){
      error<-rnorm(size)
      Y<- X%*%true_coeffs+error
      beta<-solve(t(X)%*%X+lambda*diag(n_predictors))%*%t(X)%*%Y
      result[i, ]<-beta
    }
    # Calculate bias, MSE, and RMSE
    ridge_mean_coef=apply(result,2,mean)
    var_coef=apply(result,2,var)
    bias=ridge_mean_coef-true_coeffs
    mse=bias^2+var_coef
    rmse=sqrt(mse)
    lossridge<-mean(mse)
    # Append results to the data frame
    results_df1<- rbind(results_df1, data.frame(Sample=size, Rho = rho, Bias =bias, MSE = mse, RMSE =rmse))
    loss_df1<-rbind(loss_df1, data.frame(Sample=size, Rho = rho, Loss=lossridge))
  }
}

print(results_df1)
print(loss_df1)



#STEIN RULE REGRESSION

# Initialize a list to store regression models and results
stein_results <- list()
loss_df2<-data.frame(
  Sample=numeric(),
  Rho = numeric(),
  loss=numeric()
)
results_df2 <- data.frame(
  Sample<-numeric(),
  Rho = numeric(),
  Bias = numeric(),
  MSE = numeric(),
  RMSE = numeric()
)
# Loop through each sample size
for (size in n_observation){
# Loop through each rho value
for (rho in rho_values) {
  # Set the number of observations and predictors
  correlation_matrix<-matrix(rho,nrow=n_predictors, ncol=n_predictors)
  diag(correlation_matrix)<-1
  X <- mvrnorm(n = size, mu = rep(0, n_predictors), Sigma = correlation_matrix)
  result=matrix (0,nrow=n_simulations,ncol=n_predictors)
  for(i in 1: n_simulations){
    error<-rnorm(size)
    Y<- X%*%true_coeffs+error
    beta_hat<-solve(t(X)%*%X)%*%t(X)%*%Y
    H<-X %*%solve(t(X)%*%X)%*%t(X)
    residuals<-Y-X%*%beta_hat
    sigma_squared<-var(residuals)
    shrinking_factor<-n_predictors*sigma_squared/sum(beta_hat^2)
    shrunk_beta<-(1-shrinking_factor)*c(beta_hat)
    result[i, ]<-shrunk_beta
  }
  # Calculate bias, MSE, and RMSE
  stein_mean_coef=apply(result,2,mean)
  var_coef=apply(result,2,var)
  bias=stein_mean_coef-true_coeffs
  mse=bias^2+var_coef
  rmse=sqrt(mse)
  losstein<-mean(mse)
  # Append results to the data frame
  results_df2 <- rbind(results_df2, data.frame(Sample=size, Rho = rho, SF=shrinking_factor, Bias =bias, MSE = mse, RMSE =rmse))
  # Print the results data frame
  loss_df2<-rbind(loss_df2, data.frame(Sample=size, Rho = rho, Loss=losstein))
}
}
print(results_df2)
print(loss_df2)


#baye rigde
# Initialize a list to store regression models and results
baye_ridge_results <- list()
loss_df3<-data.frame(
  Sample=numeric(),
  Rho = numeric(),
  loss=numeric()
)
results_df3 <- data.frame(
  Sample<-numeric(),
  Rho = numeric(),
  Bias = numeric(),
  MSE = numeric(),
  RMSE = numeric()
)

for(size in n_observation){
  # Loop through each rho value
  for (rho in rho_values) {
    # Set the number of observations and predictors
    correlation_matrix<-matrix(rho,nrow=n_predictors, ncol=n_predictors)
    diag(correlation_matrix)<-1
    X <- mvrnorm(n = size, mu = rep(0, n_predictors), Sigma = correlation_matrix)
    error<-rnorm(size)
    Y<- X%*%true_coeffs+error
    a=0.01
    b=0.01
    lambda=0.01
    #intial values:
    beta_hat<-solve(t(X)%*%X)%*%t(X)%*%Y
    residuals<-Y-X%*%beta_hat
    sigma_squared<-var(residuals)
    
    #Initialize matrix to store results:
    keep.new_sigma<-rep(0, size)
    keep.beta<-matrix(0, size, n_predictors)
    
    #precompute
    betaridge<-solve(t(X)%*%X+lambda*diag(n_predictors))%*%t(X)%*%Y
    resid1<-Y-X%*%betaridge
    new_sigma1<-var(resid1)
    
    #Start MCMC Sampler
    for(i in 1:size){
      
      # Update beta
      beta<-mvrnorm(1, betaridge, as.numeric(new_sigma1)*solve(t(X)%*%X))
      beta<-as.vector(beta)
      
      #update sigma
      new_sigma<-1/rgamma(1, size/2+a, sum((Y-X%*%betaridge)^2)/2+b)
      
      # store results:
      #Initialize matrix to store results:
      keep.new_sigma[i]<-new_sigma1
      keep.beta[i,]<-beta
      
    }
    # Calculate bias, MSE, and RMSE
    baye_ridge_mean_coef=apply(keep.beta,2,mean)
    var_coef=apply(keep.beta,2,var)
    bias=baye_ridge_mean_coef-true_coeffs
    mse=bias^2+var_coef
    rmse=sqrt(mse)
    lossbr<-mean(mse)
    # Append results to the data frame
    results_df3 <- rbind(results_df3, data.frame(Sample=size, Rho = rho, Bias =bias, MSE = mse, RMSE =rmse))
    loss_df3<-rbind(loss_df3, data.frame(Sample=size, Rho = rho, Loss=lossbr))
  }
}
print(results_df3)
print(loss_df3)


#baye stein
# Initialize a list to store regression models and results
baye_stein_results <- list()
loss_df4<-data.frame(
  Sample=numeric(),
  Rho = numeric(),
  loss=numeric()
)
results_df4 <- data.frame(
  Sample<-numeric(),
  Rho = numeric(),
  Bias = numeric(),
  MSE = numeric(),
  RMSE = numeric()
)

for(size in n_observation){
  # Loop through each rho value
  for (rho in rho_values) {
    # Set the number of observations and predictors
    correlation_matrix<-matrix(rho,nrow=n_predictors, ncol=n_predictors)
    diag(correlation_matrix)<-1
    X <- mvrnorm(n = size, mu = rep(0, n_predictors), Sigma = correlation_matrix)
    error<-rnorm(size)
    Y<- X%*%true_coeffs+error
    a=0.01
    b=0.01
  #intial values:
    beta_hat<-solve(t(X)%*%X)%*%t(X)%*%Y
    residuals<-Y-X%*%beta_hat
    sigma_squared<-var(residuals)
  
  #Initialize matrix to store results:
  keep.new_sigma<-rep(0, size)
  keep.beta<-matrix(0, size, n_predictors)
  
  #precompute
  shrinking_factor<-n_predictors*sigma_squared/sum(beta_hat^2)
  shrunk_beta<-(1-shrinking_factor)*c(beta_hat)
  resid2<-Y-X%*%as.matrix(shrunk_beta)
  new_sigma2<-var(resid2)
  
  #Start MCMC Sampler
  for(i in 1:size){
    
    # Update beta
    beta<-mvrnorm(1, as.matrix(shrunk_beta), as.numeric(new_sigma2)*solve(t(X)%*%X))
    beta<-as.vector(beta)
    
    #update sigma
    new_sigma<-1/rgamma(1, size/2+a, sum((Y-X%*%beta)^2)/2+b)
    
    # store results:
    #Initialize matrix to store results:
    keep.new_sigma[i]<-new_sigma2
    keep.beta[i,]<-beta

  }
  # Calculate bias, MSE, and RMSE
  baye_stein_mean_coef=apply(keep.beta,2,mean)
  var_coef=apply(keep.beta,2,var)
  bias=baye_stein_mean_coef-true_coeffs
  mse=bias^2+var_coef
  rmse=sqrt(mse)
  lossbs<-mean(mse)
  # Append results to the data frame
  results_df4 <- rbind(results_df4, data.frame(Sample=size, Rho = rho, SF=shrinking_factor, Bias =bias, MSE = mse, RMSE =rmse))
  loss_df4<-rbind(loss_df4, data.frame(Sample=size, Rho = rho, Loss=lossbs))
  # Print the results data frame
  }
}
print(results_df4)
print(loss_df4)

#Overall loss
all_loss<-cbind(loss_df1,loss_df2,loss_df3,loss_df4)
all_loss

