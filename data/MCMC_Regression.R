library("stein.thinning")
library(mcmc)
library(rstan)


mydata <- read.csv("loans_good.csv")

regression_stein_thinning <- function() {
  # Linear Regression Model
  mc <- "
    data {
      int<lower=0> N;
      int<lower=0> K;
      matrix[N, K] x;
      vector[N] y;
    }
    transformed data {
      matrix[N, K] Q_ast;
      matrix[K, K] R_ast;
      matrix[K, K] R_ast_inverse;
      Q_ast = qr_thin_Q(x) * sqrt(N - 1);
      R_ast = qr_thin_R(x) / sqrt(N - 1);
      R_ast_inverse = inverse(R_ast);
    }
    parameters {
      real alpha;
      vector[K] theta;
      real<lower=0> sigma;
    }
    model {
      y ~ normal(Q_ast * theta + alpha, sigma);
    }
    generated quantities {
      vector[K] beta;
      beta = R_ast_inverse * theta;
    }
    "
  
  # Standardize x and y in initial data
  # find a way to put an initial N(0,1) on beta in the model code
  
  data_list <- list(N = nrow(mydata),
                    K = 1,
                    x = as.matrix(mydata[, 2]),
                    y = mydata[, 1])
  
  
  fit <- stan(model_code=mc, data=data_list, iter=1000000, chains=1)
  
  
  # Extract sampled points and gradients
  smp <- extract(fit, permuted=FALSE, inc_warmup=FALSE)
  smp <- smp[,,1:3]
  scr <- t(apply(smp, 1, function(x) rstan::grad_log_prob(fit, x)))
  
  # Obtain a subset of 40 points
  idx <- thin(smp, scr, 50)
  stein_data <- smp[idx,1:3]
  return(stein_data)
}

final_data <- regression_stein_thinning()
alpha <- final_data[,1]
beta <- final_data[,2]
sigma <- final_data[,3]

hist(beta)
hist(alpha)
mu <- mean(beta)
mu2 <- mean(alpha)
v <- var(beta)
shapiro_test <- shapiro.test(beta)
mu
v
shapiro_test
