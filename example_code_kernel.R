
library(quadprog)
library(numDeriv)
library(MASS)

### simulation settings
# sample size
n_seq <- seq(100, 1000, by = 100)

source("source_functions_kernel.R")


# component of inferential interest
comp_con <- c(1, 3)

# sparse irregular measurements
L_seq <- 5:15

# measurement time distribution
beta_para <- c(1, 1.25)

# noise level
rho_err <- 0.1
sigma_err <- 0.5

# bandwidth
h <- 0.15
h_con <- 0.15

# regression coefficient functions
p <- 4
rho_x <- 0.5

beta1 <- function(t) ifelse(t < 0.5, 2*sin((6/12)*2*pi*t)-1, cos((6/12)*2*pi*(t - 0.5))) 
beta2 <- function(t) ifelse(t < 0.6, 2*sin((5/12)*2*pi*t)-1, cos((10/16)*2*pi*(t - 0.6))) 
beta3 <- function(t) 2.35 * (1/(1+exp(-5*(t-0.5))) - 0.5) 
beta4 <- function(t) sin(2*pi*t)

coef_fun <- function(t) {
  
  b1 <- beta1(t)
  b2 <- beta2(t)
  b3 <- beta3(t)
  b4 <- beta4(t)
  b <- cbind(b1, b2, b3, b4)
  
  return(b)
  
}


### simulation
time1 <- Sys.time()

sim_time_mat <- matrix(nrow = length(n_seq), ncol = 2)
for (n in n_seq) {
  
  cat('sample size:',which(n == n_seq),'of',length(n_seq), '\n')
  
  # evaluation points
  t_eval <- seq(0, 1, length.out = sqrt(n))
  
  # data generation
  set.seed(99 * n / 100)
  sample_gen <- sample_gen_fun(n, rho_x, L_seq, rho_err, sigma_err)
  
  ###  all data on (0,1)
  t_list <- sample_gen$t_list
  y_list <- sample_gen$y_list
  x_mat <- sample_gen$x_mat
  
  L_vec <- unlist(lapply(t_list, length))
  t_vec <- unlist(t_list)
  y_vec <- unlist(y_list)
  
  
  # domain of inferential interest
  c1 <- 0; c2 <- 0.5
  domain_con <- c(c1, c2)
  
  # estimation w/ all data
  time01 <- Sys.time()
  
  B_vec_hat_null <- B_vec_hat_fun(t_eval = t_eval, 
                                  t_list = t_list, 
                                  y_list = y_list, 
                                  x_mat = x_mat, 
                                  h_con = h_con, 
                                  comp_con = comp_con, 
                                  domain_con = domain_con)
  
  B_vec_hat_alt <- B_vec_hat_fun(t_eval = t_eval, 
                                 t_list = t_list, 
                                 y_list = y_list, 
                                 x_mat = x_mat, 
                                 h = h)
  time02 <- Sys.time()
  
  time_1 <- as.numeric(time02 - time01, units = 'secs')
  
  
  # domain of inferential interest
  c1 <- 0.4; c2 <- 0.5
  domain_con <- c(c1, c2)
  
  # estimation w/ all data
  time01 <- Sys.time()
  
  B_vec_hat_null <- B_vec_hat_fun(t_eval = t_eval, 
                                  t_list = t_list, 
                                  y_list = y_list, 
                                  x_mat = x_mat, 
                                  h_con = h_con, 
                                  comp_con = comp_con, 
                                  domain_con = domain_con)
  
  B_vec_hat_alt <- B_vec_hat_fun(t_eval = t_eval, 
                                 t_list = t_list, 
                                 y_list = y_list, 
                                 x_mat = x_mat, 
                                 h = h)
  time02 <- Sys.time()
  
  time_2 <- as.numeric(time02 - time01, units = 'secs')
  
  sim_time_mat[which(n == n_seq), ] <- c(time_1, time_2)
  
  print(sim_time_mat[1:which(n == n_seq), ])
}


par(mfrow = c(1,1))
par(mar=c(4.5,4.5,3,1)+0.1)
ind <- 1:length(n_seq)
plot(n_seq[ind], sim_time_mat[ind, 1], type = 'l', 
     lwd = 2, lty = 1, col = 'black', 
     cex.axis = 1.25, cex.lab = 1.5, cex.main = 1.5,
     ylim = range(c(0, sim_time_mat[ind, ])),
     xlab = 'Sample size',
     ylab = 'Computatino time (sec)',
     main = 'Kernel estimation')
lines(n_seq[ind], sim_time_mat[ind, 2],
      lwd = 4, lty = 2, col = 'blue')
grid()
legend('bottomright',
       legend = c('I = [0, 0.5]', 'I = [0.4, 0.5]'),
       lwd = 2, lty = c(1, 2), col = c('black', 'blue'),
       cex = 1.5, bty = 'n')
