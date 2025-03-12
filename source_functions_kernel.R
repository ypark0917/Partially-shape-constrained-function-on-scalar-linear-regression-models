

# function-on-scalar functinoal linear regression model
reg_fun <- function(t, x) return(c(coef_fun(t)%*%x))

# within-subject dependence with AR(1) structure
err_subj_fun <- function (t, rho_err, order = 50) {
  
  ar_mat_err <- matrix(NA, nrow = order, ncol = order)
  for (i in 1:order) {
    for (j in i:order) {
      ar_mat_err[i,j] <- ar_mat_err[j,i] <- rho_err^(abs(i-j))
    }
  }
  
  # rand_coef <- ar_mat %*% rnorm(order, 0, sigma)
  rand_coef <- mvrnorm(1, rep(0, order), ar_mat_err)
  
  err <- 0
  for (k in 1:order) {
    err <- err + rand_coef[k] * (1/sqrt(k)) * sin(2*k*pi*t)/sqrt(pi)
  }
  
  return(err)
}

# data generation
sample_gen_fun <- function (n, rho_x, L_seq, rho_err, sigma_err) {
  
  # # covariate
  x_mat <- matrix(nrow = n, ncol = p)
  ar_mat_x <- matrix(NA, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in i:p) {
      ar_mat_x[i,j] <- ar_mat_x[j,i] <- rho_x^(abs(i-j))
    }
  }
  
  L_vec <- c()
  t_list <- y_list <- list()
  
  for (i in 1:n) {
    
    # xi <- pnorm(mvrnorm(1, rep(0, p), ar_mat_x))
    xi <- (runif(p, 0, 1) + rho_x * runif(1, 0, 1)) / (1 + rho_x)

    Li <- sample(L_seq, 1)
    ti <- sort(rbeta(Li, beta_para[1], beta_para[2]))
    # ti <- sort(runif(Li, 0, 1))
    
    yi <- reg_fun(ti, xi) + err_subj_fun(ti, rho_err) + rnorm(length(ti), 0, sigma_err)
    
    L_vec[i] <- Li
    t_list[[i]] <- ti
    y_list[[i]] <- yi
    x_mat[i,] <- xi
  }
  
  return (list(t_list = t_list,
               y_list = y_list,
               L_vec = L_vec,
               x_mat = x_mat))
}

# Epanechnikov kernel
K <- function (u) (3/4) * (1 - u^2) * ifelse(abs(u) <= 1, 1, 0)

# kernel weighted averages
K_avg_fun <- function (t_eval_point, t_list, y_list, x_mat, h) {
  
  L_vec <- unlist(lapply(t_list, length))
  t_vec <- unlist(t_list)
  
  n <- length(t_list)
  N <- sum(L_vec)
  p <- ncol(x_mat)
  
  D_mat <- matrix(0, nrow = (2*p), ncol = (2*p))
  d_vec <- rep(0, 2*p)
  for (i in 1:n) {
    
    kww_sum <- matrix(0, nrow = (2*p), ncol = (2*p))
    kwy_sum <- rep(0, 2*p)
    for (l in 1:L_vec[i]) {
      w_il <- c(x_mat[i,], (t_list[[i]][l] - t_eval_point) * x_mat[i,])
      
      kww_sum <- kww_sum + K((t_list[[i]][l] - t_eval_point)/h)/h * w_il %*% t(w_il)
      kwy_sum <- kwy_sum + K((t_list[[i]][l] - t_eval_point)/h)/h * w_il * y_list[[i]][l]
    }
    
    D_mat <- D_mat + kww_sum / L_vec[i]
    d_vec <- d_vec + kwy_sum / L_vec[i]
  }
  
  return(list(D_mat = D_mat,
              d_vec = d_vec))
}


# shape-constrained (monotone) local linear kernel smoothing
B_vec_hat_fun <- function (t_eval, t_list, y_list, x_mat, h = NULL, h_con = NULL, comp_con = NULL, domain_con = NULL) {
  
  L_vec <- unlist(lapply(t_list, length))
  t_vec <- unlist(t_list)
  y_vec <- unlist(y_list)
  
  n <- length(t_list)
  N <- sum(L_vec)
  p <- ncol(x_mat)
  M <- length(t_eval)
  
  D_mat <- matrix(0, nrow = 2*p*M, ncol = 2*p*M)
  d_vec <- rep(0, 2*p*M)

  tA_mat_incr_list <- list()
  tA_mat_bdd_fit_list <- list()
  for (j in 1:p) {
    
    tA_mat_incr_list[[j]] <- matrix(0, nrow = M, ncol = 2*p*M) 
    tA_mat_bdd_fit_list[[j]] <- matrix(0, nrow = M, ncol = 2*p*M)
    
  }
  
  if (is.null(comp_con)) {
    
    # print('alt')
    for (m in 1:M) {
      
      K_avg_m <- K_avg_fun(t_eval[m], t_list, y_list, x_mat, h)
      
      ind_K_avg <- seq(1+2*p*(m-1), 2*p*m)
      
      D_mat[ind_K_avg, ind_K_avg] <- K_avg_m$D_mat
      d_vec[ind_K_avg] <- K_avg_m$d_vec
      
    }
    
    D_mat_eigen <- eigen(D_mat)
    D_mat_eigen_val <- D_mat_eigen$val
    D_mat_eigen_vec <- D_mat_eigen$vec
    
    if (sum(D_mat_eigen_val[D_mat_eigen_val < h^2]) > 0) {
      
      D_mat_eigen_val[D_mat_eigen_val < h^2] <- h^2
      # D_mat_eigen_val <- D_mat_eigen_val / sum(D_mat_eigen_val)
      
      D_mat <- D_mat_eigen_vec %*% diag(D_mat_eigen_val) %*% t(D_mat_eigen_vec)
      
    }
    
    B_vec_hat_alt <- solve(D_mat) %*% d_vec
    # B_vec_hat_alt <- D_mat_eigen_vec %*% diag(inv_D_mat_eigen_val) %*% t(D_mat_eigen_vec) %*% d_vec
    B_vec_hat_alt <- matrix(B_vec_hat_alt, nrow = 2*p, ncol = M)
    
    return (B_vec_hat_alt)
    
  } else {
    
    # print('null')
    
    for (m in 1:M) {
        
      K_avg_m <- K_avg_fun(t_eval[m], t_list, y_list, x_mat, h_con)
      
      ind_K_avg <- seq(1+2*p*(m-1), 2*p*m)
      
      D_mat[ind_K_avg, ind_K_avg] <- K_avg_m$D_mat
      d_vec[ind_K_avg] <- K_avg_m$d_vec
      
      for (j in comp_con) {
        
        if (t_eval[m] >= domain_con[1] && t_eval[m] <= domain_con[2]) {
          
          if (m < M) {
            
            # increasing
            tA_mat_incr_list[[j]][m, 2*p*(m-1) + j] <- -1
            tA_mat_incr_list[[j]][m, 2*p*m + j] <- 1
            
          } 
          
          if (m == M) {
            
            # increasing
            tA_mat_incr_list[[j]][m, 2*p*(M-2) + j] <- -1
            tA_mat_incr_list[[j]][m, 2*p*(M-1) + j] <- 1
            
          }
          
        } else {
          
          
          m_bdd <- NULL
          
          if (t_eval[m] <= domain_con[1] && t_eval[m] > domain_con[1] - h_con) {
            
            m_bdd <- min(which(t_eval >= domain_con[1]))
            
          }
          
          if (t_eval[m] >= domain_con[2] && t_eval[m] < domain_con[2] + h_con) {
            
            m_bdd <- max(which(t_eval <= domain_con[2]))
            
          }
          
          if (is.null(m_bdd) == FALSE && m_bdd != m) {
            
            # boundary fit
            tA_mat_bdd_fit_list[[j]][m, 2*p*(m-1) + j] <- -1 
            tA_mat_bdd_fit_list[[j]][m, 2*p*(m_bdd-1) + j] <- 1 
            
            tA_mat_bdd_fit_list[[j]][m, 2*p*(m-1) + p + j] <- -(t_eval[m_bdd] - t_eval[m])/2
            tA_mat_bdd_fit_list[[j]][m, 2*p*(m_bdd-1) + p + j] <- -(t_eval[m_bdd] - t_eval[m])/2 
            
          }
        }
      }
    }
    
    
    
    tA_mat_incr <- c()
    tA_mat_bdd_fit <- c()
    for (j in 1:p) {
      tA_mat_incr <- rbind(tA_mat_incr, tA_mat_incr_list[[j]])
      tA_mat_bdd_fit <- rbind(tA_mat_bdd_fit, tA_mat_bdd_fit_list[[j]])
    }
    
    A_mat_incr <- t(tA_mat_incr)
    A_mat_bdd_fit <- t(tA_mat_bdd_fit)
    
    A_mat <- cbind(A_mat_incr,
                   A_mat_bdd_fit,
                   -A_mat_bdd_fit)
    
    b_vec <- c(rep(0, nrow(A_mat_incr)),
               # rep(-1e-2, 2 * nrow(A_mat_bdd_fit)))
               rep(-h_con^2, 2 * nrow(A_mat_bdd_fit)))
    
    ind_qp <- which(apply(abs(A_mat), 2, sum) != 0)
    
    A_mat <- A_mat[,ind_qp]
    b_vec <- b_vec[ind_qp]
    
    # preventing numerical error
    D_mat_eigen <- eigen(D_mat)
    D_mat_eigen_val <- D_mat_eigen$val
    D_mat_eigen_vec <- D_mat_eigen$vec
    
    if (sum(D_mat_eigen_val[D_mat_eigen_val < h_con^2]) > 0) {
      
      D_mat_eigen_val[D_mat_eigen_val < h_con^2] <- h_con^2
      # D_mat_eigen_val <- D_mat_eigen_val / sum(D_mat_eigen_val)
      
      D_mat <- D_mat_eigen_vec %*% diag(D_mat_eigen_val) %*% t(D_mat_eigen_vec)
      
    }
    
    result_qp <- solve.QP(D_mat, d_vec, A_mat, b_vec)
    
    B_vec_hat_null <- result_qp$solution
    B_vec_hat_null <- matrix(B_vec_hat_null, nrow = 2*p, ncol = M)
    
    return(B_vec_hat_null)
    
  }
  
}


# wild bootstrap
boot_fun <- function (B = 1000, t_eval, t_list, y_list, x_mat, h, h_con, comp_con, domain_con) {
  
  L_vec <- unlist(lapply(t_list, length))
  t_vec <- unlist(t_list)
  y_vec <- unlist(y_list)
  
  n <- length(t_list)
  N <- sum(L_vec)
  p <- ncol(x_mat)
  M <- length(t_eval)
  
  # fitted values under the alternative
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
  
  y_fit_null_list <- y_fit_alt_list <- resid_alt_list <- list()
  for (i in 1:n) {
    
    # cat(i, '\n')
    
    tB_vec_fit_null <- sapply(1:p,
                              function (j) spline(t_eval, B_vec_hat_null[j,], xout = t_list[[i]])$y
    )
    
    tB_vec_fit_alt <- sapply(1:p,
                             function (j) spline(t_eval, B_vec_hat_alt[j,], xout = t_list[[i]])$y
    )
    
    y_fit_null_list[[i]] <- c(tB_vec_fit_null %*% x_mat[i,])
    y_fit_alt_list[[i]] <- c(tB_vec_fit_alt %*% x_mat[i,])
    
    resid_alt_list[[i]] <- y_fit_alt_list[[i]] - y_list[[i]]
    
  }
  
  # Mammen's wild bootstrap with moment correction
  boot_pr <- (5 + sqrt(5))/10
  
  # wild bootstrap
  boot_T0 <- boot_T1 <- boot_T2 <- boot_T3 <- boot_T4 <- c()
  for (b in 1:B) {
    
    if ((b%%10) == 1) {
      
      cat('Bootstrap sampling:',b,'of',B, '\n')
      
    }
      
    
    y_list_b <- lapply(1:n,
                       function (i) {
                         y_fit_null_list[[i]] + 
                           sample(c((1 - sqrt(5))/2, (1 + sqrt(5))/2), 
                                  size = L_vec[i], replace = TRUE, 
                                  prob = c(boot_pr, 1 - boot_pr)) * 
                           resid_alt_list[[i]]
                       })
    
    B_vec_hat_null_b <- B_vec_hat_fun(t_eval = t_eval, 
                                    t_list = t_list, 
                                    y_list = y_list_b, 
                                    x_mat = x_mat, 
                                    h_con = h_con, 
                                    comp_con = comp_con, 
                                    domain_con = domain_con)
    
    B_vec_hat_alt_b <- B_vec_hat_fun(t_eval = t_eval, 
                                   t_list = t_list, 
                                   y_list = y_list_b, 
                                   x_mat = x_mat, 
                                   h = h)
    
    T0_b <- 0
    sse_null_b <- sse_alt_b <- 0
    for (i in 1:n) {
      
      # cat(i, '\n')
      
      tB_vec_fit_null_b <- sapply(1:p,
                                function (j) spline(t_eval, B_vec_hat_null_b[j,], xout = t_list[[i]])$y
      )
      
      tB_vec_fit_alt_b <- sapply(1:p,
                               function (j) spline(t_eval, B_vec_hat_alt_b[j,], xout = t_list[[i]])$y
      )
      
      T0_b <- T0_b + mean(((tB_vec_fit_null_b - tB_vec_fit_alt_b) %*% x_mat[i,])^2)
      
      # sse_null_b <- sse_null_b + mean((y_list_b[[i]] - tB_vec_fit_null_b %*% x_mat[i,])^2)
      # sse_alt_b <- sse_alt_b + mean((y_list_b[[i]] - tB_vec_fit_alt_b %*% x_mat[i,])^2)
    }
    # F_b <- ifelse(sse_null_b - sse_alt_b > 0, sse_null_b - sse_alt_b, 0) / sse_alt_b
    
    gamma_hat <- t(x_mat) %*% x_mat / n
    ind_eval <- (t_eval >= domain_con[1] & t_eval <= domain_con[2])
    T1_b <- sum(diag(t(B_vec_hat_null_b[1:p,ind_eval] - B_vec_hat_alt_b[1:p,ind_eval]) %*%
                     gamma_hat %*%
                     (B_vec_hat_null_b[1:p,ind_eval] - B_vec_hat_alt_b[1:p,ind_eval]))) * diff(t_eval)[1]
    
    T2_b <- sum(diag(t(B_vec_hat_null_b[1:p,ind_eval] - B_vec_hat_alt_b[1:p,ind_eval]) %*%
                       (B_vec_hat_null_b[1:p,ind_eval] - B_vec_hat_alt_b[1:p,ind_eval]))) * diff(t_eval)[1]
    
    T3_b <- sum(diag(t(B_vec_hat_null_b[1:p,] - B_vec_hat_alt_b[1:p,]) %*%
                       gamma_hat %*%
                       (B_vec_hat_null_b[1:p,] - B_vec_hat_alt_b[1:p,]))) * diff(t_eval)[1]
    
    T4_b <- sum(diag(t(B_vec_hat_null_b[1:p,] - B_vec_hat_alt_b[1:p,]) %*%
                       (B_vec_hat_null_b[1:p,] - B_vec_hat_alt_b[1:p,]))) * diff(t_eval)[1]

    boot_T0[b] <- T0_b
    boot_T1[b] <- T1_b
    boot_T2[b] <- T2_b
    boot_T3[b] <- T3_b
    boot_T4[b] <- T4_b
  }
  
  return (list(boot_T0 = boot_T0,
               boot_T1 = boot_T1,
               boot_T2 = boot_T2,
               boot_T3 = boot_T3,
               boot_T4 = boot_T4))
}



# CV for bandwidth selection
bw_sel_fun <- function (h_seq = seq((c2-c1)/10, (c2-c1)/2, length.out = 10),
                        cv_len = 5,
                        t_eval, t_list, y_list, x_mat, comp_con, domain_con) {
  
  n <- length(t_list)
  
  cv_null <- cv_alt <- matrix(Inf, nrow = length(h_seq), ncol = cv_len)
  
  for (s in 1:length(h_seq)) {
    
    cat(paste0(cv_len, '-fold CV Bandwith selection: ', s, ' of ', length(h_seq), ' (bw = ',round(h_seq[s], 3), ') \n'))
    
    for (cv in 1:cv_len) {
      
      set.seed(cv * cv_len + s * length(h_seq))
      # print(cv)
      
      ind_test <- sort(sample(1:n, n/cv_len, replace = FALSE))
      ind_train <- (1:n)[-ind_test]
      
      B_vec_hat_null <- try (
        B_vec_hat_fun(t_eval = t_eval, 
                      t_list = t_list[ind_train], 
                      y_list = y_list[ind_train], 
                      x_mat = x_mat[ind_train, ], 
                      h_con = h_seq[s], 
                      comp_con = comp_con, 
                      domain_con = domain_con),
        silent = TRUE
      )
      
      B_vec_hat_alt <- try (
        B_vec_hat_fun(t_eval = t_eval,
                      t_list = t_list[ind_train],
                      y_list = y_list[ind_train],
                      x_mat = x_mat[ind_train, ],
                      h = h_seq[s]),
        silent = TRUE
      )
      
      tmp_null <- class(B_vec_hat_null) %in% c('matrix', 'array')
      tmp_alt <- class(B_vec_hat_alt) %in% c('matrix', 'array')
      
      # print(c(tmp_null, tmp_alt))
      
      if (!('FALSE' %in% c(tmp_null, tmp_alt))) {
      # if (!('FALSE' %in% tmp_null)) {
        
        pe_null <- pe_alt <- 0
        for (i in ind_test) {
          
          tB_vec_fit_null <- sapply(1:p,
                                    function (j) spline(t_eval, B_vec_hat_null[j,], xout = t_list[[i]])$y
          )
          
          tB_vec_fit_alt <- sapply(1:p,
                                   function (j) spline(t_eval, B_vec_hat_alt[j,], xout = t_list[[i]])$y
          )
          
          y_pred_null <- c(tB_vec_fit_null %*% x_mat[i,])
          y_pred_alt <- c(tB_vec_fit_alt %*% x_mat[i,])
          
          pe_null <- pe_null + mean((y_pred_null - y_list[[i]])^2)
          pe_alt <- pe_alt + mean((y_pred_alt - y_list[[i]])^2)
          
        }
        
        cv_null[s, cv] <- pe_null
        cv_alt[s, cv] <- pe_alt
        
      }
      
    }
    
  }
  
  # print(cv_null)
  # print(cv_alt)
  
  cv_null <- apply(cv_null, 1, mean)
  cv_alt <- apply(cv_alt, 1, mean)
  
  h_null <- h_seq[which.min(cv_null)[1]]
  h_alt <- h_seq[which.min(cv_alt)[1]]
  # h_alt <- h_null
  
  return (list(h_null = h_null,
               h_alt = h_alt,
               cv_null = cv_null,
               cv_alt = cv_alt)
  )
  
}

