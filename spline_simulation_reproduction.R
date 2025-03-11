library(quadprog)
library(numDeriv)
library(MASS)
library(rstudioapi)
## unconstrained spline fitting
library(tvem) 
## cone projection
library(splines2)
library(Matrix)
library(coneproj)
library(caret)

#setwd('/Users/sdy897/Library/CloudStorage/OneDrive-UniversityofTexasatSanAntonio/Research/convex-constraints/Rcode')
#setwd('C:/Users/sdy897/OneDrive - University of Texas at San Antonio/Research/convex-constraints/Rcode')
setwd('/work/sdy897/shape')

source("source/functions_1101_2023_KH.R")
source("source/cons_spline_source_V3.R")

### simulation settings
# Bootstrap size
B <- 200
# sample size
n <- 500
#kint.res <- list()
load(paste('source/opt.knot.',n,'.MC200.final.RData',sep=''))
nknot.re=NA

a1 <- 0.09 # 0.09 for type I

# domain of inferential interest
c1 <- 0; c2 <- 0.5

# contrast 
shape.const.1="increasing"
shape.const.2="increasing"

# MC size
MC <- 200
#mc <- 1

domain_con <- c(c1, c2)

# component of inferential interest
comp_con <- c(1, 3)

# sparse irregular measurements
L_seq <- 5:15 
#L_seq <- 15:30
# measurement time distribution
beta_para <- c(1, 1.25)

# noise level
rho_err <- 0.1
sigma_err <- 0.01

# regression coefficient functions
p <- 4
rho_x <- 0.5

tau1 <- 0.5
tau2 <- 0.7

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

# evaluation points
t_eval <- seq(0, 1, length.out = 101)
t_eval_re <- seq(c1, c2, length.out = length(t_eval))


### simulation
#time1 <- Sys.time()
#######################################
# MC replication
######################################
# order
ndegree=1

# knot candidates
if(length(kint.res)==0){
  # for the whole data
  
  int.set<-list()
  if(n==500){
    nlength.3=2
    nlength.4=2
    # nknot=3
    tmp1<-seq(0.15, 0.4, length=(nlength.3+2))[-c(1,(nlength.3+2))]
    
    tmp2<- seq(0.4, 0.6,length=(nlength.3+2))[-c(1,(nlength.3+2))]
    
    tmp3<- seq(0.6, 0.85,length=(nlength.3+2))[-c(1,(nlength.3+2))]
    
    for(i1 in 1:length(tmp1)){
      for(i2 in 1:length(tmp2)){
        for(i3 in 1:length(tmp3)){
          tmp<-list(c(tmp1[i1], tmp2[i2], tmp3[i3]))
          int.set <- append(int.set, tmp)
        }
      }
    }    
    # nknot=4
    tmp1<-seq(0.15, 0.35, length=(nlength.4+2))[-c(1,(nlength.4+2))]
    tmp2<- seq(0.35, 0.5, length=(nlength.4+2))[-c(1,(nlength.4+2))]
    tmp3<- seq(0.5, 0.65, length=(nlength.4+2))[-c(1,(nlength.4+2))]
    tmp4<- seq(0.65, 0.85, length=(nlength.4+2))[-c(1,(nlength.4+2))]
    for(i1 in 1:length(tmp1)){
      for(i2 in 1:length(tmp2)){
        for(i3 in 1:length(tmp3)){
          for(i4 in 1:length(tmp4)){
            tmp<-list(c(tmp1[i1], tmp2[i2], tmp3[i3], tmp4[i4]))
            int.set <- append(int.set, tmp)
          }
        }
      }
    }
  }
  
  if(n==100){
    nlength.2=3;
    nlength.3=2
    # nknot=2
    tmp1<-seq(0.2, 0.5, length=nlength.2+2)[-c(1,(nlength.2+2))]
    #tmp1 <- tmp1[-c(1,length(tmp1))]
    
    tmp2<- seq(0.5, 0.8, length=nlength.2+2)[-c(1,(nlength.2+2))]
    #tmp2<-tmp2[-c(1, length(tmp2))]
    
    for(i1 in 1:length(tmp1)){
      for(i2 in 1:length(tmp2)){
        tmp<-list(c(tmp1[i1], tmp2[i2]))
        int.set <- append(int.set, tmp)
      }
      
      int.set <- int.set[-3]
    }
    
    # nknot=3
    tmp1<-seq(0.15, 0.4, length=(nlength.3+2))[-c(1,(nlength.3+2))]
    
    tmp2<- seq(0.4, 0.6,length=(nlength.3+2))[-c(1,(nlength.3+2))]
    
    tmp3<- seq(0.6, 0.85,length=(nlength.3+2))[-c(1,(nlength.3+2))]
    
    for(i1 in 1:length(tmp1)){
      for(i2 in 1:length(tmp2)){
        for(i3 in 1:length(tmp3)){
          tmp<-list(c(tmp1[i1], tmp2[i2], tmp3[i3]))
          int.set <- append(int.set, tmp)
        }
      }
    }
  }
  
  
}

# for partial data
nlength.2=2;
nlength.1=2

int.set.re<-list()
tmp1 <- seq(c1+0.02, c2-0.02, length=nlength.1)
for(i1 in 1:length(tmp1)){
  tmp <- list(tmp1[i1])
  int.set.re <- append( int.set.re, tmp)
}

mdt = (c1+c2)/2
tmp1<-seq(c1+0.02,  mdt-0.02, length=nlength.2)
tmp2<- seq(mdt+0.02, c2-0.02, length=nlength.2)

for(i1 in 1:length(tmp1)){
  for(i2 in 1:length(tmp2)){
    tmp<-list(c(tmp1[i1], tmp2[i2]))
    int.set.re <- append( int.set.re, tmp)
  }
}

###

pval_T0_MC <- pval_T1_MC <- pval_T2_MC <- pval_T1_all_MC <- pval_T2_all_MC <- c()
pval_T0_re_MC <- pval_T1_re_MC <- pval_T2_re_MC <- c()
B_vec_hat_null_MC <- B_vec_hat_alt_MC <- list()
B_vec_hat_null_re_MC <- B_vec_hat_alt_re_MC <- list()

#nknot.res<-c(); nknot.con.res<-c()
#kint.res <- list();

for (mc in 1:MC) {
  
  # data generation
  set.seed(9999 * n + 999 * mc + 9 * p + 10)
  sample_gen <- sample_gen_fun(n, rho_x, L_seq, rho_err, sigma_err)
  
  ###  all data on (0,1)
  t_list <- sample_gen$t_list
  y_list <- sample_gen$y_list
  x_mat <- sample_gen$x_mat
  
  L_vec <- unlist(lapply(t_list, length))
  t_vec <- unlist(t_list)
  y_vec <- unlist(y_list)
  
  # estimation w/ all data
  sample_id= rep( 1:n, sample_gen$L_vec )
  x_mat1 = matrix(ncol=p)#= matrix(0, nrow=length(y_vec), ncol=2)
  for (i1 in 1:nrow(x_mat)){
    aa= x_mat[i1,]
    aa1= matrix( rep(aa, length(t_list[[i1]]) ), ncol=p, byrow = TRUE)
    x_mat1 =  rbind(x_mat1, aa1)
  }
  X = x_mat1[-1,]
  
  df= data.frame(subject= sample_id, time=t_vec, y=y_vec, X1= X[,1], X2=X[,2], X3=X[,3], X4=X[,4] )
  
  ####################
  ### unconstrained
  ####################
  
  #op.kint = op.knot.V2(df, int.set.list = int.set, ndegree=1, cv=10)
  #kint.res[[mc]] = op.kint
  op.kint = kint.res[[mc]]
  
  res= vcm_irr( df=df,ndegree=ndegree, t_eval=t_eval, const=c(shape.const.1, NA,shape.const.2,NA),var.kint = op.kint, nknot = NA,
                range1=domain_con, range2=NA, range3=domain_con, range4=NA)
  Bs_vec_hat_alt=rbind(res$unconst_teval$b1hat,res$unconst_teval$b2hat,res$unconst_teval$b3hat,res$unconst_teval$b4hat)
  
  
  ########################
  #### constrained
  ########################
  op.kint2 = op.kint
  if(c1==0){
    #op.kint[which(abs(op.kint - c1)<0.08)] = c1
    op.kint2[which(abs(op.kint2 - c2)< a1)] = c2
    op.kint.con <- sort(unique(c(op.kint2, c2)))
  }else if(c2==1){
    op.kint2[which(abs(op.kint2 - c1)< a1)] = c1
    op.kint.con <- sort(unique(c(op.kint2, c1)))
  }else{
    op.kint2[which(abs(op.kint2 - c1)< a1)] = c1
    op.kint2[which(abs(op.kint2 - c2)< a1)] = c2
    op.kint.con <- sort(unique(c(op.kint2, domain_con)))
  }
  
  # loc.id <- which(op.kint.con > c1 & op.kint.con < c2)
  # if(length(loc.id)>0){
  #   op.kint.con <- op.kint.con[-loc.id]
  # }
  
  # to make diff
  res2= vcm_irr_con( df=df,ndegree=ndegree, t_eval=t_eval, const=c(shape.const.1, NA,shape.const.2,NA),var.kint = op.kint,
                     var.kint.con = op.kint.con, nknot = NA,range1=domain_con, range2=NA, range3=domain_con, range4=NA)
  
  Bs_vec_hat_null=rbind(res2$const_teval$b1hat,res2$const_teval$b2hat,res2$const_teval$b3hat,res2$const_teval$b4hat)
  
  
  B_vec_hat_null_MC[[mc]] <- Bs_vec_hat_null
  B_vec_hat_alt_MC[[mc]] <- Bs_vec_hat_alt
  
  T0 <- 0
  sse_null <- sse_alt <- 0
  for (i in 1:n) {
    
    # cat(i, '\n')
    
    tB_vec_sfit_null <- sapply(1:p,
                               function (j) spline(t_eval, Bs_vec_hat_null[j,], xout = t_list[[i]])$y
    )
    
    tB_vec_sfit_alt <- sapply(1:p,
                              function (j) spline(t_eval, Bs_vec_hat_alt[j,], xout = t_list[[i]])$y
    )
    
    T0 <- T0 + mean(((tB_vec_sfit_null - tB_vec_sfit_alt) %*% x_mat[i,])^2)
    
    sse_null <- sse_null + mean((y_list[[i]] - c(tB_vec_sfit_null %*% x_mat[i,]))^2)
    sse_alt <- sse_alt + mean((y_list[[i]] - c(tB_vec_sfit_alt %*% x_mat[i,]))^2)
  }
  
  #F0 <- ifelse(sse_null - sse_alt > 0, sse_null - sse_alt, 0) / sse_alt
  
  gamma_hat <- t(x_mat) %*% x_mat / n
  ind_eval <- (t_eval >= domain_con[1] & t_eval <= domain_con[2])
  T1 <- sum(diag(t(Bs_vec_hat_null[1:p,ind_eval] - Bs_vec_hat_alt[1:p,ind_eval]) %*%
                   gamma_hat %*%
                   (Bs_vec_hat_null[1:p,ind_eval] - Bs_vec_hat_alt[1:p,ind_eval]))) * diff(t_eval)[1]
  
  
  T2 <- sum(diag(t(Bs_vec_hat_null[1:p,ind_eval] - Bs_vec_hat_alt[1:p,ind_eval]) %*%
                   (Bs_vec_hat_null[1:p,ind_eval] - Bs_vec_hat_alt[1:p,ind_eval]))) * diff(t_eval)[1]
  
  T1_all <- sum(diag(t(Bs_vec_hat_null[1:p,] - Bs_vec_hat_alt[1:p,]) %*%
                       gamma_hat %*%
                       (Bs_vec_hat_null[1:p,] - Bs_vec_hat_alt[1:p,]))) * diff(t_eval)[1]
  
  
  T2_all <- sum(diag(t(Bs_vec_hat_null[1:p,] - Bs_vec_hat_alt[1:p,]) %*%
                       (Bs_vec_hat_null[1:p,] - Bs_vec_hat_alt[1:p,]))) * diff(t_eval)[1]
  
  boot_test <- boot_fun_sp(B = B, 
                           df=df, 
                           ndegree=ndegree, 
                           t_eval=t_eval, 
                           t_list=t_list, 
                           y_list=y_list, 
                           x_mat=x_mat,
                           nknot=NA,
                           var.kint.con = op.kint.con,
                           var.kint = op.kint,
                           const=c(shape.const.1,NA,shape.const.2,NA), 
                           range1=domain_con, range2=NA, range3=domain_con, range4=NA, partial=FALSE)
  
  
  boot_T0 <- boot_test$boot_T0
  boot_T1 <- boot_test$boot_T1
  boot_T2 <- boot_test$boot_T2
  boot_T1_all <- boot_test$boot_T1_all
  boot_T2_all <- boot_test$boot_T2_all  
  
  
  pval_T0_MC[mc] <- mean(boot_T0 > T0)
  pval_T1_MC[mc] <- mean(boot_T1 > T1)
  pval_T2_MC[mc] <- mean(boot_T2 > T2)
  pval_T1_all_MC[mc] <- mean(boot_T1_all > T1_all)
  pval_T2_all_MC[mc] <- mean(boot_T2_all > T2_all) 
  
  ###############################
  ### partial data on (c1, c2)
  t_list_re <- lapply(1:n, 
                      function (i) (t_list[[i]])[t_list[[i]] >= c1 & t_list[[i]] <= c2]
  )
  y_list_re <- lapply(1:n, 
                      function (i) (y_list[[i]])[t_list[[i]] >= c1 & t_list[[i]] <= c2]
  )
  
  ex_ind <- lapply(1:n,
                   function (i) length(t_list_re[[i]])==0
  )
  
  t_list_re <- t_list_re[ex_ind==FALSE]
  y_list_re <- y_list_re[ex_ind==FALSE]
  x_mat_re <- x_mat[ex_ind==FALSE, ]
  
  L_vec_re <- unlist(lapply(t_list_re, length))
  t_vec_re <- unlist(t_list_re)
  y_vec_re <- unlist(y_list_re)
  n_re <- length(t_list_re)
  
  #n_re_MC[mc] <- n_re
  
  ####
  # estimation w/ partial data
  sample_id_re= rep( 1:n_re,lengths(t_list_re) )
  #head(x_mat)
  x_mat1_re = matrix(0, ncol=p)#= matrix(0, nrow=length(y_vec), ncol=2)
  for (i1 in 1:nrow(x_mat_re)){
    aa= x_mat_re[i1,]
    aa1= matrix( rep(aa, length(t_list_re[[i1]]) ), ncol=p, byrow = TRUE)
    x_mat1_re =  rbind(x_mat1_re, aa1)
  }
  
  X_re = x_mat1_re[-1,]
  
  df_re= data.frame(subject= sample_id_re, time=t_vec_re, y=y_vec_re, X1= X_re[,1], X2=X_re[,2],X3= X_re[,3], X4=X_re[,4]  )
  
  ########################
  ### unconstrained
  #######################
  op.kint.re = op.knot.V2(df=df_re, int.set.list = int.set.re, ndegree=1, cv=10)
  res_re= vcm_irr( df=df_re, ndegree=ndegree, t_eval=t_eval_re, const=c(shape.const.1,NA,shape.const.2,NA), var.kint = op.kint.re, nknot = nknot.re,
                   range1=domain_con, range2=NA, range3=domain_con, range4=NA)
  Bs_vec_hat_alt_re = rbind(res_re$unconst_teval$b1hat,res_re$unconst_teval$b2hat,res_re$unconst_teval$b3hat,res_re$unconst_teval$b4hat)
  
  ########################
  ### constrained
  #######################
  op.kint.con.re = op.knot.con.V2(df=df_re, int.set.list = int.set.re, const=c(shape.const.1,NA,shape.const.2,NA), ndegree=1, cv=10)
  res2_re= vcm_irr_con( df=df_re, ndegree=ndegree, t_eval=t_eval_re, const=c(shape.const.1, NA,shape.const.2,NA),var.kint = op.kint.con.re, var.kint.con= op.kint.con.re, nknot = nknot.re,
                        range1=domain_con, range2=NA, range3=domain_con, range4=NA)
  
  Bs_vec_hat_null_re = rbind(res2_re$const_teval$b1hat, res2_re$const_teval$b2hat, res2_re$const_teval$b3hat, res2_re$const_teval$b4hat)
  
  B_vec_hat_null_re_MC[[mc]] <- Bs_vec_hat_null_re
  B_vec_hat_alt_re_MC[[mc]] <- Bs_vec_hat_alt_re
  
  
  T0_re <- 0
  sse_null <- sse_alt <- 0
  for (i in 1:n_re) {
    
    # cat(i, '\n')
    
    tB_vec_fit_null <- sapply(1:p,
                              function (j) spline(t_eval_re, Bs_vec_hat_null_re[j,], xout = t_list_re[[i]])$y
    )
    
    tB_vec_fit_alt <- sapply(1:p,
                             function (j) spline(t_eval_re, Bs_vec_hat_alt_re[j,], xout = t_list_re[[i]])$y
    )
    
    T0_re <- T0_re + mean(((tB_vec_fit_null - tB_vec_fit_alt) %*% x_mat_re[i,])^2)
    
    sse_null <- sse_null + mean((y_list_re[[i]] - c(tB_vec_fit_null %*% x_mat_re[i,]))^2)
    sse_alt <- sse_alt + mean((y_list_re[[i]] - c(tB_vec_fit_alt %*% x_mat_re[i,]))^2)
  }
  
  
  gamma_hat <- t(x_mat_re) %*% x_mat_re / n_re
  ind_eval <- (t_eval_re >= domain_con[1] & t_eval_re <= domain_con[2])
  T1_re <- sum(diag(t(Bs_vec_hat_null_re[1:p,ind_eval] - Bs_vec_hat_alt_re[1:p,ind_eval]) %*%
                      gamma_hat %*%
                      (Bs_vec_hat_null_re[1:p,ind_eval] - Bs_vec_hat_alt_re[1:p,ind_eval]))) * diff(t_eval_re)[1]
  
  T2_re <- sum(diag(t(Bs_vec_hat_null_re[1:p,ind_eval] - Bs_vec_hat_alt_re[1:p,ind_eval]) %*%
                      (Bs_vec_hat_null_re[1:p,ind_eval] - Bs_vec_hat_alt_re[1:p,ind_eval]))) * diff(t_eval_re)[1]
  
  boot_test_re <- boot_fun_sp(B = B, 
                              df=df_re, 
                              ndegree=ndegree, 
                              t_eval=t_eval_re, 
                              t_list=t_list_re, 
                              y_list=y_list_re, 
                              x_mat=x_mat_re,
                              nknot=1,
                              var.kint.con = op.kint.con.re,
                              var.kint = op.kint.re,
                              const=c(shape.const.1,NA,shape.const.2,NA), 
                              range1=domain_con, range2=NA, range3=domain_con, range4=NA, partial=TRUE)
  
  
  boot_T0_re <- boot_test_re$boot_T0
  boot_T1_re <- boot_test_re$boot_T1
  boot_T2_re <- boot_test_re$boot_T2
  
  #boot_F <- boot_test$boot_F
  
  
  pval_T0_re_MC[mc] <- mean(boot_T0_re > T0_re)
  pval_T1_re_MC[mc] <- mean(boot_T1_re > T1_re)
  pval_T2_re_MC[mc] <- mean(boot_T2_re > T2_re) 
  
  
  
  cat('MC replication:',mc,'of',MC, ':',( mean(boot_T0 > T0)<0.05), (mean(boot_T1_all > T1_all)<0.05), (mean(boot_T1_re > T1_re)<0.05),'\n')
}

save(pval_T0_MC, pval_T0_re_MC, pval_T1_MC, pval_T1_re_MC, pval_T2_MC, pval_T2_re_MC, pval_T1_all_MC,  pval_T2_all_MC,
     B_vec_hat_null_MC, B_vec_hat_null_re_MC,
     B_vec_hat_alt_MC, B_vec_hat_alt_re_MC, file= paste('res/032824_final_V3_MC_', MC,'_B_',B,'_n_',n,'_c1_',c1*10,'_c2_',c2*10,'_a1_',a1*100,'.RData', sep = ''))
