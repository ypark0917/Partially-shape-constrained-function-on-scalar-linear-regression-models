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
setwd("C:/Users/sdy897/OneDrive - University of Texas at San Antonio/Research/convex-constraints/Rcode")

source("cons_spline_source_datv1.R")
load("dat1_schizophrenia.RData")
load("schizophrenia.rda")
dat = schizophrenia

#names(dat)
#table(dat$id) # mostly 4 or 3

range(dat$Week) # (0,6)
length(unique(dat$id)) # 437 subjects

id.placebo = unique(dat[dat$TxDrug==0, "id"])
id.treat = unique(dat[dat$TxDrug==1, "id"])

# plot(seq(0,6, length=10), rep(4,10), type="n", ylim=c(0,8), ylab="Item 79", xlab="week")
# for(i in id.placebo){
#   tmp = dat[dat$id == i,]
#   lines(tmp$Week, tmp$imps79)
# }
# for(j in id.treat){
#   tmp = dat[dat$id == j,]
#   lines(tmp$Week, tmp$imps79, col="blue", lty=3)
# }
# 
# mean(table(dat$Week)[c(1,2,4,7)]/437)
# mean(table(dat$Week)[c(3,5,6)]/437)

#####################
n = length(unique(dat$id))
B=500
p=2
a1=1.1
# c1 <- 0; c2 <- 6
#c1<-5; c2<-6
#shape.const.1="decreasing"

ndegree=1;
t_eval = seq(0,6, length=50)
tau = sort(unique(dat$Week))
comp_con <- 1

#######
df = data.frame(subject = dat$id, time = dat$Week, y = dat$imps79, X0=rep(1, nrow(dat)), X1 = dat$TxDrug)
id.subj = sort(unique(df$subject))
########

t_list<- y_list <- list()
x_mat <- c()
for (i in 1:length(id.subj)){
  tmp <- dat[dat$id == id.subj[i] , ]
  t_list[[i]] <- tmp$Week
  y_list[[i]] <- tmp$imps79
  x_mat <- rbind(x_mat, c(1,tmp$TxDrug[1]))
}

L_vec <- unlist(lapply(t_list, length))
t_vec <- unlist(t_list)
y_vec <- unlist(y_list)

#####################################
int.set<-list()

tmp1 <- seq(1.5, 5, length=5)
for(i1 in 1:length(tmp1)){
  tmp<-list(c(tmp1[i1]))
  int.set <- append(int.set, tmp)
}

nlength.2=3;
#nlength.3=2
# nknot=2
tmp1<-seq(1, 3, length=nlength.2+2)[-c(1,(nlength.2+2))]
#tmp1 <- tmp1[-c(1,length(tmp1))]

tmp2<- seq(3, 5, length=nlength.2+2)[-c(1,(nlength.2+2))]
#tmp2<-tmp2[-c(1, length(tmp2))]

for(i1 in 1:length(tmp1)){
  for(i2 in 1:length(tmp2)){
    tmp<-list(c(tmp1[i1], tmp2[i2]))
    int.set <- append(int.set, tmp)
  }
}

# # nknot=3
# tmp1<-seq(1, 2.4, length=(nlength.3+2))[-c(1,(nlength.3+2))]
# 
# tmp2<- seq(2.4, 3,length=(nlength.3+2))[-c(1,(nlength.3+2))]
# 
# tmp3<- seq(3, 5,length=(nlength.3+2))[-c(1,(nlength.3+2))]
# 
# for(i1 in 1:length(tmp1)){
#   for(i2 in 1:length(tmp2)){
#     for(i3 in 1:length(tmp3)){
#       tmp<-list(c(tmp1[i1], tmp2[i2], tmp3[i3]))
#       int.set <- append(int.set, tmp)
#     }
#   }
# }

################################
# 2/27/25 for revision - check convexity [3,6]
################################


pval.res.T1.all <- pval.res.T0 <- pval.res.T1 <- pval.res.T2 <- pval.res.T2.all<-c()
#op.kint = 1.5 # restricted to the one internal knot
op.kint = c(1,4) #c(2.5 , 3.5) # from CV
var.kint = op.kint
#intercept==FALSE

  c1<-3
  c2<-6
  shape.const.1= "convexity"
  
  domain_con <- c(c1, c2)
  ndegree <-1
  ################################
  ## unconstrained fit
  # op.kint = op.knot.dat1(df, int.set.list = int.set, ndegree=1, cv=10)
  nsubj= length(unique(df$subject))
  id.subj = sort(unique(df$subject))
  t1= sort(unique(df$time))  
  tau=sort(unique(df$time)) 
  
  cbasis0 = cSpline(t1, knots = var.kint , degree = ndegree) 
  cbasis0.all = cSpline(t_eval, knots = var.kint , degree = ndegree) 
  
  cbasis = (cbind( rep(1, nrow(cbasis0)), cbasis0))
  cbasis.all = (cbind( rep(1, nrow(cbasis0.all)), cbasis0.all))#[,-6]
  
  # Bmat = as.matrix( bdiag(as.matrix(cbasis), as.matrix(cbasis0)) ) 
  # Bmat.all = as.matrix( bdiag(cbasis.all, cbasis0.all) ) 
  
  Bmat = rbind(
    cbind(as.matrix(cbasis), matrix(0, nrow(cbasis), ncol(cbasis0))),
    cbind(matrix(0, nrow(cbasis0), ncol(cbasis)), as.matrix(cbasis0))
  )
  
  Bmat.all =rbind(
    cbind(as.matrix(cbasis.all), matrix(0, nrow(cbasis.all), ncol(cbasis0.all))),
    cbind(matrix(0, nrow(cbasis0.all), ncol(cbasis.all)), as.matrix(cbasis0.all))
  )
  Q=0;C=0
  for (i1 in id.subj){
    df_sub= df[which(df$subject==i1),]
    ni = nrow(df_sub)
    ti = df_sub$time
    Bi = Bmat  
    Xmat= as.matrix( data.frame( df_sub[, c('X0','X1')]) )
    yi = as.matrix(df_sub$y)
    wi = diag( 1/ni, ni )
    
    basis_id= which( t1 %in% ti )
    Ui = matrix(0, ncol= ncol(Bmat), nrow= nrow(Xmat))
    for (j1 in 1:ni){
      j2= basis_id[j1]
      Ui[j1,] = t(Xmat[j1,]) %*% Bmat[c(j2, nrow(cbasis)+j2), ]  # dep on p
    }
    Qi = t(Ui)%*%wi%*%Ui
    ci =  t(yi)%*%wi%*%Ui
    
    Q= Q+Qi
    C= C+ci
  }
  
  
  theta_u= solve(Q)%*% t(C) #unconstrained est 
  #### unconstrained fit
  uncost_fit = (Bmat %*% theta_u)
  unconst_b0hat=uncost_fit[1:length(tau)]
  unconst_b1hat=uncost_fit[(length(tau)+1):(2*length(tau))]
  unconst_res= data.frame( t=tau, b0hat=unconst_b0hat, b1hat=unconst_b1hat)
  
  uncost_fit.all = Bmat.all %*% theta_u
  unconst_b0hat.all=uncost_fit.all[1:length(t_eval)]
  unconst_b1hat.all=uncost_fit.all[(length(t_eval)+1):(2*length(t_eval))] 
  unconst_res.teval= data.frame( t=t_eval, b0hat=unconst_b0hat.all, b1hat=unconst_b1hat.all)
  
  ##########
  # constrained
  var.kint.con = sort(c(var.kint, 3))
  
  cbasis0.con = cSpline(t1, knots = var.kint.con , degree = ndegree)
  cbasis0.all.con = cSpline(t_eval, knots = var.kint.con , degree = ndegree) 
  
  # cbasis0 = cSpline(t1, knots = var.kint , degree = ndegree) 
  # cbasis0.all = cSpline(t_eval, knots = var.kint , degree = ndegree) 
  
  #####
  
  cbasis.con = (cbind( rep(1, nrow(cbasis0.con)), cbasis0.con))
  cbasis.all.con = (cbind( rep(1, nrow(cbasis0.all.con)), cbasis0.all.con))#[,-6]
  
  # cbasis = (cbind( rep(1, nrow(cbasis0)), cbasis0))
  # cbasis.all = (cbind( rep(1, nrow(cbasis0.all)), cbasis0.all))#[,-6]
  
  Bmat.con = rbind(
    cbind(as.matrix(cbasis), matrix(0, nrow(cbasis.con), ncol(cbasis0.con))),
    cbind(matrix(0, nrow(cbasis0), ncol(cbasis)), as.matrix(cbasis0.con))
  )
  
  Bmat.all.con =rbind(
    cbind(as.matrix(cbasis.all), matrix(0, nrow(cbasis.all.con), ncol(cbasis0.all.con))),
    cbind(matrix(0, nrow(cbasis0.all), ncol(cbasis.all)), as.matrix(cbasis0.all.con))
  )
  
  Q=0;C=0
  for (i1 in id.subj){
    df_sub= df[which(df$subject==i1),]
    ni = nrow(df_sub)
    ti = df_sub$time
    Bi = Bmat.con  
    Xmat= as.matrix( data.frame( df_sub[, c('X0','X1')]) )
    yi = as.matrix(df_sub$y)
    wi = diag( 1/ni, ni )
    
    basis_id= which( t1 %in% ti )
    Ui = matrix(0, ncol= ncol(Bmat.con), nrow= nrow(Xmat))
    for (j1 in 1:ni){
      j2= basis_id[j1]
      Ui[j1,] = t(Xmat[j1,]) %*% Bmat.con[c(j2, nrow(cbasis)+j2), ]  # dep on p
    }
    Qi = t(Ui)%*%wi%*%Ui
    ci =  t(yi)%*%wi%*%Ui
    
    Q= Q+Qi
    C= C+ci
  }
  
  
  ###########
  # # intercept
  # knot_rm0 = c(1:ncol(cbasis))
  # 
  # # p=1
  # #if (is.numeric(range1)){ 
  #   inc_id1 = which( tau>= c1 & tau<= c2)
  #   #knot_rm1= which( cbasis0[min(inc_id),] ==1   )
  #   #knot_rm2= which( cbasis0[max(inc_id), 1:ncol(cbasis0)] ==0 )
  #   
  #   knot_rm11= which( round(cbasis0[min(inc_id1),],3) ==1   )
  #   knot_rm12= which( round(cbasis0[max(inc_id1), 1:ncol(cbasis0)],3) ==0 )
  #   knot_rm1 = c(knot_rm11, knot_rm12)
  # #} else{
  # #  knot_rm1 =c(1:ncol(cbasis0))
  # #}
  
  knot_c0 = c() #c(1:ncol(cbasis))[!(c(1:ncol(cbasis)) %in% knot_rm0)]
  #knot_c1 = c(3,4,5) # c(1:ncol(cbasis0))[!(c(1:ncol(cbasis0)) %in% knot_rm1)]
  knot_c1 = c(4,5)
  
  amat = matrix(0, nrow=length(knot_c1),ncol=ncol(Bmat.con))
# convexity
  for(i1 in 1:length(knot_c1)){
    amat[i1, ncol(cbasis) + knot_c1[i1]]=1
  }
  

  ###########
  b= rep(10^-10, nrow(amat)) 
  qp.res= qprog(Q, C, amat, b )
  
  #est= Bmat.con%*%qp.res$thetahat #const est
  
  ##### constrained fit

  const_b0hat.all= (Bmat.all.con%*%qp.res$thetahat)[1:length(t_eval)] 
  const_b1hat.all= (Bmat.all.con%*%qp.res$thetahat)[(length(t_eval)+1):(2*length(t_eval))] 
  const_res.teval = data.frame( t=t_eval, b0hat=const_b0hat.all , b1hat=const_b1hat.all)
  
  Bs_vec_hat_alt=rbind(unconst_res.teval$b0hat, unconst_res.teval$b1hat)
  Bs_vec_hat_null=rbind(const_res.teval$b0hat, const_res.teval$b1hat)
  

  #Bs_vec_hat_null=rbind(res2$const$b0hat, res2$const$b1hat)
  
  par(mfrow=c(1,2))
  plot(t_eval, Bs_vec_hat_alt[1,], type="l", xlab="week", ylab=expression(beta[0](week)),lwd=2, main="intercept", ylim=c(0,7)) #ylim=range(c(Bs_vec_hat_alt[1,],Bs_vec_hat_null[1,])))
  lines(t_eval, Bs_vec_hat_null[1,], col="blue", lwd=2, lty=1)
  abline(h=0, lty=3)
  #abline(v=domain_con, col="red", lty=2)
  
  plot(t_eval, Bs_vec_hat_alt[2,], type="l", xlab="week", ylab=expression(beta[1](week)), lwd=2, main="treatment", ylim=c(-2,0.5)  ) #ylim=range(c(Bs_vec_hat_alt[2,],Bs_vec_hat_null[2,])))
  lines(t_eval, Bs_vec_hat_null[2,], col="blue", lwd=2, lty=1)
  abline(h=0, lty=3)
  abline(v=domain_con, col="red", lty=2)
  
  T0 <- 0
  sse_null <- sse_alt <- 0
  #t_eval <- tau
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
  
  #  if(intercept==TRUE){
  gamma_hat <- t(x_mat) %*% x_mat / n
  ind_eval <- (t_eval >= domain_con[1] & t_eval <= domain_con[2])
  T1 <- sum(diag(t(Bs_vec_hat_null[1:p,ind_eval] - Bs_vec_hat_alt[1:p,ind_eval]) %*% gamma_hat%*%
                   (Bs_vec_hat_null[1:p,ind_eval] - Bs_vec_hat_alt[1:p,ind_eval]))) * diff(t_eval)[1]
  
  
  T2 <-sum(diag(t(Bs_vec_hat_null[1:p,ind_eval] - Bs_vec_hat_alt[1:p,ind_eval]) %*% 
                  (Bs_vec_hat_null[1:p,ind_eval] - Bs_vec_hat_alt[1:p,ind_eval]))) * diff(t_eval)[1]
  
  T1_all <-  sum(diag(t(Bs_vec_hat_null[1:p,] - Bs_vec_hat_alt[1:p,]) %*%gamma_hat %*% 
                        (Bs_vec_hat_null[1:p,] - Bs_vec_hat_alt[1:p,]))) * diff(t_eval)[1]
  
  T2_all <- sum(diag(t(Bs_vec_hat_null[1:p,] - Bs_vec_hat_alt[1:p,]) %*% 
                       (Bs_vec_hat_null[1:p,] - Bs_vec_hat_alt[1:p,]))) * diff(t_eval)[1]
  
  #####
 
  sample_id= rep( 1:n,  unlist(lapply(t_list, length)) )
  y_fit_null_list <- y_fit_alt_list <- resid_alt_list <- list()
  for (i in 1:n) {
    
    # cat(i, '\n')
    
    tB_vec_fit_null <- sapply(1:p,
                              function (j) spline(t_eval, Bs_vec_hat_null[j,], xout = t_list[[i]])$y
    )
    
    tB_vec_fit_alt <- sapply(1:p,
                             function (j) spline(t_eval, Bs_vec_hat_alt[j,], xout = t_list[[i]])$y
    )
    
    y_fit_null_list[[i]] <- c(tB_vec_fit_null %*% x_mat[i,])
    y_fit_alt_list[[i]] <- c(tB_vec_fit_alt %*% x_mat[i,])
    
    resid_alt_list[[i]] <- y_fit_alt_list[[i]] - y_list[[i]]
    
  }
  
  
  # Mammen's wild bootstrap with moment correction
  boot_pr <- (5 + sqrt(5))/10
  
  # wild bootstrap
  boot_T0 <- boot_T1 <-boot_T2 <- boot_T1_all <- boot_T2_all <- c()
  for (bb in 1:B) {
    
    #if (b%%100==0) cat('Bootstrap sampling:',b,'of',B, '\n')
    
    y_list_b <- lapply(1:n,
                       function (i) {
                         y_fit_null_list[[i]] + 
                           sample(c((1 - sqrt(5))/2, (1 + sqrt(5))/2), 
                                  size = L_vec[i], replace = TRUE, 
                                  prob = c(boot_pr, 1 - boot_pr)) * 
                           resid_alt_list[[i]]
                       })
    
    y_vec.boot <- unlist(y_list_b)
    
    df.boot= data.frame(cbind(subject= sample_id, time=t_vec, y=y_vec.boot, X1=df[ ,-c(1:3)])) #, X1= X[,1], X2=X[,2], X3=X[,3], X4=X[,4]  )
    names(df.boot)[c(4:5)] <- c('X0','X1')
  
    
    #res.boot= vcm_irr_dat1( df=df.boot, ndegree = ndegree, t_eval=t_eval, var.kint=var.kint, nknot = NA,
    #                        const=const, range1=domain_con)
    
    Q=0;C=0
    for (i1 in sample_id){
      df.boot_sub= df.boot[which(df.boot$subject==i1),]
      ni = nrow(df.boot_sub)
      ti = df.boot_sub$time
      Bi = Bmat  
      Xmat= as.matrix( data.frame( df.boot_sub[, c('X0','X1')]) )
      yi = as.matrix(df.boot_sub$y)
      wi = diag( 1/ni, ni )
      
      basis_id= which( t1 %in% ti )
      Ui = matrix(0, ncol= ncol(Bmat), nrow= nrow(Xmat))
      for (j1 in 1:ni){
        j2= basis_id[j1]
        Ui[j1,] = t(Xmat[j1,]) %*% Bmat[c(j2, nrow(cbasis)+j2), ]  # dep on p
      }
      Qi = t(Ui)%*%wi%*%Ui
      ci =  t(yi)%*%wi%*%Ui
      
      Q= Q+Qi
      C= C+ci
    }
    
    
    theta_u= solve(Q)%*% t(C) #unconstrained est 
    #### unconstrained fit
    uncost_fit.all = Bmat.all %*% theta_u
    unconst_b0hat.all=uncost_fit.all[1:length(t_eval)]
    unconst_b1hat.all=uncost_fit.all[(length(t_eval)+1):(2*length(t_eval))] 
    unconst_res.teval= data.frame( t=t_eval, b0hat=unconst_b0hat.all, b1hat=unconst_b1hat.all)
    
    Bs_vec_hat_alt_b=rbind(unconst_res.teval$b0hat, unconst_res.teval$b1hat)
    
    
    ###########
    Q=0;C=0
    for (i1 in sample_id){
      df.boot_sub= df.boot[which(df.boot$subject==i1),]
      ni = nrow(df.boot_sub)
      ti = df.boot_sub$time
      Bi = Bmat.con  
      Xmat= as.matrix( data.frame( df.boot_sub[, c('X0','X1')]) )
      yi = as.matrix(df.boot_sub$y)
      wi = diag( 1/ni, ni )
      
      basis_id= which( t1 %in% ti )
      Ui = matrix(0, ncol= ncol(Bmat.con), nrow= nrow(Xmat))
      for (j1 in 1:ni){
        j2= basis_id[j1]
        Ui[j1,] = t(Xmat[j1,]) %*% Bmat.con[c(j2, nrow(cbasis)+j2), ]  # dep on p
      }
      Qi = t(Ui)%*%wi%*%Ui
      ci =  t(yi)%*%wi%*%Ui
      
      Q= Q+Qi
      C= C+ci
    }
    
    b= rep(10^-10, nrow(amat)) 
    qp.res= qprog(Q, C, amat, b )
    
    ##### constrained fit
 
    const_b0hat.all= (Bmat.all.con%*%qp.res$thetahat)[1:length(t_eval)] 
    const_b1hat.all= (Bmat.all.con%*%qp.res$thetahat)[(length(t_eval)+1):(2*length(t_eval))] 
    const_res.teval = data.frame( t=t_eval, b0hat=const_b0hat.all , b1hat=const_b1hat.all)
    
     Bs_vec_hat_null_b=rbind(const_res.teval$b0hat, const_res.teval$b1hat)
    
    T0_b <- 0; 
    sse_null_b <- sse_alt_b <- 0
    for (i in 1:n) {
      
      # cat(i, '\n')
      
      tB_vec_fit_null_b <- sapply(1:p,
                                  function (j) spline(t_eval, Bs_vec_hat_null_b[j,], xout = t_list[[i]])$y
      )
      
      tB_vec_fit_alt_b <- sapply(1:p,
                                 function (j) spline(t_eval, Bs_vec_hat_alt_b[j,], xout = t_list[[i]])$y
      )
      
      T0_b <- T0_b + mean(((tB_vec_fit_null_b - tB_vec_fit_alt_b) %*% x_mat[i,])^2)
      
      sse_null_b <- sse_null_b + mean((y_list_b[[i]] - tB_vec_fit_null_b %*% x_mat[i,])^2)
      sse_alt_b <- sse_alt_b + mean((y_list_b[[i]] - tB_vec_fit_alt_b %*% x_mat[i,])^2)
      
    }
    
    #    if(intercept==TRUE){
    gamma_hat <- t(x_mat) %*% x_mat / n
    ind_eval <- (t_eval >= domain_con[1] & t_eval <= domain_con[2])
    
    T1_b <- sum(diag(t(Bs_vec_hat_null_b[1:p,ind_eval] - Bs_vec_hat_alt_b[1:p,ind_eval]) %*%gamma_hat%*%
                       (Bs_vec_hat_null_b[1:p,ind_eval] - Bs_vec_hat_alt_b[1:p,ind_eval])))* diff(t_eval)[1]
    
    T2_b <-sum(diag(t(Bs_vec_hat_null_b[1:p,ind_eval] - Bs_vec_hat_alt_b[1:p,ind_eval]) %*%
                      (Bs_vec_hat_null_b[1:p,ind_eval] - Bs_vec_hat_alt_b[1:p,ind_eval])))* diff(t_eval)[1]
    
    T1_all_b <- sum(diag(t(Bs_vec_hat_null_b[1:p,] - Bs_vec_hat_alt_b[1:p,]) %*%gamma_hat%*%
                           (Bs_vec_hat_null_b[1:p,] - Bs_vec_hat_alt_b[1:p,])))* diff(t_eval)[1]
    
    T2_all_b <- sum(diag(t(Bs_vec_hat_null_b[1:p,] - Bs_vec_hat_alt_b[1:p,]) %*%
                           (Bs_vec_hat_null_b[1:p,] - Bs_vec_hat_alt_b[1:p,])))* diff(t_eval)[1]
    
   
    
    boot_T0[bb] <- T0_b
    boot_T1[bb] <- T1_b
    boot_T2[bb] <- T2_b
    boot_T1_all[bb] <- T1_all_b
    boot_T2_all[bb] <- T2_all_b
    
    if(bb %% 50 ==0) print(bb)
  }
  
  
  mean(boot_T1_all > T1_all)
  mean(boot_T0 > T0)
  mean(boot_T1 > T1)
  mean(boot_T2 > T2)
  mean(boot_T2_all > T2_all)
  
  pval.res.T1.all<- c(pval.res.T1.all, )
  pval.res.T0<- c(pval.res.T0, mean(boot_test$boot_T0 > T0))
  pval.res.T1<- c(pval.res.T1, mean(boot_test$boot_T1 > T1))
  pval.res.T2<- c(pval.res.T2, mean(boot_test$boot_T2 > T2))
  pval.res.T2.all<- c(pval.res.T2.all,mean(boot_test$boot_T2_all > T2_all))
  

  plot(t_eval, unconst_b1hat.all, type="l")
  lines(t_eval, const_b1hat.all, col="blue")
  

round(matrix(pval.res.T0, ncol=2),2)
round(matrix(pval.res.T1, ncol=2),2)
round(matrix(pval.res.T1.all, ncol=2),2)
round(matrix(pval.res.T2, ncol=2),2)
round(matrix(pval.res.T2.all, ncol=2),2)


###########################################
# for paper result on p-values for increasing/ decreasing constraint
##########################################
pval.res.T1.all <- pval.res.T0 <- pval.res.T1 <- pval.res.T2 <- pval.res.T2.all<-c()
#op.kint = 1.5 # restricted to the one internal knot
op.kint = c(1, 4) #c(2.5 , 3.5) # from CV
#intercept==FALSE
const.list<-list(c(0,6), c(0,3), c(3,6), c(4,6), c(5,6),
                 c(0,6), c(0,3), c(3,6), c(4,6), c(5,6))
shape.list<-list("decreasing","decreasing","decreasing","decreasing","decreasing",
                 "increasing","increasing","increasing","increasing","increasing")

#t_eval <- tau
for(l in 1:length(const.list)){

  c1<-const.list[[l]][1]; 
  c2<-const.list[[l]][2]
  shape.const.1=shape.list[[l]]
  
  domain_con <- c(c1, c2)
  ## unconstrained fit
  # op.kint = op.knot.dat1(df, int.set.list = int.set, ndegree=1, cv=10)
  res= vcm_irr_dat1( df=df, ndegree=1, t_eval=t_eval, const=shape.const.1, var.kint = op.kint, nknot = NA,
                     range1=domain_con)
  
  Bs_vec_hat_alt=rbind(res$unconst_teval$b0hat, res$unconst_teval$b1hat)
  #Bs_vec_hat_alt=rbind(res$unconst$b0hat, res$unconst$b1hat)
  
  ## constrained fit
  op.kint2 = op.kint
  if(c1==0 & c2 !=6){
    #op.kint[which(abs(op.kint - c1)<0.08)] = c1
    op.kint2[which(abs(op.kint2 - c2)< a1)] = c2
    op.kint.con <- sort(unique(c(op.kint2, c2)))
  }else if(c1!=0 & c2==6){
    op.kint2[which(abs(op.kint2 - c1)< a1)] = c1
    op.kint.con <- sort(unique(c(op.kint2, c1)))
  }else if (c1==0 & c2==6){
    op.kint.con <- op.kint2
  }else{
    op.kint2[which(abs(op.kint2 - c1)< a1)] = c1
    op.kint2[which(abs(op.kint2 - c2)< a1)] = c2
    op.kint.con <- sort(unique(c(op.kint2, domain_con)))
  }
  
  loc.id <- which(op.kint.con > c1 & op.kint.con < c2)
  if(length(loc.id)>0){
    op.kint.con <- op.kint.con[-loc.id]
  }

  res2= vcm_irr_con_dat1( df=df, ndegree=ndegree, t_eval=t_eval, const=shape.const.1,var.kint = op.kint ,var.kint.con = op.kint.con, nknot = NA,
                      range1=domain_con)
  Bs_vec_hat_null=rbind(res2$const_teval$b0hat, res2$const_teval$b1hat)
  #Bs_vec_hat_null=rbind(res2$const$b0hat, res2$const$b1hat)
  
  par(mfrow=c(1,2))
  plot(t_eval, Bs_vec_hat_alt[1,], type="l", xlab="week", ylab=expression(beta[0](week)),lwd=2, main="intercept", ylim=c(0,7)) #ylim=range(c(Bs_vec_hat_alt[1,],Bs_vec_hat_null[1,])))
  lines(t_eval, Bs_vec_hat_null[1,], col="blue", lwd=2, lty=1)
  abline(h=0, lty=3)
  #abline(v=domain_con, col="red", lty=2)

  plot(t_eval, Bs_vec_hat_alt[2,], type="l", xlab="week", ylab=expression(beta[1](week)), lwd=2, main="treatment", ylim=c(-2,0.5)  ) #ylim=range(c(Bs_vec_hat_alt[2,],Bs_vec_hat_null[2,])))
  lines(t_eval, Bs_vec_hat_null[2,], col="blue", lwd=2, lty=1)
  abline(h=0, lty=3)
  abline(v=domain_con, col="red", lty=2)
  
  T0 <- 0
  sse_null <- sse_alt <- 0
  #t_eval <- tau
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
  
#  if(intercept==TRUE){
    gamma_hat <- t(x_mat) %*% x_mat / n
    ind_eval <- (t_eval >= domain_con[1] & t_eval <= domain_con[2])
    T1 <- sum(diag(t(Bs_vec_hat_null[1:p,ind_eval] - Bs_vec_hat_alt[1:p,ind_eval]) %*% gamma_hat%*%
                     (Bs_vec_hat_null[1:p,ind_eval] - Bs_vec_hat_alt[1:p,ind_eval]))) * diff(t_eval)[1]
    
    
    T2 <-sum(diag(t(Bs_vec_hat_null[1:p,ind_eval] - Bs_vec_hat_alt[1:p,ind_eval]) %*% 
                    (Bs_vec_hat_null[1:p,ind_eval] - Bs_vec_hat_alt[1:p,ind_eval]))) * diff(t_eval)[1]
    
    T1_all <-  sum(diag(t(Bs_vec_hat_null[1:p,] - Bs_vec_hat_alt[1:p,]) %*%gamma_hat %*% 
                          (Bs_vec_hat_null[1:p,] - Bs_vec_hat_alt[1:p,]))) * diff(t_eval)[1]
    
    T2_all <- sum(diag(t(Bs_vec_hat_null[1:p,] - Bs_vec_hat_alt[1:p,]) %*% 
                         (Bs_vec_hat_null[1:p,] - Bs_vec_hat_alt[1:p,]))) * diff(t_eval)[1]
    
  # }else if(intercept==FALSE){
  #   gamma_hat <- t(x_mat[ ,-1]) %*% x_mat[ ,-1] / n
  #   ind_eval <- (t_eval >= domain_con[1] & t_eval <= domain_con[2])
  #   T1 <- sum((t(Bs_vec_hat_null[2:p,ind_eval] - Bs_vec_hat_alt[2:p,ind_eval]) %*%
  #                    (Bs_vec_hat_null[2:p,ind_eval] - Bs_vec_hat_alt[2:p,ind_eval]))*gamma_hat) * diff(t_eval)[1]
  #   
  #   
  #   T2 <-sum((t(Bs_vec_hat_null[2:p,ind_eval] - Bs_vec_hat_alt[2:p,ind_eval]) %*% 
  #                   (Bs_vec_hat_null[2:p,ind_eval] - Bs_vec_hat_alt[2:p,ind_eval]))) * diff(t_eval)[1]
  #   
  #   T1_all <-  sum((t(Bs_vec_hat_null[2:p,] - Bs_vec_hat_alt[2:p,]) %*%
  #                         (Bs_vec_hat_null[2:p,] - Bs_vec_hat_alt[2:p,]))*gamma_hat) * diff(t_eval)[1]
  #   
  #   T2_all <- sum((t(Bs_vec_hat_null[2:p,] - Bs_vec_hat_alt[2:p,]) %*% 
  #                        (Bs_vec_hat_null[2:p,] - Bs_vec_hat_alt[2:p,]))) * diff(t_eval)[1]
  #   
  # }
  

  #####
  boot_test <- boot_fun_sp_dat1(B = B, 
                                df=df, 
                                ndegree=ndegree, 
                                t_eval=t_eval, 
                                t_list=t_list, 
                                y_list=y_list, 
                                x_mat=x_mat,
                                nknot=NA,
                                var.kint.con=op.kint.con,
                                var.kint = op.kint,
                                const=shape.const.1, 
                                range1=domain_con
                                )
  
  pval.res.T1.all<- c(pval.res.T1.all, mean(boot_test$boot_T1_all > T1_all))
  pval.res.T0<- c(pval.res.T0, mean(boot_test$boot_T0 > T0))
  pval.res.T1<- c(pval.res.T1, mean(boot_test$boot_T1 > T1))
  pval.res.T2<- c(pval.res.T2, mean(boot_test$boot_T2 > T2))
  pval.res.T2.all<- c(pval.res.T2.all,mean(boot_test$boot_T2_all > T2_all))
  
  print(c(l, pval.res.T0[l]))
}

round(matrix(pval.res.T0, ncol=2),2)
round(matrix(pval.res.T1, ncol=2),2)
round(matrix(pval.res.T1.all, ncol=2),2)
round(matrix(pval.res.T2, ncol=2),2)
round(matrix(pval.res.T2.all, ncol=2),2)

##############################################
##############################################
# code for Figure in the draft (05/13/24)
p=2
ndegree=1
var.kint = c(1,4) # optimal from CV selection
var.kint_const = c(1,3,4) # constraint: [0,3] decreasing and [3,6] constant

# overall decreasing fit using previous df
res0= vcm_irr_con_dat1_V2( df=df, ndegree=ndegree, t_eval=t_eval, const="decreasing",var.kint = var.kint ,var.kint.con = var.kint_const, nknot = NA,
                           range1=c(0,6))
Bs_vec_hat_null0 = rbind(res0$const_teval$b0hat, res0$const_teval$b1hat)

# with constratint defined by our finding
range2=c(3,6)
load("dat1_schizophrenia.RData")
colnames(x_mat)= c( 'intercept', 'treatment' )
X_mat= cbind( subject= 1:nrow(x_mat), x_mat)

L_vec <- unlist(lapply(t_list, length))
t_vec <- unlist(t_list)
y_vec <- unlist(y_list)

n=length(L_vec)
sample_id= rep( 1:n, L_vec )

df0= data.frame(subject= sample_id, time=t_vec, y=y_vec )
df= merge(df0, X_mat)
head(df)

t_unique = sort(unique( df$time))
t_eval =  seq( min(t_unique), max(t_unique), length=50)

knots= sort( c(var.kint, min(t_unique), max(t_unique) ) )
knots_const= sort( c(var.kint_const, min(t_unique), max(t_unique) ) )

nsubj= length(unique(df$subject))
id.subj = sort(unique(df$subject))
t1= sort(unique(df$time))  
tau=sort(unique(df$time)) 

ibasis0 = iSpline(t1, knots = var.kint , degree = ndegree) 
ibasis0.all = iSpline(t_eval, knots = var.kint , degree = ndegree) 

ibasis = (cbind( rep(1,nrow(ibasis0)), ibasis0))
ibasis.all = (cbind( rep(1, nrow(ibasis0.all)), ibasis0.all)) 

###
ibasis0_const =  iSpline(t1, knots = var.kint_const , degree = ndegree) 
ibasis0_const.all = iSpline(t_eval, knots = var.kint_const , degree = ndegree) 

ibasis_const = (cbind( rep(1, nrow(ibasis0_const)), ibasis0_const))
ibasis_const.all = (cbind( rep(1, nrow(ibasis0_const.all)), ibasis0_const.all))#[,-6]

###
Bmat_const = as.matrix( bdiag(ibasis, ibasis_const ) ) 
Bmat_const.all = as.matrix( bdiag(ibasis.all, ibasis_const.all ) ) # to evaluate from t_eval grids

Bmat = as.matrix( bdiag(ibasis, ibasis) ) 
Bmat.all = as.matrix( bdiag(ibasis.all, ibasis.all ) ) # to evaluate from t_eval grids

#########
##unconst
Q=0;C=0
for (i1 in 1:nsubj){
  df_sub= df[which(df$subject==i1),]
  ni = nrow(df_sub)
  ti = df_sub$time
  Bi = Bmat  
  Xmat= as.matrix( data.frame( df_sub[,  colnames(x_mat)]) )
  yi = as.matrix(df_sub$y)
  wi = diag( 1/ni, ni )
  
  basis_id= which( t1 %in% ti )
  Ui = matrix(0, ncol= ncol(Bmat), nrow= nrow(Xmat))
  for (j1 in 1:ni){
    j2= basis_id[j1]
    Ui[j1,] = t(Xmat[j1,]) %*% Bmat[c(j2, nrow(ibasis)+j2),]  # dep on p
  }
  Qi = t(Ui)%*%wi%*%Ui
  ci =  t(yi)%*%wi%*%Ui
  
  Q= Q+Qi
  C= C+ci
}

tol=10^-10
theta_u= ginv(Q)%*% t(C) 

###########
### const

Q=0;C=0
for (i1 in 1:nsubj){
  df_sub= df[which(df$subject==i1),]
  ni = nrow(df_sub)
  ti = df_sub$time
  Bi = Bmat_const  
  Xmat= as.matrix( data.frame( df_sub[,  colnames(x_mat)]) )
  yi = as.matrix(df_sub$y)
  wi = diag( 1/ni, ni )
  
  basis_id= which( t1 %in% ti )
  Ui = matrix(0, ncol= ncol(Bmat_const), nrow= nrow(Xmat))
  for (j1 in 1:ni){
    j2= basis_id[j1]
    Ui[j1,] = t(Xmat[j1,]) %*% Bmat_const[c(j2, nrow(ibasis)+j2 ),]  # dep on p
  }
  Qi = t(Ui)%*%wi%*%Ui
  ci =  t(yi)%*%wi%*%Ui
  
  Q= Q+Qi
  C= C+ci
}

if ( is.numeric(range2) ){
  inc_id2 = which( tau >= range2[1] & tau <= range2[2])
  if(length(inc_id2)>1){
    knot_rm210= which(  ibasis_const[min(inc_id2),] ==1   )
    knot_rm211= which(  ibasis_const[max(inc_id2), 1:ncol(ibasis_const)] ==0 )
    knot_rm21 = c(knot_rm210, knot_rm211)
  } #else{
  #knot_rm21 = c(1:ncol(ibasis_const))
  #}
  knot_rm2 = unique( c(1, knot_rm21) )#1 for the intercept of sp basis 
} else{
  #knot_rm12 =vector()
  knot_rm2 = c(1:ncol(ibasis_const)) 
}   
knot_c2 = c(1:ncol(ibasis_const))[!(c(1:ncol(ibasis_const)) %in% knot_rm2)]

amat0= diag( 0, ncol(ibasis)  )
amat02_2= diag( 1, ncol(ibasis_const)  ) 
amat02_1=  diag( -1, ncol(ibasis_const)  ) 

a21= amat02_1[knot_rm2,]
a21= a21[-1,]
a22= amat02_2[-knot_rm2,]
amat02= rbind( a21, a22, -a22)

amat = as.matrix( bdiag( amat0, amat02  ) ) 

if ( length(which(apply(abs(amat),1,sum)==0)) > 0 ){
  amat = amat[-which(apply(abs(amat),1,sum)==0),]
} 

b= rep(0, nrow(amat))
qp.res= qprog(Q  , C, amat, b ) 

const_b1hat.all= (Bmat_const.all%*%qp.res$thetahat)[1:length(t_eval)]
const_b2hat.all= (Bmat_const.all%*%qp.res$thetahat)[(length(t_eval)+1):(2*length(t_eval))] 
Bs_vec_hat_null = rbind(const_b1hat.all, const_b2hat.all)


uncost_fit.all = Bmat.all %*% theta_u
unconst_b1hat.all=uncost_fit.all[1:length(t_eval)]
unconst_b2hat.all=uncost_fit.all[(length(t_eval)+1):(2*length(t_eval))] 
Bs_vec_hat_alt = rbind(unconst_b1hat.all, unconst_b2hat.all)


###############
load("051324_schizophrenia_result.RData")
###############
# 06/19/24 calculate group mean

y.mat <-matrix(NA, nrow=n, ncol=7)
for(i in 1:n){
  tmp.id <- t_list[[i]]+1
  y.mat[i, tmp.id] <- y_list[[i]]
}

y.mean.0 <- colMeans(y.mat[x_mat[,2]==0, ], na.rm=TRUE)
y.mean.1 <- colMeans(y.mat[x_mat[,2]==1, ], na.rm=TRUE)

png(paste0('figures/Rplot_intro_dat_V3.png'), width = 950, height = 410)

par(mfrow=c(1,2))

par(mar=c(4.5,6,2,1)+0.1)
n.gr=30
set.seed(65) # 324
placebo.gr <- sample(id.placebo, size=n.gr)
treat.gr <- sample(id.treat, size=n.gr)

t_list_1 <-t_list[id.subj %in% placebo.gr]
t_list_2 <-t_list[id.subj %in% treat.gr]

y_list_1 <- y_list[id.subj %in% placebo.gr]
y_list_2 <- y_list[id.subj %in% treat.gr]

plot(t_list_1[[1]], y_list_1[[1]],
     type = 'b', cex = 1.5, pch = 20, col = 'darkslategray3',
     lwd = 1, cex.lab = 1.5, cex.axis = 1.5,
     xlim = c(0,6),
     ylim = c(0.3, 7.5)  , #1.1*range(unlist(y_list[1:50])),  
     xlab = 'Week', ylab = 'Disease severity',
     yaxt="n")
axis(2,at=c(1,3,5,7),labels=c(1,3,5,7))

for (i in 2:length(t_list_1)) {
  lines(t_list_1[[i]], y_list_1[[i]],
        type = 'b', cex = 1, pch = 20, lwd = 1, lty=2, col = 'darkslategray3')
}
for (i in 2:length(t_list_2)) {
  lines(t_list_2[[i]], y_list_2[[i]],
        type = 'b', cex = 1, pch = 20, lwd = 1, lty=2, col = 'salmon2')
}
lines(c(0,1,3,6), y.mean.0[c(1,2,4,7)], lwd = 3, col = 'darkslategray3')
lines(c(0,1,3,6), y.mean.1[c(1,2,4,7)], lwd = 3, col = 'salmon2')

grid(lwd = 1.5, col="grey70")
legend("bottomleft", col=c("salmon2","darkslategray3"),legend=c("treatment","placebo"), lty=1, bg="transparent",cex=1.5, bty = "n")


plot(t_eval, Bs_vec_hat_alt[2,], cex=1.5, cex.axis = 1.5, cex.lab=1.5, type="l", xlab="Week", ylab=expression(beta[1](week)), lwd=2, main="Estimated treatment effect on disease severity", ylim=c(-1.7,0.1) , col="black", lty=2,
     cex.main=1.5) #ylim=range(c(Bs_vec_hat_alt[2,],Bs_vec_hat_null[2,])))
lines(t_eval, Bs_vec_hat_null0[2,], col="black", lwd=2, lty=3)
#lines(t_eval, Bs_vec_hat_null[2,], col="blue", lwd=2, lty=1)
lines(t_eval, Bs_vec_hat_null[2,], col="blue", lwd=2, lty=1)
grid(lwd = 1.5, col="grey70")
dev.off()
#save(t_eval, Bs_vec_hat_alt, Bs_vec_hat_null0, Bs_vec_hat_null, file="051324_schizophrenia_result.RData")

#################################3
##################################
# plot for paper with decreasing constraint on [0,6]
op.kint = c(1,4) # optimal from CV selection
var.kint = op.kint
op.kint.con = c(1,4)

res0= vcm_irr_con_dat1_V2( df=df, ndegree=ndegree, t_eval=t_eval, const="decreasing",var.kint = op.kint ,var.kint.con = op.kint.con, nknot = NA,
                        range1=c(0,6))
Bs_vec_hat_null0 = rbind(res0$const_teval$b0hat, res0$const_teval$b1hat)

lines(t_eval,Bs_vec_hat_null0[2,], col="red")

# plot for paper with constraints on [0,3]

ndegree=1
range1 = c(0,3)

nsubj= length(unique(df$subject))
id.subj = sort(unique(df$subject))
t1= sort(unique(df$time))  
tau=sort(unique(df$time)) 

ibasis0.con = iSpline(t1, knots = var.kint.con , degree = ndegree)
ibasis0.all.con = iSpline(t_eval, knots = var.kint.con , degree = ndegree) 

ibasis0 = iSpline(t1, knots = var.kint , degree = ndegree) 
ibasis0.all = iSpline(t_eval, knots = var.kint , degree = ndegree) 

#####
ibasis.con = (cbind( rep(1, nrow(ibasis0.con)), ibasis0.con))
ibasis.all.con = (cbind( rep(1, nrow(ibasis0.all.con)), ibasis0.all.con))

ibasis = (cbind( rep(1, nrow(ibasis0)), ibasis0))
ibasis.all = (cbind( rep(1, nrow(ibasis0.all)), ibasis0.all))  

Bmat.u = as.matrix( bdiag(ibasis, ibasis) ) # including intercept # put 0 if want to see zero at t=0
Bmat.all.u = as.matrix( bdiag(ibasis.all, ibasis.all) ) 

Q.u=0;C.u =0
for (i1 in id.subj){
  df_sub= df[which(df$subject==i1),]
  ni = nrow(df_sub)
  ti = df_sub$time
  Bi = Bmat.u  
  Xmat= as.matrix( data.frame( df_sub[, c('X0','X1')]) )
  yi = as.matrix(df_sub$y)
  wi = diag( 1/ni, ni )
  
  basis_id= which( t1 %in% ti )
  Ui = matrix(0, ncol= ncol(Bmat.u), nrow= nrow(Xmat))
  for (j1 in 1:ni){
    j2= basis_id[j1]
    Ui[j1,] = t(Xmat[j1,]) %*% Bmat.u[c(j2, nrow(ibasis)+j2), ]  # dep on p
  }
  Qi.u = t(Ui)%*%wi%*%Ui
  ci.u =  t(yi)%*%wi%*%Ui
  
  Q.u= Q.u+Qi.u
  C.u= C.u+ci.u
}
theta_u= solve(Q.u)%*% t(C.u)

uncost_fit.all = Bmat.all.u %*% theta_u
unconst_b0hat.all=uncost_fit.all[1:length(t_eval)]
unconst_b1hat.all=uncost_fit.all[(length(t_eval)+1):(2*length(t_eval))] 

Bs_vec_hat_alt=rbind(unconst_b0hat.all, unconst_b1hat.all)

inc_id1 = which( tau>= range1[1] & tau<= range1[2])
knot_rm11= which( round(ibasis.con[min(inc_id1),],3) ==1   )
knot_rm12= which( round(ibasis.con[max(inc_id1), 1:ncol(ibasis.con)],3) ==0 )
knot_rm1 = c(knot_rm11, knot_rm12)

knot_c1 = c(1:ncol(ibasis.con))[!(c(1:ncol(ibasis.con)) %in% knot_rm1)]

Bmat = as.matrix( bdiag(ibasis, ibasis.con) ) # including intercept # put 0 if want to see zero at t=0
Bmat.all = as.matrix( bdiag(ibasis.all, ibasis.all.con) ) 


amat = matrix(0, nrow=length(knot_c1), ncol=ncol(Bmat))
amat[1, (ncol(ibasis)+2) ] = -1
amat[2, (ncol(ibasis)+3) ] = -1
amat[3, (ncol(ibasis)+4) ] = -1

Q=0;C=0
for (i1 in id.subj){
  df_sub= df[which(df$subject==i1),]
  ni = nrow(df_sub)
  ti = df_sub$time
  Bi = Bmat  
  Xmat= as.matrix( data.frame( df_sub[, c('X0','X1')]) )
  yi = as.matrix(df_sub$y)
  wi = diag( 1/ni, ni )
  
  basis_id= which( t1 %in% ti )
  Ui = matrix(0, ncol= ncol(Bmat), nrow= nrow(Xmat))
  for (j1 in 1:ni){
    j2= basis_id[j1]
    Ui[j1,] = t(Xmat[j1,]) %*% Bmat[c(j2, nrow(ibasis)+j2), ]  # dep on p
  }
  Qi = t(Ui)%*%wi%*%Ui
  ci =  t(yi)%*%wi%*%Ui
  
  Q= Q+Qi
  C= C+ci
}

b= rep(10^-10, nrow(amat)) 
qp.res= qprog(Q, C, amat, b )
est= Bmat%*%qp.res$thetahat #const est

const_b0hat.all= (Bmat.all%*%qp.res$thetahat)[1:length(t_eval)] 
const_b1hat.all= (Bmat.all%*%qp.res$thetahat)[(length(t_eval)+1):(2*length(t_eval))] 

Bs_vec_hat_null = rbind(const_b0hat.all, const_b1hat.all)

b1.hat <- Bs_vec_hat_null[2,]
b1.hat[c(26:length(b1.hat))]<- b1.hat[26]
###############
par(mfrow=c(1,2))
par(mar=c(4.5,6,1,1)+0.1)
n.gr=30
set.seed(65) # 324
placebo.gr <- sample(id.placebo, size=n.gr)
treat.gr <- sample(id.treat, size=n.gr)

t_list_1 <-t_list[id.subj %in% placebo.gr]
t_list_2 <-t_list[id.subj %in% treat.gr]

y_list_1 <- y_list[id.subj %in% placebo.gr]
y_list_2 <- y_list[id.subj %in% treat.gr]

plot(t_list_1[[1]], y_list_1[[1]],
     type = 'b', cex = 1, pch = 20, col = 'darkslategray3',
     lwd = 1, cex.lab = 1, cex.axis = 1,
     xlim = c(0,6),
     ylim = c(0.5, 7.5)  , #1.1*range(unlist(y_list[1:50])),  
     xlab = 'Week', ylab = 'Disease severity')
for (i in 2:length(t_list_1)) {
  lines(t_list_1[[i]], y_list_1[[i]],
        type = 'b', cex = 1, pch = 20, lwd = 1, col = 'darkslategray3')
}
for (i in 2:length(t_list_2)) {
  lines(t_list_2[[i]], y_list_2[[i]],
        type = 'b', cex = 1, pch = 20, lwd = 1, col = 'gray30')
}
grid(lwd = 1.5)

plot(t_eval, Bs_vec_hat_alt[2,], type="l", xlab="week", ylab=expression(beta[1](week)), lwd=2, main="Estimated treatment effect on disease severity", ylim=c(-1.8,0.1) , col="blue", lty=2 ) #ylim=range(c(Bs_vec_hat_alt[2,],Bs_vec_hat_null[2,])))
lines(t_eval, Bs_vec_hat_null0[2,], col="green", lwd=2, lty=1)
#lines(t_eval, Bs_vec_hat_null[2,], col="blue", lwd=2, lty=1)
lines(t_eval, b1.hat, col="blue", lwd=2, lty=1)
grid(lwd = 1.5)
#abline(h=0, lty=3)
#abline(v=3, col="red", lty=2)





################
par(mfrow=c(1,2))
par(mar=c(4.5,6,1,1)+0.1)
plot(t_eval, Bs_vec_hat_alt[1,], type="l", xlab="week", ylab=expression(beta[0](week)),lwd=2, main="Intercept", ylim=c(3,6.5), col="blue", lty=2) #ylim=range(c(Bs_vec_hat_alt[1,],Bs_vec_hat_null[1,])))
lines(t_eval, Bs_vec_hat_null[1,], col="blue", lwd=2, lty=1)
abline(h=0, lty=3)
#abline(v=domain_con, col="red", lty=2)

plot(t_eval, Bs_vec_hat_alt[2,], type="l", xlab="week", ylab=expression(beta[1](week)), lwd=2, main="Treatment", ylim=c(-2,0.5) , col="blue", lty=2 ) #ylim=range(c(Bs_vec_hat_alt[2,],Bs_vec_hat_null[2,])))
lines(t_eval, Bs_vec_hat_null[2,], col="blue", lwd=2, lty=1)
abline(h=0, lty=3)
abline(v=3, col="red", lty=2)

###
library(tvem)
model3 <- tvem(data=df,
               formula=y~X1,               
               #formula=y~X2, 
               id=subject,
               time=time,
               num_knots= 1)
print(model3)
plot(model3)

names(model3)
model3$grid_fitted_coefficients

tvem_time= model3$time_grid
tvem_res2= model3$grid_fitted_coefficients$X2[,1]
tvem_res1= model3$grid_fitted_coefficients$X1[,1]

# ###
# Bs_vec_hat_alt_56 <- Bs_vec_hat_alt; Bs_vec_hat_null_56 <- Bs_vec_hat_null
# Bs_vec_hat_alt_46 <- Bs_vec_hat_alt; Bs_vec_hat_null_46 <- Bs_vec_hat_null
# Bs_vec_hat_alt_06 <- Bs_vec_hat_alt; Bs_vec_hat_null_06 <- Bs_vec_hat_null
# 
# plot(t_eval, Bs_vec_hat_alt_06[1,], type="l", xlab="week", ylab=expression(beta[1](week)), lwd=2)
# lines(t_eval, Bs_vec_hat_null_06[1,], col="blue", lwd=2, lty=1)
# lines(t_eval, Bs_vec_hat_null_46[1,], col="blue", lwd=2, lty=2)
# lines(t_eval, Bs_vec_hat_null_56[1,], col="blue", lwd=2, lty=3)
# abline(v=c(4,5,6), col="red", lty=2)
# 
# plot(t_eval, Bs_vec_hat_alt[1,], type="l", xlab="week", ylab=expression(beta[1](week)), lwd=2)
# lines(t_eval, Bs_vec_hat_null[1,], col="blue", lwd=2, lty=1)
# abline(v=domain_con, col="red", lty=2)

