#https://hbiostat.org/rmsc/long
# https://pubmed.ncbi.nlm.nih.gov/10534248/
# http://ndl.ethernet.edu.et/bitstream/123456789/88820/1/Davis2002%20-%20Statistical%20Methods%20for%20the%20Analysis%20of%20Repeated%20Measurements.pdf


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


#load("/Users/sdy897/Downloads/mixor/data/schizophrenia.rda")
#setwd('/Users/sdy897/Library/CloudStorage/OneDrive-UniversityofTexasatSanAntonio/Research/convex-constraints/Rcode')
setwd("C:/Users/sdy897/OneDrive - University of Texas at San Antonio/Research/convex-constraints/Rcode")
#library(medicaldata)
#data("cdystonia")
#dat <- as.data.frame(cdystonia)
#save(cdystonia, file="cdystonia.RData")

source("cons_spline_source_datv1.R")
load("cdystonia.RData")

dat = cdystonia
# names(dat)
table(dat$week)
# length(unique(dat$id))

# generate unique id
uid = (100*as.vector(dat$site) + as.vector(dat$id))
length(unique(uid)) # 109 subjects

dat = cbind(dat, uid)
dat$treat = as.vector(dat$treat)
uid.trt1 = unique(dat[dat$treat==1, "uid"]) # 10000 units
uid.trt2 = unique(dat[dat$treat==2, "uid"]) # 5000 units
uid.trt3 = unique(dat[dat$treat==3, "uid"]) # placebo

# plot(seq(0,16, length=10), rep(0,10), type="n", ylim=range(dat$twstrs), ylab="TWSTRS", xlab="week")
# for(i1 in uid.trt1){
#   tmp = dat[dat$uid == i1,]
#   lines(tmp$week, tmp$twstrs)
# }
# for(i2 in uid.trt2){
#   tmp = dat[dat$uid == i2,]
#   lines(tmp$week, tmp$twstrs, col="blue")
# }
# for(i3 in uid.trt3){
#   tmp = dat[dat$uid == i3,]
#   lines(tmp$week, tmp$twstrs, col="red")
# }

##############################
# treatment vs. placebo

treatment <- ifelse(dat$treat %in% c(1,2),1,0)
dat <- cbind(dat,treatment)

######################

dat$sex <- (as.vector(dat$sex)-1 )
dat$age <- as.vector(dat$age)
#####################
n = length(unique(dat$uid))
B=500
p=4
a1=1
# c1 <- 0; c2 <- 6
#c1<-0; c2<-16
#shape.const.1="decreasing"

ndegree=1;
t_eval = seq(0,16, length=50) # originally 100
#domain_con <- c(c1, c2)
comp_con <- 1

#######
df = data.frame(subject = dat$uid, time = dat$week, y = dat$twstrs, X0= rep(1, nrow(dat)), X1 = dat$treatment, X2 = dat$sex, X3 = dat$age)
id.subj = sort(unique(df$subject))
########

t_list<- y_list <- list()
x_mat <- c()
for (i in 1:length(id.subj)){
  tmp <- dat[dat$uid == id.subj[i] , ]
  t_list[[i]] <- tmp$week
  y_list[[i]] <- tmp$twstrs
  x_mat <- rbind(x_mat, c(1,tmp$treatment[1], tmp$sex[1], tmp$age[1]))
}

L_vec <- unlist(lapply(t_list, length))
t_vec <- unlist(t_list)
y_vec <- unlist(y_list)

#########
int.set<-list()
nlength.1=4
nlength.2=3;

# nknot=1
#tmp1 <- seq(5,11, length=nlength.1+2)[-c(1,(nlength.1+2))]
tmp1 <- c(2,4,8)
for(i1 in 1:length(tmp1)){
    tmp<-list(c(tmp1[i1]))
    int.set <- append(int.set, tmp)
}
# nknot=2
# tmp1<-seq(2, 8, length=nlength.2+2)[-c(1,(nlength.2+2))]
# tmp2<- seq(8, 14, length=nlength.2+2)[-c(1,(nlength.2+2))]

tmp1 <- c(2,4,6)
tmp2 <- c(8,12)
for(i1 in 1:length(tmp1)){
  for(i2 in 1:length(tmp2)){
    tmp<-list(c(tmp1[i1], tmp2[i2]))
    int.set <- append(int.set, tmp)
  }
}


#op.kint <- op.knot.dat2(df=df, int.set.list = int.set, ndegree=1, cv=10)

pval.res<-c()
op.kint = 8 # from the optimal choice

const.list<-list(c(0,4), c(4,16),c(0,5), c(5,16),
                 c(0,4), c(4,16),c(0,5), c(5,16))
shape.list<-list("decreasing","decreasing","decreasing","decreasing",
                 "increasing","increasing","increasing","increasing")

for(l in 1:length(const.list)){
  
  c1<-const.list[[l]][1]; 
  c2<-const.list[[l]][2]
  shape.const.1=shape.list[[l]]
  
  domain_con <- c(c1, c2)

  res= vcm_irr_dat2( df=df, ndegree=1, t_eval=t_eval, const=c(shape.const.1,NA,NA), var.kint = op.kint, nknot = NA,
                     range1=domain_con, range2=NA, range3=NA)
  Bs_vec_hat_alt=rbind(res$unconst_teval$b0hat, res$unconst_teval$b1hat, res$unconst_teval$b2hat, res$unconst_teval$b3hat)
  
  ## constrained fit
  op.kint2 = op.kint
  if(c1==0 & c2 !=16){
    #op.kint[which(abs(op.kint - c1)<0.08)] = c1
    op.kint2[which(abs(op.kint2 - c2)< a1)] = c2
    op.kint.con <- sort(unique(c(op.kint2, c2)))
  }else if(c1!=0 & c2==16){
    op.kint2[which(abs(op.kint2 - c1)< a1)] = c1
    op.kint.con <- sort(unique(c(op.kint2, c1)))
  }else if (c1==0 & c2==16){
    op.kint.con <- op.kint2
  }else{
    op.kint2[which(abs(op.kint2 - c1)< a1)] = c1
    op.kint2[which(abs(op.kint2 - c2)< a1)] = c2
    op.kint.con <- sort(unique(c(op.kint2, domain_con)))
  }
  
  res2= vcm_irr_con_dat2( df=df, ndegree=ndegree, t_eval=t_eval, const=c(shape.const.1,NA,NA), var.kint = op.kint, var.kint.con=op.kint.con, 
                          nknot = NA, range1=domain_con, range2=NA, range3=NA)
  Bs_vec_hat_null=rbind(res2$const_teval$b0hat, res2$const_teval$b1hat,res2$const_teval$b2hat,res2$const_teval$b3hat)
  
  
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
  
  ###
  par(mfrow=c(2,2))
  plot(t_eval, Bs_vec_hat_alt[1,], type="l", xlab="week", ylab=expression(beta[0](week)),lwd=2, main="intercept", ylim=c(0,60)) #ylim=range(c(Bs_vec_hat_alt[1,],Bs_vec_hat_null[1,])))
  lines(t_eval, Bs_vec_hat_null[1,], col="blue", lwd=2, lty=1)
  abline(h=0, lty=3)
  #abline(v=domain_con, col="red", lty=2)
  
  plot(t_eval, Bs_vec_hat_alt[2,], type="l", xlab="week", ylab=expression(beta[1](week)), lwd=2, main="treatment", ylim=c(-10,10)  ) #ylim=range(c(Bs_vec_hat_alt[2,],Bs_vec_hat_null[2,])))
  lines(t_eval, Bs_vec_hat_null[2,], col="blue", lwd=2, lty=1)
  abline(h=0, lty=3)
  abline(v=domain_con, col="red", lty=2)
  
  plot(t_eval, Bs_vec_hat_alt[3,], type="l", xlab="week", ylab=expression(beta[2](week)), lwd=2, main="sex",  ylim=c(-6,2)) #ylim=range(c(Bs_vec_hat_alt[3,],Bs_vec_hat_null[3,])))
  lines(t_eval, Bs_vec_hat_null[3,], col="blue", lwd=2, lty=1)
  abline(h=0, lty=3)
  
  plot(t_eval, Bs_vec_hat_alt[4,], type="l", xlab="week", ylab=expression(beta[2](week)), lwd=2, main="age",  ylim=c(-0.3, 0.3)) #ylim=range(c(Bs_vec_hat_alt[4,],Bs_vec_hat_null[4,])))
  lines(t_eval, Bs_vec_hat_null[4,], col="blue", lwd=2, lty=1)
  abline(h=0, lty=3)

  ######
  boot_test <- boot_fun_sp_dat2(B = B, 
                                df=df, 
                                ndegree=ndegree, 
                                t_eval=t_eval, 
                                t_list=t_list, 
                                y_list=y_list, 
                                x_mat=x_mat,
                                nknot=NA,
                                var.kint.con = op.kint.con,
                                var.kint = op.kint,
                                const=c(shape.const.1,NA,NA), 
                                range1=domain_con, range2=NA, range3=NA)
  
  
  boot_T0 <- boot_test$boot_T0
  boot_T1 <- boot_test$boot_T1
  boot_T2 <- boot_test$boot_T2
  boot_T1_all <- boot_test$boot_T1_all
  boot_T2_all <- boot_test$boot_T2_all  
  
  pval.res<- c(pval.res, mean(boot_test$boot_T1_all > T1_all))
  mean(boot_test$boot_T0 > T0)
  mean(boot_test$boot_T1 > T1)
  mean(boot_test$boot_T2 > T2)
  mean(boot_test$boot_T2_all > T2_all) 
  
  print(c(l, pval.res[l] ))
}

#########################################
# 04/22/24
# final fit under given constraints
## decreasing over [0,4] and increasing over [4,16]

op.kint = 8

# res= vcm_irr_dat2( df=df, ndegree=1, t_eval=t_eval, const=c("increasing",NA,NA), var.kint = op.kint, nknot = NA,
#                    range1=c(0,16), range2=NA, range3=NA)
# Bs_vec_hat_alt=rbind(res$unconst_teval$b0hat, res$unconst_teval$b1hat, res$unconst_teval$b2hat, res$unconst_teval$b3hat)

## constrained fit
var.kint = op.kint
var.kint.con = c(4,8)
ndegree=1
range1 = c(0,16)

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

Bmat.u = as.matrix( bdiag(ibasis, ibasis, ibasis, ibasis) ) # including intercept # put 0 if want to see zero at t=0
Bmat.all.u = as.matrix( bdiag(ibasis.all, ibasis.all, ibasis.all,ibasis.all) ) 

#Bmat = as.matrix( bdiag(ibasis, ibasis0.con, ibasis0, ibasis0) ) 
#Bmat.all = as.matrix( bdiag(ibasis.all, ibasis0.all.con, ibasis0.all,ibasis0.all) ) 
#head(df)

Q.u=0;C.u =0
for (i1 in id.subj){
  df_sub= df[which(df$subject==i1),]
  ni = nrow(df_sub)
  ti = df_sub$time
  Bi = Bmat.u  
  Xmat= as.matrix( data.frame( df_sub[, c('X0','X1','X2','X3')]) )
  yi = as.matrix(df_sub$y)
  wi = diag( 1/ni, ni )
  
  basis_id= which( t1 %in% ti )
  Ui = matrix(0, ncol= ncol(Bmat.u), nrow= nrow(Xmat))
  for (j1 in 1:ni){
    j2= basis_id[j1]
    Ui[j1,] = t(Xmat[j1,]) %*% Bmat.u[c(j2, nrow(ibasis)+j2, 2*nrow(ibasis)+j2, 3*nrow(ibasis)+j2), ]  # dep on p
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
unconst_b2hat.all=uncost_fit.all[(2*length(t_eval)+1):(3*length(t_eval))] 
unconst_b3hat.all=uncost_fit.all[(3*length(t_eval)+1):(4*length(t_eval))] 

Bs_vec_hat_alt=rbind(unconst_b0hat.all, unconst_b1hat.all, unconst_b2hat.all, unconst_b3hat.all)

inc_id1 = which( tau>= range1[1] & tau<= range1[2])
knot_rm11= which( round(ibasis.con[min(inc_id1),],3) ==1   )
knot_rm12= which( round(ibasis.con[max(inc_id1), 1:ncol(ibasis.con)],3) ==0 )
knot_rm1 = c(knot_rm11, knot_rm12)
knot_rm2 = c(1:ncol(ibasis))
knot_rm3 =  c(1:ncol(ibasis)) 


knot_c0 = c(1:ncol(ibasis))[!(c(1:ncol(ibasis)) %in% knot_rm0)]
knot_c1 = c(1:ncol(ibasis.con))[!(c(1:ncol(ibasis.con)) %in% knot_rm1)]
knot_c2 = c(1:ncol(ibasis0))[!(c(1:ncol(ibasis0)) %in% knot_rm2)]
knot_c3 = c(1:ncol(ibasis0))[!(c(1:ncol(ibasis0)) %in% knot_rm3)]

Bmat = as.matrix( bdiag(ibasis, ibasis.con, ibasis, ibasis) ) # including intercept # put 0 if want to see zero at t=0
Bmat.all = as.matrix( bdiag(ibasis.all, ibasis.all.con, ibasis.all,ibasis.all) ) 


amat = matrix(0, nrow=length(knot_c1), ncol=ncol(Bmat))
amat[1, (ncol(ibasis)+2)] = -1
amat[2, (ncol(ibasis)+3)] = -1
amat[3, (ncol(ibasis)+4)] = 1
amat[4, (ncol(ibasis)+5)] = 1


Q=0;C=0
for (i1 in id.subj){
  df_sub= df[which(df$subject==i1),]
  ni = nrow(df_sub)
  ti = df_sub$time
  Bi = Bmat  
  Xmat= as.matrix( data.frame( df_sub[, c('X0','X1','X2','X3')]) )
  yi = as.matrix(df_sub$y)
  wi = diag( 1/ni, ni )
  
  basis_id= which( t1 %in% ti )
  Ui = matrix(0, ncol= ncol(Bmat), nrow= nrow(Xmat))
  for (j1 in 1:ni){
    j2= basis_id[j1]
    Ui[j1,] = t(Xmat[j1,]) %*% Bmat[c(j2, nrow(ibasis)+j2, 2*nrow(ibasis)+j2, 3*nrow(ibasis)+j2), ]  # dep on p
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
const_b2hat.all= (Bmat.all%*%qp.res$thetahat)[(2*length(t_eval)+1):(3*length(t_eval))]
const_b3hat.all= (Bmat.all%*%qp.res$thetahat)[(3*length(t_eval)+1):(4*length(t_eval))]

Bs_vec_hat_null = rbind(const_b0hat.all, const_b1hat.all, const_b2hat.all, const_b3hat.all)

png(paste0('figures/Rplot_cdystonia_V2.png'), width = 1300*2/3, height = 400*2/3)
par(mfrow=c(1,3))
par(mar=c(5,5,2,5))
#par(mar=c(4.5,6,2,1)+0.1)
# plot(t_eval, Bs_vec_hat_alt[1,], type="l", xlab="week", ylab=expression(beta[0](week)),lwd=2, main="Intercept", ylim=c(35,50), col="blue", lty=2) #ylim=range(c(Bs_vec_hat_alt[1,],Bs_vec_hat_null[1,])))
# lines(t_eval, Bs_vec_hat_null[1,], col="blue", lwd=2, lty=1)
#abline(h=0, lty=3)
#abline(v=domain_con, col="red", lty=2)
cex.set = 1.75

plot(t_eval, Bs_vec_hat_alt[2,], type="l", xlab="week", ylab=expression(beta[1](week)), lwd=2, main="Treatment", ylim=c(-8,7), col="black", lty=2, cex.main=cex.set, cex.lab=cex.set, cex.sub=cex.set, cex.axis = cex.set ) #ylim=range(c(Bs_vec_hat_alt[2,],Bs_vec_hat_null[2,])))
lines(t_eval, Bs_vec_hat_null[2,], col="black", lwd=2, lty=1)
abline(v=4, lty=3, col="red", lwd=2)
points(c(0,var.kint.con,16), -8 + rep( 0.1,length=4) + 0.3 , pch=4, cex=2, col="blue", lwd=2 )
grid(lwd = 1.5)
abline(h=0, lty=3)
#abline(v=domain_con, col="red", lty=2)

plot(t_eval, Bs_vec_hat_alt[4,], type="l", xlab="week", ylab=expression(beta[2](week)), lwd=2, main="Age",  ylim=c(-0.15, 0.15), col="black", lty=2, cex.main=cex.set, cex.lab=cex.set, cex.axis = cex.set) #ylim=range(c(Bs_vec_hat_alt[4,],Bs_vec_hat_null[4,])))
lines(t_eval, Bs_vec_hat_null[4,], col="black", lwd=2, lty=1)
points(c(0,var.kint,16), -0.15 + rep( 0,length=3) + 0.01 , pch=4, cex=2, col="blue", lwd=2 )
grid(lwd = 1.5)
abline(h=0, lty=3)

plot(t_eval, Bs_vec_hat_alt[3,], type="l", xlab="week", ylab=expression(beta[3](week)), lwd=2, main="Sex",  ylim=c(-5.5,1), col="black", lty=2, cex.main=cex.set, cex.lab=cex.set, cex.axis = cex.set) #ylim=range(c(Bs_vec_hat_alt[3,],Bs_vec_hat_null[3,])))
lines(t_eval, Bs_vec_hat_null[3,], col="black", lwd=2, lty=1)
points(c(0,var.kint,16), -5.5 + rep( 0,length=3) + 0.2 , pch=4, cex=2, col="blue", lwd=2 )
grid(lwd = 1.5)
abline(h=0, lty=3)
dev.off()

mean(df[ df$time==0 & df$X1==0, "y"])
mean(df[ df$time==0 & df$X1==1, "y"])

#res2= vcm_irr_con_dat2( df=df, ndegree=ndegree, t_eval=t_eval, const=c(shape.const.1,NA,NA), var.kint = op.kint, var.kint.con=op.kint.con, 
#                        nknot = NA, range1=domain_con, range2=NA, range3=NA)
#Bs_vec_hat_null=rbind(res2$const_teval$b0hat, res2$const_teval$b1hat,res2$const_teval$b2hat,res2$const_teval$b3hat)



########################################3
# 2/13/25 constrained fit and confidence band

## decreasing over [0,4] and increasing over [4,16]

op.kint = 8

# res= vcm_irr_dat2( df=df, ndegree=1, t_eval=t_eval, const=c("increasing",NA,NA), var.kint = op.kint, nknot = NA,
#                    range1=c(0,16), range2=NA, range3=NA)
# Bs_vec_hat_alt=rbind(res$unconst_teval$b0hat, res$unconst_teval$b1hat, res$unconst_teval$b2hat, res$unconst_teval$b3hat)

## constrained fit
var.kint = op.kint
#var.kint.con = c(4,8)
var.kint.con0 = c(4,8)
var.kint.con = c(3.7, 4.05 , 8)
ndegree=1
range1 = c(0,16)

nsubj= length(unique(df$subject))
id.subj = sort(unique(df$subject))
t1= sort(unique(df$time))  
tau=sort(unique(df$time)) 

ibasis0.con0 = iSpline(t1, knots = var.kint.con0 , degree = ndegree) # 6 times 4
ibasis0.all.con0 = iSpline(t_eval, knots = var.kint.con0 , degree = ndegree) # 50 times 4


ibasis0.con = iSpline(t1, knots = var.kint.con , degree = ndegree) # 6 times 4
ibasis0.all.con = iSpline(t_eval, knots = var.kint.con , degree = ndegree) # 50 times 4

ibasis0 = iSpline(t1, knots = var.kint , degree = ndegree) # 6 times 3
ibasis0.all = iSpline(t_eval, knots = var.kint , degree = ndegree) # 50 times 3

#####
ibasis.con0 = (cbind( rep(1, nrow(ibasis0.con0)), ibasis0.con0))
ibasis.all.con0 = (cbind( rep(1, nrow(ibasis0.all.con0)), ibasis0.all.con0))


ibasis.con = (cbind( rep(1, nrow(ibasis0.con)), ibasis0.con))
ibasis.all.con = (cbind( rep(1, nrow(ibasis0.all.con)), ibasis0.all.con))

ibasis = (cbind( rep(1, nrow(ibasis0)), ibasis0))
ibasis.all = (cbind( rep(1, nrow(ibasis0.all)), ibasis0.all))  

Bmat.u = as.matrix( bdiag(ibasis, ibasis, ibasis, ibasis) ) # 24 times 16 # including intercept # put 0 if want to see zero at t=0
Bmat.all.u = as.matrix( bdiag(ibasis.all, ibasis.all, ibasis.all,ibasis.all) ) # 200 times 16 

#Bmat = as.matrix( bdiag(ibasis, ibasis0.con, ibasis0, ibasis0) ) 
#Bmat.all = as.matrix( bdiag(ibasis.all, ibasis0.all.con, ibasis0.all,ibasis0.all) ) 
#head(df)

# unconstrained fit for comparison
Q.u=0;C.u =0
for (i1 in id.subj){
  df_sub= df[which(df$subject==i1),]
  ni = nrow(df_sub)
  ti = df_sub$time
  Bi = Bmat.u  
  Xmat= as.matrix( data.frame( df_sub[, c('X0','X1','X2','X3')]) )
  yi = as.matrix(df_sub$y)
  wi = diag( 1/ni, ni )
  
  basis_id= which( t1 %in% ti )
  Ui = matrix(0, ncol= ncol(Bmat.u), nrow= nrow(Xmat))
  for (j1 in 1:ni){
    j2= basis_id[j1]
    Ui[j1,] = t(Xmat[j1,]) %*% Bmat.u[c(j2, nrow(ibasis)+j2, 2*nrow(ibasis)+j2, 3*nrow(ibasis)+j2), ]  # dep on p
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
unconst_b2hat.all=uncost_fit.all[(2*length(t_eval)+1):(3*length(t_eval))] 
unconst_b3hat.all=uncost_fit.all[(3*length(t_eval)+1):(4*length(t_eval))] 

Bs_vec_hat_alt=rbind(unconst_b0hat.all, unconst_b1hat.all, unconst_b2hat.all, unconst_b3hat.all)

inc_id1 = which( tau>= range1[1] & tau<= range1[2])
knot_rm11= which( round(ibasis.con[min(inc_id1),],3) ==1   )
knot_rm12= which( round(ibasis.con[max(inc_id1), 1:ncol(ibasis.con)],3) ==0 )
knot_rm0 = c(1:ncol(ibasis))
knot_rm1 = c(knot_rm11, knot_rm12)
knot_rm2 = c(1:ncol(ibasis))
knot_rm3 =  c(1:ncol(ibasis)) 


knot_c0 = c(1:ncol(ibasis))[!(c(1:ncol(ibasis)) %in% knot_rm0)]
knot_c1 = c(1:ncol(ibasis.con))[!(c(1:ncol(ibasis.con)) %in% knot_rm1)]
knot_c2 = c(1:ncol(ibasis0))[!(c(1:ncol(ibasis0)) %in% knot_rm2)]
knot_c3 = c(1:ncol(ibasis0))[!(c(1:ncol(ibasis0)) %in% knot_rm3)]


# constrained
Bmat0 = as.matrix( bdiag(ibasis, ibasis.con0, ibasis, ibasis) ) # 24 times 17 # including intercept # put 0 if want to see zero at t=0
Bmat.all0 = as.matrix( bdiag(ibasis.all, ibasis.all.con0, ibasis.all,ibasis.all) ) # 200 times 17 


Bmat = as.matrix( bdiag(ibasis, ibasis.con, ibasis, ibasis) ) # 24 times 17 # including intercept # put 0 if want to see zero at t=0
Bmat.all = as.matrix( bdiag(ibasis.all, ibasis.all.con, ibasis.all,ibasis.all) ) # 200 times 17 


amat0 = matrix(0, nrow=length(knot_c1), ncol=ncol(Bmat0))
amat0[1, (ncol(ibasis)+2)] = -1
amat0[2, (ncol(ibasis)+3)] = -1
amat0[3, (ncol(ibasis)+4)] = 1
amat0[4, (ncol(ibasis)+5)] = 1
#amat0[5, (ncol(ibasis)+6)] = 1


Q=0;C=0
for (i1 in id.subj){
  df_sub= df[which(df$subject==i1),]
  ni = nrow(df_sub)
  ti = df_sub$time
  Bi = Bmat0  
  Xmat= as.matrix( data.frame( df_sub[, c('X0','X1','X2','X3')]) )
  yi = as.matrix(df_sub$y)
  wi = diag( 1/ni, ni )
  
  basis_id= which( t1 %in% ti )
  Ui = matrix(0, ncol= ncol(Bmat0), nrow= nrow(Xmat))
  for (j1 in 1:ni){
    j2= basis_id[j1]
    Ui[j1,] = t(Xmat[j1,]) %*% Bmat0[c(j2, nrow(ibasis)+j2, 2*nrow(ibasis)+j2, 3*nrow(ibasis)+j2), ]  # dep on p
  }
  Qi = t(Ui)%*%wi%*%Ui
  ci =  t(yi)%*%wi%*%Ui
  
  Q= Q+Qi
  C= C+ci
}

b= rep(10^-10, nrow(amat0)) 
qp.res= qprog(Q, C, amat0, b )
est= Bmat0%*%qp.res$thetahat #const est

const_b0hat.all= (Bmat.all0%*%qp.res$thetahat)[1:length(t_eval)]
const_b1hat.all= (Bmat.all0%*%qp.res$thetahat)[(length(t_eval)+1):(2*length(t_eval))]
const_b2hat.all= (Bmat.all0%*%qp.res$thetahat)[(2*length(t_eval)+1):(3*length(t_eval))]
const_b3hat.all= (Bmat.all0%*%qp.res$thetahat)[(3*length(t_eval)+1):(4*length(t_eval))]

# const_b0hat.all= (Bmat%*%qp.res$thetahat)[1:6] 
# const_b1hat.all= (Bmat%*%qp.res$thetahat)[7:12] 
# const_b2hat.all= (Bmat%*%qp.res$thetahat)[13:18]
# const_b3hat.all= (Bmat%*%qp.res$thetahat)[19:24]


Bs_vec_hat_null = rbind(const_b0hat.all, const_b1hat.all, const_b2hat.all, const_b3hat.all)


#### bootstrap ###
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

amat = matrix(0, nrow=length(knot_c1), ncol=ncol(Bmat))
amat[1, (ncol(ibasis)+2)] = -1
amat[2, (ncol(ibasis)+3)] = -1
amat[3, (ncol(ibasis)+4)] = 1
amat[4, (ncol(ibasis)+5)] = 1
amat[5, (ncol(ibasis)+6)] = 1


# wild bootstrap
boot.list<-c()
for (ii in 1:B) {
  
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
  
  df.boot= data.frame(cbind(subject= uid, time=t_vec, y=y_vec.boot, X1=df[ ,-c(1:3)])) #, X1= X[,1], X2=X[,2], X3=X[,3], X4=X[,4]  )
  names(df.boot)[c(4:7)] <- c('X0','X1','X2','X3')

  
  Q=0;C=0
  for (i1 in id.subj){
    df.boot_sub= df.boot[which(df.boot$subject==i1),]
    ni = nrow(df.boot_sub)
    ti = df.boot_sub$time
    Bi = Bmat  
    Xmat= as.matrix( data.frame( df.boot_sub[, c('X0','X1','X2','X3')]) )
    yi = as.matrix(df.boot_sub$y)
    wi = diag( 1/ni, ni )
    
    basis_id= which( t1 %in% ti )
    Ui = matrix(0, ncol= ncol(Bmat), nrow= nrow(Xmat))
    for (j1 in 1:ni){
      j2= basis_id[j1]
      Ui[j1,] = t(Xmat[j1,]) %*% Bmat[c(j2, nrow(ibasis)+j2, 2*nrow(ibasis)+j2, 3*nrow(ibasis)+j2), ]  # dep on p
    }
    Qi = t(Ui)%*%wi%*%Ui
    ci =  t(yi)%*%wi%*%Ui
    
    Q= Q+Qi
    C= C+ci
  }
  
  
  b= rep(10^-10, nrow(amat)) 
  qp.res.boot= qprog(Q, C, amat, b )
  #est.boot= Bmat%*%qp.res.boot$thetahat #const est
  
  const_b0hat.all.boot= (Bmat.all%*%qp.res.boot$thetahat)[1:length(t_eval)] 
  const_b1hat.all.boot= (Bmat.all%*%qp.res.boot$thetahat)[(length(t_eval)+1):(2*length(t_eval))] 
  const_b2hat.all.boot= (Bmat.all%*%qp.res.boot$thetahat)[(2*length(t_eval)+1):(3*length(t_eval))]
  const_b3hat.all.boot= (Bmat.all%*%qp.res.boot$thetahat)[(3*length(t_eval)+1):(4*length(t_eval))]
  
  Bs_vec_hat_null_boot = rbind(const_b0hat.all.boot, const_b1hat.all.boot, 
                               const_b2hat.all.boot, const_b3hat.all.boot)
  
  plot(t_eval, Bs_vec_hat_null_boot[2,])
  abline(v=4, col="red")
  
  boot.list[[ii]] <- Bs_vec_hat_null_boot
  
  if(ii %% 100 ==0) print(ii)
}


#plot(t_eval, Bs_vec_hat_null_boot[2,])


beta.1.boot<-beta.2.boot <- beta.3.boot <- c()
for(ii in 1:B){
  tmp<- boot.list[[ii]]
  beta.1.boot <- rbind(beta.1.boot, tmp[2,])
  beta.2.boot <- rbind(beta.2.boot, tmp[3,])
  beta.3.boot <- rbind(beta.3.boot, tmp[4,])
}

beta.1.low<- apply(beta.1.boot, 2, quantile, probs = 0.025); beta.1.high<- apply(beta.1.boot, 2, quantile, probs = 0.975)
beta.2.low<- apply(beta.2.boot, 2, quantile, probs = 0.025); beta.2.high<- apply(beta.2.boot, 2, quantile, probs = 0.975)
beta.3.low<- apply(beta.3.boot, 2, quantile, probs = 0.025); beta.3.high<- apply(beta.3.boot, 2, quantile, probs = 0.975)

beta.1.low<- apply(beta.1.boot, 2, quantile, probs = 0.05); beta.1.high<- apply(beta.1.boot, 2, quantile, probs = 0.95)
beta.2.low<- apply(beta.2.boot, 2, quantile, probs = 0.05); beta.2.high<- apply(beta.2.boot, 2, quantile, probs = 0.95)
beta.3.low<- apply(beta.3.boot, 2, quantile, probs = 0.05); beta.3.high<- apply(beta.3.boot, 2, quantile, probs = 0.95)


png(paste0('figures/Rplot_cdystonia_90band_V2.png'), width = 1300*2/3, height = 350*2/3)
par(mfrow=c(1,3))
par(mar=c(5,5,2,5))
#par(mar=c(4.5,6,2,1)+0.1)
# plot(t_eval, Bs_vec_hat_alt[1,], type="l", xlab="week", ylab=expression(beta[0](week)),lwd=2, main="Intercept", ylim=c(35,50), col="blue", lty=2) #ylim=range(c(Bs_vec_hat_alt[1,],Bs_vec_hat_null[1,])))
# lines(t_eval, Bs_vec_hat_null[1,], col="blue", lwd=2, lty=1)
#abline(h=0, lty=3)
#abline(v=domain_con, col="red", lty=2)
cex.set = 1.75

plot(t_eval, Bs_vec_hat_alt[2,], type="l", xlab="week", ylab=expression(beta[1](week)), lwd=2, main="Treatment", ylim=c(-10,7), col="black", lty=2, cex.main=cex.set, cex.lab=cex.set, cex.sub=cex.set, cex.axis = cex.set ) #ylim=range(c(Bs_vec_hat_alt[2,],Bs_vec_hat_null[2,])))
# lines(t_eval, beta.1.low, col="purple")
# lines(t_eval, beta.1.high, col="purple")
polygon(c(t_eval, rev(t_eval)), c(beta.1.high, rev(beta.1.low)), col = rgb(0.7, 0.7, 0.7, 0.5), border = NA)
lines(t_eval, Bs_vec_hat_null[2,], col="black", lwd=2, lty=1)
abline(v=4, lty=3, col="red", lwd=2)
points(c(0,var.kint.con0,16), -10 + rep( 0.1,length=4) + 0.05 , pch=4, cex=2, col="blue", lwd=2 )
grid(lwd = 1.5)
abline(h=0, lty=3)
#abline(v=domain_con, col="red", lty=2)

plot(t_eval, Bs_vec_hat_alt[4,], type="l", xlab="week", ylab=expression(beta[2](week)), lwd=2, main="Age",  ylim=c(-0.15, 0.15), col="black", lty=2, cex.main=cex.set, cex.lab=cex.set, cex.axis = cex.set) #ylim=range(c(Bs_vec_hat_alt[4,],Bs_vec_hat_null[4,])))
polygon(c(t_eval, rev(t_eval)), c(beta.3.high, rev(beta.3.low)), col = rgb(0.7, 0.7, 0.7, 0.5), border = NA)
lines(t_eval, Bs_vec_hat_null[4,], col="black", lwd=2, lty=1)
#lines(t_eval, beta.3.low, col="purple")
#lines(t_eval, beta.3.high, col="purple")
points(c(0,var.kint,16), -0.15 + rep( 0,length=3) + 0.0005 , pch=4, cex=2, col="blue", lwd=2 )
grid(lwd = 1.5)
abline(h=0, lty=3)

plot(t_eval, Bs_vec_hat_alt[3,], type="l", xlab="week", ylab=expression(beta[3](week)), lwd=2, main="Sex",  ylim=c(-7,5), col="black", lty=2, cex.main=cex.set, cex.lab=cex.set, cex.axis = cex.set) #ylim=range(c(Bs_vec_hat_alt[3,],Bs_vec_hat_null[3,])))
polygon(c(t_eval, rev(t_eval)), c(beta.2.high, rev(beta.2.low)), col = rgb(0.7, 0.7, 0.7, 0.5), border = NA)
lines(t_eval, Bs_vec_hat_null[3,], col="black", lwd=2, lty=1)
#lines(t_eval, beta.2.low, col="purple")1es#
#lines(t_eval, beta.2.high, col="purple")
points(c(0,var.kint,16), -7 + rep( 0,length=3) + 0.1 , pch=4, cex=2, col="blue", lwd=2 )
grid(lwd = 1.5)
abline(h=0, lty=3)
dev.off()

mean(df[ df$time==0 & df$X1==0, "y"])
mean(df[ df$time==0 & df$X1==1, "y"])




#######
library(tvem)
model3 <- tvem(data=df,
               formula=y~X1+X2+X3,               
               #formula=y~X2, 
               id=subject,
               time=time,
               num_knots= 2)
print(model3)
plot(model3)

names(model3)
model3$grid_fitted_coefficients

tvem_time= model3$time_grid
tvem_res2= model3$grid_fitted_coefficients$X2[,1]
tvem_res1= model3$grid_fitted_coefficients$X1[,1]
