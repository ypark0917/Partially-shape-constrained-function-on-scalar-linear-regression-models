# 05/23/23
# use t_eval for defining i-splines - .all output
# Modify part for removal of knots according to domain_con

# 01/08/24 choice of knots
op.knot = function(df, nknot.set=c(1:8) , int.set1=seq(0.1, 0.4, length=10), int.set2=seq(0.6, 0.9, length=10), ndegree=1, cv=10){
  
  t1= sort(unique(df$time))  
  tau=sort(unique(df$time)) 
  cv.fold = createFolds(1:n, k = 10, list = TRUE, returnTrain = FALSE)
  
  if(is.numeric(nknot.set)){
    cv.res <- c()
    for(v in 1:cv){
      cv.err<-c()
      cv.id = unlist(cv.fold[[v]])
      df.train = df[!(df$subject %in% cv.fold[[v]]), ]
      df.test = df[which(df$subject %in% cv.fold[[v]]), ]
      nsubj= length(unique(df.train$subject))
      subj.df = unique(df.train$subject)
      
      for (nk in nknot.set){
        var.knots = seq(min(df.train$time), max(df.train$time), length= nk+2 ) 
        var.kint=var.knots[-c(1,length(var.knots))]
        
        ibasis0 = iSpline(t1, knots = var.kint , degree = ndegree) 
        ibasis0.all = iSpline(t_eval, knots = var.kint , degree = ndegree) 
        
        ibasis = (cbind( rep(1, nrow(ibasis0)), ibasis0))
        ibasis.all = (cbind( rep(1, nrow(ibasis0.all)), ibasis0.all))
        
        Bmat = as.matrix( bdiag(ibasis, ibasis, ibasis, ibasis) ) 
        Bmat.all = as.matrix( bdiag(ibasis.all, ibasis.all, ibasis.all, ibasis.all) ) # to evaluate from t_eval grids
        
        Q=0;C=0
        for (i1 in subj.df){
          df.train_sub= df.train[which(df.train$subject==i1),]
          ni = nrow(df.train_sub)
          ti = df.train_sub$time
          Bi = Bmat  
          Xmat= as.matrix( data.frame( df.train_sub[, c('X1','X2','X3','X4')]) )
          yi = as.matrix(df.train_sub$y)
          wi = diag( 1/ni, ni )
          
          basis_id= which( t1 %in% ti )
          Ui = matrix(0, ncol= ncol(Bmat), nrow= nrow(Xmat))
          for (j1 in 1:ni){
            j2= basis_id[j1]
            Ui[j1,] = t(Xmat[j1,]) %*% Bmat[c(j2, nrow(ibasis)+j2, 2*nrow(ibasis)+j2,3*nrow(ibasis)+j2),]  # dep on p
          }
          Qi = t(Ui)%*%wi%*%Ui
          ci =  t(yi)%*%wi%*%Ui
          
          Q= Q+Qi
          C= C+ci
        }
        
        theta_u= solve(Q)%*% t(C) #unconstrained est 
        uncost_fit = Bmat %*% theta_u
        unconst_b1hat=uncost_fit[1:length(tau)]
        unconst_b2hat=uncost_fit[(length(tau)+1):(2*length(tau))]
        unconst_b3hat=uncost_fit[(2*length(tau)+1):(3*length(tau))] 
        unconst_b4hat=uncost_fit[(3*length(tau)+1):(4*length(tau))]
        
        test.err <- c()
        for (i2 in cv.id){
          tmp <- df.test[df.test$subject == i2 , ]
          Y.fit.test <- tmp[1,"X1"]*unconst_b1hat + tmp[1,"X2"]*unconst_b2hat + tmp[1,"X3"]*unconst_b3hat + tmp[1,"X4"]*unconst_b4hat 
          test.err<- c(test.err, (sum(Y.fit.test[which(t1 %in% tmp$time)] - tmp$y)^2))
        }
        
        cv.err <- c(cv.err, mean(test.err)) 
      }
      cv.res <- rbind(cv.res, cv.err)
    }
    
    nknot.op <- nknot.set[which.min(colMeans(cv.res))]
    var.knots = seq(0, 1, length= nknot.op+2 ) 
    var.kint.op =var.knots[-c(1,length(var.knots))]
  }
  if(!is.numeric(nknot.set)){
    cv.res <- array(dim=c(length(int.set1), length(int.set2), cv ))
    for(v in 1:cv){
      cv.id = unlist(cv.fold[[v]])
      df.train = df[!(df$subject %in% cv.fold[[v]]), ]
      df.test = df[which(df$subject %in% cv.fold[[v]]), ]
      nsubj= length(unique(df.train$subject))
      subj.df = unique(df.train$subject)
      
      cv.err2 <- c()
      #for (nk in nknot.set){
      for (k1 in int.set1){  
        
        cv.err <-c()
        for(k2 in int.set2){
          
          #var.knots = seq(min(df.train$time), max(df.train$time), length= nk+2 ) 
          #var.kint=var.knots[-c(1,length(var.knots))]
          var.kint = c(k1, 0.5, k2)
          ibasis0 = iSpline(t1, knots = var.kint , degree = ndegree) 
          ibasis0.all = iSpline(t_eval, knots = var.kint , degree = ndegree) 
          
          ibasis = (cbind( rep(1, nrow(ibasis0)), ibasis0))
          ibasis.all = (cbind( rep(1, nrow(ibasis0.all)), ibasis0.all))
          
          Bmat = as.matrix( bdiag(ibasis, ibasis, ibasis, ibasis) ) 
          Bmat.all = as.matrix( bdiag(ibasis.all, ibasis.all, ibasis.all, ibasis.all) ) # to evaluate from t_eval grids
          
          Q=0;C=0
          for (i1 in subj.df){
            df.train_sub= df.train[which(df.train$subject==i1),]
            ni = nrow(df.train_sub)
            ti = df.train_sub$time
            Bi = Bmat  
            Xmat= as.matrix( data.frame( df.train_sub[, c('X1','X2','X3','X4')]) )
            yi = as.matrix(df.train_sub$y)
            wi = diag( 1/ni, ni )
            
            basis_id= which( t1 %in% ti )
            Ui = matrix(0, ncol= ncol(Bmat), nrow= nrow(Xmat))
            for (j1 in 1:ni){
              j2= basis_id[j1]
              Ui[j1,] = t(Xmat[j1,]) %*% Bmat[c(j2, nrow(ibasis)+j2, 2*nrow(ibasis)+j2,3*nrow(ibasis)+j2),]  # dep on p
            }
            Qi = t(Ui)%*%wi%*%Ui
            ci =  t(yi)%*%wi%*%Ui
            
            Q= Q+Qi
            C= C+ci
          }
          
          theta_u= solve(Q)%*% t(C) #unconstrained est 
          uncost_fit = Bmat %*% theta_u
          unconst_b1hat=uncost_fit[1:length(tau)]
          unconst_b2hat=uncost_fit[(length(tau)+1):(2*length(tau))]
          unconst_b3hat=uncost_fit[(2*length(tau)+1):(3*length(tau))] 
          unconst_b4hat=uncost_fit[(3*length(tau)+1):(4*length(tau))]
          
          test.err <- c()
          for (i2 in cv.id){
            tmp <- df.test[df.test$subject == i2 , ]
            Y.fit.test <- tmp[1,"X1"]*unconst_b1hat + tmp[1,"X2"]*unconst_b2hat + tmp[1,"X3"]*unconst_b3hat + tmp[1,"X4"]*unconst_b4hat 
            test.err<- c(test.err, (sum(Y.fit.test[which(t1 %in% tmp$time)] - tmp$y)^2))
          }
          
          cv.err <- c(cv.err, mean(test.err)) 
        }
        
        cv.err2<- rbind(cv.err2, cv.err)
      }
      cv.res[ , , v] <- cv.err2
    }
    
    tmp<-0
    for(i in 1:10){
      tmp<- tmp+cv.res[,,i]
    }
    
    kint.gr = which(tmp == min(tmp), arr.ind = TRUE)
    #nknot.op <- nknot.set[which.min(rowMeans(cv.res))]
    var.kint.op <- c(int.set1[kint.gr[1]], 0.5, int.set2[kint.gr[2]])
  }
  
  #return(nknot.op)
  return(var.kint.op)
}

# 01/24/24 amomg two or three internal knots
###################################
op.knot.V2 = function(df , int.set.list=NA, ndegree=1, cv=10){
  
  t1= sort(unique(df$time))  
  tau=sort(unique(df$time)) 
  cv.fold = createFolds(1:n, k = 10, list = TRUE, returnTrain = FALSE)
  
  
  cv.res <- c()
  for(v in 1:cv){
    cv.id = unlist(cv.fold[[v]])
    df.train = df[!(df$subject %in% cv.fold[[v]]), ]
    df.test = df[which(df$subject %in% cv.fold[[v]]), ]
    nsubj= length(unique(df.train$subject))
    subj.df = unique(df.train$subject)
    
    cv.err <- c()
    #for (nk in nknot.set){
    for (k1 in 1:length(int.set.list)){  
      
      
      #var.knots = seq(min(df.train$time), max(df.train$time), length= nk+2 ) 
      #var.kint=var.knots[-c(1,length(var.knots))]
      var.kint = int.set.list[[k1]]
      ibasis0 = iSpline(t1, knots = var.kint , degree = ndegree) 
      ibasis0.all = iSpline(t_eval, knots = var.kint , degree = ndegree) 
      
      ibasis = (cbind( rep(1, nrow(ibasis0)), ibasis0))
      ibasis.all = (cbind( rep(1, nrow(ibasis0.all)), ibasis0.all))
      
      Bmat = as.matrix( bdiag(ibasis, ibasis, ibasis, ibasis) ) 
      Bmat.all = as.matrix( bdiag(ibasis.all, ibasis.all, ibasis.all, ibasis.all) ) # to evaluate from t_eval grids
      
      Q=0;C=0
      for (i1 in subj.df){
        df.train_sub= df.train[which(df.train$subject==i1),]
        ni = nrow(df.train_sub)
        ti = df.train_sub$time
        Bi = Bmat  
        Xmat= as.matrix( data.frame( df.train_sub[, c('X1','X2','X3','X4')]) )
        yi = as.matrix(df.train_sub$y)
        wi = diag( 1/ni, ni )
        
        basis_id= which( t1 %in% ti )
        Ui = matrix(0, ncol= ncol(Bmat), nrow= nrow(Xmat))
        for (j1 in 1:ni){
          j2= basis_id[j1]
          Ui[j1,] = t(Xmat[j1,]) %*% Bmat[c(j2, nrow(ibasis)+j2, 2*nrow(ibasis)+j2,3*nrow(ibasis)+j2),]  # dep on p
        }
        Qi = t(Ui)%*%wi%*%Ui
        ci =  t(yi)%*%wi%*%Ui
        
        Q= Q+Qi
        C= C+ci
      }
      
      if(qr(Q)$rank != ncol(Q)){
        Q <- as.matrix(nearPD(Q)$mat)
      }
      theta_u= solve(Q)%*% t(C) #unconstrained est 
      uncost_fit = Bmat %*% theta_u
      unconst_b1hat=uncost_fit[1:length(tau)]
      unconst_b2hat=uncost_fit[(length(tau)+1):(2*length(tau))]
      unconst_b3hat=uncost_fit[(2*length(tau)+1):(3*length(tau))] 
      unconst_b4hat=uncost_fit[(3*length(tau)+1):(4*length(tau))]
      
      test.err <- c()
      for (i2 in cv.id){
        tmp <- df.test[df.test$subject == i2 , ]
        Y.fit.test <- tmp[1,"X1"]*unconst_b1hat + tmp[1,"X2"]*unconst_b2hat + tmp[1,"X3"]*unconst_b3hat + tmp[1,"X4"]*unconst_b4hat 
        test.err<- c(test.err, (sum(Y.fit.test[which(t1 %in% tmp$time)] - tmp$y)^2))
      }
      
      cv.err <- c(cv.err, mean(test.err)) 
    }
    cv.res<- rbind(cv.res, cv.err)
  }
  
  #kint.gr = which(tmp == min(tmp), arr.ind = TRUE)
  #nknot.op <- nknot.set[which.min(rowMeans(cv.res))]
  var.kint.op <- int.set.list[[which.min(colMeans(cv.res))]]
  
  
  #return(nknot.op)
  return(var.kint.op)
}


op.knot.con.V2 = function(df, int.set.list=NA , const=c(shape.const.1, shape.const.2,NA,NA), ndegree=1, cv=10){
  
  t1= sort(unique(df$time))  
  tau=sort(unique(df$time)) 
  cv.fold = createFolds(1:n, k = 10, list = TRUE, returnTrain = FALSE)
  
  cv.res <- c()
  for(v in 1:cv){
    cv.id = unlist(cv.fold[[v]])
    df.train = df[!(df$subject %in% cv.fold[[v]]), ]
    df.test = df[which(df$subject %in% cv.fold[[v]]), ]
    nsubj= length(unique(df.train$subject))
    subj.df = unique(df.train$subject)
    
    
    cv.err <- c()
    for (k1 in 1:length(int.set.list)){  
      
      var.kint = int.set.list[[k1]]
      ibasis0 = iSpline(t1, knots = var.kint , degree = ndegree) 
      ibasis0.all = iSpline(t_eval, knots = var.kint , degree = ndegree) 
      
      ibasis = (cbind( rep(1, nrow(ibasis0)), ibasis0))
      ibasis.all = (cbind( rep(1, nrow(ibasis0.all)), ibasis0.all))
      
      Bmat = as.matrix( bdiag(ibasis, ibasis, ibasis, ibasis) ) 
      Bmat.all = as.matrix( bdiag(ibasis.all, ibasis.all, ibasis.all, ibasis.all) ) # to evaluate from t_eval grids
      
      Q=0;C=0
      for (i1 in subj.df){
        df.train_sub= df.train[which(df.train$subject==i1),]
        ni = nrow(df.train_sub)
        ti = df.train_sub$time
        Bi = Bmat  
        Xmat= as.matrix( data.frame( df.train_sub[, c('X1','X2','X3','X4')]) )
        yi = as.matrix(df.train_sub$y)
        wi = diag( 1/ni, ni )
        
        basis_id= which( t1 %in% ti )
        Ui = matrix(0, ncol= ncol(Bmat), nrow= nrow(Xmat))
        for (j1 in 1:ni){
          j2= basis_id[j1]
          Ui[j1,] = t(Xmat[j1,]) %*% Bmat[c(j2, nrow(ibasis)+j2, 2*nrow(ibasis)+j2,3*nrow(ibasis)+j2),]  # dep on p
        }
        Qi = t(Ui)%*%wi%*%Ui
        ci =  t(yi)%*%wi%*%Ui
        
        Q= Q+Qi
        C= C+ci
      }
      
      if(qr(Q)$rank != ncol(Q)){
        Q <- as.matrix(nearPD(Q)$mat)
      }
      knot_c1=c(1:ncol(ibasis.all))
      knot_c2=c()
      knot_c3=c(1:ncol(ibasis.all))
      knot_c4=c()
      
      amat = matrix(0, nrow=length(knot_c1)+length(knot_c2)+length(knot_c3)+length(knot_c4),ncol=ncol(Bmat))
      # p=1
      if(!is.na(const[1])){
        if(const[1] == "increasing"){
          for(i1 in 1:length(knot_c1)){
            amat[i1, knot_c1[i1]]=1
          }
        }else if(const[1] == "decreasing"){
          for(i1 in 1:length(knot_c1)){
            amat[i1, knot_c1[i1]]=-1
          }
        } 
      }
      
      # p=2
      if(!is.na(const[2])){
        if(const[2] == "increasing"){
          for(i2 in 1:length(knot_c2)){
            amat[length(knot_c1)+i2, ncol(ibasis.all) + knot_c2[i2]]=1
          }
        }else if(const[2] == "decreasing"){
          for(i2 in 1:length(knot_c2)){
            amat[length(knot_c1)+i2,  ncol(ibasis.all) + knot_c2[i2]]=-1
          }
        }
      }
      
      # p=3
      if(!is.na(const[3])){
        if(const[3] == "increasing"){
          for(i3 in 1:length(knot_c3)){
            amat[length(knot_c1)+length(knot_c2)+i3, 2*ncol(ibasis.all) + knot_c3[i3]]=1
          }
        }else if(const[3] == "decreasing"){
          for(i3 in 1:length(knot_c3)){
            amat[length(knot_c1)+length(knot_c2)+i3,  2*ncol(ibasis.all) + knot_c3[i3]]=-1
          }
        }
      }
      
      # p=4
      if(!is.na(const[4])){
        if(const[4] == "increasing"){
          for(i4 in 1:length(knot_c4)){
            amat[length(knot_c1)+length(knot_c2)+length(knot_c3)+i4, 3*ncol(ibasis.all) + knot_c4[i4]]=1
          }
        }else if(const[4] == "decreasing"){
          for(i4 in 1:length(knot_c4)){
            amat[length(knot_c1)+length(knot_c2)+length(knot_c3)+i4,  3*ncol(ibasis.all) + knot_c4[i4]]=-1
          }
        }
      }
      
      b= rep(10^-10, nrow(amat)) 
      qp.res= qprog(Q, C, amat, b )
      
      const_b1hat= (Bmat%*%qp.res$thetahat)[1:length(tau)]
      const_b2hat= (Bmat%*%qp.res$thetahat)[(length(tau)+1):(2*length(tau))] 
      const_b3hat= (Bmat%*%qp.res$thetahat)[(2*length(tau)+1):(3*length(tau))] 
      const_b4hat= (Bmat%*%qp.res$thetahat)[(3*length(tau)+1):(4*length(tau))] 
      const_res = data.frame( t=tau, b1hat=const_b1hat, b2hat=const_b2hat, b3hat=const_b3hat, b4hat=const_b4hat )
      
      
      test.err <- c()
      for (i2 in cv.id){
        tmp <- df.test[df.test$subject == i2 , ]
        Y.fit.test <- tmp[1,"X1"]*const_b1hat + tmp[1,"X2"]*const_b2hat + tmp[1,"X3"]*const_b3hat + tmp[1,"X4"]*const_b4hat 
        test.err<- c(test.err, (sum(Y.fit.test[which(t1 %in% tmp$time)] - tmp$y)^2))
      }
      
      cv.err <- c(cv.err, mean(test.err)) 
    }
    cv.res <- rbind(cv.res, cv.err)
  }
  
  var.kint.op <- int.set.list[[which.min(colMeans(cv.res))]]
  
  return(var.kint.op)
}

#########################################################
# optimal knot selection for partial data
op.knot.re = function(df, int.set, ndegree=1, cv=10){
  
  t1= sort(unique(df$time))  
  tau=sort(unique(df$time)) 
  cv.fold = createFolds(1:n, k = 10, list = TRUE, returnTrain = FALSE)
  
  cv.res <- c()
  for(v in 1:cv){
    cv.id = unlist(cv.fold[[v]])
    df.train = df[!(df$subject %in% cv.fold[[v]]), ]
    df.test = df[which(df$subject %in% cv.fold[[v]]), ]
    nsubj= length(unique(df.train$subject))
    subj.df = unique(df.train$subject)
    
    cv.err <- c()
    #for (nk in nknot.set){
    for (k1 in int.set){  
      
      #var.knots = seq(min(df.train$time), max(df.train$time), length= nk+2 ) 
      #var.kint=var.knots[-c(1,length(var.knots))]
      var.kint = k1
      ibasis0 = iSpline(t1, knots = var.kint , degree = ndegree) 
      ibasis0.all = iSpline(t_eval, knots = var.kint , degree = ndegree) 
      
      ibasis = (cbind( rep(1, nrow(ibasis0)), ibasis0))
      ibasis.all = (cbind( rep(1, nrow(ibasis0.all)), ibasis0.all))
      
      Bmat = as.matrix( bdiag(ibasis, ibasis, ibasis, ibasis) ) 
      Bmat.all = as.matrix( bdiag(ibasis.all, ibasis.all, ibasis.all, ibasis.all) ) # to evaluate from t_eval grids
      
      Q=0;C=0
      for (i1 in subj.df){
        df.train_sub= df.train[which(df.train$subject==i1),]
        ni = nrow(df.train_sub)
        ti = df.train_sub$time
        Bi = Bmat  
        Xmat= as.matrix( data.frame( df.train_sub[, c('X1','X2','X3','X4')]) )
        yi = as.matrix(df.train_sub$y)
        wi = diag( 1/ni, ni )
        
        basis_id= which( t1 %in% ti )
        Ui = matrix(0, ncol= ncol(Bmat), nrow= nrow(Xmat))
        for (j1 in 1:ni){
          j2= basis_id[j1]
          Ui[j1,] = t(Xmat[j1,]) %*% Bmat[c(j2, nrow(ibasis)+j2, 2*nrow(ibasis)+j2,3*nrow(ibasis)+j2),]  # dep on p
        }
        Qi = t(Ui)%*%wi%*%Ui
        ci =  t(yi)%*%wi%*%Ui
        
        Q= Q+Qi
        C= C+ci
      }
      
      if(qr(Q)$rank != ncol(Q)){
        Q <- as.matrix(nearPD(Q)$mat)
      }
      theta_u= solve(Q)%*% t(C) #unconstrained est 
      uncost_fit = Bmat %*% theta_u
      unconst_b1hat=uncost_fit[1:length(tau)]
      unconst_b2hat=uncost_fit[(length(tau)+1):(2*length(tau))]
      unconst_b3hat=uncost_fit[(2*length(tau)+1):(3*length(tau))] 
      unconst_b4hat=uncost_fit[(3*length(tau)+1):(4*length(tau))]
      
      test.err <- c()
      for (i2 in cv.id){
        tmp <- df.test[df.test$subject == i2 , ]
        Y.fit.test <- tmp[1,"X1"]*unconst_b1hat + tmp[1,"X2"]*unconst_b2hat + tmp[1,"X3"]*unconst_b3hat + tmp[1,"X4"]*unconst_b4hat 
        test.err<- c(test.err, (sum(Y.fit.test[which(t1 %in% tmp$time)] - tmp$y)^2))
      }
      
      cv.err <- c(cv.err, mean(test.err)) 
    }
    cv.res <- cbind(cv.res, cv.err)
  }
  
  kint.gr = which.min(rowMeans(cv.res))
  #nknot.op <- nknot.set[which.min(rowMeans(cv.res))]
  var.kint.op <- int.set[kint.gr]
  
  #return(nknot.op)
  return(var.kint.op)
}




# 06/29/23 modify for var.kint and removal of knots if round(eval, 2) # same number of knots
# 11/15/23 modify for generalization
# 03/25/24 modify for beta1 and beta3
vcm_irr = function(df, nknot, t_eval=t_eval, ndegree=2, var.kint=NA,
                   const=c("increasing",NA,NA,NA,NA),range1=c(0,1), range2=NA, range3=NA, range4=NA){
  
  nsubj= length(unique(df$subject))
  t1= sort(unique(df$time))  
  tau=sort(unique(df$time)) 
  
  if(is.numeric(var.kint)){
    var.kint=var.kint
  } else{
    var.knots = seq(min(df$time), max(df$time), length= nknot+2 ) 
    var.kint=var.knots[-c(1,length(var.knots))]
  }
  ibasis0 = iSpline(t1, knots = var.kint , degree = ndegree) #basis ftn for one param
  ibasis0.all = iSpline(t_eval, knots = var.kint , degree = ndegree) #basis ftn for one param
  
  ibasis = (cbind( rep(1, nrow(ibasis0)), ibasis0))
  ibasis.all = (cbind( rep(1, nrow(ibasis0.all)), ibasis0.all))#[,-6]
  
  
  ###### need to fix this for the general form
  #### depeding on p # 06/23/23
  Bmat = as.matrix( bdiag(ibasis, ibasis, ibasis, ibasis) ) 
  Bmat.all = as.matrix( bdiag(ibasis.all, ibasis.all, ibasis.all, ibasis.all) ) # to evaluate from t_eval grids
  
  #matplot(ibasis0, type='l')
  
  #head(df)
  Q=0;C=0
  for (i1 in 1:nsubj){
    df_sub= df[which(df$subject==i1),]
    ni = nrow(df_sub)
    ti = df_sub$time
    Bi = Bmat  
    Xmat= as.matrix( data.frame( df_sub[, c('X1','X2','X3','X4')]) )
    yi = as.matrix(df_sub$y)
    wi = diag( 1/ni, ni )
    
    basis_id= which( t1 %in% ti )
    Ui = matrix(0, ncol= ncol(Bmat), nrow= nrow(Xmat))
    for (j1 in 1:ni){
      j2= basis_id[j1]
      Ui[j1,] = t(Xmat[j1,]) %*% Bmat[c(j2, nrow(ibasis)+j2, 2*nrow(ibasis)+j2,3*nrow(ibasis)+j2),]  # dep on p
    }
    Qi = t(Ui)%*%wi%*%Ui
    ci =  t(yi)%*%wi%*%Ui
    
    Q= Q+Qi
    C= C+ci
  }
  
  if(qr(Q)$rank != ncol(Q)){
    Q <- as.matrix(nearPD(Q)$mat)
  }
  theta_u= solve(Q)%*% t(C) #unconstrained est 
  
  ###########
  # p=1
  if (is.numeric(range1)){ 
    inc_id1 = which( tau> range1[1] & tau<range1[2])
    #knot_rm1= which( ibasis0[min(inc_id),] ==1   )
    #knot_rm2= which( ibasis0[max(inc_id), 1:ncol(ibasis0)] ==0 )
    
    knot_rm11= which( round(ibasis[min(inc_id1),],3) ==1   )
    knot_rm12= which( round(ibasis[max(inc_id1), 1:ncol(ibasis)],3) ==0 )
    knot_rm1 = c(knot_rm11, knot_rm12)
  } else{
    knot_rm1 =c(1:ncol(ibasis))
  }
  
  #p=2
  if (is.numeric(range2)){ 
    inc_id2 = which( tau> range2[1] & tau<range2[2])
    knot_rm21= which( round(ibasis[min(inc_id2),],3) ==1   )
    knot_rm22= which( round(ibasis[max(inc_id2), 1:ncol(ibasis)],3) ==0 )
    knot_rm2 =  c(knot_rm21, knot_rm22)
  } else{
    #knot_rm12 =vector()
    knot_rm2 = c(1:ncol(ibasis)) 
  }  
  
  # p=3
  if (is.numeric(range3)){ 
    inc_id3 = which( tau> range3[1] & tau<range3[2])
    knot_rm31= which( round(ibasis[min(inc_id3),],3) ==1   )
    knot_rm32= which( round(ibasis[max(inc_id3), 1:ncol(ibasis)],3) ==0 )
    knot_rm3 =  c(knot_rm31, knot_rm32)
  } else{
    knot_rm3 =  c(1:ncol(ibasis)) 
  }  
  
  # p=4
  if (is.numeric(range4)){ 
    inc_id4 = which( tau> range4[1] & tau<range4[2])
    knot_rm41= which( round(ibasis[min(inc_id4),],3) ==1   )
    knot_rm42= which( round(ibasis[max(inc_id4), 1:ncol(ibasis)],3) ==0 )
    knot_rm4 =  c(knot_rm41, knot_rm42)
  } else{
    knot_rm4 =  c(1:ncol(ibasis)) 
  }  
  
  #knot_rm= c(knot_rm1, knot_rm2, knot_rm3, knot_rm4)
  knot_c1 = c(1:ncol(ibasis))[!(c(1:ncol(ibasis)) %in% knot_rm1)]
  knot_c2 = c(1:ncol(ibasis))[!(c(1:ncol(ibasis)) %in% knot_rm2)]
  knot_c3 = c(1:ncol(ibasis))[!(c(1:ncol(ibasis)) %in% knot_rm3)]
  knot_c4 = c(1:ncol(ibasis))[!(c(1:ncol(ibasis)) %in% knot_rm4)]
  
  
  # amat0= diag( 1, ncol(Bmat) )
  # if (length(knot_rm)>1){
  #   amat= amat0[-c(knot_rm),]
  # } else{
  #   amat= amat0 
  # }
  amat = matrix(0, nrow=length(knot_c1)+length(knot_c2)+length(knot_c3)+length(knot_c4),ncol=ncol(Bmat))
  # p=1
  if(!is.na(const[1])){
    if(const[1] == "increasing"){
      for(i1 in 1:length(knot_c1)){
        amat[i1, knot_c1[i1]]=1
      }
    }else if(const[1] == "decreasing"){
      for(i1 in 1:length(knot_c1)){
        amat[i1, knot_c1[i1]]=-1
      }
    } 
  }
  
  # p=2
  if(!is.na(const[2])){
  if(const[2] == "increasing"){
    for(i2 in 1:length(knot_c2)){
      amat[length(knot_c1)+i2, ncol(ibasis) + knot_c2[i2]]=1
    }
  }else if(const[2] == "decreasing"){
    for(i2 in 1:length(knot_c2)){
      amat[length(knot_c1)+i2,  ncol(ibasis) + knot_c2[i2]]=-1
    }
  }
  }
  
  # p=3
  if(!is.na(const[3])){
    if(const[3] == "increasing"){
      for(i3 in 1:length(knot_c3)){
        amat[length(knot_c1)+length(knot_c2)+i3, 2*ncol(ibasis) + knot_c3[i3]]=1
      }
    }else if(const[3] == "decreasing"){
      for(i3 in 1:length(knot_c3)){
        amat[length(knot_c1)+length(knot_c2)+i3,  2*ncol(ibasis) + knot_c3[i3]]=-1
      }
    }
  }
  
  # p=4
  if(!is.na(const[4])){
    if(const[4] == "increasing"){
      for(i4 in 1:length(knot_c4)){
        amat[length(knot_c1)+length(knot_c2)+length(knot_c3)+i4, 3*ncol(ibasis) + knot_c4[i4]]=1
      }
    }else if(const[4] == "decreasing"){
      for(i4 in 1:length(knot_c4)){
        amat[length(knot_c1)+length(knot_c2)+length(knot_c3)+i4,  3*ncol(ibasis) + knot_c4[i4]]=-1
      }
    }
  }
  ###########
  b= rep(10^-10, nrow(amat)) 
  qp.res= qprog(Q, C, amat, b )
  
  # if(const=="decreasing"){
  #   qp.res= qprog(Q, C, -amat, b ) # for decreasing
  # } else if (const=="increasing"){
  #   qp.res= qprog(Q, C, amat, b ) # for increasing
  # }
  # 
  
  est= Bmat%*%qp.res$thetahat #const est
  
  const_b1hat= (Bmat%*%qp.res$thetahat)[1:length(tau)]
  const_b2hat= (Bmat%*%qp.res$thetahat)[(length(tau)+1):(2*length(tau))] 
  const_b3hat= (Bmat%*%qp.res$thetahat)[(2*length(tau)+1):(3*length(tau))] 
  const_b4hat= (Bmat%*%qp.res$thetahat)[(3*length(tau)+1):(4*length(tau))] 
  const_res = data.frame( t=tau, b1hat=const_b1hat, b2hat=const_b2hat, b3hat=const_b3hat, b4hat=const_b4hat )
  
  const_b1hat.all= (Bmat.all%*%qp.res$thetahat)[1:length(t_eval)]
  const_b2hat.all= (Bmat.all%*%qp.res$thetahat)[(length(t_eval)+1):(2*length(t_eval))] 
  const_b3hat.all= (Bmat.all%*%qp.res$thetahat)[(2*length(t_eval)+1):(3*length(t_eval))] 
  const_b4hat.all= (Bmat.all%*%qp.res$thetahat)[(3*length(t_eval)+1):(4*length(t_eval))] 
  const_res.teval = data.frame( t=t_eval, b1hat=const_b1hat.all, b2hat=const_b2hat.all, b3hat=const_b3hat.all, b4hat=const_b4hat.all )
  
  #uncost est
  uncost_fit = Bmat %*% theta_u
  unconst_b1hat=uncost_fit[1:length(tau)]
  unconst_b2hat=uncost_fit[(length(tau)+1):(2*length(tau))]
  unconst_b3hat=uncost_fit[(2*length(tau)+1):(3*length(tau))] 
  unconst_b4hat=uncost_fit[(3*length(tau)+1):(4*length(tau))] 
  
  unconst_res= data.frame( t=tau, b1hat=unconst_b1hat, b2hat=unconst_b2hat,
                           b3hat=unconst_b3hat, b4hat=unconst_b4hat)
  
  uncost_fit.all = Bmat.all %*% theta_u
  unconst_b1hat.all=uncost_fit.all[1:length(t_eval)]
  unconst_b2hat.all=uncost_fit.all[(length(t_eval)+1):(2*length(t_eval))] 
  unconst_b3hat.all=uncost_fit.all[(2*length(t_eval)+1):(3*length(t_eval))] 
  unconst_b4hat.all=uncost_fit.all[(3*length(t_eval)+1):(4*length(t_eval))]
  unconst_res.teval= data.frame( t=t_eval, b1hat=unconst_b1hat.all, b2hat=unconst_b2hat.all,
                                 b3hat=unconst_b3hat.all, b4hat=unconst_b4hat.all)
  
  res= list( const= const_res, const_teval=const_res.teval, 
             unconst= unconst_res, unconst_teval=unconst_res.teval)
  
  return(res)
}

# to add constrained only for beta1 and beta3
vcm_irr_con = function(df, nknot, t_eval=t_eval, ndegree=2, var.kint.con=NA, var.kint=NA,
                       const=c("increasing",NA,NA,NA,NA),range1=c(0,1), range2=NA, range3=NA, range4=NA){
  
  nsubj= length(unique(df$subject))
  t1= sort(unique(df$time))  
  tau=sort(unique(df$time)) 
  
  if(is.numeric(var.kint)){
    var.kint=var.kint
  } else{
    var.knots = seq(min(df$time), max(df$time), length= nknot+2 ) 
    var.kint=var.knots[-c(1,length(var.knots))]
    var.kint.con=var.knots[-c(1,length(var.knots))]
  }
  
  ibasis0.con = iSpline(t1, knots = var.kint.con , degree = ndegree) #basis ftn for one param
  ibasis0.all.con = iSpline(t_eval, knots = var.kint.con , degree = ndegree) #basis ftn for one param
  
   ibasis0 = iSpline(t1, knots = var.kint , degree = ndegree) #basis ftn for one param
   ibasis0.all = iSpline(t_eval, knots = var.kint , degree = ndegree) #basis ftn for one param
   
   ibasis = (cbind( rep(1, nrow(ibasis0)), ibasis0))
   ibasis.all = (cbind( rep(1, nrow(ibasis0.all)), ibasis0.all))#[,-6]

  ibasis.con = (cbind( rep(1, nrow(ibasis0.con)), ibasis0.con))
  ibasis.all.con = (cbind( rep(1, nrow(ibasis0.all.con)), ibasis0.all.con))#[,-6]
  
  
  ###### need to fix this for the general form
  #### depeding on p # 06/23/23
  ### previous version - using different knots for constrained and unconstrained beta's
  #ibasis = ibasis.con
  #ibasis.all = ibasis.all.con
  
  Bmat = as.matrix( bdiag(ibasis.con, ibasis, ibasis.con, ibasis) )
  Bmat.all = as.matrix( bdiag(ibasis.all.con, ibasis.all, ibasis.all.con, ibasis.all) ) # to evaluate from t_eval grids
  

  
  #head(df)
  Q=0;C=0
  for (i1 in 1:nsubj){
    df_sub= df[which(df$subject==i1),]
    ni = nrow(df_sub)
    ti = df_sub$time
    Bi = Bmat
    Xmat= as.matrix( data.frame( df_sub[, c('X1','X2','X3','X4')]) )
    yi = as.matrix(df_sub$y)
    wi = diag( 1/ni, ni )
    
    basis_id= which( t1 %in% ti )
    Ui = matrix(0, ncol= ncol(Bmat), nrow= nrow(Xmat))
    for (j1 in 1:ni){
      j2= basis_id[j1]
      Ui[j1,] = t(Xmat[j1,]) %*% Bmat[c(j2, nrow(ibasis)+j2, 2*nrow(ibasis)+j2,3*nrow(ibasis)+j2),]  # dep on p
    }
    Qi = t(Ui)%*%wi%*%Ui
    ci =  t(yi)%*%wi%*%Ui
    
    Q= Q+Qi
    C= C+ci
  }
  #theta_u= solve(Q)%*% t(C) #unconstrained est 
  
  ###########
  # p=1
  if (is.numeric(range1)){ 
    inc_id1 = which( tau> range1[1] & tau<range1[2])
    #knot_rm1= which( ibasis.con0[min(inc_id),] ==1   )
    #knot_rm2= which( ibasis.con0[max(inc_id), 1:ncol(ibasis.con0)] ==0 )
    
    knot_rm11= which( round(ibasis.con[min(inc_id1),],3) ==1   )
    knot_rm12= which( round(ibasis.con[max(inc_id1), 1:ncol(ibasis.con)],3) ==0 )
    knot_rm1 = c(knot_rm11, knot_rm12)
  } else{
    knot_rm1 =c(1:ncol(ibasis.con))
  }
  
  #p=2
  if (is.numeric(range2)){ 
    inc_id2 = which( tau> range2[1] & tau<range2[2])
    knot_rm21= which( round(ibasis[min(inc_id2),],3) ==1   )
    knot_rm22= which( round(ibasis[max(inc_id2), 1:ncol(ibasis)],3) ==0 )
    knot_rm2 =  c(knot_rm21, knot_rm22)
  } else{
    #knot_rm12 =vector()
    knot_rm2 = c(1:ncol(ibasis)) 
  }  
  
  # p=3
  if (is.numeric(range3)){ 
    inc_id3 = which( tau> range3[1] & tau<range3[2])
    knot_rm31= which( round(ibasis.con[min(inc_id3),],3) ==1   )
    knot_rm32= which( round(ibasis.con[max(inc_id3), 1:ncol(ibasis.con)],3) ==0 )
    knot_rm3 =  c(knot_rm31, knot_rm32)
  } else{
    knot_rm3 =  c(1:ncol(ibasis.con)) 
  }  
  
  # p=4
  if (is.numeric(range4)){ 
    inc_id4 = which( tau> range4[1] & tau<range4[2])
    knot_rm41= which( round(ibasis[min(inc_id4),],3) ==1   )
    knot_rm42= which( round(ibasis[max(inc_id4), 1:ncol(ibasis)],3) ==0 )
    knot_rm4 =  c(knot_rm41, knot_rm42)
  } else{
    knot_rm4 =  c(1:ncol(ibasis)) 
  }  
  
  #knot_rm= c(knot_rm1, knot_rm2, knot_rm3, knot_rm4)
  knot_c1 = c(1:ncol(ibasis.con))[!(c(1:ncol(ibasis.con)) %in% knot_rm1)]
  knot_c2 = c(1:ncol(ibasis))[!(c(1:ncol(ibasis)) %in% knot_rm2)]
  knot_c3 = c(1:ncol(ibasis.con))[!(c(1:ncol(ibasis.con)) %in% knot_rm3)]
  knot_c4 = c(1:ncol(ibasis))[!(c(1:ncol(ibasis)) %in% knot_rm4)]
  
  
  # amat0= diag( 1, ncol(Bmat) )
  # if (length(knot_rm)>1){
  #   amat= amat0[-c(knot_rm),]
  # } else{
  #   amat= amat0 
  # }
  amat = matrix(0, nrow=length(knot_c1)+length(knot_c2)+length(knot_c3)+length(knot_c4),ncol=ncol(Bmat))
  # p=1
  if(!is.na(const[1])){
    if(const[1] == "increasing"){
      for(i1 in 1:length(knot_c1)){
        amat[i1, knot_c1[i1]]=1
      }
    }else if(const[1] == "decreasing"){
      for(i1 in 1:length(knot_c1)){
        amat[i1, knot_c1[i1]]=-1
      }
    } 
  }
  
  # p=2
  if(!is.na(const[2])){
    if(const[2] == "increasing"){
      for(i2 in 1:length(knot_c2)){
        amat[length(knot_c1)+i2, ncol(ibasis.con) + knot_c2[i2]]=1
      }
    }else if(const[2] == "decreasing"){
      for(i2 in 1:length(knot_c2)){
        amat[length(knot_c1)+i2,  ncol(ibasis.con) + knot_c2[i2]]=-1
      }
    }
  }
  
  # p=3
  if(!is.na(const[3])){
    if(const[3] == "increasing"){
      for(i3 in 1:length(knot_c3)){
        amat[length(knot_c1)+length(knot_c2)+i3, ncol(ibasis.con) +ncol(ibasis) + knot_c3[i3]]=1
      }
    }else if(const[3] == "decreasing"){
      for(i3 in 1:length(knot_c3)){
        amat[length(knot_c1)+length(knot_c2)+i3,  ncol(ibasis.con) +ncol(ibasis) + knot_c3[i3]]=-1
      }
    }
  }
  
  # p=4
  if(!is.na(const[4])){
    if(const[4] == "increasing"){
      for(i4 in 1:length(knot_c4)){
        amat[length(knot_c1)+length(knot_c2)+length(knot_c3)+i4, 2*ncol(ibasis.con) +ncol(ibasis)+ knot_c4[i4]]=1
      }
    }else if(const[4] == "decreasing"){
      for(i4 in 1:length(knot_c4)){
        amat[length(knot_c1)+length(knot_c2)+length(knot_c3)+i4,  2*ncol(ibasis.con) +ncol(ibasis)+ knot_c4[i4]]=-1
      }
    }
  }
  ###########
  if(qr(Q)$rank != ncol(Q)){
    Q <- as.matrix(nearPD(Q)$mat)
  }
  b= rep(10^-10, nrow(amat)) 
  qp.res= qprog(Q, C, amat, b )
  
  # if(const=="decreasing"){
  #   qp.res= qprog(Q, C, -amat, b ) # for decreasing
  # } else if (const=="increasing"){
  #   qp.res= qprog(Q, C, amat, b ) # for increasing
  # }
  # 
  
  est= Bmat%*%qp.res$thetahat #const est
  
  const_b1hat= (Bmat%*%qp.res$thetahat)[1:length(tau)]
  const_b2hat= (Bmat%*%qp.res$thetahat)[(length(tau)+1):(2*length(tau))] 
  const_b3hat= (Bmat%*%qp.res$thetahat)[(2*length(tau)+1):(3*length(tau))] 
  const_b4hat= (Bmat%*%qp.res$thetahat)[(3*length(tau)+1):(4*length(tau))] 
  const_res = data.frame( t=tau, b1hat=const_b1hat, b2hat=const_b2hat, b3hat=const_b3hat, b4hat=const_b4hat )
  
  const_b1hat.all= (Bmat.all%*%qp.res$thetahat)[1:length(t_eval)]
  const_b2hat.all= (Bmat.all%*%qp.res$thetahat)[(length(t_eval)+1):(2*length(t_eval))] 
  const_b3hat.all= (Bmat.all%*%qp.res$thetahat)[(2*length(t_eval)+1):(3*length(t_eval))] 
  const_b4hat.all= (Bmat.all%*%qp.res$thetahat)[(3*length(t_eval)+1):(4*length(t_eval))] 
  const_res.teval = data.frame( t=t_eval, b1hat=const_b1hat.all, b2hat=const_b2hat.all, b3hat=const_b3hat.all, b4hat=const_b4hat.all )
  
  #uncost est
  
  res= list( const= const_res, const_teval=const_res.teval)
  
  return(res)
}


# 07/13/23 modify for CD4 data with different number of knots # not modified yet ##
vcm_irr_V2 = function(df, p=4, nknot=c(1,1,3,5), t_eval=t_eval, ndegree=2, var.kint=NA,
                      const="decreasing", range1=NA, range2=NA, range3=NA, range4=NA, tol=10^-5){
  
  n= length(unique(df$subject))
  tau=sort(unique(df$time)) 
  
  range.list <- list(range1, range2, range3, range4)
  # 1st covariate
  ibasis<-list(); ibasis_eval<-list()
  Bmat<-matrix(nrow=1, ncol=1)
  Bmat_eval<-matrix(nrow=1, ncol=1)
  for( ic in 1:p){
    var.knots.tmp <- seq(min(df$time), max(df$time), length= nknot[ic]+2 ) 
    var.kint.tmp<- var.knots.tmp[-c(1,length(var.knots.tmp))]
    ibasis[[ic]]<- cbind(1,iSpline(tau, knots = var.kint.tmp , spline.degree= ndegree))
    ibasis_eval[[ic]] <- cbind(1,iSpline(t_eval, knots = var.kint.tmp , spline.degree= ndegree))
    Bmat = as.matrix(bdiag(Bmat, ibasis[[ic]]))
    Bmat_eval = as.matrix(bdiag(Bmat_eval, ibasis_eval[[ic]]))
  }
  
  Bmat<- Bmat[-1,-1]
  Bmat_eval<- Bmat_eval[-1,-1]
  
  subid_all= unique(df$subject)
  Q=C=0
  for (i1 in 1:n){
    df_sub= df[which(df$subject==subid_all[i1] ),]
    ni = nrow(df_sub)
    ti = df_sub$time
    Bi = Bmat  
    #df_sub$intercept= rep(1, nrow(df_sub))
    Xmat= as.matrix( data.frame( df_sub[, c((ncol(df_sub)-p+1): ncol(df_sub))]) ) 
    yi = as.matrix(df_sub$y)
    
    wi = diag( 1/ni, ni )
    
    
    #basis_id= which( tau %in% ti )
    basis_id= match( ti , tau )
    Ui = matrix(0, ncol= ncol(Bmat), nrow= nrow(Xmat))
    for (j1 in 1:ni){
      j2= basis_id[j1]
      Ui[j1,] = t(Xmat[j1,]) %*% Bmat[c(j2, nrow(ibasis[[1]])+j2, 2*nrow(ibasis[[2]])+j2, 3*nrow(ibasis[[3]])+j2 ), ]
    }
    Qi = t(Ui)%*%wi%*%Ui
    ci =  t(yi)%*%wi%*%Ui
    
    Q= Q+Qi
    C= C+ci
    #print( c( i1,Qi[1,1],Q[1,1])  )
  }
  
  theta_u= solve(Q+diag(tol, ncol(Q)) )%*% t(C) #unconstrained est 
  
  #### constrained  
  knot.rm.vec <- c()
  basis.col.vec<- c(0,unlist(lapply(ibasis, ncol)))
  for (ip in 1:p){
    
    basis.tmp <- ibasis[[ip]]
    if (is.numeric(range.list[[ip]])){ 
      inc_id = which( tau> range.list[[ip]][1] & tau< range.list[[ip]][2])
      
      knot_rm1= which( round(basis.tmp[min(inc_id),],1) ==1   )
      knot_rm2= which( round(basis.tmp[max(inc_id), 1:ncol(basis.tmp)],1) ==0 )
      knot_rm12 = c(knot_rm1, knot_rm2)
    } else{
      knot_rm12 =c(1:ncol(basis.tmp))
    }
    
    knot.rm.vec <- c(knot.rm.vec, cumsum(basis.col.vec)[ip] + knot_rm12)
  }
  
  amat0= diag( 1, ncol(Bmat) )
  if (length(knot.rm.vec)>1){
    amat= amat0[-c(knot.rm.vec),]
  } else{
    amat= amat0 
  }
  ###########
  
  b= rep(10^-10, nrow(amat)) 
  
  if(const=="decreasing"){
    qp.res= qprog(Q, C, -amat, b ) # for decreasing
  } else if (const=="increasing"){
    qp.res= qprog(Q, C, amat, b ) # for increasing
  }
  
  est= Bmat%*%qp.res$thetahat #const est
  
  
  #constr fit 
  const_b1hat= (Bmat%*%qp.res$thetahat)[1:length(tau)]
  const_b2hat= (Bmat%*%qp.res$thetahat)[(length(tau)+1):(2*length(tau))] 
  const_b3hat= (Bmat%*%qp.res$thetahat)[(2*length(tau)+1):(3*length(tau))] 
  const_b4hat= (Bmat%*%qp.res$thetahat)[(3*length(tau)+1):length(est)] 
  const_res = data.frame( t=tau, b1hat=const_b1hat, b2hat=const_b2hat, b3hat=const_b3hat, b4hat=const_b4hat )
  
  
  const_b1hat_eval= (Bmat_eval%*%qp.res$thetahat)[1:length(t_eval)]
  const_b2hat_eval= (Bmat_eval%*%qp.res$thetahat)[(length(t_eval)+1):(2*length(t_eval))] 
  const_b3hat_eval= (Bmat_eval%*%qp.res$thetahat)[(2*length(t_eval)+1):(3*length(t_eval))] 
  const_b4hat_eval= (Bmat_eval%*%qp.res$thetahat)[(3*length(t_eval)+1):(4*length(t_eval))] 
  const_res_eval = data.frame( t=t_eval, b1hat=const_b1hat_eval, b2hat=const_b2hat_eval, b3hat=const_b3hat_eval, b4hat=const_b4hat_eval )
  
  
  #uncost est
  uncost_fit = Bmat %*% theta_u
  unconst_b1hat=uncost_fit[1:length(tau)]
  unconst_b2hat=uncost_fit[(length(tau)+1):(2*length(tau))]
  unconst_b3hat=uncost_fit[(2*length(tau)+1):(3*length(tau))] 
  unconst_b4hat=uncost_fit[(3*length(tau)+1):(4*length(tau))] 
  unconst_res= data.frame( t=tau, b1hat=unconst_b1hat, b2hat=unconst_b2hat, b3hat=unconst_b3hat, b4hat=unconst_b4hat)
  
  uncost_fit_eval = Bmat_eval %*% theta_u
  unconst_b1hat_eval=uncost_fit_eval[1:length(t_eval)]
  unconst_b2hat_eval=uncost_fit_eval[(length(t_eval)+1):(2*length(t_eval))] 
  unconst_b3hat_eval=uncost_fit_eval[(2*length(t_eval)+1):(3*length(t_eval))] 
  unconst_b4hat_eval=uncost_fit_eval[(3*length(t_eval)+1):(4*length(t_eval))]
  unconst_res_eval= data.frame( t=t_eval, b1hat=unconst_b1hat_eval, b2hat=unconst_b2hat_eval,b3hat=unconst_b3hat_eval, b4hat=unconst_b4hat_eval)
  
  
  res= list( const= const_res, const_eval= const_res_eval,
             unconst= unconst_res, unconst_eval= unconst_res_eval)
  
  return(res)
}


vcm = function(df, t_eval=t_eval_re, nknot=5, ndegree=3, range1=c(0,1), range2=NA ){ # add ndegree as the option
  
  ndegree= 3 
  n= length(unique(df$subject))
  nsubj=n
  
  #t1=tau
  t1= sort(unique(df$time))  
  #if(length(t1)<1){
  #    t1= sort(unique(df$time))  
  # }
  
  #if(length(tau)<1){
  tau=sort(unique(df$time));
  #  t1=sort(unique(df$time)) 
  #}
  
  
  var.knots = seq(min(df$time), max(df$time), length= nknot+2 )
  var.kint=var.knots[-c(1,length(var.knots))]
  ibasis0 = iSpline(t1, knots = var.kint , spline.degree = ndegree) #basis ftn for one param
  ibasis0.all = iSpline(t_eval, knots = var.kint , spline.degree = ndegree) #basis ftn for one param
  
  ibasis0 = cbind( rep(1, nrow(ibasis0)), ibasis0)
  ibasis=ibasis0
  
  ibasis0.all = cbind( rep(1, nrow(ibasis0.all)), ibasis0.all)
  ibasis.all=ibasis0.all
  
  ###### need to fix this for the general form
  
  Bmat = as.matrix( bdiag(ibasis0, ibasis0) ) 
  Bmat.all = as.matrix( bdiag(ibasis0.all, ibasis0.all) ) # to evaluate from t_eval grids
  #Bmat = ibasis0 
  #matplot(ibasis0, type='l')
  
  #head(df)
  Q=0;C=0
  for (i1 in 1:nsubj){
    df_sub= df[which(df$subject==i1),]
    ni = nrow(df_sub)
    ti = df_sub$time
    Bi = Bmat  
    Xmat= as.matrix( data.frame( df_sub[, c('X1','X2')]) )
    yi = as.matrix(df_sub$y)
    
    wi = diag( 1/ni, ni )
    
    Ui = matrix(0, ncol= ncol(Bmat), nrow= nrow(Xmat))
    for (j1 in 1:ni){
      Ui[j1,] = t(Xmat[j1,]) %*% Bmat[c(j1, ni+j1),]  
    }
    
    Qi = t(Ui)%*%wi%*%Ui
    ci =  t(yi)%*%wi%*%Ui
    
    
    Q= Q+Qi
    C= C+ci
  }
  
  theta_u= solve(Q)%*% t(C) #unconstrained est 
  #######
  # same result using kronecker product # by YPark
  ######
  #X.mat2=c(); #I.mat=c(); 
  #for(i in 1:nsubj){
  #  #I.mat = rbind(I.mat, diag(1, length(t1)))
  #  X.mat2 = rbind(X.mat2, cbind(x_mat[i,1]*ibasis, x_mat[i,2]*ibasis))
  #}  
  
  #theta_u=solve(t(X.mat2)%*%X.mat2)%*%t(X.mat2)%*%y_vec
  #beta1.uncon=as.vector(ibasis%*%theta_u[1:ncol(ibasis)])
  #beta2.uncon=as.vector(ibasis%*%theta_u[(ncol(ibasis)+1):(2*ncol(ibasis))])
  
  ###########
  #for const. matrix 
  if(length(range1)>1 & diff(range1)==1){
    knot_rm11==vector()
  }
  if (length(range1)>1 & diff(range1)<1){ 
    inc_id = which( tau> range1[1] & tau<range1[2])
    knot_rm1= which( ibasis0[min(inc_id),] ==1   )
    knot_rm2= which( ibasis0[max(inc_id), 1:ncol(ibasis0)] ==0 )
    knot_rm11 = c(knot_rm1, knot_rm2)
  } else{
    #knot_rm11 =vector()
    knot_rm11 =c(1:ncol(ibasis0))
  }
  
  if (length(range2)>1){ 
    inc_id2 = which( tau> range2[1] & tau<range2[2])
    knot_rm21= which( ibasis0[min(inc_id2),] ==1   )
    knot_rm22= which( ibasis0[max(inc_id2), 1:ncol(ibasis0)] ==0 ) 
    knot_rm12 = ncol(ibasis) + c(knot_rm21, knot_rm22)
  } else{
    #knot_rm12 =vector()
    knot_rm12 =ncol(ibasis)+c(1:ncol(ibasis0))
  }  
  knot_rm= c(knot_rm11, knot_rm12)
  
  amat0= diag( 1, ncol(Bmat) )
  if (length(knot_rm)>1){
    amat= amat0[-c(knot_rm),]
  } else{
    amat= amat0 
  }
  ###########
  
  b= rep(10^-10, nrow(amat)) 
  
  qp.res= qprog(Q, C, amat, b )
  
  #res$thetahat #const est of coeff 
  est= Bmat%*%qp.res$thetahat #const est
  
  const_b1hat= (Bmat%*%qp.res$thetahat)[1:length(tau)]
  const_b2hat= (Bmat%*%qp.res$thetahat)[(length(tau)+1):(2*length(tau))] 
  const_res = data.frame( t=tau, b1hat=const_b1hat, b2hat=const_b2hat )
  
  const_b1hat.all= (Bmat.all%*%qp.res$thetahat)[1:length(t_eval)]
  const_b2hat.all= (Bmat.all%*%qp.res$thetahat)[(length(t_eval)+1):(2*length(t_eval))] 
  const_res.teval = data.frame( t=t_eval, b1hat=const_b1hat.all, b2hat=const_b2hat.all )
  
  #uncost est
  uncost_fit = Bmat %*% theta_u
  unconst_b1hat=uncost_fit[1:length(tau)]
  unconst_b2hat=uncost_fit[(length(tau)+1):(2*length(tau))] 
  unconst_res= data.frame( t=tau, b1hat=unconst_b1hat, b2hat=unconst_b2hat )
  
  uncost_fit.all = Bmat.all %*% theta_u
  unconst_b1hat.all=uncost_fit.all[1:length(t_eval)]
  unconst_b2hat.all=uncost_fit.all[(length(t_eval)+1):(2*length(t_eval))] 
  unconst_res.teval= data.frame( t=t_eval, b1hat=unconst_b1hat.all, b2hat=unconst_b2hat.all )
  
  uncost_fit.all = Bmat.all %*% theta_u
  unconst_b1hat.all=uncost_fit.all[1:length(t_eval)]
  unconst_b2hat.all=uncost_fit.all[(length(t_eval)+1):(2*length(t_eval))] 
  unconst_res.teval= data.frame( t=t_eval, b1hat=unconst_b1hat.all, b2hat=unconst_b2hat.all )
  
  res= list( const= const_res, const_teval=const_res.teval, 
             unconst= unconst_res, unconst_teval=unconst_res.teval)
  
  return(res)
}


####################################################
# wild bootstrap
boot_fun_sp <- function (B = 1000, df=df, ndegree=ndegree,t_eval=t_eval, var.kint=NA, var.kint.con=NA,
                         nknot=NA, t_list=t_list, y_list=y_list, x_mat=x_mat, 
                         const=c(NA, NA,NA,NA), range1=NA, range2=NA, range3=NA, range4=NA, partial= FALSE) {
  
  L_vec <- unlist(lapply(t_list, length))
  t_vec <- unlist(t_list)
  y_vec <- unlist(y_list)
  # 
  n <- length(t_list)
  N <- sum(L_vec)
  p <- ncol(x_mat)
  M <- length(t_eval)
  
  sample_id= rep( 1:n,  unlist(lapply(t_list, length)) )
  # x_mat = as.matrix(df[,-c(1:3)])
  # 
  # x_mat1 = matrix(ncol=p)#= matrix(0, nrow=length(y_vec), ncol=2)
  # for (i1 in 1:nrow(x_mat)){
  #   aa= x_mat[i1,]
  #   aa1= matrix( rep(aa, length(t_list[[i1]]) ), ncol=p, byrow = TRUE)
  #   x_mat1 =  rbind(x_mat1, aa1)
  # }
  # X = x_mat1[-1,]
  # 
  # df= data.frame(subject= sample_id, time=t_vec, y=t_vec, X1=t_vec, X2=X[,2], X3=X[,3], X4=X[,4] )
  
  res= vcm_irr( df=df,ndegree=ndegree, t_eval=t_eval, const=const, var.kint = var.kint, nknot=nknot,
                range1=domain_con, range2=NA, range3=domain_con, range4=NA)
  
  
  #Bs_vec_hat_null=rbind(res$const_teval$b1hat,res$const_teval$b2hat,res$const_teval$b3hat,res$const_teval$b4hat)
  Bs_vec_hat_alt=rbind(res$unconst_teval$b1hat,res$unconst_teval$b2hat,res$unconst_teval$b3hat,res$unconst_teval$b4hat)
  
  if(partial==FALSE){
    res2= vcm_irr_con( df=df,ndegree=ndegree, t_eval=t_eval, const=const, var.kint =var.kint, var.kint.con=var.kint.con, nknot=nknot,range1=domain_con, range2=NA, range3=domain_con, range4=NA)
  } else if(partial==TRUE){
    res2= vcm_irr_con( df=df,ndegree=ndegree, t_eval=t_eval, const=const, var.kint =var.kint.con, var.kint.con=var.kint.con, nknot=nknot,range1=domain_con, range2=NA, range3=domain_con, range4=NA)
  }
  #res2= vcm_irr( df=df,ndegree=ndegree, t_eval=t_eval, const=const, var.kint =var.kint, nknot=nknot,range1=domain_con, range2=domain_con, range3=NA, range4=NA)
  
  Bs_vec_hat_null=rbind(res2$const_teval$b1hat,res2$const_teval$b2hat,res2$const_teval$b3hat,res2$const_teval$b4hat)
  #Bs_vec_hat_alt=rbind(res$unconst_teval$b1hat,res$unconst_teval$b2hat,res$unconst_teval$b3hat,res$unconst_teval$b4hat)
  
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
  boot_T0 <- boot_T1 <-boot_T2 <- boot_T1_all <-boot_T2_all <- c()
  for (b in 1:B) {
    
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
    
    df.boot= data.frame(cbind(subject= sample_id, time=t_vec, y=y_vec.boot, df[ ,-c(1:3)])) #, X1= X[,1], X2=X[,2], X3=X[,3], X4=X[,4]  )
    res.boot= vcm_irr( df=df.boot,t_eval=t_eval, var.kint=var.kint, nknot = nknot,
                       const=const, range1=domain_con, range2=NA, range3=domain_con, range4=NA)
    Bs_vec_hat_alt_b=rbind(res.boot$unconst_teval$b1hat,res.boot$unconst_teval$b2hat,res.boot$unconst_teval$b3hat,res.boot$unconst_teval$b4hat)
    
    if(partial==FALSE){
      res.boot2= vcm_irr_con( df=df.boot,t_eval=t_eval, var.kint=var.kint, var.kint.con=var.kint.con, nknot = nknot,
                              const=const, range1=domain_con, range2=NA, range3=domain_con, range4=NA)
    } else if (partial==TRUE){
      res.boot2= vcm_irr_con( df=df.boot,t_eval=t_eval, var.kint=var.kint.con, var.kint.con=var.kint.con, nknot = nknot,
                              const=const, range1=domain_con, range2=NA, range3=domain_con, range4=NA)
    }
    Bs_vec_hat_null_b=rbind(res.boot2$const_teval$b1hat,res.boot2$const_teval$b2hat,res.boot2$const_teval$b3hat,res.boot2$const_teval$b4hat)
    
    #B_vec_hat_b <- B_vec_hat_fun(t_eval, t_list, y_list_b, x_mat, h, comp_con, domain_con)
    #B_vec_hat_null_b <- B_vec_hat_b$B_vec_hat_null
    #B_vec_hat_alt_b <- B_vec_hat_b$B_vec_hat_alt
    
    
    T0_b <- 0
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
    
    #F_b <- ifelse(sse_null_b - sse_alt_b > 0, sse_null_b - sse_alt_b, 0) / sse_alt_b
    
    gamma_hat <- t(x_mat) %*% x_mat / n
    ind_eval <- (t_eval >= domain_con[1] & t_eval <= domain_con[2])
    
    T1_b <- sum(diag(t(Bs_vec_hat_null_b[1:p,ind_eval] - Bs_vec_hat_alt_b[1:p,ind_eval]) %*%
                       gamma_hat %*%
                       (Bs_vec_hat_null_b[1:p,ind_eval] - Bs_vec_hat_alt_b[1:p,ind_eval]))) * diff(t_eval)[1]
    
    T2_b <- sum(diag(t(Bs_vec_hat_null_b[1:p,ind_eval] - Bs_vec_hat_alt_b[1:p,ind_eval]) %*%
                       (Bs_vec_hat_null_b[1:p,ind_eval] - Bs_vec_hat_alt_b[1:p,ind_eval]))) * diff(t_eval)[1]
    
    T1_all_b <- sum(diag(t(Bs_vec_hat_null_b[1:p,] - Bs_vec_hat_alt_b[1:p,]) %*%
                           gamma_hat %*%
                           (Bs_vec_hat_null_b[1:p,] - Bs_vec_hat_alt_b[1:p,]))) * diff(t_eval)[1]
    
    T2_all_b <- sum(diag(t(Bs_vec_hat_null_b[1:p,] - Bs_vec_hat_alt_b[1:p,]) %*%
                           (Bs_vec_hat_null_b[1:p,] - Bs_vec_hat_alt_b[1:p,]))) * diff(t_eval)[1]
    
    
    boot_T0[b] <- T0_b
    boot_T1[b] <- T1_b
    boot_T2[b] <- T2_b
    boot_T1_all[b] <- T1_all_b
    boot_T2_all[b] <- T2_all_b
    #boot_F[b] <- F_b
    
  }
  
  return (list(boot_T0 = boot_T0,
               boot_T1 = boot_T1,
               boot_T2 = boot_T2,
               boot_T1_all = boot_T1_all,
               boot_T2_all = boot_T2_all
  ))
}



boot_fun_sp_V2 <- function (B = 500, nknot, var.kint, t_eval, t_list, y_list, x_mat, 
                            range1=range1, range2=range2, range3=range3, range4=range4) {
  
  L_vec <- unlist(lapply(t_list, length))
  t_vec <- unlist(t_list)
  y_vec <- unlist(y_list)
  
  n <- length(t_list)
  N <- sum(L_vec)
  p <- ncol(x_mat)
  M <- length(t_eval)
  sample_id= rep( 1:n,lengths(t_list) )
  
  x_mat1 = matrix(0, ncol=p)#= matrix(0, nrow=length(y_vec), ncol=2)
  for (i1 in 1:nrow(x_mat)){
    aa= x_mat[i1,]
    aa1= matrix( rep(aa, length(t_list[[i1]]) ), ncol=p, byrow = TRUE)
    x_mat1 =  rbind(x_mat1, aa1)
  }
  X = x_mat1[-1,]
  df= data.frame(subject= sample_id, time=t_vec, y=y_vec, 
                 X1= as.numeric(X[,1]), X2=as.numeric(X[,2]), X3=as.numeric(X[,3]), X4=as.numeric(X[,4])  )
  
  # fitted values under the alternative
  res.tmp= vcm_irr_V2( df=df, p=4, t_eval=t_eval, nknot=nknot,ndegree=ndegree, const=shape.const,var.kint = NA,
                       range1=range1, range2=range2, range3=range3, range4=range4, tol=tol)
  
  
  Bs_vec_hat_null=rbind(res.tmp$const_eval$b1hat,res.tmp$const_eval$b2hat,res.tmp$const_eval$b3hat,res.tmp$const_eval$b4hat)
  Bs_vec_hat_alt=rbind(res.tmp$unconst_eval$b1hat,res.tmp$unconst_eval$b2hat,res.tmp$unconst_eval$b3hat,res.tmp$unconst_eval$b4hat)
  
  y_fit_null_list <- y_fit_alt_list <- resid_alt_list <- list()
  
  for (i in 1:n) {
    
    # cat(i, '\n')
    
    tB_vec_fit_null <- sapply(1:p,
                              function (j) spline(t_eval, Bs_vec_hat_null[j,], xout = t_list[[i]])$y
    )
    
    tB_vec_fit_alt <- sapply(1:p,
                             function (j) spline(t_eval, Bs_vec_hat_alt[j,], xout = t_list[[i]])$y
    )
    
    y_fit_null_list[[i]] <- c(matrix(tB_vec_fit_null, ncol=p) %*% matrix(as.numeric(x_mat[i,]), ncol=1))
    y_fit_alt_list[[i]] <- c(matrix(tB_vec_fit_alt, ncol=p) %*% matrix(as.numeric(x_mat[i,]), ncol=1))
    
    resid_alt_list[[i]] <- y_fit_alt_list[[i]] - y_list[[i]]
    
  }
  
  
  # Mammen's wild bootstrap with moment correction
  boot_pr <- (5 + sqrt(5))/10
  
  # wild bootstrap
  boot_T0 <- boot_T1 <- c()
  for (b in 1:B) {
    
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
    
    df.boot= data.frame(subject= sample_id, time=t_vec, y=y_vec.boot, 
                        X1= as.numeric(X[,1]), X2=as.numeric(X[,2]), X3=as.numeric(X[,3]), X4=as.numeric(X[,4]))
    res.boot= vcm_irr_V2( df=df.boot, p=4, t_eval=t_eval, nknot=nknot, ndegree=ndegree, const=shape.const, var.kint = NA,
                          range1=range1, range2=range2, range3=range3, range4=range4, tol=tol)
    
    
    Bs_vec_hat_null_b=rbind(res.boot$const_eval$b1hat,res.boot$const_eval$b2hat,res.boot$const_eval$b3hat,res.boot$const_eval$b4hat)
    Bs_vec_hat_alt_b=rbind(res.boot$unconst_eval$b1hat,res.boot$unconst_eval$b2hat,res.boot$unconst_eval$b3hat,res.boot$unconst_eval$b4hat)
    
    
    T0_b <- 0
    for (i in 1:n) {
      
      # cat(i, '\n')
      
      tB_vec_fit_null_b <- sapply(1:p,
                                  function (j) spline(t_eval, Bs_vec_hat_null_b[j,], xout = t_list[[i]])$y
      )
      
      tB_vec_fit_alt_b <- sapply(1:p,
                                 function (j) spline(t_eval, Bs_vec_hat_alt_b[j,], xout = t_list[[i]])$y
      )
      
      T0_b <- T0_b + mean(((matrix(tB_vec_fit_null_b - tB_vec_fit_alt_b, ncol=p)) %*% matrix(as.numeric(x_mat[i,]),ncol=1))^2)
    }
    
    T1_b <- sqrt(n) * sum(apply((Bs_vec_hat_null_b[1:2,] - Bs_vec_hat_alt_b[1:2,])^2, 1, sum)) * diff(t_eval)[1]
    
    boot_T0[b] <- T0_b
    boot_T1[b] <- T1_b
    
    if(b %% 100 ==0) print(b)
  }
  
  return (list(boot_T0 = boot_T0,
               boot_T1 = boot_T1))
}
