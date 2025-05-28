mfpca.face <- function(data, newdata=NULL, pve=0.95, knots=10, lambda = 0, center=T){
  ## this is multivariate face for multivariate functional data

  require(refund)
  FPCA <- list()
  J <- length(data)
  for(j in 1:J){
    FPCA[[j]] <- fpca.face(data[[j]], knots=knots, argvals =(1:ncol(data[[j]]))/ncol(data[[j]]),
                           pve= pve,
                           lambda=lambda,
                           center=center)
  }

  ## first recover covariance functions
  C_list <- vector("list",J)


  for(j1 in 1:J){
    for(j2 in 1:J){
      if(j1==j2) temp <- diag(FPCA[[j1]]$evalues)
      if(j1!=j2) temp <- cov(FPCA[[j1]]$scores, FPCA[[j2]]$scores)
      C_list[[j1]][[j2]] <- FPCA[[j1]]$efunctions%*%temp%*%t(FPCA[[j2]]$efunctions)
    }
  }


  ##  MFPCA
  c <- unlist(lapply(data,function(x){dim(x)[2]}))
  c_all <- sum(c)
  Cs <- matrix(NA,c_all, c_all)
  for(j1 in 1:J){
    for(j2 in 1:J){
      if(j1==1) sel1 <- 1:c[1]
      if(j1>1) sel1 <- sum(c[1:(j1-1)]) + (1:c[j1])
      if(j2==1) sel2 <- 1:c[1]
      if(j2>1) sel2 <- sum(c[1:(j2-1)]) + (1:c[j2])

      Cs[sel1,sel2] <- C_list[[j1]][[j2]]
    }
  }
  Eigen <- eigen(Cs)
  d <- which.max(cumsum(Eigen$values)/sum(Eigen$values)>=pve)
  U_joint <- Eigen$vectors[,1:d]
  lambda <- Eigen$values[1:d]


  ## for prediction of newdata
  if(is.null(newdata)) newdata = data
  n_new <- dim(newdata[[1]])[1]
  m <- lapply(newdata,function(x){dim(x)[2]})

  mu_newdata <- list()
  c_newdata <- list()
  sc_newdata <- list()
  for(j in 1:J){
    mu_newdata[[j]] <- matrix(FPCA[[j]]$mu,nrow=n_new,ncol=m[[j]],byrow=TRUE)
    c_newdata[[j]] <- newdata[[j]] - mu_newdata[[j]]
    sc_newdata[[j]] <- c_newdata[[j]]%*%FPCA[[j]]$efunctions%*%t(FPCA[[j]]$efunctions)
  }


  scores <- do.call(cbind,sc_newdata)%*%U_joint
  sc_newdata_pred <- scores%*%t(U_joint)



  #prediction by  MFPCA
  newdata_pred <- list()
  for(j1 in 1:J){
    if(j1==1) sel <- 1:c[1]
    if(j1>1) sel <- sum(c[1:(j1-1)]) + (1:c[j1])
    newdata_pred[[j1]] <- sc_newdata_pred[,sel] + mu_newdata[[j1]]
  }

  ## prediction via BLUP

  sig2_list <- lapply(FPCA,function(x){mean((t(x$Yhat-x$Y)-x$mu)^2)})
  for(j in 1:J){

    if(j==1) sel <- 1:c[1]
    if(j>1) sel <- sum(c[1:(j-1)]) + (1:c[j])

    Cs[sel,sel] <- Cs[sel,sel] + diag(length(sel))*sig2_list[[j]]
  }

  Lambda <- diag(lambda)
  scores <- Lambda%*%t(U_joint)%*%solve(Cs)%*%t(do.call(cbind,c_newdata))
  scores <- t(scores)
  newdata_predu <- scores%*%t(U_joint)

  newdata_pred_blup <- list()
  for(j in 1:J) {
    if(j==1) sel <- 1:c[1]
    if(j>1) sel <- sum(c[1:(j-1)]) + (1:c[j])
    newdata_pred_blup[[j]] <- mu_newdata[[j]] + newdata_predu[,sel]

  }

  res <- list(
    data = data,
    FPCA = FPCA,
    J = J,
    C_list = C_list,
    efunctions = U_joint,
    evalues = Eigen$values[1:d],
    npc = d,
    scores = scores,
    newdata_pred = newdata_pred,
    newdata_pred_blup = newdata_pred_blup,
    c = c,
    c_all = c_all,
    scores = scores
  )
  class(res) <- "mfpca.face"
  return(res)
}
