FisherInfoInv <- function(fit,data){
  ## this function calculates Fisher information matrix used for asLRT
  subj <- unique(data$ID)
  nsubj <- length(subj)
  npc <- (ncol(data) - 3)/2
  var.est <- as.data.frame(VarCorr(fit))[,4]
  tau0 <- var.est[1]
  tau1 <- var.est[2]
  sigma2 <- var.est[3]
  FI <- matrix(0, npc + 3, npc + 3)

  for(i in 1:nsubj){
    sel <- which(data$ID==subj[i])
    ni <- length(sel)
    datai <- data[sel,]
    Sigmai <- sigma2*diag(ni) + tau0*datai$Z%x%t(datai$Z) + tau1*rep(1,ni)%*%t(rep(1,ni))
    Sigmai_inv <- solve(Sigmai)
    for(j in 1:(npc+1)){
      for(k in 1:(npc+1)){
        FI[j,k] <- FI[j,k] + 1/2*(sum(datai[, j+2]*(Sigmai_inv%*%datai[,k+2])))^2
      }
    }
    for(j in 1:(npc+1)){
      FI[npc+2, j] <- FI[npc+2, j] + 1/2*sum((Sigmai_inv%*%datai[,j+2])^2)
      FI[j, npc+2] <- FI[npc+2, j]
    }
    FI[npc+2, npc+2] <- FI[npc+2, npc+2] + 1/2*sum(Sigmai_inv^2)

    ##for subject-specific random intercept
    for(j in 1:(npc+1)){
      FI[j, npc + 3] <- FI[j, npc + 3] + 1/2*(sum(datai[, j+2]*(Sigmai_inv%*%rep(1,ni))))^2
      FI[npc + 3, j] <- FI[j, npc + 3]
    }
    FI[npc + 2, npc + 3] <- FI[npc + 2, npc + 3] + 1/2*sum((Sigmai_inv%*%rep(1,ni))^2)
    FI[npc + 3, npc + 2] <- FI[npc + 2, npc + 3]
    FI[npc + 3, npc + 3] <- FI[npc + 3, npc + 3] + 1/2*(sum(rep(1,ni)*(Sigmai_inv%*%rep(1,ni))))^2

  }
  return(solve(FI)[1:npc,1:npc])
}

Chi_bar_square <- function(lrt, R,alpha){
  #this function evaluates the chi-bar asymptotic dist. for asLRT
  npc <- nrow(R)
  RC <- diag(sqrt(1/diag(R)))%*%R%*%diag(sqrt(1/diag(R)))
  if(npc == 3){

    RPC <- rep(0,3)
    RPC[1] <- (RC[2,3] - RC[3,1]*RC[2,1])/sqrt(1-RC[3,1]^2)/sqrt(1-RC[3,1]^2)
    RPC[2] <- (RC[1,3] - RC[1,2]*RC[3,2])/sqrt(1-RC[1,2]^2)/sqrt(1-RC[3,2]^2)
    RPC[3] <- (RC[1,2] - RC[1,3]*RC[2,3])/sqrt(1-RC[1,3]^2)/sqrt(1-RC[2,3]^2)
    w <- rep(0,4)
    w[4] <- 1/(4*pi)*(2*pi - acos(RC[1,2]) - acos(RC[1,3]) - acos(RC[2,3]))
    w[3] <- 1/(4*pi)*(3*pi - acos(RPC[3]) - acos(RPC[2]) - acos(RPC[1]))
    w[2] <- 1/2 - w[4]
    w[1] <- 1/2 - w[3]
  }

  if(npc == 2){
    w <- rep(0, 3)
    w[1] <- 1/(2*pi)*acos(RC[1,2])
    w[2] <- 1/2
    w[3] <- 1/2 - w[1]
  }

  if(npc ==1){
    w <- c(1/2,1/2)
  }

  obj <- function(x){
    prob <- 0
    for(i in 0:npc){
      prob <- prob +  w[i+1]*pchisq(x,df=i)
    }
    return((prob-(1-alpha))^2)
  }
  cutoff <- nlminb(w[1],obj)$par

  pvalue <- 0
  for(i in 1:npc){
    pvalue <- pvalue + w[i+1]*(1-pchisq(lrt,df=i))
  }

  return(list("observed" = lrt, "cutoff" = cutoff,
              "pvalue" = pvalue ))
}

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


fmm_sofr <- function(data,
                     multi=1,
                     control = list(
                              pve = 0.95,
                              knots = 7,
                              p = 5,
                              m = 3,
                              center = TRUE,
                              lambda = 0)){

  # input contains data and control
  # both data and control are lists
  # multi denotes the number of multivariate functional predictor; an integer needs to be specified for multi; by default, multi = 1

  Y <- data$Y   # response vector Yij
  Z <- data$Z   # numeric variable vector Zij
  Xt <- as.matrix(data$Xt) # functional covariate Xij(t) matrix
  ID <- data$ID # subject index

  subs <- unique(ID)
  nSubj <- length(subs)
  nRep <- sapply(subs, function(x) length(which(ID==x)))

  t <- 1:(dim(Xt)[2]/multi)/(dim(Xt)[2]/multi) # time grid
  D <- length(t) # time grid length

  library(refund)
  library(lme4)
  library(nlme)
  library(arm)
  library(RLRsim)
  library(MASS)

 ###########################################################################

  if(multi>1){
    ### mfpca on Xt, here Xt needs to be reformulated into a list ###
    Xt_list <- list()
    for(i in 1:multi){
      Xt_list[[i]] <- as.matrix(Xt[,((i-1)*D+1):(i*D)])
    }

    results <- mfpca.face(Xt_list,
                       center = control$center,
                       knots = 10, pve = control$pve, lambda = control$lambda)

  }else{
    ### fpca on Xt ###
    results <- fpca.face(Xt,
                       center = control$center, argvals = t,
                       knots = control$knots, pve = control$pve,
                       p = control$p, lambda = control$lambda)
  }

  npc <- results$npc
  score <- results$scores
  ascore <- score[, 1:npc]/sqrt(D)
  efunctions <- results$efunctions*sqrt(D)
  evalues <- results$evalues/D
  bscore <- ascore*Z #interaction term

  ###########################################################################

  ### data for bonferroni testing and FLM/FMM regression ###

  ## note by Luo: the design matrix has to be of this order
  ## for the (asymptotic) LRT

  designMatrix.reg <- data.frame(Y = Y,
                                 ID = as.factor(ID),
                                 ascore = ascore,
                                 Z = Z,
                                 bscore = bscore)

  ###########################################################################

  ### equal variance testing ###
  z.sim.uni = c()
  ID.uni <- c()
  cscore <- c()

  index <- matrix(1:(nSubj*npc), nrow = npc, ncol = nSubj)
  for (i in 1:length(nRep)){
    ID.uni = c(ID.uni, c(index[,i], rep(0, nRep[i] - npc)))
  }
  z.sim.uni = c()
  # svd on random scores A_i for each subject i
  for(k in 1:nSubj){
    if(k==1){
      svd <- svd(ascore[1:nRep[1], ] %*% t(ascore[1:nRep[1], ])) # SVD on A_i
    }else{
      svd <- svd(ascore[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k]), ] %*% t(ascore[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k]), ])) #SVD on A_i
    }
    u.tra <- t(svd$v)
    u <- svd$u
    d <- (svd$d)[1:npc]
    cscore <- c(cscore,  rowSums(u.tra))
    if(k==1){
      Y[1:nRep[k]] <- u.tra %*% Y[1:nRep[k]]
      Z[1:nRep[k]] <- u.tra %*% Z[1:nRep[k]]
      ascore[1:nRep[k], ] <- rbind(u.tra[1:npc, ] %*% ascore[1:nRep[k], ],
                                   matrix(0, nrow = nRep[k] - npc,
                                          ncol = npc))
      bscore[1:nRep[k], ] <- u.tra %*% bscore[1:nRep[k], ]

    }else{
      Y[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k])] <- u.tra %*% Y[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k])]
      Z[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k])] <- u.tra %*% Z[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k])]
      ascore[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k]), ] <- rbind(u.tra[1:npc, ] %*% ascore[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k]), ],
                                                               matrix(0,
                                                                      nrow = nRep[k] - npc,
                                                                      ncol = npc))
      bscore[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k]), ] <- u.tra %*% bscore[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k]), ]
    }
    # z.sim.uni is the coefficient for the random slope to be tested
    z.sim.uni <- c(z.sim.uni, sqrt(d), rep(0, nRep[k] - npc))
  }

  ###########################################################################

  ### data for equal-variance testing ###

  designMatrix <- data.frame(Y = Y,
                             Z = Z,
                             ID = as.factor(ID),
                             ID.uni = as.factor(ID.uni),
                             ascore = ascore,
                             z.sim.uni = z.sim.uni,
                             bscore = bscore,
                             cscore = cscore)

  ###########################################################################

  ## equal-variance test ##
  additive0.sim <- paste(1:npc, collapse = " + ascore.")
  additive.sim.b <- paste(1:npc, collapse = " + bscore.")

  model.sim <- as.formula(paste("Y ~ 1 + Z + ascore.",
                                additive0.sim,
                                " + bscore.",
                                 additive.sim.b,
                                " + (0 + cscore | ID)",
                                " + (0 + Z | ID) + (0 + z.sim.uni | ID.uni)",
                                sep = ""))
  # fullReml is the model under alternative
  fullReml <- lmer(model.sim, data = designMatrix)
  # m.slope only contains the random effect to be tested
  f.slope <- as.formula(paste("Y ~ 1 + Z + ascore.",
                              additive0.sim,
                              " + bscore.",
                              additive.sim.b,
                              " + (0 + z.sim.uni | ID.uni)",
                              sep = ""))
  m.slope <- lmer(f.slope, data = designMatrix)
  # m0 is the model under the null
  f0 <- as.formula(" . ~ . - (0 + z.sim.uni | ID.uni)")
  m0 <- update(fullReml, f0)

  EqualVar.test <- exactRLRT(m.slope, fullReml, m0)
  EqualVar.pvalue <- EqualVar.test$p[1]
  ## end of equal-variance test

  ## bonferroni test ##
  additive.heter <- paste0(" + (0 + ascore.", 1:npc, " | ID)", collapse = "")
  bonf.test <- list()
  for(i in 1:npc){
    ii <- paste("ascore.", i, sep = "")
    # f.slope only contains the random effect to be tested
    f.slope <- as.formula(paste("Y ~ 1 + Z + ascore.",
                                additive0.sim,
                                " + bscore.",
                                additive.sim.b,
                                " + (0 +", ii, " | ID)",
                                sep = ""))
    m.slope <- lmer(f.slope, data = designMatrix.reg)
    # mA is the model under alternative
    mA <- as.formula(paste("Y ~ 1 + Z + ascore.",
                           additive0.sim,
                           " + bscore.",
                           additive.sim.b,
                           " + (0 +", ii, " | ID)", " + (1 | ID)",
                           "+ (0 + Z | ID)",
                           sep = ""))
    fullReml <- lmer(mA, data = designMatrix.reg)
    #m0 is model under the null
    f0 <- as.formula(paste(" . ~ . - (0 + ", ii, "| ID)"))
    m0 <- update(fullReml, f0)
    bonf.test[[i]] <- exactRLRT(m.slope, fullReml, m0)
  }
  multiTest <- sapply(bonf.test, function(x) {
    c(statistic = x$statistic[1],
      "p-value" = x$p[1])})
  # use bonferroni correctiong method to adjust p-value
  bonf.pvalue <- p.adjust(multiTest[2,], "bonferroni")
  ## end of bonferroni test

  ##  LRT (added by Luo Apr 4, 2019)
  ###########################################################################
  additive0.sim <- paste(1:npc, collapse = " + ascore.")
  fixed <- paste("Y ~ 1 + Z + ascore.",
                 additive0.sim,
                 " + bscore.",
                 additive.sim.b,
                 sep="")

  scores <- sapply(sapply(1:npc,function(x){paste("ascore.",x, sep="")}),
                   function(x){paste("(0+",x,"|ID)",sep="")})
  scores <- c(scores)
  random <- paste(scores, collapse=" + ")
  random <- paste("(0 + Z|ID) + ", random, sep="")
  random <- paste0(random, " + (1 | ID)")


  # fullReml is the model under alternative
  model.sim <- as.formula(paste(fixed, "+",random,sep=""))
  fit.full <- lmer(model.sim, data = designMatrix.reg, REML=FALSE)

  model.sim.null <- as.formula(paste(fixed,"+ (0+Z|ID)"," + (1 | ID)", sep=""))
  fit.null <- lmer(model.sim.null, data = designMatrix.reg, REML=FALSE)

  lrt <- -2*(as.numeric(logLik(fit.null))-as.numeric(logLik(fit.full)))


  R <- FisherInfoInv(fit.null,designMatrix.reg)
  LRT.test <- Chi_bar_square(lrt,R,0.05)
  LRT.pvalue <- LRT.test$pvalue
  ##end of LRT

  ###########################################################################

  ## estimation ##

    # FLM: estimation without subject-specific random effect
    noRandom.simpart <- paste(1:npc, collapse = " + ascore.")
    additive.sim.b <- paste(1:npc, collapse = " + bscore.")
    noRandom.sim <- as.formula(paste("Y ~ 1 + Z + ascore.",
                                     noRandom.simpart,
                                     " + bscore.",
                                     additive.sim.b,
                                     " + (0 + Z | ID) + (1 | ID)",
                                     sep = ""))
    FLM <- lmer(noRandom.sim, data = designMatrix.reg)

    # FMM: estimation with subject-specific random effect
    Random.simpart <- paste0(" + (0 + ascore.", 1:npc,"|ID) ", collapse="")
    Random.sim <- as.formula(paste("Y ~ 1 + Z + ascore.",
                                   noRandom.simpart,
                                   " + bscore.",
                                   additive.sim.b,
                                   " + (0 + Z | ID) + (1 | ID)",
                                   Random.simpart,
                                   sep = ""))
    FMM <- lmer(Random.sim, data = designMatrix.reg)
  ###########################################################################

  ## get fixed effect beta(t) ##

  # FLM
  fixeff1 <- fixef(FLM)
  beta1_a <- efunctions %*% as.vector(fixeff1[grep("^ascore", names(fixeff1))]) # population effect beta(t)
  beta1_b <- efunctions %*% as.vector(fixeff1[grep("^bscore", names(fixeff1))]) # delta(t)
  raneffect1 <- ranef(FLM)$ID
  FLM.sum <- summary(FLM)
  fixed.coeff1 <- FLM.sum$coefficients # fixed coefficient
  fixed.vcov1  <- FLM.sum$vcov # coefficient covariance
  # sigma is standard deviation for error term
  sigma1 <- FLM.sum$sigma
  # se_a is standard deviation for population effect beta(t)
  se1_a = apply(efunctions, 1, function(x) sqrt(x %*% as.matrix(fixed.vcov1[match(paste0("ascore.",1:npc),names(fixeff1)),match(paste0("ascore.",1:npc),names(fixeff1))]) %*% x))
  # se_b is standard deviation for delta(t)
  se1_b = apply(efunctions, 1, function(x) sqrt(x %*% as.matrix(fixed.vcov1[match(paste0("bscore.",1:npc),names(fixeff1)),match(paste0("bscore.",1:npc),names(fixeff1))]) %*% x))

  yhat1 <- predict(FLM)

  # FMM
  fixeff2 <- fixef(FMM)
  beta2_a <- efunctions %*% as.vector(fixeff2[grep("^ascore", names(fixeff2))]) # population effect beta(t)
  beta2_b <- efunctions %*% as.vector(fixeff2[grep("^bscore", names(fixeff2))]) # population effect delta(t)
  raneffect2 <- ranef(FMM)$ID
  betai_2 <- efunctions %*% t(as.matrix(raneffect2[grep("^ascore", names(raneffect2))])) # beta_i(t): subject deviation from population effect beta(t)
  FMM.sum <- summary(FMM)
  fixed.coeff2 <- FMM.sum$coefficients # fixed coefficient
  fixed.vcov2  <- FMM.sum$vcov
  sigma2 <- FMM.sum$sigma # error term se
  # se_a is standard deviation for population effect
  se2_a = apply(efunctions, 1, function(x) sqrt(x %*% as.matrix(fixed.vcov2[match(paste0("ascore.",1:npc),names(fixeff2)),match(paste0("ascore.",1:npc),names(fixeff2))]) %*% x))
  # se_b is standard deviation for delta(t)
  se2_b = apply(efunctions, 1, function(x) sqrt(x %*% as.matrix(fixed.vcov2[match(paste0("bscore.",1:npc),names(fixeff2)),match(paste0("bscore.",1:npc),names(fixeff2))]) %*% x))
  yhat2 <- predict(FMM)

  ###########################################################################

  return(list(
        'fPCA_result' = list(
                             fpca_result = results,
                             npc = npc,
                             ascore = ascore,
                             efunctions = efunctions,
                             evalues = evalues),
        'test_result' = list(
                            'equal-variance' = list(
                                                  EqualVar.test = EqualVar.test,
                                                  EqualVar.pvalue = EqualVar.pvalue),
                            'bonferroni' = list(
                                                bonf.test = bonf.test,
                                                bonf.pvalue = bonf.pvalue),
                            'asLRT' = list(
                                           LRT.test = LRT.test,
                                           LRT.pvalue = LRT.pvalue)),
        'estimation_result' = list(
                                  'FLM' = list(
                                              'fixed' = list(
                                                            'beta(t)' = beta1_a,
                                                            'delta(t)' = beta1_b,
                                                            'coefficient' = fixed.coeff1,
                                                            'vcov' = fixed.vcov1,
                                                            'sigma' = sigma1,
                                                            'se_beta(t)' = se1_a,
                                                            'se_delta(t)' = se1_b
                                                            ),
                                              'random' = list('Zi' = raneffect1),
                                              'yhat' = yhat1,
                                              'FLM.fit' = FLM),
                                  'FMM' = list(
                                              'fixed' = list(
                                                            'beta(t)' = beta2_a,
                                                            'delta(t)' = beta2_b,
                                                            'coefficient' = fixed.coeff2,
                                                            'vcov' = fixed.vcov2,
                                                            'sigma' = sigma2,
                                                            'se_beta(t)' = se2_a,
                                                            'se_delta(t)' = se2_b
                                                            ),
                                              'random' = list('Zi'=raneffect2[,1],
                                                              'beta_i(t)'=betai_2),
                                              'yhat' = yhat2,
                                              'FMM.fit' = FMM)
  )))
}
