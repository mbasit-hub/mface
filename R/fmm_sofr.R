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
