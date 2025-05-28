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
