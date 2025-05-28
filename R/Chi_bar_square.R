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
