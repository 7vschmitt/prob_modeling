# Define the Bernoulli likelihood function, p(D|theta).
# The argument theta could be a vector, not just a scalar.
likelihood <- function(i, params){
  llAlpha <- params$alpha
  llBeta <- params$beta
  llPsi <- params$psi
  llMu <- params$mu[i,]
  data_new_i <- data_new %>% filter(user == i)
  sum_k_log <- 0 # wird spaeter mit ergebnis aus k-loop multipliziert
  for (k in 1:length(ad_types)){
    prod_l <- 1 # wird sp?ter mit ergebnis aus l-loop multipliziert
    if (total_counts[i , k+1] != 0){
      first_two_lines <- numeric()
      for (l in 1:(total_counts[i , k+1])){ 
        first_line <- llMu[k] * exp(llPsi[k] * (data_new_i %>% 
                                                  filter(event == k) %>%
                                                  .[l, ] %>% 
                                                  select(purchase))) %>% 
          as.numeric() # erste Zeile
        
        sum_j <- numeric()
        for (j in 1:3){ # zweite Zeile erste summe 1:K-1
          # j schleife umgehen
          sum_index <- ifelse(data_new_i %>% 
                                filter(event == k) %>% 
                                plyr::empty(), #if is empty
                              0, #then
                              data_new_i %>% #else
                                filter(event == k) %>% 
                                filter(row_number() ==l) %>% 
                                select(j+3) %>% as.numeric()
          )
          if (sum_index > 0){
            sum_m <- numeric()
            t_l <- data_new_i %>% 
              filter(event==k) %>% 
              filter(row_number() == l) %>% 
              select(timestamp) %>% as.numeric()
            t_m <- data_new_i %>% filter(event==j) %>% 
              select(timestamp) 
            # second line:
            sum_m <- llAlpha[j,k] * exp(-llBeta[j] * (t_l - t_m[1:sum_index,]))
            
          } else {sum_m <- 0} 
          sum_j[j] <- sum(sum_m) # Summe j
        } # j Ende
        first_two_lines[l] <- first_line + sum(sum_j) # erste und zweite Zeile
      } 
    } else {
      first_two_lines <- 1
    }
    prod_l <- prod(first_two_lines)
    
    # Produkt ueber l; wenn von 1 bis 0, dann ergebnis auf 0 setzen
    sum_m <- 0
    purchase_cnt <- total_counts[i,5]
    if (purchase_cnt > 0){
      for (m in 0:(purchase_cnt-1)){ # 0 bis Anzahl Kaeufe in 120 Tagen -1
        t_m1 <- data_new_i %>% filter(event == K) %>% 
          filter(row_number() == m+1) %>% 
          select(timestamp) %>% 
          as.numeric()
        t_m <- ifelse(m == 0, 
                      0,
                      ifelse(purchase_cnt != 0, 
                             data_new_i %>% filter(event == K) %>% 
                               filter(row_number() == m) %>% 
                               select(timestamp) %>% as.numeric, 
                             0)
        )
        tmp_m <- exp(llPsi[k]*(m+1))*(t_m1 - t_m)
        sum_m <- sum_m + tmp_m # Summe m
      } # m Ende
      exponent_teil1 <- sum_m * llMu[k]
      
    } # if ende
    else {exponent_teil1 <- 0}
    browser(expr = is.na(exponent_teil1))
    exponent_teil2 <- 0
    for (j in 1:3){ # letzte Zeile erste Summe
      # j schleife umgehen?
      sum_m <- numeric()
      if (total_counts[i,j+1] > 0) {
        # letzte Zeile zweite Summe
        t_m_neu <- data_new_i %>% filter(event == j) %>% 
          select(timestamp) 
        tmp_m <- (llAlpha[j,k]/llBeta[j]) * (1-exp(-llBeta[j] * (T - t_m_neu)))
        sum_m <- sum(tmp_m) # Summe m
        #} # m Ende
        exponent_teil2 <- exponent_teil2 + sum_m
      } 
    } # j Ende
    
    third_fourth_line <- exponent_teil1 - exponent_teil2
    sum_k_log <- log(prod_l) + third_fourth_line + sum_k_log %>% as.numeric()
  }  # loop over all ad types
  return(sum_k_log)
} # Funktion Ende
# Define the prior density function. 
prior <- function(alpha = NULL, beta = NULL, psi = NULL, 
                  mu = NULL, theta_mu = NULL, Sigma_mu = NULL) {
  result <- list()
  if(!is.null(alpha)) {
    result$alpha <- foreach(i=iter(alpha), .combine = rbind) %dopar% { 
      dgamma(i, .5, .5, log = T) } %>% matrix(nrow=nrow(alpha))
  }
  if(!is.null(beta)) {
    result$beta <- foreach(i=iter(beta), .combine = c) %dopar% dgamma(i, 25, .5, log = T)
  }
  if(!is.null(psi)) {
    result$psi <- dmvnorm(psi, 
                          mean = rep(0, length(psi)), 
                          sigma = 10^4 * diag(nrow = length(psi)), log = T)
  }
  if(!is.null(mu)) {
    result$mu <- foreach(m=iter(mu, by = 'row'), .combine = rbind) %dopar% 
      #exp(dmvnorm(m, mean = theta_mu, sigma = Sigma_mu, log = T))
      dlnorm.rplus(m, theta_mu, Sigma_mu) %>% log
  }
  if(!is.null(theta_mu)) {
    result$theta_mu <- dmvnorm(theta_mu, 
                               mean = rep(0, length(theta_mu)), 
                               sigma = 10^6 * diag(nrow = length(theta_mu)), log = T)
  }
  if(!is.null(Sigma_mu)) {
    result$Sigma_mu <- diwish(Sigma_mu, S = diag(nrow = K), v = K)
  }
  return(result)
}

calc_likelihoods <- function(users, params) {
  ll <- foreach(i = iter(users), .combine = rbind) %dopar% likelihood(i, params)
  ll
}

makeSymmetric <- function(matrix) {
  return((matrix + t(matrix)) / 2)
}

proposedParams <- function(current_params, ...) {
  change_params <- list(...)
  keys <- names(change_params)
  new_params <- current_params
  for(k in keys) new_params[k] <- change_params[k]
  new_params
}

adjustVariance <- function(var, nAccept, nReject) {
  acceptRate <- nAccept / nReject
  if(acceptRate < 0.1) {
    var <- var * 0.75
  } else if(acceptRate > 0.4) {
    var <- var * 1.25
  }
  return(var)
}

# Specify the length of the trajectory, i.e., the number of jumps to try:
trajLength = 50000 # arbitrary large number

user_cnt = 25
# Initialize the vector that will store the results:
initals = list(mu=rep(exp(theta_mu_mean), user_cnt) %>% 
                 matrix(ncol = 4, byrow = T), 
               alpha = alpha_mean, 
               beta = beta_mean, 
               psi = psi_mean,
               Sigma_mu = Sigma_mu_mean,
               theta_mu = theta_mu_mean)
# Now generate the random walk. The 't' index is time or trial in the walk.
# Specify seed to reproduce same random walk:
set.seed(1899)
trajectory_file <- "mcmc_trajectory"
if(file.exists(trajectory_file)) {
  load(trajectory_file)
  startCount <- length(trajectory)
} else {
  trajectory <- list()
  # Specify where to start the trajectory:
  trajectory[[1]] <- initals
  startCount <- 1
}
# Specify the burn-in period:
burnIn = ceiling( 0.0 * trajLength ) # arbitrary number, less than trajLength
# Initialize accepted, rejected counters, just to monitor performance:
counts_file <- "mcmc_counts"
if(file.exists(counts_file)) {
  load(counts_file)
} else {
  nAlphaAccepted <- 0
  nAlphaRejected <- 0
  nBetaAccepted <- 0
  nBetaRejected <- 0
  nPsiAccepted <- 0
  nPsiRejected <- 0
  nMuAccepted <- numeric(length = user_cnt)
  nMuRejected <- numeric(length = user_cnt)
  tailAlphaAcc <- c(0, 0)
  tailBetaAcc <- c(0, 0)
  tailPsiAcc <- c(0, 0)
  tailMuAcc <- c(0, 0)
  varAlpha <- alpha_sd * 1.5
  varBeta <- beta_sd * 1.5
  varPsi <- psi_sd * 1.5
  varMu <- rep(0.001, K)
}


paste(Sys.time(), "Start MCMC sampling") %>% print
for ( ct in startCount:(trajLength-1) ) {
  # initialize current position
  currentPosition <- trajectory[[ct]]
  newPosition <- currentPosition
  # calc likelihoods
  current_lls <- foreach(i = 1:user_cnt, .combine = rbind) %dopar% {
    likelihood(i, currentPosition)
  }
  sum_ll <- sum(current_lls)
  # calc priors
  priors <- do.call(prior, currentPosition)
  
  #### draw alpha
  
  if(ct %% 10 == 0) {
    varAlpha <- adjustVariance(varAlpha, nAlphaAccepted - tailAlphaAcc[1], nAlphaRejected - tailAlphaAcc[2])
    tailAlphaAcc <- c(nAlphaAccepted, nAlphaRejected)
  }
  # alphaStar <- foreach(a=iter(c(currentPosition$alpha)), .combine = c) %dopar% {
  #   exp(rnorm(1, mean = log(a), sd = varAlpha))
  # } %>% matrix(nrow=nrow(currentPosition$alpha))
  alphaStar <- foreach(j = 1:nrow(currentPosition$alpha), .combine = rbind) %:% 
    foreach(k = 1:ncol(currentPosition$alpha)) %dopar% {
      rlnorm(1, meanlog = log(currentPosition$alpha[j,k]), 
             sdlog = sqrt(log(1 + varAlpha[j,k]^2 / currentPosition$alpha[j,k])))    
  } %>% unlist() %>% matrix(nrow=nrow(currentPosition$alpha))
  # Compute the probability of accepting the proposed jump.
  proposedAlpha <- proposedParams(currentPosition, alpha=alphaStar)
  # exp(log(proposed) - log(current)) = proposed/current
  probAlpha <- exp((sum(calc_likelihoods(1:user_cnt, proposedAlpha)) + 
                  sum(prior(alpha = alphaStar)$alpha)) - 
    (sum_ll + sum(priors$alpha)))
  probAlpha <- ifelse(is.nan(probAlpha),0 , probAlpha)
  alphaStarAccept <- min(1, probAlpha)
  
  # Generate a random uniform value from the interval [0,1] to
  # decide whether or not to accept the proposed jump.
  if (runif(1) < alphaStarAccept) {
    # accept the proposed jump
    newPosition$alpha <- alphaStar
    # increment the accepted counter, just to monitor performance
    if ( ct > burnIn ) { nAlphaAccepted <- nAlphaAccepted + 1 }
  } else {
    # increment the rejected counter, just to monitor performance
    if ( ct > burnIn ) { nAlphaRejected <- nAlphaRejected + 1 }
  }
  
  #### draw beta
  if(ct %% 10 == 0) {
    varBeta <- adjustVariance(varBeta, nBetaAccepted - tailBetaAcc[1], nBetaRejected - tailBetaAcc[2])
    tailBetaAcc <- c(nBetaAccepted, nBetaRejected)
  }
  # betaStar <- foreach(b=iter(currentPosition$beta), .combine = c) %dopar% {
  #   exp(rnorm(1, mean = log(b), sd = varBeta))
  # }
  betaStar <- foreach(j = 1:length(currentPosition$beta), .combine = c) %dopar% {
    rlnorm(1, meanlog = log(currentPosition$beta[j]), 
           sdlog = sqrt(log(1 + varBeta[j]^2/currentPosition$beta[j]^2)))
  }
  proposedBeta <- proposedParams(currentPosition, beta=betaStar)
  probBeta <- exp((sum(calc_likelihoods(1:user_cnt, proposedBeta)) + 
                 sum(prior(beta = betaStar)$beta)) - 
    (sum_ll + sum(priors$beta)))
  probBeta <- ifelse(is.nan(probBeta), 0, probBeta)
  betaStarAccept <- min(1, probBeta)
  if (runif(1) < betaStarAccept) {
    # accept the proposed jump
    newPosition$beta <- betaStar
    # increment the accepted counter, just to monitor performance
    if ( ct > burnIn ) { nBetaAccepted <- nBetaAccepted + 1 }
  } else {
    # increment the rejected counter, just to monitor performance
    if ( ct > burnIn ) { nBetaRejected <- nBetaRejected + 1 }
  }
  
  #### draw psi
  acceptRate <- nPsiAccepted / nPsiRejected
  if(ct %% 10 == 0) {
    varPsi <- adjustVariance(varPsi, nPsiAccepted - tailPsiAcc[1], nPsiRejected - tailPsiAcc[2])
    tailPsiAcc <- c(nPsiAccepted, nPsiRejected)
  }
  psiStar <- foreach(k = 1:length(currentPosition$psi), .combine = c) %dopar% {
    rnorm(1, mean = currentPosition$psi[k], sd = varPsi[k])
  }
  proposedPsi <- proposedParams(currentPosition, psi = psiStar)
  probPsi <- exp((sum(calc_likelihoods(1:user_cnt, proposedPsi)) + prior(psi = psiStar)$psi) - 
    (sum_ll + priors$psi))
  probPsi <- ifelse(is.nan(probPsi), 0, probPsi)
  psiStarAccept <- min(1, probPsi)
  if (runif(1) < psiStarAccept) {
    # accept the proposed jump
    newPosition$psi <- psiStar
    # increment the accepted counter, just to monitor performance
    if ( ct > burnIn ) { nPsiAccepted <- nPsiAccepted + 1 }
  } else {
    # increment the rejected counter, just to monitor performance
    if ( ct > burnIn ) { nPsiRejected <- nPsiRejected + 1 }
  }
  
  #### draw mu^i
  if(ct %% 10 == 0) {
    varMu <- adjustVariance(varMu, sum(nMuAccepted) - tailMuAcc[1], sum(nMuRejected) - tailMuAcc[2])
    tailMuAcc <- c(sum(nMuAccepted), sum(nMuRejected))
  }
  foreach(i = 1:user_cnt, .combine=rbind, .export = c("newPosition", "nMuAccepted", "nMuRejected")) %do% {
    muiStar <- foreach(k = 1:length(currentPosition$mu[i,]), .combine = cbind) %do% {
      rlnorm(1, meanlog = log(currentPosition$mu[i,k]), 
             sdlog = sqrt(log(1 + varMu^2 / currentPosition$mu[i,k]^2)))
    }
    propMui <- currentPosition$mu
    propMui[i, ] <- muiStar
    proposedMu <- proposedParams(currentPosition, mu = propMui)
    probMui <- exp((calc_likelihoods(i, proposedMu) + 
                  prior(mu = muiStar, 
                        theta_mu = currentPosition$theta_mu, 
                        Sigma_mu = currentPosition$Sigma_mu)$mu) - 
      (current_lls[i] + priors$mu[i]))
    probMui <- ifelse(is.nan(probMui), 0, probMui)
    muiStarAccept <- min(1, probMui)
    
    if (runif(1) < muiStarAccept) {
      # accept the proposed jump
      newPosition$mu[i, ] <- muiStar
      # increment the accepted counter, just to monitor performance
      if ( ct > burnIn ) { nMuAccepted[i] <- nMuAccepted[i] + 1 }
    } else {
      # increment the rejected counter, just to monitor performance
      if ( ct > burnIn ) { nMuRejected[i] <- nMuRejected[i] + 1 }
    }
  }
  
  #### draw theta_mu
  I <- diag(nrow = nrow(currentPosition$Sigma_mu))
  Sigma_theta_mu <- 10^6 * I
  InvSigmaMu <- solve(currentPosition$Sigma_mu)
  InvSigmaThetaMu <- solve(Sigma_theta_mu)
  #B <- solve(I %*% InvSigmaMu + InvSigmaThetaMu)
  #A <- t(B) %*% t( t(colMeans(log(currentPosition$mu))) %*% InvSigmaMu + rep(0, K) %*% InvSigmaThetaMu )
  # Gelman formula p. 71
  
  B <- solve(user_cnt * InvSigmaMu + InvSigmaThetaMu)
  A <- B %*% (user_cnt * InvSigmaMu %*% colMeans(log(currentPosition$mu)) + InvSigmaThetaMu %*% rep(0, K))
  if(!isSymmetric(B)) B <- makeSymmetric(B)
  
  thetaMuStar <- rmvnorm(1, mean = A, sigma = B) %>% as.numeric()
  newPosition$theta_mu <- thetaMuStar
  
  #### draw Sigma_mu
  S <- solve(
    foreach(i = 1:nrow(currentPosition$mu), .combine = '+') %dopar% {
      (log(currentPosition$mu[i, ]) - currentPosition$theta_mu) %*% 
        t(log(currentPosition$mu[i, ]) - currentPosition$theta_mu)
    } + diag(K)
  )
  v <- 1 + user_cnt
  SigmaMuStar <- riwish(S = S, v = v)
  # avoid assymetry through bad rounding
  if(!isSymmetric(SigmaMuStar)) SigmaMuStar <- makeSymmetric(SigmaMuStar)
  newPosition$Sigma_mu <- SigmaMuStar
  
  #### finally set new trajectory
  trajectory[[ct+1]] <- newPosition
  
  # save milestones
  if (ct %% 100 == 0) {
    paste(Sys.time(), "Reached iteration", ct) %>% print
    save(trajectory, file = trajectory_file)
    save(nAlphaAccepted,
         nAlphaRejected,
         nBetaAccepted,
         nBetaRejected,
         nPsiAccepted,
         nPsiRejected,
         nMuAccepted,
         nMuRejected, 
         tailAlphaAcc,
         tailBetaAcc,
         tailPsiAcc,
         tailMuAcc,
         varAlpha,
         varBeta,
         varPsi,
         varMu,
         file = counts_file)
  } 
}

# Extract the post-burnIn portion of the trajectory.
acceptedTraj = trajectory[ (burnIn+1) : length(trajectory) ]

# End of Metropolis algorithm.

#-----------------------------------------------------------------------
# Display the chain.

# openGraph(width=4,height=8)
# layout( matrix(1:3,nrow=3) )
# par(mar=c(3,4,2,1),mgp=c(2,0.7,0))
# 
# # Posterior histogram:
# paramInfo = plotPost( acceptedTraj$alpha[1,1] , xlim=c(0,1) , xlab=bquote(theta) , 
#                       cex.main=2.0 ,
#                       main=bquote( list( "Prpsl.SD" == .(proposalSD) ,
#                                          "Eff.Sz." == .(round(effectiveSize(acceptedTraj),1)) ) ) )
# 
# # Trajectory, a.k.a. trace plot, end of chain:
# idxToPlot = (trajLength-100):trajLength
# plot( trajectory[idxToPlot]$alpha[1,1] , idxToPlot , main="End of Chain" ,
#       xlab=bquote(alpha[1,1]) , xlim=c(0,1) , ylab="Step in Chain" ,
#       type="o" , pch=20 , col="skyblue" , cex.lab=1.5 )
# # Display proposal SD and acceptance ratio in the plot.
# text( 0.0 , trajLength , adj=c(0.0,1.1) , cex=1.75 ,
#       labels = bquote( frac(N[acc],N[pro]) == 
#                          .(signif( nAlphaAccepted/length(acceptedTraj) , 3 ))))
# 
# # Trajectory, a.k.a. trace plot, beginning of chain:
# idxToPlot = 1:100
# plot( trajectory[idxToPlot] , idxToPlot , main="Beginning of Chain" ,
#       xlab=bquote(theta) , xlim=c(0,1) , ylab="Step in Chain" ,
#       type="o" , pch=20 , col="skyblue" , cex.lab=1.5 )
# # Indicate burn in limit (might not be visible if not in range):
# if ( burnIn > 0 ) {
#   abline(h=burnIn,lty="dotted")
#   text( 0.5 , burnIn+1 , "Burn In" , adj=c(0.5,1.1) )
# }
# 
# #saveGraph( file=paste0( "MCMC" , 
# #                        "SD" , proposalSD ,
# #                        "Init" , trajectory[1] ) , type="eps" )
# 
# #------------------------------------------------------------------------
