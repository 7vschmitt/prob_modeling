source("DBDA2E-utilities.R")

want = c("mvtnorm")
have = want %in% rownames(installed.packages())
if ( any(!have) ) { install.packages( want[!have] ) }

try(library(mvtnorm))

# Specify the data, to be used in the likelihood function.
myData = c(rep(0,6),rep(1,14))

# Define the Bernoulli likelihood function, p(D|theta).
# The argument theta could be a vector, not just a scalar.
likelihood = function( theta , data ) {
  z = sum( data )
  N = length( data )
  pDataGivenTheta = theta^z * (1-theta)^(N-z)
  # The theta values passed into this function are generated at random,
  # and therefore might be inadvertently greater than 1 or less than 0.
  # The likelihood for theta > 1 or for theta < 0 is zero:
  pDataGivenTheta[ theta > 1 | theta < 0 ] = 0
  return( pDataGivenTheta )
}

# Define the prior density function. 
prior <- function(alpha = NULL, beta = NULL, psi = NULL, 
                  mu = NULL, theta_mu = NULL, Sigma_mu = NULL) {
  result <- list()
  if(!is.null(alpha)) {
    result$alpha <- foreach(i=iter(alpha), .combine = rbind) %dopar% { 
      dgamma(i, .5, .5) } %>% matrix(nrow=nrow(alpha))
  }
  if(!is.null(beta)) {
    result$beta <- foreach(i=iter(beta), .combine = c) %dopar% dgamma(i, .5, .5)
  }
  if(!is.null(psi)) {
    result$psi <- dmvnorm(psi, 
                    mean = rep(0, length(psi)), 
                    sigma = 10^4 * diag(nrow = length(psi)))
  }
  if(!is.null(mu)) {
    result$mu <- foreach(i=iter(mu, by = 'row'), .combine = rbind) %dopar% 
      exp(dmvnorm(i, mean = theta_mu, sigma = Sigma_mu))
  }
  if(!is.null(theta_mu)) {
    result$theta_mu <- dmvnorm(theta_mu, 
                        mean = rep(0, length(theta_mu)), 
                        sigma = 10^6 * diag(nrow = length(theta_mu)))
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

# Specify the length of the trajectory, i.e., the number of jumps to try:
trajLength = 50000 # arbitrary large number

user_cnt = 1200
# Initialize the vector that will store the results:
initals = list(mu=rep(c(0.05, 0.005, 0.005, 0.001), user_cnt) %>% 
                 matrix(ncol = 4, byrow = T), 
               alpha = alpha_mean, 
               beta = beta_mean, 
               psi = psi_mean,
               Sigma_mu = Sigma_mu_mean,
               theta_mu = theta_mu_mean)

trajectory <- list()
# Specify where to start the trajectory:
trajectory[[1]] <- initals
# Specify the burn-in period:
burnIn = ceiling( 0.2 * trajLength ) # arbitrary number, less than trajLength
# Initialize accepted, rejected counters, just to monitor performance:
nAlphaAccepted <- 0
nAlphaRejected <- 0
nBetaAccepted <- 0
nBetaRejected <- 0
nPsiAccepted <- 0
nPsiRejected <- 0


# Now generate the random walk. The 't' index is time or trial in the walk.
# Specify seed to reproduce same random walk:
set.seed(1899)

proposedParams <- function(current_params, ...) {
  change_params <- list(...)
  keys <- names(change_params)
  new_params <- current_params
  for(k in keys) new_params[k] <- change_params[k]
  new_params
}

for ( t in 1:(trajLength-1) ) {
  # initialize current position
  currentPosition <- trajectory[t]
  newPosition <- currentPosition
  # calc likelihoods
  current_lls <- foreach(i = 1:user_cnt, .combine = rbind) %dopar% {
    likelihood(i, currentPosition)
  }
  sum_ll <- prod(current_lls)
  # calc priors
  priors <- do.call(prior, currentPosition)
  
  #### draw alpha
  varAlpha <- 0.5 # Improvement: Adapt to have acceptance b/w 0.1 & 0.4
  alphaStar <- exp(foreach(a=iter(currentPosition$alpha), .combine = rbind) %dopar% {
    rnorm(1, mean = log(a), varAlpha) 
  } %>% matrix(nrow=nrow(currentPosition$alpha)))
  # Compute the probability of accepting the proposed jump.
  proposedAlpha <- proposedParams(currentPosition, alpha=alphaStar)
  alphaStarAccept <- min(1, 
                         sum(calc_likelihoods(1:user_cnt, proposedAlpha)) * 
                           prior(alpha = alphaStar) * prod(alphaStar) / 
                           sum_ll * priors$alpha * prod(currentPosition$alpha))
  
  # Generate a random uniform value from the interval [0,1] to
  # decide whether or not to accept the proposed jump.
  if (runif(1) < alphaStarAccept) {
    # accept the proposed jump
    newPosition <- proposedParams(newPosition, alpha=alphaStar)
    # increment the accepted counter, just to monitor performance
    if ( t > burnIn ) { nAlphaAccepted <- nAlphaAccepted + 1 }
  } else {
    # increment the rejected counter, just to monitor performance
    if ( t > burnIn ) { nAlphaRejected <- nAlphaRejected + 1 }
  }
  
  #### draw beta
  varBeta <- 0.5
  betaStar <- exp(foreach(b=iter(currentPosition$beta), .combine = c) %dopar% {
    rnorm(1, mean = log(a), varBeta) 
  })
  proposedBeta <- proposedParams(currentPosition, beta=betaStar)
  betaStarAccept <- min(1, 
                         sum(calc_likelihoods(1:user_cnt, proposedBeta)) * 
                           prior(beta = betaStar) * prod(betaStar) / 
                           sum_ll * priors$beta * prod(currentPosition$beta))
  
  if (runif(1) < betaStarAccept) {
    # accept the proposed jump
    newPosition <- proposedParams(newPosition, beta=betaStar)
    # increment the accepted counter, just to monitor performance
    if ( t > burnIn ) { nBetaAccepted <- nBetaAccepted + 1 }
  } else {
    # increment the rejected counter, just to monitor performance
    if ( t > burnIn ) { nBetaRejected <- nBetaRejected + 1 }
  }
  
  #### draw psi
  varPsi <- 0.5
  psiStar <- foreach(p = iter(currentPosition$psi), .combine = c) %dopar% {
    rnorm(1, mean = p, sd = varPsi)
  }
  proposedPsi <- proposedParams(currentPosition, psi = psiStar)
  psiStarAccept <- min(1, 
                       sum(calc_likelihoods(1:user_cnt, proposedPsi)) * 
                         prior(psi = psiStar) / 
                         sum_ll * priors$psi)
  
  if (runif(1) < psiStarAccept) {
    # accept the proposed jump
    newPosition <- proposedParams(newPosition, psi=psiStar)
    # increment the accepted counter, just to monitor performance
    if ( t > burnIn ) { nPsiAccepted <- nPsiAccepted + 1 }
  } else {
    # increment the rejected counter, just to monitor performance
    if ( t > burnIn ) { nPsiRejected <- nPsiRejected + 1 }
  }
  
  
  #### finally set new trajectory
  trajectory[ t+1 ] <- newPosition
}

# Extract the post-burnIn portion of the trajectory.
acceptedTraj = trajectory[ (burnIn+1) : length(trajectory) ]

# End of Metropolis algorithm.

#-----------------------------------------------------------------------
# Display the chain.

openGraph(width=4,height=8)
layout( matrix(1:3,nrow=3) )
par(mar=c(3,4,2,1),mgp=c(2,0.7,0))

# Posterior histogram:
paramInfo = plotPost( acceptedTraj$alpha[1,1] , xlim=c(0,1) , xlab=bquote(theta) , 
                      cex.main=2.0 ,
                      main=bquote( list( "Prpsl.SD" == .(proposalSD) ,
                                         "Eff.Sz." == .(round(effectiveSize(acceptedTraj),1)) ) ) )

# Trajectory, a.k.a. trace plot, end of chain:
idxToPlot = (trajLength-100):trajLength
plot( trajectory[idxToPlot]$alpha[1,1] , idxToPlot , main="End of Chain" ,
      xlab=bquote(alpha[1,1]) , xlim=c(0,1) , ylab="Step in Chain" ,
      type="o" , pch=20 , col="skyblue" , cex.lab=1.5 )
# Display proposal SD and acceptance ratio in the plot.
text( 0.0 , trajLength , adj=c(0.0,1.1) , cex=1.75 ,
      labels = bquote( frac(N[acc],N[pro]) == 
                         .(signif( nAlphaAccepted/length(acceptedTraj) , 3 ))))

# Trajectory, a.k.a. trace plot, beginning of chain:
idxToPlot = 1:100
plot( trajectory[idxToPlot] , idxToPlot , main="Beginning of Chain" ,
      xlab=bquote(theta) , xlim=c(0,1) , ylab="Step in Chain" ,
      type="o" , pch=20 , col="skyblue" , cex.lab=1.5 )
# Indicate burn in limit (might not be visible if not in range):
if ( burnIn > 0 ) {
  abline(h=burnIn,lty="dotted")
  text( 0.5 , burnIn+1 , "Burn In" , adj=c(0.5,1.1) )
}

#saveGraph( file=paste0( "MCMC" , 
#                        "SD" , proposalSD ,
#                        "Init" , trajectory[1] ) , type="eps" )

#------------------------------------------------------------------------
