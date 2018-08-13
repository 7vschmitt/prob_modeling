try( library(rjags) )
try( library(runjags) )
try( runjags.options( inits.warning=FALSE , rng.warning=FALSE ) )

# set default number of chains and parallelness for MCMC:
library(parallel) # for detectCores().
nCores = detectCores() 
if ( !is.finite(nCores) ) { nCores = 1 } 
if ( nCores > 4 ) { 
  nChainsDefault = 4  # because JAGS has only 4 rng's.
  runjagsMethodDefault = "parallel"
}
if ( nCores == 4 ) { 
  nChainsDefault = 3  # save 1 core for other processes.
  runjagsMethodDefault = "parallel"
}
if ( nCores < 4 ) { 
  nChainsDefault = 3 
  runjagsMethodDefault = "rjags" # NOT parallel
}

modelStr <- "
data {
C <- 10000 # constant for keeping scaled lambda < 1
}
model {
# lambda
for(i in 1:nSubj){
ones[i] <- 1
lambda[i] <- ()/C
ones[i] ~ dbern(lambda[i]) # use bernulli ones trick (Kruschke, 2015, p. 214)
}

# alpha
for(j in 1:K-1) {
for(k in 1:K) {
alpha[j, k] ~ dgamma(0.5, 0.5)
}
}

# beta
for(j in 1:K-1) {
beta[j] ~ dgamma(0.5, 0.5)
}

# psi
psi ~ dmnorm(0, 10^4 * identMat)

# mu
for(i in 1:nSubj){
mu[i] ~ exp(dmnorm(theta_mu, Sigma_mu))
}

# theta_mu
theta_mu ~ dmnorm(0, 10^6 * identMat)
# for ( varIdx in 1:K ) { theta_mu[varIdx] ~ dnorm( 0 , 10^6 ) }  # alternative

# Sigma_mu
Sigma_mu_inv ~ dwish(scaleMat)
Sigma_mu <- inverse(Sigma_mu_inv)
}
"
## Run JAGS
writeLines( modelString , con="TEMPmodel.txt" )


dataList = list(    # Put the information into a list.
  userCounts = data %>% 
    select(user, event) %>% 
    mutate(count = as.integer(1)) %>% 
    group_by(user, event) %>% 
    summarise(count = sum(count)) %>% 
    spread(key = event, value = count, fill = as.integer(0)),
  identMat = diag(1, K),
  scaleMat = diag(1, K),
  nSubj = data$user %>% unique %>% length
)
# Initialize the chains based on MLE of data.
# Option: Use single initial value for all chains:
#  thetaInit = sum(y)/length(y)
#  initsList = list( theta=thetaInit )
# Option: Use function that generates random values for each chain:
initsList = function() {
  resampledY = sample( y , replace=TRUE )
  thetaInit = sum(resampledY)/length(resampledY)
  thetaInit = 0.001+0.998*thetaInit # keep away from 0,1
  return( list( theta=thetaInit ) )
}

# Run the chains:
jagsModel = jags.model( file="TEMPmodel.txt" , data=dataList , inits=initsList , 
                        n.chains=3 , n.adapt=500 )
update( jagsModel , n.iter=500 )
codaSamples = coda.samples( jagsModel , variable.names=c("theta") ,
                            n.iter=3334 )
save( codaSamples , file=paste0(fileNameRoot,"Mcmc.Rdata") )
