library('MASS')

n <- 12000
k <- 4

aBarAlpha # ?
bBarAlpha # ?

aBarBeta # ?
bBarBeta # ?

thetaBarPsi # ?
sigmaBarPsi # ?

alpha.k <- rgamma(n = n, shape = aBarAlpha, rate = bBarAlpha) # do that for each type combination ?
beta.k <- rgamma(n = n, shape = aBarBeta, rate = bBarBeta) # do this for all K-1 types
psi <- mvrnorm(n = n, mean = thetaBarPsi, sd = sigmaBarPsi)