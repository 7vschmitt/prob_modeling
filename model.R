library('MASS')
#library(mixAK)
#install.packages("mixAK")


n <- 200 #12000
k <- 4

aBarAlpha = 0.5
bBarAlpha = 0.5

aBarBeta = 0.5
bBarBeta = 0.5

alpha.k <- rgamma(n = n, shape = aBarAlpha, rate = bBarAlpha) # do that for each type combination ?
beta.k <- rgamma(n = n, shape = aBarBeta, rate = bBarBeta) # do this for all K-1 types
psi <- mvrnorm(n = n, mean = thetaBarPsi, sd = sigmaBarPsi)


thetaBarPsi # ?
sigmaBarPsi # ?

#alpha_j_k: mutually exciting effects between different types of points
alpha_j_k <- matrix(c(2.817, 0.0860, 0.5381, 0.6167,
                      0.1614, 1.7818, 0.2055, 0.0845,
                      0.4647, 0.1270, 8.0526, 0.8384), nrow = 3, ncol = 4)
row.names(alpha_j_k) <- c("search", "display", "other")
colnames(alpha_j_k) <- c("search", "display", "other", "purchase")

# beta_j: decaying factor
beta_j <- c(34.0188, 46.8854, 51.5114)
names(beta_j) <- c("search", "display", "other")

#psi_k: influence of a purchase on the probability of future occurence of type k points
psi_k <- c(-0.5664, -0.7556, -0.6235, 0.2787)
names(psi_k) <- c("search", "display", "other", "purchase")

# theta_mu: mean over all 12,000 theta_mu_i => input paramter for mu_i
theta_mu <- c(-5.3926, -6.1027, -5.8063, -9.7704) 
names(theta_mu) <- c("search", "display", "other", "purchase")

#Sigma_mu: mean over all 12,000 Sigma_mu_i: input paramter for mu_i
Sigma_mu <- matrix(c(0.4584, -0.1197, -0.4942, 0.2256,
                     -0.1197, 0.5934, -0.3383, -0.4157,
                     -0.4942, -0.3380, 1.0014, 0.0762,
                     0.2256, -0.4157, 0.0762, 2.3914), nrow = 4, ncol = 4)
row.names(Sigma_mu) <- c("search", "display", "other", "purchase")
colnames(Sigma_mu) <- c("search", "display", "other", "purchase")


# mu: sample from l
mu <- mvrnorm(n, mu = theta_mu, Sigma = Sigma_mu)
mu <- exp(mu)

lambda_k <- matrix(nrow = n, ncol = k) # 200 * 4 matrix

N_k_i <- function(t) # for k=3
  N_k_i_init <- # initialize N_K_i with t=0 (?) and K=3
  
  t_l_ji <- function(l,j,i)  #?
    return 

cond_intens_func <- function(t){
  for (i in c(1:n)){
    for (k in c(1:k)){
      for (j in c(1:3)) # 3 = K-1
        for (l in c(1:N_k_i))
          lambda_k <- mu[i,k]*exp(psi_k[k]*N_k_i(t)[i,3]) #k=3 
        + alpha_j_k[j,k]*exp(-beta_j[j]*(t-t_l_ji(l,j,i)))
        return(lambda_k)
    }
  }
}