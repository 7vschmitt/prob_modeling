#Making Plate notation plots

setwd("C:/Users/User/PycharmProjects/prob_modeling")
source("plot_dist.R")

#############################################
#Full hierarchical model
#############################################

#level 1:
#mu|theta, Sigma ~ log MVN dependensies (2)
plot_dist(dists$log_mvn, labels = c(mean = expression(mu[i]), right_sd = expression(Sigma[mu]), left_sd = expression(Theta[mu])))


######################################################
#level 2:

#alpha_j_k ~ Gamma
plot_dist(dists$gamma, labels = c(params = expression(alpha[jk](bar(a)[alpha], bar(b)[alpha]))))
#beta_j ~ Gamma
plot_dist(dists$gamma, labels = c(params = expression(beta[j](bar(a)[beta], bar(b)[beta]))))


#psi_k ~ MVN (no dep)
plot_dist(dists$mvn_normal, labels = c(mean = expression(psi), right_sd = expression(bar(Sigma)[mu]), left_sd = expression(bar(Theta)[mu])))

#level 3:

#1. theta_mu ~ MVN
plot_dist(dists$mvn_normal, labels = c(mean = expression(Theta[mu])))
#Sigma_mu ~IW 
plot_dist(dists$inv_wishart, labels = c(mean = expression(Sigma[mu])))







