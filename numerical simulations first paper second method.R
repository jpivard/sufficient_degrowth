########## A model of consumption choices and sufficiency with environmental and social image preferences - Second simulation method ###########

library(lattice)
library(tidyverse)
library (gtools)
library(RColorBrewer)

######## INTRODUCTION ###############


rm(list=ls())

#This script is designed for numerical simulations associated to the model presented here : https://fr.overleaf.com/project/653a7f038e83eb35864f2593 
# It aims to calculate the consumption areas and quantities of different goods. 
#These goods are defined along two dimensions: green/brown, conspicuous/discreet.
#Depending on their personal preferences, agents can choose to consume one or several goods.


#Unlike in the first method, here we aim to take into account all types of consumers, regardless of the number of goods they consume.
#To do so, we will use more generic formulas to characterize preferences.
#However we still rely on rasterization for our computations.

#Matrix size is defined as a function of the step between 2 pixels.
step = 0.025
size = 1/step + 1
size

#Let's also define our preference parameters, over which we will iterate
alphas <- seq(from = 0,to = 1,by=step)
betas_rev <- seq(from = 1,to = 0,by=-step) 
betas <- seq(from = 0,to = 1,by=step) 


######## POPULATION DISTRIBUTIONS  ##############


#To start with (scenario 1), we generate a uniformly distributed population density matrix of this size, giving the proportions of individuals positioned at each point on the map. 
n1 <- size
n2 <- size
pop <- matrix (rep(1/(size*size), n1*n2), n1, n2)
pop
#For example, by choosing a step size of 0.2, we have a matrix of size 6*6, with 36 coefficients ; while it is 11*11 (and 101 coeffs) for a 0.1 step (the default choice in the first place)

#Scenario 2 : see first script



###### BORDER PARAMETRIZATION #################

#We define our "borders", delimiting the consumption zones according to the exogenous model parameters (so that we can study the different cases by changing only the following values).

#we initialize the parameters of our "borders", i.e. the prices of the four goods, individual income and marginal damage as well as impact parameters (satisfying model assumptions)
p_go = 3
p_gd = 2
p_bo = 1.99
p_bd = 1  #this price is fixed
R = 20
d_prime = 0.01
gamma_bo = 4
gamma_bd = 2
gamma_go = 1.6
gamma_gd = 0.8



##### CONSUMPTION TRADEOFFS (Generic expressions) #####


### Partial derivatives appearing in the first order conditions of the maximization problem 
### Or even start from the Lagrangian ? 







###One good

#Let's start by good GO (and then we'll write a more generic code)


