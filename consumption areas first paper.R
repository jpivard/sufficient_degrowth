########## A model of consumption choices and sufficiency with environmental and social image preferences ###########


library(lattice)
library(tidyverse)
#install.packages("rasterVis") : not working 
#library (rasterVis) 
#install.packages("gtools") 
library (gtools)
library(RColorBrewer)
#install.packages("copula")
#install.packages("rgl")
#install.packages("psych")


#library(copula)  : probl?me de package
#install.packages("gsl") : idem
#library(psych) : toujours la m?me...

######## INTRODUCTION ###############


rm(list=ls())

#This script is designed for numerical simulations associated to the model presented here : https://fr.overleaf.com/project/653a7f038e83eb35864f2593 
# It aims to calculate the consumption areas and quantities of different goods. 
#These goods are defined along two dimensions: green/brown, conspicuous/discreet.
#Depending on their personal preferences, agents can choose to consume one or several goods.




#Method used in this script : As our heterogeneous consumers are placed in a [0,1]*[0,1] grid according to their preferences (on the x-axis, environmental awareness; on the y-axis, concern for social image), 
#The idea behind this numerical method is to discretize/rasterize the square into several points/pixels (expressing the two types of preferences in percentages if the number of points is 101 on the square ; for example individuals lying at the origin being indifferent to both, while those at the top right hand corner have maximum sensitivities to both).


#Matrix size is defined as a function of the step between 2 pixels.
step = 0.1
size = 1/step + 1
size


#Let's also define our preference parameters, over which we will iterate
alphas <- seq(from = 0,to = 1,by=step)
betas_rev <- seq(from = 1,to = 0,by=-step) #reversed so that coordinates appear in the right order
betas <- seq(from = 0,to = 1,by=step) #only for graphs

######## POPULATION DISTRIBUTIONS  ##############

#Here, square matrices of dimension 'size' are defined. 


#To start with (scenario 1), we generate a uniformly distributed population density matrix of this size, giving the proportions of individuals positioned at each point on the map. 
n1 <- size
n2 <- size
pop <- matrix (rep(1/(size*size), n1*n2), n1, n2)
pop
#For example, by choosing a step size of 0.2, we have a matrix of size 6*6, with 36 coefficients ; while it is 11*11 (and 101 coeffs) for a 0.1 step (the default choice in the first place)



#Scenario 2: population concentration around mean/median values
#This case is modeled by a product of independent beta laws for each preference parameter.

pop2 <- matrix(,n1,n2) #an empty matrix
 
#Basic scenario 2: two beta laws with the two same shape parameters (bell curves) to model the concentration in the middle of the distribution
a=b=c=d=2

#Variant (2b): the shape parameters are changed (the beta parameter is now assumed to follow a uniform distribution, i.e. a beta distribution with c=d=1).
#a=b=2
#c=d=1


#Scenario 2B: Lower/laxer environmental norm (concentration on the left, beta law with shape parameters 5 and 1 on alpha, uniform on beta).
#a=5
#b=c=d=1

#Sc?nario 2C : Higher/stringent environmental norm (right-hand concentration)
#a=c=d=1
#b=5

#Sc?nario 2D : social image matters less (widespread discretion)
#a=b=c=1
#d=5

#Sc?nario 2E: social image matters more ("Red Queen Economy", Keep Up With the Joneses...)
#a=b=d=1
#c=5




###Code that computes discrete densities###


#We create our quantile vectors (in line)

alphas_vector <- t(matrix(alphas))
betas_vector <- t(matrix(betas))

#Remark : Repartition functions can be represented (the same ones, with identical shape parameters)
fdr_beta_dim1 <-pbeta(alphas, shape1 = a, shape2 =b)
plot(fdr_beta_dim1,alphas)

fdr_beta_dim2 <-pbeta(betas, shape1 = c, shape2 =d)
plot(fdr_beta_dim2,betas)


#Starting from the repartition function, we define the discrete density of a beta (one for each of the 2 dimensions).
proba_values_dim1 <-matrix(,1,size)
proba_values_dim2 <-matrix(,1,size)

for (n in 1:(1/step)-1:1) {
  densite_beta_dim1 <- function(n) {
  densite <-  pbeta(step/2 + 1-n*step,a,b) - pbeta(step/2 + 1-(n+1)*step,a,b)
  return(densite)
  }
  for (i in 1:ncol(proba_values_dim1)-1){
   proba_values_dim1[1,i+1]= densite_beta_dim1(i)  #values are iteratively stored in a matrix
  }
}

#first and last values are defined outside the loop
densite_beta_dim1(0) <- pbeta (1,a,b) - pbeta (1-step/2,a,b)
densite_beta_dim1(1/step)<- pbeta(step/2,a,b) - pbeta(0,a,b)

proba_values_dim1[1,1]=densite_beta_dim1(0)
proba_values_dim1[1,size]=densite_beta_dim1(1/step)

sum(proba_values_dim1)
#We get the values we need (summing to 1, bell-shaped density, symmetrical)


#plot(densite_beta_dim1,n)  
#This gives a continuous density (whereas finite function only for integers...) instead of a discrete one, but has the right shape
#Density should be a function of alpha, not n
#An interval is missing to have 11 probability values


#Idem for the betas (for the moment we use exactly the same probability laws)
for (n in 1:(1/step)-1:1) {
  densite_beta_dim2 <- function(n) {
    densite <-  pbeta(step/2 + 1-n*step,c,d) - pbeta(step/2 + 1-(n+1)*step,c,d)
    return(densite)
  }
  for (i in 1:ncol(proba_values_dim2)-1){
    proba_values_dim2[1,i+1]= densite_beta_dim2(i)  
  }
}

densite_beta_dim2(0) <- pbeta (1,c,d) - pbeta (1-step/2,c,d)   
densite_beta_dim2(1/step)<- pbeta(step/2,c,d) - pbeta(0,c,d)

proba_values_dim2[1,1]=densite_beta_dim2(0)
proba_values_dim2[1,size]=densite_beta_dim2(1/step)

#We multiply our two probability vectors to obtain the probabilities for the entire population.
pop2 = t(proba_values_dim2)%*%proba_values_dim1
#There is indeed a higher population concentration on median values compared to the uniform case.

sum(pop2) 
#just checking, we're good :)
        

#Scenario 3: Correlated Beta laws (positively or negatively) - come back at it later

#M?thode des copules ?
#Th?orie : https://www.r-bloggers.com/2011/10/copulas-made-easy/
#Docu :
#https://www.r-bloggers.com/2015/10/modelling-dependence-with-copulas-in-r/
#https://www.r-bloggers.com/2016/03/how-to-fit-a-copula-model-in-r-heavily-revised-part-1-basic-tools/


#3a : On commence par mod?liser une forte corr?lation positive entre nos deux types de sensibilit?s

#Premi?re ?tape : on g?n?re une beta bivari?e avec corr?lation positive entre les marginales

#require(mvtnorm)
#set.seed(235)

#S <- matrix(c(1,.8,.8,1),2,2) #Une matrice de corr?lation avec des coefficients hors diagonale proches de 1
#AB_copula <- rmvnorm(mean=c(0.5,0.5),sig=S,n=size*size) #On d?finit une normale bivari?e centr?e autour de la m?diane de l'intervalle et avec le bon nombre d'observations

#U <- pnorm(AB_copula) #On convertit en uniforme - se v?rifie avec hist(U[,1]) or hist(U[,2])

#alpha_cop <- qbeta(U[,1],2,2) 
#beta_cop <- qbeta(U[,2],2,2) 
#On choisit ici comme lois marginales deux Beta de m?me param?tres de forme, correspondant ? une concentration des valeurs vers le centre

#plot(alpha_cop,beta_cop)  #On obtient notre r?partition de population en fonction des param?tres, qu'il faut maintenant convertir en matrice de probabilit?s
#On obtient la densit? de cette distribution bivari?e ? l'aide d'une fonction du package :

#library(MASS)




############# BORDER PARAMETRIZATION ####################

#We define our "borders", delimiting the consumption zones according to the exogenous model parameters (so that we can study the different cases by changing only the following values).

#we initialize the parameters of our "borders", i.e. the prices of the four goods, individual income and marginal damage as well as impact parameters (satisfying model assumptions)
p_go = 3
p_gd= 2
p_bo = 1.99
p_bd= 1  #this price is fixed
R = 20
d_prime = 0.01
gamma_bo = 4
gamma_bd = 2
gamma_go = 1.6
gamma_gd = 0.8



######## COMPUTATION OF CONSUMPTION ZONE MATRICES (first without quantities?) ###############


#1st method (inspired from the simpler model):
#We run a nested loop on our individuals (i.e. their alpha and beta parameters, which we will convert - by normalizing - into (x,y) coordinates where the coordinates lie between 1 and size).
#For each point (i.e. pair of parameters) we look at which of the 12 constraints are satisfied, which determines the type of good consumed at that point (and the corresponding matrix among the 4 defined above, which will be filled by a "1" at the associated (x,y) point where the good is exclusively consumed).
#The idea is to create four size*size boolean matrices (true = 1 /false = 0), one for each type of consumer; and to examine for each point on the square (using nested loops) whether it satisfies the constraints corresponding to a certain type.

#We start by exclusive consumption zones only (defining the 4 corresponding matrices)
#Then we compute the mixed consumption zones matrices 

#In the end, we'll check that the sum of these four matrices gives 1s in the exclusive consumption zones and integers from 2 to 4 elsewhere.

#1.Exclusive consumption zones

#initialization of the matrices : we create empty matrices of the right dimension
mat_go_only <- matrix(, nrow = size, ncol = size)  
mat_gd_only <- matrix(, nrow = size, ncol = size)
mat_bo_only <- matrix(, nrow = size, ncol = size)
mat_bd_only <- matrix(, nrow = size, ncol = size)
test<-matrix(, nrow = size, ncol = size)

#nested loop to qualify each constraint everywhere in the map

for (beta in betas_rev)  {
  
  for (alpha in alphas) {
    
    x=size-beta/step
    y=alpha/step+1
    
    #1. Preferences for good GO 
    
    #Border f1 (coefficient takes value one if individual lies above, meaning she prefers good GO over GD)
    if ((1-(p_go/p_gd)*alpha)*beta>(((p_go/p_gd)-1)*(p_go/R)+d_prime*alpha*(gamma_go-gamma_gd*(p_go/p_gd)))) {
      bool_f1 = 1
    }
    else {bool_f1 = 0}
    
    #f2
    if ((1-(p_go/p_bo)*(1-alpha))*beta>(((p_go/p_bo)-1)*(p_go/R)+d_prime*alpha*(gamma_go-gamma_bo*(p_go/p_bo)))) {
      bool_f2 = 1
    }
    else {bool_f2 = 0}
    
    #f3
    if (beta> ((p_go-1)*(p_go/R)+d_prime*alpha*(gamma_go-gamma_bd*p_go)) ) {
      bool_f3 = 1
    }
    else {bool_f3 = 0}
    
    #2. Preferences for good GD
    
    if ((alpha-(p_gd/p_go))*beta>(((p_gd/p_go)-1)*(p_gd/R)+d_prime*alpha*(gamma_gd-gamma_go*(p_gd/p_go)))) {
      bool_g1 = 1
    }
    else {bool_g1 = 0}
    
    if ((alpha-(p_gd/p_bo)*(1-alpha))*beta>(((p_gd/p_bo)-1)*(p_gd/R)+d_prime*alpha*(gamma_gd-gamma_bo*(p_gd/p_bo)))) {
      bool_g2 = 1
    }
    else {bool_g2 = 0}
    
    if (alpha*beta> ((p_gd-1)*(p_gd/R)+d_prime*alpha*(gamma_gd-gamma_bd*p_gd)) ) {
      bool_g3 = 1
    }
    else {bool_g3 = 0}
    
    
    #3. Preferences for good BO
    
    if ((1-alpha-(p_bo/p_go))*beta>(((p_bo/p_go)-1)*(p_bo/R)+d_prime*alpha*(gamma_bo-gamma_go*(p_bo/p_go)))) {
      bool_h1 = 1
    }
    else {bool_h1 = 0}
    
    if ((1-alpha*(p_bo/p_gd))*beta>(((p_bo/p_gd)-1)*(p_bo/R)+d_prime*alpha*(gamma_bo-gamma_gd*(p_gd/p_bo)))) {
      bool_h2 = 1
    }
    else {bool_h2 = 0}
    
    if ((1-alpha)*beta> ((p_bo-1)*(p_bo/R)+d_prime*alpha*(gamma_bo-gamma_bd*p_bo)) ) {
      bool_h3 = 1
    }
    else {bool_h3 = 0}
    
    
    #4. Preferences for good BD
    
    if (beta< ((p_go-1)*(1/R)+d_prime*alpha*(gamma_go-gamma_bd*p_go)) ) {
      bool_i1 = 1
    }
    else {bool_i1 = 0}
    
    if (beta< ((p_gd-1)*(1/R)+d_prime*alpha*(gamma_gd-gamma_bd*p_gd)) ) {
      bool_i2 = 1
    }
    else {bool_i2 = 0}
    
    if ((1-alpha)*beta< ((p_bo-1)*(1/R)+d_prime*alpha*(gamma_bo-gamma_bd*p_bo)) ) {
      bool_i3 = 1
    }
    else {bool_i3 = 0}

#We fill our boolean matrices of exclusive consumption.

    mat_go_only[x,y]= bool_f1*bool_f2*bool_f3
    mat_gd_only[x,y]= bool_g1*bool_g2*bool_g3
    mat_bo_only[x,y]= bool_h1*bool_h2*bool_h3
    mat_bd_only[x,y]= bool_i1*bool_i2*bool_i3
    test[x,y]= mat_go_only[x,y]+mat_bo_only[x,y]+mat_gd_only[x,y]+mat_bd_only[x,y]        
    
  }
  
}
  
#Our exclusive consumption matrices are fine, they tie in with our graphs
#Summing the 4 matrices, we have a few zero coefficients, corresponding to zones in which several goods are consumed.



####### GRAPHICAL REPRESENTATIONS ##########

#To make better figures : https://r-graph-gallery.com/27-levelplot-with-lattice.html


#1. Exclusive consumption zones

#Good GO

mat_go_only<- mat_go_only %>% 
  as.tibble() %>% 
  as.list() %>% 
  lapply(function(x) rev(x)) %>% 
  do.call(rbind, .) %>% 
  as.matrix()

colnames(mat_go_only)<-alphas
row.names(mat_go_only)<-betas

#convert to binary scale on the plot (currently not working)
#mat1 <- raster(mat1)
#mat1 <- ratify(mat1)
#rat <- levels(mat1)[[1]]
#rat$code <- c(0,1)
#levels(mat1) <- rat

levelplot(mat_go_only, main="GO exclusive consumption zone (in pink)", xlab="alpha", ylab="beta", col.regions = cm.colors(121))   

#GD
#For now I have to launch this chunk of code several times to make it rotate so that it appears at the right place

mat_gd_only<- mat_gd_only %>% 
  as.tibble() %>% 
  as.list() %>% 
  lapply(function(x) rev(x)) %>% 
  do.call(rbind, .) %>% 
  as.matrix()

colnames(mat_gd_only)<-alphas
row.names(mat_gd_only)<-betas

levelplot(mat_gd_only, main="GD exclusive consumption zone (in dark blue)", xlab="alpha", ylab="beta")    



#BO

mat_bo_only<- mat_bo_only %>% 
  as.tibble() %>% 
  as.list() %>% 
  lapply(function(x) rev(x)) %>% 
  do.call(rbind, .) %>% 
  as.matrix()

colnames(mat_bo_only)<-alphas
row.names(mat_bo_only)<-betas


levelplot(mat_bo_only, main="BO exclusive consumption zone (in pink)", xlab="alpha", ylab="beta",col.regions = cm.colors(121))   


#BD

mat_bd_only<- mat_bd_only %>% 
  as.tibble() %>% 
  as.list() %>% 
  lapply(function(x) rev(x)) %>% 
  do.call(rbind, .) %>% 
  as.matrix()

colnames(mat_bd_only)<-alphas
row.names(mat_bd_only)<-betas

levelplot(mat_bd_only, main="BD exclusive consumption zone (in dark blue)", xlab="alpha", ylab="beta")  



######COMPUTATION OF AGGREGATED QUANTITIES########

#In each case, we first compute the consumption shares and the quantities consumed.
#And then we derive the subsequent total environmental impacts.

### 1) Exclusive consumption 
 
###### 1. Uniformly distributed population ######

#Consumption shares
pop_go_unif = pop*mat_go_only
pop_gd_unif = pop*mat_gd_only
pop_bo_unif = pop*mat_bo_only
pop_bd_unif = pop*mat_bd_only

#We express these consumption shares in percentages
excl_cons_go_scenar1 = sum(pop_go_unif)*100
excl_cons_gd_scenar1 = sum(pop_gd_unif)*100
excl_cons_bo_scenar1 = sum(pop_bo_unif)*100
excl_cons_bd_scenar1 = sum(pop_bd_unif)*100

excl_cons_go_scenar1 
excl_cons_gd_scenar1
excl_cons_bo_scenar1 
excl_cons_bd_scenar1

#We define the quantities consumed in each exclusive zone by each consumer.
indiv_quantity_go = R/p_go
indiv_quantity_gd = R/p_gd
indiv_quantity_bo = R/p_bo
indiv_quantity_bd = R

#Multiplying the population matrix coefficients by these quantities yields total quantities for each good : 
total_quantity_go = sum(indiv_quantity_go*mat_go_only)
total_quantity_gd = sum(indiv_quantity_gd*mat_gd_only) 
total_quantity_bo = sum(indiv_quantity_bo*mat_bo_only)
total_quantity_bd = sum(indiv_quantity_bd*mat_bd_only) 

#Finally, weighting these total quantities by the coefficients of environmental impact yields us a total environmental impact value of each good:
total_impact_go = gamma_go*total_quantity_go
total_impact_gd = gamma_gd*total_quantity_gd
total_impact_bo = gamma_bo*total_quantity_bo
total_impact_bd = gamma_bd*total_quantity_bd
total_impact_consumption = total_impact_go + total_impact_gd + total_impact_bo + total_impact_bd


#To study the different cases within this first scenario :
#we go to the 'parametrization' section above and directly change the parameters according to what we want to test, and then rerun the code up to this point.













###### 2. Sc?nario d'une concentration de la population autour de certaines valeurs (produit de lois betas ind?pendantes) ######

#On commence par concentration au milieu

pop_vo_scenar2 = pop2*mat_vo
pop_bo_scenar2 = pop2*mat_bo
pop_vd_scenar2 = pop2*mat_vd
pop_bd_scenar2 = pop2*mat_bd

conso_vo_scenar2 = sum(pop_vo_scenar2)*100   #ne marche plus quand les matrices sont grandes (alors que tout le reste si)...
conso_bo_scenar2 = sum(pop_bo_scenar2)*100
conso_vd_scenar2 = sum(pop_vd_scenar2)*100
conso_bd_scenar2 = sum(pop_bd_scenar2)*100

conso_vo_scenar2
conso_bo_scenar2
conso_vd_scenar2
conso_bd_scenar2

#Dans certains cas, ces changements de distribution ne changent pas la hi?rarchie entre les biens, mais uniquement les parts de consommation ; dans d'autres l'ordre est boulevers?
#On garde un d?coupage 41*41 pour les r?sultats num?riques (probl?me ? r?soudre)

#Cas de r?f?rence : l'ordre n'est pas modifi?, les parts deviennent maintenant dans l'ordre croissant  :
#0,05% pour VO (en baisse // au sc?nario uniforme) ;
#8,8% pour BO (en forte baisse) 
#22,3% pour BD (en baisse)
#68,9% pour VD (en hausse)
#On remarque que dans ce sc?nario seul VD profite de la concentration des consommateurs autour de la moyenne
#C'est le contraire pour les autres biens, surtout BO dont la zone est celle qui recoupe le plus des valeurs extr?mes dans ce premier cas (et dans ce sc?nario il y a peu de monde pour ces valeurs)


##### 2b. Concentration au milieu pour les sensibilit?s environnementales seulement (beta avec a=b=2 toujours) mais distribution uniforme des pr?f?rences pour le statut (c=d=1)


##### 2bis et ter. Concentrations sur des valeurs inf?rieures/sup?rieures de sensibililt? environnementale (et toujours uniforme sur le statut)

### 2B : Soci?t? moins pr?occup?e par l'environnement


### 2C : Soci?t? plus pr?occup?e par l'environnement


### 2D: Soci?t? moins pr?occup?e par le statut


### 2E : Economie de la Reine Rouge




