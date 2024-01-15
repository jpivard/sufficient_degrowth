########## A model of consumption choices and sufficiency with environmental and social image preferences - First simulation method ###########

# Install ggplot2 package if not already installed
#install.packages("ggplot2")


library(lattice)
library(tidyverse)
library (gtools)
library(RColorBrewer)
library(ggplot2)


######## INTRODUCTION ###############


rm(list=ls())

#This script is designed for numerical simulations associated to the model presented here : https://fr.overleaf.com/project/653a7f038e83eb35864f2593 
# It aims to calculate the consumption areas and quantities of different goods. 
#These goods are defined along two dimensions: green/brown, conspicuous/discreet.
#Depending on their personal preferences, agents can choose to consume one or several goods.




#Method used in this script : As our heterogeneous consumers are placed in a [0,1]*[0,1] grid according to their preferences (on the x-axis, environmental awareness; on the y-axis, concern for social image), 
#The idea behind this numerical method is to discretize/rasterize the square into several points/pixels (expressing the two types of preferences in percentages if the number of points is 101 on the square ; for example individuals lying at the origin being indifferent to both, while those at the top right hand corner have maximum sensitivities to both).


#Matrix size is defined as a function of the step between 2 pixels.
step = 0.025
size = 1/step + 1
size


#Let's also define our preference parameters, over which we will iterate
alphas <- seq(from = 0,to = 1,by=step)
betas_rev <- seq(from = 1,to = 0,by=-step) #reversed so that coordinates appear in the right order
betas <- seq(from = 0,to = 1,by=step) #only for graphs


######### SCENARIO CHOICE - Assign values to shape parameters according to the scenario ############

#Basic scenario 2 (2A): two beta laws with the two same shape parameters (bell curves) to model the concentration in the middle of the distribution
a=b=c=d=2

#Variant (2b): the shape parameters are changed (the beta parameter is now assumed to follow a uniform distribution, i.e. a beta distribution with c=d=1).
#a=b=2
#c=d=1


#Scenario 2B: Lower/laxer environmental norm (concentration on the left, beta law with shape parameters 5 and 1 on alpha, uniform on beta).
#a=5
#b=c=d=1

#Scenario 2C : Higher/stringent environmental norm (right-hand concentration)
#a=c=d=1
#b=5



#Scenario 2D : social image matters less (widespread discretion)
#a=b=c=1
#d=5

#Scenario 2E: social image matters more ("Red Queen Economy", Keep Up With the Joneses...)
#a=b=d=1
#c=5

######## POPULATION DISTRIBUTIONS (Part 1) ##############

#Here, square matrices of dimension 'size' are defined. 


#To start with (scenario 1), we generate a uniformly distributed population density matrix of this size, giving the proportions of individuals positioned at each point on the map. 
n1 <- size
n2 <- size
pop <- matrix (rep(1/(size*size), n1*n2), n1, n2)
pop
#For example, by choosing a step size of 0.2, we have a matrix of size 6*6, with 36 coefficients ; while it is 11*11 (and 101 coeffs) for a 0.1 step (the default choice in the first place)



#Scenario 2: population concentration around some values
#This case is modeled by a product of independent beta laws for each preference parameter.

pop2 <- matrix(,n1,n2) #an empty matrix
 

###Code that computes discrete densities for scenario 2###

#We create our quantile vectors (in line)

alphas_vector <- t(matrix(alphas))
betas_vector <- t(matrix(betas))

#Starting from the repartition function, we define the discrete density of a beta (one for each of the 2 dimensions).

densite_beta_dim1 <- function(n, a, b, step) {
  if (n == 0) {
    return(pbeta(1, a, b) - pbeta(1 - step/2, a, b))
  } 
  if (n == 1/step) {
    return(pbeta(step/2,a,b) - pbeta(0,a,b))
  } 
  else {
    densite <- pbeta(step/2 + 1 - n * step, a, b) - pbeta(step/2 + 1 - (n + 1) * step, a, b)
    return(densite)
  }
}

densite_beta_dim2 <- function(n, c, d, step) {
  if (n == 0) {
    return(pbeta(1, c, d) - pbeta(1 - step/2, c, d))
  } 
  if (n == 1/step) {
    return(pbeta(step/2,c,d) - pbeta(0,c,d))
  } 
  else {
    densite <- pbeta(step/2 + 1 - n * step, c, d) - pbeta(step/2 + 1 - (n + 1) * step, c, d)
    return(densite)
  }
}



####################################### POPULATION DISTRIBUTIONS (Part 2) ###########




# Calling the functions at n = 0 and storing the result
result_at_0_dim1 <- densite_beta_dim1(0, a, b, step)
result_at_0_dim2 <- densite_beta_dim2(0, c, d, step)


# Storing other values iteratively in a matrix 
num_iterations <- size - 1  # Change this value as needed
proba_values_dim1 <- matrix(0, nrow = 1, ncol = num_iterations + 1)
proba_values_dim2 <- matrix(0, nrow = 1, ncol = num_iterations + 1)

for (i in 1:num_iterations) {
  proba_values_dim1[1, i + 1] <- densite_beta_dim1(i, a, b, step)
  proba_values_dim2[1, i + 1] <- densite_beta_dim2(i, c, d, step)
}

proba_values_dim1[1]=result_at_0_dim1
proba_values_dim2[1]=result_at_0_dim2

sum(proba_values_dim1)
sum(proba_values_dim2)
#bell-shaped densities, symmetrical and sum to one : perfect !!!


#We multiply our two probability vectors to obtain the probabilities for the entire population.
pop2 = t(proba_values_dim2)%*%proba_values_dim1
#There is indeed a higher population concentration on median values compared to the uniform case.

sum(pop2) 
#just checking, we're good in all dimensions :)
        


#Remark : Repartition functions can be represented (the same ones, with identical shape parameters)
fdr_beta_dim1 <-pbeta(alphas, shape1 = a, shape2 =b)
plot(fdr_beta_dim1,alphas)

fdr_beta_dim2 <-pbeta(betas, shape1 = c, shape2 =d)
plot(fdr_beta_dim2,betas)


#Astuce pour pas avoir à tout relancer = faire une fonction qui demande les paramètres et renvoie tous les résultats
programme<- function() {
  
  x <- as.numeric(readline(prompt = "valeur de x ?"))
  y<- as.numeric(readline(prompt = "valeur de y ?"))
  z<- as.numeric(readline(prompt = "valeur de z ?"))
  
  fonction <- function(x,y,z){
    return(x*y + z*x)
  }
  print(fonction(x, y, z))
}



############# BORDER PARAMETRIZATION ####################

#We define our "borders", delimiting the consumption zones according to the exogenous model parameters (so that we can study the different cases by changing only the following values).
#we initialize the parameters of our "borders", i.e. the prices of the four goods, individual income and marginal damage as well as impact parameters (satisfying model assumptions)
p_go = 3
p_gd = 2
p_bo = 1.99
p_bd = 1 
R = 20
d_prime = 0.01
gamma_bo = 4
gamma_bd = 2
gamma_go = 1.6
gamma_gd = 0.8





######## COMPUTATION OF CONSUMPTION ZONE MATRICES  ###############


#1st method (inspired from the simpler model):
#We run a nested loop on our individuals (i.e. their alpha and beta parameters, which we will convert - by normalizing - into (x,y) coordinates where the coordinates lie between 1 and size).
#For each point (i.e. pair of parameters) we look at which of the 12 constraints are satisfied, which determines the type of good consumed at that point (and the corresponding matrix among the 4 defined above, which will be filled by a "1" at the associated (x,y) point where the good is exclusively consumed).
#The idea is to create four size*size boolean matrices (true = 1 /false = 0), one for each type of consumer; and to examine for each point on the square (using nested loops) whether it satisfies the constraints corresponding to a certain type.

#We start by exclusive consumption zones only (defining the 4 corresponding matrices)
#Then we compute the mixed consumption zones matrices 

#In the end, we'll check that the sum of these four matrices gives 1s in the exclusive consumption zones and integers from 2 to 4 elsewhere.


#1.Exclusive consumption zones


# Function to calculate boolean value based on condition
calculate_boolean <- function(condition) {
  if (condition) {
    return(1)
  } else {
    return(0)
  }
}

# Function to calculate preferences for each good
calculate_preferences <- function(alpha, beta, p_go, p_gd, p_bo, p_bd, R, d_prime, 
                                  gamma_go, gamma_gd, gamma_bo, gamma_bd) {
  bool_f1 <- calculate_boolean((1 - (p_go / p_gd) * alpha) * beta > (((p_go / p_gd) - 1) * (p_go / R) + d_prime * alpha * (gamma_go - gamma_gd * (p_go / p_gd))))
  bool_f2 <- calculate_boolean((1 - (p_go / p_bo) * (1 - alpha)) * beta > (((p_go / p_bo) - 1) * (p_go / R) + d_prime * alpha * (gamma_go - gamma_bo * (p_go / p_bo))))
  bool_f3 <- calculate_boolean(beta > ((p_go/p_bd - 1) * (p_go / R) + d_prime * alpha * (gamma_go - gamma_bd * p_go/p_bd)))
  
  bool_g1 <- calculate_boolean((alpha - (p_gd / p_go)) * beta > (((p_gd / p_go) - 1) * (p_gd / R) + d_prime * alpha * (gamma_gd - gamma_go * (p_gd / p_go))))
  bool_g2 <- calculate_boolean((alpha - (p_gd / p_bo) * (1 - alpha)) * beta > (((p_gd / p_bo) - 1) * (p_gd / R) + d_prime * alpha * (gamma_gd - gamma_bo * (p_gd / p_bo))))
  bool_g3 <- calculate_boolean(alpha * beta > ((p_gd/p_bd - 1) * (p_gd / R) + d_prime * alpha * (gamma_gd - gamma_bd * p_gd/p_bd)))
  
  bool_h1 <- calculate_boolean((1 - alpha - (p_bo / p_go)) * beta > (((p_bo / p_go) - 1) * (p_bo / R) + d_prime * alpha * (gamma_bo - gamma_go * (p_bo / p_go))))
  bool_h2 <- calculate_boolean((1 - alpha * (p_bo / p_gd)) * beta > (((p_bo / p_gd) - 1) * (p_bo / R) + d_prime * alpha * (gamma_bo - gamma_gd * (p_gd / p_bo))))
  bool_h3 <- calculate_boolean((1 - alpha) * beta > ((p_bo/p_bd - 1) * (p_bo / R) + d_prime * alpha * (gamma_bo - gamma_bd * p_bo/p_bd)))
  
  bool_i1 <- calculate_boolean((p_bd/p_go)*beta < ((1- p_bd/p_go) * (p_bd / R) - d_prime * alpha * (gamma_bd - gamma_go * p_bd/p_go)))
  bool_i2 <- calculate_boolean(alpha*(p_bd/p_gd)*beta < ((1- p_bd/p_gd) * (p_bd / R) - d_prime * alpha * (gamma_bd - gamma_gd * p_bd/p_bd)))
  bool_i3 <- calculate_boolean((1-alpha)*(p_bd/p_bo)*beta < ((1- p_bd/p_bo) * (p_bd / R) - d_prime * alpha * (gamma_bd - gamma_bo * p_bd/p_bo)))
  return(c(bool_f1, bool_f2, bool_f3, bool_g1, bool_g2, bool_g3, bool_h1, bool_h2, bool_h3, bool_i1, bool_i2, bool_i3))
}

# Initialize matrices
mat_go_only <- matrix(, nrow = size, ncol = size)  
mat_gd_only <- matrix(, nrow = size, ncol = size)
mat_bo_only <- matrix(, nrow = size, ncol = size)
mat_bd_only <- matrix(, nrow = size, ncol = size)
test <- matrix(, nrow = size, ncol = size)

# Nested loop to qualify each constraint everywhere in the map
for (beta in betas_rev)  {
  for (alpha in alphas) {
    x = size - beta / step
    y = alpha / step + 1
    
    # Calculate preferences for each good
    preferences <- calculate_preferences(alpha, beta, p_go, p_gd, p_bo, p_bd, R, d_prime, gamma_go, gamma_gd, gamma_bo, gamma_bd)
    
    # Fill boolean matrices of exclusive consumption
    mat_go_only[x, y] <- prod(as.numeric(preferences[1:3]))
    mat_gd_only[x, y] <- prod(as.numeric(preferences[4:6]))
    mat_bo_only[x, y] <- prod(as.numeric(preferences[7:9]))
    mat_bd_only[x, y] <- prod(as.numeric(preferences[10:12]))
    
    test[x, y] <- mat_go_only[x, y] + mat_bo_only[x, y] + mat_gd_only[x, y] + mat_bd_only[x, y]
  }
}

#Our exclusive consumption matrices are almost fine, they perfectly tie in with our graphs.
#Summing the 4 matrices, we have a few zero coefficients, corresponding to zones in which several goods are consumed.


####### GRAPHICAL REPRESENTATIONS ##########

#To make better figures : https://r-graph-gallery.com/27-levelplot-with-lattice.html


#1. Exclusive consumption zones

#Good GO

  mat1<- mat_go_only %>% 
  as.tibble() %>% 
  as.list() %>% 
  lapply(function(x) rev(x)) %>% 
  do.call(rbind, .) %>% 
  as.matrix()

colnames(mat1)<-alphas
row.names(mat1)<-betas

levelplot(mat1, main="GO exclusive consumption zone (in pink)", xlab="alpha", ylab="beta", col.regions = cm.colors(121))   

#GD
mat2<- mat_gd_only %>% 
  as.tibble() %>% 
  as.list() %>% 
  lapply(function(x) rev(x)) %>% 
  do.call(rbind, .) %>% 
  as.matrix()

colnames(mat2)<-alphas
row.names(mat2)<-betas

levelplot(mat2, main="GD exclusive consumption zone (in blue)", xlab="alpha", ylab="beta")    


#BO

mat3<- mat_bo_only %>% 
  as.tibble() %>% 
  as.list() %>% 
  lapply(function(x) rev(x)) %>% 
  do.call(rbind, .) %>% 
  as.matrix()

colnames(mat3)<-alphas
row.names(mat3)<-betas


levelplot(mat3, main="BO exclusive consumption zone (in pink)", xlab="alpha", ylab="beta",col.regions = cm.colors(121))   


#BD

mat4<- mat_bd_only %>% 
  as.tibble() %>% 
  as.list() %>% 
  lapply(function(x) rev(x)) %>% 
  do.call(rbind, .) %>% 
  as.matrix()

colnames(mat4)<-alphas
row.names(mat4)<-betas

levelplot(mat4, main="BD exclusive consumption zone (in blue)", xlab="alpha", ylab="beta")  


########### COMPUTATION OF AGGREGATED QUANTITIES IN DIFFERENT SCENARIOS########

#In each case, we first compute the consumption shares and the quantities consumed.
#And then we derive the subsequent total environmental impacts.

### 1) Exclusive consumption only
 
###### Scenario 1  : Uniformly distributed population ######

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
total_exclusive = excl_cons_go_scenar1 + excl_cons_gd_scenar1  + excl_cons_bo_scenar1 + excl_cons_bd_scenar1 

excl_cons_go_scenar1 
excl_cons_gd_scenar1
excl_cons_bo_scenar1 
excl_cons_bd_scenar1
total_exclusive



#We represent these shares in a diagram 

categories <- c("Good GO", "Good GD", "Good BO", "Good BD", "Mixed")
percentages <- c(excl_cons_go_scenar1, excl_cons_gd_scenar1, excl_cons_bo_scenar1, excl_cons_bd_scenar1, 100-total_exclusive)

# Create a data frame
data <- data.frame(Category = categories, Percentage = percentages)

# Set colors
custom_colors <- c("lightgrey", "brown", "lightgreen", "darkgreen", "blue") 

# Create a bar plot (caption to be changed according to the case tested)
plot <- ggplot(data, aes(x = Category, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity") +
  labs(title = "Distribution of consumers - Environmental tax, Scenario 1",
       x = "Categories",
       y = "Percentage") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Format y-axis as percentage
  theme_minimal()+
  scale_fill_manual(values = custom_colors)

print(plot)




###### Scenario 2 : A concentration of the population around certain values (product of independent Beta laws) ######

pop_go_scenar2 = pop2*mat_go_only
pop_gd_scenar2 = pop2*mat_gd_only
pop_bo_scenar2 = pop2*mat_bo_only
pop_bd_scenar2= pop2*mat_bd_only

excl_conso_go_scenar2 = sum(pop_go_scenar2)*100   
excl_conso_bo_scenar2 = sum(pop_bo_scenar2)*100
excl_conso_gd_scenar2 = sum(pop_gd_scenar2)*100
excl_conso_bd_scenar2 = sum(pop_bd_scenar2)*100
total_exclusive2 = excl_conso_go_scenar2 + excl_conso_gd_scenar2  + excl_conso_bo_scenar2 + excl_conso_bd_scenar2

excl_conso_go_scenar2
excl_conso_bo_scenar2
excl_conso_gd_scenar2
excl_conso_bd_scenar2
total_exclusive2

#We represent these shares in a diagram 

categories2 <- c("Good GO", "Good GD", "Good BO", "Good BD", "Mixed")
percentages2 <- c(excl_conso_go_scenar2, excl_conso_gd_scenar2, excl_conso_bo_scenar2, excl_conso_bd_scenar2, 100-total_exclusive2)  

# Create a data frame
data2 <- data.frame(Category = categories2, Percentage = percentages2)

custom_colors <- c("lightgrey", "brown", "lightgreen", "darkgreen", "blue")  


# Create a bar plot (caption to be changed according to the case/scenario tested)
plot2 <- ggplot(data2, aes(x = Category, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity") +
  labs(title = "Distribution of consumers - 20% drop in income, Scenario 2E",
       x = "Categories",
       y = "Percentage") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Format y-axis as percentage
  theme_minimal()+
  scale_fill_manual(values = custom_colors)

print(plot2)


### Common to all scenarios  : 

###We define the quantities consumed in each exclusive zone by each consumer.

indiv_quantity_go = R/p_go
indiv_quantity_gd = R/p_gd
indiv_quantity_bo = R/p_bo
indiv_quantity_bd = R/p_bd


### Multiplying the population matrix coefficients by these quantities yields total quantities for each good : 

#Scenario 1 : for a uniform distribution, it is sufficient to sum all the coefficients wihtout using the pop matrix as we know the population is equally distributed

total_quantity_go = sum(indiv_quantity_go*mat_go_only)
total_quantity_gd = sum(indiv_quantity_gd*mat_gd_only) 
total_quantity_bo = sum(indiv_quantity_bo*mat_bo_only)
total_quantity_bd = sum(indiv_quantity_bd*mat_bd_only) 
total_quantity_consumed = total_quantity_go + total_quantity_gd + total_quantity_bo + total_quantity_bd

total_quantity_go
total_quantity_gd
total_quantity_bo
total_quantity_bd
total_quantity_consumed



#Scenario 2 

#We first rescale population matrices by the size of the population (the square of the size variable) so that a bit more than 1 individual lie in each pixel in a concentrated area (max seems to be around 2.3), whereas one can find values between 0 and 1 when one gets further from the concentrated area
pop_go_scenar2_rescaled = pop_go_scenar2*size^2
pop_gd_scenar2_rescaled = pop_gd_scenar2*size^2
pop_bo_scenar2_rescaled = pop_bo_scenar2*size^2
pop_bd_scenar2_rescaled= pop_bd_scenar2*size^2


#And then we apply the same principle

total_quantity_go_scenar2 = sum(indiv_quantity_go*pop_go_scenar2_rescaled)
total_quantity_gd_scenar2 = sum(indiv_quantity_gd*pop_gd_scenar2_rescaled) 
total_quantity_bo_scenar2 = sum(indiv_quantity_bo*pop_bo_scenar2_rescaled)
total_quantity_bd_scenar2 = sum(indiv_quantity_bd*pop_bd_scenar2_rescaled) 
total_quantity_consumed_scenar2 = total_quantity_go_scenar2 + total_quantity_gd_scenar2 + total_quantity_bo_scenar2 + total_quantity_bd_scenar2

total_quantity_go_scenar2
total_quantity_gd_scenar2
total_quantity_bo_scenar2
total_quantity_bd_scenar2
total_quantity_consumed_scenar2






###Finally, in each scenario, weighting these total quantities by the coefficients of environmental impact yields us a total environmental impact value of each good:

#Scenario 1

total_impact_go = gamma_go*total_quantity_go
total_impact_gd = gamma_gd*total_quantity_gd
total_impact_bo = gamma_bo*total_quantity_bo
total_impact_bd = gamma_bd*total_quantity_bd
total_impact_consumption = total_impact_go + total_impact_gd + total_impact_bo + total_impact_bd

total_impact_go
total_impact_gd
total_impact_bo
total_impact_bd
total_impact_consumption


#Plot the impacts 

data_impacts <- data.frame(
  category = c("Good BD", "Good  BO", "Good GD ", "Good  GO"),
  values = c(total_impact_bd, total_impact_bo, total_impact_gd, total_impact_go)
)

bar_colors <- c("lightgrey", "brown", "lightgreen", "darkgreen")

ylim_values <- c(0, 50000)

barplot(data_impacts$values, names.arg = data_impacts$category,  col = bar_colors, main = "Env. impacts by good type - Environmental tax, Sce. 1", ylab = "Impacts", border = "black",ylim = ylim_values)



#Scenario 2

total_impact_go_scenar2 = gamma_go*total_quantity_go_scenar2
total_impact_gd_scenar2 = gamma_gd*total_quantity_gd_scenar2
total_impact_bo_scenar2 = gamma_bo*total_quantity_bo_scenar2
total_impact_bd_scenar2 = gamma_bd*total_quantity_bd_scenar2
total_impact_consumption_scenar2 = total_impact_go_scenar2 + total_impact_gd_scenar2 + total_impact_bo_scenar2 + total_impact_bd_scenar2

total_impact_go_scenar2
total_impact_gd_scenar2
total_impact_bo_scenar2
total_impact_bd_scenar2
total_impact_consumption_scenar2

#Plot
data_impacts <- data.frame(
  category = c("Good BD", "Good  BO", "Good GD ", "Good  GO"),
  values = c(total_impact_bd_scenar2, total_impact_bo_scenar2, total_impact_gd_scenar2, total_impact_go_scenar2)
)

ylim_values <- c(0, 50000) #to be changed according to the analysis

bar_colors <- c("lightgrey", "brown", "lightgreen", "darkgreen")

barplot(data_impacts$values, names.arg = data_impacts$category,  col = bar_colors, main = "Env.impacts by good type - Reference case, Sce. 2B", ylab = "Impacts", border = "black",ylim = ylim_values)




###As an additional sufficiency metrics, let's compute the total impact per capita

#Scenario 1 : we simply divide aggregate impact by the population inside exclusive consumption zones (here we can simply divide total impact by the percentage of exclusive pixels as we have one pixel = one person)

total_impact_percapita = total_impact_consumption/((excl_cons_go_scenar1+excl_cons_gd_scenar1+excl_cons_bo_scenar1+excl_cons_bd_scenar1)*1/100*size^2)
total_impact_percapita

#Scenario 2 : To obtain the impact per capita, here we need to know how many exclusive consumers there are taking into account the rescaled matrices 

number_exclusives = sum(pop_go_scenar2_rescaled)+sum(pop_gd_scenar2_rescaled)+sum(pop_bo_scenar2_rescaled)+sum(pop_bd_scenar2_rescaled)
number_exclusives #different from size square

total_impact_percapita_scenar2 = total_impact_consumption_scenar2/number_exclusives
total_impact_percapita_scenar2

#Compute the variation (in percentage) of impacts per capita when changing pop distribution from uniform to concentration
variation_impactpercapita = ((total_impact_percapita_scenar2-total_impact_percapita)/total_impact_percapita)*100
variation_impactpercapita



#To study the different cases within each scenario :
#we go to the 'parametrization' section above and directly change the parameters according to what we want to test, and then rerun the code (corresponding to the chosen scenario) up to this point.


########### COMPARING CASES AND SCENARIOS #########

#To ensure better comparability, let's plot per capita variables (+variations?) in the different cases, in both scenarios.
#To do this, we have to store the results for each case in different variables, that we will finally use in the code for the plot.

#1. Per capita impacts

#Scenario 1 results : store in the suited variable according to case tested.
#total_impact_percapita_unif_ref <- total_impact_percapita 
#total_impact_percapita_unif_case2 <-  total_impact_percapita 
#total_impact_percapita_unif_case3 <-  total_impact_percapita 
#total_impact_percapita_unif_case4 <-  total_impact_percapita 
#total_impact_percapita_unif_case5 <-  total_impact_percapita 
#total_impact_percapita_unif_case6 <-  total_impact_percapita 
#total_impact_percapita_unif_case7 <-  total_impact_percapita 
#total_impact_percapita_unif_degrowth <- total_impact_percapita 
#total_impact_percapita_unif_tax <- total_impact_percapita

#Scenario 2A results
#total_impact_percapita_scenar2A_ref <- total_impact_percapita_scenar2 
#total_impact_percapita_scenar2A_case2 <-  total_impact_percapita_scenar2 
#total_impact_percapita_scenar2A_case3 <-  total_impact_percapita_scenar2 
#total_impact_percapita_scenar2A_case4 <-  total_impact_percapita_scenar2 
#total_impact_percapita_scenar2A_case5 <-  total_impact_percapita_scenar2 
#total_impact_percapita_scenar2A_case6 <-  total_impact_percapita_scenar2 
#total_impact_percapita_scenar2A_case7 <-  total_impact_percapita_scenar2 
#total_impact_percapita_scenar2A_degrowth <- total_impact_percapita_scenar2 




#Scenario 2B results
#total_impact_percapita_scenar2B_ref <- total_impact_percapita_scenar2 : OK
#total_impact_percapita_scenar2B_case2 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2B_case3 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2B_case4 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2B_case5 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2B_case6 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2B_case7 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2B_degrowth <- total_impact_percapita_scenar2 : OK 


#Scenario 2C results
#total_impact_percapita_scenar2C_ref <- total_impact_percapita_scenar2 : OK
#total_impact_percapita_scenar2C_case2 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2C_case3 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2C_case4 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2C_case5 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2C_case6 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2C_case7 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2C_degrowth <- total_impact_percapita_scenar2 



#Scenario 2D results
#total_impact_percapita_scenar2D_ref <- total_impact_percapita_scenar2
#total_impact_percapita_scenar2D_case2 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2D_case3 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2D_case4 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2D_case5 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2D_case6 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2D_case7 <-  total_impact_percapita_scenar2

#Scenario 2E results
#total_impact_percapita_scenar2E_ref <- total_impact_percapita_scenar2 : OK
#total_impact_percapita_scenar2E_case2 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2E_case3 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2E_case4 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2E_case5 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2E_case6 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2E_case7 <-  total_impact_percapita_scenar2
#total_impact_percapita_scenar2E_degrowth <- total_impact_percapita_scenar2 




#Let's start with case comparison in scenario 1.

# Store per capita impacts and case categories in lists or vectors
categories <- c("Ref", "Case 2", "Case 3", "Case 4", "Case 5", "Case 6","Degrowth","Tax")
per_capita_impacts <- c(
  total_impact_percapita_unif_ref,
  total_impact_percapita_unif_case2,
  total_impact_percapita_unif_case3,
  total_impact_percapita_unif_case4,
  total_impact_percapita_unif_case5,
  total_impact_percapita_unif_case6,
  total_impact_percapita_unif_degrowth,
  total_impact_percapita_unif_tax
  
)

# Create a data frame
data <- data.frame(Category = factor(categories, levels = c("Ref", "Tax", "Degrowth", "Case 2", "Case 3", "Case 4", "Case 5", "Case 6")), Numbers = per_capita_impacts)

# Set custom colors for the gradient fill
custom_gradient_colors <- colorRampPalette(c("green", "tan", "saddlebrown"))(length(per_capita_impacts))

# Create the bar plot using ggplot2
plot_pcimpact_unif <- ggplot(data, aes(x = Category, y = Numbers, fill = Numbers)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Per capita impacts in the different cases, Uniform distribution",
    x = "Cases",
    y = "Impacts"
  ) +
  theme_minimal() +
  scale_fill_gradientn(
    colors = custom_gradient_colors,
    guide = "legend",
    name = "Impacts"
  )

# Display the plot
print(plot_pcimpact_unif)



#We then do the same in scenario 2A.

# Store per capita impacts and case categories in lists or vectors
categories <- c("Ref", "Case 2", "Case 3", "Case 4", "Case 5", "Case 6","Degrowth")
per_capita_impacts <- c(
  total_impact_percapita_scenar2A_ref,
  total_impact_percapita_scenar2A_case2,
  total_impact_percapita_scenar2A_case3,
  total_impact_percapita_scenar2A_case4,
  total_impact_percapita_scenar2A_case5,
  total_impact_percapita_scenar2A_case6,
  total_impact_percapita_scenar2A_degrowth
)

# Create a data frame
data <- data.frame(Category = factor(categories, levels = c("Ref", "Degrowth", "Case 2", "Case 3", "Case 4", "Case 5", "Case 6")), Numbers = per_capita_impacts)

# Set custom colors for the gradient fill
custom_gradient_colors <- colorRampPalette(c("green", "tan", "saddlebrown"))(length(per_capita_impacts))

# Create the bar plot using ggplot2
plot_pcimpact_scen2 <- ggplot(data, aes(x = Category, y = Numbers, fill = Numbers)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Per capita impacts in the different cases, Scenario 2A",
    x = "Cases",
    y = "Impacts"
  ) +
  theme_minimal() +
  scale_fill_gradientn(
    colors = custom_gradient_colors,
    guide = "legend",
    name = "Impacts"
  )

# Display the plot
print(plot_pcimpact_scen2)


#Let's now plot the per capita impacts in the ref case according to the scenario

# Store per capita impacts and scenario categories in lists or vectors
categories <- c("Unif (Sce. 1)", "Sce. 2A", "Sce. 2B", "Sce. 2C", "Sce. 2D", "Sce. 2E")
per_capita_impacts <- c(
  total_impact_percapita_unif_ref,
  total_impact_percapita_scenar2A_ref,
  total_impact_percapita_scenar2B_ref,
  total_impact_percapita_scenar2C_ref,
  total_impact_percapita_scenar2D_ref,
  total_impact_percapita_scenar2E_ref
)

# Create a data frame
data <- data.frame(Category = factor(categories, levels = c("Unif (Sce. 1)", "Sce. 2A", "Sce. 2B", "Sce. 2C", "Sce. 2D", "Sce. 2E")), Numbers = per_capita_impacts)

# Set custom colors for the gradient fill
custom_gradient_colors <- colorRampPalette(c("green", "tan", "saddlebrown"))(length(per_capita_impacts))

# Create the bar plot using ggplot2
plot_pcimpact_ref <- ggplot(data, aes(x = Category, y = Numbers, fill = Numbers)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Per capita impacts in the different population scenarios, Ref. case",
    x = "Population distribution scenarios",
    y = "Impacts"
  ) +
  theme_minimal() +
  scale_fill_gradientn(
    colors = custom_gradient_colors,
    guide = "legend",
    name = "Impacts"
  )

# Display the plot
print(plot_pcimpact_ref)



#We now plot the different degrowth scenarios (depending on evolution of preferences)

# Store per capita impacts and scenario categories in lists or vectors
categories <- c("Ref-Center", "Degrowth-Center", "Degrowth-Left", "Degrowth-Right", "Degrowth-Top")
per_capita_impacts <- c(
  total_impact_percapita_scenar2A_ref,
  total_impact_percapita_scenar2A_degrowth,
  total_impact_percapita_scenar2B_degrowth,
  total_impact_percapita_scenar2C_degrowth,
  total_impact_percapita_scenar2E_degrowth
)

# Create a data frame
data <- data.frame(Category = factor(categories, levels = c("Ref-Center", "Degrowth-Center", "Degrowth-Left", "Degrowth-Right", "Degrowth-Top")), Numbers = per_capita_impacts)

# Set custom colors for the gradient fill
custom_gradient_colors <- colorRampPalette(c("green", "tan", "saddlebrown"))(length(per_capita_impacts))

# Create the bar plot using ggplot2
plot_pcimpact_degrowth <- ggplot(data, aes(x = Category, y = Numbers, fill = Numbers)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Per capita impacts in different degrowth scenarios relative to baseline",
    x = "Scenarios (case-type of population distribution concentration)",
    y = "Per capita impacts"
  ) +
  theme_minimal() +
  scale_fill_gradientn(
    colors = custom_gradient_colors,
    guide = "legend",
    name = "Impacts"
  )

# Display the plot
print(plot_pcimpact_degrowth)



