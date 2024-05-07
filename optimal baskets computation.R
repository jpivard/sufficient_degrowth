########## A model of consumption choices and sufficiency with environmental and social image preferences - Global optima computations ###########


#In this file, we aim to compute the optimal consumption baskets of consumers depending on their preference parameters alpha and beta (and the exogenous parameters)
#These baskets can be composed of one to four goods, so that we take into account mixed consumption, unlike in the first simulations
#To do so, we use a loop generation approach, making it possible to determine for each couple of coordinates (alpha,beta) the composition of the basket maximizing the associated utility

rm(list=ls())

###Packages

#install.packages("rgl")
#install.packages("reshape")
#install.packages("ggplot2") 
#install.packages("plotly")  

library(graphics)
library(rgl)
library("reshape")
library("ggplot2")
library(tibble)
library("plotly")     





######We first calibrate the parameters 

R = 20
p_go = 3
p_gd = 2
p_bo = 1.99
p_bd = 1
p <- c(p_go,p_gd,p_bo,p_bd)

d_prime = 0.01
gammaGO <- 1.6
gammaGD <- 0.8
gammaBD <- 2
gammaBO <- 4
gamma <- c(gammaGO, gammaGD, gammaBO, gammaBD)

theta = 1000  #How to calibrate that?

  
##### The size of the matrices (i.e. the accuracy of the computation) must also be chosen here.
step = 0.1
size = 1/step + 1



#Let us define an optimization function that loops over all possible baskets for each pair of parameters and chooses for each the one maximizing utility.

opti_function <- function(alpha, beta, R, p, gamma, d_prime, theta) {
  X1 <- seq(0, R/p[1], by = step)
  X2 <- seq(0, R/p[2], by = step)
  X3 <- seq(0, R/p[3], by = step)
  X4 <- seq(0, R/p[4], by = step)
  
  Umax <- -10000
  xs <- numeric(4)
  
  for (i in seq_along(X1)) {
    x1 <- X1[i]
    for (j in seq_along(X2)) {
      x2 <- X2[j]
      if (p[1]*x1 + p[2]*x2 <= R) {
        for (k in seq_along(X3)) {
          x3 <- X3[k]
          if (p[1]*x1 + p[2]*x2 + p[3]*x3 <= R) {
            x4 <- (R - p[1]*x1 - p[2]*x2 - p[3]*x3)/p[4]
            u_cons <- log(x1 + x2 + x3 + x4 + .Machine$double.eps)
            u_damage <- -alpha * (theta + d_prime * (gamma[1]*x1 + gamma[2]*x2 + gamma[3]*x3 + gamma[4]*x4))
            v_status <- beta * (alpha * (x1 + x2) + (1 - alpha) * (x1 + x3))
            U <- u_cons +  u_damage +  v_status
            if (U >= Umax) {
              Umax <- U
              xs[1] <- x1
              xs[2] <- x2
              xs[3] <- x3
              xs[4] <- x4
            }
          }
        }
      }
    }
  }
  
  xopt <- xs
  yopt <- Umax
  
  return(list(xopt = xopt, yopt = yopt))
}


#Test for some specific values of preference parameters


#alpha=beta=0
opti_function(0,0,R,p,gamma,d_prime,theta)
# Seems OK : such a consumer only consumes the BD good, and we can find the corresponding max utility

#alpha=beta=1
opti_function(1,1,R,p,gamma,d_prime,theta)
#It is coherent too : only consumes GD

#What happens in the middle of the alpha distribution and bottom of the beta distrib? There should be mixed consumption according to the previous method
opti_function(1/2,1/4,R,p,gamma,d_prime,theta)
#Coherent with the math again : the point (0.5,0.25) is in the GO+GD mixed zone

#CCL : the function seems to work fine.



#Let us now loop over preference parameters and represent optimal consumption baskets for each couple of coordinates

a <- seq(0, 1, by = step)
b <- seq(0, 1, by = step)
rev_a <-  seq(1, 0, by = -step)
rev_b <-  seq(1, 0, by = -step)
A <- matrix(rep(a, length(b)), nrow = length(a), byrow = TRUE)
B <- matrix(rep(rev_b, length(a)), nrow = length(b), byrow = FALSE)

D <- matrix(0, nrow = nrow(A), ncol = ncol(A))
D1 <- matrix(0, nrow = nrow(A), ncol = ncol(A))
D2 <- matrix(0, nrow = nrow(A), ncol = ncol(A))
D3 <- matrix(0, nrow = nrow(A), ncol = ncol(A))
D4 <- matrix(0, nrow = nrow(A), ncol = ncol(A))
M <- matrix(0, nrow = nrow(A), ncol = ncol(A))

DD <- list()
for (k in 1:4) {
  DD[[k]] <- matrix(0, nrow = nrow(A), ncol = ncol(A))
}

for (i in 1:nrow(A)) {
  for (j in 1:ncol(A)) {
    alpha <- A[i, j]
    beta <- B[i, j]
    
    optibasket <- opti_function(alpha, beta, R, p, gamma, d_prime, theta)
    xopt <- optibasket$xopt
    yopt <- optibasket$yopt
    l <- which(xopt > 0)
    
    if (length(l) > 0) {
      D[i, j] <- sum(2^(l-1))
      M[i, j] <- yopt
      D1[i, j] <- xopt[1]
      D2[i, j] <- xopt[2]
      D3[i, j] <- xopt[3]
      D4[i, j] <- xopt[4]
      for (k in 1:length(l)) {   #To get the number of goods consumed depending on location in the plane
        DD[[l[k]]][i, j] <- 1
      }
    }
  }
}

#We check the matrices yield the quantities consumed of each good depending on the values of (alpha,beta) i.e. on the position in the graph
print(D1)
print(D2)
print(D3)
print(D4)
#Seems fine (at least for small dimensions)
#Except maybe BD being slightly consumed in the BO exclusive zone, why ?


#Computing the size of global solutions = Compute the sum of DD{1}, DD{2}, DD{3}, and DD{4} in order to get the number of goods consumed in each point
S <- DD[[1]] + DD[[2]] + DD[[3]] + DD[[4]]
S


#Computing the budget shares spent in each good (in percentage)

shareGO<- matrix(0, nrow = nrow(A), ncol = ncol(A))
shareGD <- matrix(0, nrow = nrow(A), ncol = ncol(A))
shareBO <- matrix(0, nrow = nrow(A), ncol = ncol(A))
shareBD <- matrix(0, nrow = nrow(A), ncol = ncol(A))

for (i in 1:nrow(shareGO)) {
  for (j in 1:ncol(shareGO)) {
   shareGO[i,j]<- (p[1]*D1[i,j]/R)*100
   shareGD[i,j]<- (p[2]*D2[i,j]/R)*100
   shareBO[i,j]<- (p[3]*D3[i,j]/R)*100
   shareBD[i,j]<- (p[4]*D4[i,j]/R)*100
  }
}


###### Graphical representations of the matrices : TO BE CONTINUED########

#1. Number of goods consumed 

#Transformation to be applied to each matrix so that the display is correct







# Define new row and column names
new_row_names <- rev_b
new_col_names <- a

# Rename rows and columns
rownames(S_plot) <- new_row_names
colnames(S_plot) <- new_col_names

# Print the matrix with renamed rows and columns
print(S_plot)




# Create the heatmap
heatmap(S, scale = "none", Rowv = NA, Colv = NA,
        col = c("red", "yellow", "blue"), 
        main = "Number of goods consumed - Reference case",
        cexRow = 0.7, cexCol = 0.7)
legend("right", legend = c("1", "2", "3"), fill = c("red", "yellow", "blue"))


plot_ly(z = S_plot, type = "heatmap")  




#3D representations of the zones where good k is consumed 
for (k in 1:4) {
  # Open a new 3D plot window
  open3d()
  # Create a 3D surface plot
  persp3d(A, B, DD[[k]], col = "lightblue", xlab = "alpha", ylab = "beta", zlab = "DD", main = paste("Domain of global solutions with good", k, ""))
}
#Zones are correct, but graphs should be improved.



#Legend : all possible baskets appearing in the figures
Z <- list(c(1, 2, 3, 4), c(1, 2, 3), c(1, 2, 4), c(1, 3, 4), c(2, 3, 4), c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4), 1, 2, 3, 4)
num <- numeric(length(Z))
for (k in seq_along(Z)) {
  if (is.list(Z[[k]])) {
    num[k] <- sum(2^(Z[[k]] - 1))
    cat(num[k], " ", Z[[k]], "\n")
  } else {
    num[k] <- 2^(Z[[k]] - 1)
    cat(num[k], " ", Z[[k]], "\n")
  }
}


#Plotting the different baskets

# Create a title string using sprintf
str <- sprintf("global solutions")

# Set the plot title and the labels
title(str)
xlab <- "alpha"
ylab <- "beta"

# Open a new plot window
plot.new()

# Create an image plot
image(a, b, D, col = heat.colors(100), xlab = xlab, ylab = ylab)

# Reverse the y-axis direction
axis(2, at = seq_along(b), labels = b)












