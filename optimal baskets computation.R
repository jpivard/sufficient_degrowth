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
  
  Umax <- -1000   #How to calibrate that?
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

alphas <- seq(0, 1, by = step)
betas <- seq(0, 1, by = step)
rev_alphas <-  seq(1, 0, by = -step)
rev_betas <-  seq(1, 0, by = -step)

a <- seq(0, 1, by = step)
b <- seq(0, 1, by = step)
rev_a <-  seq(1, 0, by = -step)
rev_b <-  seq(1, 0, by = -step)
A <- matrix(rep(a, length(b)), nrow = length(a), byrow = TRUE)
B <- matrix(rep(rev_b, length(a)), nrow = length(b), byrow = FALSE)



D <- matrix(0, nrow = size, ncol = size) 
D1 <- matrix(0, nrow = size, ncol = size) 
D2 <- matrix(0, nrow = size, ncol = size)
D3 <- matrix(0, nrow = size, ncol = size)
D4 <- matrix(0, nrow = size, ncol = size)
M <- matrix(0, nrow = size, ncol = size)



#D <- matrix(0, nrow = nrow(A), ncol = ncol(A))
#D1 <- matrix(0, nrow = nrow(A), ncol = ncol(A))
#D2 <- matrix(0, nrow = nrow(A), ncol = ncol(A))
#D3 <- matrix(0, nrow = nrow(A), ncol = ncol(A))
#D4 <- matrix(0, nrow = nrow(A), ncol = ncol(A))
#M <- matrix(0, nrow = nrow(A), ncol = ncol(A))


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
rownames(S) <- new_row_names
colnames(S) <- new_col_names

# Print the matrix with renamed rows and columns
#print(S)

# Reverse the rows of the matrix to flip the y-axis direction
S <- S[nrow(S):1, ]

# Create the heatmap
heatmap(S, scale = "none", Rowv = NA, Colv = NA,
        col = c("red", "yellow", "blue"), 
        main = "Number of goods consumed - Reference case",
        cexRow = 0.7, cexCol = 0.7)
legend("right", legend = c("1", "2", "3"), fill = c("red", "yellow", "blue"))




#2. Quantities consumed and Share of expenses in each good
#Only for good GO to start with

# Load required packages
#library(RColorBrewer)


#a) Quantities consumed

rownames(D1) <- rev_b
colnames(D1) <- a

## Trick to get the right legend

# Specify the number of decimal places for rounding
decimal_places <- 8  # Increase this value if necessary

# Specify a tolerance for considering values as equal
tolerance <- 1e-6  # Adjust this as needed

# Function to group values within a tolerance
group_values <- function(values, tolerance) {
  groups <- list()
  sorted_values <- sort(values)
  group <- c(sorted_values[1])
  
  for (i in seq(2, length(sorted_values))) {
    if (sorted_values[i] - group[length(group)] <= tolerance) {
      group <- c(group, sorted_values[i])
    } else {
      groups <- c(groups, list(group))
      group <- c(sorted_values[i])
    }
  }
  groups <- c(groups, list(group))
  
  # Extract unique group means as unique coefficients
  unique_coeffs <- sapply(groups, mean)
  
  return(unique_coeffs)
}

# Round the matrix values to the specified number of decimal places
D1_rounded <- round(D1, decimal_places)

# Get the unique coefficients using the grouping function with tolerance
unique_coeffs <- sort(group_values(D1_rounded, tolerance))

# Determine the number of unique values in the rounded matrix
num_unique_values <- length(unique_coeffs)

# Create a color palette from white to dark green based on the number of unique values
light_green_palette <- colorRampPalette(c("white", "darkgreen"))
color_palette <- light_green_palette(num_unique_values)

# Specify the desired size for the main title
cex_main <- 0.5  # Adjust this value as needed

# Set the global graphical parameter for main title size
par(cex.main = cex_main)

# Reverse the rows of the matrix to flip the y-axis direction
D1_reversed <- D1_rounded[nrow(D1_rounded):1, ]

# Plot the heatmap with the reversed matrix
heatmap(D1_reversed, scale = "none", Rowv = NA, Colv = NA,
        col = color_palette,
        main = "Quantities consumed of GO good - Reference case",
        cexRow = 0.7, cexCol = 0.7)

# Create a named vector associating each unique value in the matrix with a color
color_map <- setNames(color_palette, unique_coeffs)

# Adjust plot margins using `par()` function
# Adjust the right margin to allocate more space for the legend
par(mar = c(4, 4, 4, 4))

# Define the coordinates for the legend
# x and y values should range from 0 (bottom/left) to 1 (top/right)
# Try `x = 1.05` or higher to move the legend further right of the plot
legend_x <- 0.8  # Increase the value as needed to move the legend further right
legend_y <- 0.8  # Adjust this as needed

# Create the legend
legend(x = legend_x, y = legend_y,
       legend = unique_coeffs, fill = color_map,
       cex = 0.5, pt.cex = 1.5, # Adjust `cex` and `pt.cex` as needed
       y.intersp = 1.5) # Adjust spacing between legend items as needed






#b) Share of expenses

# Set the row and column names
rownames(shareGO) <- rev_b
colnames(shareGO) <- a

## Trick to get the right legend

# Specify the number of decimal places for rounding
decimal_places <- 8  # Increase this value if necessary

# Specify a tolerance for considering values as equal
tolerance <- 1e-6  # Adjust this as needed

# Function to group values within a tolerance
group_values <- function(values, tolerance) {
  groups <- list()
  sorted_values <- sort(values)
  group <- c(sorted_values[1])
  
  for (i in seq(2, length(sorted_values))) {
    if (sorted_values[i] - group[length(group)] <= tolerance) {
      group <- c(group, sorted_values[i])
    } else {
      groups <- c(groups, list(group))
      group <- c(sorted_values[i])
    }
  }
  groups <- c(groups, list(group))
  
  # Extract unique group means as unique coefficients
  unique_coeffs <- sapply(groups, mean)
  
  return(unique_coeffs)
}

# Round the matrix values to the specified number of decimal places
shareGO_rounded <- round(shareGO, decimal_places)

# Get the unique coefficients using the grouping function with tolerance
unique_coeffs <- sort(group_values(shareGO_rounded, tolerance))

# Determine the number of unique values in the rounded matrix
num_unique_values <- length(unique_coeffs)

# Create a color palette from white to dark green based on the number of unique values
light_green_palette <- colorRampPalette(c("white", "darkgreen"))
color_palette <- light_green_palette(num_unique_values)

# Specify the desired size for the main title
cex_main <- 0.5  # Adjust this value as needed

# Set the global graphical parameter for main title size
par(cex.main = cex_main)

# Reverse the rows of the matrix to flip the y-axis direction
shareGO_reversed <- shareGO_rounded[nrow(shareGO_rounded):1, ]

# Plot the heatmap with the reversed matrix
heatmap(shareGO_reversed, scale = "none", Rowv = NA, Colv = NA,
        col = color_palette,
        main = "Share of budget spent in GO - Reference case",
        cexRow = 0.7, cexCol = 0.7)

# Create a named vector associating each unique value in the matrix with a color
color_map <- setNames(color_palette, unique_coeffs)

# Adjust plot margins using par() function
# Adjust the right margin to allocate more space for the legend
par(mar = c(4, 4, 4, 4))

# Create the legend
# Define the coordinates for the legend
# x and y values should range from 0 (bottom/left) to 1 (top/right)
# To move the legend more to the right, increase the x value
# Try x = 1.05 or higher to move the legend further right of the plot
legend_x <- 0.8  # Increase the value as needed to move the legend further right
legend_y <- 0.8  # You can adjust this as needed

# Create the legend
legend(x = legend_x, y = legend_y,
       legend = unique_coeffs, fill = color_map,
       cex = 0.5, pt.cex = 1.5, # Adjust cex and pt.cex as needed
       y.intersp = 1.5) # Adjust spacing between legend items as needed




#Same graph for good BO 

# Set the row and column names
rownames(shareBO) <- rev_b
colnames(shareBO) <- a

## Trick to get the right legend

# Specify the number of decimal places for rounding
decimal_places <- 8  # Increase this value if necessary

# Specify a tolerance for considering values as equal
tolerance <- 1e-6  # Adjust this as needed

# Function to group values within a tolerance
group_values <- function(values, tolerance) {
  groups <- list()
  sorted_values <- sort(values)
  group <- c(sorted_values[1])
  
  for (i in seq(2, length(sorted_values))) {
    if (sorted_values[i] - group[length(group)] <= tolerance) {
      group <- c(group, sorted_values[i])
    } else {
      groups <- c(groups, list(group))
      group <- c(sorted_values[i])
    }
  }
  groups <- c(groups, list(group))
  
  # Extract unique group means as unique coefficients
  unique_coeffs <- sapply(groups, mean)
  
  return(unique_coeffs)
}

# Round the matrix values to the specified number of decimal places
shareBO_rounded <- round(shareBO, decimal_places)

# Get the unique coefficients using the grouping function with tolerance
unique_coeffs <- sort(group_values(shareBO_rounded, tolerance))

# Determine the number of unique values in the rounded matrix
num_unique_values <- length(unique_coeffs)

# Create a color palette from white to dark green based on the number of unique values
light_green_palette <- colorRampPalette(c("white", "brown"))
color_palette <- light_green_palette(num_unique_values)

# Specify the desired size for the main title
cex_main <- 0.5  # Adjust this value as needed

# Set the global graphical parameter for main title size
par(cex.main = cex_main)

# Reverse the rows of the matrix to flip the y-axis direction
shareBO_reversed <- shareBO_rounded[nrow(shareBO_rounded):1, ]

# Plot the heatmap with the reversed matrix
heatmap(shareBO_reversed, scale = "none", Rowv = NA, Colv = NA,
        col = color_palette,
        main = "Share of budget spent in BO - Reference case",
        cexRow = 0.7, cexCol = 0.7)

# Create a named vector associating each unique value in the matrix with a color
color_map <- setNames(color_palette, unique_coeffs)

# Adjust plot margins using par() function
# Adjust the right margin to allocate more space for the legend
par(mar = c(4, 4, 4, 4))

# Define the coordinates for the legend
# x and y values should range from 0 (bottom/left) to 1 (top/right)
# To move the legend more to the right, increase the x value
# Try x = 1.05 or higher to move the legend further right of the plot
legend_x <- 0.8  # Increase the value as needed to move the legend further right
legend_y <- 0.8  # You can adjust this as needed

# Create the legend
legend(x = legend_x, y = legend_y,
       legend = unique_coeffs, fill = color_map,
       cex = 0.5, pt.cex = 1.5, # Adjust cex and pt.cex as needed
       y.intersp = 1.5) # Adjust spacing between legend items as needed






# (To be improved) Let us try to generalize this code with a function enabling us to perform the same tasks for all goods

# Define a function to plot heatmaps with specified parameters
plot_heatmap <- function(matrix, row_names, col_names, title = "Heatmap", decimal_places = 8, tolerance = 1e-6,
                         cex_main = 0.5, legend_x = 1.05, legend_y = 0.5) {
  
  # Set the row and column names
  rownames(matrix) <- row_names
  colnames(matrix) <- col_names
  
  # Round the matrix values to the specified number of decimal places
  rounded_matrix <- round(matrix, decimal_places)
  
  # Get the unique coefficients using the grouping function with tolerance
  unique_coeffs <- sort(group_values(rounded_matrix, tolerance))
  
  # Determine the number of unique values in the rounded matrix
  num_unique_values <- length(unique_coeffs)
  
  # Create a color palette from white to dark green based on the number of unique values
  light_green_palette <- colorRampPalette(c("white", "darkgreen"))
  color_palette <- light_green_palette(num_unique_values)
  
  # Reverse the rows of the matrix to flip the y-axis direction
  reversed_matrix <- rounded_matrix[nrow(rounded_matrix):1, ]
  
  # Plot the heatmap with the reversed matrix
  heatmap(reversed_matrix, scale = "none", Rowv = NA, Colv = NA,
          col = color_palette,
          main = title,
          cexRow = 0.7, cexCol = 0.7)
  
  # Create a named vector associating each unique value in the matrix with a color
  color_map <- setNames(color_palette, unique_coeffs)
  
  # Create the legend with the given x and y coordinates
  legend(x = legend_x, y = legend_y,
         legend = unique_coeffs, fill = color_map,
         cex = 0.5, pt.cex = 1.5, # Adjust cex and pt.cex as needed
         y.intersp = 1.5, title = "Legend")
}

# Define a list of matrices to be plotted
matrix_list <- list(
 shareGO,
 shareGD,
 shareBO,
 shareBD
  # Add more matrices as needed
)

# Define the desired row and column names
row_names <- rev_b
col_names <- a

# Iterate through the list of matrices and plot heatmaps
for (matrix_name in names(matrix_list)) {
  matrix <- matrix_list[[matrix_name]]
  
  # Customize the title as needed
  title <- paste("Share of expenses in good", matrix_name)
  
  # Call the plot_heatmap function to plot each matrix
  plot_heatmap(matrix, row_names, col_names, title = title)
}

## Not working for the moment

















