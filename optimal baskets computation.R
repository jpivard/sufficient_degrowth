########## Numerical simulation code for a model of lifestyle choices with environmental and image preferences - Global optima and scenario computations ###########


#In this file, we aim to compute the optimal consumption baskets of consumers depending on their preference parameters alpha and beta (and the exogenous parameters)
#These baskets can be composed of one to four goods, so that we take into account mixed consumption, unlike in the first simulations
#To do so, we use a loop generation approach, making it possible to determine for each couple of coordinates (alpha,beta) the composition of the basket maximizing the associated utility
#We then deduce optimal levels of consumption of each type of good and associated pollution
#We finally compare such aggregates across a range of population distribution scenarios


####### Preliminary #####

#To clean space
rm(list=ls())

###Packages to be installed

#install.packages("rgl")
#install.packages("reshape")
#install.packages("ggplot2") 
#install.packages("plotly")  
#install.packages("gridExtra")
#install.packages("cowplot")
#install.packages("parallel")
#install.packages("doParallel")
#install.packages("foreach")
#install.packages("ggtext")
#install.packages("fields")
#install.packages("patchwork")

# Load necessary libraries

library(graphics)
library(rgl)
library("reshape")
library("ggplot2")
library(tibble)
library("plotly")     
library(RColorBrewer)
library("gridExtra")
library(cowplot)
library(grDevices)
library(parallel)
library(doParallel)
library(foreach)
library("ggtext")
library(fields)
library(grid)
library(patchwork)
library(png)
library(dplyr)

###### We first calibrate the parameters (to be changed according to the case tested) #####

R = 16
p_go = 3
p_gd = 2
p_bo = 1.9
p_bd = 1
p <- c(p_go,p_gd,p_bo,p_bd)

d_prime = 0.1
gammaGO = 0.15
gammaGD = 0.05
gammaBD = 0.2
gammaBO = 0.5
gamma <- c(gammaGO, gammaGD, gammaBO, gammaBD)


### The size of the matrices (i.e. the accuracy of the computation) must also be chosen here.
step = 0.025  #smallest possible step chosen in order to be accurate enough while having reasonable computation times with paralleled code
size = 1/step + 1


###### The optimization algorithm #######

### Let us define an optimization function that loops over all possible baskets for each pair of parameters and chooses for each the one maximizing utility.

# Utility Function Definition
utility_function <- function(x1, x2, x3, x4, alpha, beta, gamma, d_prime) {
  u_cons <- log(x1 + x2 + x3 + x4 + .Machine$double.eps)
  u_damage <- -alpha * d_prime * (gamma[1] * x1 + gamma[2] * x2 + gamma[3] * x3 + gamma[4] * x4)
  v_status <- beta * (alpha * (x1 + x2) + (1 - alpha) * (x1 + x3))
  return(u_cons + u_damage + v_status)
}

# Optimization Function
opti_function <- function(alpha, beta, R, p, gamma, d_prime, step) {
  X <- seq(0, R / min(p), by = step)
  Umax <- -Inf
  xs <- numeric(4)
  
  for (x1 in X) {
    for (x2 in X) {
      if (p[1] * x1 + p[2] * x2 <= R) {
        for (x3 in X) {
          if (p[1] * x1 + p[2] * x2 + p[3] * x3 <= R) {
            x4 <- (R - p[1] * x1 - p[2] * x2 - p[3] * x3) / p[4]
            if (x4 >= 0) {
              U <- utility_function(x1, x2, x3, x4, alpha, beta, gamma, d_prime)
              if (U > Umax) {
                Umax <- U
                xs <- c(x1, x2, x3, x4)
              }
            }
          }
        }
      }
    }
  }
  
  return(list(xopt = xs, yopt = Umax))
}

# Parallel Grid Search Function
grid_search_parallel <- function(R, p, gamma, d_prime, step) {
  alphas <- seq(0, 1, by = step)
  betas <- seq(0, 1, by = step)
  rev_betas <- rev(betas)
  size <- length(alphas)
  
  A <- matrix(rep(alphas, each = size), nrow = size)
  B <- matrix(rep(rev_betas, size), nrow = size)
  
  # Prepare storage for results
  results <- vector("list", size * size)
  
  # Create a cluster using all available cores
  cl <- makeCluster(detectCores() - 1) # Use one less core than available
  
  # Export necessary variables and functions to the cluster
  clusterExport(cl, c("A", "B", "R", "p", "gamma", "d_prime", "step", "opti_function", "utility_function"))
  
  # Define the parallelized grid search operation
  grid_search_task <- function(index) {
    i <- ((index - 1) %/% size) + 1
    j <- ((index - 1) %% size) + 1
    
    alpha <- A[i, j]
    beta <- B[i, j]
    
    optibasket <- opti_function(alpha, beta, R, p, gamma, d_prime, step)
    xopt <- optibasket$xopt
    yopt <- optibasket$yopt
    l <- which(xopt > 0)
    
    result <- list(D = 0, D1 = 0, D2 = 0, D3 = 0, D4 = 0, M = yopt, DD = rep(0, 4))
    if (length(l) > 0) {
      result$D <- sum(2^(l-1))
      result$D1 <- xopt[1]
      result$D2 <- xopt[2]
      result$D3 <- xopt[3]
      result$D4 <- xopt[4]
      for (k in l) {
        result$DD[k] <- 1
      }
    }
    
    return(list(i = i, j = j, result = result))
  }
  
  # Run the grid search in parallel
  results <- parLapply(cl, 1:(size * size), grid_search_task)
  
  # Stop the cluster
  stopCluster(cl)
  
  # Initialize result matrices
  D <- matrix(0, nrow = size, ncol = size)
  D1 <- matrix(0, nrow = size, ncol = size)
  D2 <- matrix(0, nrow = size, ncol = size)
  D3 <- matrix(0, nrow = size, ncol = size)
  D4 <- matrix(0, nrow = size, ncol = size)
  M <- matrix(0, nrow = size, ncol = size)
  DD <- lapply(1:4, function(x) matrix(0, nrow = size, ncol = size))
  
  # Populate result matrices
  for (res in results) {
    i <- res$i
    j <- res$j
    D[i, j] <- res$result$D
    D1[i, j] <- res$result$D1
    D2[i, j] <- res$result$D2
    D3[i, j] <- res$result$D3
    D4[i, j] <- res$result$D4
    M[i, j] <- res$result$M
    for (k in 1:4) {
      DD[[k]][i, j] <- res$result$DD[k]
    }
  }
  
  return(list(D = D, D1 = D1, D2 = D2, D3 = D3, D4 = D4, M = M, DD = DD, A = A, B = B, alphas = alphas, rev_betas = rev_betas))
}


results <- grid_search_parallel(R, p, gamma, d_prime, step)

# Access results and ensure A, B, alphas, and rev_betas are available globally
D <- results$D
D1 <- results$D1
D2 <- results$D2
D3 <- results$D3
D4 <- results$D4
M <- results$M
DD <- results$DD
A <- results$A
B <- results$B
alphas <- results$alphas
rev_betas <- results$rev_betas


#Compute the optimal budget shares spent in each good (in percentage) according to intrinsic preferences
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


###### Graphical representation of the quantity matrices ########


### Correcting for algorithmic approximations

#As there seems to be a progressive convergence to 100% for some points theoretically located in ostentatious exclusive zones when matrices get bigger, we first need to check utilities when more than 99% of the budget is spent on a given good and compare them with the utility of allocating the whole budget on the good
#Current utility levels (before double checking those pixels) are stored in the M matrix.
#We would like to compare those levels to utility levels that would be achieved in the theoretical case of exclusive consumption (absence of consumption for BD) for all coefficients equal to one in the dummy matrix.
#In particular pixels where consumption of BD is very close to zero could actually be zeros (artifact due to the budget allocation method)

##Step 1 : Create for each good a dummy matrix that locates the pixels that are concerned by the problem

# Initialize the matrices 'checkexclusives' and 'checkzeros' with the same dimensions as 'A'
n_rows <- nrow(A)
n_cols <- ncol(A)
checkexclusivesGO <- matrix(NA, nrow = n_rows, ncol = n_cols)
checkexclusivesGD <- matrix(NA, nrow = n_rows, ncol = n_cols)
checkexclusivesBO <- matrix(NA, nrow = n_rows, ncol = n_cols)
checkzerosBD <- matrix(NA, nrow = n_rows, ncol = n_cols)


# We first create matrices with positive coefficients only in the concerned pixels
for (i in 1:nrow(shareGO)) {
  for (j in 1:ncol(shareGO)) {
    checkexclusivesGO[i,j] <- shareGO[i,j] - 98.99
    checkexclusivesGD[i,j] <- shareGD[i,j] - 98.99
    checkexclusivesBO[i,j] <- shareBO[i,j] - 98.99
    checkzerosBD[i,j] <- 1 - shareBD[i,j]  # For BD, we need to check absence of consumption rather than exclusive consumption as we see many points with tiny shares of the budget being spent in the good
  }
}

#Convert to dummy matrices where positive coefficients are converted into 1s and the rest to zeros
create_dummy_matrix <- function(matrix_input) {
  # Create the dummy matrix
  dummy_matrix <- as.numeric(matrix_input > 0)
  # Convert the dummy_matrix to the same dimensions as the original matrix
  dummy_matrix <- matrix(dummy_matrix, nrow = nrow(matrix_input), ncol = ncol(matrix_input))
  # Return the dummy matrix
  return(dummy_matrix)
}

dummy_GOcheck <- create_dummy_matrix(checkexclusivesGO)
dummy_GDcheck <- create_dummy_matrix(checkexclusivesGD)
dummy_BOcheck <- create_dummy_matrix(checkexclusivesBO)
dummy_BDcheck <- create_dummy_matrix(checkzerosBD)
#For each good, we now have all the pixels where we need to compare utilities outside the grid search.


# Initialize matrices to track which utility value is maximized
maximizing_GO <- matrix(FALSE, nrow = nrow(A), ncol = ncol(A))
maximizing_GD <- matrix(FALSE, nrow = nrow(A), ncol = ncol(A))
maximizing_BO <- matrix(FALSE, nrow = nrow(A), ncol = ncol(A))
maximizing_BD <- matrix(FALSE, nrow = nrow(A), ncol = ncol(A))

# Function to compute utility when only one good is consumed (except D4 where it is not consumed)
compute_utility_exclusive <- function(good_index, alpha, beta, gamma, d_prime, theta, R, p) {
  x_new <- rep(0, 4)
  if (good_index != 4) {
    x_new[good_index] <- R / p[good_index]
  }
  U_exclusive <- utility_function(x_new[1], x_new[2], x_new[3], x_new[4], alpha, beta, gamma, d_prime, theta)
  return(list(U_exclusive = U_exclusive, x_new = x_new))
}

# Function to compute utility with x4 set to zero
compute_utility_no_x4 <- function(alpha, beta, gamma, d_prime, theta, R, p, D1, D2, D3) {
  x_new <- c(D1, D2, D3, 0)  # Set x4 to zero
  U_no_x4 <- utility_function(x_new[1], x_new[2], x_new[3], x_new[4], alpha, beta, gamma, d_prime, theta)
  return(list(U_no_x4 = U_no_x4, x_new = x_new))
}

# Function to compute utility with various allocations of residual revenue
compute_utility_with_allocations <- function(alpha, beta, gamma, d_prime, theta, R, p, D1, D2, D3) {
  max_utility <- -Inf
  best_allocation <- c(D1, D2, D3, 0)
  
  # Get indices of goods that have positive quantities
  positive_indices <- which(c(D1, D2, D3) > 0)
  
  # Total initial consumption expenditure
  initial_expenditure <- sum(p[1:3] * c(D1, D2, D3))
  
  # Residual revenue to allocate
  residual_revenue <- R - initial_expenditure
  
  # Generate all combinations of allocating residual revenue to positive goods
  if (length(positive_indices) > 0) {
    n <- length(positive_indices)
    step <- residual_revenue / 10000  # Change this step size as needed for more precision or speed
    for (alloc in expand.grid(rep(list(seq(0, residual_revenue, by = step)), n))) {
      if (sum(alloc) <= residual_revenue) {
        allocation <- c(D1, D2, D3)
        allocation[positive_indices] <- allocation[positive_indices] + alloc / p[positive_indices]
        utility <- utility_function(allocation[1], allocation[2], allocation[3], 0, alpha, beta, gamma, d_prime, theta)
        if (utility > max_utility) {
          max_utility <- utility
          best_allocation <- c(allocation, 0)
        }
      }
    }
  }
  
  return(list(U_best = max_utility, x_new = best_allocation))
}

# Loop through the elements of the dummy matrices
for (i in 1:nrow(dummy_GOcheck)) {
  for (j in 1:ncol(dummy_GOcheck)) {
    if (dummy_GOcheck[i, j] == 1) {
      alpha <- A[i, j]
      beta <- B[i, j]
      
      # For dummy_GOcheck
      result <- compute_utility_exclusive(1, alpha, beta, gamma, d_prime, theta, R, p)
      if (result$U_exclusive > M[i, j]) {
        M[i, j] <- result$U_exclusive
        D1[i, j] <- result$x_new[1]
        D2[i, j] <- result$x_new[2]
        D3[i, j] <- result$x_new[3]
        D4[i, j] <- result$x_new[4]
        maximizing_GO[i, j] <- TRUE
      }
    }
    
    if (dummy_GDcheck[i, j] == 1) {
      alpha <- A[i, j]
      beta <- B[i, j]
      
      # For dummy_GDcheck
      result <- compute_utility_exclusive(1, alpha, beta, gamma, d_prime, theta, R, p)
      if (result$U_exclusive > M[i, j]) {
        M[i, j] <- result$U_exclusive
        D1[i, j] <- result$x_new[1]
        D2[i, j] <- result$x_new[2]
        D3[i, j] <- result$x_new[3]
        D4[i, j] <- result$x_new[4]
        maximizing_GD[i, j] <- TRUE
      }
    }
    
    
    if (dummy_BOcheck[i, j] == 1) {
      alpha <- A[i, j]
      beta <- B[i, j]
      
      # For dummy_BOcheck
      result <- compute_utility_exclusive(3, alpha, beta, gamma, d_prime, theta, R, p)
      if (result$U_exclusive > M[i, j]) {
        M[i, j] <- result$U_exclusive
        D1[i, j] <- result$x_new[1]
        D2[i, j] <- result$x_new[2]
        D3[i, j] <- result$x_new[3]
        D4[i, j] <- result$x_new[4]
        maximizing_BO[i, j] <- TRUE
      }
    }
    
    if (dummy_BDcheck[i, j] == 1) {
      alpha <- A[i, j]
      beta <- B[i, j]
      
      # For dummy_BDcheck
      result <- compute_utility_no_x4(alpha, beta, gamma, d_prime, theta, R, p, D1[i, j], D2[i, j], D3[i, j])
      if (result$U_no_x4 > M[i, j]) {
        M[i, j] <- result$U_no_x4
        D4[i, j] <- 0
        maximizing_BD[i, j] <- TRUE
      }
    }
  }
}

# Recalculate the DD list to reflect updated quantities
for (i in 1:nrow(DD[[1]])) {
  for (j in 1:ncol(DD[[1]])) {
    DD[[1]][i, j] <- ifelse(D1[i, j] > 0, 1, 0)
    DD[[2]][i, j] <- ifelse(D2[i, j] > 0, 1, 0)
    DD[[3]][i, j] <- ifelse(D3[i, j] > 0, 1, 0)
    DD[[4]][i, j] <- ifelse(D4[i, j] > 0, 1, 0)
  }
}

# Check for D4 coefficients still below 0.01 but different from zero and update matrices again if necessary
for (i in 1:nrow(D4)) {
  for (j in 1:ncol(D4)) {
    if (D4[i, j] < 0.01 && D4[i, j] > 0) {
      alpha <- A[i, j]
      beta <- B[i, j]
      
      # Compare utilities for optimal distribution excluding D4
      result <- compute_utility_with_allocations(alpha, beta, gamma, d_prime, theta, R, p, D1[i, j], D2[i, j], D3[i, j])
      current_utility <- M[i, j]
      
      if (result$U_best > current_utility) {
        M[i, j] <- result$U_best
        D1[i, j] <- result$x_new[1]
        D2[i, j] <- result$x_new[2]
        D3[i, j] <- result$x_new[3]
        D4[i, j] <- 0
      }
    }
  }
}

# Now M contains the updated utility levels and D1-D4 contain the updated quantities
comparison_results <- list(M = M, D1 = D1, D2 = D2, D3 = D3, D4 = D4, DD = DD)

#Now recompute the budget shares 
for (i in 1:nrow(shareGO)) {
  for (j in 1:ncol(shareGO)) {
    shareGO[i,j]<- (p[1]*D1[i,j]/R)*100
    shareGD[i,j]<- (p[2]*D2[i,j]/R)*100
    shareBO[i,j]<- (p[3]*D3[i,j]/R)*100
    shareBD[i,j]<- (p[4]*D4[i,j]/R)*100
  }
}

#We now have the 'consolidated' budget share matrices after checking all points in the exclusive zones are really optimized.

### Additional correction step in the case of a tax : as prices are not round numbers anymore, the algorithm has to artificially finish allocating small remaining quantities to another good that the preferred one.
#For instance, in GD exclusive zone, between 0.3 and 0.6% can be allocated to another lifestyle (BD or GO)
#Although this does not matter at individual level, this might affect our collective aggregates artificially hence the need for correction
#We correct quantities (not spending because what enters following impact or market shares computations is quantities) so that each time a lifestyle represents less than 1% of consumer spending, this remaining 1% is allocated to the otherwise preferred lifestyle.
#Constraint : quantities should never exceed R/p however, if this is the case the remaining quantity is rounded to R/p


matrices <- list(D1 = D1, D2 = D2, D3 = D3, D4 = D4)  # Store all matrices in a list
thresholds <- R / p  # Compute quantity thresholds for each matrix

for (matrix_index in seq_along(matrices)) {
  matrix_name <- names(matrices)[matrix_index]  # Current matrix name (D1, D2, ...)
  current_matrix <- matrices[[matrix_name]]    # Current matrix
  
  for (i in 1:nrow(current_matrix)) {
    for (j in 1:ncol(current_matrix)) {
      # Step 1: Check if the value is below the threshold and reallocate if needed
      if (current_matrix[i, j] < 0.01 * thresholds[matrix_index]) {
        # Find the matrix with the max value at (i, j) among the other three
        other_matrices <- matrices[-matrix_index]  # Exclude the current matrix
        max_matrix_name <- names(which.max(sapply(other_matrices, function(m) m[i, j])))
        
        # Add the value to the max matrix and set the current value to 0
        matrices[[max_matrix_name]][i, j] <- matrices[[max_matrix_name]][i, j] + current_matrix[i, j]
        current_matrix[i, j] <- 0
      }
    }
  }
  
  # Step 2: Ensure values don't exceed the threshold or fall within 99%-threshold range
  current_matrix <- ifelse(
    current_matrix > 0.99 * thresholds[matrix_index],  # If value exceeds 99% of the threshold
    thresholds[matrix_index],                         # Set to the threshold
    pmin(current_matrix, thresholds[matrix_index])    # Otherwise, cap at the threshold
  )
  
  matrices[[matrix_name]] <- current_matrix  # Update the modified current matrix in the list
}

# Assign the updated matrices back to their respective variables
D1 <- matrices$D1
D2 <- matrices$D2
D3 <- matrices$D3
D4 <- matrices$D4

#And compute the final spending matrices
for (i in 1:nrow(shareGO)) {
  for (j in 1:ncol(shareGO)) {
    shareGO[i,j]<- (p[1]*D1[i,j]/R)*100
    shareGD[i,j]<- (p[2]*D2[i,j]/R)*100
    shareBO[i,j]<- (p[3]*D3[i,j]/R)*100
    shareBD[i,j]<- (p[4]*D4[i,j]/R)*100
  }
}



## Graphical representations of how consumers spend their budget in the different lifestyles

generate_heatmap <- function(share, lifestyle, color_palette, legend_title, add_legend = FALSE) {
  if (!is.null(share)) {
    # Set the row and column names
    rownames(share) <- rev_betas
    colnames(share) <- alphas
    
    # Reverse the rows of the matrix to flip the y-axis direction
    share_reversed <- share[nrow(share):1, ]
    
    # Calculate the maximum dimension to ensure square heatmaps
    max_dim <- max(nrow(share_reversed), ncol(share_reversed))
    
    # Plot the heatmap with the reversed matrix using image
    image(1:ncol(share_reversed), 1:nrow(share_reversed), t(share_reversed), col = color_palette, axes = FALSE, 
          main = paste("% of income spent in", lifestyle, "- BO excluded / High damage"), 
          xlab = "alpha", ylab = "beta", xlim = c(1, max_dim), ylim = c(1, max_dim), asp = 1)
    
    # Add axis labels
    axis(1, at = 1:ncol(share_reversed), labels = colnames(share_reversed))
    axis(2, at = 1:nrow(share_reversed), labels = rownames(share_reversed))
  }
  
  if (add_legend) {
    # Determine the number of colors and percentage cutoffs
    num_colors <- length(color_palette)  # Match the number of intervals
    cutoffs <- seq(0, 100, by = 100 / (num_colors - 1))  # Ensure cutoffs match num_colors
    
    # Create the legend
    legend("center", legend = as.character(cutoffs), fill = color_palette, title = legend_title, 
           cex = 0.7, pt.cex = 1.5, y.intersp = 1.5, xpd = TRUE)
  }
}

# Define color palettes for different lifestyles
num_colors <- 11  # Ensure this matches the number of intervals
color_palette_GO <- colorRampPalette(c("white", "darkgreen"))(num_colors)
color_palette_GD <- colorRampPalette(c("white", "darkgreen"))(num_colors)
color_palette_BO <- colorRampPalette(c("white", "brown"))(num_colors)
color_palette_BD <- colorRampPalette(c("white", "brown"))(num_colors)

# Set up the layout for 2x3 grid of heatmaps and legends
layout(matrix(c(1, 2, 5, 3, 4, 6), nrow = 2, ncol = 3, byrow = TRUE), widths = c(1, 1, 0.3), heights = c(1, 1))

# Set margins for the heatmaps
par(mar = c(4, 4, 4, 0))

# Generate heatmap plots for different lifestyles
generate_heatmap(shareGO, "GO", color_palette_GO, "Shares")
generate_heatmap(shareGD, "GD", color_palette_GD, "Shares")
generate_heatmap(shareBO, "BO", color_palette_BO, "Shares")
generate_heatmap(shareBD, "BD", color_palette_BD, "Shares")

# Add legends
par(mar = c(4, 0, 4, 0))  # Reduce margins for legends
plot.new()
generate_heatmap(NULL, "", color_palette_GO, "Percentages", add_legend = TRUE)
plot.new()
generate_heatmap(NULL, "", color_palette_BO, "Percentages", add_legend = TRUE)

# Reset the graphical parameters to default
par(mfrow = c(1, 1))





##Compute and plot 'market shares' of the different goods

#Start by computing total quantities consumed of each good
quantity_GO = sum(D1)
quantity_GD = sum(D2)
quantity_BO = sum(D3)
quantity_BD = sum(D4)
total_quantity = sum(D1+D2+D3+D4)

#Then infer the 'market share' of each composite good
marketshare_GO = (quantity_GO/total_quantity)*100
marketshare_GD = (quantity_GD/total_quantity)*100
marketshare_BO = (quantity_BO/total_quantity)*100
marketshare_BD = (quantity_BD/total_quantity)*100

#Compute the indicator 'market share of green goods' for each case, to study its relationship with impacts.
marketshare_green = marketshare_GO + marketshare_GD


## Plot those market shares in an histogram

consumer_types <- c("GO", "GD", "BO", "BD")
percentages <- c(marketshare_GO, marketshare_GD, marketshare_BO, marketshare_BD)

# Create a data frame in which goods are ranked from the dirtiest to the cleanest (in absolute, regardless of income effects)
data <- data.frame(Category = factor(consumer_types, levels = c("BO", "BD", "GO", "GD")), 
                   Percentage = percentages)

# Set colors
custom_colors <- c("brown", "lightgrey", "darkgreen", "lightgreen") 

# Create a bar plot (caption to be changed according to the case tested!)
plot <- ggplot(data, aes(x = Category, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
            position = position_stack(vjust = 0.5), 
            color = "black", 
            size = 5) +  # Adjust the size if needed
  labs(title = "Market shares in volume - 50% higher income (R=24), Low damage baseline, Uniform distribution of the population",
       x = "Lifestyles (from the most to the least polluting)",
       y = "Percentage of quantities consumed") +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +  # Format y-axis as percentage and normalize the scale
  theme_minimal() +
  scale_fill_manual(values = custom_colors, guide = "none")

print(plot)



#### The environmental part #####
#From quantities consumed of each good, we can infer the environmental impacts induced by each good (assuming first implicitly a uniform repartition of the population)
#Method : create for each good matrices of environmental impacts per pixel

#Initialization of matrices
P <- matrix(0, nrow = size, ncol = size) 
P1 <- matrix(0, nrow = size, ncol = size) 
P2 <- matrix(0, nrow = size, ncol = size)
P3 <- matrix(0, nrow = size, ncol = size)
P4 <- matrix(0, nrow = size, ncol = size)

#We multiply by 100 the 'total marginal damage' (i.e. accouting for heterogenous impacts across goods) to normalize the marginal damage coefficient to one in the low damage case 
#This way we obtain impacts per basket between 0 and 100, instead of having tiny impact values.
for (i in 1:nrow(P1)) {
  for (j in 1:ncol(P1)) {
    P1[i,j]<- 100*d_prime*gammaGO*D1[i,j]
    P2[i,j]<- 100*d_prime*gammaGD*D2[i,j]
    P3[i,j]<- 100*d_prime*gammaBO*D3[i,j]
    P4[i,j]<- 100*d_prime*gammaBD*D4[i,j]
  }
}

#Compute the total impact per pixel/consumer (to be compared between cases and then plotted)
P <- P1+P2+P3+P4

#And the contribution of each lifestyle to total pollution
impact_GO = sum(P1)
impact_GD = sum(P2)
impact_BO = sum(P3)
impact_BD = sum(P4)
total_impacts = (impact_GO+impact_GD+impact_BO+impact_BD)

#And finally per capita impacts
totalimpacts_percapita = total_impacts/(size^2)

##Store and plot those impacts in the different cases (we leave two in comments for illustration purposes and clean the rest)
#totalimpacts_percapita_reflowdamage_unif <- totalimpacts_percapita
#totalimpacts_percapita_refhighdamage_unif <- totalimpacts_percapita



#Secondary indicator: Plot the relative contribution of each lifestyle to pollution
contrib_GO = (impact_GO/total_impacts)*100
contrib_GD = (impact_GD/total_impacts)*100
contrib_BO = (impact_BO/total_impacts)*100
contrib_BD = (impact_BD/total_impacts)*100

consumer_types <- c("GO", "GD", "BO", "BD")
percentages <- c(contrib_GO, contrib_GD, contrib_BO, contrib_BD)

# Create a data frame
data <- data.frame(Category = consumer_types, Percentage = percentages)

# Set colors
custom_colors <- c("lightgrey", "brown", "lightgreen", "darkgreen") 

# Create a bar plot (caption to be changed according to the case tested)
plot <- ggplot(data, aes(x = Category, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity") +
  labs(title = "Contributions to environmental impacts by the different lifestyles - Low damage baseline, Uniform distribution",
       x = "Lifestyles",
       y = "Proportion of environmental impacts") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Format y-axis as percentage
  theme_minimal() +
  scale_fill_manual(values = custom_colors, guide="none")

print(plot)





#### Representing Pollution/Environmental impacts per basket per consumer type (i.e. depending on preferences)


# Define the matrix P 
rownames(P) <- rev_betas
colnames(P) <- alphas

# Reverse the rows of the matrix to flip the y-axis direction
P <- P[nrow(P):1, ]

# Determine the number of colors and create a continuous color palette
num_colors <- 100  # Use 100 colors for a smooth gradient
color_palette <- colorRampPalette(c("darkgreen", "lightgreen", "yellow", "orange", "red", "brown"))(num_colors)

# Calculate the range and create breaks based on equal intervals
min_val <- min(P)
max_val <- max(P)
breaks <- seq(min_val, max_val, length.out = num_colors + 1)

# Create a named vector associating each unique value in the matrix with a color
color_map <- cut(P, breaks = breaks, include.lowest = TRUE, labels = color_palette)

# Adjust plot margins using par() function before creating the heatmap
# Increase the bottom margin to make space for the x-axis label
par(mar = c(5,6,6,7))  # c(bottom, left, top, right) - increase the bottom margin

# Create the heatmap
heatmap(P, scale = "none", Rowv = NA, Colv = NA,
        col = color_palette,
        main = "Env. impacts of a potential consumer - BO excluded from the market, Low damage",
        cexRow = 0.7, cexCol = 0.7,
        ylab = expression(beta))

# Add the x-axis label after creating the heatmap
mtext(expression(alpha), side = 1, line = 3, las = 1)  # Adjust 'line' parameter to position the label properly

#Optional: Add the y-axis label (only when hidden by the panel)
#mtext(expression(beta), side = 2, line = 3, las = 1)

# Create the legend
legend_labels <- seq(min_val, max_val, length.out = 10)
legend_colors <- colorRampPalette(c("darkgreen", "lightgreen", "yellow", "orange", "red", "brown"))(10)
legend("right", inset = c(-0.1, 0), legend = round(legend_labels, 2), fill = legend_colors,
       title = "Impact units/pixel",  # Add a title to the legend
       cex = 0.7, pt.cex = 1,  # Adjust cex and pt.cex as needed
       y.intersp = 1.5, xpd = NA)  # Allow plotting outside the plot region and adjust spacing between legend items



#### Population scenarios definition ####

#Until here we have assumed the population to be uniformly distributed in the square.
#We now use Beta distributions to model various population concentration scenarios.
#We divide the plane into 9 squares of population concentration that we successively study 
#Careful: only select one block at a time! 
#For instance select parameter values in "Square A" (corresponding to high image and low environmental concerns) and then directly jump to "jump and start here" below
#Then come back for B, C and so on...


#Square A (North-West) : high image, low environment  
a=c=5
b=d=1

#Square B (Center-North) : high image, medium environment
a=b=c=5
d=1

#Square C (North-East) : high image, high environment
a=d=1
b=c=5

#Square D (Center-West) : medium image, low environment
a=c=d=5
b=1

#Square E (Center) : medium image, average environment
a=b=c=d=5

#Square F (Center-East) : medium image, high environment
b=c=d=5
a=1

#Square G (South-West) : low image, low environment
a=d=5
b=c=1

#Square H (Center-South) : low image, medium environment
c=1
a=b=d=5

#Square I (South-East) : low image, high environment
a=c=1
b=d=5



## Alternative scenarios (for sensitivity analysis): 

#1. Smaller concentration zones

#Square A2 : very high image, very low environment
a=c=8
b=d=1

#Square B2 (Center-North) : very high image, medium environment
a=b=c=8
d=1

#Square C2 (North-East) : very high image, very high environment
a=d=1
b=c=8

#Square D2 (Center-West) : "very medium" image, very low environment
a=c=d=8
b=1

#Square E2 (Center) :"very medium" image, medium environment - more concentrated
a=b=c=d=8

#Square F2 (Center-East) : "very medium" image, very high environment
b=c=d=8
a=1

#Square G2 (South-West) : very low image, very low environment
a=d=8
b=c=1

#Square H2 (Center-South) : very low image, medium environment
c=1
a=b=d=8

#Square I2 (South-East) : very low image, very high environment
a=c=1
b=d=8


#2. Bigger concentration zones (less 'extreme' societies)

#Square A3 : medium/high image, medium/low environment
a=c=4
b=d=2

#Square B3 (Center-North) 
a=b=c=4
d=2

#Square C3 (North-East) 
a=d=2
b=c=4

#Square D3 (Center-West) 
a=c=d=4
b=2

#Square E3 (Center) : medium image, medium environment - more dispersed
a=b=c=d=2

#Square F3 (Center-East) 
b=c=d=4
a=2

#Square G3 (South-West) 
a=d=4
b=c=2

#Square H3 (Center-South) 
c=2
a=b=d=4

#Square I3 (South-East) 
a=c=2
b=d=4




### Jump and start here with your values for each scenario

##We now compute population matrices according to the chosen scenario

#initialize  with an empty matrix
pop <- matrix(,size,size) 

##To compute discrete densities: 

#We first create our quantile vectors (in line)

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

# Calling the functions at n = 0 and storing the result
result_at_0_dim1 <- densite_beta_dim1(0, a, b, step)
result_at_0_dim2 <- densite_beta_dim2(0, c, d, step)

# Storing other values iteratively in a matrix 
num_iterations <- size - 1  
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
#bell-shaped densities, symmetrical and sum to one 

#We multiply our two probability vectors to obtain the probabilities for the entire population.
pop= t(proba_values_dim2)%*%proba_values_dim1
#There is indeed a higher population concentration on median values compared to the uniform case.

sum(pop) 
#just to check it sums to one




###Plot the population concentration in a heatmap for each scenario

# Define the matrix P with rownames and colnames
rownames(pop) <- rev_betas
colnames(pop) <- alphas

# Reverse the rows of the matrix to flip the y-axis direction for the figure - to be executed once only !
pop <- pop[nrow(pop):1, ]

# Determine the number of colors and create a continuous color palette
num_colors <- 100  # Use 100 colors for a smooth gradient
color_palette <- colorRampPalette(c("lightblue", "blue", "darkblue"))(num_colors)

# Plot the heatmap without a title
heatmap(
  pop, 
  scale = "none", 
  Rowv = NA, 
  Colv = NA,
  col = color_palette,
  main = expression("Density of population - Scenario " * (alpha[l] * "," * beta[h])* ""),  #Example for square A (that we name alpha_l, beta_h in the paper)
  cexRow = 0.7, 
  cexCol = 0.7,
  ylab = expression(beta)
)


# Add the x-axis label after creating the heatmap
mtext(expression(alpha), side = 1, line = 3, las = 1)  # Adjust 'line' parameter to position the label properly

#Put back the rows in the other direction for correct market shares computation
pop <- pop[nrow(pop):1, ]




### Application of these population scenarios to compute environmental impacts

#We first define the population of consumers of each lifestyle depending on the scenario

# Step 1 is to define a function to create a dummy matrix corresponding to consumers (or not) of the lifestyle
dummy_GO <- create_dummy_matrix(D1)
dummy_GD <- create_dummy_matrix(D2)
dummy_BO <- create_dummy_matrix(D3)
dummy_BD <- create_dummy_matrix(D4)

#Step 2 is to compute the number of people following each lifestyle depending on the population distribution
go_consumers = pop*dummy_GO
gd_consumers = pop*dummy_GD
bo_consumers = pop*dummy_BO
bd_consumers = pop*dummy_BD

#Before computing per capita impacts, we must rescale population matrices by the size of the population (the square of the size variable) so that a bit more than 1 individual lie in each pixel in a concentrated area, whereas one can find values between 0 and 1 when one gets further from the concentrated area
go_consumers_rescaled = go_consumers*size^2
gd_consumers_rescaled = gd_consumers*size^2
bo_consumers_rescaled = bo_consumers*size^2
bd_consumers_rescaled = bd_consumers*size^2

#We can then compute the quantities after accounting for the change in population distribution
#This will yield us the market shares for other population distributions
quantity_GO_nonuniform = sum(D1*go_consumers_rescaled)
quantity_GD_nonuniform = sum(D2*gd_consumers_rescaled)
quantity_BO_nonuniform = sum(D3*bo_consumers_rescaled)
quantity_BD_nonuniform = sum(D4*bd_consumers_rescaled)
total_quantity_nonuniform = quantity_GO_nonuniform+quantity_GD_nonuniform+quantity_BO_nonuniform+quantity_BD_nonuniform

#Then infer the 'market share' of each composite good, store them for figures (we only leave one for example and clean the rest)
marketshare_GO_nonuniform = (quantity_GO_nonuniform/total_quantity_nonuniform)*100
marketshare_GD_nonuniform = (quantity_GD_nonuniform/total_quantity_nonuniform)*100
marketshare_BO_nonuniform = (quantity_BO_nonuniform/total_quantity_nonuniform)*100
marketshare_BD_nonuniform = (quantity_BD_nonuniform/total_quantity_nonuniform)*100

green_marketshare_nonuniform= marketshare_GO_nonuniform + marketshare_GD_nonuniform 
#green_marketshare_scenA_lowinclowdam <- green_marketshare_nonuniform

##Then plot those market shares in an histogram (only thing to change is the title according to the scenario/case)

consumer_types <- c("GO", "GD", "BO", "BD")
percentages <- c(marketshare_GO_nonuniform, marketshare_GD_nonuniform, marketshare_BO_nonuniform, marketshare_BD_nonuniform)

# Create a data frame
data <- data.frame(Category = factor(consumer_types, levels = c("BO", "BD", "GO", "GD")), 
                   Percentage = percentages)

# Set colors
custom_colors <- c("brown", "lightgrey", "darkgreen", "lightgreen") 

# Create a bar plot (caption to be changed according to the case tested)
plot <- ggplot(data, aes(x = Category, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
            position = position_stack(vjust = 0.5), 
            color = "black", 
            size = 5) +  # Adjust the size if needed
  labs(title = "Market shares in volume - High environmental concern, Low image concern - Higher income (R=24), Low damage",
       x = "Lifestyles (from the most to the least polluting)",
       y = "Percentage of quantities consumed") +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +  # Format y-axis as percentage and normalize the scale
  theme_minimal() +
  scale_fill_manual(values = custom_colors, guide="none")

print(plot)


## Let's finally compute the associated per capita impacts.

#We can easily adapt the pollution matrices according to the population scenario
P1_nonunif = P1*go_consumers_rescaled
P2_nonunif = P2*gd_consumers_rescaled
P3_nonunif = P3*bo_consumers_rescaled
P4_nonunif = P4*bd_consumers_rescaled

#And infer the new contribution of each lifestyle to total pollution
impact_GO_nonunif = sum(P1_nonunif)
impact_GD_nonunif = sum(P2_nonunif)
impact_BO_nonunif = sum(P3_nonunif)
impact_BD_nonunif = sum(P4_nonunif)
total_impacts_nonunif = impact_GO_nonunif+impact_GD_nonunif+impact_BO_nonunif+impact_BD_nonunif

#And finally per capita impacts, that we compare to the uniform scenario
#NB: Next lines to be changed depending on cases !
totalimpacts_percapita_higherincomelowdamage_nonunif = total_impacts_nonunif/(size^2)
totalimpacts_percapita_reflowdamage_unif = totalimpacts_percapita  
variation_with_unif = ((totalimpacts_percapita_hugeincome_lowdamage_nonunif- totalimpacts_percapita)/totalimpacts_percapita)*100


#Secondary indicator: Plot the relative contribution of each lifestyle to pollution
contrib_GO_nonunif = (impact_GO_nonunif/total_impacts_nonunif)*100
contrib_GD_nonunif = (impact_GD_nonunif/total_impacts_nonunif)*100
contrib_BO_nonunif = (impact_BO_nonunif/total_impacts_nonunif)*100
contrib_BD_nonunif = (impact_BD_nonunif/total_impacts_nonunif)*100

consumer_types <- c("GO", "GD", "BO", "BD")
percentages <- c(contrib_GO_nonunif, contrib_GD_nonunif, contrib_BO_nonunif, contrib_BD_nonunif)

# Create a data frame
data <- data.frame(Category = consumer_types, Percentage = percentages)

# Set colors
custom_colors <- c("lightgrey", "brown", "lightgreen", "darkgreen") 

# Create a bar plot (caption to be changed according to the case tested)
plot <- ggplot(data, aes(x = Category, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity") +
  labs(title = "Contributions to environmental impacts - Reference case, Middle concentration",
       x = "Lifestyles",
       y = "Proportion of environmental impacts") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Format y-axis as percentage
  theme_minimal() +
  scale_fill_manual(values = custom_colors, guide="none")

print(plot)








###Compare cases across scenarios to study the interaction between socioeconomic and cultural conditions

#Store and plot per capita impacts in the different scenarios. Again we only leave one example

#totalimpacts_percapita_squareA_reflowdamage <- totalimpacts_percapita_reflowdamage_nonunif
#totalimpacts_percapita_squareB_reflowdamage <- totalimpacts_percapita_reflowdamage_nonunif
#totalimpacts_percapita_squareC_reflowdamage <- totalimpacts_percapita_reflowdamage_nonunif
#totalimpacts_percapita_squareD_reflowdamage <- totalimpacts_percapita_reflowdamage_nonunif
#totalimpacts_percapita_squareE_reflowdamage <- totalimpacts_percapita_reflowdamage_nonunif
#totalimpacts_percapita_squareF_reflowdamage <- totalimpacts_percapita_reflowdamage_nonunif
#totalimpacts_percapita_squareG_reflowdamage <- totalimpacts_percapita_reflowdamage_nonunif
#totalimpacts_percapita_squareH_reflowdamage <- totalimpacts_percapita_reflowdamage_nonunif
#totalimpacts_percapita_squareI_reflowdamage <- totalimpacts_percapita_reflowdamage_nonunif






### Codes for the figures ####

#NB: for each figure to be reproduced as in the paper, one should first store the impact numbers in variables and adapt the code accordingly (either leaving it untouched if using same variable names, or modifying accordingly)



#Let's compare cases across scenarios, zooming on different population concentration zones

##First: Uniform distribution of consumers and compare cases only

# Categories and per capita impacts
categories <- c("Baseline-Low damage", "Baseline-High damage", "Low income-Low damage", "Low income-Low high damage", "Pigovian tax-Low damage", "Pigovian tax-High damage")
per_capita_impacts <- c(
  totalimpacts_percapita_reflowdamage_unif,
  totalimpacts_percapita_refhighdamage_unif,
  totalimpacts_percapita_lowerincomelowdamage_unif,
  totalimpacts_percapita_lowincomehighdamage_unif,
  totalimpacts_percapita_taxlowdamage_unif,
  totalimpacts_percapita_taxhighdamage_unif
)

# Create a data frame
data <- data.frame(Category = factor(categories, levels = c("Baseline-Low damage", "Baseline-High damage", "Low income-Low damage", "Low income-Low high damage", "Pigovian tax-Low damage", "Pigovian tax-High damage")), 
                   Numbers = per_capita_impacts)

# Reorder the levels of Category based on Numbers in descending order
data$Category <- factor(data$Category, levels = data$Category[order(-data$Numbers)])

# Set custom colors for the gradient fill
custom_gradient_colors <- colorRampPalette(c("green", "tan", "saddlebrown"))(length(per_capita_impacts))

# Create the bar plot with ggplot2 and include demarcation
library(ggplot2)

plot_pcimpact_unif <- ggplot(data, aes(x = Category, y = Numbers, fill = Numbers)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Per capita impacts by case in two states of the world (high or low damage), Uniform distribution of preferences",
    x = "Cases",
    y = "Impacts"
  ) +
  theme_minimal() +
  scale_fill_gradientn(
    colors = custom_gradient_colors,
    guide = "legend",
    name = "Impacts"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # For better readability
    plot.margin = unit(c(1, 1, 2.5, 1), "cm") # Extend the bottom margin for annotations
  ) +
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "black", size = 1) + # Adds a dashed line for demarcation
  annotate("text", x = 1.75, y = -0.1 * max(data$Numbers), label = "High damage world", size = 4, hjust = 0.5) + # Text for the left side
  annotate("text", x = 5.25, y = -0.1 * max(data$Numbers), label = "Low damage world", size = 4, hjust = 0.5)   # Text for the right side

# Display the plot
print(plot_pcimpact_unif)


## Alternative figure: plot those impacts against the share of green goods in the corresponding case.

#Compute cumulative market shares of green lifestyles
#Here we do it manually (reading market share figures plotted above) not to rerun everything in my case, but this isn automatized for figures below
marketshare_green_lowdbase <- 37.1 + 6.1
marketshare_green_highdbase <- 43.8 + 8.6
marketshare_green_lowdtax <- 37.9 + 6.7
marketshare_green_highdtax <- 46.2 + 10.4
marketshare_green_lowdlowinc <- 35.1 + 3.9 
marketshare_green_highdlowinc <- 41.3 + 5.8

# Define axes
green_marketshares <- c(marketshare_green_lowdbase, marketshare_green_highdbase, 
                        marketshare_green_lowdtax, marketshare_green_highdtax, 
                        marketshare_green_lowdlowinc, marketshare_green_highdlowinc)

per_capita_impacts <- c(
  totalimpacts_percapita_reflowdamage_unif,
  totalimpacts_percapita_refhighdamage_unif,
  totalimpacts_percapita_taxlowdamage_unif,
  totalimpacts_percapita_taxhighdamage_unif,
  totalimpacts_percapita_lowerincomelowdamage_unif,
  totalimpacts_percapita_lowincomehighdamage_unif
)

cases <- c("Baseline-Low damage", "Baseline-High damage",  
           "Pigovian tax-Low damage", "Pigovian tax-High damage", 
           "Low income-Low damage", "Low income-High damage")

# Create a data frame
data <- data.frame(green_marketshares, per_capita_impacts, cases)

# Order data by market share (ascending)
data <- data[order(data$green_marketshares),]

# Split data into Low damage and High damage cases
low_damage <- data[grepl("Low damage", data$cases), ]
high_damage <- data[grepl("High damage", data$cases), ]

# Generate a color gradient based on impacts
custom_gradient_colors <- colorRampPalette(c("green", "tan", "saddlebrown"))(length(per_capita_impacts))

# Plot base (without lines)
plot(data$green_marketshares, data$per_capita_impacts, 
     pch=16, col=custom_gradient_colors, type="n", # "n" to avoid default line
     xlab="Cumulative Market Share of Green Lifestyles (%)", 
     ylab="Average Impact", 
     main="Average Impacts vs. Green Market Shares", 
     xlim=c(30, 70), ylim=range(per_capita_impacts), xaxt='n')

# Add separate lines
lines(low_damage$green_marketshares, low_damage$per_capita_impacts, col="green", lwd=2, type="b", pch=16)
lines(high_damage$green_marketshares, high_damage$per_capita_impacts, col="saddlebrown", lwd=2, type="b", pch=16)

# Custom x-axis labels
axis(1, at = seq(30, 70, by=5), labels = seq(30, 70, by=5), las=1)

# Add text labels to the left of the points
text(data$green_marketshares, data$per_capita_impacts, labels=data$cases, pos=2, cex=0.8, col="black")

# Add legend
legend("topright", legend=c("Low damage cases", "High damage cases"), 
       col=c("green", "saddlebrown"), lwd=2, pch=16, bty="n")



###Compare effect of the tax across scenarios. (Old figure, removed from the latest version but still interesting so that we leave it)

##Under Hyp 4A (ostentatious pollute more)

# Define custom labels with Greek symbols manually
custom_labels2 <- c(
  expression("("*alpha[m]*", "*beta[m]*")"),
  expression("("*alpha[m]*", "*beta[l]*")"),
  expression("("*alpha[l]*", "*beta[l]*")"),
  "Uniform",
  expression("("*alpha[l]*", "*beta[h]*")"),
  expression("("*alpha[h]*", "*beta[l]*")"),
  expression("("*alpha[h]*", "*beta[m]*")")
)

# Define categories and impacts
categories <- rep(c("Baseline", "Pigovian tax"), 7)  # 14 categories, 7 pairs
group_labels <- rep(1:7, each = 2)  # Group pairs to apply custom labels

# Define impact values for each bar (Baseline and Pigovian tax for each pair)
per_capita_impacts <- c(
  totalimpacts_percapita_refhighdamage_unif,
  totalimpacts_percapita_taxhighdamage_unif,
  totalimpacts_percapita_squareA_refhighdamage,
  totalimpacts_percapita_squareA_taxhighdamage,
  totalimpacts_percapita_squareG_refhighdamage,
  totalimpacts_percapita_squareG_taxhighdamage,
  totalimpacts_percapita_squareE_refhighdamage,
  totalimpacts_percapita_squareE_taxhighdamage,
  totalimpacts_percapita_squareH_refhighdamage,
  totalimpacts_percapita_squareH_taxhighdamage,
  totalimpacts_percapita_squareI_refhighdamage,
  totalimpacts_percapita_squareI_taxhighdamage,
  totalimpacts_percapita_squareF_refhighdamage,
  totalimpacts_percapita_squareF_taxhighdamage
)


# Compute percentage changes correctly
percentage_changes <- ((per_capita_impacts[seq(2, 14, by = 2)] - per_capita_impacts[seq(1, 13, by = 2)]) /
                         per_capita_impacts[seq(1, 13, by = 2)]) * 100

# Ensure `custom_labels2_str` is character, NOT expression (which ggplot may struggle with)
custom_labels2_str <- as.character(custom_labels2)

# Get absolute values of the already computed percentage changes
abs_percentage_changes <- abs(percentage_changes)

# Rank scenarios in descending order
rank_order <- order(abs_percentage_changes, decreasing = TRUE)

# Reorder labels and impact values while keeping Baseline & Pigovian tax pairs together
custom_labels2_sorted <- custom_labels2[rank_order]
sorted_indices <- as.vector(rbind(rank_order * 2 - 1, rank_order * 2))  # Keep pairs together

# Reorder the data based on sorted indices
data_sorted <- data.frame(
  Group = factor(rep(custom_labels2_sorted, each = 2), levels = custom_labels2_sorted),
  Category = rep(c("Baseline", "Pigovian tax"), 7),
  Numbers = as.numeric(per_capita_impacts[sorted_indices])
)

# Create labels for smallest bars
percentage_labels <- data_sorted %>%
  group_by(Group) %>%
  mutate(Numbers = as.numeric(Numbers)) %>%  # Ensure Numbers is numeric
  summarise(
    Category = Category[which.min(Numbers)],  # Find which bar is smaller
    Numbers = min(Numbers),  # Get the smallest bar's value
    label = sprintf("%.1f%%", percentage_changes[which(custom_labels2_str == first(Group))])  # Assign correct percentage
  ) %>%
  ungroup()  # Ungroup to prevent ggplot errors

# Debug check: Print structure of labels
str(percentage_labels)

# Create the bar plot
bar_plot <- ggplot(data_sorted, aes(x = Group, y = Numbers, fill = Category)) +  
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  # Side-by-side bars
  geom_text(data = percentage_labels, aes(x = Group, y = Numbers, label = label, group = Category),  
            position = position_dodge(width = 0.8), vjust = -0.5, color = "black", size = 7) +  # Labels above both bars
  geom_text(data = percentage_labels, aes(x = Group, y = Numbers, label = label, group = Category),  
            position = position_dodge(width = 0.8), vjust = -0.5, color = "black", size = 7) +  # Duplicate to ensure both bars
  labs(
    title = "Comparison of tax effects on impacts (sorted by absolute % change)",
    x = "Scenario (sorted)",
    y = "Per capita impacts"
  ) + 
  scale_fill_manual(values = c("Baseline" = "brown", "Pigovian tax" = "green")) +  # Categorical colors
  scale_x_discrete(labels = custom_labels2) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, margin = margin(t = 5, r = 10, b = 0, l = 10), hjust = 0.5),
    panel.grid.major.x = element_blank()
  )

print(bar_plot)

# Load and create heatmap plots for specific scenarios (excluding "Uniform" case)
heatmap_files <- c("heatmap scen E.png", "heatmap scen H.png", "heatmap scen G.png",
                   "heatmap scen A.png", "heatmap scen I.png", "heatmap scen F.png")

heatmap_plots <- list()

# Generate ggplot objects for heatmaps
for (i in 1:length(heatmap_files)) {
  heatmap_image <- readPNG(heatmap_files[i])
  heatmap_grob <- rasterGrob(as.raster(heatmap_image), interpolate = TRUE)
  
  # Create a ggplot for the heatmap
  heatmap_plot <- ggplot() + 
    annotation_custom(heatmap_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
    theme_void()
  
  heatmap_plots[[i]] <- heatmap_plot
}

# Adjust heatmap plots so they don't interfere with scales
heatmap_plots_fixed <- lapply(heatmap_plots, function(p) {
  p + scale_y_continuous(expand = c(0, 0))  # Prevents conflicts with y-scales
})

# Create a corrected heatmap grid with a blank space for "Uniform" (now in 4th position)
heatmap_grid_fixed <- wrap_plots(
  heatmap_plots_fixed[[1]], heatmap_plots_fixed[[2]], heatmap_plots_fixed[[3]],
  plot_spacer(),  # Empty space for "Uniform"
  heatmap_plots_fixed[[4]], heatmap_plots_fixed[[5]], heatmap_plots_fixed[[6]],
  ncol = 7
)

# Combine bar plot and heatmap grid correctly
final_plot <- bar_plot / heatmap_grid_fixed + plot_layout(heights = c(1, 0.25))

# Display the final plot
print(final_plot)



#Compute impact variations with tax between scenarios
ref <-  totalimpacts_percapita_squareH_refhighdamage
taxcase <- totalimpacts_percapita_squareH_taxhighdamage
tax_var = ((taxcase-ref)/ref)*100



##Test (not in the paper): Figure combining the plots and the legend heatmaps

# Define data (manually setting market shares for each scenario)
scenarios <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")
marketshare_base_A <- 1.0 + 3.8
marketshare_base_B <- 30.8 + 33.3
marketshare_base_C <- 94.2 + 4.5
marketshare_base_D <- 1.9 + 2.3
marketshare_base_E <- 40.1 + 18.6
marketshare_base_F <- 95.6 + 2.4
marketshare_base_G <- 1.3 + 0.2
marketshare_base_H <- 20.2 + 1.6
marketshare_base_I <- 52.3 + 0.5

green_marketshares <- c(
  marketshare_A, marketshare_B, marketshare_C, 
  marketshare_D, marketshare_E, marketshare_F, 
  marketshare_G, marketshare_H, marketshare_I
)

per_capita_impacts <- c(
  totalimpacts_percapita_squareA_reflowdamage,
  totalimpacts_percapita_squareB_reflowdamage,
  totalimpacts_percapita_squareC_reflowdamage,
  totalimpacts_percapita_squareD_reflowdamage,
  totalimpacts_percapita_squareE_reflowdamage,
  totalimpacts_percapita_squareF_reflowdamage,
  totalimpacts_percapita_squareG_reflowdamage,
  totalimpacts_percapita_squareH_reflowdamage,
  totalimpacts_percapita_squareI_reflowdamage
)

# Create data frame
data <- data.frame(
  Scenario = factor(scenarios, levels = scenarios), 
  MarketShare = green_marketshares, 
  Impact = per_capita_impacts
)

# Sort data by market share (ascending)
data <- data[order(data$MarketShare),]

# Define color gradient for points
custom_gradient_colors <- colorRampPalette(c("green", "tan", "saddlebrown"))(nrow(data))

# Create scatter plot with a single line connecting all points
scatter_plot <- ggplot(data, aes(x = MarketShare, y = Impact)) +
  geom_line(color = "lightblue", size = 0.8) +  # Thinner light blue line
  geom_point(aes(color = Impact), size = 3) +  # Smaller points
  scale_color_gradientn(colors = custom_gradient_colors) +  # Gradient for points
  labs(
    title = "Average Impacts vs. Green Market Shares in Different Scenarios (by Ascending Order of Green Shares), Low damage baseline",
    x = "Cumulative Market Share of Green Lifestyles (%)",
    y = "Average Impacts",
    color = "Impact Level"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

print(scatter_plot)

# Load heatmap images
heatmap_files <- c("heatmap scen G.png", "heatmap scen D.png", "heatmap scen A.png",
                   "heatmap scen H.png", "heatmap scen I.png", "heatmap scen E.png",
                   "heatmap scen B.png", "heatmap scen F.png", "heatmap scen C.png")

# Create heatmap plots for each scenario
heatmap_plots <- list()
for (i in 1:length(heatmap_files)) {
  heatmap_image <- readPNG(heatmap_files[i])  # Read the heatmap image
  heatmap_grob <- rasterGrob(as.raster(heatmap_image), interpolate = TRUE)  # Create a rasterGrob
  
  # Create a ggplot for the heatmap
  heatmap_plot <- ggplot() + 
    annotation_custom(heatmap_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
    theme_void()  # No axes or grid lines
  
  # Store the heatmap plot in the list
  heatmap_plots[[i]] <- heatmap_plot
}

# Combine scatter plot and heatmap plots using patchwork
heatmap_grid <- wrap_plots(heatmap_plots, ncol = length(heatmap_plots))

# Combine scatter plot and heatmap grid
final_plot <- scatter_plot / heatmap_grid + plot_layout(heights = c(1, 0.5))

# Display the final plot
print(final_plot)





### Histograms comparing scenarios for different cases

## 1. Low damage baseline

# Data for bar plot
categories <- c("A", "B","C","D","E","F","G","H","I")
per_capita_impacts <- c(
  totalimpacts_percapita_squareA_reflowdamage,
  totalimpacts_percapita_squareB_reflowdamage,
  totalimpacts_percapita_squareC_reflowdamage,
  totalimpacts_percapita_squareD_reflowdamage,
  totalimpacts_percapita_squareE_reflowdamage,
  totalimpacts_percapita_squareF_reflowdamage,
  totalimpacts_percapita_squareG_reflowdamage,
  totalimpacts_percapita_squareH_reflowdamage,
  totalimpacts_percapita_squareI_reflowdamage
)

# Create a data frame and order it by Numbers
data <- data.frame(
  Category = factor(categories, levels = categories),
  Numbers = per_capita_impacts
)
data <- data[order(data$Numbers, decreasing = TRUE), ]  # Order by impacts

# Set custom colors for the gradient fill
custom_gradient_colors <- colorRampPalette(c("green", "tan", "saddlebrown"))(length(per_capita_impacts))

# Create a custom legend for each bar with Greek symbols
custom_labels <- c(
  expression("("*alpha[l]*", "*beta[m]*")"),
  expression("("*alpha[l]*", "*beta[h]*")"),
  expression("("*alpha[l]*", "*beta[l]*")"),
  expression("("*alpha[m]*", "*beta[l]*")"),
  expression("("*alpha[m]*", "*beta[m]*")"),
  expression("("*alpha[m]*", "*beta[h]*")"),
  expression("("*alpha[h]*", "*beta[l]*")"),
  expression("("*alpha[h]*", "*beta[h]*")"),
  expression("("*alpha[h]*", "*beta[m]*")")
)

# Define the baseline impact value
baseline_impact <- totalimpacts_percapita_reflowdamage_unif   # Replace with actual value

# Find the approximate position where the baseline should be inserted
data$Rank <- rank(-data$Numbers, ties.method = "first")  # Rank in decreasing order
insert_pos <- sum(data$Numbers > baseline_impact) + 0.5  # Position between bars

# Create the bar plot with a custom legend and vertical line
bar_plot <- ggplot(data, aes(x = reorder(Category, -Numbers), y = Numbers, fill = Numbers)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_gradientn(colors = custom_gradient_colors) +  # Gradient fill applied here
  labs(
    title = "Average impacts in different population scenarios - Low damage baseline",
    x = "Population concentration scenario (environmental axis, image axis) - l:low, m:medium, h:high",
    y = "Average impacts"
  ) +
  scale_x_discrete(labels = custom_labels) +  # Use custom labels for the x-axis
  geom_segment(aes(x = insert_pos, xend = insert_pos, y = 0, yend = baseline_impact), 
               linetype = "dotted", color = "blue", linewidth = 1) +  # Vertical line with proper height
  annotate("text", x = insert_pos, y = baseline_impact * 1.05, label = "Uniform", 
           color = "blue", size = 3, fontface = "bold", hjust = 0.5) +  # Legend above the line
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 15)  # Adjust text size and angle for better readability
  )

print(bar_plot)

# Load heatmap images (put them manually in the right order, or change the code!)
heatmap_files <- c("heatmap scen D.png", "heatmap scen A.png", "heatmap scen G.png",
                   "heatmap scen H.png", "heatmap scen E.png", "heatmap scen B.png",
                   "heatmap scen I.png", "heatmap scen C.png", "heatmap scen F.png")


# Create heatmap plots for each scenario
for (i in 1:length(heatmap_files)) {
  heatmap_image <- readPNG(heatmap_files[i])  # Read the heatmap image
  heatmap_grob <- rasterGrob(as.raster(heatmap_image), interpolate = TRUE)  # Create a rasterGrob
  
  # Create a ggplot for the heatmap
  heatmap_plot <- ggplot() + 
    annotation_custom(heatmap_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
    theme_void()  # No axes or grid lines
  
  # Store the heatmap plot in the list
  heatmap_plots[[i]] <- heatmap_plot
}

# Combine bar plot and heatmap plots using patchwork
# Create a plot grid with reduced spacing
heatmap_grid <- wrap_plots(heatmap_plots, ncol = length(heatmap_plots))

# Combine bar plot and heatmap grid with less spacing
final_plot <- bar_plot / heatmap_grid + plot_layout(heights = c(1, 0.5))  # Adjust heights as needed

# Display the combined plot
print(final_plot)



##2. Do the same with R=24 (50% higher income)

# Data for bar plot
categories <- c("A", "B","C","D","E","F","G","H","I")
per_capita_impacts <- c(
  totalimpacts_percapita_squareA_higherincomelowdamage,
  totalimpacts_percapita_squareB_higherincomelowdamage,
  totalimpacts_percapita_squareC_higherincomelowdamage,
  totalimpacts_percapita_squareD_higherincomelowdamage,
  totalimpacts_percapita_squareE_higherincomelowdamage,
  totalimpacts_percapita_squareF_higherincomelowdamage,
  totalimpacts_percapita_squareG_higherincomelowdamage,
  totalimpacts_percapita_squareH_higherincomelowdamage,
  totalimpacts_percapita_squareI_higherincomelowdamage
)

# Create a data frame and order it by Numbers
data <- data.frame(
  Category = factor(categories, levels = categories),
  Numbers = per_capita_impacts
)
data <- data[order(data$Numbers, decreasing = TRUE), ]  # Order by impacts

# Set custom colors for the gradient fill
custom_gradient_colors <- colorRampPalette(c("green", "tan", "saddlebrown"))(length(per_capita_impacts))

# Create a custom legend for each bar with Greek symbols
custom_labels <- c(
  expression("("*alpha[l]*", "*beta[m]*")"),
  expression("("*alpha[l]*", "*beta[h]*")"),
  expression("("*alpha[l]*", "*beta[l]*")"),
  expression("("*alpha[m]*", "*beta[l]*")"),
  expression("("*alpha[m]*", "*beta[m]*")"),
  expression("("*alpha[m]*", "*beta[h]*")"),
  expression("("*alpha[h]*", "*beta[l]*")"),
  expression("("*alpha[h]*", "*beta[h]*")"),
  expression("("*alpha[h]*", "*beta[m]*")")
)

# Define the baseline impact value
baseline_impact <- totalimpacts_percapita_higherincomelowdamage_unif   # Replace with actual value

# Find the approximate position where the baseline should be inserted
data$Rank <- rank(-data$Numbers, ties.method = "first")  # Rank in decreasing order
insert_pos <- sum(data$Numbers > baseline_impact) + 0.5  # Position between bars

# Create the bar plot with a custom legend and vertical line
bar_plot <- ggplot(data, aes(x = reorder(Category, -Numbers), y = Numbers, fill = Numbers)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_gradientn(colors = custom_gradient_colors) +  # Gradient fill applied here
  labs(
    title = "Average impacts in different population scenarios - R=24, low damage",
    x = "Population concentration scenario (environmental axis, image axis) - l:low, m:medium, h:high",
    y = "Average impacts"
  ) +
  scale_x_discrete(labels = custom_labels) +  # Use custom labels for the x-axis
  geom_segment(aes(x = insert_pos, xend = insert_pos, y = 0, yend = baseline_impact), 
               linetype = "dotted", color = "blue", linewidth = 1) +  # Vertical line with proper height
  annotate("text", x = insert_pos, y = baseline_impact * 1.05, label = "Uniform", 
           color = "blue", size = 3, fontface = "bold", hjust = 0.5) +  # Legend above the line
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 15)  # Adjust text size and angle for better readability
  )

print(bar_plot)

#Load heatmap images: 
#Looking at our impacts vector, we can see that the order remains the same so we don't have to load them again

# Combine bar plot and heatmap grid with less spacing
final_plot <- bar_plot / heatmap_grid + plot_layout(heights = c(1, 0.5))  # Adjust heights as needed

# Display the combined plot
print(final_plot)






##3. Same for high damage baseline

# Data for bar plot
categories <- c("A", "B","C","D","E","F","G","H","I")
per_capita_impacts <- c(
  totalimpacts_percapita_squareA_refhighdamage,
  totalimpacts_percapita_squareB_refhighdamage,
  totalimpacts_percapita_squareC_refhighdamage,
  totalimpacts_percapita_squareD_refhighdamage,
  totalimpacts_percapita_squareE_refhighdamage,
  totalimpacts_percapita_squareF_refhighdamage,
  totalimpacts_percapita_squareG_refhighdamage,
  totalimpacts_percapita_squareH_refhighdamage,
  totalimpacts_percapita_squareI_refhighdamage
)

# Create a data frame and order it by Numbers
data <- data.frame(
  Category = factor(categories, levels = categories),
  Numbers = per_capita_impacts
)
data <- data[order(data$Numbers, decreasing = TRUE), ]  # Order by impacts

# Set custom colors for the gradient fill
custom_gradient_colors <- colorRampPalette(c("green", "tan", "saddlebrown"))(length(per_capita_impacts))

# Create a custom legend for each bar with Greek symbols
custom_labels <- c(
  expression("("*alpha[l]*", "*beta[m]*")"),
  expression("("*alpha[l]*", "*beta[h]*")"),
  expression("("*alpha[l]*", "*beta[l]*")"),
  expression("("*alpha[m]*", "*beta[l]*")"),
  expression("("*alpha[m]*", "*beta[h]*")"),
  expression("("*alpha[m]*", "*beta[m]*")"),
  expression("("*alpha[h]*", "*beta[l]*")"),
  expression("("*alpha[h]*", "*beta[h]*")"),
  expression("("*alpha[h]*", "*beta[m]*")")
)

# Define the baseline impact value
baseline_impact <- totalimpacts_percapita_refhighdamage_unif   # Replace with actual value

# Find the approximate position where the baseline should be inserted
data$Rank <- rank(-data$Numbers, ties.method = "first")  # Rank in decreasing order
insert_pos <- sum(data$Numbers > baseline_impact) + 0.5  # Position between bars

# Create the bar plot with a custom legend and vertical line
bar_plot <- ggplot(data, aes(x = reorder(Category, -Numbers), y = Numbers, fill = Numbers)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_gradientn(colors = custom_gradient_colors) +  # Gradient fill applied here
  labs(
    title = "Average impacts in different population scenarios - High damage baseline",
    x = "Population concentration scenario (environmental axis, image axis) - l:low, m:medium, h:high",
    y = "Average impacts"
  ) +
  scale_x_discrete(labels = custom_labels) +  # Use custom labels for the x-axis
  geom_segment(aes(x = insert_pos, xend = insert_pos, y = 0, yend = baseline_impact), 
               linetype = "dotted", color = "blue", linewidth = 1) +  # Vertical line with proper height
  annotate("text", x = insert_pos, y = baseline_impact * 1.05, label = "Uniform", 
           color = "blue", size = 3, fontface = "bold", hjust = 0.5) +  # Legend above the line
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 15)  # Adjust text size and angle for better readability
  )

print(bar_plot)

# Load heatmap images : here the order changes a bit so we reload the heatmaps
heatmap_files_2 <- c("heatmap scen D.png", "heatmap scen A.png", "heatmap scen G.png",
                   "heatmap scen H.png", "heatmap scen B.png", "heatmap scen E.png",
                   "heatmap scen I.png", "heatmap scen C.png", "heatmap scen F.png")

# Create heatmap plots for each scenario
for (i in 1:length(heatmap_files_2)) {
  heatmap_image_2 <- readPNG(heatmap_files_2[i])  # Read the heatmap image
  heatmap_grob_2 <- rasterGrob(as.raster(heatmap_image_2), interpolate = TRUE)  # Create a rasterGrob
  
  # Create a ggplot for the heatmap
  heatmap_plot_2 <- ggplot() + 
    annotation_custom(heatmap_grob_2, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
    theme_void()  # No axes or grid lines
  
  # Store the heatmap plot in the list
  heatmap_plots_2[[i]] <- heatmap_plot_2
}

# Combine bar plot and heatmap plots using patchwork
# Create a plot grid with reduced spacing
heatmap_grid_2 <- wrap_plots(heatmap_plots_2, ncol = length(heatmap_plots_2))

# Combine bar plot and heatmap grid with less spacing
final_plot <- bar_plot / heatmap_grid_2 + plot_layout(heights = c(1, 0.5))  # Adjust heights as needed

# Display the combined plot
print(final_plot)




##Same with other scenarios (here lower dispersion of preferences)

# Data for bar plot
categories <- c("A", "B","C","D","E","F","G","H","I")
per_capita_impacts <- c(
  totalimpacts_percapita_squareA2_reflowdamage,
  totalimpacts_percapita_squareB2_reflowdamage,
  totalimpacts_percapita_squareC2_reflowdamage,
  totalimpacts_percapita_squareD2_reflowdamage,
  totalimpacts_percapita_squareE2_reflowdamage,
  totalimpacts_percapita_squareF2_reflowdamage,
  totalimpacts_percapita_squareG2_reflowdamage,
  totalimpacts_percapita_squareH2_reflowdamage,
  totalimpacts_percapita_squareI2_reflowdamage
)

# Create a data frame and order it by Numbers
data <- data.frame(
  Category = factor(categories, levels = categories),
  Numbers = per_capita_impacts
)
data <- data[order(data$Numbers, decreasing = TRUE), ]  # Order by impacts

# Set custom colors for the gradient fill
custom_gradient_colors <- colorRampPalette(c("green", "tan", "saddlebrown"))(length(per_capita_impacts))

# Create a custom legend for each bar with Greek symbols
custom_labels <- c(
  expression("("*alpha[l]*", "*beta[m]*")"),
  expression("("*alpha[l]*", "*beta[h]*")"),
  expression("("*alpha[l]*", "*beta[l]*")"),
  expression("("*alpha[m]*", "*beta[l]*")"),
  expression("("*alpha[m]*", "*beta[m]*")"),
  expression("("*alpha[m]*", "*beta[h]*")"),
  expression("("*alpha[h]*", "*beta[l]*")"),
  expression("("*alpha[h]*", "*beta[h]*")"),
  expression("("*alpha[h]*", "*beta[m]*")")
)

# Define the baseline impact value
baseline_impact <- totalimpacts_percapita_reflowdamage_unif   # Replace with actual value

# Find the approximate position where the baseline should be inserted
data$Rank <- rank(-data$Numbers, ties.method = "first")  # Rank in decreasing order
insert_pos <- sum(data$Numbers > baseline_impact) + 0.5  # Position between bars

# Create the bar plot with a custom legend and vertical line
bar_plot <- ggplot(data, aes(x = reorder(Category, -Numbers), y = Numbers, fill = Numbers)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_gradientn(colors = custom_gradient_colors) +  # Gradient fill applied here
  labs(
    title = "Average impacts in scenarios with less preference dispersion - Low damage baseline",
    x = "Population concentration scenario (environmental axis, image axis) - l:low, m:medium, h:high",
    y = "Average impacts"
  ) +
  scale_x_discrete(labels = custom_labels) +  # Use custom labels for the x-axis
  geom_segment(aes(x = insert_pos, xend = insert_pos, y = 0, yend = baseline_impact), 
               linetype = "dotted", color = "blue", linewidth = 1) +  # Vertical line with proper height
  annotate("text", x = insert_pos, y = baseline_impact * 1.05, label = "Uniform", 
           color = "blue", size = 3, fontface = "bold", hjust = 0.5) +  # Legend above the line
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 15)  # Adjust text size and angle for better readability
  )

print(bar_plot)


# Combine bar plot and heatmap grid (same order) with less spacing
final_plot <- bar_plot / heatmap_grid + plot_layout(heights = c(1, 0.5))  # Adjust heights as needed

# Display the combined plot
print(final_plot)





##Same under Assumption 4B

# Data for bar plot
categories <- c("A", "B","C","D","E","F","G","H","I")
per_capita_impacts_ref4B <- c(
  totalimpacts_percapita_squareA_discrmorepollhighdamage,
  totalimpacts_percapita_squareB_discrmorepollhighdamage,
  totalimpacts_percapita_squareC_discrmorepollhighdamage,
  totalimpacts_percapita_squareD_discrmorepollhighdamage,
  totalimpacts_percapita_squareE_discrmorepollhighdamage,
  totalimpacts_percapita_squareF_discrmorepollhighdamage,
  totalimpacts_percapita_squareG_discrmorepollhighdamage,
  totalimpacts_percapita_squareH_discrmorepollhighdamage,
  totalimpacts_percapita_squareI_discrmorepollhighdamage
)

# Create a data frame and order it by Numbers
data <- data.frame(
  Category = factor(categories, levels = categories),
  Numbers = per_capita_impacts_ref4B
)
data <- data[order(data$Numbers, decreasing = TRUE), ]  # Order by impacts

# Set custom colors for the gradient fill
custom_gradient_colors <- colorRampPalette(c("green", "tan", "saddlebrown"))(length(per_capita_impacts))

# Create a custom legend for each bar with Greek symbols
custom_labels <- c(
  expression("("*alpha[l]*", "*beta[m]*")"),
  expression("("*alpha[l]*", "*beta[l]*")"),
  expression("("*alpha[l]*", "*beta[h]*")"),
  expression("("*alpha[m]*", "*beta[l]*")"),
  expression("("*alpha[m]*", "*beta[m]*")"),
  expression("("*alpha[m]*", "*beta[h]*")"),
  expression("("*alpha[h]*", "*beta[l]*")"),
  expression("("*alpha[h]*", "*beta[h]*")"),
  expression("("*alpha[h]*", "*beta[m]*")")
)

# Define the baseline impact value
baseline_impact <- totalimpacts_percapita_discrmorepollhighdamage_unif   # Replace with actual value

# Find the approximate position where the baseline should be inserted
data$Rank <- rank(-data$Numbers, ties.method = "first")  # Rank in decreasing order
insert_pos <- sum(data$Numbers > baseline_impact) + 0.5  # Position between bars

# Create the bar plot with a custom legend and vertical line
bar_plot_4B <- ggplot(data, aes(x = reorder(Category, -Numbers), y = Numbers, fill = Numbers)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_gradientn(colors = custom_gradient_colors) +  # Gradient fill applied here
  labs(
    title = "Average impacts in different population scenarios under Assumption 4B - High damage baseline",
    x = "Population concentration scenario (environmental axis, image axis) - l:low, m:medium, h:high",
    y = "Average impacts"
  ) +
  scale_x_discrete(labels = custom_labels) +  # Use custom labels for the x-axis
  geom_segment(aes(x = insert_pos, xend = insert_pos, y = 0, yend = baseline_impact), 
               linetype = "dotted", color = "blue", linewidth = 1) +  # Vertical line with proper height
  annotate("text", x = insert_pos, y = baseline_impact * 1.05, label = "Uniform", 
           color = "blue", size = 3, fontface = "bold", hjust = 0.5) +  # Legend above the line
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 15)  # Adjust text size and angle for better readability
  )

print(bar_plot_4B)

# Load heatmap images: the order changes again
heatmap_files_3 <- c("heatmap scen D.png", "heatmap scen G.png", "heatmap scen A.png",
                   "heatmap scen H.png", "heatmap scen E.png", "heatmap scen B.png",
                   "heatmap scen I.png", "heatmap scen C.png", "heatmap scen F.png")

# Create heatmap plots for each scenario
for (i in 1:length(heatmap_files_3)) {
  heatmap_image_3 <- readPNG(heatmap_files_3[i])  # Read the heatmap image
  heatmap_grob_3 <- rasterGrob(as.raster(heatmap_image_3), interpolate = TRUE)  # Create a rasterGrob
  
  # Create a ggplot for the heatmap
  heatmap_plot_3 <- ggplot() + 
    annotation_custom(heatmap_grob_3, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
    theme_void()  # No axes or grid lines
  
  # Store the heatmap plot in the list
  heatmap_plots_3[[i]] <- heatmap_plot_3
}

# Combine bar plot and heatmap plots using patchwork
# Create a plot grid with reduced spacing
heatmap_grid_3 <- wrap_plots(heatmap_plots_3, ncol = length(heatmap_plots_3))

# Combine bar plot and heatmap grid with less spacing
final_plot_4B <- bar_plot_4B / heatmap_grid_3 + plot_layout(heights = c(1, 0.5))  # Adjust heights as needed

# Display the combined plot
print(final_plot_4B)



##Same with R=48

# Data for bar plot
categories <- c("A", "B","C","D","E","F","G","H","I")
per_capita_impacts <- c(
  totalimpacts_percapita_squareA_hugeincome_lowdamage,
  totalimpacts_percapita_squareB_hugeincome_lowdamage,
  totalimpacts_percapita_squareC_hugeincome_lowdamage,
  totalimpacts_percapita_squareD_hugeincome_lowdamage,
  totalimpacts_percapita_squareE_hugeincome_lowdamage,
  totalimpacts_percapita_squareF_hugeincome_lowdamage,
  totalimpacts_percapita_squareG_hugeincome_lowdamage,
  totalimpacts_percapita_squareH_hugeincome_lowdamage,
  totalimpacts_percapita_squareI_hugeincome_lowdamage
)

# Create a data frame and order it by Numbers
data <- data.frame(
  Category = factor(categories, levels = categories),
  Numbers = per_capita_impacts
)
data <- data[order(data$Numbers, decreasing = TRUE), ]  # Order by impacts

# Set custom colors for the gradient fill
custom_gradient_colors <- colorRampPalette(c("green", "tan", "saddlebrown"))(length(per_capita_impacts))

# Create a custom legend for each bar with Greek symbols
custom_labels <- c(
  expression("("*alpha[l]*", "*beta[h]*")"),
  expression("("*alpha[l]*", "*beta[m]*")"),
  expression("("*alpha[l]*", "*beta[l]*")"),
  expression("("*alpha[m]*", "*beta[h]*")"),
  expression("("*alpha[h]*", "*beta[m]*")"),
  expression("("*alpha[m]*", "*beta[m]*")"),
  expression("("*alpha[m]*", "*beta[l]*")"),
  expression("("*alpha[h]*", "*beta[l]*")"),
  expression("("*alpha[h]*", "*beta[h]*")")
)

# Create the bar plot with a custom legend
library(ggplot2)
bar_plot <- ggplot(data, aes(x = reorder(Category, -Numbers), y = Numbers, fill = Numbers)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_gradientn(colors = custom_gradient_colors) +  # Gradient fill applied here
  labs(
    title = "Average impacts in different population scenarios - R=48, Low damage",
    x = "Population concentration scenario (environmental axis, image axis) - l:low, m:medium, h:high",
    y = "Average impacts"
  ) +
  scale_x_discrete(labels = custom_labels) +  # Use custom labels for the x-axis
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 15)  # Adjust text size and angle for better readability
  )


# Display the bar plot
print(bar_plot)


# Create heatmap plots for each scenario (same order as high damage baseline)

# Combine bar plot and heatmap grid with less spacing
final_plot <- bar_plot / heatmap_grid_2 + plot_layout(heights = c(1, 0.5))  # Adjust heights as needed

# Display the combined plot
print(final_plot)




##Sensitivity analysis (alternative scenarios), higher damage reference case.

#More extreme scenarios

categories <- c("Uniform", "A2", "B2", "C2", "D2", "E2", "F2", "G2", "H2", "I2")

per_capita_impacts <- c(
  totalimpacts_percapita_refhighdamage_unif,
  totalimpacts_percapita_squareA2_refhighdamage,
  totalimpacts_percapita_squareB2_refhighdamage,
  totalimpacts_percapita_squareC2_refhighdamage,
  totalimpacts_percapita_squareD2_refhighdamage,
  totalimpacts_percapita_squareE2_refhighdamage,
  totalimpacts_percapita_squareF2_refhighdamage,
  totalimpacts_percapita_squareG2_refhighdamage,
  totalimpacts_percapita_squareH2_refhighdamage,
  totalimpacts_percapita_squareI2_refhighdamage
)

# Create a data frame
data <- data.frame(
  Category = factor(categories, levels = categories),
  Numbers = per_capita_impacts
)

# Reorder the Category factor based on descending Numbers
data <- data[order(data$Numbers, decreasing = TRUE), ]
data$Category <- factor(data$Category, levels = data$Category)

# Set custom colors for the gradient fill
custom_gradient_colors <- colorRampPalette(c("green", "tan", "saddlebrown"))(length(per_capita_impacts))

# Define custom labels with Greek letters and subscripts for each pair
custom_labels <- c(
  expression("("* alpha[l] * "," * beta[m] * ")" * scriptstyle("2")),
  expression("("* alpha[l] * "," * beta[h] * ")" * scriptstyle("2")),
  expression("("* alpha[l] * "," * beta[l] * ")" * scriptstyle("2")),
  expression("Uniform"),
  expression("("* alpha[m] * "," * beta[h] * ")" * scriptstyle("2")),
  expression("("* alpha[m] * "," * beta[m] * ")" * scriptstyle("2")),
  expression("("* alpha[m] * "," * beta[l] * ")" * scriptstyle("2")),
  expression("("* alpha[h] * "," * beta[l] * ")" * scriptstyle("2")),
  expression("("* alpha[h] * "," * beta[h] * ")" * scriptstyle("2")),
  expression("("* alpha[h] * "," * beta[m] * ")" * scriptstyle("2"))
)


# Create the bar plot using ggplot2
plot_pcimpact_refhighd_sens <- ggplot(data, aes(x = Category, y = Numbers, fill = Numbers)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(
    title = "Per capita impacts with more 'extreme' population scenarios, high damage baseline",
    x = "Population concentration scenario (environmental axis, image axis)",
    y = "Per capita impacts"
  ) +
  scale_fill_gradientn(
    colors = custom_gradient_colors,
    guide = "legend",
    name = "Impacts"
  ) +
  scale_x_discrete(labels = custom_labels) +  # Use custom Greek labels for x-axis
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.75)
  )

# Display the plot
print(plot_pcimpact_refhighd_sens)



#Less extreme scenarios
categories <- c("Uniform", "A3", "B3", "C3", "D3", "E3", "F3", "G3", "H3", "I3")

per_capita_impacts <- c(
  totalimpacts_percapita_refhighdamage_unif,
  totalimpacts_percapita_squareA3_refhighdamage,
  totalimpacts_percapita_squareB3_refhighdamage,
  totalimpacts_percapita_squareC3_refhighdamage,
  totalimpacts_percapita_squareD3_refhighdamage,
  totalimpacts_percapita_squareE3_refhighdamage,
  totalimpacts_percapita_squareF3_refhighdamage,
  totalimpacts_percapita_squareG3_refhighdamage,
  totalimpacts_percapita_squareH3_refhighdamage,
  totalimpacts_percapita_squareI3_refhighdamage
)

# Create a data frame
data <- data.frame(
  Category = factor(categories, levels = categories),
  Numbers = per_capita_impacts
)

# Reorder the Category factor based on descending Numbers
data <- data[order(data$Numbers, decreasing = TRUE), ]
data$Category <- factor(data$Category, levels = data$Category)

# Set custom colors for the gradient fill
custom_gradient_colors <- colorRampPalette(c("green", "tan", "saddlebrown"))(length(per_capita_impacts))

# Define custom labels with Greek letters and subscripts for each pair
custom_labels <- c(
  expression("("* alpha[l] * "," * beta[m] * ")" * scriptstyle("3")),
  expression("("* alpha[l] * "," * beta[h] * ")" * scriptstyle("3")),
  expression("("* alpha[l] * "," * beta[l] * ")" * scriptstyle("3")),
  expression("Uniform"),
  expression("("* alpha[m] * "," * beta[m] * ")" * scriptstyle("3")),
  expression("("* alpha[m] * "," * beta[h] * ")" * scriptstyle("3")),
  expression("("* alpha[h] * "," * beta[h] * ")" * scriptstyle("3")),
  expression("("* alpha[h] * "," * beta[m] * ")" * scriptstyle("3")),
  expression("("* alpha[m] * "," * beta[l] * ")" * scriptstyle("3")),
  expression("("* alpha[h] * "," * beta[l] * ")" * scriptstyle("3"))
)


# Create the bar plot using ggplot2
plot_pcimpact_refhighd_sens2 <- ggplot(data, aes(x = Category, y = Numbers, fill = Numbers)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(
    title = "Per capita impacts with less 'extreme' population scenarios, high damage baseline",
    x = "Population concentration scenario (environmental axis, image axis)",
    y = "Per capita impacts"
  ) +
  scale_fill_gradientn(
    colors = custom_gradient_colors,
    guide = "legend",
    name = "Impacts"
  ) +
  scale_x_discrete(labels = custom_labels) +  # Use custom Greek labels for x-axis
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.75)
  )

# Display the plot
print(plot_pcimpact_refhighd_sens2)







### Figure to compare tax and nudges by scenarios in the low damage baseline

# Define scenarios
scenarios <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")

# Define impacts for each scenario (Baseline, Tax, Low Income)
impacts_baseline <- c(
  totalimpacts_percapita_squareA_reflowdamage,
  totalimpacts_percapita_squareB_reflowdamage,
  totalimpacts_percapita_squareC_reflowdamage,
  totalimpacts_percapita_squareD_reflowdamage,
  totalimpacts_percapita_squareE_reflowdamage,
  totalimpacts_percapita_squareF_reflowdamage,
  totalimpacts_percapita_squareG_reflowdamage,
  totalimpacts_percapita_squareH_reflowdamage,
  totalimpacts_percapita_squareI_reflowdamage
)

impacts_tax <- c(
  totalimpacts_percapita_squareA_taxlowdamage,
  totalimpacts_percapita_squareB_taxlowdamage,
  totalimpacts_percapita_squareC_taxlowdamage,
  totalimpacts_percapita_squareD_taxlowdamage,
  totalimpacts_percapita_squareE_taxlowdamage,
  totalimpacts_percapita_squareF_taxlowdamage,
  totalimpacts_percapita_squareG_taxlowdamage,
  totalimpacts_percapita_squareH_taxlowdamage,
  totalimpacts_percapita_squareI_taxlowdamage
)

impacts_low_income <- c(
  totalimpacts_percapita_squareA_lowerincomelowdamage,
  totalimpacts_percapita_squareB_lowerincomelowdamage,
  totalimpacts_percapita_squareC_lowerincomelowdamage,
  totalimpacts_percapita_squareD_lowerincomelowdamage,
  totalimpacts_percapita_squareE_lowerincomelowdamage,
  totalimpacts_percapita_squareF_lowerincomelowdamage,
  totalimpacts_percapita_squareG_lowerincomelowdamage,
  totalimpacts_percapita_squareH_lowerincomelowdamage,
  totalimpacts_percapita_squareI_lowerincomelowdamage
)

# Rank scenarios by descending order of baseline impacts
ordered_indices <- order(impacts_baseline, decreasing = TRUE)
ordered_scenarios <- scenarios[ordered_indices]
ordered_impacts_baseline <- impacts_baseline[ordered_indices]
ordered_impacts_tax <- impacts_tax[ordered_indices]
ordered_impacts_low_income <- impacts_low_income[ordered_indices]


# Create a data frame respecting the sorted order, aligning the impact values to the sorted scenarios
data <- data.frame(
  Scenario = ordered_scenarios,
  Impact_Baseline = ordered_impacts_baseline,
  Impact_Tax = ordered_impacts_tax,
  Impact_Low_Income = ordered_impacts_low_income
)

# Transform the data to long format for ggplot (so that each scenario's impacts are stacked)
data_long <- tidyr::pivot_longer(data, 
                                 cols = starts_with("Impact"), 
                                 names_to = "Case", 
                                 values_to = "Impact")

# Adjust the "Case" variable to represent the different scenarios' impacts
data_long$Case <- factor(data_long$Case, levels = c("Impact_Baseline", "Impact_Tax", "Impact_Low_Income"))

# Compute the minimum Impact for each Scenario to use as the yend for vertical lines
min_impacts <- data_long %>%
  group_by(Scenario) %>%
  summarize(min_impact = min(Impact))

# Merge min_impacts back into data_long to get the min_impact for each row
data_long <- left_join(data_long, min_impacts, by = "Scenario")


# Create the scatter plot with vertical lines stopping at the lowest point
scatter_plot <- ggplot(data_long, aes(x = Scenario, y = Impact, color = Impact, shape = Case)) +
  geom_point(size = 3) +  # Add points
  geom_segment(aes(xend = Scenario, yend = min_impact), color = "lightblue", size = 1) +  # Vertical lines stopping at the lowest impact
  scale_color_gradientn(colors = c("green", "tan", "saddlebrown")) +  # Color gradient for Impact values
  scale_shape_manual(values = c(16, 17, 15), name = "Case") +  # Different shapes for B, T, I
  labs(
    title = "Impacts by Cases in Different Population Scenarios - Low Damage Baseline",
    x = "Population Scenarios ranked by descending order of Baseline Impacts",
    y = "Average Impacts"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove text (letters) above x-axis
    axis.title.x = element_text(margin = margin(t = 10)),  # Adjust space between x-axis title and labels
    legend.position = "right",  # Move legend to the right
    axis.ticks.x = element_blank()  # Remove ticks on x-axis
  ) +
  scale_x_discrete(limits = ordered_scenarios)  # Order the x-axis according to the ordered_scenarios

# Display the plot
print(scatter_plot)


# Load heatmap images (agin the order is that of the  low damage baseline)

# Create heatmap plots for each scenario
heatmap_plots <- list()
for (i in seq_along(heatmap_files)) {
  heatmap_image <- readPNG(heatmap_files[i])  # Read the heatmap image
  heatmap_grob <- rasterGrob(as.raster(heatmap_image), interpolate = TRUE)  # Create a rasterGrob
  
  # Create a ggplot for the heatmap
  heatmap_plot <- ggplot() + 
    annotation_custom(heatmap_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
    theme_void()  # No axes or grid lines
  
  heatmap_plots[[i]] <- heatmap_plot
}

# Combine scatter plot and heatmap plots using patchwork
library(patchwork)
heatmap_grid <- wrap_plots(heatmap_plots, ncol = length(heatmap_plots))

# Combine scatter plot and heatmap grid
final_plot <- scatter_plot / heatmap_grid + plot_layout(heights = c(1, 0.5))

# Display the final plot
print(final_plot)




##Same for high damage levels

# Define scenarios
scenarios <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")

# Define impacts for each scenario (Baseline, Tax, Low Income)
impacts_baseline_highd <- c(
  totalimpacts_percapita_squareA_refhighdamage,
  totalimpacts_percapita_squareB_refhighdamage,
  totalimpacts_percapita_squareC_refhighdamage,
  totalimpacts_percapita_squareD_refhighdamage,
  totalimpacts_percapita_squareE_refhighdamage,
  totalimpacts_percapita_squareF_refhighdamage,
  totalimpacts_percapita_squareG_refhighdamage,
  totalimpacts_percapita_squareH_refhighdamage,
  totalimpacts_percapita_squareI_refhighdamage
)

impacts_tax_highd <- c(
  totalimpacts_percapita_squareA_taxhighdamage,
  totalimpacts_percapita_squareB_taxhighdamage,
  totalimpacts_percapita_squareC_taxhighdamage,
  totalimpacts_percapita_squareD_taxhighdamage,
  totalimpacts_percapita_squareE_taxhighdamage,
  totalimpacts_percapita_squareF_taxhighdamage,
  totalimpacts_percapita_squareG_taxhighdamage,
  totalimpacts_percapita_squareH_taxhighdamage,
  totalimpacts_percapita_squareI_taxhighdamage
)

impacts_low_income_highd <- c(
  totalimpacts_percapita_squareA_lowincomehighdamage,
  totalimpacts_percapita_squareB_lowincomehighdamage,
  totalimpacts_percapita_squareC_lowincomehighdamage,
  totalimpacts_percapita_squareD_lowincomehighdamage,
  totalimpacts_percapita_squareE_lowincomehighdamage,
  totalimpacts_percapita_squareF_lowincomehighdamage,
  totalimpacts_percapita_squareG_lowincomehighdamage,
  totalimpacts_percapita_squareH_lowincomehighdamage,
  totalimpacts_percapita_squareI_lowincomehighdamage
)

# Rank scenarios by descending order of baseline impacts
ordered_indices <- order(impacts_baseline_highd, decreasing = TRUE)
ordered_scenarios_highd <- scenarios[ordered_indices]
ordered_impacts_baseline_highd <- impacts_baseline_highd[ordered_indices]
ordered_impacts_tax_highd  <- impacts_tax_highd[ordered_indices]
ordered_impacts_low_income_highd  <- impacts_low_income_highd[ordered_indices]

# Create a data frame respecting the sorted order, aligning the impact values to the sorted scenarios
data <- data.frame(
  Scenario = ordered_scenarios_highd,
  Impact_Baseline = ordered_impacts_baseline_highd,
  Impact_Tax = ordered_impacts_tax_highd,
  Impact_Low_Income = ordered_impacts_low_income_highd 
)

# Transform the data to long format for ggplot (so that each scenario's impacts are stacked)
data_long <- tidyr::pivot_longer(data, 
                                 cols = starts_with("Impact"), 
                                 names_to = "Case", 
                                 values_to = "Impact")

# Adjust the "Case" variable to represent the different scenarios' impacts
data_long$Case <- factor(data_long$Case, levels = c("Impact_Baseline", "Impact_Tax", "Impact_Low_Income"))

# Compute the minimum Impact for each Scenario to use as the yend for vertical lines
min_impacts <- data_long %>%
  group_by(Scenario) %>%
  summarize(min_impact = min(Impact))

# Merge min_impacts back into data_long to get the min_impact for each row
data_long <- left_join(data_long, min_impacts, by = "Scenario")


# Create the scatter plot with vertical lines stopping at the lowest point
scatter_plot_highd <- ggplot(data_long, aes(x = Scenario, y = Impact, color = Impact, shape = Case)) +
  geom_point(size = 3) +  # Add points
  geom_segment(aes(xend = Scenario, yend = min_impact), color = "lightblue", size = 1) +  # Vertical lines stopping at the lowest impact
  scale_color_gradientn(colors = c("green", "tan", "saddlebrown")) +  # Color gradient for Impact values
  scale_shape_manual(values = c(16, 17, 15), name = "Case") +  # Different shapes for B, T, I
  labs(
    title = "Impacts by Cases in Different Population Scenarios - High Damage Baseline",
    x = "Population Scenarios ranked by descending order of Baseline Impacts",
    y = "Average Impacts"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove text (letters) above x-axis
    axis.title.x = element_text(margin = margin(t = 10)),  # Adjust space between x-axis title and labels
    legend.position = "right",  # Move legend to the right
    axis.ticks.x = element_blank()  # Remove ticks on x-axis
  ) +
  scale_x_discrete(limits = ordered_scenarios_highd)  # Order the x-axis according to the ordered_scenarios

# Display the plot
print(scatter_plot_highd)


# Load heatmap images (setting the right order manually to avoid problems, even if the same as in previous high damage baseline)
heatmap_files <- c("heatmap scen D.png", "heatmap scen A.png", "heatmap scen G.png",
                   "heatmap scen H.png", "heatmap scen B.png", "heatmap scen E.png",
                   "heatmap scen I.png", "heatmap scen C.png", "heatmap scen F.png")

# Create heatmap plots for each scenario
heatmap_plots <- list()
for (i in seq_along(heatmap_files)) {
  heatmap_image <- readPNG(heatmap_files[i])  # Read the heatmap image
  heatmap_grob <- rasterGrob(as.raster(heatmap_image), interpolate = TRUE)  # Create a rasterGrob
  
  # Create a ggplot for the heatmap
  heatmap_plot <- ggplot() + 
    annotation_custom(heatmap_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
    theme_void()  # No axes or grid lines
  
  heatmap_plots[[i]] <- heatmap_plot
}

# Combine scatter plot and heatmap plots using patchwork
library(patchwork)
heatmap_grid <- wrap_plots(heatmap_plots, ncol = length(heatmap_plots))

# Combine scatter plot and heatmap grid
final_plot_highd <- scatter_plot_highd / heatmap_grid + plot_layout(heights = c(1, 0.5))

# Display the final plot
print(final_plot_highd)


