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
#install.packages("gridExtra")
#install.packages("cowplot")

#install.packages("parallel")
#install.packages("doParallel")
#install.packages("foreach")


library(graphics)
library(rgl)
library("reshape")
library("ggplot2")
library(tibble)
library("plotly")     
library(RColorBrewer)
library("gridExtra")
library(cowplot)

library(parallel)
library(doParallel)
library(foreach)


# Detect the number of available cores
num_cores <- detectCores() - 1  # It's a good practice to leave one core free

# Set up the parallel backend to use multiple cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)



###### We first calibrate the parameters #####

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

theta = 1000  

##### The size of the matrices (i.e. the accuracy of the computation) must also be chosen here.
step = 0.1
size = 1/step + 1



### Let us define an optimization function that loops over all possible baskets for each pair of parameters and chooses for each the one maximizing utility.

opti_function <- function(alpha, beta, R, p, gamma, d_prime, theta) {
  X1 <- seq(0, R/p[1], by = step)
  X2 <- seq(0, R/p[2], by = step)
  X3 <- seq(0, R/p[3], by = step)
  X4 <- seq(0, R/p[4], by = step)
  
  Umax <- -1000  
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
rev_alphas <- seq(1, 0, by = -step)
rev_betas <- seq(1, 0, by = -step)
A <- matrix(rep(alphas, length(betas)), nrow = length(alphas), byrow = TRUE)
B <- matrix(rep(rev_betas, length(alphas)), nrow = length(betas), byrow = FALSE)

size <- nrow(A)
D <- matrix(0, nrow = size, ncol = size)
D1 <- matrix(0, nrow = size, ncol = size)
D2 <- matrix(0, nrow = size, ncol = size)
D3 <- matrix(0, nrow = size, ncol = size)
D4 <- matrix(0, nrow = size, ncol = size)
M <- matrix(0, nrow = size, ncol = size)

DD <- list()
for (k in 1:4) {
  DD[[k]] <- matrix(0, nrow = size, ncol = size)
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
#OK (99% in the ostentatious exclusive zones because of limited accuracy due to small dimensions being used)





###### Graphical representations of the quantity matrices ########

#1. Number of goods consumed 


#Computing the size of global solutions = Compute the sum of DD{1}, DD{2}, DD{3}, and DD{4} in order to get the number of goods consumed in each point
S <- DD[[1]] + DD[[2]] + DD[[3]] + DD[[4]]


# Rename rows and columns
rownames(S) <- rev_betas
colnames(S) <- alphas


# Reverse the rows of the matrix to flip the y-axis direction
S <- S[nrow(S):1, ]

# Create the heatmap
heatmap(S, scale = "none", Rowv = NA, Colv = NA,
        col = c("red", "yellow", "blue"), 
        main = "Number of goods consumed - Baseline",
        cexRow = 0.7, cexCol = 0.7)
legend("right", legend = c("1", "2", "3"), fill = c("red", "yellow", "blue"))




#2. Graphical representations of how consumers spend their budget in the different lifestyles

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
          main = paste("Share of total income spent in", lifestyle, "lifestyle"), 
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



#3. Compute and plot 'market shares' of the different goods

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

#And plot those market shares in an histogram

consumer_types <- c("GO", "GD", "BO", "BD")
percentages <- c(marketshare_GO, marketshare_GD, marketshare_BO, marketshare_BD)

# Create a data frame
data <- data.frame(Category = consumer_types, Percentage = percentages)

# Set colors
custom_colors <- c("lightgrey", "brown", "lightgreen", "darkgreen") 

# Create a bar plot (caption to be changed according to the case tested)
plot <- ggplot(data, aes(x = Category, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity") +
  labs(title = "Market shares - Reference case, Uniformly distributed population",
       x = "Lifestyles",
       y = "Percentage of quantities consumed") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Format y-axis as percentage
  theme_minimal() +
  scale_fill_manual(values = custom_colors, guide="none")

print(plot)


#### The environmental part #####

#From quantities consumed of each good, we can infer the environmental impacts induced by each good (assuming first implicitly a uniform repartition of the population)

P <- matrix(0, nrow = size, ncol = size) 
P1 <- matrix(0, nrow = size, ncol = size) 
P2 <- matrix(0, nrow = size, ncol = size)
P3 <- matrix(0, nrow = size, ncol = size)
P4 <- matrix(0, nrow = size, ncol = size)

#Create for each good matrices of environmental impacts per pixel

for (i in 1:nrow(P1)) {
  for (j in 1:ncol(P1)) {
    P1[i,j]<- gammaGO*D1[i,j]
    P2[i,j]<- gammaGD*D2[i,j]
    P3[i,j]<- gammaBO*D3[i,j]
    P4[i,j]<- gammaBD*D4[i,j]
  }
}

#Compute the total impact per pixel/consumer (to be compared between cases and then plotted)
P <- P1+P2+P3+P4


#And the contribution of each lifestyle to total pollution
impact_GO = sum(P1)
impact_GD = sum(P2)
impact_BO = sum(P3)
impact_BD = sum(P4)
total_impacts = impact_GO+impact_GD+impact_BO+impact_BD

contrib_GO = (impact_GO/total_impacts)*100
contrib_GD = (impact_GD/total_impacts)*100
contrib_BO = (impact_BO/total_impacts)*100
contrib_BD = (impact_BD/total_impacts)*100

#Plot the relative contribution of each lifestyle to pollution

consumer_types <- c("GO", "GD", "BO", "BD")
percentages <- c(contrib_GO, contrib_GD, contrib_BO, contrib_BD)

# Create a data frame
data <- data.frame(Category = consumer_types, Percentage = percentages)

# Set colors
custom_colors <- c("lightgrey", "brown", "lightgreen", "darkgreen") 

# Create a bar plot (caption to be changed according to the case tested)
plot <- ggplot(data, aes(x = Category, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity") +
  labs(title = "Contributions to environmental impacts - Reference case, Uniform distribution",
       x = "Lifestyles",
       y = "Proportion of environmental impacts") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Format y-axis as percentage
  theme_minimal() +
  scale_fill_manual(values = custom_colors, guide="none")

print(plot)





#And finally per capita impacts
totalimpacts_percapita = total_impacts/(size^2)


#Store and plot per capita impacts in the different cases

#totalimpacts_percapita_ref_unif <- totalimpacts_percapita 
#totalimpacts_percapita_ref_case2 <- totalimpacts_percapita
#totalimpacts_percapita_ref_case3 <- totalimpacts_percapita
#totalimpacts_percapita_ref_case4 <- totalimpacts_percapita
#totalimpacts_percapita_ref_case5 <- totalimpacts_percapita

categories <- c("Ref", "GO more exp.", "BO cheaper", "Income doubled", "Damage doubled")
per_capita_impacts <- c(
  totalimpacts_percapita_ref_unif,
  totalimpacts_percapita_ref_case2,
  totalimpacts_percapita_ref_case3,
  totalimpacts_percapita_ref_case4,
  totalimpacts_percapita_ref_case5)
#total_impact_percapita_unif_degrowth,
#total_impact_percapita_unif_tax

# Create a data frame
data <- data.frame(Category = factor(categories, levels = c("Ref", "GO more exp.", "BO cheaper", "Income doubled", "Damage doubled")), Numbers = per_capita_impacts)

# Reorder the levels of Category based on Numbers
data$Category <- factor(data$Category, levels = data$Category[order(data$Numbers)])

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
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # For better readability

# Display the plot
print(plot_pcimpact_unif)





# Pollution/Environmental impacts per consumer type 

# Rename rows and columns
rownames(P) <- rev_betas
colnames(P) <- alphas

# Reverse the rows of the matrix to flip the y-axis direction
P <- P[nrow(P):1, ]

# Determine the number of colors and intervals
num_colors <- 6  # Adjust as needed
color_palette_impacts <- c("darkgreen", "lightgreen", "yellow", "orange", "red", "brown")

# Calculate the range and create breaks based on equal intervals
min_val <- min(P)
max_val <- max(P)
breaks <- seq(min_val, max_val, length.out = num_colors + 1)

# Round the break values to the next highest integer
breaks_rounded <- ceiling(breaks)

# Generate labels for the legend based on the intervals of the rounded breaks
cutoff_labels <- paste0("[", breaks_rounded[-(num_colors + 1)], ", ", breaks_rounded[-1], "]")

# Create a named vector associating each unique value in the matrix with a color
color_map <- cut(P, breaks = breaks, include.lowest = TRUE, labels = color_palette_impacts)

# Specify the desired size for the main title
cex_main <- 0.8  # Adjust this value as needed

# Adjust plot margins using par() function before creating the heatmap
# Increase the bottom margin to make space for the x-axis label
par(mar = c(5, 4, 4, 8))  # c(bottom, left, top, right) - increase the right margin

# Create the heatmap
heatmap(P, scale = "none", Rowv = NA, Colv = NA,
        col = color_palette_impacts,
        main = "Environmental impacts (in units) - BO 10% cheaper",
        cexRow = 0.7, cexCol = 0.7,
        ylab = "beta")

# Add the x-axis label after creating the heatmap
mtext("alpha", side = 1, line = 3, las = 1)  # Adjust 'line' parameter to position the label properly

# Create the legend with inset parameter specifying the distance from the margins
legend("right", inset = c(-0.2, 0), legend = cutoff_labels, fill = color_palette_impacts,
       title = "Impacts/consumer",  # Add a title to the legend
       cex = 0.7, pt.cex = 1.5,  # Adjust cex and pt.cex as needed
       y.intersp = 1.5, xpd = NA)  # Allow plotting outside the plot region and adjust spacing between legend items







#### Population scenarios  ####

#Until here we have assumed the population to be uniformly distributed in the square.
#We now use Beta distributions to model various population concentration scenarios.


##SCENARIO CHOICE - Assign values to shape parameters according to the scenario 

#Basic scenario 2 (2A): two beta laws with the two same shape parameters (bell curves) to model the concentration in the middle of the distribution
a=b=c=d=2

#Variant (2b) - Only concentrated for environmental sensitivity while uniform for status: the shape parameters are changed (the beta parameter is now assumed to follow a uniform distribution, i.e. a beta distribution with c=d=1).
#a=b=2
#c=d=1

#Scenario 2B: Lower/laxer environmental norm (concentration on the left, beta law with shape parameters 5 and 1 on alpha, uniform on beta).
#a=5
#b=c=d=1

#Scenario 2C : Higher/stringent environmental norm (right-hand concentration)
#a=c=d=1
#b=5


#initialize with an empty matrix
pop <- matrix(,size,size) 

#To compute discrete densities : 

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
#just to check it sums to one, and it does !



## Application of these population scenarios to compute environmental impacts


#We first define the population of consumers of each lifestyle depending on the scenario


# Step 1 : Define a function to create a dummy matrix corresponding to consumers (or not) of the lifestyle
create_dummy_matrix <- function(matrix_input) {
  # Create the dummy matrix
  dummy_matrix <- as.numeric(matrix_input > 0)
  # Convert the dummy_matrix to the same dimensions as the original matrix
  dummy_matrix <- matrix(dummy_matrix, nrow = nrow(matrix_input), ncol = ncol(matrix_input))
  # Return the dummy matrix
  return(dummy_matrix)
}

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


#Then infer the 'market share' of each composite good
marketshare_GO_nonuniform = (quantity_GO_nonuniform/total_quantity_nonuniform)*100
marketshare_GD_nonuniform = (quantity_GD_nonuniform/total_quantity_nonuniform)*100
marketshare_BO_nonuniform = (quantity_BO_nonuniform/total_quantity_nonuniform)*100
marketshare_BD_nonuniform = (quantity_BD_nonuniform/total_quantity_nonuniform)*100


#And plot those market shares in an histogram

consumer_types <- c("GO", "GD", "BO", "BD")
percentages <- c(marketshare_GO_nonuniform, marketshare_GD_nonuniform, marketshare_BO_nonuniform, marketshare_BD_nonuniform)

# Create a data frame
data <- data.frame(Category = consumer_types, Percentage = percentages)

# Set colors
custom_colors <- c("lightgrey", "brown", "lightgreen", "darkgreen") 

# Create a bar plot (caption to be changed according to the case tested)
plot <- ggplot(data, aes(x = Category, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity") +
  labs(title = "Market shares - Reference case, Concentration in the middle",
       x = "Lifestyles",
       y = "Percentage of quantities consumed") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Format y-axis as percentage
  theme_minimal() +
  scale_fill_manual(values = custom_colors, guide="none")

print(plot)


## Let's finally compute the associated per capita impacts.

#We can easily adapt the pollution matrices according to the pop scenario
P1_nonunif = P1*go_consumers_rescaled
P2_nonunif = P2*gd_consumers_rescaled
P3_nonunif = P3*bo_consumers_rescaled
P4_nonunif = P2*bd_consumers_rescaled

#And infer the new contribution of each lifestyle to total pollution
impact_GO_nonunif = sum(P1_nonunif)
impact_GD_nonunif = sum(P2_nonunif)
impact_BO_nonunif = sum(P3_nonunif)
impact_BD_nonunif = sum(P4_nonunif)
total_impacts_nonunif = impact_GO_nonunif+impact_GD_nonunif+impact_BO_nonunif+impact_BD_nonunif

contrib_GO_nonunif = (impact_GO_nonunif/total_impacts_nonunif)*100
contrib_GD_nonunif = (impact_GD_nonunif/total_impacts_nonunif)*100
contrib_BO_nonunif = (impact_BO_nonunif/total_impacts_nonunif)*100
contrib_BD_nonunif = (impact_BD_nonunif/total_impacts_nonunif)*100

#Plot the relative contribution of each lifestyle to pollution

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



#And finally per capita impacts, that we compare with the uniform scenario
totalimpacts_percapita_nonunif = total_impacts_nonunif/(size^2)
variation_with_unif = ((totalimpacts_percapita_nonunif-totalimpacts_percapita)/totalimpacts_percapita)*100


#Store and plot per capita impacts in the different scenarios 
#A FAIRE avec un pas de 0.05





