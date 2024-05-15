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

library(graphics)
library(rgl)
library("reshape")
library("ggplot2")
library(tibble)
library("plotly")     
library(RColorBrewer)
library("gridExtra")
library(cowplot)




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
step = 0.05
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
rownames(S) <- rev_b
colnames(S) <- a


# Reverse the rows of the matrix to flip the y-axis direction
S <- S[nrow(S):1, ]

# Create the heatmap
heatmap(S, scale = "none", Rowv = NA, Colv = NA,
        col = c("red", "yellow", "blue"), 
        main = "Number of goods consumed - Baseline",
        cexRow = 0.7, cexCol = 0.7)
legend("right", legend = c("1", "2", "3"), fill = c("red", "yellow", "blue"))




#2. Graphical representations of how consumers spend their budget in the different lifestyles

generate_heatmap <- function(share, lifestyle, color_palette, legend_title) {
  # Set the row and column names
  rownames(share) <- rev_b
  colnames(share) <- a
  
  # Determine the number of colors and percentage cutoffs
  num_colors <- 5  # Adjust as needed
  cutoffs <- seq(0, 100, by = 100 / (num_colors - 1))  # Ensure cutoffs match num_colors
  
  # Specify the desired size for the main title
  cex_main <- 0.8  # Adjust this value as needed
  
  # Set the global graphical parameter for main title size
  par(cex.main = cex_main)
  
  # Reverse the rows of the matrix to flip the y-axis direction
  share_reversed <- share[nrow(share):1, ]
  
  # Plot the heatmap with the reversed matrix
  heatmap(share_reversed, scale = "none", Rowv = NA, Colv = NA,
          col = color_palette,
          main = paste("Share of total income spent in", lifestyle, "lifestyle (in percentage) - BO 10% cheaper"),
          cexRow = 0.7, cexCol = 0.7,
          ylab = "beta")  # Add labels for x and y axes
  
  # Rotate the x-axis label horizontally
  mtext("alpha", side = 1, line = 2, las = 1)
  
  # Create breaks for color mapping
  breaks <- seq(0, 100, length.out = num_colors + 1)  # num_colors + 1 intervals
  
  # Create a named vector associating each unique value in the matrix with a color
  color_map <- cut(share_reversed, breaks = breaks, include.lowest = TRUE,
                   labels = color_palette)
  
  # Adjust plot margins using par() function
  par(mar = c(4, 4, 4, 6))  # Increased right margin for better label alignment
  
  # Create the legend with inset parameter specifying the distance from the margins
  legend("right", inset = c(0, 0), legend = as.character(cutoffs), fill = color_palette,
         title = legend_title,  # Add a title to the legend
         cex = 0.7, pt.cex = 1.5,  # Adjust cex and pt.cex as needed
         y.intersp = 1.5, xpd = TRUE)  # Adjust spacing between legend items as needed
}


# Define color palettes for different lifestyles
num_colors <- 5  # Ensure this matches the number of intervals
color_palette_GO <- colorRampPalette(c("white", "darkgreen"))(num_colors)
color_palette_GD <- colorRampPalette(c("white", "darkgreen"))(num_colors)
color_palette_BO <- colorRampPalette(c("white", "brown"))(num_colors)
color_palette_BD <- colorRampPalette(c("white", "brown"))(num_colors)

# Generate heatmap plots for different lifestyles
heatmap_GO <- generate_heatmap(shareGO, "GO", color_palette_GO, "Shares")
heatmap_GD <- generate_heatmap(shareGD, "GD", color_palette_GD, "Shares")
heatmap_BO <- generate_heatmap(shareBO, "BO", color_palette_BO, "Shares")
heatmap_BD <- generate_heatmap(shareBD, "BD", color_palette_BD, "Shares")

# Arrange the heatmaps in a grid
heatmap_grid <- plot_grid(heatmap_GO, heatmap_GD, heatmap_BO, heatmap_BD, ncol = 2)

# Display the grid
print(heatmap_grid)










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
  labs(title = "Market shares - BO 10% cheaper, Uniformly distributed population",
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
  labs(title = "Contributions to environmental impacts - BO 10% cheaper, Uniform distribution",
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
rownames(P) <- rev_b
colnames(P) <- a

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







#### Population scenarios (later) ####






