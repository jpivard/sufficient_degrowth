########## Simulation numérique modèle simple de consommation version stage ###########


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


#library(copula)  : problème de package
#install.packages("gsl") : idem
#library(psych) : toujours la même...

######## INTRODUCTION ###############


rm(list=ls())

#Ce script est destiné au calcul des aires de consommation des différents biens définis dans un modèle simple (présenté ici https://www.overleaf.com/project/64a43062cf37016e7480691b) 
#de consommation définie selon deux dimensions  : vert/brun, ostentatoire/discret.

#Comme nos consommateurs hétérogènes sont placés dans un repère [0,1]*[0,1] en fonction de leurs préférences (en abscisses la sensibilité environnementale, en ordonnées le souci du statut social), 
#l'idée de cette méthode numérique va être de discrétiser le carré en plusieurs points/pixels (exprimant les deux types de préférences en pourcentages si le nombre de points est 101 au carré par exemple; les individus à l'origine du repère étant indifférents aux deux alors que ceux en haut à droite ont des sensibilités maximales aux deux).


#On définit la taille des matrices en fonction du pas entre les pixels.
step = 0.1
size = 1/step + 1
size

#ainsi que les jeux de paramètres
alphas <- seq(from = 0,to = 1,by=step)
betas_rev <- seq(from = 1,to = 0,by=-step) #on change le sens pour que les coordonnées apparaissent dans le bon ordre
betas <- seq(from = 0,to = 1,by=step) #juste pour les graphiques

######## DISTRIBUTIONS DE POPULATION ##############

#On définit des matrices carrées de dimension notée 'size'. 

#Pour commencer (scénario 1), on génère selon une distribution uniforme une matrice de densité de population de cette taille donnant les proportions d'individus positionnés à chaque point du repère. 
n1 <- size
n2 <- size
pop <- matrix (rep(1/(size*size), n1*n2), n1, n2)
pop
#Par ex, en choisissant un pas de 0.2, on a donc une matrice de taille 6*6, avec 36 coefficients.


#Scénario 2 : concentration de la population autour des valeurs moyennes/médianes
#On modélise ce cas par un produit de lois beta indépendantes sur chaque paramètre

pop2 <- matrix(,n1,n2) #on crée une matrice vide
 
#Scénario 2 de base: deux lois beta de même paramètres de forme 2 (courbes en cloche)pour modéliser la concentration au milieu
a=b=c=d=2

#Variante (2b): on change les paramètres de forme (on suppose maintenant le paramètre beta suivant une loi uniforme i.e. une loi beta avec c=d=1)
#a=b=2
#c=d=1


#Scénario 2B: norme environnementale plus BASSE (concentration à gauche, loi beta de paramètres de forme 5 et 1 sur alpha, uniforme sur beta)
#Même méthode, on change juste les paramètres de forme de nos betas indépendantes
#a=5
#b=c=d=1

#Scénario 2C : norme environnementale plus haute (concentration à droite)
#a=c=d=1
#b=5

#Scénario 2D : norme de statut plus basse (discrétion généralisée)
#a=b=c=1
#d=5

#Scénario 2E: norme de statut plus élevée ("Economie de la Reine Rouge", Keep Up With the Joneses...)
#a=b=d=1
#c=5




#Code pour calculer nos densités discrètes#


#On créé nos vecteurs de quantiles (en ligne)

alphas_vector <- t(matrix(alphas))
betas_vector <- t(matrix(betas))

#On peut représenter les fonctions de répartition (les mêmes en prenant les paramètres de forme identiques)

# fdr_beta_dim1 <-pbeta(alphas, shape1 = a, shape2 =b)
# plot(fdr_beta_dim1,alphas)
# 
# fdr_beta_dim2 <-pbeta(betas, shape1 = c, shape2 =d)
# plot(fdr_beta_dim2,betas)


#Partant de la fonction de répartition, on définit la densité discrète d'une beta (une pour chacune des 2 dimensions)

proba_values_dim1 <-matrix(,1,size)
proba_values_dim2 <-matrix(,1,size)


for (n in 1:(1/step)-1:1) {
  densite_beta_dim1 <- function(n) {
  densite <-  pbeta(step/2 + 1-n*step,a,b) - pbeta(step/2 + 1-(n+1)*step,a,b)
  return(densite)
  }
  for (i in 1:ncol(proba_values_dim1)-1){
   proba_values_dim1[1,i+1]= densite_beta_dim1(i)  #on stocke les valeurs dans une matrice
  }
}

densite_beta_dim1(0) <- pbeta (1-step,a,b) - pbeta (1-step/2,a,b)
densite_beta_dim1(1/step)<- pbeta(step/2,a,b) - pbeta(0,a,b)

proba_values_dim1[1,1]=densite_beta_dim1(0)
proba_values_dim1[1,size]=densite_beta_dim1(1/step)


sum(proba_values_dim1)
#On obtient les valeurs qu'il faut (somme à 1, densité en "escalier en cloche", symétrique)


#plot(densite_beta_dim1,n)  #donne une densité continue (alors que fonction définie seulement pour des entiers...) au lieu d'une discrète, mais qui a la bonne forme
#il faudrait que la densité soit fonction d'alpha et pas de n
#il manque un intervalle pour avoir 11 valeurs de probas


#Idem pour les betas (pour l'instant on prend exactement les mêmes lois)
for (n in 1:(1/step)-1:1) {
  densite_beta_dim2 <- function(n) {
    densite <-  pbeta(step/2 + 1-n*step,c,d) - pbeta(step/2 + 1-(n+1)*step,c,d)
    return(densite)
  }
  for (i in 1:ncol(proba_values_dim2)-1){
    proba_values_dim2[1,i+1]= densite_beta_dim2(i)  #on stocke les valeurs dans une matrice
  }
}

densite_beta_dim2(0) <- pbeta (1-step,c,d) - pbeta (1-step/2,c,d)
densite_beta_dim2(1/step)<- pbeta(step/2,c,d) - pbeta(0,c,d)

proba_values_dim2[1,1]=densite_beta_dim2(0)
proba_values_dim2[1,size]=densite_beta_dim2(1/step)

#On multiplie nos deux vecteurs de probabilités pour obtenir les probabilités sur l'ensemble de la population
pop2 = t(proba_values_dim2)%*%proba_values_dim1
#On a bien une plus forte concentration sur les valeurs médianes que dans le cas uniforme

sum(pop2) 
#juste pour vérifier, c'est bon !!!
        





#Scénario 3 : Lois Beta corrélées entre elles (positivement ou négativement)

#Méthode des copules ?
#Théorie : https://www.r-bloggers.com/2011/10/copulas-made-easy/
#Docu du package et méthodes : 
#https://www.r-bloggers.com/2015/10/modelling-dependence-with-copulas-in-r/
#https://www.r-bloggers.com/2016/03/how-to-fit-a-copula-model-in-r-heavily-revised-part-1-basic-tools/



#3a : On commence par modéliser une forte corrélation positive entre nos deux types de sensibilités

#Première étape : on génère une beta bivariée avec corrélation positive entre les marginales

require(mvtnorm)
set.seed(235)

S <- matrix(c(1,.8,.8,1),2,2) #Une matrice de corrélation avec des coefficients hors diagonale proches de 1
AB_copula <- rmvnorm(mean=c(0.5,0.5),sig=S,n=size*size) #On définit une normale bivariée centrée autour de la médiane de l'intervalle et avec le bon nombre d'observations

U <- pnorm(AB_copula) #On convertit en uniforme - se vérifie avec hist(U[,1]) or hist(U[,2])

alpha_cop <- qbeta(U[,1],2,2) 
beta_cop <- qbeta(U[,2],2,2) 
#On choisit ici comme lois marginales deux Beta de même paramètres de forme, correspondant à une concentration des valeurs vers le centre

#plot(alpha_cop,beta_cop)  #On obtient notre répartition de population en fonction des paramètres, qu'il faut maintenant convertir en matrice de probabilités
#On obtient la densité de cette distribution bivariée à l'aide d'une fonction du package :

library(MASS)


















############# PARAMETRISATION DES FRONTIERES ####################

#Définissons nos frontières délimitant les zones de consommation en fonction des paramètres de prix et des dérivées (afin de pouvoir étudier les différents cas en changeant seulement les valeurs suivantes).

#on initialise les paramètres de nos "frontières", à savoir les prix des quatre biens, l'utilité et les dommages marginaux (en respectant les hypothèses du modèle)
p_vo = 2.9
p_bo = 2
p_vd= 3/2
p_bd= 1  #même valeur dans tous les cas
u_prime = 0.5
d_prime = 0.2
#u_prime forcément plus grand que d_prime

#On fait une boucle imbriquée sur nos individus (ie leurs paramètres alpha et beta, que l'on va convertir - en normalisant- en coordonnées (x,y) où les coordonnées sont comprises entre 1 et la taille "size").
#Pour chaque point (ie paire de paramètres) on regarde quels contraintes sont satisfaites, ce qui détermine le type de bien consommé en ce point (et la matrice correspondante parmi les 4 définies plus haut qui va être remplie par un "1" au point (x,y) associé).

#L'idée va être de créer quatre matrices de booléens size*size (vrai = 1 /faux = 0), une par type de consommateur ; et d'examiner pour chaque point du carré (à l'aide de boucles imbriquées)  s'il satisfait les contraintes correspondant à un certain type.
#Dans un second temps, on vérifiera que la somme de ces quatre matrices donne la matrice 1 de taille size*size



#### BROUILLON 



#Méthode 1 : boucles sur alpha et beta

#data <- matrix(, nrow = size, ncol = size)  #on créé une matrice vide


#alphas <- seq(from = 0,to = 1,by=step)
#betas <- seq(from = 0,to = 1,by=step)

#for (beta in betas) {
  
#  for (alpha in alphas){
    #    print(c(alpha,1-beta))
    
#   }
#}
# on a bien les bons couples de paramètres dans l'ordre souhaité (avec (0,0) en bas à gauche comme dans un repère si on arrive à remplir la matrice par le haut))




#Autre méthode

#data1 <- data
#data1 <- sapply(1:size, function(i) sapply(1:size, function(j) paste((j-1)*step,(size-i)*step,sep = ", ")))
#data1 <- t(data1)

#On a bien notre matrice de coordonnées (alpha,beta) avec les couples placés comme on le souhaite!!!




#On veut créer une matrice de booléens avec uniquement des 1 (True) et des 0 (False) selon que (alpha,beta) satisfont certaines contraintes.


#for (i in 1:nrow(data1)){
  
#  for (j in 1:ncol(data1)){
    
#    if ((1-i)*j > ((p_vo-p_vd)*(u_prime-d_prime*i)))  {
      
#      data1[i,j]= 1}
    
#    else  {data1[i,j]= 0}
    
#  }
#}

#On a bien une matrice de booléens mais je pense qu'il ne teste pas les coordonnées (alpha,beta) converties mais les coordonnées (i,j) entre 1 et 6
#Trop difficile avec une matrice où chaque coefficient est un couple de coordonnées ?








#Méthode 3 : intégrer les booléens à la boucle et convertir les coordonnées


#data <- matrix(, nrow = size, ncol = size)  #on créé une matrice vide
#data2 <- data

#alphas <- seq(from = 0,to = 1,by=step)
#betas_rev <- seq(from = 1,to = 0,by=-step) #on change le sens pour que les coordonnées apparaissent dans le bon ordre

#for (beta in betas_rev) {
  
#  for (alpha in alphas){
    
   #print(c(alpha,beta))  #les paramètres apparaissent dans le bon ordre
   
#   x=size-beta/step
#    y=alpha/step+1
#    print(c(x,y))         #les coordonnées sont bonnes (correspond à celles d'une matrice "classique")
    
#    if ((1-alpha)*beta>((p_vo-p_vd)*(u_prime-d_prime*alpha))) {
#      bool_f = 1
#    }
#    else {bool_f = 0}
    
    
#    data2[x,y]= bool_f
      
#  }
    
#}

#ça marche !!!!!!! 





######## CALCUL DES MATRICES DE ZONES DE CONSOMMATION ###############

#On généralise pour trouver nos zones de consommation de chacun des quatre biens.

mat_vo <- matrix(, nrow = size, ncol = size)  #on créé 4 matrices vides
mat_vd <- matrix(, nrow = size, ncol = size)
mat_bo <- matrix(, nrow = size, ncol = size)
mat_bd <- matrix(, nrow = size, ncol = size)
test<-matrix(, nrow = size, ncol = size)

alphas <- seq(from = 0,to = 1,by=step)
betas_rev <- seq(from = 1,to = 0,by=-step) #on change le sens pour que les coordonnées apparaissent dans le bon ordre
betas <- seq(from = 0,to = 1,by=step)

for (beta in betas_rev)  {
  
  for (alpha in alphas) {
    
    #alpha=0.29  #tester avec 0.28, 0.29 et 0.3 par exemple pour montrer le problème : la colonne 30 est ignorée par R...
    x=size-beta/step
    y=alpha/step+1
    
  
    if ((1-alpha)*beta>((p_vo-p_vd)*(u_prime-d_prime*alpha))) {
      bool_f = 1
    }
    else {bool_f = 0}
    
    if (((1-alpha)*beta>((p_bo-p_bd)*(u_prime-d_prime*alpha)))) {
      bool_g = 1
    }
    else {bool_g = 0}
    
    if (alpha*beta> -alpha*d_prime + ((p_vo-p_bo)*(u_prime-d_prime*alpha))) {
      bool_h = 1
    }
    else {bool_h = 0}
    
    if (alpha*beta>-alpha*d_prime + ((p_vd-p_bd)*(u_prime-d_prime*alpha))) {
      bool_k = 1
    }
    else {bool_k = 0}
    
    if (beta> -d_prime*alpha + ((p_vo-p_bd)*(u_prime-d_prime*alpha))) {
      bool_i = 1
    }
    else {bool_i = 0}
    
    if ((1-2*alpha)*beta>(d_prime*alpha + (p_bo-p_vd)*(u_prime-d_prime*alpha))) {
      bool_j = 1
    }
    else {bool_j = 0}
    
    mat_vo[x,y]= bool_f*bool_h*bool_i
    mat_bo[x,y]= bool_g*bool_j*(1-bool_h)
    mat_vd[x,y]= bool_k*(1-bool_f)*(1-bool_j)
    mat_bd[x,y]= (1-bool_g)*(1-bool_k)*(1-bool_i)
    test[x,y]= mat_vo[x,y]+mat_bo[x,y]+mat_vd[x,y]+mat_bd[x,y]          #Vérifions que la somme de nos matrices donne bien la matrice 1. 
    
  }
  
}
  
#Cas des individus sur les frontières : la moitié consomme l'un des biens, l'autre moitié l'autre bien
#Le code tourne mais le problème persiste voire s'aggrave (plus de NA, des 2...) 

# for (beta in betas_rev)    {
#   
#   for (alpha in alphas) {
#     
#     x=size-beta/step
#     y=alpha/step+1
#     
#     if ((1-alpha)*beta==((p_vo-p_vd)*(u_prime-d_prime*alpha))) {
#       bool_frontiere_f = 1
#     }
#     else {bool_frontiere_f = 0}
#     
#     if ((1-alpha)*beta == ((p_bo-p_bd)*(u_prime-d_prime*alpha))) {
#       bool_frontiere_g = 1
#     }
#     else {bool_frontiere_g = 0}
#     
#     if (alpha*beta== -alpha*d_prime + (p_vo-p_bo)*(u_prime-d_prime*alpha)) {
#       bool_frontiere_h = 1
#     }
#     else {bool_frontiere_h = 0}
#     
#     if (alpha*beta==-alpha*d_prime + (p_vd-p_bd)*(u_prime-d_prime*alpha)) {
#       bool_frontiere_k = 1
#     }
#     else {bool_frontiere_k = 0}
#     
#     if (beta== -d_prime*alpha + (p_vo-p_bd)*(u_prime-d_prime*alpha)) {
#       bool_frontiere_i = 1
#     }
#     else {bool_frontiere_i = 0}
#     
#     if ((1-2*alpha)*beta == d_prime*alpha + (p_bo-p_vd)*(u_prime-d_prime*alpha)) {
#       bool_frontiere_j = 1
#     }
#     else {bool_frontiere_j = 0}
#     
#     mat_vo[x,y]= bool_f*bool_h*bool_i 
#     mat_bo[x,y]= bool_g*bool_j*(1-bool_h) + bool_frontiere_h
#     mat_vd[x,y]= bool_k*(1-bool_f)*(1-bool_j)+bool_frontiere_f+bool_frontiere_g
#     mat_bd[x,y]= (1-bool_g)*(1-bool_k)*(1-bool_i)+bool_frontiere_g+bool_frontiere_k+bool_frontiere_i
#     test[x,y]= mat_vo[x,y]+mat_bo[x,y]+mat_vd[x,y]+mat_bd[x,y]          #Vérifions que la somme de nos matrices donne bien la matrice 1. 
#     
#   }
#   
# }

#La matrice somme des 4 matrices de booléens est bien la matrice 1 
#Après correction d'une erreur, c'est toujours le cas même pour des matrices plus grandes et d'autres valeurs des prix !!!
#Mais quand la dimension grandit beaucoup, des NA apparaissent sur les frontières...




####### REPRESENTATIONS GRAPHIQUES DES ZONES ##########

#Exporter les plots dans le cas 101*101 (step=0,01) 
#Essayer d'améliorer encore les couleurs et leur échelle : https://r-graph-gallery.com/27-levelplot-with-lattice.html


#VO

mat1<- mat_vo %>% 
  as.tibble() %>% 
  as.list() %>% 
  lapply(function(x) rev(x)) %>% 
  do.call(rbind, .) %>% 
  as.matrix()

colnames(mat1)<-alphas
row.names(mat1)<-betas

#convert to binary scale on the plot (currently not working)
#mat1 <- raster(mat1)
#mat1 <- ratify(mat1)
#rat <- levels(mat1)[[1]]
#rat$code <- c(0,1)
#levels(mat1) <- rat

levelplot(mat1, main="zone de consommation de VO (en violet)", xlab="alpha", ylab="beta", col.regions = cm.colors(121))   



#BO

mat2<- mat_bo %>% 
  as.tibble() %>% 
  as.list() %>% 
  lapply(function(x) rev(x)) %>% 
  do.call(rbind, .) %>% 
  as.matrix()

colnames(mat2)<-alphas
row.names(mat2)<-betas


levelplot(mat2, main="zone de consommation de BO (en violet)", xlab="alpha", ylab="beta",col.regions = cm.colors(121))   

#VD

mat3<- mat_vd %>% 
  as.tibble() %>% 
  as.list() %>% 
  lapply(function(x) rev(x)) %>% 
  do.call(rbind, .) %>% 
  as.matrix()

colnames(mat3)<-alphas
row.names(mat3)<-betas

levelplot(mat3, main="zone de consommation de VD (en bleu ciel)", xlab="alpha", ylab="beta")    


#BD

mat4<- mat_bd %>% 
  as.tibble() %>% 
  as.list() %>% 
  lapply(function(x) rev(x)) %>% 
  do.call(rbind, .) %>% 
  as.matrix()

colnames(mat4)<-alphas
row.names(mat4)<-betas

levelplot(mat4, main="zone de consommation de BD (en bleu ciel)", xlab="alpha", ylab="beta")  








######CALCUL DES PARTS DE CONSO POUR CHAQUE BIEN########
 
###### 1. Scénario d'une distribution uniforme de la population ######

pop_vo_unif = pop*mat_vo
pop_bo_unif = pop*mat_bo
pop_vd_unif = pop*mat_vd
pop_bd_unif = pop*mat_bd

conso_vo_scenar1 = sum(pop_vo_unif)    #ne marche plus quand les matrices sont trop grandes (alors que tout le reste si) car on a des NA sur les frontières...
conso_bo_scenar1 = sum(pop_bo_unif)
conso_vd_scenar1 = sum(pop_vd_unif)
conso_bd_scenar1 = sum(pop_bd_unif)

#en pourcentages
conso_vo_scenar1*100
conso_bo_scenar1*100
conso_vd_scenar1*100
conso_bd_scenar1*100


#Ainsi, avec une distribution uniforme de la population (scénario 1), on obtient :


#Dans le cas de référence, avec un découpage (11*11) du carré, la part de marché de VO tend vers zéro (0,008), alors que BO est à 15%, BD à 29% et VD à 55%. 
#Un découpage un peu plus fin (41*41) donne les parts suivantes (dans le même ordre) : 0.4% ; 14,2%; 26,7%; 58,7%.

#Pour la suite, revenons au découpage 41*41 (le plus fin pour lequel je puisse obtenir les pourcentages pour l'instant) :

#Cas 1b (avec le prix de VO égal à trois) : VO est chassé du marché ; les autres ont peu bougé (14,4, 26,7 et 58,9%)

#Cas 2 (prix dans l'ordre = 9/4;2;7/4) : BD à 34%; VO à 32%; VD à 28% et BO à 6%

#Cas 3 (prix dans l'ordre = 5/2;2;3/2) : VD à 53%; BD à 27%; BO à 11%, VO à 10%

#Cas 4 (cas 1 mais avec doublement du revenu) : VD à 41%; VO à 34%; BO à 17% et BD à 8%

#Cas 5 (cas 1 mais avec doublement des dommages marginaux): VD à 64%; BD à 16%; BO à 13%; VO à 7%



###### 2. Scénario d'une concentration de la population autour de certaines valeurs (produit de lois betas indépendantes) ######

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

#Dans certains cas, ces changements de distribution ne changent pas la hiérarchie entre les biens, mais uniquement les parts de consommation ; dans d'autres l'ordre est bouleversé
#On garde un découpage 41*41 pour les résultats numériques (problème à résoudre)

#Cas de référence : l'ordre n'est pas modifié, les parts deviennent maintenant dans l'ordre croissant  :
#0,05% pour VO (en baisse // au scénario uniforme) ;
#8,8% pour BO (en forte baisse) 
#22,3% pour BD (en baisse)
#68,9% pour VD (en hausse)
#On remarque que dans ce scénario seul VD profite de la concentration des consommateurs autour de la moyenne
#C'est le contraire pour les autres biens, surtout BO dont la zone est celle qui recoupe le plus des valeurs extrêmes dans ce premier cas (et dans ce scénario il y a peu de monde pour ces valeurs)

#Deuxième cas : ordre modifié, le passage au scénario 2 est favorable au bien VO qui devient le plus consommé avec 42%
#Devant BD avec 30% (en baisse), VD avec 26% (en baisse) et BO à seulement 2% (forte baisse)

#Troisième cas : ordre modifié, le passage au scénario 2 fait tomber BO à la dernière place à 6% (forte baisse, comme dans le cas précédent)
#Sinon VD est toujours le plus consommé mais progresse à 63%, devant BD en baisse à 22% et VO à peu près stable à 9%

#Quatrième cas : ordre modifié, le passage au scénario 2 est de nouveau favorable au bien VO (le plus cher) qui progresse fortement, à 46%
#Il est suivi des biens VD (39%, en baisse), BO (12%, en baisse) et BD à 3%

#Cinquième cas : l'ordre n'est pas modifié, seules les parts ont changé, au profit du bien VD qui progresse à 78%
# BD, BO et VO descendent respectivement à 10%, 8% et 4%


##### 2b. Concentration au milieu pour les sensibilités environnementales seulement (beta avec a=b=2 toujours) mais distribution uniforme des préférences pour le statut (c=d=1)

#Cas de référence : par rapport au scénario précédent, l'ordre n'est pas modifié, mais les parts évoluent sensiblement, avec respectivement 63,6%, 26%, 10% et 0,4%

#Cas 2 : Ordre pas modifié, mais BD progresse à 35% alors que VD descend à 21%

#Cas 3 : Idem pour l'ordre, VD redescend à 55%, devant BD à 26%, VO à 13% et BO à 6%

#Cas 4: Ordre pas modifié, mais VO baisse un peu à 45%, VD stable à 39,5%, BO baisse à 9,5% et BD monte à 6%

#Cas 5 : VD toujours largement en tête à 71% mais l'ordre entre les 3 autres n'est plus aussi clair :  BD semble un peu devant 11% puis les biens ostentatoires suivent autour de 9%


##### 2bis et ter. Concentrations sur des valeurs inférieures/supérieures de sensibililté environnementale (et toujours uniforme sur le statut)

### 2B : Société moins préoccupée par l'environnement

#Cas 1 : VD en baisse à 51, devant BD en hausse à 39 ; les autres sont stables

#Cas 2 : BD en forte hausse à 51 passe en tête, devant VO à 25 et BO à 23, VD 

#Cas 3: BD en forte hausse à 51, devant BO à 33,5, VO à 8,5 et VD à 7

#Cas 4 : BO à 39, VO à 34, BD à 22 et VD à 5

#Cas 5 : BD 43; BO 31; VO 16 et VD 10

### 2C : Société plus préoccupée par l'environnement

#Cas 1 : 97 pour VD ! 3 pour BD et les autres sont presque inexistants

#Cas 2 : 86 pour VD, VO à 14, BD et BO presque à 0

#Cas 3 : 95,5 VD; 3 BD; 1,5 VO ; BO toujours presque à 0

#Cas 4 : 78 VD; 22 VO

#Cas 5 : VD 98; VO 1,5; BD 0,2; BO 0,1


### 2D: Société moins préoccupée par le statut

#Cas 1 : BD 58; VD 42 et BO 0,5

#Cas 2 : VD 62; BD 31 et VO 7

#Cas 3 : mêmes parts que dans le premier cas

#4 : VD 66; BD 25; BO 6 et VO 3

#5 : 66,5 VD; 33 BD; 0,5 BO et VO exclu du marché


### 2E : Economie de la Reine Rouge

#1 : 66 VD; 31 BO; 2 BD et 1 VO

#2 : 78 VO; 22 VD

#3 : 48,5 VD ; 29 VO; 20,5 BO et 2 pour BD

#4 : 67 VO; 19 BO et 14 VD

#5 : 51 VD; 26 BO; 22 VO et 1 pour BD


