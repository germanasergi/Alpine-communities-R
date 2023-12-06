setwd("C://UNI/R")

#load useful libraries
library (ade4)
library(vegan)
library(tidyverse)
library(heatmaply)
library(randomForest)
library(gplots)
library(ggplot2)
library(ggpubr)
library("GGally")
library("tidyr")
library(pastecs)
library(psych)
library(NbClust)
library(dendextend)
library(partykit)

data (aravo)

spe<-aravo[[1]]
env<-aravo[[2]]
traits<-aravo[[3]]

##change of variables
#Form
env$Form <- as.numeric(env$Form)

#Zoological disturbance
ZoogD_mapping <- c("no"=0, "some"=1, "high"=2)
env$ZoogD <- ZoogD_mapping[env$ZoogD]

# Center and scale = standardize all the variables
env_std <- decostand(env, "standardize")
apply(env_std, 2, mean)	# check if means = 0
apply(env_std, 2, sd)
#heatmaply(env_std)

# Center and scale = standardize all the speices
spe_std <- decostand(spe, "standardize")

#DATA RESEMBLANCE

##Q-mode dissimilarity and distance measures for species

#modification of image() to include row and column axis labels and to construct a heatmap of species
par(mfrow=c(1,1))
image(as.matrix(spe)) #original
image.real <- function(mat) { 
  mat <- t(mat)[,nrow(mat):1]
  image(mat, axes = FALSE, col = hcl.colors(15, palette="viridis"))
  axis(1, at = seq(0, 1, length = nrow(mat)), labels = rownames(mat), las=2)
  axis(2, at = seq(0, 1, length = ncol(mat)), labels = colnames(mat))
  box() 
}
image.real(as.matrix(spe)) #modified
heatmap.2(as.matrix(spe), Colv = FALSE, Rowv =FALSE, dendrogram="none", trace="none")

#percentage difference (Bray-Curtis) dissimilarity matrix on raw species
spe.db <- vegdist(spe)
image.real(as.matrix(spe.db))

#custom function coldiss (Francois Gillet)
coldiss(dune.db, nc = 16, diag = TRUE)

#percentage difference (Bray-Curtis) dissimilarity matrix on log-transformed abundances
spe.dbln <- vegdist(log1p(spe))
image.real(as.matrix(spe.dbln))

# Compute matrix of chord distances among sites
spe.norm <- decostand(spe, "normalize")
spe.ch <- vegdist(spe.norm, "euc")


#Hellinger distance matrix - one type of euclidian distance
spe.hel <- decostand(spe, "hel")
spe.dh <- dist(spe.hel)
image.real(as.matrix(spe.dh))

#Jaccard dissimilarity matrix for binary data
spe.dj <- vegdist(spe, "jac", binary = TRUE)
image.real(as.matrix(spe.dj))


#UNSUPERVISED CLASSIFICATION

# Compute single linkage agglomerative clustering
spe.ch.single <- hclust(spe.ch, method = "single")

# Plot the dendrogram using default options
plot(spe.ch.single, main = "Chord - Single linkage")


# Compute and plot complete-linkage agglomerative clustering
spe.ch.complete <- hclust(spe.ch, method = "complete")
plot(spe.ch.complete, main = "Chord - Complete linkage")


# Compute UPGMA clustering
spe.ch.UPGMA <- hclust(spe.ch, method = "average")
plot(spe.ch.UPGMA, main = "Chord - UPGMA")


# Compute UPGMC clustering
spe.ch.centroid <- hclust(spe.ch, method = "centroid")
plot(spe.ch.centroid, main = "Chord - Centroid")


# Compute Ward's minimum variance clustering
#when there are branches pointing up (ex 17-18-19) means that there is more dissimilarity
spe.ch.ward <- hclust(spe.ch, method = "ward.D2")
plot(spe.ch.ward,  main = "Chord - Ward")


#k-means clustering, you can change the number 5 to change the number of groups for splitting
spe.ch.k<-  kmeans(spe.ch, centers=5)
spe.ch.k


# Cophenetic correlations, gives coefficient of original dissimilarity of data
# Single linkage clustering
spe.ch.single.coph <- cophenetic(spe.ch.single)
cor(spe.ch, spe.ch.single.coph)


# Complete linkage clustering
spe.ch.comp.coph <- cophenetic(spe.ch.complete)
cor(spe.ch, spe.ch.comp.coph)


# Average clustering
spe.ch.UPGMA.coph <- cophenetic(spe.ch.UPGMA)
cor(spe.ch, spe.ch.UPGMA.coph)


# Ward clustering
spe.ch.ward.coph <- cophenetic(spe.ch.ward)
cor(spe.ch, spe.ch.ward.coph)

##find optimal number of clusters

#looking for 'ch' maximum distance, best number of clusters and best partition
Nb.UPGMA<-NbClust(spe, diss=spe.ch, distance = NULL, min.nc=3, max.nc=30, 
                  method = "average", index = "ch")
Nb.UPGMA
plot(3:30,Nb.UPGMA$All.index, xlab="Number of clusters", ylab="Calinski and Harabasz index")
abline(v=5, col="red", lty=2)

plot(spe.ch, spe.ch.ward.coph,
     xlab = "Chord distance",
     ylab = "Cophenetic distance",
     asp = 1, xlim = c(0, sqrt(2)),
     ylim = c(0, sqrt(2)),
     main = c("Single linkage", paste("Cophenetic correlation =", round(cor(spe.ch, spe.ch.ward.coph), 3))))
abline(0, 1) #addig line to the plot (incidence, slope)
lines(lowess(spe.ch, spe.ch.ward.coph, f = 0.1), col = "red", lwd=3)


#SUPERVISED CLASSIFICATION
spe.cw.g <- cutree(spe.ch.ward, 4)
spe.cw.g

D<-cbind(spe.cw.g,env,spe)

class.groups=ctree(as.factor(spe.cw.g)~ PhysD+Slope+ZoogD+Aspect+Form+Snow, data = D)
plot(class.groups)
?ctree

class.snow<-ctree(ZoogD~PhysD+Slope+Snow+Aspect+Form, data = D)
plot(class.snow)

rf <- randomForest(as.factor(spe.cw.g)~., env, ntree=500, mtry=, importance=TRUE,
                 na.action=na.omit,do.trace=100,proximity=T)
plot(rf)

#partial plot
par(mfrow=c(2,3))
partialPlot(rf, env, Aspect)
partialPlot(rf, env, Slope)
partialPlot(rf, env, Form)
partialPlot(rf, env, PhysD)
partialPlot(rf, env, ZoogD)
partialPlot(rf, env, Snow)



### ORDINATION

##SPECIES
#Compute CA
spe.ca <- cca(spe)
summary(spe.ca)		# default scaling 2
summary(spe.ca, scaling = 1)

# CA biplots
par(mfrow = c(1, 2))
# Scaling 1: sites are centroids of species
plot(spe.ca, 
     scaling = 1, 
     main = "CA abundances - biplot scaling 1"
)
# Scaling 2 (default): species are centroids of sites
plot(spe.ca, main = "CA abundances - biplot scaling 2")

#to remove the horseshoe we use DCA
spe.DCA<-decorana(spe, iweigh=0, iresc=4, ira=0, mk=26, short=0, before=NULL, after=NULL)
spe.DCA
plot(spe.DCA, disp="sites")
plot(spe.DCA, scaling = 1, main = "PCA - scaling 1")

# Adjusted R^2 retrieved from the rda object
(R2adj <- RsquareAdj(spe.rda)$adj.r.squared)

#try nmds or to filter before cca
spe.nmds <- metaMDS(spe, distance = "bray")
spe.nmds$stress
plot(spe.nmds, type = "t",main = paste("NMDS Bray Curtis; Stress =",round(spe.nmds$stress, 3)))

par(mfrow = c(1, 2))
stressplot(spe.nmds, main = "Shepard plot")
gof <- goodness(spe.nmds)
gof
plot(spe.nmds, type = "t", main = "Goodness of fit")
points(spe.nmds, display = "sites", cex = gof * 300) #to see the goodness of fit if the sites


#ENV VARIABLES
#Compute RDA
#to compute rda we need non numerical variable
#this is to substitute the ZoogD variable with 3 variables, one for each value
encoded_vars <- model.matrix(~ ZoogD - 1, data = D)
env_vars_processed <- cbind(env[, -which(names(env) %in% "ZoogD")], encoded_vars)

env.pca <- rda(env_vars_processed, scale = TRUE, row.wt = spe.ca$lw) #ADDED THE WEIGHT HERE
env.pca
summary(env.pca)
par(mfrow = c(1, 2))
biplot(env.pca, main = "PCA - scaling 2")
biplot(env.pca, scaling = 1, main = "PCA - scaling 1")

#RDA constrained by Hel. matrix
spe.rda <- rda(spe.hel ~ ., env) #"." means constrained for all the env par in env3
summary(spe.rda)

# Triplots of the rda results (lc site scores)
# Site scores as linear combinations of the environmental variables
par(mfrow=c(1,3))
plot(spe.rda, scaling = 1,   display = c("lc"), main = "RDA - sites")
plot(spe.rda, scaling = 1,   display = c("sp"), main = "RDA - species")
plot(spe.rda, scaling = 1,   display = c("cn"), main = "RDA - constraints")

spe.good <- goodness(spe.rda)
sel.sp <- which(spe.good[, 2] >= 0.6)
anova(spe.rda, permutations = how(nperm = 999))
vif.cca(spe.rda)

mod0 <- rda(spe.hel ~ 1, data = env)
spe.rda.all <- rda(spe.hel ~ ., data = env)
step.forward <- ordistep(mod0, scope = formula(spe.rda.all), direction = "forward", permutations = how(nperm = 999))
RsquareAdj(step.forward)


#TRAITS
#Compute RDA
traits.pca <- rda(traits, scale = TRUE, row.wt = spe.ca$cw) #ADDED THE WEIGHT HERE
traits.pca
summary(traits.pca)
par(mfrow = c(1, 2))
biplot(traits.pca, main = "PCA - scaling 2")
biplot(traits.pca, scaling = 1, main = "PCA - scaling 1")
arrows(0, 0,  traits.pca[, 1], traits.pca[, 2])
text(traits.pca, disp = "species",scaling=1)

#plot env and traits together
par(mfrow = c(1, 2))
biplot(env.pca, scaling = 1, main = "PCA on Environmental variables")
biplot(traits.pca, scaling = 1, main = "PCA on Traits")
#custom_pink <- rgb(255, 170, 100, maxColorValue = 255)
text(traits.pca, disp = "species",scaling=1, col = "pink")








##PAPER

#rlq
afcL.aravo <- dudi.coa(aravo$spe, scannf = FALSE)
acpR.aravo <- dudi.hillsmith(aravo$env, row.w = afcL.aravo$lw,
                             scannf = FALSE)
acpQ.aravo <- dudi.pca(aravo$traits, row.w = afcL.aravo$cw,
                       scannf = FALSE)
rlq.aravo <- rlq(acpR.aravo, afcL.aravo, acpQ.aravo,
                 scannf = FALSE)
plot(rlq.aravo)

par(mfrow = c(1, 3))
s.arrow(rlq.aravo$l1)
s.arrow(rlq.aravo$c1)
s.label(rlq.aravo$lQ, boxes = FALSE)

#4th corner
nrepet <- 999
four.comb.aravo <- fourthcorner(aravo$env, aravo$spe,
                                aravo$traits, modeltype = 6, p.adjust.method.G = "none",
                                p.adjust.method.D = "none", nrepet = nrepet)
plot(four.comb.aravo)


#LDA
library(MASS)
env.lda <- lda(Snow ~ ., as.data.frame(env))
plot(env.lda, 
     col=c("red","green","blue")[env$Snow])