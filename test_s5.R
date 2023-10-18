library (ade4)
library(tidyverse)
library(heatmaply)

setwd("C:/Users/fanny/Desktop/stats R/project")
data (aravo)
spe <- aravo[[1]]
env <- aravo[[2]]
traits <- aravo[[3]]

env$Form <- as.numeric(env$Form)
env$ZoogD <- ifelse(env$ZoogD == "no", 0, env$ZoogD)
env$ZoogD <- ifelse(env$ZoogD == "some", 1, env$ZoogD)
env$ZoogD <- ifelse(env$ZoogD == "high", 2, env$ZoogD)

#standardize
env_std <- scale(env)
#env_std[,c("Aspect", "Slope", "PhysD","Snow")] <- scale(env[, c("Aspect", "Slope", "PhysD","Snow")])

heatmaply(env_std)
