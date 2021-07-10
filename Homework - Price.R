####################
##  1. Libraries  ##
####################

library(fractal)
detach("package:fractal", unload = TRUE)
library(nonlinearTseries)


##############
#  2. Input  #
##############

path <- "C:/Users/thxsg/OneDrive - International Hellenic University/1. IHU (2021)/5. Timeseries Forecasting/2. Non-Linear/2. Homework/EP_Italy_CentralSouth3.csv"
clm <- 5


############################
##  3. Setting Variables  ##
############################

#x <- read.table("clipboard" , sep = "\t", header = TRUE)
x <- read.csv(path)

z <- x[[clm]]
len <- length(z)


####################
##  4. Solutions  ##
####################

#######
#  a  #

plot.ts(z)

#######
#  b  #

d <-timeLag(
  z,
  technique = "ami",
  selection.method = "first.minimum",
  lag.max = 15,
  do.plot = TRUE,
  main = "Mutual Information - First Minimum")
d

d2 <- timeLag(
  z,
  technique = "acf",
  selection.method = "first.e.decay",
  lag.max = 15,
  do.plot = TRUE,
  main = "Autocorrelation - 1/e")

d3 <- timeLag(
  z,
  technique = "acf",
  selection.method = "first.zero",
  lag.max = 15,
  do.plot = TRUE,
  main = "Autocorrelation - 0")


#######
#  c  #

qq <- as.data.frame(buildTakens(z, embedding.dim = 2, time.lag = d))
plot(qq[, 1], qq[, 2],
     pch = ".",
     main = "Phase Portrait - X, Y",
     xlab = "X",
     ylab = "Y")


#######
#  d  #

m <- estimateEmbeddingDim(
  z,
  number.points = len,
  time.lag = d,
  max.embedding.dim = 15
)

correl <- corrDim(
  z,
  min.embedding.dim = 1,
  max.embedding.dim = 2*m,
  time.lag = d,
  min.radius = 0.001,
  max.radius = max(z),
  corr.order = 2,
  n.points.radius = 30,
  theiler.window = 3*d2,
  do.plot = FALSE,
  number.boxes = NULL)

plot(correl, type = "l", log = "xy", xlab = "lnr", ylab = "lnC(r)")
estimate(correl)
plotLocalScalingExp(correl, type = "l", main = "Correlation Integrals")

detach("package:nonlinearTseries", unload = TRUE)
library(fractal)
corrDim(z, dimension = 2*m, tlag=d, olag=3*d2, resolution=2)


#######
#  e  #

detach("package:fractal", unload = TRUE)
library(nonlinearTseries)

plotLocalScalingExp(correl, type = "l")

# Correlation dimensions:
c0 <- c(0.7, 0.9, 1.08, 1.2, 1.4, 1.55, 1.68, 1.79, 1.88, 1.97,
        2.05, 2.1, 2.16, 2.225, 2.28, 2.33, 2.38, 2.42, 2.44, 2.45)

plot(1:length(c0), c0,
     main = "Correlation Dim VS Embedding Dim",
     xlab = "Embedding Dimensions",
     ylab = "Correlation Dimensions",
     pch = 15)

fractal_dim <- estimate(correl); cat("fractal_dim =", fractal_dim, "\n")
min_dim <- ceiling(fractal_dim); cat("min_dim =", min_dim, "\n")
essential_embedding_dim <- m; cat("essential_embedding_dim =", essential_embedding_dim, "\n")


#######
#  f  #

detach("package:nonlinearTseries", unload = TRUE)
library(fractal)

lyap <- lyapunov(z,
               tlag=d,
               dimension=m,
               local.dimension = min_dim,
               olag=3*d2,
               metric= Inf,
               scale=NULL)
plot(lyap, type = "l")

corrDim(z, dimension = m, tlag=d, olag=3*d2, resolution=2)

#######
#  g  #

detach("package:fractal", unload = TRUE)
library(nonlinearTseries)

entropy <- sampleEntropy(correl,do.plot = TRUE)
estimate(entropy)




#####################################
##                B                ##
#####################################


library(nonlinearTseries)


z2 <- c()
for (i in (1:(len/4))){
  z2 <- c(z2, (z[4*(i-1) + 1] + z[4*(i-1) + 2] + z[4*(i-1) + 3] + z[4*(i-1) + 4])/4)
}
z <- z2
len <- length(z2)

#######
#  a  #

plot.ts(z)

#######
#  b  #

d <-timeLag(
  z,
  technique = "ami",
  selection.method = "first.minimum",
  lag.max = 15,
  do.plot = TRUE,
  main = "Mutual Information - First Minimum")
d

d2 <- timeLag(
  z,
  technique = "acf",
  selection.method = "first.e.decay",
  lag.max = 15,
  do.plot = TRUE,
  main = "Autocorrelation - 1/e")

d3 <- timeLag(
  z,
  technique = "acf",
  selection.method = "first.zero",
  lag.max = 15,
  do.plot = TRUE,
  main = "Autocorrelation - 0")


#######
#  c  #

qq <- as.data.frame(buildTakens(z, embedding.dim = 2, time.lag = d))
plot(qq[, 1], qq[, 2],
     pch = ".",
     main = "Phase Portrait - X, Y",
     xlab = "X",
     ylab = "Y")


#######
#  d  #

m <- estimateEmbeddingDim(
  z,
  number.points = len,
  time.lag = d,
  max.embedding.dim = 15
)

correl <- corrDim(
  z,
  min.embedding.dim = 1,
  max.embedding.dim = 2*m,
  time.lag = d,
  min.radius = 0.001,
  max.radius = max(z),
  corr.order = 2,
  n.points.radius = 30,
  theiler.window = 3*d2,
  do.plot = FALSE,
  number.boxes = NULL)

plot(correl, type = "l", log = "xy", xlab = "lnr", ylab = "lnC(r)")
estimate(correl)
plotLocalScalingExp(correl, type = "l", main = "Correlation Integrals")

detach("package:nonlinearTseries", unload = TRUE)
library(fractal)
corrDim(z, dimension = 2*m, tlag=d, olag=3*d2, resolution=2)


#######
#  e  #

detach("package:fractal", unload = TRUE)
library(nonlinearTseries)

plotLocalScalingExp(correl, type = "l")

# Correlation dimensions:
c0 <- c(0.5, 0.75, 0.9, 1.1, 1.2, 1.28, 1.36, 1.43, 1.5,
        1.56, 1.62, 1.67, 1.71, 1.74, 1.77, 1.79, 1.81, 1.82)

plot(1:length(c0), c0,
     main = "Correlation Dim VS Embedding Dim",
     xlab = "Embedding Dimensions",
     ylab = "Correlation Dimensions",
     pch = 16)

fractal_dim <- estimate(correl); cat("fractal_dim =", fractal_dim, "\n")
min_dim <- ceiling(fractal_dim); cat("min_dim =", min_dim, "\n")
essential_embedding_dim <- m; cat("essential_embedding_dim =", essential_embedding_dim, "\n")


#######
#  f  #

detach("package:nonlinearTseries", unload = TRUE)
library(fractal)

lyap <- lyapunov(z,
                 tlag=d,
                 dimension=m,
                 local.dimension = 3,
                 olag=3*d2,
                 metric= Inf,
                 scale=NULL)
plot(lyap, type = "l")

corrDim(z, dimension = m, tlag=d, olag=3*d2, resolution=2)


#######
#  g  #

detach("package:fractal", unload = TRUE)
library(nonlinearTseries)

entropy <- sampleEntropy(correl,do.plot = TRUE)
estimate(entropy)
