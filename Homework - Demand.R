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

x <- read.table("clipboard" , sep = "\t", header = TRUE)
#x <- read.csv(path)

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
c0 <- c(0.996, 1.999, 2.186, 2.518, 2.428, 2.508, 2.404, 2.656, 2.551, 2.611,
        2.681)

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
c0 <- c(1.008, 1.965, 2.715, 2.388, 2.109, 2.291, 2.35, 2.523, 2.389,
        2.388, 2.381, 2.37, 2.479)


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
