####################
##  1. Libraries  ##
####################

library(fractal)
detach("package:fractal", unload = TRUE)
library(nonlinearTseries)


##############
#  2. Input  #
##############

path <- "C:/Users/thxsg/Downloads/w9.csv"
clm <- 1


############################
##  3. Setting Variables  ##
############################

#x <- read.table("clipboard" , sep = "\t", header = TRUE)
x <- read.csv(path, header = FALSE)

z <- x[[clm]]
len <- length(z)



# EXTRA: Time series
plot.ts(z)



####################
##  4. Solutions  ##
####################

#######
#  i  #

d <-timeLag(
  z,
  technique = "ami",
  selection.method = "first.minimum",
  lag.max = 15,
  do.plot = TRUE,
  main = "Mutual Information - First Minimum")


d2 <- timeLag(
  z,
  technique = "acf",
  selection.method = "first.e.decay",
  lag.max = 15,
  do.plot = TRUE,
  main = "Autocorrelation - 1/e")


paste("We select d =", d, ", as our proper time delay, because Mutual Information can identify Non-Linear Correlations as well:")


########
#  ii  #
paste("For Theiler's window we need to find the plateau in the space-time plot",
      "but for simplicity we can just set it equal to 3 times the time lag based on first 1/e value of the Autocorrelation graph")

m <- 5

correl <- corrDim(
  z,
  min.embedding.dim = 1,
  max.embedding.dim = 2*m,
  time.lag = d,
  min.radius = 0.001,
  max.radius = max(z)/2,
  corr.order = 2,
  n.points.radius = 30,
  theiler.window = 3*d2,
  do.plot = FALSE,
  number.boxes = NULL)

estimate(correl)
plotLocalScalingExp(correl, type = "l")

# Correlation Integrals
detach("package:nonlinearTseries", unload = TRUE)
library(fractal)
corrDim(z, dimension = 2*m, tlag=d, olag=3*d2, resolution=2)

# Dimensions
detach("package:fractal", unload = TRUE)
library(nonlinearTseries)
fractal_dim <- estimate(correl); cat("fractal_dim =", fractal_dim, "\n")
min_dim <- ceiling(fractal_dim); cat("min_dim =", min_dim, "\n")


#########
#  iii  #

detach("package:nonlinearTseries", unload = TRUE)
library(fractal)

lyap <- lyapunov(z,
                 tlag=d,
                 dimension = min_dim,
                 local.dimension = 3,
                 olag=3*d2,
                 metric= Inf,
                 scale=NULL)

plot(lyap, type = "l")

cat("Lyapunov Exponents:", "\n", "??1 = -0.11", "\n", "??2 = 0", "\n", "??3 = 0.02", "\n")
cat("Sum < 0, so its a chaotic system", "\n")


########
#  iV  #

detach("package:fractal", unload = TRUE)
library(nonlinearTseries)

entropy <- sampleEntropy(correl,do.plot = TRUE)
estimate(entropy)[min_dim]


#######
#  V  #

nonLinearPrediction(
  z,
  embedding.dim = min_dim,
  time.lag = 3*d2,
  prediction.step = 1,
  radius = 10^(-10),
  radius.increment = 10^(-3)
)

# Extra: Phase Portrait
qq <- as.data.frame(buildTakens(z, embedding.dim = 2, time.lag = d))
plot(qq[, 1], qq[, 2],
     pch = ".",
     main = "Phase Portrait - X, Y",
     xlab = "X",
     ylab = "Y")