#############
# LIBRARIES #
#############

library(itsmr)
library(stats)



#############
#   INPUT   # 
#############

X <- uspop


#############
# SOLUTIONS #
#############

#####
# 1 #

plot.ts(X)


#####
# 2 #

# A) CLASSICAL DECOMPOSITION #

periodicity <- 1                                                         # INPUT

# Let X = s(seasonality) + m(trend) + Y(random noise)
s <- season(X, periodicity)
mY <- ts(X-s)
poly_rank <- 2                                                           # INPUT
m <- trend(mY,poly_rank)
Y <- ts(mY-m)

# B) DIFFERENCING METHOD #
Y2 <- diff(X, lag = periodicity)
Y2 <- diff(Y2, lag = poly_rank - periodicity)

#####
# 3 #

plot.ts(X)
plot.ts(Y)
plot.ts(Y2)


#####
# 4 #

# Selecting the Classical Decomposition Model
Y2 <- Y

p <- 2                                                                   # INPUT
q <- 1                                                                   # INPUT

###  AR Models  ###
# Yule-Walker
Yule_Walker <- yw(Y2, p)

# Burg
Burg <- burg(Y2, p)

###  General Models  ###
# ARMA
ARMA <- arma(Y2, p, q)

# Auto-fit
Autofit <- autofit(Y2, p = 1:p, q = 1:q)


#####
# 5 #

# Yule-Walker
print(Yule_Walker)

# Burg
print(Burg)

# ARMA
print(ARMA)

# Auto-fit
print(Autofit)


#####
# 6 #

c <- c("Yule_Walker", "Burg", "ARMA", "Autofit")
minmodel <- c[1]
minAICC <- Yule_Walker$aicc

for (i in c){
  model <- get(paste(i))
  if (minAICC <= model$aicc){
    #pass
  }
  else {
    minAICC <- model$aicc
    minmodel <- i}
}
cat("The model with the minimum AICC is the: ", minmodel, "\n", "AICC = ", minAICC, "\n", sep = "")

for (i in minmodel){
  model <- get(paste(i))
  phis <- c()
  phis2 <- c()
  thetas <- c()
  thetas2 <- c()
  for (j in (model[1])){
    phis <- c(phis, j)}
  for (k in (model[2])){
    thetas <- c(thetas, k)}
  for (j in (1:(length(phis)))){
    phis2 <- c(phis2, "(", phis[j], ") X_t-", j, " + ")}
  if(i == c[1] | i == c[2]){
    thetas2 <- 0}
  else{
    for (k in (1:(length(thetas)))){
      thetas2 <- c(thetas2, "(", thetas[k], ") Z_t-", k, " + ")}}
  cat("Model: ", i, "\n", sep = "")
  cat("X_t + ", phis2," = Z_t + ", thetas2, "\n", "{Zt} ~ WN(0, ", model$sigma2, ")", "\n", "\n", sep = "")
}

