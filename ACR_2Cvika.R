#' #Advanced Methods of Time Series Analysis Applied to Quarterly Estimates of Unemployment Rate#

#' ## Introduction ##
#' The chosen source of data is the Labour Force Survey (LFS) quarterly estimates of unemployment rate
#' in the UK since March 1971, up to March 2018.
#' 
#'---------------------------------------------------------------------------------------
#' 
#' ## 1. Elementary Modeling by an AR Process ## 
#' We begin by extracting the data from a downloaded file

setwd("D:/Martin Data/TimeSeriesAnalysis/ACR_2Repo")
dat <- read.table("LFS_unemployment.csv", header=F, skip=6, sep = ",", as.is = T)
dat <- dat[50:239,]
names(dat) <- c("quartile", "unemp")
dat$time <- dat$quartile[grepl('\\d{4}\\sQ\\d', dat$quartile)]

data.frame(head(dat$time)) # how the quarterly data looks

months <- (as.numeric(gsub("\\d{4}\\sQ(\\d)","\\1", dat$time)) - 1) * 3
data.frame(head(months))
years <- as.numeric(gsub("(\\d{4})\\sQ\\d","\\1", dat$time))
data.frame(head(years))
years <- years + months / 12
data.frame(head(years)) # how the time data should look

dat$time <- years

#' ### 1.1: Data Plot ###

```{r initDataPlot, echo=T, fig.width=10, fig.height=5}
plot(x=dat$time, y=dat$unemp, type="l", lwd=2, xlab="year", ylab = "unemployment [%]", main="Unemployment Rate (quarterly)")

#' Now even thought we are working with annual data there should not be any seasonal components or trend
#' since the results do not depend on periodic observable phenomena, but rather the complex economic 
#' situation over multiple decades. Also the data may include exponentially decaying decrease in unemployment,
#' but only after year 1980, which would suggest a regime-switching stochastic process. 
#' 
#' ### 1.2: Test Part and Evaluation Part of the Time Series ###
#' 
#' Now we separate the time series into test part, where a suitable model of a stochastic process will
#' be found, and the evaluation part, where predictions given by such model are evaluated. Since our
#' dataset contains quarterly data, we choose the length of the evaluation part of the time series as 
#' L = k * 4 where k is an arbitrary (and sufficiently small) positive integer.

k = 5 # how many years
( L = k * 4 )
N = length(dat$unemp)
unempseries <- list(unempts = ts(dat$unemp))

unempseries$test = window(unempseries$unempts, start=1, end=N - L, extend=T)
unempseries$eval = window(unempseries$unempts, start=N - L + 1, end=N, extend=T)

```{r testEvalPlot, echo=T, fig.width=9, fig.height=4}
par(mfrow=c(1,1))
plot(x=dat$time, y=dat$unemp, main="Test and Evaluation Part of the Quarterly Series", xlab="year", ylab="unemployment [%]",type="l")
lines(head(dat$time, n = N - L), unempseries$test, col="red", lwd=2)
lines(tail(dat$time, n = L), unempseries$eval, col="blue", lwd=2)
legend("topright", legend=c("test","eval"), col=c("red","blue"), lty=1, lwd=2, cex=0.8)

#' ### 1.3: Mean, Variance, ACF, and PACF of the Test Part ###

stat <- summary(as.numeric(dat$unemp)[1:(N - L)])
stddev <- var(unempseries$test, na.rm = TRUE)
v <- as.numeric(c(stat[1], stat[6], stat[4], stat[3], stddev))
names(v)<-c("Min.","Max.","Mean","Median","Variance")
v

```{r ACFplot1, fig.width=9, fig.height=4}
par(mfrow=c(1,2))
acf(as.numeric(unempseries$test), lag.max=N, main="unemp. ACF")
acf(as.numeric(unempseries$test), lag.max=N, type="partial",main="unemp. PACF")

#' As we mentioned in section 1.1, the underlying process which gave rise to the observed results is aperiodic,
#' yet it is undoubtedly a process with memory. Unemployment rate strongly depends (aside from other important aspects)
#' on its own history which might extend generations into the past. The results are, however, significantly influenced
#' by external phenomena, such as the global economic crisis in late 2000's.
#'
#' ### 1.4: Finding a Suitable AR Model ###
#' 
#' Since the economic situation and the job market remembers its past, we choose a simple AR(p) process with parameter p
#' corresponding to the number of steps after which the process still "remembers" its past. 
#' The inbuilt ar() function automatically finds the model with the lowest AIC (Akaike's Information Criterion). 
#' And by plotting the $aic parameter we obtain differences AIC_min - AIC_k for all models.

kmax = 30 # set a maximum order

model <- list()
model$ar<- ar(as.numeric(unempseries$test), order.max=kmax)

```{r AICplot1, fig.width=10, fig.height=4}
par(mfrow=c(1,2))
plot(0:(length(model$ar$aic)-1), xlab="p", model$ar$aic, ylab="dAIC", main="Differences in AIC")
tmp <- sapply(1:kmax, function(x) ar(as.numeric(unempseries$test), aic=F, order.max=x)$var.pred)
plot( tmp, xlab="p", ylab="sigma^2", main="Residual variances")

#' As we can see in the figures, the lowest variance of residues corresponds to an AR(3) process
#' with coefficients:

rbind(coef=model$ar$ar, se=sqrt(diag(model$ar$asy.var.coef))) 

( p <- model$ar$order )

#' Unfortunately, the `ar` function does not return fitted values, thus we need to model the time series via the
#' `arima` function using the AR order `p` from the previous result.

suppressMessages(library(forecast))
( model$AR <- arima(as.numeric(unempseries$test), order=c(p, 0, 0)) )
unempseries$ARfit <- fitted(model$AR)
unempseries$ARresid <- residuals(model$AR)

```{r FitARPlot, fig.width=10, fig.height=5}
par(mfrow=c(1,2))
plot(x=dat$time[1:(N - L)], y=unempseries$test, 
     main="Fitted Values of the Test Part", 
     xlab="year", ylab="unemployment [%]",type="p")
lines(x=dat$time[1:(N - L)], y=unempseries$test, col="blue", lwd=2)
legend("topright", legend=c("AR"), col=c("blue"), lty=1, lwd=2, cex=0.8)
plot(x=dat$time[1:(N - L)], y=unempseries$ARresid, main="AR Residuals", xlab="year", ylab="unemployment [%]", type="l")

#' The given AR model seems to fit the time series very well, which may be due to its low oscillation rate.
#' 
#' ### 1.5: 1-Step Predictions Over the Evaluation Part ###
#' 

tmp <- Arima(as.numeric(dat$unemp), model = model$AR)
model$ARpred <- window(fitted(tmp), start = N - L + 1)

( RMSE <- sqrt(mean((as.numeric(unempseries$eval) - as.numeric(model$ARpred))^2)) )

```{r PredictionsARPlot, fig.width=10, fig.height=6}
par(mfrow=c(1,1))
plot(x=dat$time[(N - L + 1):N], y=dat$unemp[(N - L + 1):N], 
     main=paste("Prediction AR(",p,"), RMSE = ",round(RMSE, digits = 4)), 
     xlab="year", ylab="unemployment [%]",type="p", pch=20)
lines(x=dat$time[(N - L + 1):N], y=model$ARpred, lwd=2.5, col="blue")

#' ### 1.6.: Conclusion ###
#' 
#' Since it has very low rate of local oscillation, but does not have an easily predictable systematic pattern,
#' the analyzed unemployment rate time series seems to be well-estimated by an AR(3) process with low prediction errors.
#' However, as we mentioned earlier we might be dealing with a 'regime-switching' process. 
#' Further analysis will be carried out in the next chapter.
#' 
#'  
#'---------------------------------------------------------------------------------------
#' 
#' ## 2. Finding the Parameters of a SETAR Model ## 
#' 
#' As we mentioned, the unemployment time series might be a result of a regime-switching process. Naturally, the behavior
#' of the unemployment rate in a given country should depend on the current economic situation. The change in the local economy
#' can be described via a set of "thresholds" which determine whether the stochastic process changes its regime. The regime of
#' a stochastic process is defined as a unique ARMA or any other linear stochastic process with unique parameters. We begin 
#' by finding the parameters of a Self-Exciting Threshold Autoregressive (SETAR) process, that is: a process whose regime 
#' is described by a random variable determined by the very process itself, more specifically its history of up to d steps behind,
#' which in an essence means that the process "influences its regime" up to d time steps into the future.
#' 
#' For the purposes of this analysis we consider only 2 regimes, namely the regime of "job crisis" when the unemployment rate may 
#' fluctuate or drop more wildly compared to the regime of "job stability" when the unemployment rate stabilizes or grows.
#' 
#' ### 2.1: Implementation of Useful Functions ###
#' 
#' First we define a, so called, "indicator function" which essentially returns a boolean value from a given input process
#'  value `x` and threshold value `c`:

Indicator <- function(x, c) ifelse(x > c, 1, 0)

#' Afterwards, we define the 'basis' function for a single regime

Yt <- function(x, t, p) c(1, x[(t - 1):(t - p)])

#' which can then be used in the "basis" for two regimes:

Xt <- function(x, t, p, d, c, z = x) {  
  # z is the threshold variable 
  I <- Indicator(z[t - d], c)
  Y <- Yt(x, t, p)
  c((1 - I) * Y, I * Y)
}

#' We can test the function on a given subset of the unemployment time series. Due to the fact that the
#' examined time series is rather 'smooth', for further use, we will examine its differences:

xt <- diff(as.numeric(dat$unemp))

d <- 1; t <- 3;
xt[(t - d): t]
Xt(xt, t, p, d, c = 0)

t <- 100;
xt[(t - d): t]
Xt(xt, t, p, d, c = 0)

#' As we can see, in the first case, with time series values crossing zero from above, the latter half of the
#' coefficient vector gets expressed, corresponding to the series assuming the second regime.
#' 
#' Then we need a function defining a deterministic skeleton of the model:

SkeletonSETAR <- function(x, t, p, d, c, theta, z = x) theta %*% Xt(x, t, p, d, c, z)

#' where `theta` corresponds to the parameter vector, for example:

SkeletonSETAR(xt, t, p, d, c = 0, c(1, 1, 1, 1, 1, 1, 1, 1))

#' and the last group of functions we need for the upcoming procedure are functions for the information criteria
#' of a SETAR model:

# Akaike
AIC_SETAR <- function(orders, regimeDataCount, resVariances) {
  sum(regimeDataCount * log(resVariances) + 2 * (orders + 1))
}

# Bayesian
BIC_SETAR <- function(orders, regimeDataCount, resVariances) {
  sum(regimeDataCount * log(resVariances) + log(regimeDataCount) * (orders + 1))
} 

#' and test it out:

AIC_SETAR(c(2, 2), c(10, 10), c(0.5, 0.7))
BIC_SETAR(c(2, 2), c(10, 10), c(0.5, 0.7))

#'
#' ### 2.2: The Estimation of Parameters of a SETAR Model  ###
#' 
#' Given a dataset `x` and parameters `p` (AR order), `d` (SETAR delay), and the threshold `c` we find the coefficients
#' of a SETAR model with these parameters by performing a multivariate linear regression. The coefficient vector`PhiParams`
#' is the vector of unknowns of a linear system with matrix `X` and a right-hand-side vector `y` 
#' given by the time series. Although for higher values of `p` the inversion of matrix `X'X` (with dimensions 
#' `(2*p + 2)`x`(2*p + 2)`) might be computationally demanding, we will determine the covariance matrix, i.e.: `(X'X)^(-1)`
#' using a function `inv` from the `matlib` package:

suppressMessages(library(matlib))

EstimSETAR <- function(x, p, d, c) {
  resultModel <- list()
  resultModel$p = p; resultModel$d = d; resultModel$c = c;
  resultModel$data = x;  n = length(x);  resultModel$n = n;
  k <- max(p, d)
  
  X <- as.matrix(apply(as.matrix((k + 1):n), MARGIN=1, function(t) Xt(x, t, p, d, c) ))
  y <- as.matrix(x[(k + 1):n])
  
  A = crossprod(t(X), t(X));  b = crossprod(t(X), y)

  if (abs(det(A)) > 0.000001) {
    resultModel$PhiParams <- solve(A, b) # solving (X'X)*phi = X'y
    resultModel$PhiStErrors <- sqrt(diag(inv(A)) / n)  # standard errors
    skel <- crossprod(X, resultModel$PhiParams)
    resultModel$residuals <- (y - skel)
    resultModel$resSigmaSq <- 1 / (n - k) * sum(resultModel$residuals ^ 2)
    
    resultModel$n1 <- sum(apply(y, MARGIN = 1, function(x) (1 - Indicator(x, c))))
    resultModel$n2 <- sum(apply(y, MARGIN = 1, function(x) Indicator(x, c)))
    
    resultModel$resSigmaSq1 <- sum(apply(y, MARGIN = 1, function(x) ifelse((1 - Indicator(x, c)),(x - skel)^2, 0))) / resultModel$n1
    resultModel$resSigmaSq2 <- sum(apply(y, MARGIN = 1, function(x) ifelse(Indicator(x, c),(x - skel)^2, 0))) / resultModel$n2
    
    return(resultModel)
  } else {
    return(NA)
  }
}

#' and now we test the function for suitable parameters:

str(EstimSETAR(xt, 2, 1, c=0))


#' It should be noted that for some values of `p` and `d` the indices of arrays in the algorithms might 
#' get out of range. For that reason we implement exceptions for the outputs of `EstimSETAR` in the following
#' algorithm.
#' 
#' ### 2.3: SETAR Parameter Estimation Procedure  ###
#' 
#' To answer the question: 'how does one find the right parameters `p`, `d` and `c` for their desired SETAR model?',
#' we implement the following procedure:
#' 
#' 

pmax <- 7 # set maximum order p
# limit the c parameter by the 7.5-th and 92.5 percentile
cmin <- as.numeric(quantile(xt, 0.075)); cmax <- as.numeric(quantile(xt, 0.925));
h = (cmax - cmin) / 100 # determine the step by which c should be iterated

models <- list()
modelColumns <- list()
for (p in 1:pmax) {
  for (d in 1:p) {
    pdModels <- list()
    for (c in seq(cmin, cmax, h)) {
      tmp <- EstimSETAR(xt, p, d, c) # try to run the function
      # then test whether it returns`NA` as a result
      if (!as.logical(sum(is.na(tmp))) ) {
        pdModels[[length(pdModels) + 1]] <- tmp
      }
    }
    sigmas <- as.numeric(lapply(pdModels, function(m) m$resSigmaSq))
    orders <- order(sigmas)
    min_sigma_model <- pdModels[[ orders[1] ]]
    # only the model whose parameter c gives the lowest residual square sum is chosen
    min_sigma_model$AIC <- AIC_SETAR(c(p, p), c(min_sigma_model$n1, min_sigma_model$n2), c(min_sigma_model$resSigmaSq1, min_sigma_model$resSigmaSq2))
    min_sigma_model$BIC <- BIC_SETAR(c(p, p), c(min_sigma_model$n1, min_sigma_model$n2), c(min_sigma_model$resSigmaSq1, min_sigma_model$resSigmaSq2))
    models[[length(models) + 1]] <- min_sigma_model
    modelColumns[[length(modelColumns) + 1]] <- c(
      p, d, min_sigma_model$c,
      min_sigma_model$n1, min_sigma_model$n2,
      min_sigma_model$AIC, min_sigma_model$BIC,
      min_sigma_model$resSigmaSq)
  }
}
modelColumns <- data.frame(matrix(unlist(modelColumns), nrow=length(modelColumns), byrow=T))
names(modelColumns) <- c(
  "p", "d", "c",
  "n1", "n2", "AIC", "BIC",
  "resSigmaSq"
)

head(modelColumns, n=12)


#' Now we have a set of models in their original order. To find the best suitable model, we choose
#' 12 models with the lowest BIC (Bayesian Information Criterion):

BICs <- sapply(models, function(m) m$BIC)
orders <- order(BICs)
modelColumns <- modelColumns[orders,]

head(modelColumns, n=12)

#' and we can also include errors of the estimated regression coefficients:

modelCoeffErrors <- list()
for(i in 1:12) {
  p <- models[[ orders[i] ]]$p
  d <- models[[ orders[i] ]]$d
  c <- models[[ orders[i] ]]$c
  key <- paste(p, d, round(c, digits=4), sep="/")
  modelCoeffErrors[[key]] <- rbind(t(models[[ orders[i] ]]$PhiParams), t(models[[ orders[i] ]]$PhiStErrors))
  row.names(modelCoeffErrors[[key]]) <- t(c("Phi", "stdError"))
}
modelCoeffErrors

#' We can now visualize the results of the top 3 models:

```{r SETARTop3Plot, fig.width=9, fig.height=4}
plotNmax <- 75
par(mfrow=c(1,2))
for (i in 1:3) {
  model <- models[[orders[i]]]
  SetarFit <- xt - append(matrix(0., ncol=model$p), model$residuals)
  m <- length(model$residuals)

  plot(x=dat$time[1:plotNmax], y=xt[1:plotNmax],
       main=paste("SETAR(",model$p,",",model$d,",",round(model$c, digits=4),")"), xlab="year", ylab="%")
  lines(x=dat$time[1:plotNmax], y=SetarFit[1:plotNmax], col="blue",lwd=2)
  legend("topleft", legend=c("fitted SETAR"), col=c("blue"), lty=1, lwd=2, cex=0.8)
  
  plot(x=dat$time[1:plotNmax], y=model$residuals[1:plotNmax], type="l", main=paste("SETAR(",model$p,",",model$d,",",round(model$c, digits=4),") Residuals"),
       xlab="year", ylab="%")
}

#' The results suggest that the best SETAR models have a threshold `c` quite close to zero and more-or-less the same
#' RSS. The very first `SETAR(1,1, 0.005)`, having a significanly lower `BIC` (Bayesian Information Criterion), has slightly
#' larger RSS than the models that come after it. To verify the correctness of our procedure we will need to compare it
#' with inbuilt functions from a verified library.
#' 
#' ### 2.4: SETAR Equilibrium Simulations ###
#'  
#'  It is also essential to find out whether the skeletons of the selected SETAR models have some equilibria. 
#'  We are going to find out by simulating trajectories of the top 3 models selected in the previous section:

par(mfrow=c(1,1))

```{r SETARSkeletTop3Data, fig.width=9, fig.height=4}
plot(x=dat$time[1:plotNmax], y=xt[1:plotNmax],
     main="Skeletons of Chosen Models With the Provided Data", xlab="year", ylab="%")
for (i in 1:3) {
  s <- array()
  k <- max(models[[orders[i]]]$p, models[[orders[i]]]$d)
  for (t in (k + 1):plotNmax) {
    tmp <- try(SkeletonSETAR(xt, t, models[[orders[i]]]$p, models[[orders[i]]]$d, models[[orders[i]]]$c, t(models[[orders[i]]]$PhiParams)),
               silent=T)
    if (class(tmp) == "try-error") {
      s[t - k] <- NA
    } else {
      s[t - k] <- tmp
    }
  }
  lines(x=dat$time[(k + 1):plotNmax], y=s[1:(plotNmax - k)], col=c("blue","purple","red")[i],lwd=2)
}
legend("topleft", legend=c(paste("SETAR(",models[[orders[1] ]]$p,",",models[[orders[1] ]]$d,",",round(models[[orders[1] ]]$c, digits=4),")"),
                           paste("SETAR(",models[[orders[2] ]]$p,",",models[[orders[2] ]]$d,",",round(models[[orders[2] ]]$c, digits=4),")"),
                           paste("SETAR(",models[[orders[3] ]]$p,",",models[[orders[3] ]]$d,",",round(models[[orders[3] ]]$c, digits=4),")")), 
       col=c("blue","purple","red"), lty=1, lwd=2, cex=0.8)


```{r SETARsimulations, fig.width=9, fig.height=4}
xmax <- 5 * max(xt); xmin <- min(xt) - 0.5 * xmax;
xmax0 <- max(xt); xmin0 <- min(xt)
par(mfrow=c(1,2))
for (i in 1:3) {
  equilib_sims <- list()
  p <- models[[orders[i] ]]$p; d <- models[[orders[i] ]]$d; c <- models[[orders[i] ]]$c;
  sigmaSq <- models[[orders[i] ]]$resSigmaSq
  k <- max(p, d)
  n_offsets <- 10
  for (j in 0:n_offsets) {
    x0 <- (xmax - xmin) / n_offsets * j - 0.5 * (xmax - xmin)
    equilib_sim <- array()
    equilib_sim[1] <- x0
    for (t in 2:plotNmax) {
      if (t < (k + 1)) {
        equilib_sim[t] <- equilib_sim[t - 1] + rnorm(1, 0, sqrt(sigmaSq))
      } else {
        equilib_sim[t] <- SkeletonSETAR(equilib_sim, t, p, d, c, t(models[[orders[i]]]$PhiParams)) + rnorm(1, 0, sqrt(sigmaSq))
      }
    }
    equilib_sims[[j + 1]] <- equilib_sim
    if (j < 1) {
      plot(x=dat$time[1:plotNmax], y=equilib_sim, type="l", col="gray", ylim=c(1.5 * xmin, 1.5 * xmax), 
           main=paste("SETAR(",models[[orders[i] ]]$p,",",models[[orders[i] ]]$d,",",round(models[[orders[i] ]]$c, digits=4),")"),
           xlab="year", ylab="%")
    } else {
      lines(x=dat$time[1:plotNmax], y=equilib_sim, col="gray")
    }
  }
  for (j in 0:n_offsets) {
    epsilon <- (xmax - xmin) * 0.025
    x0 <- c + (n_offsets / 2 - j) * epsilon #small initial perturbations from the threshold value
    equilib_sim <- array()
    equilib_sim[1] <- x0
    for (t in 2:plotNmax) {
      if (t < (k + 1)) {
        equilib_sim[t] <- equilib_sim[t - 1] + rnorm(1, 0, sqrt(sigmaSq))
      } else {
        equilib_sim[t] <- SkeletonSETAR(equilib_sim, t, p, d, c, t(models[[orders[i]]]$PhiParams)) + rnorm(1, 0, sqrt(sigmaSq))
      }
    }
    equilib_sims[[j + 1]] <- equilib_sim
    lines(x=dat$time[1:plotNmax], y=equilib_sim, lwd=2)
    lines(x=c(dat$time[1], dat$time[plotNmax]), y=c(c, c), col="green", lty="dashed", lwd=2)
  }
  # comparison of the original data with the simulation of the original time series
  
  x0 <- xt[1]
  equilib_sim <- array()
  equilib_sim[1] <- x0
  for (t in 2:plotNmax) {
    if (t < (k + 1)) {
      equilib_sim[t] <- equilib_sim[t - 1] + rnorm(1, 0, sqrt(sigmaSq))
    } else {
      equilib_sim[t] <- SkeletonSETAR(equilib_sim, t, p, d, c, t(models[[orders[i]]]$PhiParams)) + rnorm(1, 0, sqrt(sigmaSq))
    }
  }
  simMax <- max(equilib_sim, xmax0); simMin <- min(equilib_sim, xmin0);
  plot(x=dat$time[1:plotNmax], y=xt[1:plotNmax], type="l", lwd=1.5, ylim=c(simMin, 1.1 * simMax), 
       main=paste("Simulation SETAR(",models[[orders[i] ]]$p,",",models[[orders[i] ]]$d,",",round(models[[orders[i] ]]$c, digits=4),")"),
       xlab="year", ylab="%")
  lines(x=dat$time[1:plotNmax], y=equilib_sim, col="purple", lwd=2, )
  lines(x=c(dat$time[1], dat$time[plotNmax]), y=c(c, c), col="green", lty="dashed", lwd=2)
  legend("topleft", legend=c("original ts","simulated SETAR","threshold"), col=c("black","purple","green"), lty=c(1,1,2), lwd=2, cex=0.8)
}

#' The trajectories of all of the first three models seem to gravitate toward `0` significantly fast (or alternatively:
#' towards their threshold values which are close to zero as well). The relatively low oscillation rate of the original 
#' time series suggests that the differences of this time series will, at most, fluctuate around 0. The change between 
#' the 'high' and 'low' regimes does not seem very significant, at leat on the larger scale. The validity of the model 
#' will be tested in chapter 3.
#' 
#' ### 2.5: Comparison Of the Results With Inbuilt Functions ###
#' 
#' To verify the correctness of our methods we proceed to construct the top 3 SETAR models by plugging their parameters into 
#' inbuilt functions:
#' 

suppressMessages(library(tsDyn))

#' Testing a function which selects an orders automatically:
mmax <- 2
```{r selectSetar1, fig.width=10, fig.height=4}
par(mfrow=c(1,1))
( result1 <- selectSETAR(xt, m=mmax, thDelay=0:(mmax-1), criterion="BIC", same.lags=T, trim=0.1)  )
#' the estimated thDelay corresponds to d-1. 
tmp <- EstimSETAR(xt, 2, 1, -0.1)
str(tmp)

#' Note that we set `thDelay=0:(mmax-1)` instead of `1:mmax`. `selectSETAR` uses `thDelay = 0` for step `d=1` delay
#' correspondence: `x(t-d) < c` or `x(t-d) > c`.
#' the resulting BIC's are different, possibly due to the package using a different formula
mmax <- pmax
```{r selectSetar2, fig.width=10, fig.height=10}
par(mfrow=c(1,1))
( result2 <- selectSETAR(xt, m=mmax, thDelay=0:(mmax-1), criterion="BIC", same.lags=T, trim=0.1) )

#' Setting higher `mmax`, the function returns a list of models similar to the one given by our procedure in section 2.3.
#' We can also compare the accuracy of the computation of the regression coefficients in our `estimSETAR` method, with
#' for example: `setar()` function (from `tsDyn` library as well):

setars <- list()
coeffComparison <- list()
resSigmaComparison <- list()
n <- length(xt)
for (i in 1:3) {
  p <- models[[ orders[i] ]]$p
  d <- models[[ orders[i] ]]$d
  c <- models[[ orders[i] ]]$c
  setars[[i]] <- setar(xt, m=p)
  k <- max(p, d)
  resSigmaComparison[[i]] <- c(
    (1 / (n - k) * sum(setars[[i]]$residuals ^ 2)), 
    models[[ orders[i] ]]$resSigmaSq
    )
  inbuiltParams <- t(setars[[i]]$coefficients)
  key <- paste(p, d, round(c, digits=4), sep=" / ")
  coeffComparison[[key]] <- rbind(
    inbuiltParams, 
    t(append(models[[orders[i]]]$PhiParams, models[[orders[i]]]$c))
  )
  row.names(coeffComparison[[key]]) <- t(c("inbuilt", "custom"))
}
coeffComparison

# comparing RSS

resSigmaComparison <- data.frame(matrix(unlist(resSigmaComparison), nrow=3, byrow=T))
colnames(resSigmaComparison) <- c("inbuilt", "custom")
row.names(resSigmaComparison) <- t(paste("rss",1:3))
resSigmaComparison

#' Without specifying the threshold value, the inbuilt `setar` function finds threshold values quite close
#' to those of our custom procedures. The AR order `p` for both regimes, however, has to be specified in advance.
#' The comparison of the model coefficients suggests that our custom method was more-or-less accurate. 
#'
#' ### 2.6: Conclusion ###
#'
#'  The results of the SETAR Parameter Estimation Procedure in section 2.3 show that the 3 best 2-regime SETAR
#'  models are `SETAR(1,1,0.005)`, `SETAR(2,2,0.1077)`, and `SETAR(2,1,0.4)`. The first model with the lowest 
#'  `BIC` (Bayesian Information Criterion) has the most accurate estimation of its 4 regression parameters, 
#'  with the highest residual square sum. The first model seems to have a stable equilibrium at their threshold values.
#'  
#'  --------------------------------------------------------------------------------------------------------------
#'  
#' ## 3: Tests of Linearity/Nonlinearity of SETAR models ##
#'
#'  We need to make sure a non-linear model (SETAR, for example) is really suitable for describing the process. 
#'  In order to find out, we test the null hypothesis that a linear model is more suitable than a non-linear one. 
#'  In the case of a 2-regime model we are looking for, so called, nuisance parameters, i.e.: `H0: Phi1 = Phi2` where
#'  `Phi1` and `Phi2` are the parameters of the low and the high regime respectively.
#'  
#' ### 3.1: Hansen's Conditions ###
#' 
#' Hansen proposed three conditions to test whether a SETAR model can be tested for linearity using the so called
#' Likelihood-Ratio (LR) test:

Hansen <- function(d, c, Phi) {
  p <- (length(Phi)/2) - 1  
  Phi <- do.call(rbind, split(Phi,rep(1:2,each=(p + 1))))  #separate regimes into rows
  c1 <- !isTRUE(all.equal( 0, apply(Phi[,c(1,1 + d),drop=F],2, diff) %*% c(1,c) ))  # (p10-p20)+(p1d-p2d)*c <= 0
  c2 <- all(apply(Phi[,-c(1,1 + d),drop=F], 2, function(x) !identical(0, diff(x))))  # p1j neq p2j, j notin {0,d}
  c3 <- all(apply(Phi[,-1,drop=F], 1, function(x) sum(abs(x))) < 1)  # sum_j|pij| < 1 forall i=1,2
  c(cond1=c3, cond2=c2, cond3=c3)
}

#' If all three are satisfied the model can be tested using the LR test:

hansenResults <- matrix(NA, ncol=3, nrow=12)
keys <- array()
for (i in 1:12) {
  p <- models[[ orders[i] ]]$p
  d <- models[[ orders[i] ]]$d
  c <- models[[ orders[i] ]]$c
  keys[i] <- paste(p, d, round(c, digits=4), sep=" / ")
  hansenResults[i,] <- Hansen(d, c, models[[ orders[i] ]]$PhiParams)
}
row.names(hansenResults) <- keys
colnames(hansenResults) <- c("cond1","cond2","cond3")

hansenResults

#' It appears that only the first and the fifth model can be tested using the LR test. The rest will have to be 
#' assessed using the Lagrange Multiplier (LM) test.
#' 
#' ### 3.2: LR and LM Tests ###
#' 
#' In this section we formulate the basic procedures for the LR (Likelihood Ratio), and LM (Lagrange Multiplier) tests:

LRtest <- function(x, p, var, alpha=0.05) {
  tmp <- ar(x, aic=F, order.max=pmax, method = "ols")
  tmp <- tmp$var.pred  # linear model residual variance
  testat <- length(x)*(tmp-var)/tmp  # test statistic
  CDF <- Vectorize( function(t) {  # test statistic CDF
    fun <- function(t) 1 + sqrt(t/(2*pi))*exp(-t/8) + 1.5*exp(t)*pnorm(-1.5*sqrt(t)) - (t+5)*pnorm(-sqrt(t)/2)/2
    if(abs(t)>300 || is.infinite(t)) return(sign(t))
    if(t >=0 ) fun(t) else 1-fun(-t)
  })
  if(alpha==0.05) critval <- 7.68727553 # for alpha=2.5%: CV=11.03329250
  else critval <- uniroot(function(x) CDF(x) - (1-alpha), c(-1000,1000))$root
  c(TS=testat, CV=critval, p_value=1-CDF(testat))  # (test statistics, critical value, p-value)
}

LRtest(xt, models[[ orders[1] ]]$p, models[[ orders[1] ]]$resSigmaSq)

#' --------------------------------------------------------------------------------------------------------

suppressMessages(library(dynlm))

LMtest <- function(x, p, d, alpha = 0.05) {
  names(p) <- NULL   # prevent from passing (accidental and needless) name to result
  x <- as.ts(x)  # if x is not a ts object, by chance
  model1 <- dynlm(x ~ L(x,1:p))  # requires dynlm package (it can be implemented withou dynlm, see model2)
  y <- model1$residuals
  tmp <- c(  # a list of shifted time series
    list(y),
    lapply(1:p, function(i) stats::lag(x, -i)), 
    lapply(1:p, function(i) stats::lag(x, -i)*stats::lag(x,-d)),
    list(stats::lag(x,-d)^3)
  )
  tmp <- do.call(function(...) ts.intersect(..., dframe=T), tmp)
  names(tmp) <- c("y", paste0("x",1:p), paste0("xd",1:p), "xd^3")
  model2 <- lm(y ~ ., data = tmp)  # cannot be done with the dynlm package
  z <- model2$residuals
  testat <- (length(x)-p) * (sum(y^2)/sum(z^2) - 1)
  c(TS=testat, CV=qchisq(1-alpha, df=p+1), p_value=1-pchisq(testat, df=p+1))
}

LMtest(xt, models[[ orders[1] ]]$p, models[[ orders[1] ]]$d)

#' We can easily automate the testing procedure in the following loop:
alpha = 0.05
results <- list()
nonlinear <- list()
nonLinCount <- 0
for (i in 1:12) {
  p <- models[[ orders[i] ]]$p; d <- models[[ orders[i] ]]$d; c <- models[[ orders[i] ]]$c;
  hansenResult <- Hansen(d, c, models[[ orders[i] ]]$PhiParams)
  if (FALSE %in% hansenResult) {
    hansenResult <- FALSE
    testResult <- LMtest(xt, p, d)
  } else {
    hansenResult <- TRUE
    testResult <- LRtest(xt, p, models[[ orders[i] ]]$resSigmaSq)
  }
  if (testResult[3] < alpha) {
    nonLinCount <- nonLinCount + 1
    nonlinear[[nonLinCount]] <- models[[ orders[i] ]]
  }
  results[[i]] <- append(cbind(p, d, c, hansenResult), testResult)
}
results <- data.frame(matrix(unlist(results), nrow=12, byrow=T))
colnames(results) <- c("p", "d", "c", "Hansen Cond.", "TS", "CV", "p-value")
row.names(results) <- orders[1:12]
results[,4] <- as.logical(results[,4])

results

#' ### 3.3 Modified LR Test Via Boostrapping ###
#' 
#' The proposed LR test has a significant drawback in the fact that it can only be done when Hansen's 
#' conditions are satisfied. This is due to the fact that we do not know the distribution of the resulting 
#' F-statistic. According to Hansen (1996), however, the distribution of a bootstrapped statistic F* converges
#' weakly in probability to the distribution of F, so that repeated bootstrap draws from F* can be used to
#' approximate the asymptotic distribution of F.
#' 
#' To obtain F* we will make regressions of the time series `xt` and the multi-regime matrix `Xt`,
#' where `y_star` will be the vector of dependent variables - a series of random values generated from N(0,1).
#' 

```{r bootstrap1, fig.width=8, fig.height=4}

bootstrapSigmaSq1 <- function(x, y_star) {
  tilde_model <- lm(y_star ~ x, data=data.frame(x))
  return (sum((tilde_model$residuals - x)^2) / length(x))
}

bootstrapSigmaSq2 <- function(x, y_star, p, d, c) {
  n <- length(x); k <- max(p, d);
  Yc <- matrix(0., ncol = (2 * p + 2), nrow = (2 * p + 2)) # covariance matrix
  Y <- matrix(0., ncol = (2 * p + 2), nrow = n) # matrix of basis vectors
  Xc <- matrix(0., nrow = (2 * p + 2)) # rhs vector
  i <- 1 
  for (t in (k + 1):n) {
    XT <- Xt(x, t, p, d, c); Y[i,] <- XT;
    Yc <- Yc + (XT %o% XT) # accumulating elements of the cov matrix
    Xc <- Xc + XT * y_star[t] # using random draws as rhs
    i <- (i + 1)
  }
  det <- det(Yc)
  if (det > -0.00001 && det < 0.00001) {
    return(NA) # return NA if Yc is (almost) singular
  } else {
    params <- inv(Yc) %*% Xc
    z <- array()
    for (t in 1:(n - k)) {
      z[t] <- crossprod(Y[t,], params)
    }
   residuals <- (y_star[1:(n - k)] - z)
    return(sum(residuals ^ 2)/ (n - k) )
  }
}

n <- length(xt)
n_draws <- 1000 # of bootstrap draws per threshold c
model_p_values <- list();

for (i in 1:12) {
  p <- models[[ orders[i] ]]$p
  d <- models[[ orders[i] ]]$d
  c <- models[[ orders[i] ]]$c
  sigmaSq <- models[[ orders[i] ]]$resSigmaSq
  Fstat <- LRtest(xt, p, sigmaSq)[1]
  draw_count <- 0
  for (j in 1:n_draws) {
    y_star <- rnorm(n)
    sigmaSqTilde <- bootstrapSigmaSq1(xt, y_star)
    sigmaSqHat <- bootstrapSigmaSq2(xt, y_star, p, d, c)
    FstatDraw <- (n * (sigmaSqTilde - sigmaSqHat) / sigmaSqHat)
    if (FstatDraw > Fstat) {
      draw_count <- draw_count + 1
    }
  }
  p_val <- draw_count / n_draws
  model_p_values[[i]] <- c(paste("(",p,",",d,",",round(c, digits=4),")"), p_val, 
                           ifelse(p_val < 0.1 && p_val >= 0.05, ".",
                                  ifelse(p_val < 0.05 && p_val >= 0.01, "*",
                                         ifelse(p_val < 0.01, "**", ""))))
}

bootstrap_results <- data.frame(matrix(unlist(model_p_values), nrow=12, ncol=3, byrow=T))
colnames(bootstrap_results) <- c("model","P-value","")
bootstrap_results

#' ### 3.4 Visualisation of Non-Linear Models ###
#' 
#' From the results of the previous procedure, we will visualize the models for which the linearity
#' null-hypothesis was rejected based on the LR and LM tests:

```{r trueSETARs, fig.width=10, fig.height=4}
par(mfrow=c(1,2))
xmax <- max(xt); xmin <- min(xt) - 0.5 * xmax;
xmax <- xmax + 0.2 * (xmax - xmin)
for (i in 1:nonLinCount) {
  p <- nonlinear[[i]]$p
  d <- nonlinear[[i]]$d
  c <- nonlinear[[i]]$c
  x0 <- as.numeric(unempseries$test)[1]
  fitted <- xt[1:(N - L)] - nonlinear[[i]]$residuals[1:(N - L)]
  #diffsum <- sapply(1:(N - L), function(j) { x0 + sum(fitted[1:j]) })
  plot(x=dat$time[1:(N - L)], y=xt[1:(N - L)], ylim=c(xmin, xmax),
       main=paste("SETAR(",p,",",d,",",round(c, digits=4),")"), 
       xlab="year", ylab="diff(unemp) [%]",type="p")
  lines(x=dat$time[1:(N - L)], y=fitted, lwd=2, col="royalblue2")
  lines(x=c(dat$time[1], dat$time[(N - L)]), y=c(c, c), col="brown3", lty="dashed", lwd=2)
  legend("topleft", legend=c("fitted SETAR","threshold"), col=c("royalblue2","brown3"), lty=c(1,2), lwd=2, cex=0.8)
  
  plot(x=dat$time[1:(N - L)], y=nonlinear[[i]]$residuals[1:(N - L)], type="l", lwd=1.5,
       main=paste("SETAR(",p,",",d,",",round(c, digits=4),") residuals"), xlab="year", ylab="diff(unemp) [%]", ylim=c(xmin, xmax))
  print( c(resSigmaSq = nonlinear[[i]]$resSigmaSq) )
}

#' ### 3.5 Conclusion ###
#' 
#' Since the differences in the unemployment rate have been used, we show the threshold value as well as the fitted values of the 
#' models in the same plot. The switching between the high and the low regimes can is clearly visible for all the selected models
#' (perhaps, except the second one, with its threshold value quite close to zero ). It is not yet clear whether another regime
#' should be present in the stochastic process. This will be assessed in the following chapter.
#' 
#' ## 4. 3-Regime SETARs and Diagnostic Tests of SETAR Models ##
#' 
#' The next step in the analysis using SETAR models is verifying whether 2 regimes suffice. If they do not, we will have 
#' to consider the possibility that a third regime needs to be added. In that case, we need to write methods for such model
#' 
#' ### 4.1 Useful Functions ###
#' 

# the indicator function for 3 regimes:
Indicator3 <- function(x, c) {
  tmp <- rep(F,3)
  tmp[findInterval(x, c, left.open = T) + 1] <- T
  tmp
}

Indicator3(-4, c(-1,1))
Indicator3(0, c(-1,1))
Indicator3(4, c(-1,1))

# SETAR3 basis vector
Yt <- function(x, t, p) c(1, x[(t - 1):(t - p)])

# SETAR3 skeleton
Xt <- function(x, t, p, d, c, z = x) {  
  # z is the threshold variable 
  I <- Indicator3(z[t - d], c)
  Y <- Yt(x, t, p)
  c(I[1] * Y, I[2] * Y, I[3] * Y)
}

suppressMessages(library(matlib))

# covariance matrix of the 3 regime SETAR
CovMat3 <- function(x, p, d, c) {
  n <- length(x)
  Yc <- matrix(0., ncol = (3 * p + 3), nrow = (3 * p + 3)) # this will become the covariance matrix
  k <- max(p, d)
  for (t in (k + 1):n) {
    XT <- Xt(x, t, p, d, c)
    Yc <- Yc + (XT %o% XT)
  }
  det <- det(Yc)
  if (det > -0.00001 && det < 0.00001) {
    return(NA)
  } else {
    return(inv(Yc))
  }
}

CovMat3(xt, p=2, d=1, c=c(-0.1, 0.2))

# skeleton of a model with regression coeffs: theta
SkeletonSETAR3 <- function(x, t, p, d, c, theta, z = x) c(theta %*% CovMat3(x, p, d, c)[t,]) 
SkeletonSETAR3(xt, t=3, p=2, d=1, c=c(-0.1, 0.2), rep(1, 3*3))

#' 
#' To find out, whether the third regime should be added, we need to test for the independence of residuals:
#'
#' ### 4.2 The Brock-Dechert-Scheinkman (BDS) Test ###
#' 



#' 
#' ## 5. Predictions via SETAR Models and Their Evaluation ##
#' 
#' ## 6. Tests for Non-Linearity of STAR Models ###
#' 
#' 
