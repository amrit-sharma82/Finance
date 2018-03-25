rm(list = ls())

vanillaBlack <- function(K, F, sigmaSqrtT, isCall) {
  logFonK <- log(F/K)
  oneOverSigmaSqrtT <- 1.0 / sigmaSqrtT
  halfSigmaSquaredT <- 0.5 * sigmaSqrtT * sigmaSqrtT
  dplus <- (logFonK + halfSigmaSquaredT) * oneOverSigmaSqrtT
  dminus <- (logFonK - halfSigmaSquaredT) * oneOverSigmaSqrtT

  #prices <- ifelse(isCall, F * pnorm(dplus) - K * pnorm(dminus), K * pnorm(-dminus) - F * pnorm(-dplus))

  putPrices <- K * pnorm(-dminus) - F * pnorm(-dplus)
  prices <- ifelse(isCall, F - K + putPrices, putPrices)
  return(prices)
}

volSwapIntegrand <- function(K, F, sigmaSqrtT) {
  isCall <- (K > F)
  strikeTo3On2 <- K * sqrt(K)
  blackOption <- vanillaBlack(K, F, sigmaSqrtT, isCall)
  value <- blackOption / strikeTo3On2
  return(value)
}

volSwap <- function(sigma, T, volOfVar) {
  var <- sigma * sigma
  volOfVarSqrtT <- volOfVar * sqrt(T)
  putIntegrandObj <- integrate(volSwapIntegrand, lower = 0.0, upper = var, F = var, sigmaSqrtT = volOfVarSqrtT)
  callIntegrandObj <- integrate(volSwapIntegrand, lower = var, upper = Inf, F = var, sigmaSqrtT = volOfVarSqrtT)
  putIntegrand <- putIntegrandObj$value
  callIntegrand <- callIntegrandObj$value
  integrand <- putIntegrand + callIntegrand
  fairVolSwap <- sqrt(var) - 0.25 * integrand
  vol <- sqrt(var)
  volSwapDF <- data.frame(fairVolSwap, vol, putIntegrand, callIntegrand, integrand)
  return(volSwapDF)
}

T <- 1.0
sigma <- 0.2
volOfVar <- 0.5
bumpSize <- 0.01

base <- volSwap(sigma, T, volOfVar)
bumped <- volSwap(sigma + bumpSize, T, volOfVar)

vega <- (bumped$fairVolSwap - base$fairVolSwap) / bumpSize
vsVega <- (bumped$vol - base$vol) / bumpSize
put_vega <- (bumped$putIntegrand - base$putIntegrand) / bumpSize
call_vega <- (bumped$callIntegrand - base$callIntegrand) / bumpSize

paste("vega = ", vega)
paste("put integrand vega = ", put_vega)
paste("call integrand vega = ", call_vega)
paste("vs fair strike vega = ", vsVega)
paste("assembled vega = ", vsVega - 0.25 * (put_vega + call_vega))

