#####################################################
## Project title:                                  ##
## Value-at-Risk and Expected Shortfall Analysis   ##
## Using GARCH Models                              ##
#####################################################


# Setting
library(quantmod)
library(rugarch)

# Define the date range for the analysis
start_date <- "2016-04-01"
end_date   <- "2025-03-31"

# Define the stock symbols
symbols <- c("QCOM", "AAPL")
getSymbols(symbols, src = "yahoo", from = start_date, to = end_date)

# Extract the Adjusted Close prices for each stock using the Ad() function
QCOM_adj <- Ad(QCOM)
AAPL_adj <- Ad(AAPL)

# Calculate daily log returns for each stock
QCOM_ret <- log(1 + dailyReturn(QCOM_adj, type = "arithmetic"))
AAPL_ret <- log(1 + dailyReturn(AAPL_adj, type = "arithmetic"))

# Merge the returns into a single xts object for convenience
returns <- merge(QCOM_ret, AAPL_ret)
colnames(returns) <- c("QCOM", "AAPL")

# Display the first few rows of the returns data
head(returns)

# Set the portfolio position
position <- 1000000  # $1 million

# Calculating the 1-day 99% VaR and expected shortfall for
# for my position ($1 million) for the next trading day 
# using RiskMetrics method.
start_date <- "2016-04-01"
end_date   <- "2025-03-31"
getSymbols("QCOM", src = "yahoo", from = start_date, to = end_date)
QCOM_adj <- Ad(QCOM)
qcom_logret <- log(1 + dailyReturn(QCOM_adj, type = "arithmetic"))

plot(qcom_logret, type = 'l', main = "QCOM Daily Log Returns")
mean(qcom_logret)

rm_spec_norm <- ugarchspec(
  variance.model = list(model = "iGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm",
  fixed.pars = list(omega = 0)
)

qcom_rm <- ugarchfit(rm_spec_norm, data = qcom_logret)
qcom_rm
coef(qcom_rm)

lambda <- coef(qcom_rm)["beta1"]
alpha1 <- coef(qcom_rm)["alpha1"]

tail(qcom_logret)
tail(qcom_rm@fit$sigma)

fc <- ugarchforecast(qcom_rm, n.ahead = 1)
sigma1 <- fc@forecast$sigmaFor[1]
sigma1

qnorm(0.99) * sigma1             # 1% VaR (per $1 of investment)
dnorm(qnorm(0.99)) / 0.01        # Normal distribution constant for 1% ES
dnorm(qnorm(0.99)) / 0.01 * sigma1  # 1% ES (per $1)

position <- 1000000
VaR_1day <- qnorm(0.99) * sigma1 * position
ES_1day <- dnorm(qnorm(0.99)) / 0.01 * sigma1 * position

VaR_10day <- VaR_1day * sqrt(10)
ES_10day <- ES_1day * sqrt(10)

cat("Lambda (beta1):", lambda, "\n")
cat("1-day 99% VaR: $", round(VaR_1day, 2), "\n")
cat("1-day 99% ES:  $", round(ES_1day, 2), "\n")
cat("10-day 99% VaR: $", round(VaR_10day, 2), "\n")
cat("10-day 99% ES:  $", round(ES_10day, 2), "\n")

# we will now consider a GARCH(1,1) model with Gaussian 
# innovations, without an intercept term in the mean
# equation.
# 
library(quantmod)
library(fGarch)

start_date <- "2016-04-01"
end_date <- "2025-03-31"
getSymbols("QCOM", src = "yahoo", from = start_date, to = end_date)
QCOM_adj <- Ad(QCOM)

qcom_ret <- log(1 + dailyReturn(QCOM_adj, type = "arithmetic"))

nqcom <- -qcom_ret

m1 <- garchFit(~garch(1,1), data = nqcom, trace = FALSE)

pred1 <- predict(m1, n.ahead = 10)

mu_1d <- pred1$meanForecast[1]
sigma_1d <- pred1$standardDeviation[1]

mu_10d <- sum(pred1$meanForecast[1:10])
sigma_10d <- sqrt(sum(pred1$standardDeviation[1:10]^2))

VaR_1d <- mu_1d + qnorm(0.99) * sigma_1d
VaR_10d <- mu_10d + qnorm(0.99) * sigma_10d

position <- 1000000
VaR_1d_usd <- VaR_1d * position
VaR_10d_usd <- VaR_10d * position

cat("===== GARCH(1,1) with ARMA(2,0) (Normal) =====\n")
cat("1-day 99% VaR:  $", round(VaR_1d_usd, 2), "\n")
cat("10-day 99% VaR: $", round(VaR_10d_usd, 2), "\n")

# replacing the underlying distribution for innovations
# with student t-distribution
# Load required packages
library(quantmod)
library(rugarch)

start_date <- "2016-04-01"  
end_date <- "2025-03-31"
getSymbols("QCOM", src = "yahoo", from = start_date, to = end_date)
QCOM_adj <- Ad(QCOM)
qcom_ret <- log(1 + dailyReturn(QCOM_adj, type = "arithmetic"))

nqcom <- -qcom_ret

spec_t <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model = "std"
)

fit_t <- ugarchfit(spec = spec_t, data = nqcom)
fit_t

forecast_t <- ugarchforecast(fit_t, n.ahead = 10)
sigma_t <- forecast_t@forecast$sigmaFor
mu_t <- rep(0, 10)  # No intercept → mean = 0
nu <- coef(fit_t)["shape"]
q_t_99 <- qt(0.99, df = nu)

VaR_1d_t <- q_t_99 / sqrt(nu / (nu - 2)) * sigma_t[1]
VaR_10d_t <- q_t_99 / sqrt(nu / (nu - 2)) * sqrt(sum(sigma_t[1:10]^2))

position <- 1000000
VaR_1d_t_usd <- VaR_1d_t * position
VaR_10d_t_usd <- VaR_10d_t * position

cat("===== GARCH(1,1), Student-t, No Intercept =====\n")
cat("Degrees of Freedom (ν):", round(nu, 2), "\n")
cat("1-day 99% VaR:  $", round(VaR_1d_t_usd, 2), "\n")
cat("10-day 99% VaR: $", round(VaR_10d_t_usd, 2), "\n")




# other garch model analysis for AAPL

library(quantmod)
library(rugarch)
library(quantreg)

# 1. Load data
start_date <- "2016-04-01"
end_date   <- "2025-03-31"
# Fetch Apple and Qualcomm together
getSymbols(c("AAPL","QCOM"), src = "yahoo",
           from = start_date, to = end_date)

# Fetch S&P 500 separately so we get a clean object
SP500_xts <- getSymbols("^GSPC", src = "yahoo",
                        from = start_date, to = end_date,
                        auto.assign = FALSE)


AAPL_adj <- Ad(AAPL)
QCOM_adj <- Ad(QCOM)
GSPC_adj <- Ad(SP500_xts)


AAPL_ret <- log(1 + dailyReturn(AAPL_adj, type = "arithmetic"))
QCOM_ret <- log(1 + dailyReturn(QCOM_adj, type = "arithmetic"))
GSPC_ret <- log(1 + dailyReturn(GSPC_adj, type = "arithmetic"))
head(QCOM_ret)
position <- 1e6  # $1 million

## fitting IGARCH(1,1), μ = 0, α₀ = 0, normal innovations

spec_igarch <- ugarchspec(
  variance.model    = list(model = "iGARCH", garchOrder = c(1,1)),
  mean.model        = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model= "norm",
  fixed.pars        = list(omega = 0)     # only ω here
)
fit_igarch <- ugarchfit(spec_igarch, data = AAPL_ret)




# Forecast volatility 10 days ahead
fc_igarch <- ugarchforecast(fit_igarch, n.ahead = 10)
σ1_a  <- fc_igarch@forecast$sigmaFor[1]
σ10_a <- sqrt(sum(fc_igarch@forecast$sigmaFor^2))

VaR1_a  <- qnorm(0.99) * σ1_a  * position
ES1_a   <- (dnorm(qnorm(0.99)) / 0.01) * σ1_a  * position
VaR10_a <- qnorm(0.99) * σ10_a * position
ES10_a   <- (dnorm(qnorm(0.99)) / 0.01) * σ10_a  * position

# 10‑day VaR of combined AAPL + QCOM using sample covariance
rets_2      <- merge(AAPL_ret, QCOM_ret, join = "inner")
weights     <- c(1e6, 1e6)
cov_mat     <- cov(rets_2)
var_port1   <- t(weights) %*% cov_mat %*% weights
var_port10  <- var_port1 * 10
VaR10_port  <- qnorm(0.99) * sqrt(var_port10)

cat(
    "1‑day 99% VaR (AAPL):", round(VaR1_a, 2), "\n",
    "1‑day 99% ES  (AAPL):", round(ES1_a,  2), "\n",
    "10‑day 99% VaR (AAPL):", round(VaR10_a, 2), "\n",
    "10‑day 99% ES  (AAPL):", round(ES10_a,  2), "\n",
    "10‑day 99% VaR (Portfolio):", round(VaR10_port, 2), "\n\n")

## GJRGARCH(1,1), Normal Innovations

library(rugarch)

# 1. Specify and fit the GJRGARCH(1,1) model with no mean intercept
spec_gjr <- ugarchspec(
  variance.model    = list(model = "gjrGARCH", garchOrder = c(1,1)),
  mean.model        = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model= "norm"
)
fit_gjr <- ugarchfit(spec_gjr, data = AAPL_ret)

# 2. Write down the fitted model:
#    σ_t² = ω + α₁ (ε_(t-1))² + γ₁ (ε_(t-1))²·I(ε_(t-1) < 0) + β₁ (σ_(t-1))²
#    Estimated parameters (from coef(fit_gjr)):
print(coef(fit_gjr))
#    ω    = 1.474418e-05
#    α₁   = 2.164411e-02
#    γ₁   = 8.572375e-01
#    β₁   = 1.639202e-01 
#    Mean equation has no intercept term.

# 3. Test for serial correlation in standardized residuals
std_resid <- residuals(fit_gjr, standardize = TRUE)
box <- Box.test(std_resid, type = "Ljung-Box")
print(box)
#    Box-Ljung test on std. residuals:
#    X-squared = 0.48176, df = 1, p-value = 0.4876
#    → p-value > 0.05, no evidence of autocorrelation in residuals

# 4. Diagnostic plots
#    (i) ACF of standardized residuals
acf(std_resid,
    main = "ACF of Standardized Residuals")

#    (ii) ACF of squared standardized residuals
acf(std_resid^2,
    main = "ACF of Squared Standardized Residuals")

#    (iii) QQ‑plot of standardized residuals against Normal
qqnorm(std_resid,
       main = "QQ-Plot of Standardized Residuals")
qqline(std_resid)

# Comments on diagnostics:
#  - ACF(std_resid): All autocorrelation bars lie well within the confidence bounds.
#.                  This confirms the Box–Ljung test (p ≈ 0.49): there is no significant 
#                   serial correlation remaining in the standardized residuals.
#  - ACF(std_resid^2): Almost all autocorrelations of the squared residuals sit within 
#                   the bands, indicating the GJRGARCH model has largely captured the 
#.                  ARCH effects.
#                   A couple of very small spikes (e.g. around lag 25) are barely outside 
#                   the bounds—suggesting only minimal leftover volatility clustering, 
#                   which may be inconsequential for VaR estimation.
#  - QQ‑plot: The bulk of the points fall close to the 45° line, showing a generally good 
#                   fit to the Normal distribution in the center.
#                   In the tails, however, the points deviate markedly: 
#                   the left tail is slightly heavier than normal, and the right tail shows
#                   more extreme outliers.
#                   This indicates some residual kurtosis and slight asymmetry remain, 
#                   despite using GJRGARCH with normal innovations.

# 5. Recompute 1‑day and 10‑day VaR and ES for completeness
fc_gjr <- ugarchforecast(fit_gjr, n.ahead = 10)
σ1_b   <- fc_gjr@forecast$sigmaFor[1]
σ10_b  <- sqrt(sum(fc_gjr@forecast$sigmaFor^2))

VaR1_b  <- qnorm(0.99) * σ1_b  * position
ES1_b   <- (dnorm(qnorm(0.99)) / 0.01) * σ1_b  * position
VaR10_b <- qnorm(0.99) * σ10_b * position
ES10_b  <- (dnorm(qnorm(0.99)) / 0.01) * σ10_b * position

cat(
    "1‑day 99% VaR (AAPL):", round(VaR1_b,  2), "\n",
    "1‑day 99% ES  (AAPL):", round(ES1_b,   2), "\n",
    "10‑day 99% VaR:",       round(VaR10_b, 2), "\n",
    "10‑day 99% ES:",        round(ES10_b,  2), "\n")


## Quantile regression on S&P 500 returns

# 1) Fetch S&P 500 log‐returns
SP500_xts <- getSymbols("^GSPC",
                        src         = "yahoo",
                        from        = start_date,
                        to          = end_date,
                        auto.assign = FALSE)
GSPC_ret <- log(1 + dailyReturn(Ad(SP500_xts), type = "arithmetic"))

# 2) Inner‑merge AAPL and S&P returns so they align perfectly
merged_rets <- na.omit( merge(AAPL_ret, GSPC_ret, join = "inner") )
#    column 1 = AAPL_ret; column 2 = GSPC_ret

# 3) Pull out pure numeric vectors and name them loss & sp
loss_vec <- -as.numeric( merged_rets[,1] )  # losses = –AAPL return
sp_vec   <-  as.numeric( merged_rets[,2] )  # S&P return

# 4) Build the regression data.frame
df_q <- data.frame(
  loss = loss_vec,
  sp   = sp_vec
)

# 5) Fit the 99% quantile regression
library(quantreg)
fit_q <- rq(loss ~ sp, data = df_q, tau = 0.99)

# 6) Display coefficients and compute 1‑day VaR
cat("Quantile regression coefficients:\n")
print(coef(fit_q))

last_sp <- tail(df_q$sp, 1)
VaR1_c  <- predict(fit_q, newdata = data.frame(sp = last_sp)) * position
cat("1‑day 99% VaR via quantile regression:", round(VaR1_c, 2), "\n")


## Model comparison and choice

# Among the methods, the GJRGARCH(1,1) VaR is preferred, as it captures both volatility
# clustering and leverage effects.  It also produces a more conservative (higher) 99%
# VaR compared to the IGARCH and regression‑quantile approaches.



# Disclaimer
# this was an assignment from my timeseries class in 
# National University of Singapore, in QF4205. This 
# assignment has been done with the help of my professor
# as well as my teammate. 