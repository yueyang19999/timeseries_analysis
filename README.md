# Value-at-Risk and Expected Shortfall Analysis Using GARCH Models

This project applies **Value-at-Risk (VaR)** and **Expected Shortfall (ES)** analysis to equity positions in **Qualcomm (QCOM)** and **Apple (AAPL)**.  
The study explores multiple approaches, including RiskMetrics, GARCH-type models, and quantile regression, to estimate tail risk under different statistical assumptions.

---

## Objectives
- Compute **1-day and 10-day 99% VaR and ES** for $1M positions.
- Compare the performance of different volatility models:
  - RiskMetrics (IGARCH with fixed parameters)
  - GARCH(1,1) with Gaussian innovations
  - GARCH(1,1) with Student-t innovations
  - IGARCH(1,1) for Apple
  - GJRGARCH(1,1) with leverage effects
  - Quantile regression based on S&P 500 log returns
- Assess which approach produces the most appropriate risk measure.

---

## Tools & Libraries
- **R**  
- [`quantmod`](https://cran.r-project.org/package=quantmod) – financial data retrieval & returns  
- [`rugarch`](https://cran.r-project.org/package=rugarch) – GARCH family models  
- [`fGarch`](https://cran.r-project.org/package=fGarch) – alternative GARCH fitting  
- [`quantreg`](https://cran.r-project.org/package=quantreg) – quantile regression  

Data source: **Yahoo Finance** (QCOM, AAPL, ^GSPC).

---

## Workflow
1. **Data Collection**  
   - Download adjusted close prices for QCOM, AAPL, and S&P 500 (2016–2025).  
   - Compute daily log returns.  

2. **RiskMetrics Method (IGARCH)**  
   - Estimate conditional volatility with fixed `ω = 0`.  
   - Compute 1-day and 10-day 99% VaR & ES.  

3. **GARCH Models**
   - **Gaussian innovations** (no mean intercept).  
   - **Student-t innovations** (capture fat tails).  
   - Compute VaR/ES at 1-day and 10-day horizons.  

4. **Apple Stock Extensions**
   - IGARCH(1,1) with normal innovations.  
   - GJRGARCH(1,1) with leverage effects.  
   - Portfolio VaR combining AAPL & QCOM.  

5. **Quantile Regression Approach**
   - Estimate 99% quantile of AAPL losses conditional on S&P 500 returns.  

6. **Model Comparison**
   - Evaluate adequacy via diagnostics (residuals, QQ-plots).  
   - Compare VaR magnitudes across methods.  
   - Recommend GJRGARCH(1,1) as the most robust due to leverage effect modeling.

---

## Key Results (Illustrative)
- **RiskMetrics (IGARCH)**  
  - Produces relatively low VaR estimates because it assumes normally distributed returns and does not capture fat tails.  
  - Likely underestimates tail risk for both QCOM and AAPL.  

- **GARCH(1,1) with Gaussian innovations**  
  - Improves volatility modeling by capturing clustering.  
  - Still too light-tailed — extreme events are underestimated.  

- **GARCH(1,1) with Student-t innovations**  
  - Heavier tails → higher VaR values than Gaussian.  
  - Provides more conservative estimates, better matching real market shocks.  

- **IGARCH(1,1) for Apple**  
  - Captures persistence of volatility but can exaggerate long memory.  
  - Works as a baseline but less flexible compared to GJR-GARCH.  

- **GJRGARCH(1,1)**  
  - Captures **leverage effects** (asymmetric volatility: bad news → higher volatility than good news).  
  - Produces the most realistic VaR estimates, especially during drawdowns.  

- **Quantile Regression (S&P 500 as factor)**  
  - Provides a factor-based conditional VaR estimate.  
  - Useful for stress-testing against market moves, but less conservative than GJR-GARCH.  

** Overall:**  
The **GJRGARCH(1,1)** model with leverage effects gave the most robust and conservative 99% VaR estimates. RiskMetrics consistently underestimated extreme risks, while Student-t GARCH provided a strong improvement by incorporating fat tails.  


---

## How to Run
1. Clone this repo.
2. Open the R script in **RStudio** or run from R console.  
3. Install required packages:
   ```R
   install.packages(c("quantmod", "rugarch", "fGarch", "quantreg"))
