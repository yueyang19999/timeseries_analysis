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
- 1-day and 10-day VaR/ES estimates for QCOM and AAPL.  
- GJRGARCH(1,1) produced the most conservative VaR, capturing volatility clustering + asymmetry.  
- RiskMetrics underestimates tail risk.  
- Quantile regression provides an alternative, but less conservative, measure.

---

## How to Run
1. Clone this repo.
2. Open the R script in **RStudio** or run from R console.  
3. Install required packages:
   ```R
   install.packages(c("quantmod", "rugarch", "fGarch", "quantreg"))
