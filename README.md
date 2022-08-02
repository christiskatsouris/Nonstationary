# Nonstationary

# I. Structural Break Testing 

The notion of "nonstationarity" in time series analysis and time series econometrics has various interpreations which depends on the exact modelling environment as well as the assumptions we impose on model parameters and moment conditions of innovation sequences. Thus, we provide some examples related to nonstationarity in terms of the presence of parameter instability for commonly used models in empirical studies. For both examples, we assume that the model parameters are such that the underline stochastic processes are stationarity and ergodic. The following examples are implemented in R.    

## Example 1: Stationary first-order autoregressive model  

$$y_t = \rho y_{t-1} + e_t, \ \ \ 0< \rho < 1, \ \ \ t = 1,...,n, \ \ \ e_t \sim N (0,1).$$

Assume that there is a single break-point at an unknown location within the full sample. Then, the R code below can be employed to determine the location of the break-point (if one exists) using standard F-tests. 

```R

install.packages("strucchange")
library("strucchange")

# insert dataset
mydata <- read.table("series.txt")

# obtain subseries
ciss     <- mydata$V1
money    <- mydata$V2
finance  <- mydata$V3
bond     <- mydata$V4
equity   <- mydata$V5
exchange <- mydata$V6

# take the first differences
d1 <- diff(ciss)
d2 <- diff(money)
d3 <- diff(finance)
d4 <- diff(bond)
d5 <- diff(equity)
d6 <- diff(exchange)

nfci <- d5

# F-TESTS AND DETECTION OF BREAK for model trend and intercept
nfci.new  <- nfci[1:(n-1)]
nfci.diff <- diff(nfci)

fs.nfci <- Fstats(nfci.new ~ nfci.diff)
plot(fs.nfci, main="F-statistics for CFSI")
breakpoints(fs.nfci)
lines(breakpoints(fs.nfci))

bp.nfci <- breakpoints(nfci.new ~ nfci.diff)
summary(bp.nfci)
plot(bp.nfci)

# Exctracting the breakpoint
bp <- fs.nfci$breakpoint

```

## Example 2: GARCH model (conditional heteroscedasticity)  

Consider the following Garch(1,1) model

$$\sigma^2_t = \alpha_0 + \alpha_1 \epsilon^2_{t-1} + \beta_1 \sigma^2_{t-1}, \ \ \ t = 1,...,n, \ \ \ \epsilon_t = \sigma_t Z_t, \ \ Z_t \sim N (0,1).$$

### Remarks

The main feature of multiplicative noise models is that it enable us to capture the dependence structure of financial log-return series. In particular, the Garch(1,1) process is defined as above, where the conditional heteroscedasticity (changing variance over time) is captured by the second equation such that Zt is a sequence of identically distributed (i.i.d) random variables with zero mean and unit variance.  Furthermore, the first component represents the local conditional standard deviation of the process and the two sequences are assumed to be independent. Note also, that we consider only stationary Garch processes which implies that a1 + b1 < 1. 

The R code below fits a Garch(1,1) model and obtains the estimated innovation sequences which can be used to test for structural breaks using standard retrospective tests. 

```R

library("fGarch")

# GARCH volatility model 
specs1 <- garchSpec(model=list(alpha=alpha1, omega=omega1, beta=beta1)) 
sigma1 <- garchSim(spec = specs1, n = n)

est.garch <- garchFit(formula = ~ garch(1,1), data = sigma1[1:n], include.mean = FALSE, trace = F)
est.par   <- est.garch@fit$par[1:3]
innov     <- dat[1:n]/est.garch@sigma.t
y <- (innov)^2
x <- c(0, y[1:(n-1)])

# MODEL ESTIMATION 
gspec.ru <- ugarchspec(mean.model=list(armaOrder=c(1,0)), distribution="norm" )

garch.estim1 <- ugarchfit(gspec.ru, nfci)
garch.estim2 <- ugarchfit(gspec.ru, risk)
garch.estim3 <- ugarchfit(gspec.ru, credit)
garch.estim4 <- ugarchfit(gspec.ru, leverage)

garch.estim1
persistence(garch.estim1)

values1 <- residuals(garch.estim1, standardize=FALSE)

```

## Example 3: CUSUM test statistic  

In this Example we present an application of a sequential monitoring framework which can be employed for testing for structural breaks in time series using the OLS-CUSUM and OLS-MOSUM statistics. 

```R

### CUSUM TEST ###

# Step 1: Use the sub-sample of the times series with upper bound the estimated break-point

x <- as.vector(series[1:(bp-1)])
y <- as.vector(series[2:bp])
n <- length(nfci)

nfci.cusum  <- mefp(y ~ x, type="OLS-CUSUM", border=newborder3)

# Step 2: Use the time series of the full sample 

x <-as.vector(nfci[1:(n-1)])
y <-as.vector(nfci[2:n])

nfci.cusum <- monitor(nfci.cusum)
stat <- round(nfci.cusum $statistic, digit=4)
m    <- nfci.cusum $breakpoint  

### MOSUM TEST ###

# Step 1: Use the sub-sample of the times series with upper bound the estimated break-point

x <- as.vector(series[1:(bp-1)])
y <- as.vector(series[2:bp])
n <- length(nfci)

nfci.mosum <- mefp(y ~ x, type="OLS-MOSUM", border=newborder3, h=0.5, alpha=0.05, period=2)

# Step 2: Use the time series of the full sample 

x <- as.vector(nfci[1:(n-1)])
y <- as.vector(nfci[2:n])

nfci.mosum <- monitor(nfci.mosum)
stat <- round(nfci.mosum$statistic, digit=4)
m    <- nfci.mosum$breakpoint  

mos.bound <- zoo( c(rep(NA,(bp-1)), newborder3(bp:n)), index(nfci) )

# Ploting 
plot( zoo(c(nfci.mosum$efpprocess, nfci.mosum$process), index(nfci)), 
      ylim = c(-5, 5), 
      xlab = "Time", ylab = "empirical fluctuation process",
      main = "OLS-MOSUM Process with Boundary 3")
      
abline(0, 0)
lines(mos.bound, col = 2)
lines(-mos.bound, col = 2)

```
## Assignment 1

Using one of the following Statistical Software: [R](https://www.r-project.org/), [Matlab](https://uk.mathworks.com/help/matlab/getting-started-with-matlab.html) or [Stata](https://www.stata.com/bookstore/getting-started-windows/) prepare a short empirical study based on a suitable time series dataset. The main focus is the use of both the first-order autoregressive AR(1) model as well as the p-th order autoregressive conditionally heteroscedastic model ARCH(p) and GARCH(p,q) for time-series data. The analysis should include: (i) model fitting and residual analysis, such as testing for the presence of autocorrelation in the residual series and testing for the presence of Arch/Garch effects, and (ii) testing for structural breaks using both retrospective and sequential break-point tests. 

> Many of the econometric tools and state-of-the-art methodologies we introduce in this course are commonly used by government bodies, regulatory  authorities and   global banks for monitoring financial stability and economic conditions. Specifically, Early Warning Models (EWMs) can be classified into two broad categories, i.e., the models used to assess the performance of the banking sector as well as the models which are used to assess the overall economic and financial conditions for signalling banking, currency or economic crises. 

## References

(a) Restrospective Testing for Structural Break 
- Chow, G. C. (1960). Tests of equality between sets of coefficients in two linear regressions. Econometrica: Journal of the Econometric Society, 591-605.
- Toyoda, T. (1974). Use of the Chow test under heteroscedasticity. Econometrica: Journal of the Econometric Society, 601-608.

(b) Testing for Multiple Structural Breaks
- Andreou, E., & Ghysels, E. (2002). Detecting multiple breaks in financial market volatility dynamics. Journal of Applied Econometrics, 17(5), 579-600
- Bai, J., & Perron, P. (1998). Estimating and testing linear models with multiple structural changes. Econometrica, 47-78.

(b) Sequential Monitoring for Structural Breaks
- Chu, C. S. J., Stinchcombe, M., & White, H. (1996). Monitoring structural change. Econometrica: Journal of the Econometric Society, 1045-1065.
- Horváth, L., Hušková, M., Kokoszka, P., & Steinebach, J. (2004). Monitoring changes in linear models. Journal of statistical Planning and Inference, 126(1), 225-251.
- Zeileis, A., Leisch, F., Kleiber, C., & Hornik, K. (2005). Monitoring structural change in dynamic econometric models. Journal of Applied Econometrics, 20(1), 99-121.
- Zeileis, A., Leisch, F., Hornik, K., & Kleiber, C. (2002). strucchange: An R package for testing for structural change in linear regression models. Journal of Statistical Software, 7, 1-38.
- Katsouris, C. (2017). Sequential Break-Point Detection in Stationary Time Series: An Application to Monitoring Economic Indicators. [arXiv preprint:2112.06889](https://arxiv.org/abs/2112.06889).

(c) Garch and volatility modelling
- Bollerslev, T. (1986). Generalized autoregressive conditional heteroskedasticity. Journal of Econometrics, 31(3), 307-327.
- Babikir, A., Gupta, R., Mwabutwa, C., & Owusu-Sekyere, E. (2012). Structural breaks and GARCH models of stock return volatility: The case of South Africa. Economic Modelling, 29(6), 2435-2443.
- Stefan Richter, Weining Wang and Wei Biao Wu (2021) Testing for parameter change epochs in GARCH time series. The Econometrics Journal.
- Katsouris, C. (2021). Forecast Evaluation in Large Cross-Sections of Realized Volatility. [arXiv preprint:2112.04887](https://arxiv.org/abs/2112.04887).

(d) Testing for autocorrelation in time series models
- Henshaw Jr, R. C. (1966). Testing single-equation least squares regression models for autocorrelated disturbances. Econometrica: Journal of the Econometric Society, 646-660.
- Wallis, K. F. (1972). Testing for fourth order autocorrelation in quarterly regression equations. Econometrica: Journal of the Econometric Society, 617-636.

# Reading List

- Gourieroux, C., & Jasiak, J. (2018). Financial econometrics. In Financial Econometrics. Princeton University Press.
- Tsay, R. S. (2005). Analysis of financial time series. John wiley & sons. (Chapter 3: Conditional Heteroscedastic Models)
- Bauwens, L., Lubrano, M., & Richard, J. F. (2000). Bayesian inference in dynamic econometric models. OUP Oxford. (Chapter 7: Heteroscedasticity and Arch)

# II. Nonlinear Time-Series Models

Modelling nonlinearity (e.g., "changing-regime" dynamics or threshold-type non-linearity) in time series econometrics is commonly captured using threshold autoregressive specifications and threshold regression models. A typology of these models is provided by the nature of the switching function F(.) as well as by the nature of the switching variable which can be either time index (structural change) or a continuous variable (threshold variable). Furthermore, threshold regressions can be also employed for jointly modelling nonstationarity and nonlinearity, although we consider this case as a more advanced application of the current framework.    

## Application 1: Threshold Autoregression Model

$$y_t = 
\begin{cases}
\alpha_1 + \beta_1^{\top} y_{t-1} + u_{t} , \ \ q_t \leq \gamma \\
\alpha_2 + \beta_2^{\top} y_{t-1} + u_{t} , \ \ q_t > \gamma
\end{cases}$$

```R

install.packages("chngpt")
library("chngpt")

# Fitting Threshold Regression Models

```

## Assignment 2

Using one of the following Statistical Software: [R](https://www.r-project.org/), [Matlab](https://uk.mathworks.com/help/matlab/getting-started-with-matlab.html) or [Stata](https://www.stata.com/bookstore/getting-started-windows/) prepare a short empirical study based on a suitable time series dataset. In your dataset include the [Economic Policy Uncertainy](https://www.policyuncertainty.com/) (EPU) index as one of the variables. Using a dependent variable of your choice, a set of control variables and EPU as the threshold variable, investigate an economic application by fitting an appropriate threshold model. Describe in details the research question and the modelling approach. Then, using your dataset evaluate the performance of the model and explain the econometric intuition of your findings. Hint: The research question of your study should be tackling an economic theory problem that supports the use of a threshold regression. Some examples include: poverty traps, technological traps, financial development traps, liquidity traps or energy traps.  

### Remarks

Although for the purpose of this course, we focus on the implementation of the threshold-model based on time series data, threshold models are widely applicable in various settings. According to Lee, Seo and Shin (2011), the threshold effect can be specified as a change-point due to an unknown threshold in a particular covariate. In economics and sociology, racial segregation can be modeled as a threshold effect. For example, one can investigate the existence of race-based tipping in neighborhoods. Moreover, examples from economic empirical studies include Durlauf and Johnson (1995) who argue that cross-country growth models with multiple equilibria can exhibit threshold effects. In empirical finance, Peseran and Pick (2007) argue that the effect of financial contagion can be described as a discontinuous threshold effect, hence testing for threshold effects implies testing for the presence of financial contagion. In biostatistics, dose-response models are typically specified with some unknown threshold parameters. Lastly, in Epidemiology, logistic regressions with unknown change points are used to model the relationship between the continuous exposure variable and disease risk.   

## References

- Che, X., & Jiang, M. (2021). Economic Policy Uncertainty, Financial Expenditure and Energy Poverty: Evidence Based on a Panel Threshold Model. Sustainability, 13(21), 11594.
- Dagenais, M. G. (1969). A threshold regression model. Econometrica: Journal of Econometric Society, 193-203.
- Davies, R. B. (1987). Hypothesis testing when a nuisance parameter is present only under the alternative. Biometrika, 74(1), 33-43.
- Durlauf, S. N., & Johnson, P. A. (1995). Multiple regimes and cross‐country growth behaviour. Journal of Applied Econometrics, 10(4), 365-384.
- Fong, Y., Huang, Y., Gilbert, P. B., & Permar, S. R. (2017). chngpt: Threshold regression model estimation and inference. BMC bioinformatics, 18(1), 1-7.
- Galvao Jr, A. F., Montes-Rojas, G., & Olmo, J. (2013). A panel data test for poverty traps. Applied Economics, 45(14), 1943-1952.
- Gonzalo, J., & Pitarakis, J. Y. (2002). Estimation and model selection based inference in single and multiple threshold models. Journal of Econometrics, 110(2), 319-352.
- K. S. Chan and H. Tong (1990). On Likelihood Ratio Tests for Threshold Autoregression. Journal of the Royal Statistical Society, Series B, Methodological, 52, 469-476.  
- Lee, S., Seo, M. H., & Shin, Y. (2011). Testing for threshold effects in regression models. Journal of the American Statistical Association, 106(493), 220-231.
- Peng, F., Cebula, R. J., Foley, M., & Zhan, K. (2016). Estimation of the liquidity trap using a panel threshold model. Applied Economics Letters, 23(16), 1134-1137.
- Pesaran, M. H., & Pick, A. (2007). Econometric issues in the analysis of contagion. Journal of Economic Dynamics and Control, 31(4), 1245-1277.
- Yu, P. (2012). Likelihood estimation and inference in threshold regression. Journal of Econometrics, 167(1), 274-294.


## Application 2: Threshold Cointegration Model (Advanced Topics)

A threshold cointegration model allows to investigate nonlinearities which exhibit time series nonstationarity and cointegration dynamics. A suitable modelling framework to capture such dynamics is proposed by Hansen and Seo (2002). Furthermore, the threshold cointegration model allows to capture short-run dynamics and therefore is commonly used in international finance empirical studies to investigate related hypotheses, such as the twin deficit hypothesis. 

```R

# Fitting Threshold Cointegration Models


```

### Remarks

According to Ahmad, Aworinde and Martin (2015): 

> Persistent fiscal and current account deficits are a major policy concern, irrespective of whether the country is affected is developed or developing.  This is because large fiscal deficits may lead to crowding-out of private investment if they cause interest rates to rise. Similarly, a large current account deficit could lead to a decline in competitiveness, a transfer of wealth to foreign nationals and a depletion of foreign exchange reserves, possibly triggering a currency crisis. From the traditional open-economy macroeconomic perspective, there are three main reasons to expect a positive relationship between the fiscal deficit and the current account deficit, the “Twin Deficit” hypothesis. First, an increase in the fiscal deficit may induce
an increase in the interest rate that causes capital inflows and an appreciation of the exchange rate, with unfavourable effects on the current account. Second, an increase in the fiscal deficit may lead to an increase in the demand for imports, causing a worsening of the current account. And third, a worsening of the current account deficit will reduce tax revenue and thus increase the fiscal deficit. In contrast to the traditional Keynesian view, the Ricardian equivalence hypothesis of
Barro (1974, 1989) argues that the fiscal deficits and the current account deficits are unrelated. Kim and Roubini (2008) argue for a negative relationship, a “Twin Divergence” hypothesis.



## References

- Ahmad, A. H., Aworinde, O. B., & Martin, C. (2015). Threshold cointegration and the short-run dynamics of twin deficit hypothesis in African countries. The Journal of Economic Asymmetries, 12(2), 80-91.
- Hansen, B. E., & Seo, B. (2002). Testing for two-regime threshold cointegration in vector error-correction models. Journal of econometrics, 110(2), 293-318.

# Reading List
- Tsay, R. S. (2005). Analysis of financial time series. John wiley & sons. (Chapter 4: Nonlinear Models and Their Applications)
- Bauwens, L., Lubrano, M., & Richard, J. F. (2000). Bayesian inference in dynamic econometric models. OUP Oxford. (Chapter 8: Non-linear Time Series Models)

# Remarks

1. In this teaching page, we present some examples for financial time series modelling using linear and nonlinear time series regression models. We consider the notion of "nonstationarity" as the presence of time-variation in model parameters, which requires to employ structural break tests to identify the presence of breaksm while the notion of "nonlinearity" corresponds to changing regime dynamics. Therefore, unless otherwise specified with the term "nonstationary" we mean processes with non-constant moments, and so we consider processes whose sum of absolute autocovariances is finite.   

2. A different stream of literature, namely "nonstationary time series econometrics" focuses on the development of asymptotic theory and inference methods for nonstationary time series processes and regression models (e.g., using local-to-unity asymptotics). In particular, that case implies that the second moments of underline stochastic processes can be unbounded (e.g., presence of unit roots). The presentation of the particular applications is beyond the scope of the current teaching page.   

3. In practise the econometric analyses of macroeconomic time series versus financial time series often require different concepts and tools. Although there is an overlap of methodologies, usually the analysis of financial time series which includes stock returns and/or financial and economic indicators corresponds to stationary stochastic processes. On the other hand, the analysis of macroeconomic time series involves the use of unit root and cointegration theory due to the fact that the underline stochastic processes are considered to be nonstationary (not necessarily because of the presence of structural breaks, although this is an additional aspect of consideration). 

# Disclaimer

The author declares no conflicts of interest.

The proposed Course Syllabus is currently under development and has not been officially undergone quality checks. All rights reserved.

Any errors or omissions are the responsibility of the author.

# How to Cite a Website

See: https://www.mendeley.com/guides/web-citation-guide/
