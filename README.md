# Nonstationary

# I. Structural Break Testing 

The notion of "nonstationarity" in time series analysis and time series econometrics has various interpreations which depends on the exact modelling environment as well as the assumptions we impose on model parameters and moment conditions of innovation sequences. Thus, we provide some examples related to nonstationarity in terms of the presence of parameter instability for commonly used models in empirical studies. For both examples, we assume that the model parameters are such that the underline stochastic processes are stationarity and ergodic. The following examples are implemented in R.    

## Example 1: Stationary first-order autoregressive model  

$$y_t = \rho y_{t-1} + e_t, \ \ \ 0< \rho < 1, \ \ \ t = 1,...,n, \ \ \ e_t \sim N (0,1).$$


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

$$\sigma^2_t = \alpha_0 + \alpha_1 \epsilon^2_{t-1} + \beta_1 \sigma^2_{t-1}, \ \ \ t = 1,...,n, \ \ \ \epsilon_t = \sigma_t z_t, \ \ z_t \sim N (0,1).$$

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

## References

- Chow, G. C. (1960). Tests of equality between sets of coefficients in two linear regressions. Econometrica: Journal of the Econometric Society, 591-605.
- Toyoda, T. (1974). Use of the Chow test under heteroscedasticity. Econometrica: Journal of the Econometric Society, 601-608.
- Bollerslev, T. (1986). Generalized autoregressive conditional heteroskedasticity. Journal of Econometrics, 31(3), 307-327.
- Chu, C. S. J., Stinchcombe, M., & White, H. (1996). Monitoring structural change. Econometrica: Journal of the Econometric Society, 1045-1065.
- Bai, J., & Perron, P. (1998). Estimating and testing linear models with multiple structural changes. Econometrica, 47-78.
- Andreou, E., & Ghysels, E. (2002). Detecting multiple breaks in financial market volatility dynamics. Journal of Applied Econometrics, 17(5), 579-600.
- Zeileis, A., Leisch, F., Kleiber, C., & Hornik, K. (2005). Monitoring structural change in dynamic econometric models. Journal of Applied Econometrics, 20(1), 99-121.
- Zeileis, A., Leisch, F., Hornik, K., & Kleiber, C. (2002). strucchange: An R package for testing for structural change in linear regression models. Journal of Statistical Software, 7, 1-38.


## Further Reading

- Gourieroux, C., & Jasiak, J. (2018). Financial econometrics. In Financial Econometrics. Princeton University Press.
- Tsay, R. S. (2005). Analysis of financial time series. John wiley & sons. (Chapter 3: Conditional Heteroscedastic Models)
- Bauwens, L., Lubrano, M., & Richard, J. F. (2000). Bayesian inference in dynamic econometric models. OUP Oxford. (Chapter 7: Heteroscedasticity and Arch)

# II. Nonstationarity and Nonlinearity

Modelling jointly nonstationarity and nonlinearity (e.g., "changing-regime" dynamics or threshold-type non-linearity) in time series econometrics is commonly captured using threshold autoregressive specifications and threshold regression models. A typology of these models is provided by the nature of the switching function F(.) as well as by the nature of the switching variable which can be either time index (structural change) or a continuous variable (threshold variable).   

## Example 1: Threshold Autoregression Model

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

## References

- Dagenais, M. G. (1969). A threshold regression model. Econometrica: Journal of Econometric Society, 193-203.
- K. S. Chan and H. Tong (1990). On Likelihood Ratio Tests for Threshold Autoregression. Journal of the Royal Statistical Society, Series B, Methodological, 52, 469-476.  
- Gonzalo, J., & Pitarakis, J. Y. (2002). Estimation and model selection based inference in single and multiple threshold models. Journal of Econometrics, 110(2), 319-352.
- Yu, P. (2012). Likelihood estimation and inference in threshold regression. Journal of Econometrics, 167(1), 274-294.
- Fong, Y., Huang, Y., Gilbert, P. B., & Permar, S. R. (2017). chngpt: Threshold regression model estimation and inference. BMC bioinformatics, 18(1), 1-7.

## Further Reading
- Tsay, R. S. (2005). Analysis of financial time series. John wiley & sons. (Chapter 4: Nonlinear Models and Their Applications)
- Bauwens, L., Lubrano, M., & Richard, J. F. (2000). Bayesian inference in dynamic econometric models. OUP Oxford. (Chapter 8: Non-linear Time Series Models)

# Remarks

1. In this teaching page, we present some examples for financial time series modelling using linear and non-linear time series models. We consider the notion of "nonstationarity" as the presence of time-variation in model parameters, which requires to employ structural break tests to identify the presence of breaks. 

2. A different stream of literature, namely "nonstationary time series econometrics" focuses on the development of asymptotic theory and inference methods for nonstationry time series models (e.g., using local-to-unity asymptotics). The presentation of the particular applications is beyond the scope of the current teaching page.   
