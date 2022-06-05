# Nonstationary

# I. Structural Break Testing 

The notion of "nonstationarity" in time series analysis and time series econometrics has various interpreations which depends on the exact modelling environment as well as the assumptions we impose on model parameters and moment conditions of innovation sequences. Here, we provide some examples related to nonstationarity in terms of the presence of parameter instability in commonly used models in empirical studies. Notice that for both examples, we assume that the model parameters are such that the underline stochastic processes are stationarity and ergodic.    

## Example 1: Stationary first-order autoregressive model  

$$y_t = \rho y_{t-1} + e_t, \ \ \ 0< \rho < 1, \ \ \ t = 1,...,n, \ \ \ e_t \sim N (0,1).$$


```R
library("strucchange")

# insert dataset
mydata<-read.table("series.txt")

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

# Take series up to the breakpoint
x <- as.vector(series[1:(bp-1)])
y <- as.vector(series[2:bp])
  
series.cusum <- mefp(y ~ x, type=type, alpha=0.05, period=2, h=0.5, border=newboundary)
  
# take all the series and monitor
x <- as.vector(series[1:(n-1)])
y <- as.vector(series[2:n])
  
series.cusum <- monitor(series.cusum)
stat <- round(series.cusum$statistic, digit=4)
m    <- series.cusum$breakpoint  

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

values1<-residuals(garch.estim1, standardize=FALSE)

```

## Example 3: CUSUM test statistic  

```R

```

## References

- Andreou, E., & Ghysels, E. (2002). Detecting multiple breaks in financial market volatility dynamics. Journal of Applied Econometrics, 17(5), 579-600.
- Bai, J., & Perron, P. (1998). Estimating and testing linear models with multiple structural changes. Econometrica, 47-78.
- Bollerslev, T. (1986). Generalized autoregressive conditional heteroskedasticity. Journal of econometrics, 31(3), 307-327.
- Zeileis, A., Leisch, F., Kleiber, C., & Hornik, K. (2005). Monitoring structural change in dynamic econometric models. Journal of Applied Econometrics, 20(1), 99-121.
- Zeileis, A., Leisch, F., Hornik, K., & Kleiber, C. (2002). strucchange: An R package for testing for structural change in linear regression models. Journal of statistical software, 7, 1-38.




