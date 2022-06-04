# Nonstationary

# Structural Break Testing 

## Example 1: Stationary first-order autoregressive model  

$$y_t = \rho y_{t-1} + e_t, t = 1,...,n, e_t \sim N (0,1).$$

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
```

## Example 2: GARCH model (conditional heteroscedasticity)  

```R

# GARCH volatility model 
specs1 <- garchSpec(model=list(alpha=alpha1, omega=omega1, beta=beta1)) 
sigma1 <- garchSim(spec = specs1, n = n)

est.garch <- garchFit(formula = ~ garch(1,1), data = sigma1[1:n], include.mean = FALSE, trace = F)
est.par   <- est.garch@fit$par[1:3]
innov     <- dat[1:n]/est.garch@sigma.t
y <- (innov)^2
x <- c(0, y[1:(n-1)])


```

# References

- Bai, J., & Perron, P. (1998). Estimating and testing linear models with multiple structural changes. Econometrica, 47-78.
- Bollerslev, T. (1986). Generalized autoregressive conditional heteroskedasticity. Journal of econometrics, 31(3), 307-327.
- Zeileis, A., Leisch, F., Kleiber, C., & Hornik, K. (2005). Monitoring structural change in dynamic econometric models. Journal of Applied Econometrics, 20(1), 99-121.




