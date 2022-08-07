# Nonstationary

A teaching page presenting various topics related to Time Series Econometrics Using R (Drafted July 2022).

### Course Overview:

The main philosophy with this course is to combine traditional statistical modelling methodologies with modern econometric specifications suitable for both cross-sectional and time series data. Emphasis with this course is to introduce some important economic and finance applications such as the monitoring of financial time series for structural breaks as well as the modelling of nonlinear dynamics in time series. Furthermore, we introduce state-of-the-art techniques and programming capabilities with R for each topic covered.

<img src="https://github.com/christiskatsouris/Nonstationary/blob/main/Data/time_series_plot.jpg" width="825"/>

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

<img src="https://github.com/christiskatsouris/Nonstationary/blob/main/Data/retrospective_tests.jpg" width="825"/>

## Example 2: GARCH model (conditional heteroscedasticity)  

Consider the following Garch(1,1) model

$$\sigma^2_t = \alpha_0 + \alpha_1 \epsilon^2_{t-1} + \beta_1 \sigma^2_{t-1}, \ \ \ t = 1,...,n, \ \ \ \epsilon_t = \sigma_t Z_t, \ \ Z_t \sim N (0,1).$$

### Remarks

The main feature of multiplicative noise models is that it enable us to capture the dependence structure of financial log-return series. In particular, the Garch(1,1) process is defined as above, where the conditional heteroscedasticity (changing variance over time) is captured by the second equation such that $Z_t$ is a sequence of identically distributed (i.i.d) random variables with zero mean and unit variance.  Furthermore, the first component represents the local conditional standard deviation of the process and the two sequences are assumed to be independent. Note also, that we consider only stationary Garch processes which implies that $a_1 + b_1 < 1$. 

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

> Many of the econometric tools and state-of-the-art methodologies we introduce in this course are commonly used by government bodies, regulatory  authorities and   global banks for monitoring financial stability and economic conditions. Specifically, Early Warning Models (EWMs) can be classified into two broad categories, i.e., the models used to assess the performance of the banking sector as well as the models which are used to assess the overall economic and financial conditions for signalling banking, currency or economic crises. 


```R

# R Output Example: Fitting Garch model on financial time series

> garch.estim1

*---------------------------------*
*          GARCH Model Fit        *
*---------------------------------*

Conditional Variance Dynamics 	
-----------------------------------
GARCH Model	: sGARCH(1,1)
Mean Model	: ARFIMA(1,0,0)
Distribution	: norm 

Optimal Parameters
------------------------------------
        Estimate  Std. Error  t value Pr(>|t|)
mu     -0.002042    0.007188 -0.28414 0.776304
ar1     0.387700    0.076242  5.08512 0.000000
omega   0.000318    0.000206  1.54072 0.123384
alpha1  0.242776    0.065444  3.70966 0.000208
beta1   0.756224    0.059930 12.61846 0.000000

Robust Standard Errors:
        Estimate  Std. Error  t value Pr(>|t|)
mu     -0.002042    0.006184 -0.33023 0.741226
ar1     0.387700    0.080344  4.82548 0.000001
omega   0.000318    0.000464  0.68491 0.493400
alpha1  0.242776    0.157520  1.54124 0.123259
beta1   0.756224    0.166300  4.54734 0.000005

LogLikelihood : 203.2973 

Information Criteria
------------------------------------
                    
Akaike       -1.9159
Bayes        -1.8354
Shibata      -1.9170
Hannan-Quinn -1.8834

Nyblom stability test
------------------------------------
Joint Statistic:  0.795
Individual Statistics:              
mu     0.05877
ar1    0.18327
omega  0.06785
alpha1 0.14477
beta1  0.05210

Asymptotic Critical Values (10% 5% 1%)
Joint Statistic:     	 1.28 1.47 1.88
Individual Statistic:	 0.35 0.47 0.75

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

## Remarks

Notice that the additional assumption of autocorrelated errors can be imposed when modeling financial time series, however a typical assumption in these models with autocorrelated errors is that the variance of the errors is constant over the sampling points (e.g., in discrete time models). 

<img src="https://github.com/christiskatsouris/Nonstationary/blob/main/Data/monitoring_series.jpg" width="725"/>

<img src="https://github.com/christiskatsouris/Nonstationary/blob/main/Data/monitoring_diff_series.jpg" width="725"/>

<img src="https://github.com/christiskatsouris/Nonstationary/blob/main/Data/monitoring_diff_series_p2.jpg" width="725"/>

<img src="https://github.com/christiskatsouris/Nonstationary/blob/main/Data/monitor1.png" width="545"/>

<img src="https://github.com/christiskatsouris/Nonstationary/blob/main/Data/monitoring_residuals.jpg" width="620"/>

Obtaining useful insights and information regarding the economic and financial behaviour of individuals (such as economic agents, investors, traders) requires the development of robust statistical procedures for monitoring time series observations. In particular, these state-of-the-art econometric methodologies can be employed as a signaling mechanism during periods of possible financial turbulence and thus preventing the occurrence of a financial crisis. 

## References

(a) Restrospective Testing for Structural Break 
- Chow, G. C. (1960). Tests of equality between sets of coefficients in two linear regressions. Econometrica: Journal of the Econometric Society, 591-605.
- Toyoda, T. (1974). Use of the Chow test under heteroscedasticity. Econometrica: Journal of the Econometric Society, 601-608.
- Zeileis, A. (2005). A unified approach to structural change tests based on ML scores, F statistics, and OLS residuals. Econometric Reviews, 24(4), 445-466.

(b) Testing for Multiple Structural Breaks
- Andreou, E., & Ghysels, E. (2002). Detecting multiple breaks in financial market volatility dynamics. Journal of Applied Econometrics, 17(5), 579-600
- Bai, J., & Perron, P. (1998). Estimating and testing linear models with multiple structural changes. Econometrica, 47-78.

(b) Sequential Monitoring for Structural Breaks
- Chu, C. S. J., Stinchcombe, M., & White, H. (1996). Monitoring structural change. Econometrica: Journal of the Econometric Society, 1045-1065.
- Chevallier, J. (2011). Detecting instability in the volatility of carbon prices. Energy Economics, 33(1), 99-110.
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

# II. Nonlinear Time-Series Regression Models

Modelling nonlinearity (e.g., "changing-regime" dynamics or threshold-type non-linearity) in time series econometrics is commonly captured using threshold autoregressive specifications and threshold regression models. A typology of these models is provided by the nature of the switching function $F(.)$ as well as by the nature of the switching variable which can be either time index (structural change) or a continuous variable (threshold variable). Furthermore, threshold regressions can be also employed for jointly modelling nonstationarity and nonlinearity, although we consider this case as a more advanced application of the current framework.    

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

- Baker, S. R., Bloom, N., & Davis, S. J. (2016). Measuring economic policy uncertainty. The Quarterly Journal of Economics, 131(4), 1593-1636.
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

Generally, the framework of the cointegration model allows to capture various dynamics and thus investigate international finance and international macroeconomic puzzles such as the Purchasing Power Parity (PPP) condition. In particular, the PPP theory states that in the long run, the exchange rate between two countries is determined by their relative price levels (see, Chapter 6 in Engel, Mark and West (2007) and Deardorff (1979) for related definitions). The PPP theory puzzle provides a suitable economic theory along with time series that exhibit cointegrated dynamics and therefore allows for the implementation of different econometric techniques such as testing for unit roots (ADF), testing for unit roots at the conditional quantile distributions, nonlinear unit root testing (KPPS), testing for the presence of structural breaks and tests for cointegration (see, [Nelson and Plosser (1982)](https://www.sciencedirect.com/science/article/pii/0304393282900125) and Corbae and Ouliaris (1988)). Specifically, exchange rates and considered to be time series which exhibit persistence characteristics and thus the aformentioned testing methodologies are sufficient statistics. In addition, many empirical studies in the literature have demonstrated the existence of strong evidence of nonlinear mean reversion of real exchange rates towards a stable long-run equilibrium and thus validating the long-run PPP hypothesis for certain economic regions.  

Furthermore, cointegration modeling can be used to test linear rational expectation specifications such as present-value relations. For instance, present-value models arise from agents'intertemporal allocation of funds for consumption or investment. During a recession current incomes are low relative to expected permanent incomes and the permanent income hypothesis predicts that savings will also be low. Therefore, to ensure that at the equilibrium market conditions the expected return on savings must increase to induce an increase in the funds allocated to investment. On the other hand, the particular economic theory might not always be applicable, especially in cases of nonlinear equilibrium relationships. For example, show that when the expected rate of return varies over time the present-value model does not generally imply the existence of a stationary relationship between the integrated forcing variable and the endogenous variables in levels (see, [Timmermann (1995)](https://www.jstor.org/stable/2284940#metadata_info_tab_contents)).    

A threshold cointegration model can be employed to investigate nonlinearities which exhibit time series nonstationarity and cointegration dynamics simultaneously. In particular, the threshold cointegration model allows to capture short-run dynamics and therefore is commonly used in international finance empirical studies to investigate hypotheses such as the "twin deficit" hypothesis (“Twin Divergence”). A suitable modelling framework to capture the dynamics under the "twin deficit" hypothesis is proposed by [Hansen and Seo (2002)](https://www.sciencedirect.com/science/article/pii/S0304407602000970) as briefly explained below. Similarly, the particular modeling approach can be also employed for investigating the PPP hypothesis under regime-specific dynamics (see, [Gouveia and Rodrigues (2004)](https://www.tandfonline.com/doi/abs/10.1080/0266476032000148984)).


```R

# Fitting Threshold Cointegration Models


```

### Remarks

According to Ahmad, Aworinde and Martin (2015): 

> Persistent fiscal and current account deficits are a major policy concern, irrespective of whether the country is affected is developed or developing.  This is because large fiscal deficits may lead to crowding-out of private investment if they cause interest rates to rise. Similarly, a large current account deficit could lead to a decline in competitiveness, a transfer of wealth to foreign nationals and a depletion of foreign exchange reserves, possibly triggering a currency crisis. From the traditional open-economy macroeconomic perspective, there are three main reasons to expect a positive relationship between the fiscal deficit and the current account deficit, the “Twin Deficit” hypothesis. First, an increase in the fiscal deficit may induce
an increase in the interest rate that causes capital inflows and an appreciation of the exchange rate, with unfavourable effects on the current account. Second, an increase in the fiscal deficit may lead to an increase in the demand for imports, causing a worsening of the current account. And third, a worsening of the current account deficit will reduce tax revenue and thus increase the fiscal deficit. In contrast to the traditional Keynesian view, the Ricardian equivalence hypothesis of
Barro (1974, 1989) argues that the fiscal deficits and the current account deficits are unrelated. Kim and Roubini (2008) argue for a negative relationship, a “Twin Divergence” hypothesis.


## References

- Ayinla, A. S. (2018). Why Segregating Cointegration Test?. American Journal of Applied Mathematics and Statistics, 6(4), 121-125.
- Ahmad, A. H., Aworinde, O. B., & Martin, C. (2015). Threshold cointegration and the short-run dynamics of twin deficit hypothesis in African countries. The Journal of Economic Asymmetries, 12(2), 80-91.
- Balke, N. S., & Fomby, T. B. (1997). Threshold cointegration. International economic review, 627-645.
- Corbae, D., & Ouliaris, S. (1988). Cointegration and tests of purchasing power parity. The Review of Economics and Statistics, 508-511.
- Dhrymes, P. J. (1969). Efficient estimation of distributed lags with autocorrelated errors. International Economic Review, 10(1), 47-67.
- Deardorff, A. V. (1979). One-way arbitrage and its implications for the foreign exchange markets. Journal of Political Economy, 87(2), 351-364.
- Engel, C., Mark, N. C., West, K. D. (2007). Exchange rate models are not as bad as you think. NBER macroeconomics annual, 22, 381-473.
- Gonzalo, J., & Pitarakis, J. Y. (2006). Threshold effects in cointegrating relationships. Oxford Bulletin of Economics and Statistics, 68, 813-833.
- Gouveia, P., & Rodrigues, P. (2004). Threshold cointegration and the PPP hypothesis. Journal of Applied Statistics, 31(1), 115-127.
- Hansen, B. E., & Seo, B. (2002). Testing for two-regime threshold cointegration in vector error-correction models. Journal of Econometrics, 110(2), 293-318.
- Kwiatkowski, D., Phillips, P. C., Schmidt, P., & Shin, Y. (1992). Testing the null hypothesis of stationarity against the alternative of a unit root: How sure are we that economic time series have a unit root?. Journal of econometrics, 54(1-3), 159-178.
- Li, J., & Lee, J. (2010). ADL tests for threshold cointegration. Journal of Time Series Analysis, 31(4), 241-254.
- Nelson, C. R., & Plosser, C. R. (1982). Trends and random walks in macroeconmic time series: some evidence and implications. Journal of monetary economics, 10(2), 139-162.
- Prakash-Canjels, G., & Taylor, A. M. (1997). Measuring market integration: a model of arbitrage with an econometric application to the gold standard, 1879-1913.
- Shleifer, A., & Vishny, R. W. (1997). The limits of arbitrage. The Journal of finance, 52(1), 35-55.
- Timmermann, A. (1995). Cointegration tests of present value models with a time‐varying discount factor. Journal of Applied Econometrics, 10(1), 17-31.

# Reading List
- Tsay, R. S. (2005). Analysis of financial time series. John wiley & sons. (Chapter 4: Nonlinear Models and Their Applications)
- Bauwens, L., Lubrano, M., & Richard, J. F. (2000). Bayesian inference in dynamic econometric models. OUP Oxford. (Chapter 8: Non-linear Time Series Models)
- Hamilton, J. D. (1994). Time Series Analysis. Princeton University Press. (Chapter 19: Cointegration)
- Dhrymes, Phoebus J. (1971). Distributed lags: Problems of Estimation and Formulation. 

# Comments

1. In this teaching page, we present some examples for financial time series modelling using linear and nonlinear time series regression models. We consider the notion of "nonstationarity" as the presence of time-variation in model parameters, which requires to employ structural break tests to identify the presence of breaksm while the notion of "nonlinearity" corresponds to changing regime dynamics. Therefore, unless otherwise specified with the term "nonstationary" we mean processes with non-constant moments, and so we consider processes whose sum of absolute autocovariances is finite.   

2. A different stream of literature, namely "nonstationary time series econometrics" focuses on the development of asymptotic theory and inference methods for nonstationary time series processes and regression models (e.g., using local-to-unity asymptotics). In particular, that case implies that the second moments of underline stochastic processes can be unbounded (e.g., presence of unit roots). The presentation of the particular applications is beyond the scope of the current teaching page.   

3. In practise the econometric analyses of macroeconomic time series versus financial time series often require different concepts and tools. Although there is an overlap of methodologies, usually the analysis of financial time series which includes stock returns and/or financial and economic indicators corresponds to stationary stochastic processes. On the other hand, the analysis of macroeconomic time series involves the use of unit root and cointegration theory due to the fact that the underline stochastic processes are considered to be nonstationary (not necessarily because of the presence of structural breaks, although this is an additional aspect of consideration). 


# Learning Outcomes

1. Understand the basic properties of time series regression models under different econometric assumptions. 
2. Understand the statistical procedures for implementing structural break tests under different modelling environments. 
3. Be able to apply econometric techniques in order to uncover structural breaks, unit roots and nonlinear dynamics in time series.  
4. Be able to relate empirical results from econometric methodologies with commonly testable hypotheses of international finance/macroeconomics.    
5. Be able to use Statistical/Econometric Programming Software such as [R](https://www.r-project.org/), [Matlab](https://uk.mathworks.com/help/matlab/getting-started-with-matlab.html) or [Stata](https://www.stata.com/bookstore/getting-started-windows/).

# Disclaimer

The author (Christis G. Katsouris) declares no conflicts of interest.

The proposed Course Syllabus is currently under development and has not been officially undergone quality checks. All rights reserved.

Any errors or omissions are the responsibility of the author.

Any views and opinions expressed herein are those of the author. The author accepts no liability for any loss or damage a person suffers because that person has directly or indirectly relied on any information found on this website.

# Acknowledgments

The author greatfully acknowledges financial support from Graduate Teaching Assistantships at the School of Economic, Social and Political Sciences of the University of Southampton as well as funding from Research Grants of interdisciplinary Centers of research excellence based at the University of Cyprus (UCY) as well as at the University College London (UCL). Furthermore, the author gratefully acknowledges financial support from the Vice-Chancellor's PhD Scholarship of the University of Southampton, for the duration of the academic years 2018 to 2021.

If you are interested to collaborate on any of the topics discussed in this teaching page, don't hesitate to contact me at christiskatsouris@gmail.com

# How to Cite a Website

See: https://www.mendeley.com/guides/web-citation-guide/
