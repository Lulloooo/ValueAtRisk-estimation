#############################################################################################################################
#                          VALUE AT RISK AND EXPECTED SHORTFALL
#***************************************************************************************************************************#
# Program = VaR-code.R
# Programmer = Luca Albertini
# Date first logger: 2024-02-04
#
#            Description: This programs estimates 1-day and multi-days value at risk and expected shortfall
#                         both for a single assets and for a multi-assets portfolio. It does so through
#                         different methods.
#
#            Input files: St.Louis Fed website; Yahoo finance
#            Temp files: none
#            Output files: none
#
#**************************************************************************************************************************#

#set the working directory
setwd("~/Personal/CaseStudies/ValueAtRisk-estimation")
#install the library
install.packages("quantmod")
library(quantmod)
#chose the variable that must be analyzed (the S&P500 index in this case)
sp500 <- getSymbols("SP500",src = "FRED", auto.assign = FALSE)
#drop the NA
sp500 <- na.omit(sp500)
#set the time-frame we are interested into (2014-2024)
sp500 <- sp500["2014-03-17/2024-03-18"]
#give the name idx (index) to the data within the dataset
names(sp500) <- "idx"
#display the first 3 lines
head(sp500, 3)
#display the last 3 lines
tail(sp500, 3)
#calculate the logarithm return
library(magrittr)
logidx <- sp500$idx %>% 
  log() %>% 
  diff()
#compute the daily returns without the first NA
logidx <- logidx[-1]
#display the first 3 rows with just 6 decimals (round command)
round(head(logidx, 3),6) #interpeting: the number says how much stcks fall/rise in % (values include dividens)
#find the index from the logidx
idx <-  exp(logidx)-1
round(head(idx, 3),6) #show them
#### Weekly LogIndex value  (can be done also wit monthly, quarterly etc)
logidxweek <- apply.weekly(logidx, sum) #sum tells how to add return
round(head(logidxweek, 3),6)
##weekly returns
idxweek <-  exp(logidxweek)-1
#display them
round(head(idxweek, 3),6)
#compute the mean
mu <- round(mean(logidx),8) #this is mu
#compute the standard deviation
sig <- round(sd(logidx),8) #this is sig
#set the investment
inv <- 1000
#set the alpha
alpha <- 0.05


################ VaR COMPUTING METHODS

##### NORMALITY METHOD -> assume that logindex are normally distributed
#compute the quantile (or the VaR) of our normal distribution 
VaRnorm <- qnorm(0.05, mu,sig) #this command works just for normal distributions
#compute the actual amount of money that can be lost
lossVaRnorm <- inv * (exp(VaRnorm)-1)
#calculate the epxected shortfall
esnorm <- mu-sig*dnorm(qnorm(0.05,0,1),0,1)/0.05
#calculate the actual loss
lossEsnorm <- inv * (exp(esnorm)-1)


#HISTORICAL METHOD
#function to order xts object (https://rdrr.io/rforge/xts/src/R/sort.xts.R)
sort.xts <-
  function(x, decreasing=FALSE, MARGIN=1, ...)
  {
    if(NCOL(x) > 1) {
      as.matrix(x)[order(x[,MARGIN],decreasing=decreasing,...),]
    } else as.matrix(x)[order(x,decreasing=decreasing,...),]
  }
#order the logreturns
logidxSort <- sort.xts(logidx)
#compute the VaR
VaRhist <- quantile(logidxSort, alpha)
round(VaRhist,6)
#compute the ES
EShist <- mean(logidx[logidx<VaRhist]) #consider just number that are less than the VaR
round(EShist,6)


## ESTIMATION METHOD
#normal distribution and real data
#set the RNGVersion (this is not a must, it is just to use rounding sample)
RNGversion("3.5.1")
#set the seed value (it is used to generate pseudo random number, i.e sequence of numbers that look random)
set.seed(123789) 
#why to set the seed? because in this way when generating random number the program will give me always the same random number
#as such the code become redproducible and the same results can be obtained
#if the seed is not set, every time the code is run it will displays different random number -> result will be different everytime I run the program 
#Note that in real world, the seed is not done and the number are always left random (i.e the seed is not set)
#get 100000 random numbers from the normal distribution with mean mu and sd sig
idxvecrandom <- rnorm(100000, mu, sig) #store that numbers in idxvec
#if you want to generate the sample from the ACTUAL return
idxvecreal <- sample(as.vector(logidx), 100000, replace = TRUE) #put replace = TRUE so R know to make a sample with replacements
#compute the VaR (as a quantile) with random data and assuming normality
VaRrandom <- quantile(idxvecrandom, alpha)
round(VaRrandom,6)
#compute the ES with random data assuming normality
ESrandom <- mean(idxvecrandom[idxvecrandom<VaRrandom]) #consider just number that are less than the VaR
round(ESrandom,6)
#compute the VaR (as a quantile) with real data
VaRreal <- quantile(idxvecreal, alpha)
round(VaRreal,6)
#compute the ES with real data
ESreal <- mean(idxvecreal[idxvecreal<VaRreal]) #consider just number that are less than the VaR
round(ESreal,6)

#check if the series is NORMALLY distributed -> kurtosis, skewness and jarque-bera test for normality
library(moments)
#transfom logidx into a vector
vidx <- as.vector(logidx)
#kurtosis
round(kurtosis(vidx), 2)
#skewness
round(skewness(vidx),2)
#normality test
jarque.test(vidx)

#T-STUDENT ESTIMATION
#using the student-t estimate the VaR and ES (0.95)
# estimate the student t parameters
library(MASS)
t_fit <- (fitdistr(vidx, "t"))
round(t_fit$estimate,6)
#m = mean, s = standard deviation (or rescaling parameter),  df = degree of freedom
#set the seed
set.seed(123789)
#call the library
library(metRology)
#set the standardized student-t parameters values (mean, sd and degree of freedom)
vidx2 <- rt.scaled(100000,mean=t_fit$estimate[1],sd=t_fit$estimate[2],df=t_fit$estimate[3])
#estimate VaR
VaR_t <- quantile(vidx2, alpha)
round(VaR_t, 6)
#estimate the ES
ES_t <- mean(vidx[vidx<VaR_t])
round (ES_t, 6)

####################### 10 DAYS TIME HORIZON
###### student-t method
#set alpha if not done
#alpha <- 0-05
#set the seed
set.seed(123789)
#create a vector with 100000 zeroes 
vidx3 <- rep(0,100000)
#build the loop to execute it 10 time (1:10) -> 
#each time the loop generates 100,000 outcomes from the estimated (rescaled) t-distribution
#and then add it to the vector vidx3 (for 10 times)
for (i in 1:10) {
  vidx3 <- vidx3+rt.scaled(100000, mean=t_fit$estimate[1], sd=t_fit$estimate[2], df=t_fit$estimate[3])
}
#at the end, the vidx3 is composed by 100,000 simulated 10-day log returns
#compute the VaR
VaR10 <- round(quantile(vidx,alpha),6)
VaR10
#compute the ES
ES10 <- round(mean(vidx[vidx<VaR10]),6)
ES10

####### IID method
#set the seed
#set.seed(123789)
#build the variable
vidx4 <- rep(0,100000)
#the loop is as before, but instad of using student-t the real data distirbution is used
for (i in 1:10) {
  vidx4 <- vidx4+sample(as.vector(logidx), 100000, replace=TRUE)
} 
#compute the VaR
VaR_IID <- quantile(vidx4,alpha)
round(VaR_IID, 6)
#compute the ES
ES_IID <- round(mean(vidx4[vidx4<VaR_IID]),6)
round(ES_IID, 6)

########## block simulation
vidx5 <- rep(0, 100000)
#set a vector that contains the 1day returns
rdat <- as.vector(logidx)
#set a vector starting from 1 and ending at number of rdat last night
posn <- seq(from=1, to=length(rdat)-9, by=1)
#sample the posn randomly and store these random position in a vector
rpos <- sample(posn,100000, replace = TRUE)
#build the loop
#(first step: dd the 1-day return at the starting position of the vector and immediately after that move the positioning forward by one
# second step:add the one day return from the next day to the RVEC vector. At the end of all ten passes of the loop, I have some ten day returns over ten consecutive days. )
for (i in 1:10){
  vidx5 <- vidx5 + rdat[rpos]  rpos <- rpos+1
}
#compute the VaR
VaRblock <- quantile(vidx5,alpha)
round(VaRblock, 6)
#compute the ES
ESblock <- mean(vidx5[vidx5<VaRblock])
round(ESblock, 6)

#### VOLATILITY CLUSTERING AND GARCH
library(rugarch)
#acf of the log return
acf(logidx, lag = 100)
#volatility cluserting (acf of the absolute log returns)
acf(abs(logidx), lag = 100)
#set the parameters for the GARCH model with rescaled-t assumption
garch_t <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), #variance model: GARCH(1,1)
                      mean.model = list(armaOrder = c(0,0), include.mean = TRUE), #mean model: cosntant term equals the mean
                      distribution.model = "std") #distirbution model: rescaled t distribution
#estimate the GARCH model
fitGarch_t <- ugarchfit(spec = garch_t, data = logidx)
fitGarch_t
#Mu = a sub 0 in the mean eqaution ; omega = alpha_0 in the variance eqaution ; alpha1 = alpha_1 in variance equation ; Beta1 = Beta_1 in the variance equation
#save some outputs (dataframe with logreturns, fitted value of daily sd, fitted vaule of Epsilons in the distribution equation)
save_t <- cbind(logidx,fitGarch_t@fit$sigma, fitGarch_t@fit$z)
#change column names
colnames(save_t) <- c("logidx", "s", "z")
#save the estimated parameters
parmt <- round(fitGarch_t@fit$coef,6)
parmt@shape
## test z
mean(save_t$z) #should be 0
sd(save_t$z) #should be 1
skewness(save_t$z) #should be 0
kurtosis(save_t$z) #should be equal kp
#compute the perfect kurtosis
shape <- getElement(parmt, "shape")
kp <- 3 + (6/(getElement(parmt, "shape")-4))
#check for serial correlation and volatility clustering in the Z
acf(save_t$z, lag=100)
acf(abs(save_t$z), lag=100)
#compute the VaR and ES with the GARCH model
#set seed value
set.seed(123789)
#simulate 1-day outcomes
boot <- ugarchboot(fitGarch_t, #model to use
                    method= "Partial", #method of simulation - ignore parameter uncertainty
                    sampling="raw", #draw from standardized residuals (fitted values of the epsilons)
                    n.ahead=1, #how many period to simulate: 1-day ahead
                    n.bootpred=100000, #how many simulations -number of simulated outcomes
                    solver="solnp")
#save the outcomes in the vector
vidx_GARCH <- boot@fseries
#compute the VaR
VaRGarch <- round(quantile(vidx_GARCH,0.05),6)
VaRGarch
#compute the ES
ESGarch <- round(mean(vidx_GARCH[vidx_GARCH<VaRGarch]),6)
ESGarch

## ROLLING 1-DAY VaR (for 2023, with relative plot)
#get the position of 2022-31-12 in logidx
n2022 <- length(logidx["2014-03-19/2022-12-31"])
#filter logidx to end at the end of 2023
logidx2023 <- logidx["2014-03-19/2023-12-31"]
#let the VaR roll from the 2022 to up to now
roll_garch <- ugarchroll(spec=garch_t, #garch model
           data=logidx2023, #data
           n.ahead=1, #1-day ahead var
           forecast.length=1,
           n.start=n2022, #position of last day of 2022
           refit.every=1, #re-estimate for every obs added
           refit.window = "recursive",
           calculate.VaR = TRUE, # compute the VaR
           VaR.alpha = c(0.05), #VaR alpha quantile
           keep.coef = TRUE
           )
#display the results
report(roll_garch, type = "VaR", VaR.alpha = 0.05, conf.level = 0.95)
#plot the VaR chart along with the actual returns
plot(roll_garch, type = "VaR", VaR.alpha = 0.05, conf.level = 0.95)
#save the rolling VaR and the actual value in a df
VaRroll <- as.data.frame(roll_garch, which = "VaR")
#let the "date" be a column
VaRroll <- cbind(rownames(VaRroll), VaRroll)
rownames(VaRroll) <-  NULL
colnames(VaRroll) <- c("date", "var", "actual")
#transform col "date" into a date objects
VaRroll$date <- as.Date(VaRroll$date)
#plot the VaR and actual log ret
library(ggplot2)
#set the colors
colors <- c("VaR" = "red", "log-return" = "darkblue")
VaRplot <- ggplot(data = VaRroll, aes(x=date)) +
  geom_line(aes(y=var, color = "VaR")) +
  geom_line(aes(y=actual, color = "log-return"), alpha = 0.6) +
  labs ( x = "month", y = "log return", color = "Legend", title = "1-day 95% VaR for 2023") +
  scale_color_manual(values = colors) +
  theme(panel.grid = element_line(color = "lightgrey"),
        panel.background = element_blank())
VaRplot
#count how many time the VaR was > than the real
sum(VaRroll[,2]>VaRroll[,3])
#do the same thing by first creating a column with "yes" if var > actual and then count the "yes
Varroll2mut <- dplyr::mutate(VaRroll2, worse = ifelse((var < actual), "no", "yes"))
sum(Varroll2mut$worse == "yes")
#give the same results
#compute the % of times the real was worse than the VaR
round((sum(VaRroll[,2]>VaRroll[,3])/250)*100, 2)

### e-GARCH
#set the parameters for the GARCH model with rescaled-t assumption
egarch_t <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1,1)), #variance model: GARCH(1,1)
                       mean.model = list(armaOrder = c(0,0), include.mean = TRUE), #mean model: cosntant term equals the mean
                       distribution.model = "std") #distirbution model: rescaled t distribution
## Rolling 1-day VaR for 2023 (with relative plot)
#let the VaR roll from the 2022 to up to now
roll_egarch <- ugarchroll(spec=egarch_t, #garch model
                          data=logidx2023, #data
                          n.ahead=1, #1-day ahead var
                          forecast.length=1,
                          n.start=n2022, #position of last day of 2022
                          refit.every=1, #re-estimate for every obs added
                          refit.window = "recursive",
                          calculate.VaR = TRUE, # compute the VaR
                          VaR.alpha = c(0.05), #VaR alpha quantile
                          keep.coef = TRUE
)
#display the results
report(roll_egarch, type = "VaR", VaR.alpha = 0.05, conf.level = 0.95)
#save the rolling VaR and the actual value in a df
VaReroll <- as.data.frame(roll_egarch, which = "VaR")
#let the "date" be a column
VaReroll <- cbind(rownames(VaReroll), VaReroll)
rownames(VaReroll) <-  NULL
colnames(VaReroll) <- c("date", "var", "actual")
#transform col "date" into a date objects
VaReroll$date <- as.Date(VaReroll$date)
#plot the VaR and actual log ret
library(ggplot2)
#set the colors
colors <- c("eGARCH-VaR" = "red", "log-return" = "darkblue")
eVaRplot <- ggplot(data = VaReroll, aes(x=date)) +
  geom_line(aes(y=var, color = "eGARCH-VaR")) +
  geom_line(aes(y=actual, color = "log-return"), alpha = 0.6) +
  labs ( x = "month", y = "log return", color = "Legend", title = "1-day 95% eGARCH-VaR for 2023") +
  scale_color_manual(values = colors) +
  theme(panel.grid = element_line(color = "lightgrey"),
        panel.background = element_blank())
eVaRplot
#compute the % of times the real was worse than the VaR
sum(VaReroll[,2]>VaReroll[,3])
round((sum(VaReroll[,2]>VaReroll[,3])/250)*100, 2)


####### VAR FOR A MULTI-ASSETS PORTFOLIO
library(PerformanceAnalytics)
##get the portfolio ready
#set the tickers name
tickers = c("AMZN", "AAPL", "TSLA")
#set the weights
values <- c(30000, 50000, 20000)
w1 = values/sum(values)
#weights as a list
w = c(0.3, 0.5, 0.2)
#set the dates
start <- "2020-01-01"
end <- "2024-03-01"
# get the prices (saved as the name in "")
getSymbols("AAPL", from = start, to = end, warnings = FALSE, auto.assign = TRUE)
getSymbols("AMZN", from = start, to = end, warnings = FALSE, auto.assign = TRUE)
getSymbols("TSLA", from = start, to = end, warnings = FALSE, auto.assign = TRUE)
#let's look to the first rows 
head(AAPL)
head(AMZN)
head(TSLA)
#make a df for all the data
logretstocks <- cbind(AAPL$AAPL.Adjusted, AMZN$AMZN.Adjusted, TSLA$TSLA.Adjusted)
#rename the columns
colnames(logretstocks) <- c("AAPL", "AMZN", "TSLA")
#compute daily returns
AAPLretlog <- CalculateReturns(logretstocks$AAPL, method="log")
AMZNretlog <- CalculateReturns(logretstocks$AMZN, method="log")
TSLAretlog <- CalculateReturns(logretstocks$TSLA, method="log")
#remove the first NA
AAPLretlog <- AAPLretlog[-1]
AMZNretlog <- AMZNretlog[-1]
TSLAretlog <- TSLAretlog[-1]
# add all the logret in the same database
logretport <- cbind(AAPLretlog, AMZNretlog, TSLAretlog)
#summarize the portfolio
summary(logretport)
#get the correlation matrix
cor(logretport)
#get the var-cor matrix
cov(logretport)
# a concrete look: Scatter plot, combinations of stocks
library(graphics)
pairs(~AAPL + AMZN + TSLA, data = logretport, col = "darkblue",  upper.panel = NULL)

##NORMALITY METHOD
#compute the VaR with the Gaussian method
VaRpornorm <- PerformanceAnalytics::VaR(logretport, weights = w, p=0.95, portfolio_method = "component", method = "gaussian")
VaRpornorm
#compute the ES with the Gaussian method
ESpornorm <- PerformanceAnalytics::ETL(logretport, weights = w, p=0.95, portfolio_method = "component", method = "gaussian")
ESpornorm
#compute the amount in $$
#considering the actual investment, compute the actual loss in $
loosVaRpnorm = invpor * (1- exp(VaRpornorm$VaR))
lossESpnorm = invpor * (1-exp(ESpornorm$ES))

##HISTORICAL METHOD
#compute the VaR with the historical method
VaRporhist <- PerformanceAnalytics::VaR(logretport, weights = w, p=0.95, portfolio_method = "component", method = "historical")
VaRporhist
#compute the ES with the historical method
ESporhist <- PerformanceAnalytics::ETL(logretport, weights = w, p=0.95, portfolio_method = "component", method = "historical")
ESporhist
#considering the actual investment, compute the actual loss in $
loosVaRphist = invpor * (1- exp(VaRporhist$hVaR))
lossESphist = invpor * (1-exp(ESporhist$`-r_exceed/c_exceed`))


###### PORTFOLIO RETURNS APPROACH
#single stocks' simple returns
AAPLretsimp <- CalculateReturns(logretstocks$AAPL, method="discrete")
AMZNretsimp <- CalculateReturns(logretstocks$AMZN, method="discrete")
TSLAretsimp <- CalculateReturns(logretstocks$TSLA, method="discrete")
#remove the first NA
AAPLretsimp <- AAPLretsimp[-1]
AMZNretsimp <- AMZNretsimp[-1]
TSLAretsimp <- TSLAretsimp[-1]
#create a df with all the simple returns
simpretport <- cbind(AAPLretsimp, AMZNretsimp, TSLAretsimp)
#sanity check for the weights
library(tidyverse)
tibble(w, tickers)
#set the various weights as variables
w1 <- w[1]
w2 <- w[2]
w3 <- w[3]
#compute the daily return of the portfolio assuming daily re-balancing (this keeps the weight fixed)
pretd <- (AAPLretsimp*w1) +
  (AMZNretsimp*w2) +
  (TSLAretsimp*w3)
#renaming the column
names(pretd) <- "PortfolioReturn"
#check with a formula for port returns
Return.portfolio(simpretport, weights = w, rebalance_on = "day")
#results are the same 
#compute the log-return of the portfolio
logRetPor <- log(1+pretd$PortfolioReturn)
#rename the column as logPortReturn
names(logRetPor)[names(logRetPor) == "PortfolioReturn"] <- "logPortReturn"
#plot the return
plot(logRetPor, col = "darkblue")
#Normality test on returns
#test for normality
jarque.test(as.vector(logRetPor)) #not normal
#not normally distributed
skewness(as.vector(logRetPor)) #-0.3172239
kurtosis(as.vector(logRetPor)) #5.22476
#compute the t for the return
library(MASS)
t_por <- (fitdistr(logRetPor, "t"))
round(t_por$estimate,6)
#m = mean, s = standard deviation (or rescaling parameter),  df = degree of freedom
#set the seed
set.seed(123789)
#call the library
library(metRology)
#set the standardized student-t parameters values (mean, sd and degree of freedom)
logRetPor2 <- rt.scaled(100000,mean=t_por$estimate[1],sd=t_por$estimate[2],df=t_por$estimate[3])
#estimate VaR
VaR_port <- quantile(logRetPor2, alpha)
round(VaR_port, 6)
#estimate the ES
ES_port <- mean(logRetPor2[logRetPor2<VaR_port])
round (ES_port, 6)
#check for volatility clustering
acf(abs(logRetPor), lag=100, main = "Portfolio absolute Log return acf") #yes
#and serial correlation
acf(logRetPor, lag=100, main = "Portfolio Log return acf") #no

#GARCH model
#set the parameters for the GARCH model with rescaled-t assumption
garch_por <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), #variance model: GARCH(1,1)
                      mean.model = list(armaOrder = c(0,0), include.mean = TRUE), #mean model: cosntant term equals the mean
                      distribution.model = "std") #distirbution model: rescaled t distribution
#estimate the GARCH model
fitGarch_por <- ugarchfit(spec = garch_por, data = logRetPor)
fitGarch_por
#Mu = a sub 0 in the mean eqaution ; omega = alpha_0 in the variance eqaution ; alpha1 = alpha_1 in variance equation ; Beta1 = Beta_1 in the variance equation
#save some outputs (dataframe with logreturns, fitted value of daily sd, fitted vaule of Epsilons in the distribution equation)
save_tpor <- cbind(logRetPor,fitGarch_por@fit$sigma, fitGarch_por@fit$z)
#change column names
colnames(save_tpor) <- c("logRetPor", "s", "z")
#save the estimated parameters
parmtpor <- round(fitGarch_por@fit$coef,6)
parmtpor["shape"]
## test z
mean(save_tpor$z) #should be 0
sd(save_tpor$z) #should be 1
skewness(save_tpor$z) #should be 0
kurtosis(save_tpor$z) #should be equal kp
#compute the perfect kurtosis
shape <- getElement(parmtpor, "shape")
kp <- 3 + (6/(getElement(parmtpor, "shape")-4))
#check for serial correlation and volatility clustering in the Z
acf(save_tpor$z, lag=100, main = "Z's Log return acf")
acf(abs(save_tpor$z), lag=100, main = "Z's absolute Log return acf")
#compute the VaR and ES with the GARCH model
#set seed value
set.seed(123789)
#simulate 1-day outcomes
bootpor <- ugarchboot(fitGarch_por, #model to use
                   method= "Partial", #method of simulation - ignore parameter uncertainty
                   sampling="raw", #draw from standardized residuals (fitted values of the epsilons)
                   n.ahead=1, #how many period to simulate: 1-day ahead
                   n.bootpred=100000, #how many simulations -number of simulated outcomes
                   solver="solnp")
#save the outcomes in the vector
retpor_GARCH <- bootpor@fseries
#compute the VaR
VaRGarchpor <- round(quantile(retpor_GARCH,0.05),6)
VaRGarchpor
#compute the ES
ESGarchpor <- round(mean(retpor_GARCH[retpor_GARCH<VaRGarchpor]),6)
ESGarchpor
#compute the actual loss
VaRGarchpor_act <- invpor * (exp(VaRGarchpor)-1)
ESGarchpor_act <- invpor * (exp(ESGarchpor)-1)
#show them
VaRGarchpor_act
ESGarchpor_act
## Rolling 1-day VaR for 2023 and up to now (with relative plot)
#get the position of 2022-31-12 in logidx
por2022 <- length(logRetPor["2014-03-19/2022-12-31"])
#let the VaR roll from the 2022 to up to now
roll_garchpor <- ugarchroll(spec=garch_por, #garch model
                         data=logRetPor, #data
                         n.ahead=1, #1-day ahead var
                         forecast.length=1,
                         n.start=por2022, #position of last day of 2022
                         refit.every=1, #re-estimate for every obs added
                         refit.window = "recursive",
                         calculate.VaR = TRUE, # compute the VaR
                         VaR.alpha = c(0.05), #VaR alpha quantile
                         keep.coef = TRUE
)
#display the results
report(roll_garchpor, type = "VaR", VaR.alpha = 0.05, conf.level = 0.95)
#plot the VaR chart along with the actual returns
plot(roll_garchpor, type = "VaR", VaR.alpha = 0.05, conf.level = 0.95)
#save the rolling VaR and the actual value in a df
VaRrollpor <- as.data.frame(roll_garchpor, which = "VaR")
#let the "date" be a column
VaRrollpor <- cbind(rownames(VaRrollpor), VaRrollpor)
rownames(VaRrollpor) <-  NULL
colnames(VaRrollpor) <- c("date", "var", "actual")
#transform col "date" into a date objects
VaRrollpor$date <- as.Date(VaRrollpor$date)
#plot the VaR and actual log ret
library(ggplot2)
#set the colors
colors <- c("VaR" = "red", "log-return" = "darkblue")
VaRporplot <- ggplot(data = VaRrollpor, aes(x=date)) +
  geom_line(aes(y=var, color = "VaR")) +
  geom_line(aes(y=actual, color = "Portfolio-log-return"), alpha = 0.6) +
  labs ( x = "month", y = "log return", color = "Legend", title = "1-day 95% VaR for 2023") +
  scale_color_manual(values = colors) +
  theme(panel.grid = element_line(color = "lightgrey"),
        panel.background = element_blank())
VaRporplot
#count how many time the VaR was > than the real
sum(VaRrollpor[,2]>VaRrollpor[,3])
#compute the % of times the real was worse than the VaR
round((sum(VaRrollpor[,2]>VaRrollpor[,3])/307)*100, 2)


## FILTERED HISTORICAL SIMULATION
#load the library
library(quarks)
#compute the filtered historical simulation  with GARCH(1,1)
FHS <- fhs(logRetPor, p = 0.95, 
    model = "GARCH", 
    nboot=100000, 
    lambda = 0.94,
    variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
    mean.model=list(armaOrder = c(0,0)))
#show the VaR
FHS$VaR
#show the ES
FHS$ES
#actual loss
#considering the actual investment, compute the actual loss in $
loosVaRfhs = invpor * (1- exp(FHS$VaR))
loosVaRfhs
lossESfhs = invpor * (1-exp(FHS$ES))
lossESfhs

