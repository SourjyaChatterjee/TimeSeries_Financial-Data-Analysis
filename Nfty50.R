library(PerformanceAnalytics)
library(astsa)
library(itsmr)
library(lubridate)
library(zoo)
library(randtests)
library(forecast)
library(urca)
library(aTSA)
library(ggplot2)
library(tsoutliers)
library("gridExtra")
library(parallel)
library(rugarch)
library("readxl")
library(nlts)


##### we will use data from 23.03.2020 to 14.09.2021 to train our model 
#####and then will try to forcast and compare the prices of the next 10 days.

#nfty50 <- read.csv('NIFTY 50_Data.csv')
nfty50 <- read.csv('NIFTY 50_Data.csv')

nfty50 <- nfty50[seq(1,375),c(1,5)]

nfty50[,1] <- dmy(nfty50[,1])

par(mfrow = c(1,1))
plot(rev(nfty50$Close), ylab="Stock Prices",main="Figure : Closing prices of the stocks",type = 'l')
nfty50
tso <- zoo(rev(nfty50$Close), rev(nfty50$Date))
tso
plot(tso,ylab="Stock Prices",xlab ='Time', main="Figure : Closing prices of the stocks",type = 'l')
####################
#####log return#####
####################
Return=CalculateReturns(tso, method = 'log')
plot(Return,main='Log-Return of NIFTY',xlab='Time')
Return
plot(Return^2,main='Squre of  Log-Return of NIFTY',xlab='Time')
nfty50_return = Return
#####################################
#### Augmented Dickey Fuller Test ###
#####################################

summary(ur.df(na.omit(as.vector(nfty50_return))))

######Skewness Kurtois ##############
ggplot(aes(as.vector(na.omit(Return))), data=na.omit(Return)) + geom_histogram(bins = 100,col='black',fill='red') + ggtitle('Distribution of Log-Returns') + labs(x= 'Log-Return')
skewness((as.vector(na.omit(Return)))); kurtosis((as.vector(na.omit(Return))))

############## QQ Plot ##############
ggplot(data=nfty50, aes(sample = as.vector(Return))) +
  stat_qq() +
  stat_qq_line(col='red') + ggtitle('QQ plot of Nifty Returns')+  labs(x= 'Normal Distribution Quantile', y='Nifty Data Quantile')

#######  Normality Test ##############

shapiro.test(na.omit(as.vector(Return)))

############## iid test ##############
Box.test(na.omit(as.vector(Return)), type = "Ljung-Box")
######################################

#Absolute Return or Squared of Return are auto correlated.

a<- ggAcf(abs(na.omit(as.vector(Return))), col='red',main='Acf of Absolute Return of NIFTY')
p<- ggPacf(abs(na.omit(as.vector(Return))),col='steelblue',main='PAcf of Absolute Return of NIFTY')
grid.arrange(a,p, ncol = 2, nrow = 1)


c <- ggAcf(na.omit(as.vector(Return))^2, lag.max = 40, col='red', main='ACF of squared Return Values')
d <- ggPacf(na.omit(as.vector(Return))^2,lag.max = 40, col='steelblue',main= 'PACF of squared Return Values')
grid.arrange(c,d, ncol = 2, nrow = 1)

###############################
#### Volatility Clustering ####
###############################
chart.RollingPerformance(na.omit(Return),width = 22,FUN = 'sd.annualized',scale=252, main = 'Rolling 1 month Volatility')

############## Remarks  #############
#Stylized Facts of Financial Data
#Distributions of Returns is not normal.
#Absence of significant auto correlation in returns.
#Slowly decreasing auto correlation in squared or absolute returns.
#Volatility clustering.
#####################################

##################################################################################
################################# GARCH Model ####################################
##################################################################################

##################################################################################
## Model Selection ########
##################################################################################

aic = c()
bic = c()
aicc = c()

mean_first = c()
mean_second = c()
var_first = c()
var_second = c()

for (i in 0:4){
  for (j in 0:4){
    for (k in 0:4){
      for (l in 0:4){
        try({
          
          n50_garch <- ugarchspec(mean.model = list(armaOrder=c(i,j)),variance.model = list(model = 'sGARCH', 
                                                                                            garchOrder = c(k,l)), distribution = 'std')
          fit_garch <- ugarchfit(spec = n50_garch, data= na.omit(Return))
          
          aic = c(aic,infocriteria(fit_garch)[1])
          bic = c(bic,infocriteria(fit_garch)[2])
          mean_first = c(mean_first,i)
          mean_second = c(mean_second,j)
          var_first = c(var_first,k)
          var_second = c(var_second,l)
          
          
          aic_score = infocriteria(fit_garch)[1]
          n = length(Return)
          m = i+j+k+l
          
          aicc_score = aic_score +  ((2*m^2+2*m)/(n-m-1))
          aicc = c(aicc,aicc_score)
        })
      }
    }
  }
}

length(mean_first)
length(mean_second)
var_first
var_second
length(aic)
length(aicc)

usd_result <- data.frame(mean_first,mean_second,var_first,var_second,aic,bic,aicc)
usd_result
colnames(usd_result) <- c("1st mean Parameter","2nd mean Parameter","1st variance Parameter","2nd variance Parameter","AIC Scores","BIC Scores","AICc Scores")
write.csv(usd_result,'usd_model.csv')

which(usd_result$`AIC Scores` == min(usd_result$`AIC Scores`))
which(usd_result$`BIC Scores` == min(usd_result$`BIC Scores`))
which(usd_result$`AICc Scores` == min(usd_result$`AICc Scores`))

#so, our final model is arma(3,2)-garch(4,1) for aic
n50_garch_1 <- ugarchspec(mean.model = list(armaOrder=c(0,0)),variance.model = list(model = 'sGARCH', 
                                                                                    garchOrder = c(1,1)),distribution = 'std')
fit_garch_1 <- ugarchfit(spec = n50_garch_1, data= na.omit(Return))

coef(fit_garch_1)

#so, our final model is arma(0,0)-garch(1,1) for aicc
n50_garch_2 <- ugarchspec(mean.model = list(armaOrder=c(0,0)),variance.model = list(model = 'sGARCH', 
                                                                                    garchOrder = c(1,1)),distribution = 'std')
fit_garch_2 <- ugarchfit(spec = n50_garch_2, data= na.omit(Return))

coef(fit_garch_2)

#################################
#### Convergence of the Model ###
#################################

print(convergence(fit_garch_1))   # The model converge
print(convergence(fit_garch_2))   # The model converge

#################################
######### Forecasting ###########
#################################

for_cast1 <-ugarchforecast(fit_garch_1,data=tso,n.ahead=10)
for_cast1

for_cast2 <-ugarchforecast(fit_garch_2,data=tso,n.ahead=10)
for_cast2

##################################
######### Rolling Forecast #######
##################################

##### 1st Model #####
fit_roll <- ugarchfit(n50_garch_1, data= na.omit(Return),out.sample =10)
fore_roll <- ugarchforecast(fit_roll, n.ahead=10, n.roll=10)
fore_roll
par(mfrow=c(1,2))
plot(fore_roll,which=1)
plot(fore_roll,which=2)

##### 2nd Model #####
fit_roll <- ugarchfit(n50_garch_2, data= na.omit(Return),out.sample =10)
fore_roll <- ugarchforecast(fit_roll, n.ahead=10, n.roll=10)
fore_roll
par(mfrow=c(1,2))
plot(fore_roll,which=1)
plot(fore_roll,which=2)

############################################################
####### Compute MSE for the Prediction of the Model ########
############################################################

a2 <- fore_roll@forecast
pred2 <- a2$seriesFor[,1]
sum((Return[(length(Return)-9):length(Return)] - pred2)^2)/10
