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
library(FinTS)

##### we will use data from 23.03.2020 to 12.10.2021 to train our model 
#####and then will try to forcast and compare the prices of the next 10 days.

gold <- read_excel('gold.xlsx') # Locatoon of the DataSet

nfty50 <- gold[,c(1,2)]

nfty50$Date <- mdy(nfty50$Date)
nfty50

par(mfrow = c(1,1))
plot(nfty50, ylab="Gold Prices",type = 'l')
tso <- zoo(rev(nfty50$Price), rev(nfty50$Date))
plot(as.Date(index(tso)),coredata(tso), ylab="Gold Prices", xlab = "Time",main="Prices of Gold in INR",type = 'l')
####################
#####log return#####
####################
Return=CalculateReturns(tso, method = 'log')
plot(Return,main='Return of Gold',xlab='Time')
Return
plot(Return^2,main='Sqd. Return of Gold',xlab='Time')
gold_return = Return
#####################################
#### Augmented Dickey Fuller Test ###
#####################################

summary(ur.df(na.omit(gold_return)))

######Skewness Kurtois ##############
ggplot(aes(as.vector(na.omit(Return))), data=na.omit(Return)) + geom_histogram(bins = 100,col='black',fill='red') + ggtitle('Distribution of Log-Returns') + labs(x= 'Log-Return')
skewness((as.vector(na.omit(Return)))); kurtosis((as.vector(na.omit(Return))))
############## QQ Plot ##############
ggplot(data=nfty50, aes(sample = as.vector(Return))) +
  stat_qq() +
  stat_qq_line(col='red') + ggtitle('QQ plot of Gold Returns')+  labs(x= 'Normal Distribution Quantile', y='Nifty Data Quantile')
#######  Normality Test ##############
#jarque.bera.test(na.omit(as.vector(Return)))
shapiro.test(na.omit(as.vector(Return)))
############## iid test ##############
Box.test(na.omit(as.vector(Return)), type = "Ljung-Box")
######################################

#Absolute Return or Squared of Return are auto correlated.

a<- ggAcf(abs(na.omit(as.vector(Return))), col='red',main='Acf of Absolute Return of Gold')
p<- ggPacf(abs(na.omit(as.vector(Return))),col='steelblue',main='PAcf of Absolute Return of Gold')
grid.arrange(a,p, ncol = 2, nrow = 1)


c <- ggAcf(na.omit(as.vector(Return))^2, lag.max = 40, col='red', main='ACF of squared Return Values')
d <- ggPacf(na.omit(as.vector(Return))^2,lag.max = 40, col='steelblue',main= 'PACF of squared Return Values')
grid.arrange(c,d, ncol = 2, nrow = 1)

###############################
#### Volatility Clustering ####
###############################
chart.RollingPerformance(na.omit(Return^2),width = 22,FUN = 'sd.annualized',scale=252, main = 'Rolling 1 month Volatility')

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

gold_result <- data.frame(mean_first,mean_second,var_first,var_second,aic,bic,aicc)
gold_result
colnames(gold_result) <- c("1st mean Parameter","2nd mean Parameter","1st variance Parameter","2nd variance Parameter","AIC Scores","BIC Scores","AICc Scores")
write.csv(gold_result,'Gold_model.csv')


which(gold_result$`AIC Scores` == min(gold_result$`AIC Scores`))
which(gold_result$`BIC Scores` == min(gold_result$`BIC Scores`))
which(gold_result$`AICc Scores` == min(gold_result$`AICc Scores`))
#####################################################################################
#####################################################################################

#so, our final model is arma(0,0)-garch(1,1) according to aicc & bic score
n50_garch_1 <- ugarchspec(mean.model = list(armaOrder=c(0,0)),variance.model = list(model = 'sGARCH', 
                                                                                  garchOrder = c(0,1)), distribution = 'std')
fit_garch_1 <- ugarchfit(spec = n50_garch_1, data= na.omit(gold_return))

coef(fit_garch_1)


#so, our final model is arma(4,4)-garch(1,1) according to aic score
n50_garch_2 <- ugarchspec(mean.model = list(armaOrder=c(4,4)),variance.model = list(model = 'sGARCH', 
                                                                                    garchOrder = c(1,1)),distribution = 'std')
fit_garch_2 <- ugarchfit(spec = n50_garch_2, data= na.omit(gold_return))

fit_garch_2
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

##### 2nd Model ####

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
