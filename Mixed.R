library("readxl")
library(ggplot2)
library(corrplot)

nfty50 <- read.csv('NIFTY 50_Data.csv')
nfty50 <- rev(nfty50[seq(1,375),c(1,5)])

gold <- read_excel('gold.xlsx')
gold <- rev(gold[,c(1,2)])

usd <- read_excel('usd.xlsx')
usd <- rev(usd[,c(1,2)])

silver <- read_excel('silver.xlsx')
silver <- rev(silver[,c(1,2)])

######################################
######### Correlation ################
######################################

nfty50_return # Log-return of Nfty50 stock price
gold_return # Log-return of Gold Price
usd_return # Log-return of USD
silver_return # Log-return of Silver

library(corrplot)

data = data = cbind(as.numeric(as.vector(na.omit(nfty50$Close))),as.numeric(as.vector(na.omit(gold$Price))),as.numeric(as.vector(na.omit(usd$Price))),as.numeric(as.vector(na.omit(silver$Price))))

data = data.frame(data)
colnames(data) <- c("Nfty50","Gold","USD","Silver")
corrplot(cor(data),        # Correlation matrix
         method = "shade", # Correlation plot method
         type = "full",    # Correlation plot style (also "upper" and "lower")
         diag = TRUE,      # If TRUE (default), adds the diagonal
         tl.col = "black", # Labels color
         bg = "white",     # Background color
         title = "",       # Main title
         col = NULL)       # Color palette

library(psych)

corPlot(data, cex = 1.2)
