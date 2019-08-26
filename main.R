######### Setup #########
rm(list = ls())
getwd()
source('dataloader.R')

# Cria fatores para ordenação dos Regimes adotados 
# de acordo com cada Resolução
dados <- dados %>% 
  mutate(regime_fac  = factor(regime, levels = unique(regime),
                          labels = paste('Regime', 1:10), ordered = T))
attach(dados)

## Package loading 
if(!require(ggplot2)) {install.packages('ggplot2'); require('ggplot2')}
if(!require(tidyquant)) {install.packages('tidyquant'); require('tidyquant')}
if(!require(ggpmisc)) {install.packages('ggpmisc'); require('ggpmisc')}
theme_set(theme_bw())

## Análise Descritiva
p <- ggplot(dados, aes(Data1, IPCA)) + geom_line() +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  xlab("") + ylab("IPCA") 

p + stat_smooth(color = "#FC4E07", fill = "#FC4E07",
  method = "loess")

ggplot(dados, aes(x=regime_fac, y=IPCA, color=regime)) +
  geom_boxplot() +
  geom_jitter()

## Previsão SARIMA
library(tseries)
library(forecast)

## Define série temporal como Log-IPCA
yt <- ts(c(IPCA, 4.94, 4.66, 3.37, 3.22), start = c(1999, 12), end = c(2019, 7), frequency = 12)
zt <- window(yt, start = c(1999, 12), end = c(2018, 7), frequency = 12 )
tseries::adf.test(zt) 
auto.arima(zt)

## Teste de Dickey Fuller
# Rejeita-se a hipótese nula de que os dados não são estacionários

if(!require(ggfortify)) {install.packages('ggfortify'); require('ggfortify')}
if(!require(gridExtra)) {install.packages('gridExtra'); require('gridExtra')}
if(!require(forecast)) {install.packages('forecast'); require('forecast')}
if(!require(DescTools)) {install.packages('DescTools'); require('DescTools')}
source('functions.R')

fit <- arima(zt, order = c(1,1,0), seasonal = c(2,0,1), optim.control = list(maxit = 1000))
summary(fit)

tsdiag(fit)
pred <- forecast::forecast(fit, h = 12, level = .95)
 
date <-c(timetk::tk_index(pred$x), timetk::tk_index(pred$mean))
plot(date, yt, pch = 20, bty = 'n', xlab = '', ylab = 'IPCA')
lines(date, c(pred$fitted, pred$mean), col = 2, lwd = 2)
abline(v = date[224], lty = 2, lwd = 2)

pto.parada = 224
graf_previsao(yt,
               as.vector(pred$fitted) , 
               pto.parada,
               main = "prev_12passos_sarima.pdf",
               interv = 0,
               LS_prev = c( rep(NA, pto.parada), pred$upper),
               LI_prev = c( rep(NA, pto.parada), pred$lower)
              )





