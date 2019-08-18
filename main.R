######### Setup #########
rm(list = ls()); getwd()
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

## Define série temporal como Log-IPCA
y <- IPCA
lambda <- BoxCox.lambda(y, method = 'guerrero')

z <- ts((y^(-lambda) - 1)/(-lambda), start = c(1999, 1), end = c(2019, 3), frequency = 12)
plot(z)

## Teste de Dickey Fuller
# Rejeita-se a hipótese nula de que os dados não são estacionários
adf.test(z) 

if(!require(ggfortify)) {install.packages('ggfortify'); require('ggfortify')}
if(!require(gridExtra)) {install.packages('gridExtra'); require('gridExtra')}
if(!require(forecast)) {install.packages('forecast'); require('forecast')}
if(!require(DescTools)) {install.packages('DescTools'); require('DescTools')}
source('functions.R')
fac_facp(z, top="IPCA(1999-2019)", lag=60, titulo="")
auto.arima(z)

fit <- arima(z, order = c(1,1,2), seasonal = c(2,0,1), optim.control = list(maxit = 1000))
# summary(fit)

ljungbox(fit, 6)
norm_test(fit, titulo = '')
shapiro.test(fit$residuals)
hist(fit$residuals, freq = F)
curve(dnorm(x), add = T)
previsao(fit, 4, alpha =.95)

plot(fit$residuals)
n <- length(z)
abline(h = c(-1.96/sqrt(n),1.96/sqrt(n)))

previsao = forecast::forecast(fit, h = 4, level = .95)  
tt <- function(x, lambda){
  (1-lambda*x)^(-1/lambda)
}
previsao$lower <-  tt(previsao$lower, lambda)
previsao$upper <-  tt(previsao$upper, lambda)

previsao$mean <-  tt(previsao$mean, lambda)
previsao$x <-  tt(previsao$x, lambda)
previsao

autoplot(previsao)
