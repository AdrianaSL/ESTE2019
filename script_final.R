
options(scipen = 9999)

### ================================================================
### ================================================================
### Reading packages
### ================================================================
### ================================================================

library(readxl)
library(Metrics)
library(dlm)
library(readxl)
library(expm)
library(dplyr)
library(metRology)
library(tseries)
library(forecast)
if(!require(ggfortify)) {install.packages('ggfortify'); require('ggfortify')}
if(!require(gridExtra)) {install.packages('gridExtra'); require('gridExtra')}
if(!require(forecast)) {install.packages('forecast'); require('forecast')}
if(!require(DescTools)) {install.packages('DescTools'); require('DescTools')}

### ================================================================
### ================================================================
### Dataset 
### ================================================================
### ================================================================

#1- Setup #########
getwd()

# 2- Dados Históricos #########
## Dados mensais do Sistema IBGE de Recuperação Automática (SIDRA)
dados <- read.csv2("ipca_201903SerieHist.csv", header = T, dec = ",", sep = ";")
dados <- dados %>% filter( ano > 1995 )

dados$Data2 <- paste0(tolower(dados$mes),"/",dados$ano)

aux <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
dados$Data1 <- as.Date(paste0(rep("01", 279),
                              c(rep(aux, 23), aux[1:3]),
                              dados$ano), '%d%m%Y')

dados <- dados %>% filter( as.Date(dados$Data1) > "1999-11-01" )

y <- dados$IPCA 

y = c(y, 4.94, 4.66, 3.37, 3.22)

### ================================================================
### ================================================================
### Loading scripts
### ================================================================
### ================================================================


source("function3.R")
source("Harmonicos selecionados.R")


### ================================================================
### ================================================================
### Modeling
### ================================================================
### ================================================================


### ================================================================
### MODELO DINAMICO ================================================
### ================================================================


### Definindo F (vetor coluna de covariaveis), G (matriz de evolucao) e D (matriz do inverso do desconto)

G = dlm::bdiag(G1,G2[[1]],G2[[2]])               #G1 (tendencia), G2[[1]] (harmonico 1 -> p = 12) e G2[[2]] (harmonico 2 -> p = 6)
FF = matrix(c(1,0,1,0,1,0))                          
F.t  = matrix(FF,nrow=length(FF),ncol=length(y))
D1 <- matrix(1/0.7,2,2)                          #D1 (desconto de 0,8)
D2 <- matrix(1/1,2,2)                            #D2 (desconto de 1, i.e., nao tem desconto - estatico)
D = dlm::bdiag(D1,D2,D2)                         #D1 (desconto na tendencia), D2 (desconto nos dois harmonicos)

### Definindo as quantidades iniciais ======
n<-nrow(F.t)
m0 <- rep(0,n)
C0 <- diag(100000,n,n) 
n0 <- 0.1
d0 <- .001

# ### Definir ponto de parada =====
# pto.parada <- length(y) - 12 # previsao a partir de pto.parada

### Modelo com intervencao =====

# intervencao nos meses: "jun/2000", "jun/2001", "jun/2002", "jun/2003", "jun/2004", "jun/2005", "jun/2006", "jun/2017"

int_date = c("jun/2000", "jun/2001", "jun/2002", "jun/2003", "jun/2004", "jun/2005", "jun/2006", "jun/2017")
for( i in 1:length(int_date)){
  print( dados[which(dados$Data2==int_date[i]),] )
}

int = 12*c(1:6,17) + 7
int = c(7, int)

sum( dados$Data2[int] == int_date ) == 8


# modelos: usando horizonte de previsao ( 3, 6, 9 e 12 meses )
resultados = list()
pto.parada = NULL

for( i in 1:4 ){
  
  pto.parada[i] <- length(y) - (3*i) # previsao a partir de pto.parada
  
  resultados[[i]] <- modelo.2ordem.var.desc_prev.k_intervencao(
    y, m0, C0, n0, d0, Wt = NULL, F.t, G, D,
    h0 = pto.parada[i], int = int, desc = 100
  )
  
}

### o que vamos usar eh o resultado do modelo com horizonte de previsao de 12 passos a frente
result = resultados[[4]] 
pto.parada = pto.parada[4]

### ================================================================
### MODELO SARIMA ==================================================
### ================================================================

## Define série temporal como Log-IPCA
yt <- ts(y, start = c(1999, 12), end = c(2019, 7), frequency = 12)
zt <- window(yt, start = c(1999, 12), end = c(2018, 7), frequency = 12 )


tseries::adf.test(zt) 
## Teste de Dickey Fuller
# Rejeita-se a hipótese nula de que os dados não são estacionários

auto.arima(zt)
# modelo escolhido
# ARIMA(1,1,0)(2,0,1)[12] 

fit <- arima(zt, order = c(1,1,0), seasonal = c(2,0,1), optim.control = list(maxit = 1000))
summary(fit)
tsdiag(fit)

### predicao
pred_sarima <- forecast::forecast(fit, h = 12, level = .95)

### ================================================================
### ================================================================
### Diagnostico do modelo
### ================================================================
### ================================================================

# ### diagnostico do modelo -> k = 3, 6, 9, 12 =====
# 
# pred = list()
# real = list()
# pred_power = list()
# mape = list()
# 
# for( i in 1:4 ){
#   
#   pred[[i]] = resultados[[i]]$mt[1,][(pto.parada[i]+1):length(y)]
#   real[[i]] = y[(pto.parada[i]+1):length(y)]
#   
#   pred_power[[i]] = QPS( pred[[i]], real[[i]] )
#   mape[[i]] = MAPE( pred[[i]], real[[i]] )*100
# }
# 
# pred_in = list()
# real_in = list()
# pred_power_in = list()
# mape_in = list()
# 
# for( i in 1:4 ){
#   
#   pred_in[[i]] = resultados[[i]]$mt[1,][1:pto.parada[i]]
#   real_in[[i]] = y[1:pto.parada[i]]
#   
#   pred_power_in[[i]] = QPS( pred_in[[i]], real_in[[i]] )
#   mape_in[[i]] = MAPE( pred_in[[i]], real_in[[i]] )*100
# }


# avaliando k = 3, 6, 9, 12 passos a frente, em um modelo que faz 12 passos a frente


k = c(3, 6, 9, 12)

### erro de previsao -> dinamico ==========
pred = list()
real = list()
pred_power = list()
mape = list()

for( i in 1:4 ){
  
  pred[[i]] = result$ft[(pto.parada+1):(pto.parada + k[i])]
  real[[i]] = y[(pto.parada+1):(pto.parada + k[i])]
  
  pred_power[[i]] = QPS( pred[[i]], real[[i]] )
  mape[[i]] = MAPE( pred[[i]], real[[i]] )*100
  
  names(pred_power)[i] = paste0("k = ",3*i)
  names(mape)[i] = paste0("k = ",3*i)
  
}

mape
# 
# $`k = 3`
# [1] 2.633223
# 
# $`k = 6`
# [1] 17.13106
# 
# $`k = 9`
# [1] 18.85146
# 
# $`k = 12`
# [1] 31.04113

pred_power

# > pred_power
# $`k = 3`
# [1] 0.03418701
# 
# $`k = 6`
# [1] 1.551302
# 
# $`k = 9`
# [1] 1.729843
# 
# $`k = 12`
# [1] 4.439597


### erro de estimacao -> dinamico ==========
pred_in = NULL
real_in = NULL
pred_power_in = NULL
mape_in = NULL


pred_in = result$ft[ 1:pto.parada ]
real_in = y[ 1:pto.parada ]

pred_power_in = QPS( pred_in, real_in )
mape_in = MAPE( pred_in, real_in )*100

mape_in

# > mape_in
# [1] 6.738137

pred_power_in

# > pred_power_in
# [1] 1.488532


### erro de previsao -> SARIMA ==========

pred1 = list()
real1 = list()
pred_power1 = list()
mape1 = list()

for( i in 1:4 ){
  
  pred1[[i]] = pred_sarima$mean[ 1:k[i] ]
  real1[[i]] = y[(pto.parada+1):(pto.parada + k[i])]
  
  pred_power1[[i]] = QPS( pred1[[i]], real1[[i]] )
  mape1[[i]] = MAPE( pred1[[i]], real1[[i]] )*100
  
  names(pred_power1)[i] = paste0("k = ",3*i)
  names(mape1)[i] = paste0("k = ",3*i)
  
}

mape1

# > mape1
# $`k = 3`
# [1] 7.687257
# 
# $`k = 6`
# [1] 25.35886
# 
# $`k = 9`
# [1] 33.39693
# 
# $`k = 12`
# [1] 44.80587

pred_power1
# > pred_power1
# $`k = 3`
# [1] 0.2339105
# 
# $`k = 6`
# [1] 2.989593
# 
# $`k = 9`
# [1] 5.160667
# 
# $`k = 12`
# [1] 7.944424

### erro de estimacao -> SARIMA ==========
pred_in1 = NULL
real_in1 = NULL
pred_power_in1 = NULL
mape_in1 = NULL


pred_in1 = pred_sarima$fitted
real_in1 = y[ 1:pto.parada ]

pred_power_in1 = QPS( pred_in1, real_in1 )
mape_in1 = MAPE( pred_in1, real_in1 )*100

mape_in1
# > mape_in1
# [1] 3.676957

pred_power_in1
# > pred_power_in1
# [1] 0.2002873

### ================================================================
### ================================================================
### Graphics
### ================================================================
### ================================================================

# LS_prev = list()
# LI_prev = list()
# 
# for( i in 1:4){
#   LS_prev[[i]] <- qt.scaled(0.975,resultados[[i]]$nt[pto.parada,], resultados[[i]]$ft, sqrt(resultados[[i]]$Qt))
#   LI_prev[[i]] <- qt.scaled(0.025,resultados[[i]]$nt[pto.parada,], resultados[[i]]$ft, sqrt(resultados[[i]]$Qt))
#   
#   title = paste0("prev_3passos_interv_k=",3*i,"_3.pdf")
#   graf_previsao( resultados[[i]], 
#                  pto.parada[i],
#                  main = title,
#                  interv = int,
#                  LS_prev = LS_prev[[i]],
#                  LI_prev = LI_prev[[i]])
# }


### DINAMICO ==============================================

LS_prev <- qt.scaled(0.975,result$nt[pto.parada,], result$ft, sqrt(result$Qt))
LI_prev <- qt.scaled(0.025,result$nt[pto.parada,], result$ft, sqrt(result$Qt))

graf_estimado(y = y[1:pto.parada],
              y_hat = result$ft[1:pto.parada],
              inic = 1,
              main = "estimado_dlm_interv.pdf",
              interv = int,
              IC = T,
              LS_prev = LS_prev[1:pto.parada] ,
              LI_prev = LI_prev[1:pto.parada])

graf_previsao( y = y,
               y_hat  = result$ft,
               inic = 6,
               pto.parada = pto.parada,
               main = "previsao_dlm_interv.pdf",
               interv = int,  
               LS_prev = LS_prev,
               LI_prev = LI_prev)

graf_previsao( y = y,
               y_hat = result$ft,
               inic = 200,
               pto.parada = pto.parada,
               main = "previsao_dlm_interv_zoom200.pdf",
               interv = int,  
               LS_prev = LS_prev,
               LI_prev = LI_prev)


### SARIMA =================================================

graf_estimado(y = y[1:pto.parada],
              y_hat = as.vector(pred_sarima$fitted),
              inic = 1,
              main = "estimado_sarima.pdf",
              interv = NULL,
              IC = FALSE)

graf_previsao( y = y,
               y_hat  = c(pred_sarima$fitted, pred_sarima$mean),
               inic = 6,
               pto.parada = pto.parada,
               main = "previsao_sarima.pdf",
               interv = NULL,
               IC = T,
               LS_prev = c( rep(NA, pto.parada), pred_sarima$upper),
               LI_prev = c( rep(NA, pto.parada), pred_sarima$lower))

graf_previsao( y = y,
               y_hat  = c(pred_sarima$fitted, pred_sarima$mean),
               inic = 200,
               pto.parada = pto.parada,
               main = "previsao_sarima_zoom200.pdf",
               interv = NULL,
               IC = T,
               LS_prev = c( rep(NA, pto.parada), pred_sarima$upper),
               LI_prev = c( rep(NA, pto.parada), pred_sarima$lower))


### DINAMICO X SARIMA ==============================================

# com intervalo de previsao para todas as observacoes no modelo dinamico
graf_dlmXsarima(y = y,
                y_hat1 = result$ft,
                y_hat2 = c(pred_sarima$fitted, pred_sarima$mean),
                name_yhat1 = "dinâmico",
                name_yhat2 = "SARIMA",
                inic = 200,
                main = "previsao_dlmXsarima_zoom200.pdf",
                interv = int,
                IC = T,
                LS_prev1 = LS_prev,
                LI_prev1 = LI_prev,
                LS_prev2 = c( rep(NA, pto.parada), pred_sarima$upper),
                LI_prev2 = c( rep(NA, pto.parada), pred_sarima$lower))

# com intervalo de previsao somente a partir do ponto de parada
graf_dlmXsarima(y = y,
                y_hat1 = result$ft,
                y_hat2 = c(pred_sarima$fitted, pred_sarima$mean),
                name_yhat1 = "dinâmico",
                name_yhat2 = "SARIMA",
                inic = 200,
                main = "previsao_dlmXsarima_zoom200_2.pdf",
                interv = int,
                IC = T,
                LS_prev1 = c( rep(NA, pto.parada), LS_prev[ (pto.parada+1):length(LS_prev) ] ),
                LI_prev1 = c( rep(NA, pto.parada), LI_prev[ (pto.parada+1):length(LI_prev) ] ),
                LS_prev2 = c( rep(NA, pto.parada), pred_sarima$upper),
                LI_prev2 = c( rep(NA, pto.parada), pred_sarima$lower))







