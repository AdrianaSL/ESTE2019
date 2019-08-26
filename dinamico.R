
options(scipen = 9999)

### ================================================================
### Reading packages
### ================================================================

library(readxl)
library(Metrics)
library(dlm)
library(readxl)
library(expm)
library(dplyr)
library(metRology)

### ================================================================
### Dataset 
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
### Loading scripts
### ================================================================

source("function3.R")
source("Harmonicos selecionados.R")

### ================================================================
### Modeling
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

### diagnostico do modelo =====
# avaliando k = 3, 6, 9, 12 passos a frente, em um modelo que faz 12 passos a frente

k = c(3, 6, 9, 12)

## erro de previsao
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

## erro de estimacao
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

### ================================================================
### Graphics
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


  
  
  graf_dlmXsarima(y = y,
                  y_hat1 = result$ft,
                  y_hat2 = ,
                  name_yhat1 = "dinamico",
                  name_yhat2 = "SARIMA",
                  inic = 200,
                  main = "previsao_dlmXsarima_zoom200",
                  interv = int,
                  IC = T,
                  LS_prev1 = LS_prev,
                  LI_prev1 = LI_prev,
                  LS_prev2 = NULL,
                  LI_prev2 = NULL)





