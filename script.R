options(scipen = 9999)

library(readxl)
library(Metrics)
library(dlm)
library(readxl)
library(expm)
library(dplyr)
setwd("C:\\Users\\adria\\Google Drive\\Estudos e Congressos\\ESTE2019")

dados = read.csv2("ipca_201903SerieHist.csv", header = T, dec = ",", sep = ";")
dados = dados %>% filter( ano > 1995 )
y <- dados$IPCA 
dados$Data2 = paste0(tolower(dados$mes),"/",dados$ano)
# dados$Data = as.Date.character(dados$Data2, "%m/%y")
# dados$Data2 = as.Date.character(dados$Data, "%m/%y")


source("function3.R")
source("Harmonicos selecionados.R")





#Definindo F,G e D

G.mer = dlm::bdiag(G1,G2[[1]],G2[[2]])
FF = matrix(c(1,0,1,0,1,0)) 
F.mer  = matrix(FF,nrow=length(FF),ncol=length(y))
D1 <- matrix(1/0.8,2,2)
D2 <- matrix(1/1,2,2)
D.mer = dlm::bdiag(D1,D2,D2)
D.mer_estatico = dlm::bdiag(D1,D1,D1)


F.t<-F.mer
G<-G.mer

# Definindo as quantidades iniciais
n<-nrow(F.t)
m0 <- rep(0,n)
C0 <- diag(100000,n,n) 
n0 <- 0.1; d0 <- .001

#n0 / d0
#(n0 / 2)/((d0/2) ^ 2)




pto.parada <- length(y) - 12 #prev a partir de pto.parada

#com desconto - dinamico
resultados <- modelo.2ordem.var.desc_prev.k(y, m0, C0, n0, d0, Wt = NULL, F.t, G, D.mer,h0=pto.parada,int=NULL)
grafico(resultados, harmonic = 2, d1 = 1/D1[1,1], d2 = 1/D2[1,1], main = "mercado_componentes_din.pdf")
graf_previsao( resultados, pto.parada, main = "mercado_prev_uni.pdf")

#sem desconto - estatico
resultados_estat <- modelo.2ordem.var.desc_prev.k(y, m0, C0, n0, d0, Wt = NULL, F.t, G, D.mer_estatico,h0=pto.parada,int=NULL)
grafico(resultados_estat, harmonic = 2, d1 = 1/D1[1,1], d2 = 1/D2[1,1], main = "mercado_componentes_estat.pdf")
graf_previsao( resultados_estat, pto.parada, main = "mercado_prev_uni_estat.pdf")




### Com intervencao ==========================================================================

## Intervencao no mes 110
int = 12*c(3:10,21) + 6

resultados2 <- modelo.2ordem.var.desc_prev.k_intervencao(y, m0, C0, n0, d0, Wt = NULL, F.t, G, D.mer,h0=pto.parada,int=int, desc = 100)
graf_previsao( resultados2, pto.parada, main = "mercado_prev_interv108.pdf", interv = 108)












grid_tendencia <- c(.8, .85, .9, .95, 1) # grid de valores para tendência
grid_saz <- c(.9, .95, 1) # grid de valores para a sazonalidade


# Tabela para comparação de preditivas
tab_pred <- matrix(NA, nrow = length(grid_saz), ncol = length(grid_tendencia))
rownames(tab_pred) <- as.character(grid_saz)
colnames(tab_pred) <- as.character(grid_tendencia)
tab_rmse <- tab_pred
tab_mape <- tab_pred
#tab_pred_norm <- tab_pred
k="Mercado Total"

for(i in 1:length(grid_saz)){
  for(j in 1:length(grid_tendencia)){
    tab_pred[i,j] <- criterios(y = y,k=k, sazonalidade = grid_saz[i], tendencia = grid_tendencia[j] )$preditiva
    tab_rmse[i,j] <- criterios(y = y,k=k,sazonalidade = grid_saz[i], tendencia = grid_tendencia[j] )$rmse
    tab_mape[i,j] <- criterios(y = y,k=k,sazonalidade = grid_saz[i], tendencia = grid_tendencia[j] )$mape
    #tab_pred_norm[i,j] <- criterios(sazonalidade = grid_saz[i], tendencia = grid_tendencia[j] )$preditiva_norm
  }
}

select_min(tab_rmse)
select_min(tab_mape)
select_max(tab_pred)

