library(Metrics)
library(dlm)
library(readxl)
library(expm)

### Referencia

# http://hedibert.org/wp-content/uploads/2016/07/2016-schmidt-lopes.pdf


### Modelo Dinamico para previsao k passos
modelo.2ordem.var.desc_prev.k <- function(y, m0, C0, n0, d0, Wt = NULL, F.t, G, D,h0){
  ## consideramos a matriz de covarancia Vt desconhecida
  
  # y: variavel resposta
  # m0 : vetor de medias [ theta0 ~ N( m0, C0 ) ]
  # C0 : matriz de variancia [ theta0 ~ N( m0, C0 ) ]
  # d0, n0 : para calcular S0 (entra no filtro de kalman - descobrir como definir tais valores)
  # Wt: se Wt == NULL (default) consideramos Wt (estrutura de covariancia entre os componentes de theta_t) eh desconhecida
  # F.t : vetor de covariaveis [ se nao tem covariaveis F.t = ( 1 0 0 ... 0) ]
  # G : matrix de evolucao
  # D : matrix do inverso do desconto [ D = 1/delta ] 
  # h0 ???
  
  
  D[D==0]<-1
  
  ## Definindo os objetos
  
  n <- nrow(F.t)     # numero de linhas de F.t (equivale ao p da teoria)
  r <-1              # ???
  T.t <- length(y)   # tamanho da amostra
  S0 <- d0/n0        # ???
  
  mt <- matrix(0,nrow=n,ncol=T.t)               
  Ct <- array(rep(diag(n),T.t),dim=c(n,n,T.t))
  Rt <- array(rep(diag(n),T.t),dim=c(n,n,T.t))
  ft <- matrix(0,nrow=T.t,ncol=r)
  at <- matrix(0,nrow=n,ncol=T.t)
  Qt <- matrix(0,nrow=T.t,ncol=r)
  et <- matrix(0,nrow=T.t,ncol=r)
  At <- matrix(0,nrow=n,ncol=T.t)
  LSt <- matrix(0,nrow=T.t,ncol=r)
  LIt <- matrix(0,nrow=T.t,ncol=r)
  dt <- matrix(0,nrow=T.t,ncol=r)
  nt <- matrix(0,nrow=T.t,ncol=r)
  St <- matrix(0,nrow=T.t,ncol=r)
  
  aux=0
  
  ## Priori em t=1
  at[,1] <- G%*%m0                # evoluing the state [ at = Gt * m_{t-1} ]
  Rt[,,1] <- D*(G%*%C0%*%(t(G)))  # evoluing the state [ Rt = Gt * C_{t-1} * G't + Wt] Neste caso, como estamos supondo Wt desconhecido vamos usar o fator de desconto ] 
  # Logo, [Rt = (Gt * C_{t-1} * G't)/delta], onde delta eh o fator de desconto
  # temos, D = 1/delta
  
  ## No tempo 1
  t=1
  ft[t,] <- t(F.t[,t])%*%at[,t]                     # prediction the observation [ ft = F't * at ]
  Qt[t,] <- t(F.t[,t])%*%Rt[,,t]%*%F.t[,t] + S0  # prediction the observation [ Qt = 1 + (F't * Rt * Ft) ]
  LSt[t] <- ft[t,] + qt(0.975, n0)*sqrt(Qt[t,])     # prediction the observation [ yt|D_{t-1} ~ t_{n_{t-1}}( ft, S_{t-1}*Qt ) ]
  LIt[t] <- ft[t,] + qt(0.025, n0)*sqrt(Qt[t,])     # prediction the observation [ yt|D_{t-1} ~ t_{n_{t-1}}( ft, S_{t-1}*Qt ) ]
  
  
  At[,t] <- Rt[,,t]%*%F.t[,t]%*%solve(Qt[t,])   # NAO SEI 
  et[t,] <- y[t]-ft[t,]                         # [ yt - ft ]
  
  nt[t] <- n0+1                           
  dt[t] <- d0 + S0*(et[t,]^2)/Qt[t,]    # NAO SEI
  St[t] <- dt[t]/nt[t]                  # NAO SEI
  
  mt[,t] <- at[,t] + At[,t]*et[t,]                           # 
  Ct[,,t] <- (St[t]/S0)*(Rt[,,t]-At[,t]%*%t(At[,t])*Qt[t,])
  
  
  for(t in 2:T.t){			
    #if(t %in% int){
    #  D[1,1]<-D[1,1]/10
    #}
    
    ft[t,] <- t(F.t[,t])%*%at[,t]
    Qt[t,] <- t(F.t[,t])%*%Rt[,,t]%*%F.t[,t]+S0
    LSt[t] <- ft[t,] + qt(0.975, n0)*sqrt(Qt[t,]) 
    LIt[t] <- ft[t,] + qt(0.025, n0)*sqrt(Qt[t,]) 
    
    # Posteriori em t=1
    
    At[,t] <- Rt[,,t]%*%F.t[,t]%*%solve(Qt[t,])
    et[t,] <- y[t]-ft[t,]
    
    nt[t] <- n0+1
    dt[t] <- d0+S0*(et[t,]^2)/Qt[t,]
    St[t] <- dt[t]/nt[t]
    
    mt[,t] <- at[,t]+At[,t]*et[t,]
    Ct[,,t] <- (St[t]/S0)*(Rt[,,t]-At[,t]%*%t(At[,t])*Qt[t,])
    
    
    if(t <= h0){
      at[,t] <- G%*%mt[,t-1]
      Rt[,,t] <- D*(G%*%Ct[,,t-1]%*%(t(G)))
      
      # PrevisÄo k passo-a-frente
      ft[t,] <- t(F.t[,t])%*%at[,t]
      Qt[t,] <- t(F.t[,t])%*%Rt[,,t]%*%F.t[,t] + St[t-1]
      
      # Posteriori em t
      At[,t] <- Rt[,,t]%*%F.t[,t]%*%solve(Qt[t,])
      et[t,] <- y[t]-ft[t,]
      
      nt[t] <- nt[t-1]+1
      dt[t] <- dt[t-1]+St[t-1]*(et[t,]^2)*solve(Qt[t,])
      St[t] <- dt[t]/nt[t]
      
      mt[,t] <- at[,t]+At[,t]*et[t,]
      Ct[,,t] <- (St[t]/St[t-1])*(Rt[,,t]-At[,t]%*%t(At[,t])*Qt[t,])
      
    }
    
    else{ #evolucao a partir de h0
      
      at[,t] <- (G%^%(t-h0))%*%mt[,h0-1]   
      
      
      if(t == h0+1){
        Pt<-G%*%Ct[,,h0-1]%*%t(G)
        Wt <- D*Pt*t(D)
      }
      aux=aux+Wt
      Rt[,,t] <- Pt + aux
      
      #Rt[,,t] <- Pt + Wt
      
      
      # PrevisÄo k passo-a-frente
      ft[t,] <- t(F.t[,t])%*%at[,t]
      Qt[t,] <- t(F.t[,t])%*%Rt[,,t]%*%F.t[,t] + St[h0-1]
      
      # Posteriori em t
      At[,t] <- Rt[,,t]%*%F.t[,t]%*%solve(Qt[t,])
      et[t,] <- y[t]-ft[t,]
      
      nt[t] <- nt[h0-1]+1
      dt[t] <- dt[h0-1]+St[h0-1]*(et[t,]^2)*solve(Qt[t,])
      St[t] <- dt[h0]/nt[h0]
      
      mt[,t] <- at[,t]+At[,t]*et[t,]
      Ct[,,t] <- (St[t]/St[t-1])*(Rt[,,t]-At[,t]%*%t(At[,t])*Qt[t,])
    }  
    
    
  }
  
  result <- list(mt,Ct,ft,Qt,et,nt,Rt)
  names(result) <- c("mt", "Ct", "ft", "Qt", "et", "nt","Rt")
  return(result)
}

### Modelo Dinamico para previsao k passos com intervencao
modelo.2ordem.var.desc_prev.k_intervencao <- function(y, m0, C0, n0, d0, Wt = NULL, F.t, G, D,h0,int=NULL,desc = 10){
  D[D==0]<-1
  #Definindo os objetos
  n <- nrow(F.t) ; r <-1
  T.t <- length(y)
  S0 <- d0/n0
  mt <- matrix(0,nrow=n,ncol=T.t)
  Ct <- array(rep(diag(n),T.t),dim=c(n,n,T.t))
  Rt <- array(rep(diag(n),T.t),dim=c(n,n,T.t))
  ft <- matrix(0,nrow=T.t,ncol=r)
  at <- matrix(0,nrow=n,ncol=T.t)
  Qt <- matrix(0,nrow=T.t,ncol=r)
  et <- matrix(0,nrow=T.t,ncol=r)
  At <- matrix(0,nrow=n,ncol=T.t)
  LSt <- matrix(0,nrow=T.t,ncol=r)
  LIt <- matrix(0,nrow=T.t,ncol=r)
  dt <- matrix(0,nrow=T.t,ncol=r)
  nt <- matrix(0,nrow=T.t,ncol=r)
  St <- matrix(0,nrow=T.t,ncol=r)
  aux=0
  
  #Priori em t=1
  at[,1] <- G%*%m0
  Rt[,,1] <- D*(G%*%C0%*%(t(G)))    
  
  t=1
  ft[t,] <- t(F.t[,t])%*%at[,t]
  Qt[t,] <- t(F.t[,t])%*%Rt[,,t]%*%F.t[,t]+S0
  LSt[t] <- ft[t,] + qt(0.975, n0)*sqrt(Qt[t,]) 
  LIt[t] <- ft[t,] + qt(0.025, n0)*sqrt(Qt[t,]) 
  
  
  At[,t] <- Rt[,,t]%*%F.t[,t]%*%solve(Qt[t,])
  et[t,] <- y[t]-ft[t,]
  
  nt[t] <- n0+1
  dt[t] <- d0+S0*(et[t,]^2)/Qt[t,]
  St[t] <- dt[t]/nt[t]
  
  mt[,t] <- at[,t]+At[,t]*et[t,]
  Ct[,,t] <- (St[t]/S0)*(Rt[,,t]-At[,t]%*%t(At[,t])*Qt[t,])
  
  
  for(t in 2:T.t){
    
    if(t %in% int){
      Dant = D
      D[1:2,1:2]<-D[1:2,1:2]*desc
    }
    
    ft[t,] <- t(F.t[,t])%*%at[,t]
    Qt[t,] <- t(F.t[,t])%*%Rt[,,t]%*%F.t[,t]+S0
    LSt[t] <- ft[t,] + qt(0.975, n0)*sqrt(Qt[t,]) 
    LIt[t] <- ft[t,] + qt(0.025, n0)*sqrt(Qt[t,]) 
    
    # Posteriori em t=1
    
    At[,t] <- Rt[,,t]%*%F.t[,t]%*%solve(Qt[t,])
    et[t,] <- y[t]-ft[t,]
    
    nt[t] <- n0+1
    dt[t] <- d0+S0*(et[t,]^2)/Qt[t,]
    St[t] <- dt[t]/nt[t]
    
    mt[,t] <- at[,t]+At[,t]*et[t,]
    Ct[,,t] <- (St[t]/S0)*(Rt[,,t]-At[,t]%*%t(At[,t])*Qt[t,])
    
    
    if(t <= h0){
      at[,t] <- G%*%mt[,t-1]
      Rt[,,t] <- D*(G%*%Ct[,,t-1]%*%(t(G)))
      
      # PrevisÄo k passo-a-frente
      ft[t,] <- t(F.t[,t])%*%at[,t]
      Qt[t,] <- t(F.t[,t])%*%Rt[,,t]%*%F.t[,t] + St[t-1]
      
      # Posteriori em t
      At[,t] <- Rt[,,t]%*%F.t[,t]%*%solve(Qt[t,])
      et[t,] <- y[t]-ft[t,]
      
      nt[t] <- nt[t-1]+1
      dt[t] <- dt[t-1]+St[t-1]*(et[t,]^2)*solve(Qt[t,])
      St[t] <- dt[t]/nt[t]
      
      mt[,t] <- at[,t]+At[,t]*et[t,]
      Ct[,,t] <- (St[t]/St[t-1])*(Rt[,,t]-At[,t]%*%t(At[,t])*Qt[t,])
      
    }
    
    else{ #evolucao a partir de h0
      
      at[,t] <- (G%^%(t-h0))%*%mt[,h0-1]   
      
      
      if(t == h0+1){
        Pt<-G%*%Ct[,,h0-1]%*%t(G)
        Wt <- D*Pt*t(D)
      }
      aux=aux+Wt
      Rt[,,t] <- Pt + aux
      
      #Rt[,,t] <- Pt + Wt
      
      
      # Previsao k passo-a-frente
      ft[t,] <- t(F.t[,t])%*%at[,t]
      Qt[t,] <- t(F.t[,t])%*%Rt[,,t]%*%F.t[,t] + St[h0-1]
      
      # Posteriori em t
      At[,t] <- Rt[,,t]%*%F.t[,t]%*%solve(Qt[t,])
      et[t,] <- y[t]-ft[t,]
      
      nt[t] <- nt[h0-1]+1
      dt[t] <- dt[h0-1]+St[h0-1]*(et[t,]^2)*solve(Qt[t,])
      St[t] <- dt[h0]/nt[h0]
      
      mt[,t] <- at[,t]+At[,t]*et[t,]
      Ct[,,t] <- (St[t]/St[t-1])*(Rt[,,t]-At[,t]%*%t(At[,t])*Qt[t,])
    }  
    
    if(t %in% int){D = Dant}
    
    
  }
  
  result <- list(mt,Ct,ft,Qt,et,nt,Rt)
  names(result) <- c("mt", "Ct", "ft", "Qt", "et", "nt","Rt")
  return(result)
}

### Fun??o para criacao dos graficos =================================

grafico = function( resultados, harmonic = 1, d1, d2, main = " "){
  
  if( !require('metRology') ){ install.packages("metRology"); require('metRology')  }
  
  pdf(main,width=15,height=7)  
  
  
  if( harmonic == 2 ){
    #   par(mfrow=c(2,2), cex.lab = 1.6, cex.axis = 1.5, cex.main = 2)
    par(mfrow=c(2,2),cex.lab=1.6,cex.axis=1.6,lab=c(15,6,5),mar=c(2,10,6,2),cex.main=1.8)
    
  }
  
  if( harmonic == 1 ){
    l <- rbind(c(1,1,2,2),
               c(0,3,3,0))
    layout(l)
    par(cex.lab=1.6,cex.axis=1.6,lab=c(15,6,5),mar=c(2,10,6,2),cex.main=1.8)
    layout.show(max(l))
    
  }
  
  
  p = nrow( resultados$mt )
  inic = 2*p + 1
  
  
  
  # Estimacoo do nivel
  LS_nivel <- qt.scaled(0.975,resultados$nt[pto.parada,],resultados$mt[1,],sqrt(resultados$Ct[1,1,]))
  LI_nivel <- qt.scaled(0.025,resultados$nt[pto.parada,],resultados$mt[1,],sqrt(resultados$Ct[1,1,]))
  # LS_nivel <- qnorm(0.975,resultados$mt[1,],sqrt(resultados$Ct[1,1,]))
  # LI_nivel <- qnorm(0.025,resultados$mt[1,],sqrt(resultados$Ct[1,1,]))
  plot(as.numeric(y),pch = 20, ylab = " ", xlab = " ", main = "Estima??o do Nivel",
       xlim = c( inic, length(y) ), lwd = 5,  xaxt = "n", bty = "n", las = 1)
  # axis(1, at = c(seq(1,length(y),12),length(y)), labels = format(as.Date(dados$Data)[c(seq(1,length(y),12),length(y))], '%b-%y') )
  axis(1, at = c(seq(1,length(y),12),length(y)), labels = dados$Data[c(seq(1,length(y),12),length(y))] )
  polygon(c(c(1:length(y)), rev(c(1:length(y)))), c(LI_nivel, rev(LS_nivel)), col= adjustcolor("darkgray", alpha.f = 0.5), border = NA)
  # axis(1, at = seq(1,length(y),12), labels = 2007:2018 )
  #polygon(x = c(1:length(y), rev(1:length(y))),
  #        y = c(LI_nivel, rev(LS_nivel)), col = "blue", border = NA)
  lines(resultados$mt[1,], col = 2, lwd=3)
  lines(LS_nivel, col = "darkgray",lwd=2)
  lines(LI_nivel, col = "darkgray",lwd=2)
  
  # Fator de crescimento
  LS_tendencia <- qt.scaled(0.975,resultados$nt[pto.parada,],resultados$mt[2,],sqrt(resultados$Ct[2,2,]))
  LI_tendencia <- qt.scaled(0.025,resultados$nt[pto.parada,],resultados$mt[2,],sqrt(resultados$Ct[2,2,]))
  # LS_tendencia <- qnorm(0.975,resultados$mt[2,],sqrt(resultados$Ct[2,2,]))
  # LI_tendencia <- qnorm(0.025,resultados$mt[2,],sqrt(resultados$Ct[2,2,]))
  plot(resultados$mt[2,], type = "l", col =2, main="Estima??o do Fator de Crescimento", las = 1,
       xlim = c( inic, length(y) ),  ylab = " ", xlab = " ", xaxt = "n", bty = "n", lwd=3)
  # axis(1, at = c(seq(1,length(y),12),length(y)), labels = format(as.Date(dados$Data)[c(seq(1,length(y),12),length(y))], '%b-%y') )
  axis(1, at = c(seq(1,length(y),12),length(y)), labels = dados$Data[c(seq(1,length(y),12),length(y))] )
  polygon(c(c(1:length(y)), rev(c(1:length(y)))), c(LI_tendencia, rev(LS_tendencia)), 
          col= adjustcolor("darkgray", alpha.f = 0.5), border = NA)
  lines(LS_tendencia, col = "darkgray",lwd=2)
  lines(LI_tendencia, col = "darkgray",lwd=2)
  legend("topright",paste0("Desconto \n NT: ", round(d1,2), "\n S: ", round(d2,2)  ), bty = "n", cex = 1.3)
  
  
  # Sazonal
  LS_sazonal <- qt.scaled(0.975,resultados$nt[pto.parada,],resultados$mt[3,],sqrt(resultados$Ct[3,3,]))
  LI_sazonal <- qt.scaled(0.025,resultados$nt[pto.parada,],resultados$mt[3,],sqrt(resultados$Ct[3,3,]))
  # LS_sazonal <- qnorm(0.975,resultados$mt[3,],sqrt(resultados$Ct[3,3,]))
  # LI_sazonal <- qnorm(0.025,resultados$mt[3,],sqrt(resultados$Ct[3,3,]))
  plot(resultados$mt[3,], type = "l", col =2, main = "Estima??o da Sazonalidade - Harm?nico 1", las = 1, 
       xlim = c( inic, length(y) ),  ylab = " ", xlab = " ", xaxt = "n", bty = "n", lwd=3)
  # axis(1, at = c(seq(1,length(y),12),length(y)), labels = format(as.Date(dados$Data)[c(seq(1,length(y),12),length(y))], '%b-%y') )
  axis(1, at = c(seq(1,length(y),12),length(y)), labels = dados$Data[c(seq(1,length(y),12),length(y))] )
  polygon(c(c(1:length(y)), rev(c(1:length(y)))), c(LI_sazonal, rev(LS_sazonal)), col= adjustcolor("darkgray", alpha.f = 0.5), border = NA)
  lines(LS_sazonal, col = "darkgray",lwd=2)
  lines(LI_sazonal, col = "darkgray",lwd=2)
  
  if( harmonic == 2){
    LS_sazonal <- qt.scaled(0.975,resultados$nt[pto.parada,],resultados$mt[5,],sqrt(resultados$Ct[5,5,]))
    LI_sazonal <- qt.scaled(0.025,resultados$nt[pto.parada,],resultados$mt[5,],sqrt(resultados$Ct[5,5,]))
    # LS_sazonal <- qnorm(0.975,resultados$mt[3,],sqrt(resultados$Ct[3,3,]))
    # LI_sazonal <- qnorm(0.025,resultados$mt[3,],sqrt(resultados$Ct[3,3,]))
    plot(resultados$mt[5,], type = "l", col =2, main = "Estima??o da Sazonalidade - Harm?nico 2", las = 1, 
         xlim = c( inic, length(y) ),  ylab = " ", xlab = " ", xaxt = "n", bty = "n", lwd=3)
    # axis(1, at = c(seq(1,length(y),12),length(y)), labels = format(as.Date(dados$Data)[c(seq(1,length(y),12),length(y))], '%b-%y') )
    axis(1, at = c(seq(1,length(y),12),length(y)), labels = dados$Data[c(seq(1,length(y),12),length(y))] )
    polygon(c(c(1:length(y)), rev(c(1:length(y)))), c(LI_sazonal, rev(LS_sazonal)), col= adjustcolor("darkgray", alpha.f = 0.5), border = NA)
    lines(LS_sazonal, col = "darkgray",lwd=2)
    lines(LI_sazonal, col = "darkgray",lwd=2)
  }
  
  
  # mtext(main, side = 3, line = -21, outer = TRUE)
  dev.off()
  
}

graf_previsao = function( resultados, pto.parada, main = " ", interv = 0 ){
  
  pdf(main,width=15,height=7)  
  
  par(mfrow=c(1,1),  cex.lab=1.6,cex.axis=1.6,lab=c(15,6,5),cex.main=1.8, mar=c(2,10,6,2))
  
  p = nrow( resultados$mt )
  inic = 2*p + 1
  
  # Previsao
  LS_prev <- qt.scaled(0.975,resultados$nt[pto.parada,], resultados$ft, sqrt(resultados$Qt))
  LI_prev <- qt.scaled(0.025,resultados$nt[pto.parada,], resultados$ft, sqrt(resultados$Qt))
  plot(as.numeric(y), pch = 20, ylab = " ", xlab = " ", las = 1,
       xlim = c( inic, length(y)+10 ), ylim = c(min(y) - 5, max(y) + 5),
       lwd = 2 ,  bty = "n", cex.axis = 0.1)
  axis(1, at = c( seq(1,length(y),12) + 1, length(y) ),
       labels = dados$Data2[c( seq(1,length(y), 12 ) + 1, length(y) )]
       )
  axis(2, at = seq( 0, 20, 5 ), labels = seq( 0, 20, 5 ), las = 1 )
  polygon( c( c(1:length(y)), rev( c(1:length(y)) ) ),
           c( LI_prev, rev(LS_prev) ),
           col= adjustcolor("darkgray", alpha.f = 0.5),
           border = NA)
  # abline( v = c( seq( 1,length(y), 12),length(y) ),
   #       col = "grey" , lwd = .5 )
  abline( v = pto.parada, lty = 2, lwd = 3 )
  abline( v = interv , lwd = 1, lty = 2, col = "blue")
  lines( resultados$ft, col = 2, lwd = 3 ) 
  lines( LS_prev, col = "darkgray", lwd = 2 )
  lines( LI_prev, col = "darkgray", lwd = 2 )
  
  dev.off()
  
}


### Fun??o para tabelas com medidas comparativas ======================

criterios <-  function(y,k, sazonalidade, tendencia){ #fun??o que calcula preditiva dado diferentes valores de desconto
  # D
  
  D1 <- matrix(1/tendencia,2,2)
  D2 <- matrix(1/sazonalidade,2,2)
  D <- dlm::bdiag(D1, D2)
  for (i in 2:(h-1)) {
    D <- dlm::bdiag(D, D2)
  } 
  D <- dlm::bdiag(D, 1/sazonalidade)
  
  if(k=="Residencial"||k=="Mercado Total"){
    D = dlm::bdiag(D1,D2)
  }else{
    D = dlm::bdiag(D1,D2,D2)
  }
  resultados <- modelo.2ordem.var.desc_prev.k(y, m0, C0, n0, d0, Wt = NULL, F.t, G, D,h0=pto.parada,int=NULL)
  
  # Log - Verossimilhança preditiva 
  p<-nrow(resultados$mt)
  v <- (resultados$nt[-c(1:(2*p+1))] - 1)
  
  preditiva <- sum(log((v + ((y[-c(1:(2*p+1))]-resultados$ft[-c(1:(2*p+1))])^2)/resultados$Qt[-c(1:(2*p+1))])^(-(v+1)/2)))
  # preditiva <- sum(log(sqrt(resultados$Qt[-c(1:(2*p+1))])*(dt(y[-c(1:(2*p+1))],df = resultados$nt[-c(1:(2*p+1))]-1))
  #                      + resultados$ft[-c(1:(2*p+1))]))
  # 
  preditiva_norm <- sum(dnorm(resultados$ft[-c(1:(2*p+1))], resultados$Qt[-c(1:(2*p+1))], log = T))
  #preditiva <- sum(dt.scaled(y[-c(1:(2*p+1))],resultados$nt[-c(1:(2*p+1))]-1,
  #resultados$ft[-c(1:(2*p+1))], sqrt(resultados$Qt[-c(1:(2*p+1))]), log = TRUE))
  
  actual = y[-c(1:(2*p+1))]
  predicted = resultados$ft[-c(1:(2*p+1))]
  
  # rrmse: relative root mean square error
  rmse <- rmse(actual = y[-c(1:(2*p+1))], predicted = resultados$ft[-c(1:(2*p+1))])
  rmse<-(rmse/mean(actual))*100
  # mape (Mean Absolute Percent Error)
  #mape <- mape(actual = y[-c(1:(2*p+1))], predicted = resultados$ft[-c(1:(2*p+1))])*100
  
  mape<-mean(abs((actual-predicted)/actual))*100
  
  return(list(preditiva = preditiva, preditiva_norm = preditiva_norm, rmse = rmse, mape = mape))
}


### Fun??es para destacar min e max de tabela =========================
select_min<-function(tab){
  par(mfrow=c(1,1))
  tab<-round(tab,2)
  x = 1:ncol(tab)
  y = 1:nrow(tab)
  centers <- expand.grid(y,x)
  
  image(x,y,t(tab), col = "white",
        xaxt = 'n', 
        yaxt = 'n', 
        xlab = 'tend?ncia', 
        ylab = 'sazonalidade')
  
  
  color.picker <- function(x){
    if(x == min(tab)){return("red")}
    else {return("black")}
    
  }
  
  text.cols <- sapply(tab, color.picker)
  text(centers[,2],centers[,1], c(tab),  col= text.cols)
  abline(h=y + 0.5)
  abline(v=x + 0.5)
  
  mtext(attributes(tab)$dimnames[[2]], at=1:ncol(tab), padj = -1)
  mtext(attributes(tab)$dimnames[[1]], at=1:nrow(tab), side = 2, las = 1, adj = 1.2)
  
}

select_max<-function(tab){
  par(mfrow=c(1,1))
  tab<-round(tab,2)
  x = 1:ncol(tab)
  y = 1:nrow(tab)
  centers <- expand.grid(y,x)
  
  image(x,y,t(tab), col = "white",
        xaxt = 'n', 
        yaxt = 'n', 
        xlab = 'tend?ncia', 
        ylab = 'sazonalidade')
  
  
  color.picker <- function(x){
    if(x == max(tab)){return("red")}
    else {return("black")}
    
  }
  
  text.cols <- sapply(tab, color.picker)
  text(centers[,2],centers[,1], c(tab),  col= text.cols)
  abline(h=y + 0.5)
  abline(v=x + 0.5)
  
  mtext(attributes(tab)$dimnames[[2]], at=1:ncol(tab), padj = -1)
  mtext(attributes(tab)$dimnames[[1]], at=1:nrow(tab), side = 2, las = 1, adj = 1.2)
  
}


### Calcular o predictive power =======================================

QPS =  function( pred, real ){
 return( 2 * sum( (real - pred)^2 ) / length( pred ) )
}

MAPE = function(pred, real){
  
   out = abs(real - pred)/real
   final = mean(out)
  
  return(final)
}
