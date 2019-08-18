# Autocorrelacao e autocorrelacao parcial
fac_facp <-  function(x, titulo, lag, top){
  fac <- acf(x, plot = F, lag.max = lag)
  facp <- pacf(x, plot = F, lag.max = lag)
  h1 <- autoplot(fac) + 
    xlab('k') + 
    ylab(expression(rho(k))) + 
    ggtitle(paste("Autocorrelação",titulo))
  h2 <- autoplot(facp) + 
    xlab('k') + 
    ylab(expression(phi[kk])) + 
    ggtitle(paste("Autocorrelação Parcial",titulo))
  grid.arrange(h1, h2, ncol = 2, top=top) #facv amponta necessidade de diferencição
}
# Grafico Ljung Box
ljungbox <- function(ajuste, ordem){
  res <- ajuste$residuals/sqrt(ajuste$sigma2) #residuo padronizado
  pval = numeric(20-ordem)
  for (i in (ordem+1):20){
    estat <- Box.test(res, i,type = "Ljung-Box")$statistic
    pval[i-ordem] = pchisq(estat, i - ordem,lower.tail = FALSE)}
  ljung_box <- ggplot(data.frame(x=(ordem+1):20,y= pval),aes(x,y)) +
    geom_point()+xlab('H') + 
    ylab(expression(Q[H])) +
    ggtitle('Teste de Ljung-Box') +
    ylim(c(0,1)) +
    geom_hline(yintercept = 0.05, colour = 'blue',lty = 2)
  
  ljung_box
  }


# Avaliação de Normalidade
norm_test <- function(ajuste, titulo ){
  res = ajuste$residuals/sqrt(ajuste$sigma2)  #residuo padronizado
  #QQnorm
  dados = data.frame(res = as.numeric(res))
  qqplot=ggplot(dados, aes(sample = res)) + 
    stat_qq() +
    geom_abline(intercept = mean(res),slope =  sd(res)) +
    xlab('Quantis teóricos')+ylab('Quantis amostrais') +
    ggtitle(paste('QQnorm',titulo))
  #Teste de normalidade de Jarque-Bera
  print(DescTools::JarqueBeraTest(res))
 autoplot(qqplot)
}

# Previsão
previsao <- function(ajuste, n, alpha){
  previsao = forecast::forecast(ajuste, h = n, level = alpha)
  print(previsao$mean)
  autoplot(previsao)+
    ggtitle('Previsão')+
    xlab('t')+
    ylab(expression(x[t])) +
    # geom_line(data=IPCA) +
    scale_color_manual("Dataset",values = c("p1" = "darkgreen", "p2" = "red"))
}
