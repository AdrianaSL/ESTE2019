
p=12
mod2 <- dlmModTrig(s = p, dV = 5.1420, dW = 0)
w  <- 2*pi/p
h = p/2
h.aux = seq(1,h)
FF = matrix(c(1,0,mod2$FF)) 
F  = matrix(FF,nrow=length(FF),ncol=length(y))
# Definindo a matriz G
G1 <- matrix(c(1,0,1,1),2,2)
G2 = list(0)
for (i in 1:(h-1)) {
  G2[[i]] <- matrix(c(cos(w*h.aux[i]),-sin(w*h.aux[i]),
                      sin(w*h.aux[i]),cos(w*(h.aux[i]))),2,2)
}
G <- dlm::bdiag(G1, G2[[1]])
for (i in 2:(h-1)) {
  G <- dlm::bdiag(G, G2[[i]])
}
G = dlm::bdiag(G,-1)
# D
D1 <- matrix(1/1,2,2)
D2 <- matrix(1/0.95,2,2)
D <- dlm::bdiag(D1, D2)
for (i in 2:(h-1)) {
  D <- dlm::bdiag(D, D2)
} 
D <- dlm::bdiag(D, 1/0.95)

## Residencial - Vamos usar harmonico 2 baseado no teste F
G.res = dlm::bdiag(G1,G2[[2]])
FF = matrix(c(1,0,1,0)) 
F.res  = matrix(FF,nrow=length(FF),ncol=length(y))
D.final = dlm::bdiag(D1,D2)


## iNDUSTRIAL - OPCAO1 - Vamos usar harmonico 3 baseado no teste F
G.ind = dlm::bdiag(G1,G2[[1]],G2[[2]],-1)
FF = matrix(c(1,0,1,0,1,0,1)) 
F.ind  = matrix(FF,nrow=length(FF),ncol=length(y))
D.ind = dlm::bdiag(D1,D2,D2,1/0.95)

## iNDUSTRIAL - OPCAO2 - Vamos usar 2 harmonicoS baseado no teste F
G.ind = dlm::bdiag(G1,G2[[1]],G2[[2]])
FF = matrix(c(1,0,1,0,1,0)) 
F.ind  = matrix(FF,nrow=length(FF),ncol=length(y))
D.ind = dlm::bdiag(D1,D2,D2)

## iNDUSTRIAL - OPCAO3 - Vamos usar 2 harmonicoS baseado no teste F
G.ind = dlm::bdiag(G1,G2[[1]])
FF = matrix(c(1,0,1,0)) 
F.ind  = matrix(FF,nrow=length(FF),ncol=length(y))
D.ind = dlm::bdiag(D1,D2)


## cOMERCIAL - Vamos usar harmonico 2 baseado no teste F
G.com = dlm::bdiag(G1,G2[[1]],G2[[2]])
FF = matrix(c(1,0,1,0,1,0)) 
F.com  = matrix(FF,nrow=length(FF),ncol=length(y))
D.com = dlm::bdiag(D1,D2,D2)


## mercado total - Vamos usar harmonico 2 baseado no teste F
G.mer = dlm::bdiag(G1,G2[[2]])
FF = matrix(c(1,0,1,0)) 
F.mer  = matrix(FF,nrow=length(FF),ncol=length(y))
D.mer = dlm::bdiag(D1,D2)



## mercado total - comercial - residencial - industrial - Vamos usar harmonico 2 baseado no teste F
G.mer2 = dlm::bdiag(G1,G2[[1]],G2[[2]])
FF = matrix(c(1,0,1,0,1,0)) 
F.mer2  = matrix(FF,nrow=length(FF),ncol=length(y))
D.mer2 = dlm::bdiag(D1,D2,D2)
