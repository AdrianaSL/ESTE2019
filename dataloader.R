#1- Setup #########
getwd()
ifelse(!require(dplyr), install.packages('dplyr'), require(dplyr))

# 2- Dados Históricos #########
## Dados mensais do Sistema IBGE de Recuperação Automática (SIDRA)
dados <- read.csv2("ipca_201903SerieHist.csv", header = T, dec = ",", sep = ";")
dados <- dados %>% filter( ano > 1995 )
y <- dados$IPCA 

message("y é a variável reposta! \n=> 
        'Inflação Anual Efetiva medida mensalmente pelo IPCA'")
dados$Data2 <- paste0(tolower(dados$mes),"/",dados$ano)

aux <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
dados$Data1 <- as.Date(paste0(rep("01", 279),
                              c(rep(aux, 23), aux[1:3]),
                              dados$ano), '%d%m%Y')

# 3- Resoluções #########

res <- data.frame(
  data = 
c("1999-06-30",  "1999-06-30",  "1999-06-30",  "2000-06-28",
  "2001-06-28",  "2002-06-27",  "2002-06-27",  "2003-06-25",
  "2003-06-25",  "2004-06-30",  "2005-06-23",  "2006-06-29",
  "2007-06-26",  "2008-07-01",  "2009-06-30",  "2010-06-22",
  "2011-06-30",  "2012-06-28",  "2013-06-28",  "2014-06-25",
  "2015-06-25",  "2016-06-30",  "2017-06-29",  "2017-07-29"),
meta =
  c(8,6,4,3.5,3.25,4,3.75,5.5,rep(4.5, 14), 4.25, 4),
banda_i =
  c(6, 4, 2, 1.5, 1.25, 1.5, 1.5, 3, 2, rep(2.5, 11), 3, 3, 2.75, 2.5),
banda_s = 
  c(10, 8, 6, 5.5, 5.25, 6.5, 6.25, 8, 7, rep(6.5, 11), 6, 6, 5.75, 5.5),
ipca_efetivo = 
  c(8.94, 5.97, 7.67, 12.53, 9.3, 9.3, 7.6, 7.6, 5.69, 3.14, 4.46, 5.9, 4.31, 5.91, 6.5, 5.84, 5.91, 6.41, 10.67, 6.29, 2.95, 3.75, NA, NA),
ano = c(1999:2002, 2003, 2003, 2004, 2004, 2005:2020)
)
res <- res[-c(5,7),]

dados <- merge(dados, res, by.x = 'ano', by.y = 'ano', all.x = T)
dados$regime <- paste("Inf=", dados$banda_i, "Meta=", m=dados$meta, "Sup=", dados$banda_s)
rm(res,y, aux)

dados <- dados %>% filter( ano > 1998 )

