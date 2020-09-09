# ALEXANDRE GANDINI
# questao 5, lista 3, Computacional PPGEST-UFRGS 2020.

# limpa variaveis na memoria:
rm(list=ls())
par(mfrow=c(1,1)) # limpa layout

#####
# 1) Simule valores Y para a cadeia de Markov oculta:
#####

# Temos 4 estados, cada um com probabilidade 1/4:
pis = rep(1/4, 4)

# as medias sao zero para cada estado:
mus = rep(0, 4)

# desvios padrao de Y escolhidos arbitrariamente para cada estado:
sigmas = c(0.1, 0.5, 1.0, 2.0)

# matriz de transicao, probabilidades de ir do estado i para o estado j escolhidas arbitrariamente:
P = matrix(nrow=length(pis), ncol=length(pis))
P[1,] = c(0.7, 0.3, 0.0, 0.0)
P[2,] = c(0.25, 0.25, 0.25, 0.25)
P[3,] = c(0.3, 0.3, 0.3, 0.1)
P[4,] = c(0.0, 0.1, 0.4, 0.5)

# montando um gerador de normais(0, sigma2) a partir de uniformes, via Box-Muller:

gera_normais = function(n, desvio_padrao)
  {
  # n eh a quantidade de amostras da normal(0,1) que queremos gerar:

  # gera normais(0,1) via transformada de Box-Muller:
  unif1 = runif(n)
  unif2 = runif(n)
  
  N1 = (-2 * log(unif1))^(1/2) * cos(2 * pi * unif2)
  
  N1 = N1 * desvio_padrao

  return(N1)
  }

# testando se o gerador de normais(0,sigma2) funciona:
normal = gera_normais(n=100, desvio_padrao=5)
hist(normal)
# Como o p-valor do teste de KS foi alto, nao rejeitamos h0: as duas amostras vem da mesma distribuicao.
# Ou seja, o gerador de normais a partir de uniformes(0,1) estah funcionando OK.
ks.test(normal, "pnorm", sd=5)

# Agora finalmente simulando os dados Y do problema:

# inicializa variaveis:
estados = numeric() # estados ocultos
estados[1] = 1
Y = numeric() # log-retornos
Y[1] = mus[estados[1]]

tail(1:10,1)

qtd_dias = 100

for (i in 1:qtd_dias)
  {
  estado_atual = tail(estados,1)

  # escolhe proximo estado com base nas probabilidades da matriz de transicao P:
  proximo_estado = sample(1:4, size=1, replace=FALSE, prob = P[estado_atual,])
  estados = c(estados, proximo_estado)

  # gera Y com base no desvio padrao do estado em questao:
  Y[i] = gera_normais(n=1, desvio_padrao = sigmas[proximo_estado])
}

# Visualizando os pontos Y gerados:
layout(1:2)
hist(Y, breaks=20)
plot(Y)

#####
# 2)  Utilize a funcao indicada no script de aula para estimar os parametros do modelo.
#####

mod <- depmix(Y ~ 1, family = gaussian(), nstates=length(pis), ntimes=length(Y))
fm2 <- fit(mod, verbose=FALSE)
summary(fm2)

# comparando a matriz de transicao estimada pelo depmix com a matriz P usada na geracao dos dados:
print(P)
# parece que o modelo nao conseguiu identificar os verdadeiros parametros, ao menos
# para essa quantidade de dias (100).
