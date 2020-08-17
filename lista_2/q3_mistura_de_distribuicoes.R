# ALEXANDRE GANDINI
# questao 3, lista 2, Computacional PPGEST-UFRGS 2020.

# MISTURA DE DISTRIBUICOES NORMAIS MULTIVARIADAS

##########
# LETRA A)
##########

# utilizando pacote pra gerar dados simulados de misturas de normais:
install.packages("mvtnorm")
library(mvtnorm)

# limpa variaveis na memoria:
rm(list=ls())
par(mfrow=c(1,1)) # limpa layout

m = 3
n = 40
mu1 = as.vector( c(0,0,0) )
mu2 = as.vector( c(3,3,3) )
mus = c(mu1,mu2)
variancias = diag(m)

# Nao utilizado: randomiza quantitade a ser utilizada proveniente da distribuicao 1:
# quantidade_dist_1 = rbinom(1, n, p1)

# proporcao de elementos de cada distribuicao:
proporcao_dist1 = 0.7
proporcao_dist2 = 1 - proporcao_dist1
proporcoes = c(proporcao_dist1, proporcao_dist2)

# matriz de dados:
dist1 = rmvnorm(proporcao_dist1 * n, mean = mu1, sigma = diag(length(mu1)) )
dist2 = rmvnorm(n - dim(dist1)[1], mean = mu2, sigma = diag(length(mu1)) )
x = as.matrix(rbind(dist1, dist2))

# visualizando as misturas misturas:
layout(1:3)
bins = 10
hist(x[,1], col=rgb(1,0,0,.5), breaks=bins)
hist(x[,2], col=rgb(0,0,1,.5), breaks=bins)
hist(x[,3], col=rgb(0,1,0,.5), breaks=bins)

##########
# LETRA B)
##########

# funcao de log verosimilhanca:
log_verossimilhanca = function(x, mus, proporcoes){
  mu1 = mus[1:3]
  mu2 = mus[4:6]
  # usando dmvnorm para obter a densidade da normal multivariada:
  return( sum( log( proporcoes[1] * mvtnorm::dmvnorm(x, mean=mu1, sigma=variancias) + proporcoes[2] * mvtnorm::dmvnorm(x, mean=mu2, sigma=variancias) ) ) )
}

# chutes iniciais:
chute1 = c(-1,0,1)
chute2 = c(2,3,4)
chutes = c(chute1, chute2)

# teste da funcao:
log_verossimilhanca(x, chutes, proporcoes)

n.iter = 500

# Matriz teta, pra guardar a evolucao dos parametros, jah coloca o valor inicial dos chutes ali:
teta = matrix(0, nrow=n.iter, ncol=length(chutes))
teta[1,] = chutes

# Vetor f.eval para avaliar o retorno da verossimilhanca para o parametro do passo atual:
f.eval = rep(0, n.iter)
f.eval[1] = log_verossimilhanca(x, chutes, proporcoes)

# Otimizacao estocastica:
for(i in 2:n.iter){
  temp = 1/log(1+i) # temperatura

  # tamanho do passo:
  d = 0.2 * sqrt(temp)
  e = runif(length(chutes) ,-d, d)
  
  # soma passo aa posicao anterior:
  te1 = teta[(i-1),] + e
  
  variacao_da_log_verossimilhanca = log_verossimilhanca(x, te1, proporcoes) - log_verossimilhanca(x, teta[(i-1),], proporcoes)
  
  # se o passo levar aa um aumento da verossimilhanca, rho vai ser grande:
  rho = min(exp(variacao_da_log_verossimilhanca/temp), 1)
  
  # decide se vai dar o passo, dah o passo com probabilidade rho:
  u = runif(1)
  teta[i,] = (u <= rho) * te1 + (u > rho) * teta[(i-1),]
  
  f.eval[i] = log_verossimilhanca(x, teta[(i-1),], proporcoes)
}

# teta que resultou na maior verossimilhanca:
pos <- which.max(f.eval)
print(teta[pos,])

# Comparando com Simulated Annealing via funcao optim(), por padrao MINIMIZA, como queremos MAXIMIZAR a log_verossimilhanca utilizamos o parametro fnscale=-1:
print(optim(par=chutes, fn=log_verossimilhanca, x=x, control=list(fnscale = -1), proporcoes=proporcoes, method='SANN')$par)

# CONCLUSAO:
# nesse exercicio, simulated annealing eh bastante sensivel ao chute inicial
# tamanho do passo (d = 0.2 * sqrt(temp)) tambem eh importante.

##########
# LETRA C)
##########

# EM:

# chutes iniciais:
chute1 = c(-1,0,1)
chute2 = c(2,3,4)
chutes = c(chute1, chute2)

# Matriz teta, pra guardar a evolucao dos parametros, jah coloca o valor inicial dos chutes ali:
teta = matrix(0, nrow=n.iter, ncol=length(chutes))
teta[1,] = chutes

# criterio de convergencia:
cc_max = 0.0001
cc = 10
conta = 1
n.iter = 500

t0 = chutes

while(cc>cc_max & conta < n.iter){
  #etapa E
  
  mu1 = t0[1:3]
  mu2 = t0[4:6]
  
  pi = (apply(dnorm(t(x) - mu1), 2, prod)) / (apply(dnorm(t(x) - mu1), 2, prod) + apply(dnorm(t(x) - mu2), 2 ,prod) )
  
  #etapa M

  t1 = rep(0, length(chutes))  
  t1[1] = sum(pi * x[,1]) / sum(pi)
  t1[2] = sum(pi * x[,2]) / sum(pi)
  t1[3] = sum(pi * x[,3]) / sum(pi)
  t1[4] = sum((1-pi) * x[,1]) / sum(1-pi)
  t1[5] = sum((1-pi) * x[,2]) / sum(1-pi)
  t1[6] = sum((1-pi) * x[,3]) / sum(1-pi)

  conta = conta + 1
  cc = sum((t0-t1)^2)
  t0 = t1
  # guarda o teta:
  teta[conta,]=t0
}

head(teta)

teta_final = teta[conta,]
print('resultado final EM:')
print(teta_final)
print(paste('quantidade de iterações:',conta))

# classificação através dos pis da última iteração:
veio_da_distribuicao_1 = pi > 0.5
proporcao_distribuicao_1 = sum(veio_da_distribuicao_1) / length(veio_da_distribuicao_1)
print(proporcao_distribuicao_1)

##########
# LETRA D)
##########

# optim() padrao, utilizando Nelder-Mead, por padrao MINIMIZA, como queremos MAXIMIZAR a log_verossimilhanca utilizamos o parametro fnscale=-1:
optim(par=chutes, fn=log_verossimilhanca, x=x, control=list(fnscale = -1), proporcoes=proporcoes)$par

##########
# LETRA G)
##########

dados <- read.table(file="C:/Users/alega/Documents/Mestrado_stats/computacional/listas/lista_2/DadosGeneticos.txt", header=TRUE)
print(dim(dados))

###############################
# via kmeans implementado no R:
###############################
resposta = kmeans(dados, centers=2)
print('Clusters:')
print(resposta[1])
# aqui o kmeans conseguiu separar as populações européias das não européias.

###############
# agora via EM:
###############
x = dados

# chutes iniciais:
chutes = rnorm(200,0,2)
t0 = chutes

# Matriz teta, pra guardar a evolucao dos parametros, jah coloca o valor inicial dos chutes ali:
teta = matrix(0, nrow=n.iter, ncol=length(chutes))
teta[1,] = chutes

# criterio de convergencia:
cc_max = 0.0001
cc = 10
conta = 1
n.iter = 500

t0 = chutes

while((cc > cc_max) & (conta < n.iter)){
  #etapa E
  
  mu1 = t0[1:100]
  mu2 = t0[101:200]
  
  pi = (apply(dnorm(t(x) - mu1), 2, prod)) / (apply(dnorm(t(x) - mu2), 2, prod) + apply(dnorm(t(x) - mu2), 2, prod))
  
  #etapa M
  for(j in 1:100){
    t1[j] = sum(pi * x[,j]) / (sum(pi))
    t1[100+j] = sum((1 - pi) * x[,j]) / (sum(1-pi))
  } 
  
  #
  conta = conta + 1
  cc = sum((t0-t1)^2)
  t0 = t1
  teta[conta,] = t0
}

# EM tambem conseguiu separar a populacao europeia das demais:
which(pi>0.5)
which(pi<0.5)
