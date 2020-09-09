# ALEXANDRE GANDINI
# questao 1, lista 3, Computacional PPGEST-UFRGS 2020.

# limpa variaveis na memoria:
rm(list=ls())
par(mfrow=c(1,1)) # limpa layout

# para testar o ajuste das distribuicoes geradas:
install.packages('EnvStats')
library(EnvStats)
install.packages('evd')
library(evd)

#####
# a) Pareto:
#####

RE = 1000
unifs = runif(RE)

alpha = 1 # parametro da distribuicao pareto escolhido arbitrariamente
beta = 1 # parametro da distribuicao pareto escolhido arbitrariamente

# gera pareto:
pareto = beta / ( (1 - unifs)^(1/alpha) )

# teste:
ks.test(pareto, "ppareto", beta, alpha)
# hipotese nula: as duas amostras vem da mesma distribuicao.
# como o p-value foi alto, nao rejeitamos a hipotese nula.

#####
# b) Gumbel Padrão:
#####

RE = 1000
unifs = runif(RE)

# gera gumbel:
gumbel = -log( -log(unifs) )

# teste:
ks.test(gumbel, "pgumbel")
# hipotese nula: as duas amostras vem da mesma distribuicao.
# como o p-value foi alto, nao rejeitamos a hipotese nula.

#####
# c) Distribuicao F:
#####

RE = 1000

# gera normais(0,1) via transformada de Box-Muller:
N1 = matrix(nrow=RE, ncol=RE)
N2 = matrix(nrow=RE, ncol=RE)

i=1
for (i in (1:RE))
{
  unif1 = runif(RE)
  unif2 = runif(RE)
  
  N1[,i] = (-2 * log(unif1))^(1/2) * cos(2 * pi * unif2)
  N2[,i] = (-2 * log(unif1))^(1/2) * sin(2 * pi * unif2)
}

# teste de uma das normais:
ks.test(N1[,1], "pnorm")
ks.test(N2[,1], "pnorm")
# normais passaram no teste (p-value alto nao rejeita hipotese nula de que vem da mesma distribuicao)

# gera qui-quadrados (soma de normais ao quadrado):
Z1 = rowSums(N1^2)
Z2 = rowSums(N2^2)

# teste das qui-quadrado:
ks.test(Z1, "pchisq", RE)
ks.test(Z2, "pchisq", RE)
# qui-quadrados passaram no teste (p-value alto nao rejeita hipotese nula de que vem da mesma distribuicao)

# agora finalmente gerando a distribuicao F:

df1 = RE # parametro da distribuicao F escolhido arbitrariamente
df2 = RE # parametro da distribuicao F escolhido arbitrariamente

F = (Z1/df1) / (Z2/df2)

# teste:
ks.test(F, "pf", df1, df2)
# hipotese nula: as duas amostras vem da mesma distribuicao.
# como o p-value foi alto, nao rejeitamos a hipotese nula.

#####
# d) Binomial negativa:
#####

RE = 1000

# funcao massa de probabilidade da binomial negativa:
pmf_negative_binomial = function(x, r, p)
{
  return( choose((r+x-1),x) * p^r * (1-p)^x )
}

# parametros da binomial negativa escolhidos arbitrariamente:
p = 0.4
r = 14

# limite superior de probabilidade para definir os intervalos:
limite_superior = 0.99

# calcula as probabilidades ateh o totalizar limite superior:
prob = c(0)
x = 1 # x inicial
while( sum(prob) < limite_superior)
  {
  prob = c(prob, pmf_negative_binomial(x,r,p))
  x = x+1
}

# agora define intervalos (probabilidade acumulada):

ints = cumsum(prob)
ints = ints[ints <= limite_superior]
ints = c(ints ,1)

# possiveis valores de x ateh atingir o limite superior:
x = 1:length(ints) 

binomial_negativa = numeric(RE)

# ve em qual intervalo a uniforme[0,1] caiu, essa eh a nossa variavel gerada de interesse:
for (i in 1:RE){
  unif <- runif(1)
  posicao = sum(unif>ints)
  binomial_negativa[i] = x[posicao]
}

binomial_negativa
hist(binomial_negativa, breaks = RE/20)

#####
# e) Discreta bivariada:
#####

RE = 1000

# funcao massa de probabilidade conjunta:
pmf_conjunta = function(variaveis, n, p)
{
  x = as.numeric(variaveis[1])
  y = as.numeric(variaveis[2])
  return( choose(n,x) * choose(n,y) * ( ( x^y * (n - x)^(n - y) * p^x * (1-p)^(n - x) ) / (n^n) ) )
}

# parametro escolhido arbitrariamente:
p = 0.4
n = 12

# todas as possiveis combinacoes de valores de x e de y ateh n:

x_range = 0:n
y_range = 0:n
vars = expand.grid(x_range,y_range)

# calcula as probabilidades:
prob = c(0)

for ( i in 1:dim(vars)[1])
{
  prob = c(prob, pmf_conjunta(vars[i,],n,p))
  i = i + 1
}

# agora define intervalos (probabilidade acumulada):

ints = cumsum(prob)
ints = c(ints ,1)

bivariada_simulada = matrix(0, ncol=2, nrow=RE)

# ve em qual intervalo a uniforme[0,1] caiu, essa eh a nossa variavel gerada de interesse:
for (i in 1:RE){
  unif <- runif(1)
  posicao = sum(unif>ints)
  bivariada_simulada[i,1] = as.numeric(vars[posicao,1])
  bivariada_simulada[i,2] = as.numeric(vars[posicao,2])
}

# resultado final:
print(bivariada_simulada)

# histogramas:
layout(1:2)
hist(bivariada_simulada[,1])
hist(bivariada_simulada[,2])
