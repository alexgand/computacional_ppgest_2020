# ALEXANDRE GANDINI
# questao 2, lista 3, Computacional PPGEST-UFRGS 2020.

# limpa variaveis na memoria:
rm(list=ls())
par(mfrow=c(1,1)) # limpa layout

#####
# b) Gerar Uniformes:
#####

RE = 50

# uniforme[0,1] original:
unif = runif(RE)

# U(a - 1/2; a + 1/2):
# "a" escolhido arbitrariamente:
a = 5
X = a-(1/2) + unif
# teste, como o p-valor eh alto, nao rejeitamos a hipotese de que as duas amostras vem da mesma distribuicao:
ks.test(X, punif, a - 1/2, a + 1/2)

# U(0, b):
# "b" escolhido arbitrariamente:
b = 5
X = b * unif
# teste, como o p-valor eh alto, nao rejeitamos a hipotese de que as duas amostras vem da mesma distribuicao:
ks.test(X , punif, 0, b)

# U(a,b):
# "a" e "b" escolhidos arbitrariamente:
a = 5
b = 10
X = a + (b - a) * unif
# teste, como o p-valor eh alto, nao rejeitamos a hipotese de que as duas amostras vem da mesma distribuicao:
ks.test(X, punif, a, b)
