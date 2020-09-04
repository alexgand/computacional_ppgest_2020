# ALEXANDRE GANDINI
# questao 3, lista 3, Computacional PPGEST-UFRGS 2020.

# limpa variaveis na memoria:
rm(list=ls())
par(mfrow=c(1,1)) # limpa layout

#####
# a) Estimar E(X) e calcular erro padrao da estimativa:
#####

RE = 10000

gera_menor_n = function(RE)
{
  result = numeric()
  for (i in 1:RE)
  {
    conta = 0
    soma = 0
    while (soma < 1)
    {soma = soma + runif(n=1)
    conta = conta + 1
    }
    result[i] = conta
  }
  return(result)
}

menor_n = gera_menor_n(RE=RE)

print(paste('E(X) estimado:', mean(menor_n) ))
print(paste('Erro padrao estimado:', sqrt( var(menor_n) / RE )  ))


#####
# b) Estimar P(X >= 10) e calcular erro padrao da estimativa:
#####

shape1 = 1.6
shape2 = 12

rbeta(n=10,shape1=shape1, shape2=shape2)

betas = matrix(nrow=RE, ncol=9)
for (i in 1:9)
{
  betas[,i] = rbeta(n=RE, shape1=shape1, shape2=shape2 )
}

fx = 1

hx = rowSums(betas) <= 1

gx = apply(dbeta(betas, shape1=shape1, shape2=shape2), 1, prod)

estimativa = hx * fx / gx

print(paste('P(X >= 10) estimado:', mean(estimativa) ))
print(paste('Erro padrao estimado:', sqrt( var(estimativa) / RE )  ))


#####
# c) Estimar P(X = 10) e calcular erro padrao da estimativa:
#####

shape1 = 1.6
shape2 = 12

rbeta(n=10, shape1=shape1, shape2=shape2)

betas = matrix(nrow=RE, ncol=10)
for (i in 1:10)
{
  betas[,i] = rbeta(n=RE, shape1=shape1, shape2=shape2 )
}

fx = 1

hx = (rowSums(betas[,1:9] < 1) & (rowSums(betas) > 1) )

gx = apply(dbeta(betas, shape1=shape1, shape2=shape2), 1, prod)

estimativa = hx * fx / gx

print(paste('P(X = 10) estimado:', mean(estimativa) ))
print(paste('Erro padrao estimado:', sqrt( var(estimativa) / RE )  ))
