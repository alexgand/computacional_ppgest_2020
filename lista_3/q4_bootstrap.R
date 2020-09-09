# ALEXANDRE GANDINI
# questao 4, lista 3, Computacional PPGEST-UFRGS 2020.

# limpa variaveis na memoria:
rm(list=ls())
par(mfrow=c(1,1)) # limpa layout

# amostra original:
X = c(6.2,5.1,7.6,2.5,3.5,9.4,4.1,6.3,3.0,0.8)
Y = c(6.9,5.1,7.5,11.1,10.9,4.2,10.5,6.8,12.3,14.3)

#####
# a) IC de 95% para beta0 e beta1, estimadores de minimos quadrados da regressao Y = beta0 + beta1 * X:
#####

# numero de reamostragens bootstrap:
B = 1000

beta0 = numeric(B)
beta1 = numeric(B)

for(i in 1:B)
  {
  # Amostra bootstrap, com reposicao, da amostra original:
  Xboot = sample(X, replace=T)
  Yboot = sample(Y, replace=T)
  # adiciona coluna de 1s, pra estimar os betas de forma matricial:
  X_matriz = cbind(1, Xboot)
  # transforma Y em matriz:
  Y_matriz = matrix(Yboot,ncol=1)
  # estimadores por minimos quadrados:
  betas_hat = solve(t(X_matriz) %*% X_matriz) %*% t(X_matriz) %*% Y_matriz
  # tambem poderia ser estimado via lm do R:
  # betas_hat <- coef(lm(Yboot ~ Xboot))
  beta0[i] = betas_hat[1]
  beta1[i] = betas_hat[2]

}

hist(beta0)
hist(beta1)

mean(beta0)
print(paste("Intervalo de confianças para beta0:", quantile(beta0, c(0.025,0.975) ) ) )
mean(beta1)
print(paste("Intervalo de confianças para beta1:", quantile(beta1, c(0.025,0.975) ) ) )


#####
# b) IC de 95% para beta0 e beta1, usando bootstrap parametrico com erro ~ Normal(0, sigma2):
#####

# utiliza como parametros os estimadores de minimos quadrados da amostra original e sigma2 variancia dos residuos.

# agora utilizando lm do R:
modelo = lm(Y ~ X)
betas = modelo$coefficients
sigma = sigma(modelo)

beta0boot = numeric()
beta1boot = numeric()

for (i in 1:B)
  {
  # Y gerado para o bootstrap: X amostra bootstrap dos dados, coeficientes estimados acima, mais erro normal.
  amostra_X = sample(X, replace=T)
  Y_boot = betas[1] + betas[2] * amostra_X + rnorm( length(X), mean = 0, sd = sigma)
  betas_boot <- coef(lm(Y_boot ~ X))
  beta0boot[i] <- betas_boot[1]
  beta1boot[i] <- betas_boot[2]
  }

hist(beta0boot)
hist(beta1boot)

mean(beta0)
print(paste("Intervalo de confianças para beta0:", quantile(beta0boot, c(0.025,0.975) ) ) )
mean(beta1)
print(paste("Intervalo de confianças para beta1:", quantile(beta1boot, c(0.025,0.975) ) ) )


#####
# c) Teste de permutacao para verificar beta1 = 0. Eh possivel testar beta0 = 0 via teste de permutacao?
#####

# Se beta1 = 0, entao Y = beta0 + erro, ou seja, X nao deveria ter influencia alguma sobre Y.

# estimacao de beta1 por minimos quadrados usando a amostra original:
betas = coef(lm(Y ~ X))
beta1 = betas[2]

# permuntando Y (trocando a ordem dos dados) em B replicacoes, e mantendo a ordem dos dados de X original,
# estamos arruinando a relacao de dependencia entre os pares (X,Y), e assim podemos estimar um
# p-valor para beta1:

teste = numeric(B)

for (i in 1:B)
  {
  Y_perm <- sample(Y)
  teste[i] = coef(lm(Y_perm ~ X))[2]
  }

# histograma de beta1, via teste de permutacao:
hist(teste,50)

# Estimando o  p-valor, sob H0, do nosso teste de permutacao ter produzido beta1 tao ou mais extremo quanto o beta1 estimado com os dados originais:
mean(teste <= beta1)

# P-valor baixo: em nenhuma vez o teste de permutacao gerou um beta1 tao extremo quanto o estimado com o modelo completo,
# entao rejeitamos a hipotese nula h0: beta1 = 0.

# Eh possivel testar beta0 = 0?
# Segundo a logica do teste de permutacao acima, nao eh adequado realizar o teste dessa forma,
# pois se beta0 = 0, o modelo fica Y = beta1 * X + erro, ou seja, existe influencia
# de X sobre Y.