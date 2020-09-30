# ALEXANDRE GANDINI
# questao 3, lista 4, Computacional PPGEST-UFRGS 2020.

# limpa variaveis na memoria:
rm(list=ls())
par(mfrow=c(1,1)) # limpa layout

library(coda)
library(mcmcse)

# funcao proporcional aa densidade:
funcao_p=function(x)
  {
  return( (x[1] + x[2] + x[3] + x[4] + x[5] )^5 )
  }

# propostas:
proposta=function(xi){
  
  x_old = 10000*xi[1] + 1000*xi[2] + 100*xi[3] + 10*xi[4] + 1*xi[5]
  
  # passos de uma uniforme discreta [0:100]:
  x_new = x_old + sample(0:100,1)
  
  if(x_new < 1) x_new = sample(1:9999,1)
  if(x_new > 9999) x_new = sample(1:9999,1)

  # pra pegar os digitos:
  x_new = sprintf("%05d", x_new)
  res = as.numeric(unlist(strsplit(x_new, "")))
  return(res)
  
}

# aceitacao da proposta:
prob_aceita=function(x_old,x_new)
  {
  min(funcao_p(x_new)/funcao_p(x_old), 1)
  }

# funcao principal:
runmcmc=function(start,n_MCMC=5000)
  {
  MCMC_d=matrix(ncol=n_MCMC, nrow=5)     #declara a matriz que conterÃ¡ os valores dos dÃ�gitos ao longo do MCMC
  #cada coluna eh uma ieraÃ§Ã£o do MCMC
  
  p=vector()                 #feclara vetor que vai guardar resultado da funÃ§Ã£o densidade para cada observaÃ§Ã£o
  
  MCMC_d[,1]=start       #escolhe valores iniciais
  
  for(it in (1:(n_MCMC-1))){    ### inicia o algoritmo
    x_old = MCMC_d[,it]
    x_new = proposta(x_old)                ## cria uma proposta
    
    if(runif(1)<prob_aceita(x_old,x_new)){
      MCMC_d[,it+1]=x_new                         #com probabilidade prob_aceita registra o valor da proposta no MCMV 
    }else{
      MCMC_d[,it+1]=x_old                         #caso contrario repete o valor antigo 
    }
    
  }
  return(MCMC_d)
}


# roda o MCMC:
MCMC_d = runmcmc(start=c(0,0,0,0,0))


# montando o numero:
X <- MCMC_d[1,]*10000 + MCMC_d[2,]*1000 + 100*MCMC_d[3,] + 10*MCMC_d[4,] + MCMC_d[5,]
hist(X,nclass=100)                         ## Olhar a distribuição de X
hist(X,nclass=500)
mean(X,na.rm=T)
sd(X,na.rm=T)

MCMC_coda <- NULL
MCMC_coda <- rbind(MCMC_d,MCMC_d[1,]*10000 + MCMC_d[2,]*1000 + 100*MCMC_d[3,] + 10*MCMC_d[4,] + MCMC_d[5,],
                   (MCMC_d[1,]+MCMC_d[2,]+MCMC_d[3,]+MCMC_d[4,]+MCMC_d[5,])^5)

row.names(MCMC_coda) <- c("10k","1k","centena","dezena","unidade",'x','f(x)')

x <- coda::mcmc(t(MCMC_coda))   #cria elemento coda.mcmc
summary(x, na.rm = TRUE)
plot(x)

accept_rate <- 1 - coda::rejectionRate(x)
accept_rate

coda::effectiveSize(x)
coda::autocorr.plot(x)

# no exercicio da lista, foram necessarias mais iteracoes do MCMC do que no exemplo visto em aula.
