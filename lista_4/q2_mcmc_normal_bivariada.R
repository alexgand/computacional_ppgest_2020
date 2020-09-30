# ALEXANDRE GANDINI
# questao 2, lista 4, Computacional PPGEST-UFRGS 2020.

# limpa variaveis na memoria:
rm(list=ls())
par(mfrow=c(1,1)) # limpa layout

# proporcionalidade (do enunciado):
funcao_p = function(x)
  {
  return( exp(-(x[1]^2/2)-(x[2]^2/2)) )
  }

# Distribuicao de propostas (DISCRETAS): de -2 a 2:
proposta = function(x)
  {
  x1 = x[1] + sample(-2:2,1)
  x2 = x[2] + sample(-2:2,1)
  return(c(x1,x2))
}

# valores iniciais:
chute_inicial = c(1,1)

# probabilidade de aceitacao:
prob_aceita = function(x_new, x_old)
  {
  return( min(funcao_p(x_new)/funcao_p(x_old), 1) )
  }

# funcao principal para rodar o mcmc:
MCMC_run = function(chute_inicial, n_MCMC)
  {
  
  MCMC_d = matrix(ncol=n_MCMC, nrow=2)
  
  # o primeiro valor vai ser o chute inicial:
  MCMC_d[,1] = chute_inicial
  
  it=1
  while(it < (n_MCMC-1)){
    
    x_old = MCMC_d[,it]
    x_new = proposta(x_old) 
    
    if (runif(1) < prob_aceita(x_new,x_old))
      {
      MCMC_d[,it+1] = x_new 
      }
    else
    {
      MCMC_d[,it+1] = x_old 
    }
    
    it = it+1
  }
  return(MCMC_d)
}

# roda o MCMC e gera amostras da distribuicao:
n_MCMC = 5000
cadeia = MCMC_run(chute_inicial=chute_inicial, n_MCMC=n_MCMC)

print('media:')
rowMeans(cadeia, na.rm = TRUE)
print('variancia:')
apply(cadeia, 1, var, na.rm = TRUE)

# resultados:
require(coda)
rownames(cadeia) = c('x1', 'x2')
result = mcmc(na.omit(t(cadeia)))

summary(result)
plot(result)

print('tamanho efetivo da amostra:')
effectiveSize(result)

# marginais:

totals_x1 = table(cadeia[1,])
prob_x1 = table(cadeia[1,])/sum(totals_x1)

totals_x2 = table(cadeia[2,])
prob_x2 = table(cadeia[2,])/sum(totals_x2)

barplot(prob_x1,main="Marginal x1",ylab="Frequency",xlab="X1")
barplot(prob_x2,main="Marginal x2",ylab="Frequency",xlab="X2")
