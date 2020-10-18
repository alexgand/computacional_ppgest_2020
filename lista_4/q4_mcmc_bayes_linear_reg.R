# ALEXANDRE GANDINI
# questao 4, lista 4, Computacional PPGEST-UFRGS 2020.

# limpa variaveis na memoria:
rm(list=ls())
par(mfrow=c(1,1)) # limpa layout

library(coda)
library(mcmcse)

m0=0; m1=0; t0=1/200; t1=1/200; a=1; b=0.2

x=runif(50,1,9)
y=5-0.5*x +rnorm(50)

# a) Qual a posteriori para esse modelo?

posteriori=function(v){
  ans=-Inf
  if(v[3]>0){
    ans= sum(dnorm(y,v[1]+v[2]*x,sqrt(1/v[3]),log=T) +
               dnorm(v[1],0,sqrt(1/t0),log=T)+dnorm(v[2],0,sqrt(1/t1),log=T) + dgamma(v[3],a,b,log=T))}
  ans}

# b) metropolis hastings
proposta=function(teta_old,m=0.2){      
  p=sample(c(1,2,3),1)
  c(ifelse(1==p,  runif(1,teta_old[1]-m,teta_old[1]+m), teta_old[1]), # proposta sugerida pelo 
    ifelse(2==p,  runif(1,teta_old[2]-m,teta_old[2]+m), teta_old[2]), # exercicio
    ifelse(3==p,  runif(1,teta_old[3]-m,teta_old[3]+m), teta_old[3]))}

prob_aceita=function(teta_old,teta_new){ # probabilidade de aceitacao - como estamos trabalhando com o 
  min(exp(posteriori(teta_new)-posteriori(teta_old)), 1)} #log temos que trabalhar com a exp

n_MCMC=5000  #n�mero de itera��es do MCMC
MCMC_teta=matrix(NA,nrow=n_MCMC,ncol=3)
MCMC_teta[1,]=c(0,0,1) #valor inicial 


for(it in 1:(n_MCMC-1)){    ### inicia o algoritmo
  teta_old=MCMC_teta[it,]
  teta_new=proposta(teta_old) ## cria uma proposta
  
  if(runif(1)<prob_aceita(teta_old,teta_new)){
    MCMC_teta[it+1,]=teta_new #com probabilidade prob_aceita registra o valor da proposta no MCMV 
  }else{
    MCMC_teta[it+1,]=teta_old #caso contrario repete o valor antigo 
  }}


require("coda")

X=mcmc(MCMC_teta)
plot(X)
summary(X)
(1 - rejectionRate(X)) # taxa de aceitacao
ef1 = effectiveSize(X) # tamanho efetivo

# c) avalie m e convergencia

ESSs=matrix(NA,ncol=4,nrow=50)
valor.m=seq(0.1,5,0.1) # vou testar o m de 0.1 ate 5

for(j in 1:50){
  proposta=function(teta_old,m=valor.m[j]){    # gera uma proposta para teta com dist uniforme em um intervalo de [teta-0.2, teta+o.2]  
    p=sample(c(1,2,3),1)
    c(ifelse(1==p,  runif(1,teta_old[1]-m,teta_old[1]+m), teta_old[1]),
      ifelse(2==p,  runif(1,teta_old[2]-m,teta_old[2]+m), teta_old[2]),
      abs( ifelse(3==p,  runif(1,teta_old[3]-m,teta_old[3]+m), teta_old[3])))}
  
  prob_aceita=function(teta_old,teta_new){
    min( exp(posteriori(teta_new)-posteriori(teta_old)), 1)}
  
  n_MCMC=5000  #n�mero de itera��es do MCMC
  MCMC_teta=matrix(NA,nrow=n_MCMC,ncol=3) #declara o vetor que conter� os valores de teta 
  MCMC_teta[1,]=c(0,0,1)        #escolhe valores iniciais
  
  for(it in 1:(n_MCMC-1)){    ### inicia o algoritmo
    teta_old=MCMC_teta[it,]
    teta_new=proposta(teta_old)                ## cria uma proposta
    
    if(runif(1)<prob_aceita(teta_old,teta_new)){
      MCMC_teta[it+1,]=teta_new    #com probabilidade prob_aceita registra o valor da proposta no MCMV 
    }else{
      MCMC_teta[it+1,]=teta_old    #caso contrario repete o valor antigo 
    }}
  
  X=mcmc(MCMC_teta)
  ESSs[j,1]= mean(1 - rejectionRate(X))
  ESSs[j,2:4]=effectiveSize(X)}

plot(valor.m,ESSs[,1],main="Tx med aceit", type = "l"); abline(h=0.23) #ideal=0.23

#d) trocar o qij pra uma normal (i,s�)

# aqui s� muda a proposta que antes era uma uniforme e agora vai passar a ser uma normal 

proposta=function(teta_old,m=0.2){    # gera uma proposta para teta com dist uniforme em um intervalo de [teta-0.2, teta+o.2]  
  p=sample(c(1,2,3),1)
  c(ifelse(1==p, rnorm(1,teta_old[1],1), teta_old[1]),
    ifelse(2==p, rnorm(1,teta_old[2],1), teta_old[2]),
    ifelse(3==p, rnorm(1,teta_old[3],1), teta_old[3]))}

prob_aceita=function(teta_old,teta_new){ # probabilidade de aceitacao - como estamos trabalhando com o 
  min(exp(posteriori(teta_new)-posteriori(teta_old)), 1)} #log temos que trabalhar com a exp

n_MCMC=5000  #n�mero de itera��es do MCMC
MCMC_teta=matrix(NA,nrow=n_MCMC,ncol=3)
MCMC_teta[1,]=c(0,0,1) 


for(it in 1:(n_MCMC-1)){    ### inicia o algoritmo
  teta_old=MCMC_teta[it,]
  teta_new=proposta(teta_old) ## cria uma proposta
  
  if(runif(1)<prob_aceita(teta_old,teta_new)){
    MCMC_teta[it+1,]=teta_new #com probabilidade prob_aceita registra o valor da proposta no MCMV 
  }else{
    MCMC_teta[it+1,]=teta_old #caso contrario repete o valor antigo 
  }}

X=mcmc(MCMC_teta)
plot(X)
summary(X)
(1 - rejectionRate(X)) # taxa de aceitacao
ef2 = effectiveSize(X) # tamanho efetivo

# e) 

# f) Qual das op�c~oes de proposta �e melhor em termos do n�umero de itera�c~oes at�e obter o
# mesmo tamanho amostral efetivo.

# Analisando os dois resultados temos que o qij da uniforme � mais eficiente para var1 
# e var2 porem pra var3 o da normal � melhor
ef1; ef2

# g) Quais as estimativas de m�edia `a posteriori?
apply(MCMC_teta[-c(1:100),], 2,mean) #estimativas

# h) metodo mais efetivo - gibbs ou metropolis

# eu acho que o metodo do gibbs � mais eficiente porque a convergencia dele � mais rapida que o de 
# metropolis
