# ALEXANDRE GANDINI
# questao 4, lista 2, Computacional PPGEST-UFRGS 2020.
######################################################

# POISSON SPORTS RANKING

# Dados: matriz X, onde Xij = gols/pontos do time i no time j
# X eh modelada como uma Poisson(exp{oi - dj}), onde:
# oi = potencial ofensivo do time o
# dj = potencial defensivo do time j

##########################
# Inicializacao dos dados:
##########################

# limpa variaveis na memoria:
rm(list=ls())
par(mfrow=c(1,1)) # limpa layout

# utilizando o mesmo conjunto de dados do exemplo da aula: a entrada ij indica o numero de gols que i marcou em j
X=matrix(c(0,2,5,6,3,1,3,3,1,2,
           1,0,2,1,0,0,1,1,6,0,
           2,1,0,1,0,0,2,0,1,5,
           0,2,2,0,1,3,4,1,1,1,
           2,2,3,1,0,4,4,2,3,1,
           2,1,0,1,0,0,1,0,1,3,
           1,0,2,1,0,0,0,1,2,5,
           1,0,2,1,0,2,1,0,0,0,
           0,2,2,2,0,1,3,3,0,2,
           0,0,2,1,0,1,1,0,0,0),
         byrow=T, ncol=10, nrow=10)

# numero de times:
n = 10

# inicializa vetores para os parametros:
o = numeric(n)  
p = numeric(n)
t = c(o,p)

# chute inicial:
# t0 eh um vetor onde as n primeiras entradas sao os potenciais ofensivos, oi, e as n ultimas entradas sao os potenciais defensivos, di.
t0 = c(rep(1,n),rep(0.35,n))

######################
# log verossimilhanca:
######################

# menos log verossimilhanca, pois queremos maximiza-la:
logL <- function(t){
  temp = 0
  for(i in 1:n){
    for(j in 1:n){
      temp = temp - exp(t[i]-t[n+j]) + X[i,j] * (t[i]-t[n+j])
    }
  }
  temp = -temp
  return (temp)
}

######################################################
# otimizacao utilizando optim(), por default optim() minimiza. Estamos minimizando a -logL, portanto estamos maximizando a logL:
######################################################

resultado_optim <- optim(t0,logL)$par
o_optim = resultado_optim[1:n]
d_optim = resultado_optim[(n+1) : (2*n)]

# guarda resultados para comparacao:
resultado_final = data.frame(resultado_optim)

#####################################
# Otimizacao em blocos vista em aula:
#####################################

## Definir valores iniciais para O e D
d=rep(1,10)
o=rep(1,10)


d_old=rep(0,10)
o_old=rep(0,10)         


O=matrix(ncol=10,nrow=n)  # vamos montar matrizes para guardar o historico da otimizacao
D=matrix(ncol=10,nrow=n)
# salva valores iniciais
O[,1]=o
D[,1]=d


# Define criterio de convergencia
epslon=0.1
i=1  # contador de iteracoes


while(sum((d-d_old)^2)+sum((o-o_old)^2)>epslon){   # criterio de convergencia
  o_old=o
  o=log(rowSums(X)/sum(exp(-d)))   # registra o valor antigo, e atualiza todos os o`s
  d_old=d
  d=-log(colSums(X)/sum(exp(o)))  # registra o valor antigo, e atualiza todos os d`s
  
  #Salva os valores na matriz  
  i=i+1
  O[,i]=o
  D[,i]=d
}

# olhar os resultados:
print(O)
print(D)

# guarda resultados para comparacao:
resultado_final['em_blocos'] = c(o,d)

######################################
# agora pelo METODO DE NEWTON-RAPHSON:
######################################

# Essa implementacao estah com algum problema, soh funciona com criterio de convergencia muito alto, senao nao inverte o determinante.

gradiente_da_logL <- function(x, o, d){
  resultado = matrix(0, ncol=n, nrow=2)
  # derivada em relacao aa o:
  resultado[1,] = sum(exp(o-d)) + colSums(-x)
  # derivada em relacao aa d:
  resultado[2,] = - sum(exp(o-d)) + rowSums(x)
  return(resultado)
}

# derivada segunda em relacao aa oo:
derivada_segunda_da_logL_oo <- function(x, o, d){
  resultado = matrix(0, ncol=n, nrow=n)
  diag(resultado) = sum(exp(o-d))
  return(resultado)
}

# derivada segunda em relacao aa dd:
derivada_segunda_da_logL_dd <- function(x, o, d){
  resultado = matrix(0, ncol=n, nrow=n)
  diag(resultado) = -exp(d) * sum(exp(o))
  return(resultado)
}

# a derivadas segunda em relacao aa od eh igual aa derivada segunda em relacao aa do:
derivada_segunda_da_logL_od_do <- function(x, o, d){
  resultado = matrix( - exp(o-d) , ncol=n, nrow=n)
  return(resultado)
}

# maximo criterio de convergencia:
cc_max = 3
cc = 10

otimizador_newton_rapson <- function(cc_max, cc, t0, x)
{
  # separando os os e ds:
  o = t0[1:n]
  d = t0[(n+1) : (2*n)]
  
  conta = 0
  
  while((cc > cc_max) & (conta < 400))
  {
    
    # agora a hessiana eh uma matriz 4x4 em que cada elemento eh uma outra matriz
    # procedimento para inverter matriz 4x4: multiplica pela inversa do determinante da matriz original a matriz original com elementos a e d trocados e b e c com sinal trocado:
    
    determinante = derivada_segunda_da_logL_oo(X, o, d) %*% derivada_segunda_da_logL_dd(X, o, d) - derivada_segunda_da_logL_od_do(X, o, d) %*% derivada_segunda_da_logL_od_do(X, o, d)
    
    inversa_do_determinante = solve(determinante)

    # invertendo a hessiana (multiplica pela inversa do determinante da matriz original a matriz original com elementos a e d trocados e b e c com sinal trocado):
    hessiana_inv_11 = derivada_segunda_da_logL_dd(X, o, d) %*% inversa_do_determinante
    hessiana_inv_22 = derivada_segunda_da_logL_oo(X, o, d) %*% inversa_do_determinante
    # componente [1,2] e [2,1] sao iguais:
    hessiana_inv_12 = hessiana_inv_21 = - derivada_segunda_da_logL_od_do(X, o, d) %*% inversa_do_determinante

    o_atualizado = o - (hessiana_inv_11 %*% gradiente_da_logL(X,o,d)[1,] + hessiana_inv_12 %*% gradiente_da_logL(X,o,d)[2,] + hessiana_inv_21 %*% gradiente_da_logL(X,o,d)[1,] + hessiana_inv_22 %*% gradiente_da_logL(X,o,d)[2,])
    
    d_atualizado = d - (hessiana_inv_11 %*% gradiente_da_logL(X,o,d)[1,] + hessiana_inv_12 %*% gradiente_da_logL(X,o,d)[2,] + hessiana_inv_21 %*% gradiente_da_logL(X,o,d)[1,] + hessiana_inv_22 %*% gradiente_da_logL(X,o,d)[2,])

    # controle do loop:
    cc = ( sum((d_atualizado-d)^2) + sum((o_atualizado-o)^2) )
    
    o = o_atualizado
    d = d_atualizado
    
    conta = conta + 1
    print(paste('iteracao:', conta))
    print('ofense:')
    print(o)
    print('defense:')
    print(d)
  }
  
  return( list(o,d,conta) )
}

resposta = otimizador_newton_rapson(cc_max, cc, t0, X)
print(paste('quantidade de iteracoes ateh convergencia:',resposta[3]))
print('Parametros obtidos com NEWTON-RAPSON:')
print('ofense:')
print(resposta[1])
print('defense:')
print(resposta[2])

# guarda resultados para comparacao:
resultado_final['newton'] = c(as.vector(resposta[[1]]),as.vector(resposta[[2]]))

#########################################
# agora NEWTON-RAPHSON implementado no R:
#########################################

resultado_final['nlm'] = nlm(logL, p=t0)$estimate

####################
# resultados finais:
####################

print('Comparacao entre os parametros o e de para todos os metodos:')
print(resultado_final)
