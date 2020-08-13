# ALEXANDRE GANDINI
# questao 2, lista 2, Computacional PPGEST-UFRGS 2020.

##################
# NEWTON-RAPHSON NORMAL - APLICACAO REGRESSAO NAO LINEAR:
##################

# gera dados simulados:
set.seed(0)
n = 200
variancia = 1
erro = rnorm(n=n, mean=0, sd=sqrt(variancia))
x = runif(n=n, min=2, max=40)
beta0 = 60
beta1 = -0.05
y = beta0 * exp(beta1 * x) + erro

# visualizando Y:
plot(x,y)

# quero encontrar b0 e b1 que minimizam a soma dos erros ao quadrado:

# chute inicial:
x0 = rbind(50,-0.06)

# máximo critério de convergência:
cc_max = 0.0001
cc = 10

# gradiente da f (derivadas em relacao aa b0 e aa b1)
gradiente_da_f <- function(x0,x,y)
{
  b0 = x0[1]
  b1 = x0[2]
  grad = rbind(0,0)
  
  mu = b0 * exp(b1*x)
  
  # derivadas primeiras em relação à b0 e à b1
  grad[1] =  - sum( ( y - mu ) * exp(b1 * x) )
  grad[2] =  - sum( ( y - mu ) * mu * x) 

  return(grad)
}


# hessiana, derivadas segundas da f
hessiana <- function(x0,x,y)
{
  b0 = x0[1]
  b1 = x0[2]
  hess = matrix(0, nrow=length(x0), ncol=length(x0) )
  
  hess[1,1] = sum( exp(2*b1*x))
  
  hess[2,2] =  - sum( (y*x^2*b0*exp(b1*x) ) - (2*b0^2*x^2*exp(2*b1*x) ))
  hess[1,2] = hess[2,1] = - sum( (y*x*exp(b1*x)) - (2 * exp(2*b1*x) * x * b0 ) )
  
  return(hess)
}

otimizador_newton_rapson <- function(cc_max, cc, x0, x, y)
{
  conta = 0

  while((cc > cc_max) & (conta < 200))
  {

    x1 = x0 - solve( hessiana(x0,x,y) ) %*% gradiente_da_f(x0,x,y)
    
    # controle do loop:
    cc = sum((x1-x0)^2)
    x0 = x1
    conta = conta + 1
    print(paste('iteracao:', conta))
    print(paste('x0:', x0))
  }
  
  return(x0)
}

resposta = otimizador_newton_rapson(cc_max, cc, x0, x, y)
print('Parametros obtidos com NEWTON-RAPSON:')
print(resposta)

# conclusao: NEWTON RAPSON AQUI só funciona com chute inicial PRÓXIMO ao ponto verdadeiro.


################
# Agora com GAUSS - NEWTON:
################

# gera dados simulados:
set.seed(0)
n = 200
variancia = 1
erro = rnorm(n=n, mean=0, sd=sqrt(variancia))
x = runif(n=n, min=2, max=40)
beta0 = 60
beta1 = -0.05
y = beta0 * exp(beta1 * x) + erro

# quero encontrar b0 e b1 que minimizam a soma dos erros ao quadrado:

# chute inicial:
x0 = rbind(50,-0.07)

# máximo critério de convergência:
cc_max = 0.001
cc = 10

gradiente_da_f <- function(x0,x,y)
{
  b0 = x0[1]
  b1 = x0[2]
  grad = rbind(0,0)
  
  mu = b0 * exp(b1*x)

  # derivadas primeiras em relação à b0 e à b1
  grad[1] =  - sum( y - mu ) * sum(exp(b1 * x))
  grad[2] =  - sum( y - mu ) * sum(mu * x)
  return(grad)
}


hessiana_aproximada <- function(x0,x,y)
{
  b0 = x0[1]
  b1 = x0[2]
  
  grad_mu = c(0,0)
  
  # derivadas primeiras em relação à b0 e à b1
  grad_mu[1] = sum(exp(b1 * x))
  grad_mu[2] = sum(b0 * exp(b1*x) * x)
  
  # aproximacao da hessiana:
  hess = grad_mu %*% t(grad_mu)
  
  # gambiarra pra forcar os componentes [1,2] e [2,1] serem exatamente iguais pra nao dar erro na inversao da matriz:
  elemento_igual = round(hess[1,2])
  hess[1,2] = hess[2,1] = elemento_igual
  
  return(hess)
}


otimizador_gauss_newton <- function(cc_max, cc, x0, x, y)
{
  conta = 0
  
  while((cc > cc_max) & (conta < 20))
  {
    
    # x1 = x0 - solve( hessiana_aproximada(x0,x,y) ) %*% gradiente(x0,x,y)
    hess_aprox = hessiana_aproximada(x0,x,y)

    x1 = x0 - solve( hess_aprox ) %*% gradiente_da_f(x0,x,y)

    # controle do loop:
    cc = sum((x1-x0)^2)
    x0 = x1
    conta = conta + 1
    print(paste('iteracao:', conta))
    print(paste('x0:', x0))
  }
  
  return(x0)
}

resposta = otimizador_gauss_newton(cc_max, cc, x0, x, y)
print('Parametros obtidos com GAUSS-NEWTON:')
print(resposta)
