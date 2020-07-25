##################
# OTIMIZADOR NEWTON RAPSON PARA FUNCOES DE APENAS DUAS VARIÁVEIS (x,y), COM RESTRIÇÕES.
##################

# funcao (exercício 5 da lista 1):
f0 = function(x,y,lambda) x^4 + y^2 + 4*x*y + 
f1x = function(x,y,lambda) 4*x^3 + 4*y -2*lambda*x
f1y = function(x,y,lambda) 2*y + 4*x - 2*lambda*y
f1lambda = function(x,y,lambda) -x^2 - y^2 + 1
f2xx = function(x,y,lambda) 12*x^2 - 2*lambda
f2yy = function(x,y,lambda) 2*y + 4*x - 2*lambda*y
f2lambdalambda = function(x,y,lambda) 0
f2xy = function(x,y,lambda) 4
f2yx = function(x,y,lambda) 4
f2xlambda = function(x,y,lambda) -2*x
f2ylambda = function(x,y,lambda) -2*y
f2lambdax = function(x,y,lambda) -2*x
f2lambday = function(x,y,lambda) -2*y

# chute inicial (x,y,lambda):
x0 = c(2,2,2)
# x0 = c(-2,-2,-2)

# máximo critério de convergência:
cc_max = 0.0000000001
cc = 10

otimizador_newton_rapson <- function(f0, f1x, f1y, f1lambda, f2xx, f2yy, f2lambdalambda, f2xy, f2yx, f2xlambda, f2ylambda, f2lambdax, f2lambday, cc_max, x0)
{
  # otimiza apenas funcoes de DUAS variáveis (x,y)
  # f0 é a função a ser otimizada
  # f1x é a primeira derivada de f0 em relação à x
  # f1y é a primeira derivada de f0 em relação à y
  # f1lambda é a primeira derivada de f0 em relação à lambda
  # f2xx é a segunda derivada de f0 em relação à xx
  # f2xy é a segunda derivada de f0 em relação à xy
  # f2yy é a segunda derivada de f0 em relação à xy
  # f2yx é a segunda derivada de f0 em relação à yx
  # cc é um escalar com o máximo do critério de convergência
  # X0 é vetor com os chutes iniciais
  
  conta = 0
  
  tamanho = length(x0)
  
  while(cc > cc_max)
  {
    
    # cria Hessiana:
    hessiana = matrix(0, ncol=length(x0), nrow=length(x0))

    x = x0[1]
    y = x0[2]
    lambda = x0[3]
    
    hessiana[1,1] = f2xx(x,y,lambda)
    hessiana[2,2] = f2yy(x,y,lambda)
    hessiana[3,3] = f2lambdalambda(x,y,lambda)
    
    hessiana[1,2] = f2xy(x,y,lambda)
    hessiana[2,1] = f2yy(x,y,lambda)
    
    hessiana[1,3] = f2xlambda(x,y,lambda)
    hessiana[3,1] = f2lambdax(x,y,lambda)
    
    hessiana[2,3] = f2ylambda(x,y,lambda)
    hessiana[3,2] = f2lambday(x,y,lambda)
    
    # cria gradiente da f:
    gradiente = c( f1x(x,y,lambda), f1y(x,y,lambda), f1lambda(x,y,lambda) )
    
    # aplica newton-rapson:
    # solve inverte matriz, %*% é multiplicação matricial.
    x1 = x0 - solve(hessiana) %*% gradiente
    
    # controle do loop:
    cc = sum((x1-x0)^2)
    x0 = x1
    conta = conta + 1
    print(paste('iteracao:', conta))
    print(paste('x0:', x0))
  }
  
  return(x0)
}

resposta = otimizador_newton_rapson(f0, f1x, f1y, f1lambda, f2xx, f2yy, f2lambdalambda, f2xy, f2yx, f2xlambda, f2ylambda, f2lambdax, f2lambday, cc_max, x0)
print('Pontos encontrados (x, y, lambda):')
print(resposta)