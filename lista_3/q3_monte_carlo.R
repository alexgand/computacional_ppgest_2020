# ALEXANDRE GANDINI
# questao 3, lista 3, Computacional PPGEST-UFRGS 2020.

# limpa variaveis na memoria:
rm(list=ls())
par(mfrow=c(1,1)) # limpa layout

#####
# a) Estimar E(X) e calcular erro padrao da estimativa:
#####

RE = 10000

# funcao pra gerar o menor n tal que a soma de n uniformes[0,1] seja menor do que 1:
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

# Se formos utilizar o metodo da letra a), praticamente nao temos amostras de X >= 10, que eh um
# caso que requer que os valores das uniformes geradas sejam bem pequenos, jah que X eh a
# quantidade de uniformes somadas tal que sua soma nao ultrapasse 1.

# Uma alternativa eh utilizar importance sampling com uma funcao de importancia como a beta, por exemplo,
# que, escolhendo parametros apropriados, gera valores pequenos com alta probabilidade.

# Assim:
# f(x) = 1 # densidade da uniforme[0,1]
# h(x) = I[ (U1 + U2 + ... + U9) < 1 ] # indicadora de que a soma de 9 uniformes sejam menores do que 1, pois queremos estimar P(X >= 10).
# g(x) = produtorio das densidades de 9 betas com parametros que dao probabilidade alta para valores pequenos.

# parametros da beta escolhidos:
shape1 = 1.6
shape2 = 12

# exemplos de valores gerados com essa beta, sao pequenos:
rbeta(n=10,shape1=shape1, shape2=shape2)

# gera 9 realizacoes dessa beta, de tamanho RE:
betas = matrix(nrow=RE, ncol=9)
for (i in 1:9)
{
  betas[,i] = rbeta(n=RE, shape1=shape1, shape2=shape2 )
}

# densidade da uniforme:
fx = 1

# indicadora da soma das variaveis serem menores do que 1:
hx = rowSums(betas) <= 1

gx = apply(dbeta(betas, shape1=shape1, shape2=shape2), 1, prod)

estimativa = hx * fx / gx

print(paste('P(X >= 10) estimado:', mean(estimativa) ))
print(paste('Erro padrao estimado:', sqrt( var(estimativa) / RE )  ))


#####
# c) Estimar P(X = 10) e calcular erro padrao da estimativa:
#####

# Agora, na funcao h(x), alem da indicadora de que (U1 + ... + U9) < 1, tambem
# queremos que (U1 + ... + U10) > 1, assim cumprimos o requisito que que X eh
# exatamente igual aa 10.

# Assim:
# f(x) = 1 # densidade da uniforme[0,1]
# h(x) = I[ ((U1 + ... + U9) < 1) & ((U1 + ... + U10) > 1)] # indicadora de que a soma de 9 uniformes sejam menores do que 1, e ao mesmo tempo, de que a soma de 10 uniformes seja maior do que 1.
# g(x) = produtorio das densidades de 10 betas com parametros que dao probabilidade alta para valores pequenos.

# parametros da beta escolhidos:
shape1 = 1.6
shape2 = 12

# exemplos de valores gerados com essa beta, sao pequenos:
rbeta(n=10, shape1=shape1, shape2=shape2)

# gera 10 realizacoes dessa beta, de tamanho RE:
betas = matrix(nrow=RE, ncol=10)
for (i in 1:10)
{
  betas[,i] = rbeta(n=RE, shape1=shape1, shape2=shape2 )
}

# densidade da uniforme:
fx = 1

# indicadora:
hx = (rowSums(betas[,1:9] < 1) & (rowSums(betas) > 1) )

gx = apply(dbeta(betas, shape1=shape1, shape2=shape2), 1, prod)

estimativa = hx * fx / gx

print(paste('P(X = 10) estimado:', mean(estimativa) ))
print(paste('Erro padrao estimado:', sqrt( var(estimativa) / RE )  ))
