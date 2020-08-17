# ALEXANDRE GANDINI
# questao 1, lista 2, Computacional PPGEST-UFRGS 2020.

# O algoritmo descrito no problema foi implementado da seguinte forma:
  
# Inicializa a media dos dois grupos randomicamente, dentro da amplitude de valores do dataset

# Itera ateh convergencia:

# 1. Calcula as distancias euclidianas de cada ponto ateh as medias dos dois grupos
# 2. Atribui a cada ponto o grupo para o qual apresenta a menor distancia
# 3. Recalcula as medias dos grupos com base nas novas atribuicoes feitas no item anterior

# No caso, foi utilizado como criterio de convergencia para finalizar o algoritmo quando 
# nao houve mais alteracao de nenhuma atribuicao de grupo para qualquer ponto, de uma interacao
# para outra.

# implementacao o K-means para apenas dois grupos:

# limpa variaveis na memoria:
rm(list=ls())
par(mfrow=c(1,1)) # limpa layout

n = 100 # quantidade de pontos de cada grupo
mean_grupo_A = 1
mean_grupo_B = 3
desvio_padrao = 1

# gerando dados simulados para dois grupos (normais bivariadas com medias diferente):
XA = rnorm(n, mean=mean_grupo_A, sd=desvio_padrao)
YA = rnorm(n, mean=mean_grupo_A, sd=desvio_padrao)
XB = rnorm(n, mean=mean_grupo_B, sd=desvio_padrao)
YB = rnorm(n, mean=mean_grupo_B, sd=desvio_padrao)

# visualizacao dos grupos:
plot(XA,YA, col=4, xlim=c(-2,6), ylim=c(-2,6), pch = 19, xlab='X',ylab='Y')
points(XB,YB, col=2, pch = 19)

# criando dataset unico:
dataset = matrix(nrow = length(XA) + length(XB), ncol=2)
dataset[,1] = rbind(XA,XB)
dataset[,2] = rbind(YA,YB)

# chute inicial da media dos grupos:
g1 = c(0,0)
g1[1] = runif(n=1, min=min(XA), max=max(XA))
g1[2] = runif(n=1, min=min(YA), max=max(YA))

g2 = c(0,0)
g2[1] = runif(n=1, min=min(XB), max=max(XB))
g2[2] = runif(n=1, min=min(YB), max=max(YB))

# visualizando as medias iniciais:
points(g1[1], g1[2], col=3, pch=8, lwd=10)
text(g1[1], g1[2], labels=0,cex=0.9, font=2)
points(g2[1], g2[2], col=5, pch=8, lwd=10)
text(g2[1], g2[2], labels=0,cex=0.9, font=2)

# inicializa variavel pra guardar os grupos:
grupos = rep(0,dim(dataset)[1])
new_grupos = rep(1,dim(dataset)[1])
conta = 0
parar = FALSE

# continua iterando ateh que nao tenha mais troca de grupos:
while ((parar == FALSE) )
{
  # calculando as distancias de cada ponto pras médias iniciais:
  distancias = matrix(0, nrow=dim(dataset)[1], ncol=2)
  distancias[,1] = sqrt(rowSums((dataset - g1)**2))
  distancias[,2] = sqrt(rowSums((dataset - g2)**2))
  
  # atribui grupo a media mais próxima:
  new_grupos = apply(distancias, 1, which.min)
  
  # controle do loop:
  if (all(new_grupos == grupos))  {
    parar = TRUE}  else { grupos = new_grupos }
  
  # recalcula as medias:
  novas_medias = aggregate(cbind(dataset,grupos),by=list(grupos),FUN='mean')
  g1 = c(novas_medias[1,2], novas_medias[1,3])
  g2 = c(novas_medias[2,2], novas_medias[2,3])

  # visualiza as novas medias:
  points(g1[1], g1[2], col=3, pch=8, lwd=10)
  text(g1[1], g1[2], labels=conta,cex=0.9, font=2)
  points(g2[1], g2[2], col=5, pch=8, lwd=10)
  text(g2[1], g2[2], labels=conta,cex=0.9, font=2)

  print(paste('iteração:',conta))
  print(paste('Média g1:',g1))
  print(paste('Média g2:',g2))
    
  conta = conta + 1
  
  }

print('grupos:')
print(grupos)

print(paste('medias dos grupos apos',conta,'interacoes:'))
print('grupo 1:')
print(g1)
print('grupo 2:')
print(g2)