# inicialização de variáveis:
n_times = 1000
#sizes = c(10,100,1000,10000)
sizes = c(10,100,1000)
matrix_positiva = c()

# funcao para testar se a matriz é positivo definida, por dois métodos: autovaloes e lpm - leading principal minors:

testa_matriz_positivo_definida <- function(matrix, metodo)
{
  
  if (metodo == 'autovalores')
    {
    # autovalores - se todos os autovalores forem positivos, a matriz é positivo-definida:
    ev <- eigen(matrix)
    autovalores <- ev$values
    if ( sum(autovalores >= 0) == length(autovalores) )
    {return(TRUE)} else
    {return(FALSE)}
    }
  
  if (metodo == 'lpm')
    {
    # via leading principal minors - se todos os determinantes forem positivos, a matriz é positivo-definida:
    dets = c()
    for (i in seq(1, size))
      if (i == 1)
      { dets = c(dets, matrix[1:i,1:i] > 0) } else
      { dets = c(dets, det(matrix[1:i,1:i]) > 0) }
    if ( sum(dets == TRUE) == length(dets) )
      {return(TRUE)} else
      {return(FALSE)}
    }
}

# roda os testes conforme solicitado no exercício 3 da lista 1:

for (size in sizes)
{
  time0 = proc.time()
  for (run in seq(n_times))
    {
    # inicializa a matrix com normais padrao
    matrix = matrix( rnorm(n=size*size, mean=0, sd=1), size, size)
    
    # teste com matriz identidade:
    # matrix = diag(10)
    
    # torna a matriz simetrica, mantendo a distribuicao:
    # soma triangulo inferior com a transposta do triangulo inferior excluída a diagonal:
    matrix[lower.tri(matrix)] <- t(matrix)[lower.tri(matrix)]

    # testa matriz com a funcao definida acima:
    # resposta = testa_matriz_positivo_definida(matrix, metodo='autovalores')
    resposta = testa_matriz_positivo_definida(matrix, metodo='lpm')
    matrix_positiva = c(matrix_positiva, resposta)
    }

  time1 = proc.time()
  print(paste("size:",size))
  print(paste("proporção de matrizes positivo-definidas:",sum(matrix_positiva) / length(matrix_positiva)))
  print(paste("Tempo decorrido em segundos:"))
  print(time1-time0)
}
