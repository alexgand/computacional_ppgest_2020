import numpy as np
import pandas as pd
from tqdm import tqdm

n_times = 1000

sizes = [10,100,1000,10000]

matrix_positiva = []

for size in tqdm(sizes):
    for run in tqdm(range(n_times)):
        # inicializa a matrix com normais padrao
        matrix = np.random.normal(size=(size,size))

        # torna a matriz simetrica, mantendo a distribuicao:
        # soma triangula inferior com a transposta do triangulo inferior excluída a diagonal:
        # https://stackoverflow.com/questions/10806790/generating-symmetric-matrices-in-numpy
        sym_matrix = np.tril(matrix) + np.tril(matrix, -1).T

        # teste, transforma ela em positiva AT*A:
        # sym_matrix = np.matmul(sym_matrix.T, sym_matrix)
        # multiplica diagonal principal por numero grande:
        # sym_matrix = sym_matrix * ((np.identity(size) * 10000) + 1)
        # teste matriz com elementos positivos:
        # sym_matrix = np.ones((3,3))

        # via autovalores:
        # autovalores = np.linalg.eigvals(sym_matrix)
        # matrix_positiva.append(all(autovalores > 0))

        # via leading principal minors:
        # dets = []
        # for i in range(1, size+1):
        #     dets.append(np.linalg.det(sym_matrix[:i,:i]) > 0)
        # matrix_positiva.append(all(np.array(dets) > 0))
            
        # via decomposição de cholesky:
        # https://stackoverflow.com/questions/16266720/find-out-if-matrix-is-positive-definite-with-numpy
        try:
            np.linalg.cholesky(sym_matrix)
            matrix_positiva.append(True)
            print('ACHEI uma POSITIVA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        except :
            matrix_positiva.append(False)
    
    print('size:',size,'proporção de matrizes positivo-definidas:',str(np.array(matrix_positiva).sum() / len(matrix_positiva)))