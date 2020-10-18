# Implementação dos métodos computacionais do artigo 'Bootstrapping for Significance of Compact Clusters in Multidimensional Datasets',
# Maitra et al. (2012)
# Trabalho submetido como requisito parcial para a conclusão da disciplina de Estatística Computacional do Mestrado em Estatística do PPGEst UFRGS.
# outubro de 2020.
# Os nomes das variaveis do codigo abaixo seguem a mesma notacao utilizada no artigo.
# Autor: Alexandre Gandini
# repositorio: https://github.com/alexgand/computacional_ppgest_2020

# bibliotecas essenciais:
%matplotlib
import numpy as np
import pandas as pd
from pandas import DataFrame, Series
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import datasets
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

# fixa seed para permitir reproducibilidade:
np.random.seed(seed=0)
random_state = 0

#######
# DADOS:
#######

# a) para aplicacao com dados 'reais', utilizamos o Wine dataset:
# Nesse caso, apos gerar os dados abaixo, rodar o codigo somente a partir da linha 54 ateh 227, para evitar as replicados da simulacao de monte carlo.

wine = datasets.load_wine()
X = wine.data
y = wine.target
data = DataFrame(X)
# padronizacao das variaveis para ficarem em uma mesma escala:
scaler = StandardScaler()
data = DataFrame(scaler.fit_transform(data), index=data.index, columns=data.columns)
# visualizacao dos dados apenas ateh a quarta variavel, pra nao poluir tanto o grafico:
plot = sns.pairplot(data.iloc[:,0:4])
cols = data.columns

# ou:

# b) com dados simulados:

# quantidade de replicacoes:
reps = 1000
acertou = []

for rep in tqdm(range(reps)):

    # alterar parametros abaixo: tamanho da amostra, quantidade de grupos, quantidade de variaveis, e o desvio padrao de cada grupo:
    X, y = make_blobs(n_samples=200, centers=2, n_features=2, cluster_std=0.4)
    # visualizacao dos dados gerados:
    plt.scatter(X[:,0], X[:,1]) 
    data = DataFrame(X, columns = range(X.shape[1]))
    # padronizacao das variaveis para ficarem em uma mesma escala:
    scaler = StandardScaler()
    data = DataFrame(scaler.fit_transform(data), index=data.index, columns=data.columns)
    plt.scatter(data.iloc[:,0], data.iloc[:,1]) 
    cols = data.columns

    #######

    # quantidade de grupos a testar:
    n_clusters = [1,2,3,4,5]
    # n_clusters = [1,2,3]

    # quantidade de replicacoes da amostra via bootstrap:
    n = 50

    # inicializacao de variaveis:
    wk = {}
    XI = {}

    # para cada possivel quantidade de grupos:
    for n_cluster in tqdm(n_clusters):

        # encontra os centros dos grupos via K-Means:
        clusterer = KMeans(n_clusters=n_cluster)
        clusterer.fit(data)

        # visualizacao com os grupos:
        # plt.scatter(data.iloc[:,0],data.iloc[:,1],c=clusterer.labels_)

        mu = clusterer.cluster_centers_
        tol = 1e-10
        mu[np.abs(mu) < tol] = 0.0
        mu = DataFrame(mu,columns=data.columns)
        mu = mu[cols]
        mu.index.name = 'grupo'

        data['grupo'] = clusterer.labels_
        std = data.groupby('grupo').std()

        mapping_mu = dict()
        for i in mu.index:
            mapping_mu[i] = mu.loc[i].values

        # medias = DataFrame(data['grupo'].map(mapping_mu).tolist(), columns=['mu' + str(x) for x in mu.columns])
        medias = DataFrame(data['grupo'].map(mapping_mu).tolist(), columns=mu.columns)

        erro = data[cols] - medias[cols]

        sigma_hat = 0
        wk[n_cluster] = 0
        for i in range(len(data)):
            sigma_hat = sigma_hat + np.matmul(erro.loc[i].T, erro.loc[i])
            wk[n_cluster] = wk[n_cluster] + np.matmul(erro.loc[i].T, erro.loc[i])

        sigma_hat = 1 / (len(data) * len(cols) ) * sigma_hat

        mu_scalled = mu / sigma_hat

        mapping_mu_scalled = dict()
        for i in mu_scalled.index:
            mapping_mu_scalled[i] = mu_scalled.loc[i].values

        # medias = DataFrame(data['grupo'].map(mapping_mu_scalled).tolist(), columns=['mu' + str(x) for x in mu.columns])
        mu_df_scalled = DataFrame(data['grupo'].map(mapping_mu_scalled).tolist(), columns=mu_scalled.columns)

        ei_hat = data[cols] - mu_df_scalled

        #######
        # Agora gera as amostras de bootstrap dos dados:
        #######

        XI[n_cluster] = []

        for i in range(n):

            zi = np.random.multivariate_normal(mean=np.zeros(shape=(mu.shape[1])), cov=np.identity(n=mu.shape[1]), size=len(data))

            wi = DataFrame()
            for col in DataFrame(zi):
                wi[col] = DataFrame(zi)[col] / Series(np.linalg.norm(zi, axis=1))

            # checando se sao vetores unitarios:
            # np.linalg.norm(zi, axis=1)

            data_perm = data.sample(n=len(data))

            mu_df_scalled = DataFrame(data_perm['grupo'].map(mapping_mu_scalled).tolist(), columns=mu_scalled.columns)

            ei_hat = data_perm[cols] - mu_df_scalled.values

            ei_star = DataFrame()
            for col in wi:
                ei_star[col] = Series(np.linalg.norm(ei_hat, axis=1)) * wi[col]

            x_star = DataFrame(data_perm['grupo'].map(mapping_mu).tolist(), columns=mu.columns)

            to_add = sigma_hat * ei_star

            # to_add.columns = x_star.columns
            to_add.columns = x_star.columns

            x_star = x_star + to_add

            # plt.figure(i)
            # plt.scatter(x_star[2],x_star[3],c=i+data_perm['grupo'])
            # plt.show()

            XI[n_cluster].append(x_star)

    ######
    # Agora gera a estatistica de teste e o seu valor-p, para estimar a quantidade verdadeira de grupos:
    ######

    wk_resampled = {}

    for n_cluster in tqdm(XI):
        
        wk_resampled[n_cluster] = {}

        counter = 0
        
        for data_resampled in XI[n_cluster]:

            clusterer = KMeans(n_clusters=n_cluster)
            clusterer.fit(data_resampled)

            # plt.scatter(data_resampled[2],data_resampled[3],c=clusterer.labels_)

            mu = clusterer.cluster_centers_
            tol = 1e-10
            mu[np.abs(mu) < tol] = 0.0
            mu = DataFrame(mu,columns=data_resampled.columns)
            mu = mu[cols]
            mu.index.name = 'grupo'

            data_resampled['grupo'] = clusterer.labels_
            std = data_resampled.groupby('grupo').std()

            # mapping_mu = dict()
            # for i in mu.index:
            #     mapping_mu[i] = mu.loc[i].values

            mu_scalled = mu / sigma_hat

            mapping_mu_scalled = dict()
            for i in mu_scalled.index:
                mapping_mu_scalled[i] = mu_scalled.loc[i].values
            
            # medias = DataFrame(data_resampled['grupo'].map(mapping_mu).tolist(), columns=['mu' + str(x) for x in mu.columns])
            # medias = DataFrame(data_resampled['grupo'].map(mapping_mu).tolist(), columns=mu.columns)
            mu_df_scalled = DataFrame(data_resampled['grupo'].map(mapping_mu_scalled).tolist(), columns=mu_scalled.columns)

            erro = data_resampled[cols] - medias[cols]
            erro = data_resampled[cols] - mu_df_scalled[cols]

            # sigma_hat = 0
            wk_resampled[n_cluster][counter] = 0
            for i in range(len(data_resampled)):
                # sigma_hat = sigma_hat + np.matmul(erro.loc[i].T, erro.loc[i])
                wk_resampled[n_cluster][counter] = wk_resampled[n_cluster][counter] + np.matmul(erro.loc[i].T, erro.loc[i])
            
            counter += 1

    skk = Series(wk).diff()
    skk_resampled = DataFrame(wk_resampled).T.diff()
    valor_p = (skk_resampled > skk).sum(axis=1) / n

    print('valor-p para cada quantidade de grupos testada:')
    print(valor_p)

    sim = valor_p.loc[1] < valor_p.loc[2]
    sim = sim and (valor_p.loc[1] < 0.05)
    sim = sim and (valor_p.loc[2] > 0.05)

    acertou.append(sim)

print('Proporcao de vezes em que o procedimento acertou a quantidade correta de grupos:')
print(Series(acertou).sum() / len(acertou))
