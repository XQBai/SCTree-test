# author: Xiangqi Bai <xqbai@amss.ac.cn>

import warnings
import sys
import time
import pandas as pd
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    import wishbone

from sklearn.neighbors import NearestNeighbors
from scipy import stats
import numpy as np
import math
import random
import scipy.sparse as sp
from copy import deepcopy

def test(input_data, num_way, startcell, data_type, n=2, verbose = True):

    Data = deepcopy(input_data)
    '''

    :param Data: Dataframe of cells X expression genes
    :param num_way: The dimension of matrix D which is that subsampling size
    :param startcell: piror start cell or a random cell in data
    :param data_type: 'sc-seq'
    :param n: the number of k in k-nearest graph equals to kmin*n
    :param verbose: if the startcell is known, verbose = True, otherwise, verbose = False
    :return: p-value and SNR （signal-to-noise ratio）
    '''
    # Start cell index
    start = np.where(Data.index == startcell)[0]
    if len(start) == 0:
        raise RuntimeError('Start cell %s not found in data. Please rerun with correct start cell' % startcell)

    if not (isinstance(Data, pd.DataFrame)):
        raise TypeError('data must be of type or DataFrame')

    if not data_type in ['sc-seq']:
        raise RuntimeError('data_type must be sc-seq')

    if num_way > Data.shape[0]:
        raise RuntimeError('num_waypoints parameter is higher than the number of cells in the dataset. \
                    Please select a smaller number')

    start_time = time.process_time()


    #If the dimension is less than 20, PCA reduction dimensionality
    if Data.shape[1] > 20:
        pca_data = pca_reduce(Data)[0]
        com_list = list(range(pca_data.shape[1]))
    else:
        pca_data = Data
        com_list = list(range(Data.shape[1]))
    #Construct the knn graph
    k, knn = _knn(pca_data, com_list, n)

    if verbose == True:
        s = start
    else:
        dist = sp.csgraph.dijkstra(knn, directed=False, indices=start[0])
        s = np.where(dist == np.max(dist))[0]

    #Randomly select some points to construct discordance matrix DC
    subsample_points = random.sample(list(range(0, len(Data))), num_way)

    #Compute the discordance matrix DC
    DC = discordance_matrix(knn, subsample_points, s[0])

    # Transform the discordance matrix to Wigner-like matrix
    DC1, SNR = Normalization(pca_data, DC, subsample_points, s[0])
    DC1 = np.nan_to_num(DC1)
    DC1[np.isinf(DC1)] = 0
    DC3 = Permutation(DC1)

    eigs, eigvector = np.linalg.eig(DC3)


    v = list()
    # Bootstrap the eigenvalues
    for i in range(len(subsample_points)):
        v.extend(random.sample(list(eigs), 1))

    # normalize eigenvalues
    v = np.asarray(v)
    v = v / (math.sqrt(len(subsample_points)) * 2)

    #Non-parametric hypothesis test

    p_value = stats.kstest(v, 'semicircular')[1]

    print('p-value: %.2f  ' %p_value + 'SNR: %.2f' % SNR)
    print('Computed in：%.2f seconds' % (time.process_time() - start_time))

    return p_value, SNR


def from_csv(csv_file, data_type, cell_axis=0, normalize=True):
    '''
    csv_file: input data with n cells and p dimensions，each row represents the expressions of genes for a cell
    data_type: 'sc-seq'
    '''

    if not data_type in ['sc-seq']:
        raise RuntimeError('data_type must be sc-seq')

    # Read csv_file
    data = pd.read_csv(csv_file, sep=None, header=0, index_col=0, engine='python')

    # Dataframe for sc-data which first column is the cells
    if cell_axis != 0:
        data = data.transpose()

    # Normalize the single cell RNA-seq data
    if data_type == 'sc-seq' and normalize == True:
        data = normalize_scseq(data)

    return data


def normalize_scseq(data):
    '''
    Normalize single cell RNA-seq data: Divide each cell by its molecule count
    and multiply counts of cells by the median of the molecule counts

    This part code is referenced from Wishbone algorithm wb.py
    '''

    if not (isinstance(data, pd.DataFrame)):
        raise TypeError('data must be of type or Dataframe')

    # Compute the sum of gene expression for one cell
    molecule_counts = data.sum(axis=1)
    data = data.div(molecule_counts, axis=0) \
        .mul(np.median(molecule_counts), axis=0)

    # check that none of the genes are empty; if so remove them
    nonzero_genes = data.sum(axis=0) != 0
    data = data.ix[:, nonzero_genes].astype(np.float32)

    return data


def pca_reduce(data):


    """
    Principal component analysis of the data.

    This part code is referenced from Wishbone algorithm wb.py

    """
    print('PCA for dataset...')
    X = data.values
    # Make sure data is zero mean
    X = np.subtract(X, np.amin(X))
    X = np.divide(X, np.amax(X))

    # Compute covariance matrix
    if (X.shape[1] < X.shape[0]):
        C = np.cov(X, rowvar=0)

    # if N>D, we better use this matrix for the eigendecomposition
    else:
        C = np.multiply((1 / X.shape[0]), np.dot(X, X.T))

    # Perform eigendecomposition of C
    C[np.where(np.isnan(C))] = 0
    C[np.where(np.isinf(C))] = 0
    l, M = np.linalg.eig(C)

    # select PCA components
    #n_components = 100
    n_components = np.min(np.where(np.abs(np.diff(np.diff(l))) < 0.0001)[0] + 1)

    # Sort eigenvectors in descending order
    ind = np.argsort(l)[::-1]
    l = l[ind]
    if n_components < 1:
        n_components = np.where(np.cumsum(np.divide(l, np.sum(l)), axis=0) >= n_components)[0][0] + 1
        print('Embedding into ' + str(n_components) + ' dimensions.')
    if n_components > M.shape[1]:
        n_components = M.shape[1]
        print('Target dimensionality reduced to ' + str(n_components) + '.')

    M = M[:, ind[:n_components]]
    #ACR = np.sum(l[:n_components]) / np.sum(l)
    l = l[:n_components]
    # Apply mapping on the data
    if X.shape[1] >= X.shape[0]:
        M = np.multiply(np.dot(X.T, M), (1 / np.sqrt(X.shape[0] * l)).T)

    loadings = pd.DataFrame(data=M, index=data.columns)
    #l = pd.DataFrame(l)

    data -= np.min(np.ravel(data))
    data /= np.max(np.ravel(data))

    pca = pd.DataFrame(np.dot(data, loadings), index=data.index)

    return pca, l


def _knn(pca, com_list, n, metric='euclidean'):

    '''
    Using the PCA components to construct k-nearest neighbor graph
    '''
    print('Buliding knn graph...')
    kmin = 2
    pca = pca.ix[:, com_list]
    nbrs = NearestNeighbors(n_neighbors=kmin, metric=metric).fit(pca)
    knn = nbrs.kneighbors_graph(pca, mode='distance')
    connected_components = sp.csgraph.connected_components(knn, directed=False)[0]
    while connected_components is not 1:
        kmin = kmin + 2
        nbrs = NearestNeighbors(n_neighbors=kmin, metric=metric).fit(pca)
        knn = nbrs.kneighbors_graph(pca, mode='distance')
        connected_components = sp.csgraph.connected_components(knn, directed=False)[0]

    nbrs = NearestNeighbors(n_neighbors= n * kmin, metric=metric).fit(pca)
    knn = nbrs.kneighbors_graph(pca, mode='distance')
    knn = np.transpose(knn)
    knn.setdiag(0)
    knn.eliminate_zeros()
    k = n * kmin

    return k, knn

def discordance_matrix(knn, waypoints, startcell):

    print('Construct distance matrix DC...')
    Gromov_transform = np.zeros((len(waypoints), len(waypoints)))
    DC = np.zeros((len(waypoints), len(waypoints)))

    # Find initial trajectory
    dist = sp.csgraph.dijkstra(knn, directed=False, indices= startcell)

    # Find shortest distance from waypoints to all points
    assert isinstance(waypoints, object)
    dist_w = sp.csgraph.dijkstra(knn, directed=False, indices=waypoints)

    for i in range(len(waypoints)):
        for j in range(len(waypoints)):
            # Gromov-Farris transformation
            Gromov_transform[i, j] = 1/2 * (dist_w[i, waypoints[j]] - dist[waypoints[i]] - dist[waypoints[j]])

            # Discordance distance matrix
            DC[i, j] = 2 * np.abs(np.min([dist[waypoints[i]], dist[waypoints[j]]]) + Gromov_transform[i, j])

    return DC


# Standardize the independent matrix
def Normalization(pca_data, DC, waypoints, startcell):

    print('Normalize the discordance matrix DC...')
    n = len(DC[0, :])

    # Estimate the signal population and noise population in discordance matrix
    com_list = list(range(pca_data.shape[1]))
    __console__ = sys.stdout
    sys.stdout = open('/dev/null', 'w')
    k, knn = _knn(pca_data, com_list, n = 2)
    res = wishbone.core.wishbone(pca_data.ix[:, com_list].values, startcell, k=2 * k, l=2 * k, num_waypoints=100)
    branches = res['Branches']

    t = list(np.where(branches == 1)[0])
    b1 = list(np.where(branches == 2)[0])
    b2 = list(np.where(branches == 3)[0])

    w1 = [i for i in t if i in waypoints]
    w2 = [i for i in b1 if i in waypoints]
    w3 = [i for i in b2 if i in waypoints]
    w = w1 + w2 + w3

    D_w = discordance_matrix(knn, w, startcell)
    sys.stdout = __console__

    D_signal = [D_w[i, j] for i in range(len(w1), len(w1 + w2)) for j in range(len(w1 + w2), len(w))]
    D_noise1 = [D_w[i, j] for i in range(len(w1)) for j in range(len(w))]
    D_noise2 = [D_w[i, j] for i in range(len(w1), len(w1 + w2)) for j in range(len(w1 + w2))]
    D_noise3 = [D_w[i, j] for i in range(len(w1 + w2), len(w)) for j in range(len(w1))]
    D_noise4 = [D_w[i, j] for i in range(len(w1 + w2), len(w)) for j in range(len(w1 + w2), len(w))]
    D_noise = D_noise1 + D_noise2 + D_noise3 + D_noise4

    D_noise = np.nan_to_num(D_noise)
    D_noise[np.isinf(D_noise)] = 0
    D_signal = np.nan_to_num(D_signal)
    D_signal[np.isinf(D_signal)] = 0
    mu = np.mean(D_noise)
    sigma = np.std(D_noise)

    # Estimate the signal to noise ratio
    SNR = np.mean(D_signal) / sigma

    # Compute the normalization of matrix D
    DC1 = np.zeros((n, n))
    for i in range(n):
        for j in range(i):
            DC1[i, j] = (DC[i, j] - mu) / sigma

    return DC1, SNR

# Make matrix's elements independent
def Permutation(D1):

    print('Permuate the discordance matrix DC...')
    # Exact the lower triangular elements of the discordance matrix D
    lower_triangular = list()
    n = len(D1[0, :])
    for i in range(n):
        for j in range(i):
            lower_triangular.append(D1[i, j])

    # Disrupt the order of the list
    random.shuffle(lower_triangular)

    # Make a lower triangular matrix by out of order list
    DC2 = np.zeros((n, n))
    for i in range(n):
        for j in range(i):
            DC2[i, j] = lower_triangular.pop()

    # Make the matrix symmetrical
    rand_N = np.random.randn(len(DC2[0, :]))
    DC3 = DC2 + DC2.T
    for i in range(len(DC2[0, :])):
        DC3[i, i] = rand_N[i]

    return DC3




