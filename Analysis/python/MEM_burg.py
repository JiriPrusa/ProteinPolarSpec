"""
file              CollombBurg.py

author/translator Ernesto P. Adorio,Ph.D.
                  UPDEPP (UP Clark, Pampanga)
                  ernesto.adorio@gmail.com

Version           0.0.1 jun 11, 2010 # first release.

References        Burg's Method, Algorithm and Recursion, pp. 9-11
"""
from math import floor
import numpy as np


def burg_AR(m, X):
    """
    Based on Collomb's C++ code, pp. 10-11
    Burgs Method, algorithm and recursion
      m - number of lags in autoregressive model.
      x  - data vector to approximate.
    """
    L, N = np.shape(X)
    N = N - 1
    coeffs = np.zeros(m)

    # initialize Ak
    Ak = np.zeros(m + 1)
    mu_arr = np.zeros(m)
    Ak[0] = 1.0
    # initialize f and b.
    f = X.copy()
    b = X.copy()
    # Initialize Dk
    Dk = 0.0
    for i in range(L):
        for j in range(N + 1):
            Dk += 2.0 * f[i, j] ** 2
        Dk -= (f[i, 0] ** 2) + (b[i, N] ** 2)

    sigma = np.mean(np.power(X, 2))

    # Burg recursion
    for k in range(m):
        # compute mu
        mu = 0.0
        for i in range(L):
            for n in range(N - k):
                mu += f[i, n + k + 1] * b[i, n]
        mu *= -2.0 / Dk
        mu_arr[k] = mu

        # update Ak
        maxn = floor((k + 1) / 2 + 1)  # rounds down
        for n in range(maxn):
            t1 = Ak[n] + mu * Ak[k + 1 - n]
            t2 = Ak[k + 1 - n] + mu * Ak[n]
            Ak[n] = t1
            Ak[k + 1 - n] = t2
        # update f and b
        for i in range(L):
            for n in range(N - k):
                t1 = f[i, n + k + 1] + mu * b[i, n]
                t2 = b[i, n] + mu * f[i, n + k + 1]
                f[i, n + k + 1] = t1
                b[i, n] = t2

        # update Dk
        t1 = 0
        t2 = 0
        for i in range(L):
            t1 += f[i, k + 1] ** 2
            t2 += b[i, N - k - 1] ** 2
        # Dk = ( 1.0 - mu ** 2) * Dk - (f[ k + 1 ] ** 2) - (b[ N - k - 1 ] ** 2)
        Dk = (1.0 - mu ** 2) * Dk - t1 - t2
        sigma = sigma * (1 - Ak[-1] ** 2)

    # assign coefficients.
    coeffs = Ak[1:]

    return coeffs, mu_arr, sigma  # model coefficients, reflexion coefficients, variance of dataset

def cic(res, data_length):
    i = np.arange(len(res))
    i[0] = 1
    vi = 1 / ((data_length + 1) - i)
    fic_tail = 3 * np.cumsum(vi)
    fsic_tail = np.cumprod((1 + vi) / (1 - vi)) - 1
    cic_tail = np.amax((fic_tail, fsic_tail), 0)
    cic = np.log(res) + cic_tail
    pe_est = res * (fsic_tail + 1)
    return cic, pe_est
