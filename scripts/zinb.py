"""
This script is not used.
This is my first attempt at the ZINB.
"""


import numpy as np
import scipy.optimize
import scipy.stats


def nll_zinb(params, y):
    pi, n, p = params
    # pi = probability that the 0 values are extraneous
    # p = success probability, getting a transposon
    # n = number of hits/successes
    # y = number of failures before n successes, the TA count
    # ZINB (Equation 2)
    n_reads = pi + (1.0-pi)*scipy.stats.nbinom.pmf(y[y==0], n, p)
    y_reads = (1.0-pi)*scipy.stats.nbinom.pmf(y[y>0], n, p)
    # Compute the negative log-likelihood
    """ https://stackoverflow.com/questions/5124376/convert-nan-value-to-zero/5124409 """ 
    n_reads = np.nan_to_num(np.log(n_reads+np.finfo(float).eps))
    y_reads = np.nan_to_num(np.log(y_reads+np.finfo(float).eps))
    nll = -(np.sum(n_reads) + np.sum(y_reads))
    print(nll)
    return nll


def zinb_fit(data, debug=False):
    # Scipy function requires integers
    # rint converts to hte nearest integer
    data = np.rint(np.array(data)).astype(np.int)
    # Inital values of parameters to estimate
    pi = 0.001  # 10% of zeros are extraneous
    n = np.mean(data[data>0])  # guess around 10 sucesses before the TA count
    p = 0.5  # 5% probability of a transposon, helps keep the probability higher for very high counts
    params = [pi, n, p]
    eps = 1e-3  # using zero gives warnings, so this is a small epsilon value
    #         pi=[0,1]      n>0          p=[0,1]
    bounds = [(eps, 1-eps), (eps, None), (eps, 1-eps)]
    results = scipy.optimize.minimize(nll_zinb, params, args=(data,), bounds=bounds)
    pi, n, p = results.x
    mu = n*(1-p)/p
    var = n*(1-p)/(p*p)
    if debug:
        print("data :", data); print("pi :", pi); print("n :", n)
        print("p :", p); print("mu :", mu); print("var :", var)
    return params


if __name__ == "__main__":

    data = np.array( [20, 0, 0, 0, 0, 0] )
    # data = np.array( [2540, 10170] )

    params = zinb_fit(data, debug=True)
    
    print(np.mean(data[data>0]))

