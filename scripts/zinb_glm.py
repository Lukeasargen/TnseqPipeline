"""
Zero-Inflated Binomial Regression

Citation: Subramaniyam, S., DeJesus, M.A., Zaveri, A. et al. Statistical analysis of variability in TnSeq data across conditions using zero-inflated negative binomial regression. BMC Bioinformatics 20, 603 (2019). https://doi.org/10.1186/s12859-019-3156-z

I first made just the ZINB model. It is in a separate script.
I then added the GLM into this script.

I first implemented the formula from the paper (Equation 2).
For the probability mass function, I used scipy.
https://docs.scipy.org/doc/scipy/reference/tutorial/stats/discrete_nbinom.html
I used variable n and p to match the scipy docs.
The k in formula is the TA counts. This is the failures before 
n successes.

Then, I read the wikipedia to see how to estimate the parameters.
https://en.wikipedia.org/wiki/Negative_binomial_distribution#Maximum_likelihood_estimation
I modified the function to now output the negative log likelihood.
This is minimized to get the maximum likelihood estimates.
The reason why you take the log is because it turns multiplication
into summation. Assuming independence of the data points:
    p(x1, x2, x3...|θ) = p(x1|θ)*p(x2|θ)*p(x3|θ)*...
Under a log this becomes:
    log( π(p(θ|xi)) ) = Σ log(p(θ|xi))
Also, log is monotonically increasing, so you can compare the likelihood
between models after the log and the sign of the difference stays the same.

We want to fit the distribution to minimize the log probability for the
sum of all k values, this will give us estimates for
the n and p values, which are used to estimate the mean and variance.

I used the scipy formulas for mean and variance since my n and p are indentical:
    mu = n*(1-p)/p
    var = n*(1-p)/(p*p)

Generalized Linear Model
The next step is to implement GLM regression.

First I needed to make a binary design matrix.
This is essentially a one-hot encoded matrix. The wikipedia is really good.
https://en.wikipedia.org/wiki/Design_matrix#Multiple_regression

Second, I created a way to pack and unpack the paramters in a
1D list. Scipy optimize requires a 1D list of paramters.
pg is the first element since it is always a single value.
The rest of the list is length 2k, where k is the number of conditions.
This is because both mu and pi are estimated for each condition.

This is how I understand the parameters:
    y = TA counts
        This is the number of trials in the zinb distribution.
    pg = coefficient for dispersion parameter
        This is the number of successes in the zinb distribution.
        All condidtions are fit to the same pg. The dispersion is
        assumed to be the same for all conditions.
        As the stated in the paper, use n=np.exp(pg) to get the value
        for the dispersion. The formula ensures that n is positive and
        it allows pg to be any value (so no bounds).
    ag = coefficient for mean parameter
        mu is the mean of the zinb distribution and is calculated
        by mu = np.exp(Xg.dot(ag)). Same as pg, the mean is always
        positive and the expontential means ag can be any value.
        The difference from pg is a value for ag is fit for each
        condition. This is important for the statistical significance test.
    yg = coefficient for extraneous probability
        pi is the probability of the zero counts being extraneous and
        is calculated by pi = Xg.dot( np.exp(yg)/(1+np.exp(yg)) ).
        This function exp(x)/(1+exp(x)) is the sigmoid function and
        it domain is all values and has range of (0,1).
    n, p, and pi are calculated from these parameters and used in the zinb equation.

Statistical significance
Likelihood ratio test (LRT)
The fastest way to understand this is to read the wikipedia.
https://en.wikipedia.org/wiki/Likelihood-ratio_test

Step 1: fit 2 versions of the model.
m0) Null hypothesis model. This model assume all data comes from the
same distribution. This is implement by setting the condition to 0 for
all data points. This is a constrained version of the full model
m1) Alternative hypothesis model. This model fits a distribution for
each condition. This is the more complex model.
In short:
    condition-dependent ZINB model (m1)
    condition-independent ZINB model (m0) or null model

Step 2: compute the log-likelihood ratio.
It's easier to look at the pseudo-code than explain in text.
l1 = sum(log(m1(data, params1, conditions)))
l0 = sum(log(m1(data, params0, null_conditions)))
llr = -2*(l0-l1)  # log-likelihood ratio

Slight change is that I used the nll function from before, so
the llr is calculated by:
llr = 2*(l0-l1)

Step 3: find the p-value
As stated in the paper, the likelihood ratio (llr) is distributed
chi-squared. The degrees of freedom is the difference in the number
of parameters.
Chi-sqaured distribution is the only assumption of this test and it
is usually good, even for small sample sizes.

I used the scipy chi-sqaured pmf function to get the p-value.
https://docs.scipy.org/doc/scipy/reference/tutorial/stats/continuous_chi2.html

llr is essentailly a chi-sqaured value.
The test rejects the null if the llr is smaller than the 
chi-squared percentile (choosen by us, the analyst) for
the given degrees of freedom.

"""
import time

import numpy as np
from numpy.lib.financial import pv
from scipy.optimize import minimize
from  scipy.stats import nbinom, norm, chi2

# import warnings
# warnings.filterwarnings("ignore")


def num_conditions(conditions):
    return len(np.unique(conditions))


def design_matrix(conditions):
    """ Returns a binary design matrix of size m X k where
        m is the number of samples and
        k is the number conditions.
    """
    # Create the design matrix, essentially it's one hot encoding across the rows
    c = np.array(conditions)
    Xg = np.zeros((c.size, c.max()+1))
    Xg[np.arange(c.size), c] = 1
    return Xg


def nll_glm(params, y, Xg, Dg=None, nzmean=None, diversity=None, dist="nb"):
    """ Negative log-likelihood of the ZINB-GLM model """
    eps = 1e-20  # A little epsilon to avoid errors
    pg = params[0]
    ag, yg = np.split(params[1:], 2)
    # Convert to compatible shapes, column vectors
    ag = ag.reshape(-1,1)
    yg = yg.reshape(-1,1)
    n_reads = y[y==0].reshape(-1,1)
    y_reads = y[y>0].reshape(-1,1)
    # Use the formulas to get distribution parameters
    n = np.exp(pg)  # Dispersion parameters
    # Mean is estimated using ag if there is no nzmean arg
    mu = np.exp(Xg.dot(ag))
    # Multiply mean offsets, ag estimates a ratio of the expected mean over non-zero mean
    if nzmean is not None: mu = np.multiply(mu, Dg.dot(nzmean.reshape(-1,1)))
    # Extraneous is estimated by yg
    yg = Xg.dot(yg)
    # Add extraneous offset, yg estimates the change from the average diversity
    if diversity is not None:
        temp = Dg.dot(diversity.reshape(-1,1))
        temp[temp<=0] = eps
        temp[temp>=1] = 1-eps
        yg += np.log(temp/(1-temp))  # logit function, hopefully this is never 0 or 1
    pi = 1/(1+np.exp(-yg))  # Sigmoid, inverse of logit function
    p = n/(mu+n)

    # ZINB (Equation 2)
    if dist=="nb":
        # Negative Binomial
        n_reads = pi[y==0] + (1.0-pi[y==0])*nbinom.pmf(n_reads, n, p[y==0])
        y_reads = (1.0-pi[y>0])*nbinom.pmf(y_reads, n, p[y>0])
    elif dist=="norm":
        # Normal
        # "Normal Approximation to the Negative Binomial" by statisticsmatt
        # https://www.youtube.com/watch?v=JhmmbgLLVkQ
        var = n*(1-p)/(p*p)
        sd = np.sqrt(var)
        n_reads = pi[y==0] + (1.0-pi[y==0])*norm.pdf(n_reads, mu[y==0], sd[y==0])
        y_reads = (1.0-pi[y>0])*norm.pdf(y_reads, mu[y>0], sd[y>0])
    # Compute the negative log-likelihood
    """ https://stackoverflow.com/questions/5124376/convert-nan-value-to-zero/5124409 """ 
    # The stackoverflow answer was ok, but sometimes values are still 0
    # I added a little epsilon and everything seems to work better.
    n_reads = np.nan_to_num(np.log(n_reads + eps))
    y_reads = np.nan_to_num(np.log(y_reads + eps))
    nll = -(np.sum(n_reads) + np.sum(y_reads))
    return nll


def glm_fit(data, Xg, Dg=None, nzmean=None, diversity=None, dist="nb", debug=False):
    """ Fits a single gene to a ZINB distribution.
        There are some hard-code initial values in this function.
    """
    # Scipy function requires integers
    # np.rint converts to the nearest integer
    data = np.rint(np.array(data)).astype(np.int)
    k = Xg.shape[1]  # k conditions, columns in design matrix
    # Pack the parameters into a list, required for scipy.optimize
    # Initialization the mean as 1 assuming that the nzmean arg is supplied
    init_mean = 1
    # If there is no nzmean arg, then use the non-zero mean of the data
    if nzmean is None: init_mean = np.mean(data[data>0])
    pg = np.log(init_mean)  # coefficient for r/n
    ag = np.log(init_mean)  # coefficient for finding mean
    # Initialization extraneous probability
    init_pi = 0.5  # initial extraneous probability as with offsets
    # If there is no diversity arg, then use a low probability
    if diversity is None: init_pi = 0.000001
    yg = np.log(init_pi/(1-init_pi))  # coefficient for pi
    params = np.array([pg] + k*[ag] + k*[yg])
    # Using scipy.minimize throws erros without using bounds
    # These are some upper and lower bounds that I found to not give errors
    eps = 1e-7  # using zero gives warnings, so this is a small epsilon value
    upper_n = np.log(1e12)  # highest dispersion is 1e12, avoids overflows in gamln in the nbinom.pmf
    upper_exp = np.log(np.finfo(float).max) - 10  # largest value that np.exp() can use, minus 10 avoids most overflows
    upper_pi = 1.0-np.finfo(float).eps  # highest probabilty possible, 1-epsilon
    upper_sigmoid = np.log(upper_pi/(1-upper_pi))  # highest value of sigmoid
    #          n>0               pi(sigmoid)                 mu(exp)
    bounds = [(eps, upper_n)] + k*[(None, upper_sigmoid)] + k*[(None, upper_exp)]
    method = "L-BFGS-B"  # Recommended for bounded optimization
    options = {"ftol":1e-9, "gtol": 1e-8, "eps": 1e-8}  # TODO : tried lowering these, no much change, rescaling helped more
    results = minimize(nll_glm, params, args=(data, Xg, Dg, nzmean, diversity, dist), bounds=bounds,
                        method=method, options=options)
    # Unpack the paramters and copmute the estimated distribution parameters
    params = results.x
    pg = params[0]
    ag, yg = np.split(params[1:], 2)
    n = np.exp(pg)
    mu = np.exp(ag)
    if nzmean is not None: mu = np.multiply(Xg.dot(mu), Dg.dot(nzmean)).reshape(len(nzmean),-1)[:,0]
    if diversity is not None: x = (Xg.dot(yg) + Dg.dot(diversity)).reshape(len(diversity),-1)[:,0]
    pi = 1/(1+np.exp(-yg))
    p = n/(mu+n)  # From the paper
    # var = n*(1-p)/(p*p)  # From scipy docs
    if debug:
        print("n :", n); print("p :", p)
        print("diversity :", diversity); print("pi :", pi)
        print("ag :", ag, np.exp(ag))
        print("nzmean :", nzmean); print("mu :", mu)
    return params


def zinb_glm_llr(data, conditions, nzmean=None, diversity=None, dist="nb", rescale=0, debug=False):
    """ Compute the Likelihood ratio test on the ZINB-GLM model """
    conditions1 = np.array(conditions)
    conditions0 = [0]*len(conditions1)
    if debug: print(" data :", np.rint(np.array(data, dtype=np.float)))
    if rescale>0:
        data = (rescale/np.mean(data[data>0])) * data
        if debug: print("scale :", data)
    data = np.rint(np.array(data, dtype=np.float)).astype(np.int)
    # Difference in parameters is the degrees of freedom
    # 1 dispersion, mu and pi for all conditions, minus the 3 for the condition-indepented model
    delta_df = 1 + 2*num_conditions(conditions) - 3
    # Create the observation design matrix
    Dg = None  # Default is no observation matrix, 
    if nzmean is not None or diversity is not None:
        # It is a binary matrix for the number of samples
        # Since this test runs per gene, each gene has the same number of TA
        # site in all conditions and that information along with the number of nzmeans
        # can be used to create an observations matrix.
        x = nzmean if nzmean is not None else diversity
        y = [i for i in range(len(x)) for _ in range(int(len(conditions)/len(x)))]
        Dg = design_matrix(y)
    # Create the design matric for the null and alternative models
    Xg0 = design_matrix(conditions0)
    Xg1 = design_matrix(conditions1)
    params0 = glm_fit(data, Xg0, Dg, nzmean, diversity, dist, debug=debug)  # Null Hypothesis
    params1 = glm_fit(data, Xg1, Dg, nzmean, diversity, dist, debug=debug)  # Alternative Hypothesis
    # Compute the log-likelihood
    l0 = nll_glm(params0, data, Xg0, Dg, nzmean, diversity, dist)
    l1 = nll_glm(params1, data, Xg1, Dg, nzmean, diversity, dist)
    # Compute log-likelihood ratio
    llr = 2*(l0-l1)
    # llr is approximately chi-squared
    # pvalue is the probability that we get this llr or higher for the given degrees of freedom
    pvalue = chi2.pdf(x=llr, df=delta_df)
    if debug:
        print("l0={:e}. l1={:e}. llr={:e}. ".format(l0, l1, llr))
        print("df={}. pvalue={:e}\n".format(delta_df, pvalue))
    return pvalue


def zinb_test(data, conditions, nzmean=None, diversity=None, dist="nb", rescale=0, debug=False):
    """ Perform the ZINB-GLM test with rescaling """
    # Here is the likelihood ratio test (LRT)
    pvalue = zinb_glm_llr(data, conditions, dist, rescale, debug=debug)
    """
    If the test fails, try again with rescaled data.
    These values are arbitrary scaling factors.
    I assume this is an ok operation because when normalizing the
    values are scaled by multiplication, so why not do it here.
    My first thought is magnitude shouldn't matter, but I'm not
    sure since this a
    """
    scale = rescale
    while pvalue in [0, 0.5] and scale < 1000:
        scale += 50
        pvalue = zinb_glm_llr(data, conditions, nzmean, diversity, dist, rescale=scale, debug=debug)
    if debug: print("Scale :", scale)
    return pvalue


def zinb_glm_llr_full_test(data, conditions, nzmean=None, diversity=None,
    dist="nb", rescale=0, debug=False):
    """ Do the Likelihood Ratio Test on the ZINB-GLM for each gene. """
    # Iterate over every row, length of first axis of data
    num_genes = len(data)
    # Temporary array to fill
    pvalues = np.zeros((num_genes, 1))
    # pvalues = np.full((num_genes, 1), np.nan)
    t0 = time.time()  # Start time
    for i in range(num_genes):
        # Helpful time estimate
        if (i+1) % 10 == 0:
            duration = time.time()-t0
            remaining = duration/(i+1) * (num_genes-i+1)
            # TODO : time_to_string
            print("gene {}/{}. {:.1f} genes/second. elapsed={:.2f}. remaining={:.2f}.".format(i+1, num_genes, (i+1)/duration, duration, remaining))

        row_data = data[i]  # Select a single gene data

        # TODO : remove data less than 2 observations
        # Skip if all values are the same, which includes all zero counts
        if np.all(row_data == row_data[0]):
            pvalues[i,:] = np.nan
            continue

        if np.count_nonzero(row_data)<4:
            pvalues[i,:] = np.nan
            continue

        # Here is the likelihood ratio test (LRT)
        pvalue = zinb_glm_llr(row_data, conditions, nzmean, diversity, dist, rescale, debug=debug)

        # Here is the test with rescaling
        # pvalue = zinb_test(row_data, conditions, nzmean, diversity, dist, rescale, debug=debug)

        # Sort the pvalue in the array
        pvalues[i,:] = pvalue

        # leave early for debugging
        if i > 5 and debug:
            break

    return pvalues


if __name__ == "__main__":

    debug = True
    ta_sites = 4
    saturation = 0.5
    num_genes = 2
    min_count = 1
    max_count = 10
    factor = 100
    rescale = 0
    dist = "nb"

    conditions = [0]*ta_sites + [1]*ta_sites
    # conditions = [0, 0, 0, 0, 0, 1, 1, 1]
    # print(design_matrix(conditions))
    size = (num_genes, len(conditions))
    data = np.random.randint(min_count, max_count, size=size)
    data += np.random.randint(min_count, factor*max_count, size=size)
    # Mask based on an approximate array saturation
    mask = np.random.choice([0, 1], size=(num_genes, len(conditions)), p=[1-saturation, saturation])
    data = mask*data
    data = data + data*conditions*1.0
    data = data.astype(np.int)

    # conditions = [0, 0, 0, 0, 1, 1, 1, 1]
    # data = np.array( [
    #     [    0,  6959,  3473,     0,  2859, 23979, 22260,     0],
    #     [    0,     0,  1191,     0, 21981,  3570,  4227,     0],
    #     [    0,     0,  8641,  9077, 13916, 10370,     0,     0],
    # ] )

    nzmean = None
    diversity = None
    split = np.split(data, len(conditions)/ta_sites, axis=1)
    nzmean = np.array([s[s>0].mean() for s in split])
    diversity = np.array([np.count_nonzero(s)/s.size for s in split])

    pvalues = zinb_glm_llr_full_test(data, conditions, nzmean=nzmean, diversity=diversity,
        dist=dist, rescale=rescale, debug=debug)
    pvalues = np.squeeze(pvalues)

    sig_bool = pvalues<0.05

    no_test = np.sum(np.isnan(pvalues))
    print("No Test = {} ({:.2f}%)".format(no_test, 100*no_test/len(pvalues)))

    fails = np.sum(pvalues==0)
    print("Fails = {} ({:.2f}%)".format(fails, 100*fails/len(pvalues)))

    # exit()

    for i in range(len(data)):
        # print(data[i], sig_bool[i], pvalues[i])
        print(i, sig_bool[i], pvalues[i])
        if pvalues[i]==0:
            print(data[i])
            s1, s2 = np.split(data[i], 2)
            print("split    mean :", np.mean(s1), np.mean(s2))
            print("split nz mean :", np.mean(s1[s1>0]), np.mean(s2[s2>0]))
            print("mean={:.2f}. non-zero mean={:.2f}".format(np.mean(data[i]), np.mean(data[i][data[i]>0])))
