"""
Zero-Inflated Binomial Regression [Subramaniyam et al., 2019]

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


def nll_glm(params, y, conditions, dist="nb"):
    """ Negative log-likelihood of the ZINB-GLM model """
    pg = params[0]
    ag, yg = np.split(params[1:], 2)
    Xg = design_matrix(conditions)
    # Convert to compatible shapes, column vectors
    ag = ag.reshape(-1,1)
    yg = yg.reshape(-1,1)
    n_reads = y[y==0].reshape(-1,1)
    y_reads = y[y>0].reshape(-1,1)
    # Use the formulas to get distribution parameters
    n = np.exp(pg)
    mu = np.exp(Xg.dot(ag))
    pi = Xg.dot( np.exp(yg)/(1+np.exp(yg)) )  # This is another form of sigmoid
    # pi = Xg.dot( 1/(1+np.exp(-yg)) )  # Sigmoid, equivalent to above
    p = n/(mu+n)
    # print("n={:.60f}".format(n))
    # print("p, mu, ag :", p, mu, ag)
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
    eps = 1e-25
    n_reads = np.nan_to_num(np.log(n_reads + eps))
    y_reads = np.nan_to_num(np.log(y_reads + eps))
    nll = -(np.sum(n_reads) + np.sum(y_reads))
    # print("{:.50f}".format(nll))
    return nll


def glm_fit(data, conditions, dist="nb", debug=False):
    """ Fits a single gene to a ZINB distribution.
        There are some hard-code initial values in this function.
    """
    # Scipy function requires integers
    # rint converts to hte nearest integer
    data = np.rint(np.array(data)).astype(np.int)
    k = num_conditions(conditions)  # k conditions
    # Pack the parameters into a list, required for scipy.optimize
    # Initialization the mean and dispersion as the non-zero mean
    nz_mean = np.mean(data[data>0])
    pg = np.log(nz_mean)  # coefficient for r/n
    ag = np.log(nz_mean)  # coefficient for finding mean
    init_pi = 0.001  # extraneous probability
    yg = np.log(init_pi/(1-init_pi))  # coefficient for pi
    params = np.array([pg] + k*[ag] + k*[yg])
    eps = 1e-7  # using zero gives warnings, so this is a small epsilon value
    upper_n = np.log(1e12)  # highest dispersion is 1e12, avoids overflows in gamln in the nbinom.pmf
    upper_exp = np.log(np.finfo(float).max) - 10  # largest value that np.exp() can use, minus 10 avoids most overflows
    y = 1.0-np.finfo(float).eps  # highest probabilty possible
    upper_sigmoid = np.log(y/(1-y))  # highest value of sigmoid
    #          n>0               pi(sigmoid)                 mu(exp)
    bounds = [(eps, upper_n)] + k*[(None, upper_sigmoid)] + k*[(None, upper_exp)]
    results = minimize(nll_glm, params, args=(data, conditions, dist), bounds=bounds)
    # Unpack the paramters and copmute the estimated distribution parameters
    params = results.x
    pg = params[0]
    ag, yg = np.split(params[1:], 2)
    n = np.exp(pg)
    mu = np.exp(ag)
    pi = np.exp(yg)/(1+np.exp(yg))
    p = n/(mu+n)  # From the paper
    var = n*(1-p)/(p*p)  # From scipy docs
    if debug:
        print("pi :", pi); print("n :", n)
        print("p :", p); print("mu :", mu); print("var :", var)
    return params


def zinb_glm_llr_test(data, conditions, dist, debug=False):
    """ Compute the Likelihood ratio test on the ZINB-GLM model """
    conditions1 = np.array(conditions)
    conditions0 = [0]*len(conditions1)
    if debug: print("data :", data); print("cond :", conditions)
    # Difference in parameters is the degrees of freedom
    # 1 dispersion, mu and pi for all conditions, minus the 3 for the condition-indepented model
    delta_df = 1 + 2*num_conditions(conditions) - 3
    params0 = glm_fit(data, conditions0, dist, debug=debug)  # Null Hypothesis
    params1 = glm_fit(data, conditions1, dist, debug=debug)  # Alternative Hypothesis
    # Compute the log-likelihood
    l0 = nll_glm(params0, data, conditions0, dist)
    l1 = nll_glm(params1, data, conditions1, dist)
    # Compute log-likelihood ratio
    llr = 2*(l0-l1)
    # llr is approximately chi-squared
    # pvalue is the probability that we get this llr or higher for the given degrees of freedom
    pvalue = chi2.pdf(x=llr, df=delta_df)
    if debug:
        print("mean={:.2f}. non-zero mean={:.2f}".format(data.mean(), data[data>0].mean()))
        print("l0={:.3f}. l1={:.3f}. llr={:.3f}. ".format(l0, l1, llr))
        print("df={}. pvalue={:.8f}\n".format(delta_df, pvalue))
    return pvalue


def zinb_glm_llr_full_test(data, conditions, debug=False):
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
        
        if debug:
            print("mean :", np.mean(data))
            print("nzmean :", np.mean(data[data>0]))

        # Here is the likelihood ratio test (LRT)
        pvalue = zinb_glm_llr_test(row_data, conditions, dist="nb", debug=debug)
        pvalue = zinb_glm_llr_test(row_data, conditions, dist="norm", debug=debug)

        # TODO : check for errors in llr and pvalue

        # Sort the pvalue in the array
        pvalues[i,:] = pvalue

        # leave early for debugging
        if i > 100 and debug:
            break

    return pvalues


if __name__ == "__main__":

    debug = True
    ta_sites = 2
    saturation = 1.0
    conditions = [0]*ta_sites + [1]*ta_sites
    # print(design_matrix(conditions))

    num_genes = 2
    min_count = 0
    max_count = 10
    factor = 100
    sig = 0.05

    data = np.random.randint(min_count, max_count, size=(num_genes, len(conditions)))
    data += np.random.randint(min_count, factor*max_count, size=(num_genes, len(conditions)))
    # Mask based on an approximate array saturation
    mask = np.random.choice([0, 1], size=(num_genes, len(conditions)), p=[1-saturation, saturation])
    data = mask*data
    data = data + data*conditions*1.0

    # conditions = [0, 0, 1, 1]
    # data = np.array( [
    #     # [11400, 11200],
    #     # [142100, 142200],
    #     # [102400, 108200],
    #     # [0, 0],
    #     [200, 250, 300, 350],
    #     [2000, 2500, 3000, 3500],
    #     [20000, 25000, 30000, 35000],
    # ] )

    pvalues = zinb_glm_llr_full_test(data, conditions, debug=debug)
    pvalues = np.squeeze(pvalues)

    sig_bool = pvalues<sig

    for i in range(len(data)):
        print(data[i], sig_bool[i], pvalues[i])

    print("Genes :", num_genes)
    sig_genes = np.sum(sig_bool!=(pvalues==0))
    print("Significant genes : {} ({:.2f}%)".format(sig_genes, 100*sig_genes/num_genes))
    print("Nan pvalues : {}".format(np.sum(np.isnan(pvalues))))
    print("Zero pvalues : {}".format(np.sum(pvalues==0)))


    mask = np.isfinite(pvalues)  # update these qvals
    qvalues = np.full(pvalues.shape, np.nan)

    import statsmodels.stats.multitest
    rejected, qvalues = statsmodels.stats.multitest.fdrcorrection(pvalues, alpha=0.05)

    print("pvalues :", pvalues)
    print("sig p:", sig_bool)
    print("qvalues :", qvalues)
    print("sig q:", rejected)

