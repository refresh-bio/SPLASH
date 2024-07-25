'''
Implementation of OASIS statistical test(https://www.pnas.org/doi/10.1073/pnas.2304671121)
Tavor Z. Baharav, David Tse, and Julia Salzman

OASIS (Optimized Adaptive Statistic for Inferring Structure) utilizes a linear test-statistic, enabling the computation of closed form P-value bounds, exact asymptotic ones, and interpretable rejection of the null. It is implemented and used in SPLASH: https://github.com/refresh-bio/SPLASH.

The key method in this file is OASIS_pvalue (at the bottom of this file), which returns the p-value of the OASIS test on a given contingency table. The method can be used to compute the p-value using either the asymptotic or finite sample p-value (asymptotic flag), and can return the optimizing row and column embeddings (return_f_c flag).

Please reach out on github with any questions, or directly to tavorb@mit.edu.
'''

import numpy as np 
import scipy.stats



###### 1: utility functions

### split contingency table into train and test data
def splitCounts(mat,downSampleFrac = .5): #slight modification of https://stackoverflow.com/questions/11818215/subsample-a-matrix-python
    keys, counts = zip(*[
    ((i,j), mat[i,j])
        for i in range(mat.shape[0])
        for j in range(mat.shape[1])
        if mat[i,j] > 0
    ])
    # Make the cumulative counts array
    counts = np.array(counts, dtype=np.int64)
    sum_counts = np.cumsum(counts)

    # Decide how many counts to include in the sample
    frac_select = downSampleFrac
    count_select = int(sum_counts[-1] * frac_select)

    # Choose unique counts
    ind_select = sorted(np.random.choice(range(sum_counts[-1]), count_select,replace=False))

    # A vector to hold the new counts
    out_counts = np.zeros(counts.shape, dtype=np.int64)

    # Perform basically the merge step of merge-sort, finding where
    # the counts land in the cumulative array
    i = 0
    j = 0
    while i<len(sum_counts) and j<len(ind_select):
        if ind_select[j] < sum_counts[i]:
            j += 1
            out_counts[i] += 1
        else:
            i += 1

    # Rebuild the matrix using the `keys` list from before
    out_mat = np.zeros(mat.shape, dtype=np.int64)
    for i in range(len(out_counts)):
        out_mat[keys[i]] = out_counts[i]
        
    return out_mat

### simple wrapper of the above to downsample matrix columnwise, as opposed to overall
def splitCountsColwise(mat,downSampleFrac = .5): #slight modification of https://stackoverflow.com/questions/11818215/subsample-a-matrix-python
    out_mat = np.zeros_like(mat)
    for j in range(mat.shape[1]):
        if mat[:,j].sum()>0:
            out_mat[:,j] = splitCounts(np.reshape(mat[:,j],(mat.shape[0],-1)),downSampleFrac).flatten()
        
    return out_mat

### shift and scale an input vector to be in the range [minval, maxval]
def normalizevec(x,minval=0,maxval=1):
    if x.max()==x.min():
        return np.zeros_like(x)
    x = np.array(x)
    x01= (x-x.min())/(x.max()-x.min())
    return x01*(maxval-minval)+minval




######### 2: c,f generation (row and column embeddings for contingency table)

### starting at c, run alternating maximization
def altMaximize(X,c):
    #### if clustering put all in same cluster, perturb
    if np.all(c==c[0]):
        c[0] = -1*c[0]
        
    nj = X.sum(axis=0)
    njinvSqrt = 1.0/np.maximum(1,np.sqrt(nj)) ### avoid divide by 0 errors
    njinvSqrt[nj==0]=0
    
    Xtild = (X - 1.0/X.sum()*np.outer(X@np.ones(X.shape[1]), X.T@np.ones(X.shape[0]))) @ np.diag(njinvSqrt)
    
    Sold = 0
    i=0
    while True:
        ### find optimal f for fixed c
        f = np.sign(Xtild @ c)
        f1 = (f+1)/2 ### to rescale f to be [0,1] valued
        f2 = (1-f)/2
        f = f1
        if np.abs(f2@Xtild@c) > np.abs(f1@Xtild@c):
            f = f2
        
        ### find optimal c for fixed f
        c = Xtild.T @ f
        if np.linalg.norm(c)>0:
            c /= np.linalg.norm(c)
        
        ### compute objective value, if fixed, stop
        S = f @ Xtild @ c
        if S==Sold: ### will terminate once fOpt is fixed over 2 iterations
            break
        Sold = S
        i+=1
        if i>50:
            c = np.zeros_like(c)
            f=np.zeros_like(f)
            S=0
            break
    return c,f,np.abs(S)



### find locally-optimal (unconstrained) c and f from random initialization
### optimizes finite sample p-value bound
def generate_cf_finite_optimized(X, randSeed=0,numRandInits=10):
    np.random.seed(randSeed) ### random initialization and extension
    
    relevantTargs = X.sum(axis=1)>0
    relevantSamples = X.sum(axis=0)>0

    nrows,ncols=X.shape
    
    if relevantTargs.sum()<2 or relevantSamples.sum()<2:
        return np.zeros(ncols),np.zeros(nrows)

    X = X[np.ix_(relevantTargs,relevantSamples)]
    
    Sbase=0
    fMax=0
    cMax=0
    for _ in range(numRandInits):
        c = np.random.choice([-1,1],size=X.shape[1])
        c,f,S = altMaximize(X,c)
        if S > Sbase:
            fMax = f
            cMax = c
            Sbase = S
            

    ## extend to targets and samples that didn't occur previously
    fElong = np.random.choice([0,1],size=nrows)
    fElong[relevantTargs] = fMax
    fOpt = fElong
    
    cElong = np.zeros(ncols)
    cElong[np.arange(ncols)[relevantSamples]]=cMax ### fancy indexing
    cOpt = cElong
    
    return cOpt,fOpt



######## optimal col,row embeddings c,f for asymptotic p-value
def generate_cf_asymp_optimized(X):
    c = np.ones(X.shape[1])
    
    zeroIdxs = X.sum(axis=1)==0
    X = X[~zeroIdxs]

    zeroCols = X.sum(axis=0)==0
    X = X[:,~zeroCols]

    ### empirical probability dist p over rows of X
    p = X.sum(axis=1)/X.sum()
    Xtild = (X-np.outer(X.sum(axis=1),X.sum(axis=0))/X.sum())@np.diag(1.0/np.sqrt(X.sum(axis=0)))    

    A = np.diag(1.0/p)@Xtild @ Xtild.T

    ### set f to be principal eigenvector of A
    eigvals,eigvecs = np.linalg.eig(A)
    f = eigvecs[:,np.argmax(eigvals)]

    ## retain only real part of f
    f = np.real(f)
    
    c = Xtild.T@f
    c/= np.linalg.norm(c)

    fOld = f.copy()
    ### map the entries of fOld to f, with the zeroIdxs zeroed out
    f = np.zeros(len(zeroIdxs))
    f[~zeroIdxs] = fOld

    cOld = c.copy()
    c = np.zeros(len(zeroCols))
    c[~zeroCols] = cOld

    return c,f



######## 3. p-value computation

### asymptotic pvalue
### OASIS statistic is asymptotically normal with variance totalVar under the null
###   allowing is to provide asymptotically valid p-values
def testPval_asymp(X,cOpt,fOpt):
    if (cOpt==0).all():
        return 1
    
    cOpt = np.nan_to_num(cOpt,0)
    fOpt = np.nan_to_num(fOpt,0)


    nj = X.sum(axis=0)
    njinvSqrt = 1.0/np.maximum(1,np.sqrt(nj))
    njinvSqrt[nj==0]=0
    
    ### compute p value
    S = fOpt @ (X-X@np.outer(np.ones(X.shape[1]),nj)/X.sum())@(cOpt*njinvSqrt)

    M=X.sum()
    
    muhat = (fOpt@X).sum()/M
    
    varF = (fOpt-muhat)**2 @ X.sum(axis=1)/X.sum()
    totalVar = varF * (np.linalg.norm(cOpt)**2 - (cOpt@np.sqrt(nj))**2/M)
    
    if totalVar<=0:
        return 1 ## numerical error / issue
    
    normalizedTestStat = S/np.sqrt(totalVar)
    pval = 2*scipy.stats.norm.cdf(-np.abs(normalizedTestStat))
                
    return pval


### finite-sample valid p-value bound for fixed c and f on contingency table
def testPval_finite(X,cOpt,fOpt):
    
    cOpt = np.nan_to_num(cOpt,0)
    fOpt = np.nan_to_num(fOpt,0)

    if np.all(cOpt == cOpt[0]) or np.all(fOpt==fOpt[0]):
        return 1

    if fOpt.max()-fOpt.min()>1:
        fOpt /= (fOpt.max()-fOpt.min())    

    assert(fOpt.max()-fOpt.min()<=1)
        
    nj = X.sum(axis=0)
    njinvSqrt = 1.0/np.maximum(1,np.sqrt(nj))
    njinvSqrt[nj==0]=0
    
    ### compute test statistic
    S = fOpt @ (X-X@np.outer(np.ones(X.shape[1]),nj)/X.sum())@(cOpt*njinvSqrt)

    M=X.sum()
    
    denom = (np.linalg.norm(cOpt)**2 - (cOpt@np.sqrt(nj))**2/M)
    pval = 2*np.exp(-2*S**2/denom)

    return min(np.nan_to_num(pval,1),1)

### compute chi2 p-value for contingency table
def computeChi2Test(X):
    if len(X.shape)==1:
        return 1
    X = X[X.sum(axis=1)>0]
    X = X[:,X.sum(axis=0)>0]
    _,pv,_,_= scipy.stats.contingency.chi2_contingency(X)
    return pv


### binary effect size from OASIS paper
def effectSize_bin(X,c,f):
    if (c>0).sum()==0 or (c<0).sum()==0:
        return 0

    return np.abs(f@X@(c>0) / (X@(c>0)).sum() - f@X@(c<0) / (X@(c<0)).sum())



def OASIS_pvalue(X, numSplits=5, trainFrac=.25, asymptotic=False, return_f_c=False):
    """
    Computes the p-value using the OASIS method.

    Parameters:
    X (numpy.ndarray): The input count matrix, IxJ.
    numSplits (int, optional): The number of train/test splits for computing optimized f,c on. Default is 5
    trainFrac (float, optional): The fraction of data to be used for training. Default is 0.25.
    asymptotic (bool, optional): Whether to use the asymptotic p-value. Default is False (uses finite sample p-value).
    return_f_c (bool, optional): Whether to return the row and column embeddings. Default is False.

    Returns:
    float: The minimum p-value computed across all splits, multiplied by the number of splits (Bonferroni correction). 
           The returned value is capped at 1.
    f (numpy.ndarray): The I-dimensional row embedding vector (only returned if return_f_c is True).
    c (numpy.ndarray): The J-dimensional column embedding vector (only returned if return_f_c is True).
    """
    I,J = X.shape
    min_pval = 1
    cOpt,fOpt = (np.zeros(J),np.zeros(I))
    f_c_gen_method = generate_cf_asymp_optimized if asymptotic else generate_cf_finite_optimized
    pval_test_method = testPval_asymp if asymptotic else testPval_finite
    
    for i in range(numSplits):
        np.random.seed(i)
        Xtrain = splitCountsColwise(X,trainFrac)
        c,f = f_c_gen_method(Xtrain)
        pval = pval_test_method(X-Xtrain,c,f)
        if pval < min_pval:
            min_pval = pval
            cOpt = c
            fOpt = f
    
    if return_f_c:
        return min(1,numSplits*min_pval), fOpt, cOpt

    return min(1,numSplits*min_pval)