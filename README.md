# MBCA - MATLAB software
Multi-Batch Covariance Adjustment:

   [newY,delta_,Qstat,pval] = MBCA(Y,F,batchlabel,depth) performs a batch
   adjustment such that covariances among batches are equalized as well as the mean.

   For example, [newY] = MBCA(Y,F,batchlabel) generates an adjusted data
   set 'newY' substuting a input data set 'Y'.

   Input: 

          1) 'Y' is (p by n) whole data matrix (p > n) with K batches ; 
           e.g., Y = [Y1 Y2] where Y1,Y2 are (p by n1),(p y n2) matrix.
          
          2) 'F' is (q by n) factor matrix (q < n) ; q is the number of
          factors reduced from p dimension.
          Note : If 'F' is empty (F=[]), the MBCA will conduct 'k-means cluster'
          to generate q by n matrix. 
             
          3) 'batchlabel' is a vector of group labels  (ex. 111222333).

          4) 'depth' is a searching depth of sparse level in covariance
          estimates.
             [a](default) 'depth = 2' extends the searching range by two steps from a common delta.
             [b] depth should be integer (e.g., 1,2,3,4... but do not recommend too high number due
             to computational efficiency)
         
   Output: 

           1) 'newY' is a (p by n) adjusted data matrix substituting Y
           
           2) 'delta_' is a (K by 1) vector of sparse level for each
           batch.

           3) 'Qstat' is a test statistic for the equality of HD
           covariance matricess [2] after batch adjustment; Smaller value
           indicates better homogeneity among covariance matrices.                

           4) 'pval' is p-value of the Qstat for the null hypothesis of
           equal variance covariance for Gaussian Data.
        
   Required m-file function: covtestQsvd.m

   References:

   [1] JA Lee, KK dobbin and J Ahn. (2014), Covariance adjustment for batch effect in gene expression data, Statistics in Medicine.   
   [2] Fan et al. (2008), High dimensional covariance matrix estimation using a factor model, Journal of Econometrics.
   [3] Srivastava and Yanagihara (2010), Testing the equality of several covariance matrices with fewer observations than the dimension, Journal
   of Multivariate Analysis.
