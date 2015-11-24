%%
%% Usage:
%%
%% [ClusterID, ClusterProb] = gmix_cluster (I_chan, Mu, Sigma[, Priors[, pvalue[, impute]]])
%%
%%
%% This function determines which component of a mutivariate Gaussian mixture 
%% distribution each observation (row) in I_chan is most likely to belong to, 
%% and if asked, returns the probability of membership to each cluster. It is
%% little more than a wrapper to gmix_pdf().
%% 
%%
%% INPUTS:
%%
%% I_chan       - multichannel column vector (i.e., rows==observations,
%%                columns==channels). Missing channels are allowed in each
%%                observation, and designated using NaNs as fill values
%%                (see documentation for 'impute' argument below).
%%
%% Mu           - mean vectors for GMM...each mean is a 1Xnchan row vector,
%%                and GMM components are stored along the 3rd dimension.
%%
%% Sigma        - (optional) covariance matrices for GMM; each covariance is
%%                a nchanXnchan matrix, and GMM components are stored along
%%                the 3rd dimension...if not passed, assume uniform circular
%%                GMMs (i.e., common variance, diagonal matrices) with unit
%%                variance (the magnitude of the variance only matters for
%%                calculating probabilities).
%%
%% Priors       - (optional) prior probabilities of GMM component membership;
%%                each prior, or weight, is a 1X1 scalar, and GMM components
%%                are stored along the 3rd dimension...if not passed, or if
%%                Priors is an empty matrix, do not apply any priors to the
%%                probabilities (or rather set tham all to 1).
%%
%% pvalue       - (optional) specifies the area under the Chi-squared curve
%%                that is considered valid. If not set, consider all calculated
%%                probabilities, no matter how small/insignficiant. If it is
%%                specified, it is used to determine a threshold Mahalanobis
%%                distance (MD) for each observation for each component. If
%%                MD_max is exceeded, the probability is simply set to zero.
%%                If pvalue is a scalar, apply same threshold for each GMM
%%                component; if it is a vector equal in length to the number
%%                of components, apply each element to each component. The
%%                pvalue must be between 0 and 1.
%%
%% impute       - (optional) if not set or false, treat any observation (i.e,
%%                any row of I_chan) that has even a single NaN as bad; this
%%                sets ClusterID and ClusterProb all to zero. If true, impute
%%                missing values using available values and GMM parameters.
%%                If all values in a single row of I_chan are NaN, there is
%%                not a reasonable default value, so set ClusterID and all
%%                of ClusterProb to zero.
%%
%%
%% OUTPUTS:
%%
%% ClusterID    - GMM component each observation belongs to...this integer
%%                falls in the range 1:size(Mu|Sigma|Priors,3), and is an
%%                index into these multi-component statistical parameters.
%%
%% ClusterProb  - probability of belonging to each component of the GMM
%%                distribution...this is a nobsX1Xncomp matrix in order 
%%                to avoid confusion with nobsxnchan multichannel I_chan,
%%                and to be consistent with GMM parameters Mu, Sigma, and
%%                Priors.
%%

%%
%% REFERENCES:
%%
%% McLachlan, G., and Peel, D. (2000), Finite Mixture Models, Wiley and Sons.
%% Calinon, S. (2009), Robot Programming by Demonstration: A Probabilistic
%%   Approach, CRC Press, Taylor & Francis Group, Boca Raton, FL.
%% Fraley, C., and Raftery, A. E. (2002), Model-based clustering, discriminant
%%   analysis, and density estimation, J. American Stat. Assoc., v. 97.458,
%%   p. 611-631. (see also many references cited within).
%% Hunt, L., and Jorgensen, M. (2003), Mixture model clustering for mixed data
%%   with missing information, Comp. Stats. and Data Analysis. v. 41, pp. 429-440.
%%
%%
%% AUTHOR(s):     E. Joshua Rigler (EJR)
%%
%%
%% CHANGELOG:     18 Oct, 2012 - First public release (EJR)
%%                02 Nov, 2012 - splitting the clustering and probability
%%                               calculation into two functions; this function
%%                               is now mostly a wrapper for gmix_pdf, and
%%                               does little more than process inputs and 
%%                               find the max probability from gmix_pdf.
%%
%%


function [ClusterID, ClusterProb] = gmix_cluster (I_chan, Mu, Sigma, Priors, pvalue, impute)
   
   % process inputs
   
   % check minimal arguments
   if nargin < 2
      error('gmix_cluster: at least two input arguments are required');
   endif
   
   
   % if at least two arguments
   if nargin > 1
       
       % basic check of I_chan dimensions
      if ndims(I_chan) > 2
         error('gmix_cluster: I_chan dimensions are not valid');
      endif
      
      nobs = size(I_chan,1);
      nchan = size(I_chan,2);
      
      
      % basic check of Mu dimension
      if size(Mu,1) ~= 1 || ...
         size(Mu,2) ~= nchan
         error('gmix_cluster: I_chan and Mu are not dimensionally consistent');
      endif
      
      ncomp = size(Mu,3);
      
      
      % generate default covariance matrix if Sigma doesn't exist
      if ~exist('Sigma', 'var')
         Sigma = repmat(eye(nchan), [1,1,ncomp]);
      endif
      
      
      % generate default priors
      if ~exist('Priors', 'var')
         Priors = ones(1,1,ncomp);
      endif
      
      
      % generate default pvalue
      if ~exist('pvalue', 'var')
         pvalue = ones(ncomp, 1);
      endif
      
      
      % generate default impute
      if ~exist('impute', 'var')
         impute = false;
      endif

   endif
   
   
   % if at least three arguments
   if nargin > 2
      % pass input directly to gmix_pdf
   endif
   
   
   % if at least four arguments
   if nargin > 3
      % pass input directly to gmix_pdf
   endif
   
   
   % if at least five arguments
   if nargin > 4
      % pass input directly to gmix_pdf
   endif
   
   
   % if at least six arguments
   if nargin > 5
      % pass input directly to gmix_pdf
   endif
   
   
   
   %
   % begin algorithm
   %
   
   % call gmix_pdf to do the heavy lifting
   [~, ClusterProb] = gmix_pdf(I_chan, Mu, Sigma, Priors, pvalue, impute);
   
   
   % The assigned cluster will be the index along the third dimension in
   % ClusterProbs that corresponds to the maximum probability density
   [P_tmp, ClusterID] = max(ClusterProb, [], 3);
   
   
   % if the maximum probability is zero, all probabilities must have been
   % zero, which causes max() to return an index of 1...not desirable; 
   % if maximum probabilty is a NaN, then all probabilities must have
   % been a NaN, which indicates that at least one value in the feature
   % vector was a NaN, and therefore unclassifiable unless imputation
   % was used...if imputation *was* used, and maximum probability is a
   % NaN, then there must not have been any valid values to begin with.
   ClusterID(P_tmp == 0 | isnan(P_tmp)) = 0;
   ClusterProb(P_tmp == 0 | isnan(P_tmp), 1, :) = 0;
   
   
   % max() does not return integer-valued indices
   ClusterID = int16(ClusterID);
   
      
endfunction




