%%
%% Usage:
%%
%% [Mu, Sigma, Prior, dist] = gmix_merge (Mu1, Sigma1, Prior1, Mu2, Sigma2, Prior2, metric);
%%
%% This function combines/merges multichannel Gaussian Mixture Model (GMM)
%% statistics. It does so by finding the component of the GMM described by
%% Mu1 and Sigma1 that is "nearest" to each element of the GMM described by
%% MU2 and Sigma2, then merging the two according to their parameters and
%% respective weights. This means that the outputs always match Mu1's and 
%% Sigma1's dimensions. 
%% 
%% NOTE: once a nearest neighbor is found for GMM2 inside GMM1, it is NOT
%%       removed from the pool of GMM1's components because to do so would
%%       imply that the components of the GMMs were ordered in some way, 
%%       which they almost certainly are not. Therefore, it is entirely
%%       possible that multiple GMM2 components get merged with the same
%%       GMM1 components. All comparisons are made using the components
%%       passed into the gmix_merge function, NOT the merged components.
%%
%% FIXME: the small amount of code from gauss_merge should be included in this
%%         function, and gauss_merge should be re-implemented as a wrapper to
%%         this function, NOT the other way round as is currently the case.
%%
%%
%% INPUTS:
%%
%% Mu1          - mean vectors for GMM...each mean is a 1Xnchan row vector,
%%                and GMM components are stored along the 3rd dimension.
%%
%% Sigma1       - covariance matrices for GMM; each covariance is
%%                a nchanXnchan matrix, and GMM components are stored along
%%                the 3rd dimension.
%%
%% Prior1       - prior probabilities of GMM component membership;
%%                each prior, or weight, is a 1X1 scalar, and GMM components
%%                are stored along the 3rd dimension.
%%                FIXME: allow Priors to be two-element vectors specifying
%%                       upper and lower bounds for the GCU algorithm...
%%
%% Mu2          - mean vectors for GMM...each mean is a 1Xnchan row vector,
%%                and GMM components are stored along the 3rd dimension.
%%
%% Sigma2       - covariance matrices for GMM; each covariance is
%%                a nchanXnchan matrix, and GMM components are stored along
%%                the 3rd dimension.
%%
%% Prior2       - prior probabilities of GMM component membership;
%%                each prior, or weight, is a 1X1 scalar, and GMM components
%%                are stored along the 3rd dimension.
%%                FIXME: allow Priors to be two-element vectors specifying
%%                       upper and lower bounds for the GCU algorithm...
%%
%% metric       - (optional) "metric" to calculate divergence between GMM
%%                components...the default value is is 1; see gauss_dist docs
%%                for other options.
%%
%%
%% OUTPUTS:
%%
%% Mu           - mean vectors for GMM...each mean is a 1Xnchan row vector,
%%                and GMM components are stored along the 3rd dimension. Mu's
%%                dimensions will always match Mu1.
%%
%% Sigma        - covariance matrices for GMM; each covariance is
%%                a nchanXnchan matrix, and GMM components are stored along
%%                the 3rd dimension. Sigma's dimensions will always match
%%                Sigma1.
%%
%% Prior        - prior probabilities of GMM component membership;
%%                each prior, or weight, is a 1X1 scalar, and GMM components
%%                are stored along the 3rd dimension. Prior's dimensions will
%%                always match Prior1.
%%                FIXME: if input Priors were two-element vectors specifying
%%                       upper and lower bounds for the GCU algorithm, the
%%                       output prior should also be a two-element vector, but
%%                       one specifying the GCU-optimized Prior1 and Prior2;
%%
%% dist         - returns the matrix of distances between each component in
%%                GMM2 with each component in GMM1...only useful for dignostics.
%%

%%
%% REFERENCES:
%%
%% Reece, S., and Roberts, S. (2010), Generalized Covariance Union: a unified
%%   approach to hypothesis merging in tracking, IEEE Trans. Aerospace Elec.
%%   Systems, v. 46, n. 1, p. 207-221.
%%
%%
%%
%% AUTHOR(s):     E. Joshua Rigler (EJR)
%%
%%
%% CHANGELOG:     28 Oct, 2012 - First public release (EJR)
%%                19 Dec, 2012 - Moved distance calculations to a separate
%%                               function - gauss_distance.m

function [Mu, Sigma, Prior, dist] = gmix_merge (Mu1, Sigma1, Prior1, Mu2, Sigma2, Prior2, metric);
   
   % check inputs
   if nargin < 6
      error('gmix_merge: at least six input arguments are required');
   endif
   
   
   % if only 6 inputs
   if nargin == 6
      % default "metric"
      metric = 1;
   endif
   
   
   
   % check for valid Mu1 (Mu2 checked implicitly if first 2 dims are same size as Mu1)
   if ndims(Mu1) > 3 || ...
      size(Mu1,1) ~= 1
      error('gmix_merge: Mu1 must be a 1XnchanXncomp row vector');
   else
      nchan = size(Mu1,2);
   endif
   
   
   % check for valid Sigma1 (Sigma2 checked implicitly if first 2 dims are same size as Sigma2)
   if ndims(Sigma1) > 3 || ...
      size(Sigma1,1) ~= size(Sigma1,2) || ...
      size(Sigma1,2) ~= nchan
      error('gmix_merge: Sigma1 must be a nchanXnchan square matrix');
   endif
   
   
   % check that first GMM is sized identically to second GMM (dims 1 & 2 only)
   if ~(size(Mu1,1) == size(Mu2,1)) || ... % not sure this is necessary
      ~(size(Mu1,2) == size(Mu2,2)) || ...
      ~(size(Sigma1,1) == size(Sigma2,1)) || ...
      ~(size(Sigma1,2) == size(Sigma2,2))
      error('gmix_merge: multivariate Gaussians must have identical dimensions to merge');
   endif
   
   ncomp1 = size(Mu1,3);
   ncomp2 = size(Mu2,3);
   
   
   % check Priors
   if size(Prior1,3) ~= ncomp1 
      error('gmix_merge: Prior1 is not a valid weights array');
   endif
    if size(Prior2,3) ~= ncomp2 
      error('gmix_merge: Prior2 is not a valid weights array');
   endif
  
   % placeholder for eventually allowing u/l bounds to be passed for GCU
   if size(Prior1,1) > 1 || size(Prior1,2) > 1
      
      if size(Prior1,1) == 2 || size(Prior1,2) == 2
         error('gmix_merge: non-scalar Prior1 not supported yet');
      else
         error('gmix_merge: Prior1 not a valid weighting bounds vector');
      endif
      
   endif
   
   if size(Prior2,1) > 1 || size(Prior2,2) > 1
      
      if size(Prior2,1) == 2 || size(Prior2,2) == 2
         error('gmix_merge: non-scalar Prior2 not supported yet');
      else
         error('gmix_merge: Prior2 not a valid weighting bounds vector');
      endif
      
   endif
   
   
   
   %
   % begin actual algorithm
   %
  
   % calculate distance from GMM2 components to GMM1 components
   for k2=1:ncomp2
      for k1=1:ncomp1
         
         % calculate distance from comp2(k2) to comp1(k1)
         
         dist(k1,k2) = gauss_distance(Mu1(1,:,k1), Sigma1(:,:,k1), ...
                                      Mu2(1,:,k1), Sigma2(:,:,k1), metric);
                  
      endfor % for k1=1:ncomp1
   endfor % for k2=1:ncomp2
   
      
   
   % merge GMM2 component with nearest GMM1 component outside previous loop
   % because we wish the distances to be wrt the original distributions, NOT
   % modified distributions. ...is this really what we want? -EJR 10/2012
   Mu = Mu1;
   Sigma = Sigma1;
   Prior = Prior1;
   for k2=1:ncomp2
      
      [~, k1] = min(dist(:,k2));
      [Mu_tmp, Sigma_tmp, Prior_tmp] = gauss_merge(Mu(1,:,k1),  Sigma(:,:,k1),  Prior(:,:,k1), ...
                                                   Mu2(1,:,k2), Sigma2(:,:,k2), Prior2(:,:,k2));
      
      % copy temporary values into Mu
      Mu(1,:,k1) = Mu_tmp;
      Sigma(:,:,k1) = Sigma_tmp;
      Prior(:,:,k1) = Prior_tmp;
      
   endfor % for k2=1:ncomp2
   
   
endfunction
