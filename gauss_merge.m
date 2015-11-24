%%
%% Usage:
%%
%% [Mu, Sigma, Prior] = gauss_merge (Mu1, Sigma1, Prior1, Mu2, Sigma2, Prior2);
%%
%% This function combines/merges multichannel Gaussian statistics. Currently it 
%% does so according to a standard mixture reduction (SMR) technique reviewed  
%% in Reece & Roberts (2010). 
%% FIXME: implement optional GCU merging, as described in the same paper, so 
%%        that it utilizes two-element Priors as the upper and lower bounds.
%%
%%
%% INPUTS:
%%
%% Mu1          - 1Xnchan mean vector
%%
%% Sigma1       - nchanXnchan covariance matrix
%%
%% Prior1       - scalar weight...both priors will be normalized such that 
%%                they sum to unity, so it is OK to pass unscaled weights,
%%                even observation counts, to this function as Priors.
%%                FIXME: allow Priors to be two-element vectors specifying
%%                       upper and lower bounds for the GCU algorithm
%%
%% Mu2          - 1xnchan mean vector
%%
%% Sigma2       - nchanXnchan covriance matrix
%%
%% Prior2       - scalar weight...both priors will be normalized such that 
%%                they sum to unity, so it is OK to pass unscaled weights,
%%                even observation counts, to this function as Priors.
%%                FIXME: allow Priors to be two-element vectors specifying
%%                       upper and lower bounds for the GCU algorithm
%%
%%
%% OUTPUTS:
%%
%% Mu           - 1Xnchan mean vector
%%
%% Sigma        - nchanXnchan covariance matrix
%%
%% Prior        - simple sum of Prior1 and Prior2.
%%                FIXME: if input Priors were two-element vectors specifying
%%                       upper and lower bounds for the GCU algorithm, the
%%                       output prior should also be a two-element vector, but
%%                       one specifying the GCU-optimized Prior1 and Prior2.
%%

%%
%% REFERENCES:
%%
%% Reece, S., and Roberts, S. (2010), Generalized Covariance Union: a unified
%%   approach to hypothesis merging in tracking, IEEE Trans. Aerospace Elec.
%%   Systems, v. 46, n. 1, p. 207-221.
%%
%%
%% AUTHOR(s):     E. Joshua Rigler (EJR)
%%
%%
%% CHANGELOG:     28 Oct, 2012 - First public release (EJR)
%%

function [Mu, Sigma, Prior, S_adj1, S_adj2] = gauss_merge(Mu1, Sigma1, Prior1, Mu2, Sigma2, Prior2);
   
   % keep pre-processing simple, only proceed if all inputs are provided
   if nargin ~= 6
      error('gauss_merge: six input arguments are required');
   endif
   
   
   % check for valid Mu1 (Mu2 is checked implicitly if it is same size as Mu1)
   if ndims(Mu1) > 2 || ...
      size(Mu1,1) ~= 1
      error('gauss_merge: Mu1 must be a 1Xnchan row vector');
   else
      nchan = size(Mu1,2);
   endif
   
   
   % check for valid Sigma1 (Sigma2 is checked implicitly if it is same size as Sigma2)
   if ndims(Sigma1) > 2 || ...
      size(Sigma1,1) ~= size(Sigma1,2) || ...
      size(Sigma1,2) ~= nchan
      error('gauss_merge: Sigma1 must be a nchanXnchan square matrix');
   endif
   
   
   % check that first Gaussian is sized identically to second Gaussian
   if ~all(size(Mu1) == size(Mu2)) || ...
      ~all(size(Sigma1) == size(Sigma2))
      error('gauss_merge: multivariate Gaussians must have identical dimensions to merge');
   endif
   
   
   % check Priors
   if ~isscalar(Prior1)
      
      if numel(Prior1) == 2
         error('gauss_merge: non-scalar Prior1 not supported yet');
      else
         error('gauss_merge: Prior1 not a valid weighting bounds vector');
      endif
      
   endif
   
   if ~isscalar(Prior2)
      
      if numel(Prior2) == 2
         error('gauss_merge: non-scalar Prior2 not supported yet');
      else
         error('gauss_merge: Prior2 not a valid weighting bounds vector');
      endif
      
   endif
   
   
   %
   % And let the real games begin!!!
   %
   
   
   % Pre-allocate matrices for outputs
   Prior = Prior1 * 0;
   Mu = Mu1 * 0;
   Sigma = Sigma1 * 0;
   
   % Prior is just the sum of Prior1 and Prior2 if using SMR
   Prior = Prior1 + Prior2;
   
   % Mu is weighted average of Mu1 and Mu2 using Prior1 and Prior2 as weights
   Mu = (Mu1 * Prior1 + Mu2 * Prior2 ) / Prior;

   
   % Sigma is weighted average of Sigma1 and Sigma2 after they are adjusted by 
   % deviations of Mu1 and Mu2 from the new Mu
   S_adj1 = Sigma1 + (Mu - Mu1)' * (Mu - Mu1);
   S_adj2 = Sigma2 + (Mu - Mu2)' * (Mu - Mu2);
   Sigma = (S_adj1 * Prior1 + S_adj2 * Prior2) / Prior;

   
endfunction
