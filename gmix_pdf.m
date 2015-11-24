%%
%% Usage:
%%
%% [gmmProb, compProbs] = gmix_pdf (I_chan, Mu, Sigma[, Priors[, pvalue[, impute]]])
%%
%%
%% This function calculates the probability density of a Gaussian Mixture Model 
%% (GMM) given multichannel inputs and the GMM parameters. If requested, it 
%% returns the probability densities for the individual components of the GMM.
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
%%                GMMs (i.e., common variance, diagonal matrices) with unity
%%                variance (the magnitude of the variance only matters for
%%                calculating probabilities, and the user would probably not
%%                be interested in probabilities if they didn't pass Sigma).
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
%%                sets GmixProb and CompProbs all to zero. If true, impute
%%                missing values using available values and GMM parameters.
%%                If all values in a single row of I_chan are NaN, there is
%%                no reasonable default value, so set all probabilities to zero.
%%
%%
%% OUTPUTS:
%%
%% gmmProb      - total (cumulative) probability across all components.
%%
%% compProbs    - probability of belonging to each component of the GMM
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
%% CHANGELOG:     02 Nov, 2012 - First public release (EJR)
%%                28 Feb, 2013 - replaced calls to det() and inv() with
%%                               pseudo-determinant and pseudo-inverse, using
%%                               single-value decomposition
%%


function [gmmProb, compProbs] = gmix_pdf (I_chan, Mu, Sigma, Priors, pvalue, impute)
   
   % process inputs
   
   % check minimal arguments
   if nargin < 2
      error('gmix_pdf: at least two input arguments are required');
   endif
   
   
   % if at least two arguments
   if nargin > 1
      
      % basic check of I_chan dimensions
      if ndims(I_chan) > 2
         error('gmix_pdf: I_chan dimensions are not valid');
      endif
      
      nobs = size(I_chan,1);
      nchan = size(I_chan,2);
      
      
      % basic check of Mu dimension
      if size(Mu,1) ~= 1 || ...
         size(Mu,2) ~= nchan
         error('gmix_pdf: I_chan and Mu are not dimensionally consistent');
      endif
      
      ncomp = size(Mu,3);
      
      
      % generate default covariance matrix if Sigma doesn't exist
      if ~exist('Sigma', 'var') || ...
         isempty(Sigma)
         Sigma = repmat(eye(nchan), [1,1,ncomp]);
      endif
      
      
      % generate default priors
      no_priors = false;
      if ~exist('Priors', 'var') || ...
         isempty(Priors)
         % if no Priors is passed in, assume equal, non-normalized Priors
         no_priors = true;
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
      
      % check covariance dimensions
      if size(Sigma,3) ~= ncomp || ...
         size(Sigma,2) ~= nchan || ...
         size(Sigma,1) ~= size(Sigma,2)
         error('gmix_pdf: Sigma is not a valid covariance matrix array');
      endif
      
      
      % check if covariances are pos-def?
      
      
   endif
   
   
   % if at least four arguments
   if nargin > 3
      
      
      % check priors dimensions
      if no_priors
         
         % create priors, but don't normalize
         Priors = ones(1,1,ncomp);
         
      elseif size(Priors,3) ~= ncomp || ...
              size(Priors,1) ~= 1 || ...
              size(Priors,2) ~= 1
         
         error('gmix_pdf: Priors is not a valid array of weights');
      
      else
         
         % rescale any weights passed in to sum to unity
         % Not sure if this is appropriate for gmix_pdf...it might be for higher
         % level functions that call gmix_pdf to normalize weights -EJR 2/2013
         %Priors = Priors ./ sum(Priors,3);
      
      endif
      
            
   endif
   
   
   % if at least five arguments
   if nargin > 4
      
      % check pvalue dimensions
      if isscalar(pvalue)
         
         pvalue = ones(ncomp,1) * pvalue;
      
      elseif ~(isvector(pvalue) && numel(pvalue) == ncomp)
         
         error('gmix_pdf: invalid pvalue specified');
         
      endif
      
      % check pvalue value(s)
      if any(pvalue > 1) || any(pvalue < 0)
         
         error('gmix_pdf: pvalue must fall between zero and 1');
         
      endif
      
   endif
   
   
   % if at least six arguments
   if nargin > 5
      
      % hmm...nothing really to do here
      
   endif
   
   
   
   % initialize outputs compProbs and gmmProb
   compProbs = zeros(nobs,1,ncomp);
   gmmProb = zeros(nobs,1);
      
   % it is not unlikely that zero-filled GMM components get passed to gmix_pdf;
   % if so, we skip these, and assign zero probability to this component; note
   % that this is not actually checking for *valid* components (e.g., pos-def
   % covariance matrices), but simply whether a covariance and/or priors consist
   % entirely of zeros (zero-filled means are perfectly acceptable).
   good_comps = squeeze(Priors > 0) & ...
                any(reshape(Sigma,nchan*nchan,ncomp))';
   ngcomp = sum(good_comps);
   
   % if ngcomp == 0, just return 
   if ngcomp == 0
      return;
   endif
   
   % adjust comonent parameters
   Mu = Mu(1,:,good_comps);
   Sigma = Sigma(:,:,good_comps);
   Priors = Priors(1,1,good_comps);
   
   
   
   
%% NOT SURE WHAT THIS WAS TRYING TO ACCOMPLISH...DELETE ONCE VALIDATED -EJR 2/2013   
%   % it is possible that this routine gets called with zero-sized arguments; 
%   % this should not return an error, but return properly-sized outputs
%   if nobs == 0 || nchan == 0 || ncomp == 0
%      gmmProb = zeros(nobs,0);
%      compProbs = zeros(nobs, ncomp, 0);
%      return;
%   endif

   
   % impute missing values if function was called with the 'impute' option.
   if impute
      I_chan = gmix_impute(I_chan, Mu, Sigma, Priors);
   endif
   
      
   %
   % Begin actual algorithm
   %
      
   % this ugly-looking vectorized code speeds things up by a factor of >100
   % relative to the old looping method, giving probabilities identical to the
   % old to a tolderance of ~1e-11.
   for k=1:ngcomp

      % first, calculate deviation from component means
      dev = I_chan - repmat(Mu(1,:,k), nobs, 1);      

      
%      % this gets used after the loop, but it doesn't hurt to accumulate 
%      % it here, inside the loop
%      root_D_comp(:,1,k) = sqrt(det(Sigma(:,:,k)));
%      
%      % calculate Mahalanobis distances efficiently
%      dev_invC = dev * inv(Sigma(:,:,k));
%      md = sum(dev_invC .* dev, 2);
      
      
      % Use SVD to generate pseudo-determininant|inverse instead of det() and 
      % inv(); we could just use pinv() and eig() functions, but there would
      % be redundant operations
      [U,S,V] = svd(Sigma(:,:,k));
      tol = nchan * max(diag(S)) * eps; % default tolerance used by pinv()
      Svec = diag(S);
      Svec(Svec < tol) = 0; % set singular values to zero
      pseudoDet = prod(Svec(Svec > 0));
      Svec(Svec > 0) = 1 ./ Svec(Svec > 0);
      pseudoInv = V * diag(Svec) * U';
      
      % to be used outside this loop
      root_D_comp(:,1,k) = sqrt(pseudoDet);
      
      % efficiently calculate Mahalanobis distance for each observation
      dev_invC = dev * pseudoInv;
      md = sum(dev_invC .* dev, 2);
      
      
      % apply md threshold based on pvalue
      md(md > chi2inv(pvalue(k), nchan)) = Inf;
      
      
      % copy md into temp array
      md_comp(:,1,k) = md;
            
   endfor % for k=1:ngcomp
   
   
   % finish calculating probabilities outside the loop...not sure this is
   % necessary with the major code vectorization above, since there are now
   % only ncgomp loops anyway -EJR
   P_comp_tmp = repmat( 1 ./ (sqrt((2*pi)^nchan) * root_D_comp), nobs, 1) .* exp(-0.5 * md_comp);
   
   
   % apply prior probabilities, and copy to appropriate columns in compProbs
   compProbs(:,1,good_comps) = P_comp_tmp .* repmat(Priors, [nobs, 1, 1]);
   
   
   % finally, sum up the component probabilities for each observation
   gmmProb = sum(compProbs,3);
   
   
endfunction
