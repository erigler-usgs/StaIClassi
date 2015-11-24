%%
%% Usage:
%%
%% [Mu, Sigma, Priors, BIC] = gmix_fit (I_chan, [ncomp|Mu0[, ...
%%                                               Sigma0[, ...
%%                                               Priors0[, ...
%%                                               opt]]]]);
%%
%%
%% This function generates multichannel, multi-component statistics for a
%% Gaussian Mixture Model (GMM). It employs an expectation maximization (EM)
%% algorithm with a convergence criterion based on change in log-likelihood.
%%
%% INPUTS:
%%
%% I_chan      - multichannel column vector (i.e., rows==observations,
%%               columns==channels). Missing channels are allowed in each
%%               observation, and designated using NaNs as fill values; 
%%               the missing channels are replaced/imputed with their
%%               expected value, conditioned on the available values.
%%
%% ncomp|Mu0   - if scalar, assume GMM with ncomp components; otherwise
%%               second argument is an initial guess for component mean
%%               vectors, where each mean is a 1xnchan row vector, and
%%               components are stored along the 3rd dimension...if not
%%               provided as input, initialize Mu via kmeans, using a
%%               random value consistent with I_chan as initial guess
%%               for each cluster center.
%%               NOTE: since we overload ncomp|Mu0 based on whether the
%%                     second argument is a scalar, there is an ambiguity
%%                     that prevents the user from specifying Mu0 for a 
%%                     univariate, single-component GMM...this is not a 
%%                     problem, since doing so would be incredibly more 
%%                     complicated than just using mean() and var().
%%
%% Sigma0      - third argument is initial guess for component covariance
%%               matrices, where each covariance is a nchanXnchan matrix,
%%               and components are stored along the 3rd dimension...if
%%               not provided as input, initialize using kmeans (if Mu0 
%%               was provided, but Sigma0 was not, Mu0 is used to hold cluster 
%%               centers fixed during kmeans).
%%
%% Priors0     - fourth argument is initial guess for component priors,
%%               where each prior is a 1x1 scalar, and components are
%%               stored along the 3rd dimension...if not provided as
%%               input, initialize with equal priors.      
%% 
%% opt         - the fifth argument is a structure holding various parameters
%%               to control the iterative EM fitting algorithm. Parameters used
%%               by gmix_fit are described here, while opt will also be passed
%%               to other routines called by gmix_fit, like gmix_impute.
%%
%%               opt.llr      - threshold in abs(1-llnew/llold) below which the
%%                              algorithm is considered to converge...'ll*' is
%%                              cumulative log-likelihood of solution;
%%                              default = 1e-6
%%               opt.iter     - number of EM iterations to perform; if opt.llr
%%                              is also set, this specifies the maximum number
%%                              of iterations;
%%                              default = none
%%               opt.???      - 
%%
%%
%% OUTPUTS:
%%
%% Mu          - optimal set of mean vectors...if not requested as ouput,
%%               assume Mu==Mu0 remained fixed.
%%
%% Sigma       - optimal set of covariance matrices...if not requested as
%%               output, assume Sigma==Sigma0 remained fixed.
%%
%% Priors      - optimal set of priors, or weights...if not requested as
%%               ouput, assume Priors==Priors0 remained fixed.
%%
%% BIC         - Bayesian Information Criterion determined from log-likelihood.
%%               BIC describes how well the GMM fit the data, with appropriate 
%%               penalties for over-fitting.
%%
%%
%% NOTE:   One may choose to fix Mu, Sigma, or Priors at their initial values
%%         by *not* requesting them as output. For example:
%%
%%         >> [~, ~, Priors] = gmix_fit(Ichan, Mu0, Sigma0);
%%
%%         ...might be called if the user had prior knowledge of the component
%%         means and covariances (i.e., they assumed known distributions), but 
%%         required component priors (i.e., they desired the mixing ratios).
%%         To ensure that the user either has prior knowledge of the model
%%         components, or will extract such knowledge from this function, the
%%         function does not allow empty output arguments if corresponding
%%         input arguments are not provided. So, the following would be an
%%         invalid call to this function:
%%
%%         >> [Mu, ~, Priors] = gmix_fit(Ichan, ncomp);
%%
%%         ...because the user would never have knowledge of covariance
%%         matrices used inside the function to obtain Mu and Priors as
%%         output.
%%
%% FIXME: This is an inflexible and confusing way to fix parameters. Modify
%%         to accept mask versions of Mu0, Sigma0, and Priors0 that work in
%%         a similar manner to my state space model algorithms (which badly
%%         need to be resurrected and cleaned up anyway). -EJR 10/2012
%%

%%
%% REFERENCES:
%%
%% Calinon, S. (2009), Robot Programming by Demonstration: A Probabilistic
%%   Approach, CRC Press, Taylor & Francis Group, Boca Raton, FL.
%% Gupta, M. R., and Chen, Y. (2010), Theory and use of the EM algorithm, 
%%   Foundations and Trends in Signal Processing, v. 4.3, p. 223-296.
%% Hunt, L., and Jorgensen, M. (2003), Mixture model clustering for mixed data
%%   with missing information, Comp. Stats. and Data Analysis. v. 41, pp. 429-440.
%% Schneider, T. (2001), Analysis of incomplete climate data: estimation of mean
%%   values and covariance matrices and imputation of missing data, Journal of
%%   Climate, v. 14, pp. 853-871.
%%
%%
%% AUTHOR(s):     E. Joshua Rigler (EJR)
%%
%%
%% CHANGELOG:     18 Oct, 2012 - First public release (EJR)
%%                02 Nov, 2012 - very slight modification to call gmix_pdf
%%                               directly instead of gmix_cluster, after
%%                               splitting the clustering and probability
%%                               calculation into two functions.
%%                14 Nov, 2012 - check that Sigma is pos-def covariance
%%                30 Jan. 2013 - minor mods to accommodate and exploit changes 
%%                               in gmix_impute; log-likelihoods now calculated 
%%                               based on imputed data, not just available data.
%%                07 Mar, 2013 - minor mod to address change in gmix_impute()
%%                               that no longer normalizes the component probs.;
%%                               changed back to only considering the available
%%                               data probabilities (as opposed to imputed data
%%                               probabilities) when checking for convergence...
%%                               this is somewhat at odds with Hunt & Jorgensen
%%                               2003, but seems to give convergence properties
%%                               more consistent with what EM theory predicts.
%%                24 May, 2013 - changed back to considering complete-data probs
%%                               after re-reading H&J-2003 yet again, and peeking
%%                               at some other EM implementations; changed default
%%                               llr to 1e-3, then noted that a more sophisticated
%%                               convergence criterion may be desirable, but did
%%                               not implement one.
%%



function [Mu, Sigma, Priors, BIC] = gmix_fit (I_chan, Mu0, Sigma0, Priors0, opt)
	
   %
   % Process input arguments
   %
   
   
   % check for minimal inputs
   if nargin < 1
      error('gmix_fit: At least one input argument is required.');
   endif


   
   % if at least 1 input
   if nargin > 0
      
      % basic check of I_chan dimensions
      if ndims(I_chan) > 2
         error('gmix_fit: I_chan dimensions are not valid');
      endif
      
      nobs = size(I_chan,1);
      nchan = size(I_chan,2);
      
      do_kmeans = true;
      ncomp = 1;
      
      llr = 1e-3;
      iter = Inf;
      
      if ~exist('opt','var')
         opt = struct();
      endif
      
   endif
   
   
   
   % if at least 2 inputs
   if nargin > 1
      
      if isscalar(Mu0) && ...
         (nargin == 2 || (isempty(Sigma0) && isempty(Priors0) ) )
         
         % copy scalar Mu0 into ncomp
         ncomp = Mu0;
         
         % delete Mu0 (will be reconstructed later)
         clear Mu0;
         
      else
         
         ncomp = size(Mu0, 3);
         
      endif
      
      
      % check if Mu0 still exists, and is a valid mean vector array
      if exist('Mu0', 'var') && size(Mu0,1) > 1
         error('gmix_fit: Mu0 is not a valid mean vector array');
      endif
      
      
      % check if Mu0 still exists, and for dimensional consistency
      if exist('Mu0', 'var') && size(I_chan,2) ~= size(Mu0,2)
         error('gmix_fit: I_chan and Mu0 are not dimensionally consistent');
      endif
       
      
   endif
   
   
   
   % if at least 3 inputs
   if nargin > 2
      
      
      % check if Sigma0 is a valid covariance matrix array
      if isempty(Sigma0)
         % do nothing, this simply allows us to pass empty arguments
      elseif size(Mu0,2) ~= size(Sigma0,2) || ...
         size(Mu0,3) ~= size(Sigma0,3) || ...
         size(Sigma0,1) ~= size(Sigma0,2) || ...
         ~isdefinite(Sigma0) % lastly, check for positive definiteness
         
         error('gmix_fit: Sigma0 is not a valid covariance matrix array');
         
      else
      
         % if Sigma0 is not empty and valid, do NOT do kmeans clustering
         do_kmeans = false;
      
      endif
      
   endif
   
   
   
   % if at least 4 inputs
   if nargin > 3
            
      % check if Priors0 is a valid weights array
      if isempty(Priors0)
         % do nothing, this simply allows us to pass empty arguments
      elseif size(Mu0,3) ~= size(Priors0,3) || ...
         size(Priors0,1) > 1 || ...
         size(Priors0,2) > 1
         
         error('gmix_fit: Priors0 is not a valid weights array');
      endif
      
   endif
   
   
   % if 5 inputs
   if nargin > 4
      
      if isstruct(opt)
         
         if isfield(opt, 'llr') && isfield(opt, 'iter')
            llr = opt.llr;
            iter = opt.iter;
         elseif isfield(opt,'llr')
            llr = opt.llr;
            % iter is already set to Inf
         elseif isfield(opt,'iter')
            % if iter is passed with no llr, assume the caller wants a fixed
            % number of iterations, regardless of whether llr threshold is met
            llr = 0;
            iter = opt.iter
         endif
         
      else
         error('gmix_fit: opt must be a structure');
      endif
      
   endif
   
   
   
   %
   % "process" output arguments (i.e., determine what to do based on what was requested)
   %
      
   % check to see if we update Mu
   if ~isargout(1)
      
      if ~exist('Mu0', 'var')
         error('gmix_fit: Mu0 must be specified if it is to remain fixed');
      endif
      
      fix_Mu = true;
      
   else
      
      fix_Mu = false;
   
   endif
   
   
   % check to see if we update Sigma
   if ~isargout(2)
      
      if ~exist('Sigma0', 'var')
         error('gmix_fit: Sigma0 must be specified if it is to remain fixed');
      endif
      
      fix_Sigma = true;
      
   else
      
      fix_Sigma = false;
   
   endif
   
   
   % check to see if we update Priors
   if ~isargout(3)
      
      if ~exist('Priors0', 'var')
         error('gmix_fit: Priors0 must be specified if it is to remain fixed');
      endif
      
      fix_Priors = true;
      
   else
      
      fix_Priors = false;
   
   endif


# It turns out the following assertion is not true, and in fact, it is
# entirely appropriate to replace whole-row-NaNs with the Priors-weighted
# means, at least if/when one assumes that there is uniform *conditional*
# probability of belonging to a GMM component in the absence of any data.
# ...hmmm, having second (or is it third?) thoughts about this, go ahead
#    and remove whole-row-NaNs before proceeding, and adjust nobs, until
#    we understand this better -EJR 11/2012

   % this funtion will fit a GMM, even with incomplete observations, by
   % replacing NaNs with their most likely values conditioned on the data
   % that *are* available (i.e., "imputation"), and adjusting the covariance
   % to account for the associated increase in uncertainty. It can NOT
   % handle observations that are entirely absent (i.e., an entire row of
   % NaNs in I_chan). Issue a warning, remove any all-NaN rows, and 
   % adjust nobs accordingly.
   NaNrows = find(all(isnan(I_chan), 2));
   if numel(NaNrows) > 0
      I_chan(NaNrows,:) = [];
      nobs = nobs - numel(NaNrows);
      warning("gmix_fit: %d completely missing observations removed", numel(NaNrows));
   endif
   
   
   %
   % if they were not passed in, generate Mu0, Sigma0, and Priors0
   % using a very simple, fixed-iteration K-means clustering algorithm
   % FIXME: this should be replaced with a call to a (sub?)function to
   %         improve code readabilty
   %

   if do_kmeans
      
      % initialize Mu0 if it isn't already
      if ~exist('Mu0', 'var')
         
         % randomly sample for Mu0
         minVar = min(I_chan);
         maxVar = max(I_chan);
         rangeVar = maxVar - minVar;
         for k=1:ncomp
            Mu0(1,:,k) = rand(1,nchan) .* rangeVar + minVar;
         endfor
         
         fix_Mu0 = false;
         
      else
         
         fix_Mu0 = true;
      
      endif
      

      % start by standardizing data array and means
      % (see Chapter 7, by Tanioka and Yadohisa, in Challenges at the
      %  Interface of Data Analysis, Computer Science, and Optimization,
      %  Springer-Verlag, 2012)
      I_chan_std = I_chan ./ repmat(max(I_chan) - min(I_chan), nobs, 1);
      
      for k=1:ncomp
         cent(k,:) = Mu0(1,:,k) ./ max(I_chan) - min(I_chan);
      endfor % for k=1:ncomp
      
      
      [jj, ii] = meshgrid(1:ncomp, 1:nobs); jj = jj(:); ii = ii(:);
      
      
      for l=1:25
         
         % calculate distances
         dk = I_chan_std(ii,:) - cent(jj,:);
         
         % substitute zero for any NaNs in dk...this *should* result in a
         % valid cluster assignment that ignores the missing channels;
         % Note: the Euclidian distances will not be comparable from obs.
         %       to obs., but they are from component to component, which
         %       is all we care about for cluster assignment.
         dk(isnan(dk)) = 0;
         
         d = reshape(norm(dk, "rows"), nobs, ncomp);
         
         
         % classify
         [~, idx] = min(d, [], 2);
         
         
         % deal with singular and empty clusters
         ncluster = sum(repmat(idx,[1,ncomp]) == repmat([1:ncomp], [rows(idx), 1]));
         if any(ncluster < 2)
         %if numel(unique(idx)) ~= ncomp
            warning("gmix_fit: empty cluster in kmeans initialization...moving centroid");
            remaining_idx = 1:rows(d);
            for m=setdiff([1:ncomp], find(ncluster >= 2))
               
               cent(m,:) = sum(I_chan_std(idx==m,:))/3;
               
               [~, mmax] = max(d(remaining_idx, m));

               cent(m,:) = cent(m,:) + I_chan_std(mmax,:)/(2+sum(idx==m)); % replace centroid with a "distant" observation
               idx(mmax) = m; % reset idx
               remaining_idx(mmax) = []; % remove "distant" observation from remaining options
               
               % always need at least two obs. in a cluster if (co)variance is
               % to be calculated...repeat previous steps
               [~, mmax] = max(d(remaining_idx, m));
               cent(m,:) = cent(m,:) + I_chan_std(mmax,:)/(1+sum(idx==m)); % replace centroid with a "distant" observation
               idx(mmax) = m; % reset idx               
               remaining_idx(mmax) = []; % remove "distant" observation from remaining options
               
            endfor % endfor m=setdiff([1:ncomp], unique(idx))
                        
         end % endif numel(unique(idx)) ~= ncomp
         
         if fix_Mu0
            
            % if fixed Mu0, nothing will change in subsequent loops
            break;
            
         else
         
            % update centroid
            for k=1:ncomp
               cent(k,:) = mean(I_chan_std(idx==k,:), 1);
            end % endfor k=1:nbStates (update means)
         
         endif
      
      end % endfor l=1:25 (kmeans algorithm for centroids)
      
      
      % compute Mu0 and Sigma0
      for k=1:ncomp
         
         Mu0(1,:,k) = mean(I_chan(idx==k,:),1);
         
         % it is apparently still possible to get a cluster where one of the
         % channels is all NaNs...might be nice to fix this properly one day
         % Note: replacing "cluster mean" with population mean
         if any(isnan(Mu0(1,:,k))(:) )
            Mu0(1,isnan(Mu0(1,:,k)),k) = mean(I_chan(:, isnan(Mu0(1,:,k)) ) );
         endif
            
         % ...and of course, if there was not a single observation for an
         % entire column of I_chan, the mean from Octave's NaN package will
         % return a NaN, which must be fixed before proceeding.
         % Note: replacing "population mean" with 0
         if any(isnan(Mu0(1,:,k))(:) )
            Mu0(1,isnan(Mu0(1,:,k)),k) = 0;
         endif
         
         % turns out covm in NaN package, while very clever about calculating
         % individual covariances from sparse data, cannot guarantee a positive-
         % definite matrix, which is essential to the rest of this algorithm. 
         % A safe, if not very accurate, approach is to just calculate Sigma0 
         % by replacing missing values with their means; 
         
         % the following simply replaces missing values with the corresponding mean
         %I_murep = I_chan(idx==k,:);
         %I_murep(isnan(I_murep)) = repmat(Mu0, [nobs,1])(isnan(I_murep));
         %Sigma0(:,:,k) = covm(I_murep, 'D1'); 
         
         
         % A less robust, but possibly more accurate approach is to find the
         % "nearest" valid covariance matrix using something like Rebonato and 
         % Jackel (1999) or Higham (2002); all such algorithms are designed for
         % correlation matrices, not covariance matrices...simple enough, we 
         % just need to properly scale by univariate variances
         
         % Use corrcoef from Octave's NaN package to return pair-wise correlations 
         % Note: this corrcoef returns NaNs (expected with missing data), and very
         %       occasionally, Infs (not expected); just convert these correlations
         %       to zero, even though this is an arbitrary assumption.
         CorrPair = corrcoef(I_chan(idx==k, :));
         CorrPair(isnan(CorrPair) | isinf(CorrPair)) = 0; 
         
         % find "nearest" covariance by finding nearest correlation matrix to invalid SigNaN
         
%         % experimented with Higham-2002's own ML code, downloaded from:
%         %   http://nickhigham.wordpress.com/2013/02/13/the-nearest-correlation-matrix/
%         % it doesn't even return PSD matrix...this is somehow related to iterative
%         % solver that simply forces correlation matrix diagonal values to equal 1...not acceptable
%         nctol = 100 * eps * norm(Sigma0, 'fro'); % arbitrary tolerance borrowed from isdefinite()
%         CorrNear = nearcorr(CorrZero, nctol, [], [], [], [], 0); % Higham's code
         
         % use "simpler" Rebonato&Jackel-1999 technique (i.e., spectral decomposition);
         [V, D] = eig(CorrPair);
         T = diag(sum(V .* V .* repmat(max(diag(D), 0)', nchan, 1), 2).^(-1)); 
         B = T.^0.5 * V * diag(max(diag(D), 0)).^0.5;
         CorrNear = B * B';
         
         % convert correlation matrix and standard devations to covariance
         sd = std(I_chan(idx==k));
         sd(isnan(sd)) = 0; % convert NaNs to zeros (somewhat arbitrary assumption)
         SigNear = diag(sd) * CorrNear * diag(sd);
         
         % copy SigNear into initial guess for Sigma
         Sigma0(:,:,k) = SigNear;
         
         if isdefinite(Sigma0(:,:,k)) < 0; % this should never happen
            keyboard
         endif
         
         Priors0(1,1,k) = sum(idx==k);
         
      end % endfor k=1:ncomp (compute covariances)

      
   endif % if do_kmeans
   
   
   
   % define Priors0 if it was not passed in, or created via kmeans
   if ~exist('Priors0', 'var')
      
      % assume equal weights if none are passed with Mu0 and Sigma0
      Priors0 = ones(1,1,ncomp) ./ sum(ones(1,1,ncomp));
            
   endif
   
   
   % priors should sum to unity for EM algorithm...we rescale them back to give
   % the original sum before exiting
   PriorSum = sum(Priors0, 3);
   Priors0 = Priors0 ./ PriorSum;
   
   
   %
   % finally, EM algorithm to fit GMM parameters
   %
   
   % initialize some necessary values
   
   loglik_old = -realmax;
   nbstep = 0; % this just for debugging
   
   Mu = Mu0;
   Sigma = Sigma0;
   Priors = Priors0;
   
   
   % mask of NaNs for repeated imputation below
   I_chan_miss = zeros(size(I_chan));
   I_chan_miss(isnan(I_chan)) = NaN; % this is a missing channel mask
   
      
   % a somewhat matrix-optimized EM implementation
   
   %% cannot use "inf" in for-loop, must use a while-loop, and increment our own
   %% own counter
   %for i=1:iter
   i=0;
   
   while i < iter
      
      % increment counter
      i++;
      
      % E-step %%%%%%%%%%
      
      % (re)impute missing values
      [I_chan, Pix, I_chan_comps, S_resid] = gmix_impute(I_chan + I_chan_miss, ...
                                                         Mu, Sigma, Priors, opt);
      
      % normalize Pix
      PixNorm = Pix ./ (repmat(sum(Pix, 3), [1, 1, ncomp]) + realmin); % never divide by zero
      
      % accumulate posterior probs
      E = sum(PixNorm);
      
      
      % M-step %%%%%%%%%%
      
      % loop over components
      for k=1:ncomp
         
         % update priors
         if ~fix_Priors
            
            % this is an average weight of components
            Priors(k) = E(k) / nobs;
         
         endif
         
         
         % update means
         if ~fix_Mu
            
            % this is a weighted average of component expected values
            Mu(:,:,k) = [I_chan_comps(:,:,k)' * PixNorm(:,1,k)]' / E(k);
            
         endif
         
         
         % update covariance
         if ~fix_Sigma
                        
            % this is a weighted average of component covariances
            dev = I_chan_comps(:,:,k) - repmat(Mu(1,:,k), [nobs, 1]);
            Sigma(:,:,k) = ([repmat(PixNorm(:,1,k)', [nchan,1]) .*  dev'] * dev) / E(k);
            
            % the previous vectorized approach to weighting the deviations can
            % produce asymmetric covariances due to round-off error...symmetrize
            Sigma(:,:,k) = (Sigma(:,:,k) + Sigma(:,:,k)') / 2;
            
            % adjust for imputation error
            % NOTE: Hunt&Jorgensen did NOT calculate weighted data covariances 
            %       and weighted covariance adjustments separately; rather they
            %       calculated a weighted average of the adjusted covariances
            %       for each observation...I confirmed that this doesn't matter, 
            %       and our approach allows for some vectorization efficiency
            % FIXME: should probably ensure all the S_resid's are also symmetric, 
            %        but this must be done inside gmix_impute()
            Sigma(:,:,k) = Sigma(:,:,k) + permute(sum(repmat(PixNorm(:,1,k), ...
                                                             [1,nchan,nchan]) .* ...
                                                      S_resid(:,:,:,k), 1), ...
                                                  [2,3,1]) / E(k);
            
                        
         endif
         
      endfor % for k=1:ncomp
      
      
      % stopping criterion %%%%%%%%%%
            
      % BIG NOTE TO SELF...
      % After multiple back-and-forths, I am now convinced that the log-liklihood
      % that is supposed to get maximized by EM is the *complete-data* log-liklihood, 
      % not just the observed data log-liklihood. This is what Hunt&Jorgensen write,
      % and it is what is done in the PMTK3 toolbox (which does all StaIClassi
      % does and more...wish I'd discovered it sooner).
      %
      % The problem was that I kept thinking it was NOT increasing because the
      % log-liklihood ratio from iteration i to i+1 didn't necessarily always
      % decrease. THE TWO ARE NOT EQUIVALENT! I don't fully understand why this
      % is the case, but once I looked closely at loglik, and not just llr,
      % it became clear that loglik would always increase, even when llr might
      % might not decrease from iteration i to i+1. All this means is that the
      % error surface is not smooth, which should not be a surprise.
      %
      % So, the problem now is to determine the ideal convergence criterion.
      % You wouldn't think this would be difficult, but a survey of various
      % EM (and other iterative solvers), suggests that there is no concensus
      % on the best way to do this.
      
      % try a slighly modified version of the original complete-data loglik
      F = gmix_pdf(I_chan, Mu, Sigma, Priors);
      loglik = sum(log(F(F>0)));
      
      
      % stop process depending on the rate of increase of the log likelihood
      if abs((loglik/loglik_old)-1) < llr
         break;
      end % endif (abs((loglik/loglik_old)-1) < llr
      
      
      printf('iter: %d; Loglik=%e; abs(llr-1)=%1.12e; ncomp=%d\n', i, loglik, abs((loglik/loglik_old)-1), ncomp);

      
      loglik_old = loglik;
      
      
      fflush(stdout);
      
   %endfor; % for i=1:iter
   endwhile; % while i<iter
   
   
   % calculate BIC (i.e., -2*log(L) + k*log(n))
   BIC = -2*loglik + (~fix_Mu * numel(Mu) + ...
                      ~fix_Sigma * numel(Sigma) + ...
                      ~fix_Priors * numel(Priors)) * log(nobs);
   
   
   % scale Priors so that its sum equals PriorSum determined above
   Priors = Priors * PriorSum;
   
endfunction





