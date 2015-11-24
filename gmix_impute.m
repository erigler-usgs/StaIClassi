%%
%% Usage:
%%
%% [I_chan_out, ProbComp, I_chanComp, SigmaResidComp] = gmix_impute (I_chan_in, ...
%%                                                                   Mu[, ...
%%                                                                   Sigma[, ...
%%                                                                   Priors[, ...
%%                                                                   opt]]]])
%%
%%
%% This function uses a mutivariate Gaussian mixture model (GMM) to 'impute'
%% missing values in the multivariate data set I_chan using linear regression
%% 
%% 
%% INPUTS:
%%
%% I_chan_in    - multichannel column vector (i.e., rows==observations,
%%                columns==channels). Missing values specified with NaN
%%                values.
%%
%% Mu           - mean vectors for GMM...each mean is a 1Xnchan row vector,
%%                and GMM components are stored along the 3rd dimensioon.
%%
%% Sigma        - (optional) covariance matrices for GMM; each covariance is
%%                a nchanXnchan matrix, and GMM components are stored along
%%                the 3rd dimension...if not passed, assume uniform circular
%%                GMMs (i.e., common variance, diagonal matrices) with unity
%%                variance.
%%                NOTE: the default value is consistent with gmix_cluster, 
%%                where doing this is essentially equivalent to a step in
%%                kmeans clustering algorithm; it's not clear if there is
%%                an equally compelling reason to do it here. -EJR 10/2012
%%
%% Priors       - (optional) prior probabilities of GMM component membership;
%%                each prior, or weight, is a 1X1 scalar, and GMM components
%%                are stored along the 3rd dimension...if not passed, assume
%%                equal priors.
%%
%% opt          - (optional) structure holding options for imputation; some of
%%                these fieldnames and values are for routines adapted from
%%                T. Schneider's Imputation toolbox, which are included for
%%                testing and intercomparisons with our default fixed ridge
%%                parameter regularization):
%%                * opt.regpar      - fixed ridge parameter...if 'mridge' is set
%%                                    via opt.regress, this will execute more
%%                                    slowly, since the Schneider algorithm will 
%%                                    be called instead of our faster, matrix-
%%                                    optimized algorithm...setting both is only
%%                                    useful for debugging/testing
%%                * opt.regress     - == 'mridge', optimize and apply scalar
%%                                    ridge parameter to each regression
%%                                    == 'iridge', optimize and apply vector
%%                                    ridge parameters to each regression
%%                                    == 'ttls', not implemented yet
%%                * opt.relvar_res  - used to construct lower bound on optimized
%%                                    ridge parameter(s)
%%                * opt.minvarfrac  - used to construct upper bound on
%%                                    optimized ridge parameter(s)
%%                * opt.regavail    - if true, available obs. will be (re)imputed
%%                                    using regularized regression coefficients,
%%                                    which leads to a smoother solution, but 
%%                                    which violates certain assumptions that 
%%                                    guarantee convergence of the EM algorithm 
%%                                    (and probably other nonlinear solvers).
%%                                    If false, do not regress available data
%%                                    along with missing data (this is default); 
%%                                    NOTE: Schneider's optimization routines
%%                                          for the ridge parameter(s) cannot
%%                                          handle regularized regression of
%%                                          available obs. without significant 
%%                                          modifications, if at all; so, this
%%                                          option is only useful if fixed ridge
%%                                          parameters are specified via regpar.
%%                                          
%%
%%                FIXME: it might be nice to allow sampling the distributions,
%%                        a la Di Zio et al.; this may require some creative
%%                        thinking about how to handle imp_option, and will 
%%                        almost certainly prove to be slower, so we leave it
%%                        as a FIXME for now.
%%
%%
%% OUTPUTS:
%%
%% I_chan_out   - multichannel column vector with missing values imputed using
%%                conditional expectations
%%
%% ProbComp     - probabilities of available data for each component and obs.;
%%                if these probabilities are normalized so that their sum across
%%                components equals 1, they can be multiplied by I_chanComp
%%                to obtain I_chan_out.
%%
%% I_chanComp   - unweighted component means
%%
%% SigmaResidComp    - unweighted component covariances
%%

%%
%% REFERENCES:
%%
%% Calinon, S. (2009), Robot Programming by Demonstration: A Probabilistic
%%   Approach, CRC Press, Taylor & Francis Group, Boca Raton, FL.
%% Di Zio, M., Guarnera, U., and Luzi, Orietta (2007), Imputation through finite 
%%   Gaussian mixture models, Computational Statistics and Data Analysis, v.51,
%%   p. 5305-5316.
%% Hunt, L., and Jorgensen, M. (2003), Mixture model clustering for mixed data
%%   with missing information, Comp. Stats. and Data Analysis. v. 41, pp. 429-440.
%% Reece, S., and Roberts, S. (2010), Generalized Covariance Union: a unified
%%   approach to hypothesis merging in tracking, IEEE Trans. Aerospace Elec.
%%   Systems, v. 46, n. 1, p. 207-221.
%% Schneider, T. (2001), Analysis of incomplete climate data: estimation of mean
%%   values and covariance matrices and imputation of missing data, Journal of
%%   Climate, v. 14, pp. 853-871.
%% Sung, H. G. (2004), Gaussian mixture regression and classification, Ph.D.
%%   Thesis, Rice Univ., Houston, TX.
%%
%%
%% AUTHOR(s):     E. Joshua Rigler (EJR)
%%
%%
%% CHANGELOG:     25 Oct, 2012 - First public release (EJR)
%%                02 Nov, 2012 - very slight modification to call gmix_pdf
%%                               directly instead of gmix_cluster, after
%%                               splitting the clustering and probability
%%                               calculation into two functions.
%%                13 Nov, 2012 - added options allowing the user to regularize
%%                               regression with either a simple scalar ridge
%%                               parameter, or by calling one of the routines
%%                               adapted from Tapio Schneider's Imputation ML
%%                               toolbox.
%%                28 Jan, 2013 - modified outputs so that Sigma_resid returns
%%                               a single covariance for each observation, and
%%                               ProbComp is normalized
%%                30 Jan, 2013 - modified outputs to remove Sigma_resid, since
%%                               this is something that should be caluclated
%%                               elsewhere; added outputs for component means
%%                               and covariances for each observation.
%%                28 Feb, 2013 - replaced calls to inv() with pinv()
%%                07 Mar, 2013 - ProbComp is no longer normalized...not sure why
%%                               this ever seemed like a good idea.
%%

function [I_chan_out, ProbComp, I_chanComp, SigmaResidComp] = gmix_impute (I_chan_in, ...
                                                                       Mu, ...
                                                                       Sigma, ...
                                                                       Priors, ...
                                                                       opt);
   
   %
   % process inputs
   %
   
   % check minimial arguments
   if nargin < 2
      error('gmix_impute: at least two input arguments are required');
   endif
   
   
   % if at least 2 input arguments
   if nargin > 1
      
      % basic check of I_chan dimensions
      if ndims(I_chan_in) > 2
         error('gmix_impute: I_chan dimensions are not valid');
      endif
      
      nobs = size(I_chan_in,1);
      nchan = size(I_chan_in,2);
      
      
      % basic check of Mu dimension
      if size(Mu,1) ~= 1 || ...
         size(Mu,2) ~= nchan
         error('gmix_impute: I_chan and Mu are not dimensionally consistent');
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
      
      
      % generate default imputation options
      ridge_param = 0;
      regavail = false;
      if ~exist('opt', 'var')
         opt = struct();   % create empty structure
      endif
      
   endif
   
   
   % if at least three input arguments
   if nargin > 2
      
      % check covariance dimensions
      if size(Sigma,3) ~= ncomp || ...
         size(Sigma,2) ~= nchan || ...
         size(Sigma,1) ~= size(Sigma,2)
         error('gmix_impute: Sigma is not a valid covariance matrix array');
      endif
      
      
      % check if covariances are pos-def?
      
      
   endif
   

   % if at least 4 input arguments
   if nargin > 3
      
      % check priors dimensions
      if numel(Priors) == 0
         
         % turns out we may want to pass an empty Priors if imp_option
         % is passed, and we do NOT want the priors to all equal unity
         % (as opposed to sum to unity).
         
         % create priors, but don't normalize
         Priors = ones(1,1,ncomp);
         
      elseif size(Priors,3) ~= ncomp || ...
         size(Priors,1) ~= 1 || ...
         size(Priors,2) ~= 1
         error('gmix_impute: Priors is not a valid array of weights');
      else
         % make certain weigts sum to unity,
         Priors = Priors ./ sum(Priors,3);
      endif
      
   endif
   
   
   % if at least five arguments
   if nargin > 4
      
      if isstruct(opt)
         
         if isfield(opt, 'regavail') && opt.regavail
            regavail = true;
         endif
         
         if isfield(opt, 'regpar')
            ridge_param = opt.regpar;
         endif
         
         if isfield(opt, 'regress')
            % T. Schneider's routines cannot self-regress...it's not clear if
            % this is a mathematical limitation, or simply an algorithmic one
            regavail = false;
         endif
         
      else
         
         error('gmix_impute: opt must be a structure');
     
      endif
      
      
      
   endif

   
   %
   % begin actual algorithm
   %
   
   %
   % this bit of coding jujitsu speeds up imputation considerably when there
   % are groupings of missing channel configurations because the matrix
   % inversion required to get the regression coefficients only needs to be
   % performed once per grouping
   %
   
   chan_miss = isnan(I_chan_in);
   chan_avail = ~chan_miss;
   chan_unique = unique(chan_avail, 'rows');
   
   % initialize output array, starting with imputed values of zero
   % NOTE: Schneider (2001) leaves available data points alone (i.e., he only
   %       imputes *missing* data). It may be desirable to impute available
   %       data, which is another way to say we allow for uncertainty in the 
   %       available observations. 
   
   if isargout(1)
      I_chan_out = I_chan_in;
      I_chan_out(chan_miss) = 0; 
      if regavail
         I_chan_out(chan_avail) = 0;
      endif
   endif
   
   if isargout(2)
      % initialize ProbComp
      ProbComp = zeros(nobs, 1, ncomp);
   endif
   
   if isargout(3)
      I_chanComp = I_chan_in;
      I_chanComp(chan_miss) = 0;
      if regavail
         I_chanComp(chan_avail) = 0;
      endif
      I_chanComp = repmat(I_chanComp, [1,1,ncomp]);
   endif
   
   if isargout(4)
      SigmaResidComp = zeros(nobs, nchan, nchan, ncomp);
   endif
   
   
%%%%   % initialize imputation residual covariance matrices for each component
%%%%   Sigma_resid = zeros(nobs, size(Sigma,1), size(Sigma,2));
   
   
   % loop over unique channel configurations, thereby minimizing number of 
   % matrix inversions.
   
   for u=1:size(chan_unique,1)            
      
      %
      % (re)determine indices necessary to process current set of 
      % missing/available channels
      %
      rimpute = find(ismember(chan_avail, chan_unique(u,:), 'rows'));
      cavail = find(chan_unique(u,:));
      cmiss = find(~chan_unique(u,:));
      
      
      %
      % (re)initialize temporary arrays for current set of missing/available channels
      %
      condexp_miss = zeros(numel(rimpute), numel(cmiss), ncomp);
      condexp_avail = zeros(numel(rimpute), numel(cavail), ncomp);
      condcov_miss = zeros(numel(cmiss), numel(cmiss), ncomp);
      condcov_avail = zeros(numel(cavail), numel(cavail), ncomp);
       
       
      % calculate component posterior probabilities for available obs.
      if numel(cavail) == 0
         
         % a simple and practical approach is to just assign NaNs...
         I_chan_out(rimpute,:) = NaN;
         
         % ...one might also assign the "expected value" for this distribution,
         %    which is the weighted sum of the means of the GMM, but it's not 
         %    clear this is a reasonable default
         
         % skip to the end of the u loop, since nothing else needs to be done
         continue;
      
      else
         
         % calculate probabilities for *available* data
         [~, Pix] = gmix_pdf(I_chan_in(rimpute,cavail), Mu(1, cavail, :), ...
                             Sigma(cavail, cavail, :), Priors(1, 1, :));
         
      endif
      
      
      % normalize Pix components to sum to unity
      PixNorm = Pix ./ (repmat(sum(Pix, 3), [1, 1, ncomp]) + realmin); % never divide by zero
      
      
      if any(isnan(PixNorm(:)))
         keyboard
      endif
      
      % cannot seem to avoid looping over components
      for k=1:ncomp
         
         % generate covariance submatrices
         S_avail = Sigma(cavail, cavail, k);
         S_avail_miss = Sigma(cavail, cmiss, k);
         S_miss = Sigma(cmiss, cmiss, k);
         
         
         % apply ridge parameter and calculate regression coeffients
         if ~isfield(opt,'regress')
            
            
            % use avail covariance and avail-missing cross-covariance matrices
            % to generate regression coefficients
%            Bmiss = inv(S_avail + ridge_param^2 * diag(diag(S_avail))) * ...
%                    S_avail_miss;
            % pinv is slower, but more accurate, especially for poorly 
            % conditioned covariance matrices
            Bmiss = pinv(S_avail + ridge_param^2 * diag(diag(S_avail))) * ...
                    S_avail_miss;
            
            
            if regavail
               % will equal identity matrix if ridge parameter is zero
%               Bavail = inv(S_avail + ridge_param^2 * diag(diag(S_avail))) * ...
%                        S_avail; 
               % pinv is slower, but more accurate, especially for poorly 
               % conditioned covariance matrices
               Bavail = inv(S_avail + ridge_param^2 * diag(diag(S_avail))) * ...
                        S_avail; 
            endif
                 
            S_imp_miss = Sigma(cmiss, cmiss, k) - S_avail_miss' * Bmiss;
            
            if regavail
               % will equal zero if Bavail is identity matrix
               S_imp_avail = Sigma(cavail, cavail, k) - S_avail' * Bavail;
            endif
            
            
         else
            
            switch opt.regress
               case 'mridge'
                  
                  
                  [Bmiss, S_imp_miss] = mridge(S_avail, ...
                                               S_miss, ...
                                               S_avail_miss, ...
                                               nobs, opt);
                  
                  if regavail
                     [Bavail, S_imp_avail] = mridge(S_avail, ...
                                                    S_miss, ...
                                                    S_avail, ...
                                                    nobs, opt);
                  endif
                  
                  
               case 'iridge'
                  
                  [Bmiss, S_imp_miss] = iridge(S_avail, ...
                                               S_miss, ...
                                               S_avail_miss, ...
                                               nobs, opt);
                  
                  if regavail
                     [Bavail, S_imp_avail] = iridge(S_avail, ...
                                                    S_miss, ...
                                                    S_avail, ...
                                                    nobs, opt);
                  endif
                  
                  
               otherwise
               
               error('gmix_impute: unrecognized regression method in imp_option');
               
            endswitch
            
         endif
         
         
         %
         % calculate conditional expectations
         %
         mean_avail = repmat(Mu(1, cavail, k), [numel(rimpute), 1]);
         mean_miss = repmat(Mu(1, cmiss, k), [numel(rimpute), 1]);
         dev_avail = I_chan_in(rimpute, cavail) - mean_avail;
         
         
         % zero-mean regression
         dev_miss = dev_avail * Bmiss;
         
         if regavail
            dev_avail = dev_avail * Bavail;
         endif
         
         
         % unweighted conditional expectations per component
         condexp_miss(:,:,k) = dev_miss + mean_miss;
         
         if regavail
            condexp_avail(:,:,k) = dev_avail + mean_avail;
         endif
         
         
         % store unweighted conditional covariances for later
         condcov_miss(:,:,k) = S_imp_miss;
         
         if regavail
            condcov_avail(:,:,k) = S_imp_avail;
         endif
         
         
      endfor % for k=1:ncomp (#1)
      
      
      
      %
      % many of our outputs have a tendency to grow very large, so we only copy 
      % temporary values to them if they are actually requested
      %
      
      if isargout(1)
         % sum up weigted conditional means; copy to output array
         I_chan_out(rimpute, cmiss) = sum(condexp_miss .* ...
                                          repmat(PixNorm, [1, numel(cmiss), 1]), 3);
         if regavail
            I_chan_out(rimpute, cavail) = sum(condexp_avail .* ...
                                              repmat(PixNorm, [1, numel(cavail), 1]), 3);
         endif
      endif
      
      
      if isargout(2)
         % copy weights to output array
         ProbComp(rimpute, 1, :) = Pix;
      endif
      
      
      if isargout(3)
      
         % copy unweighted component conditional means to output array
         I_chanComp(rimpute, cmiss, :) = condexp_miss;
         if regavail
            I_chanComp(rimpute, cavail, :) = condexp_avail;
         endif
      endif
      
      
      if isargout(4)
      
         % copy unweighted component residual covariances
         SigmaResidComp(rimpute, cmiss, cmiss, :) = repmat(permute(condcov_miss, [4,1,2,3]), ...
                                                           [numel(rimpute), 1, 1, 1]);
         if regavail
            SigmaResidComp(rimpute, cavail, cavail, :) = repmat(permute(condcov_avail, [4,1,2,3]), ...
                                                                [numel(rimpute), 1, 1, 1]);
         endif
      endif
      
      
      
%%%%         %
%%%%         % calculate the imputed residual covariance matrix for each component
%%%%         % (this is useful, but is beyond the scope of this function; also,
%%%%         %  this is how Calinon-2009 did things, which may not be the most
%%%%         %  appropriate...commenting out for now, should delete in the future)
%%%%         %
%%%%         
%%%%         
%%%%         Sigma_resid(cmiss, cmiss, k) =  Sigma_resid(cmiss, cmiss, k) + ...
%%%%                                         sum(repmat(reshape(PixNorm(:, 1, k) .* ...
%%%%                                                            PixNorm(:, 1, k), ...
%%%%                                                            [1,1,numel(rimpute)]), ...
%%%%                                                    [numel(cmiss), numel(cmiss)]) .* ...
%%%%                                             repmat(S_imp_miss, [1, 1, numel(rimpute)]), 3);
%%%%         
%%%%         if regavail
%%%%            % this only adds zero unless the ridge parameter(s) are non-zero
%%%%            Sigma_resid(cavail, cavail, k) =  Sigma_resid(cavail, cavail, k) + ...
%%%%                                              sum(repmat(reshape(PixNorm(:, 1, k) .* ...
%%%%                                                                 PixNorm(:, 1, k), ...
%%%%                                                                 [1,1,numel(rimpute)]), ...
%%%%                                                         [numel(cavail), numel(cavail)]) .* ...
%%%%                                                  repmat(S_imp_avail, [1, 1, numel(rimpute)]), 3);
%%%%         endif
%%%%

      
   endfor % endfor u=1:size(chan_unique,1)
   
   
%%%%   %% NOTE: both Hunt and Jorgensen (2003), and Di Zio et al (2007) seem to 
%%%%   %%       suggest that one should normalize by 1/(nobs*Priors(k,t+1)), when
%%%%   %%       using this in an EM-based GMM estimator, but don't say why. So, we
%%%%   %%       will normalize Sigma_resid using the weights passed into this so as
%%%%   %%       to keep it self-consistent. If this breaks an EM estimator, the EM
%%%%   %%       estimator should be modified to rescale Sigma_resid accoring to its
%%%%   %%       updated Priors.
%%%%   for k=1:ncomp
%%%%      Sigma_resid(:,:,k) = Sigma_resid(:,:,k) * 1/(nobs*Priors(1,1,k));
%%%%   endfor
   
   
   
   
endfunction







%%
%% Subfunctions adapted from Tapio Schneider's Imputation toolbox. These should
%% probably be separate files in this directory, but until/unless I choose to
%% rewrite these as per StaIClassi standards, I'll just place them here to
%% ensure that any modifications that get made are updated along with the GMM
%% imputation code. -EJR 11/2012
%%


function [V, d, r] = peigs(A, rmax)
%PEIGS   Finds positive eigenvalues and corresponding eigenvectors. 
%  
%    [V, d, r] = PEIGS(A, rmax) determines the positive eigenvalues d
%    and corresponding eigenvectors V of a matrix A. The input
%    parameter rmax is an upper bound on the number of positive
%    eigenvalues of A, and the output r is an estimate of the actual
%    number of positive eigenvalues of A. The eigenvalues are returned
%    as the vector d(1:r). The eigenvectors are returned as the
%    columns of the matrix V(:, 1:r).
%
%    PEIGS calls the function EIGS to compute the first rmax
%    eigenpairs of A by Arnoldi iterations. Eigenpairs corresponding
%    eigenvalues that are nearly zero or less than zero are
%    subsequently discarded.
%
%    PEIGS is only efficient if A is strongly rank-deficient with
%    only a small number of positive eigenvalues, so that rmax <<
%    min(size(A)). If A has full rank, EIG should be used in place of
%    PEIGS.
%
%    See also: EIGS.

  error(nargchk(2, 2, nargin))                    % check number of input arguments 
  [m, n]  = size(A);
 
  if rmax > min(m,n)
    rmax  = min(m,n);                             % rank cannot exceed size of A 
  end
      
  % get first rmax eigenvectors of A
  %warning off
  %eigsopt.disp = 0;                               % do not display eigenvalues
  % eigs in Octave depends on ARPACK, which is not reliable
  if exist("OCTAVE_VERSION","builtin")
      [V, d] = eig(full(A));
      d=diag(d);
      [d, di] = sort(d);
      V = V(:,di);
      d=d(end-rmax+1:end);
      V = V(:,end-rmax+1:end);
  else
      [V, d]       = eigs(A, rmax, 'lm'); %, eigsopt);
  end
  %warning on
  
  % output of eigs differs in different Matlab versions
  if prod(size(d)) > rmax
    d          = diag(d);                         % ensure d is vector
  end

  % ensure that eigenvalues are monotonically decreasing
  [d, I]       = sort(d, 'descend');
  V            = V(:, I);
  
  % estimate number of positive eigenvalues of A
  d_min        = max(d) * max(m,n) * eps; 
  r            = sum(d > d_min);
				 
  % discard eigenpairs with eigenvalues that are close to or less than zero
  d            = d(1:r);
  V            = V(:, 1:r);
  d            = d(:);				  % ensure d is column vector

endfunction



function g = gcvfctn(h, d, fc2, trS0, dof0)
%GCVFCTN    Evaluate object function for generalized cross-validation.
%
%   GCVFCTN(h, d, fc2, trS0, dof0) returns the function values of the
%   generalized cross-validation object function
%
%                     trace [ S0 + F' * diag(g.^2) * F ]
%              G(h) = ---------------------------------- 
%                           ( dof0 + sum(g) )^2
%
%   where g = h^2 ./ (d + h^2) = 1 - d.^2 ./ (d + h^2). The argument h
%   of the GCV function is the regularization parameter, and d is a
%   column vector of eigenvalues (see GCVRIDGE for the meaning of the
%   other symbols above). GCVFCTN is an auxiliary routine that is
%   called by GCVRIDGE. The input arguments are defined in GCVRIDGE:
%
%        h:  regularization parameter,
%        d:  column vector of eigenvalues of cov(X),
%      fc2:  row sum of squared Fourier coefficients, fc2=sum(F.^2, 2),
%     trS0:  trace(S0) = Frobenius norm of generic part of residual matrix,
%     dof0:  degrees of freedom in estimate of residual covariance
%            matrix when regularization parameter is set to zero
  
%   Adapted from GCVFUN in Per Christian Hansen's REGUTOOLS Toolbox.

  filfac = (h^2) ./ (d + h^2);            
  g      = ( sum(filfac.^2 .* fc2) + trS0 ) / (dof0 + sum(filfac))^2;
endfunction



function h_opt = gcvridge(F, d, trS0, n, r, trSmin, options)
%GCVRIDGE   Finds minimum of GCV function for ridge regression.
%
%   GCVRIDGE(F, d, trS0, n, r, trSmin, OPTIONS) finds the
%   regularization parameter h that minimizes the generalized
%   cross-validation function
%
%                         trace S_h
%                 G(h) = ----------- 
%                          T(h)^2
%
%   of the linear regression model Y = X*B + E. The data matrices X
%   and Y are assumed to have n rows, and the matrix Y of dependent
%   variables can have multiple columns. The matrix S_h is the second
%   moment matrix S_h = E_h'*E_h/n of the residuals E_h = Y - X*B_h,
%   where B_h is, for a given regularization parameter h, the
%   regularized estimate of the regression coefficients,
% 
%                B_h = inv(X'*X + n h^2*I) * X'*Y.
%
%   The residual second second moment matrix S_h can be represented
%   as
%
%                S_h = S0 + F' * diag(g.^2) * F
%
%   where g = h^2 ./ (d + h^2) = 1 - d.^2 ./ (d + h^2) and d is a
%   column vector of eigenvalues of X'*X/n. The matrix F is the matrix
%   of Fourier coefficients. In terms of a singular value
%   decomposition of the rescaled data matrix n^(-1/2) * X = U *
%   diag(sqrt(d)) * V', the matrix of Fourier coefficients F can be
%   expressed as F = n^(-1/2) * U' * Y. In terms of the eigenvectors V
%   and eigenvalues d of X'*X/n, the Fourier coefficients are F =
%   diag(1./sqrt(d)) * V' * X' * Y/n. The matrix S0 is that part of
%   the residual second moment matrix that does not depend on the
%   regularization parameter: S0 = Y'*Y/n - F'*F.
% 
%   As input arguments, GCVRIDGE requires:
%        F:  the matrix of Fourier coefficients,
%        d:  column vector of eigenvalues of X'*X/n,
%     trS0:  trace(S0) = trace of generic part of residual 2nd moment matrix,
%        n:  number of degrees of freedom for estimation of 2nd moments,
%        r:  number of nonzero eigenvalues of X'*X/n,
%   trSmin:  minimum of trace(S_h) to construct approximate lower bound
%            on regularization parameter h (to prevent GCV from choosing
%            too small a regularization parameter).
%
%   The vector d of nonzero eigenvalues of X'*X/n is assumed to be
%   ordered such that the first r elements of d are nonzero and ordered 
%   from largest to smallest.
%  
%   The input structure OPTIONS contains optional parameters for the
%   algorithm:
%
%     Field name           Parameter                                  Default
%
%     OPTIONS.minvarfrac Minimum fraction of total variation in X     0
%                        that must be retained in the 
%                        regularization. From the parameter 
%                        OPTIONS.minvarfrac, an approximate upper 
%                        bound for the regularization parameter is
%                        constructed. The default value 
%                        OPTIONS.minvarfrac = 0  corresponds to no
%                        upper bound for the regularization parameter.   
  
%   References:
%   GCVRIDGE is adapted from GCV in Per Christian Hansen's REGUTOOLS
%       toolbox:
%   P.C. Hansen, "Regularization Tools: A Matlab package for
%       analysis and solution of discrete ill-posed problems,"
%       Numer. Algorithms, 6 (1994), 1--35.
%
%   see also: 
%   G. Wahba, "Spline Models for Observational Data",
%       CBMS_NSF Regional Conference Series in Applied Mathematics,
%       SIAM, 1990, chapter 4.

  error(nargchk(6, 7, nargin))     % check number of input arguments 

  % sizes of input matrices
  d      = d(:);                   % make sure that d is column vector
  if length(d) < r
    error('All nonzero eigenvalues must be given.')
  end

  % ================           process options        ======================
  if nargin == 6 | isempty(options)
    fopts       = [];
  else
    fopts       = fieldnames(options);
  end
    
  minvarfrac    = 0;
  if strmatch('minvarfrac', fopts)
    minvarfrac = options.minvarfrac;
    if minvarfrac < 0 | minvarfrac > 1
      error('OPTIONS.minvarfrac must be in [0,1].')
    end
  end
  % ========================================================================
  
  p      = size(F, 1);
  if p < r
    error(['F must have at least as many rows as there are nonzero' ...
	   ' eigenvalues d.']) 
  end
  % row sum of squared Fourier coefficients
  fc2    = sum(F.^2, 2);
  
  % accuracy of regularization parameter 
  h_tol  = .2/sqrt(n);        
  
  % heuristic upper bound on regularization parameter
  varfrac = cumsum(d)/sum(d);
  if minvarfrac > min(varfrac)
    d_max           = interp1(varfrac, d, minvarfrac, 'linear');
    h_max           = sqrt( d_max );
  else            
    h_max           = sqrt( max(d) ) / h_tol;    
  end
    
  % heuristic lower bound on regularization parameter
  if trS0 > trSmin
    % squared residual norm is greater than a priori bound for all 
    % regularization parameters
    h_min         = sqrt(eps);
  else
    % find squared residual norms of truncated SVD solutions
    rtsvd         = zeros(r, 1);
    rtsvd(r)      = trS0;
    for j = r-1: -1: 1
      rtsvd(j)    = rtsvd(j+1) + fc2(j+1);
    end
    % take regularization parameter equal to square root of eigenvalue 
    % that corresponds to TSVD truncation level with residual norm closest 
    % to a priori bound trSmin
    [dummy, rmin] = min(abs(rtsvd - trSmin));
    h_min         = sqrt( max( d(rmin), min(d)/n ) );
  end

  if h_min < h_max
    % find minimizer of GCV function
    minopt = optimset('TolX', h_tol, 'Display', 'off');
    % fminbnd in Octave doesn't know what to do with the passed parameters 
    % (ML docs don't describe this behaviro either)
    %h_opt  = fminbnd('gcvfctn', h_min, h_max, minopt, d(1:r), fc2(1:r), trS0, n-r);
    h_opt = fminbnd(@(h) gcvfctn(h, d(1:r), fc2(1:r), trS0, n-r), h_min, h_max, minopt);
  else
    warning(['Upper bound on regularization parameter smaller' ...
	     ' than lower bound.'])
    h_opt  = h_min; 
  end

endfunction



function [B, S, h, peff] = mridge(Cxx, Cyy, Cxy, dof, options);
%MRIDGE  Multiple ridge regression with generalized cross-validation.
%
%   [B, S, h, peff] = MRIDGE(Cxx, Cyy, Cxy, dof, OPTIONS) returns a
%   regularized estimate B = Mxx_h Cxy of the coefficient matrix in
%   the multivariate multiple regression model Y = X*B + noise(S). The
%   matrix Mxx_h is the regularized inverse of the covariance matrix
%   Cxx,
%
%             Mxx_h = inv(Cxx + h^2 * I).
%
%   The matrix Cxx is an estimate of the covariance matrix of the
%   independent variables X, Cyy is an estimate of the covariance
%   matrix of the dependent variables Y, and Cxy is an estimate of the
%   cross-covariance matrix of the independent variables X and the
%   dependent variables Y. The scalar dof is the number of degrees of
%   freedom that were available for the estimation of the covariance
%   matrices.
%
%   The input structure OPTIONS contains optional parameters for the
%   algorithm:
%
%     Field name         Parameter                                   Default
%
%     OPTIONS.regpar     Regularization parameter h. If regpar       not set
%                        is set, the scalar OPTIONS.regpar is
%                        taken as the regularization parameter h. 
%                        If OPTIONS.regpar is not set (default), 
%                        the regularization parameter h is selected 
%                        as the minimizer of the generalized 
%                        cross-validation (GCV) function. The output
%                        variable h then contains the selected 
%                        regularization parameter.
%  
%     OPTIONS.relvar_res Minimum relative variance of residuals.       5e-2
%                        From the parameter OPTIONS.relvar_res, a
%                        lower bound for the regularization parameter
%                        is constructed, in order to prevent GCV from
%                        erroneously choosing too small a 
%                        regularization parameter (see GCVRIDGE).
%
%   The OPTIONS structure is also passed to GCVRIDGE.
%
%   MRIDGE returns the ridge estimate B of the matrix of regression
%   coefficients. Also returned are an estimate S of the residual
%   covariance matrix, the regularization parameter h, and the scalar
%   peff, an estimate of the effective number of adjustable
%   parameters in B.  
%  
%   MRIDGE computes the estimates of the coefficient matrix and of the
%   residual covariance matrix from the covariance matrices Cxx, Cyy,
%   and Cxy by solving the regularized normal equations. The normal
%   equations are solved via an eigendecomposition of the covariance
%   matrix Cxx. However, if the data matrices X and Y are directly
%   available, a method based on a direct factorization of the data
%   matrices will usually be more efficient and more accurate.
%
%   See also: IRIDGE, GCVRIDGE.

  error(nargchk(4, 5, nargin))     % check number of input arguments 
  
  px           = size(Cxx, 1);
  py           = size(Cyy, 1);
  if size(Cxx, 2) ~= px | size(Cyy, 2) ~= py | any(size(Cxy) ~= [px, py]) 
    error('Incompatible sizes of covariance matrices.')
  end
  
  % ==============           process options        ========================
  if nargin == 4 | isempty(options)
    options    = [];
    fopts      = [];
  else
    fopts      = fieldnames(options);
  end
    
  if strmatch('regpar', fopts)
    h          = options.regpar; 
    h_given    = 1;
  else
    h_given    = 0;
  end
  
  if strmatch('relvar_res', fopts)
    relvar_res = options.relvar_res; 
  else
    relvar_res = 5e-2;
  end
  % =================           end options        =========================
  
  if nargout > 1
    S_out      = 1;
  else
    S_out      = 0;
  end
    
  % eigendecomposition of Cxx
  rmax         = min(dof, px);     % maximum possible rank of Cxx
  [V, d, r]    = peigs(Cxx, rmax);
  
  % Fourier coefficients. (The following expression for the Fourier
  % coefficients is only correct if Cxx = X'*X and Cxy = X'*Y for
  % some, possibly scaled and augmented, data matrices X and Y; for
  % general Cxx and Cxy, all eigenvectors V of Cxx must be included,
  % not just those belonging to nonzero eigenvalues.)
  F            = repmat(ones(r, 1)./sqrt(d), 1, px) .* V' * Cxy;

  % Part of residual covariance matrix that does not depend on the
  % regularization parameter h:
  if (dof > r) 
    S0         = Cyy - F'*F;
  else
    S0         = sparse(py, py);
  end
  
  if ~h_given
    % approximate minimum squared residual
    trSmin     = relvar_res * trace(Cyy);
    
    % find regularization parameter that minimizes the GCV object function
    h          = gcvridge(F, d, trace(S0), dof, r, trSmin, options);
  end

  
  % get matrix of regression coefficients
  B            = V * (repmat(sqrt(d)./(d + h^2), 1, py) .* F);
  
  if S_out
    % get estimate of covariance matrix of residuals
    S          = S0 + F' * (repmat(h^4./(d + h^2).^2, 1, py) .* F);
  end
  
  if nargout == 4
    % effective number of adjusted parameters: peff = trace(Mxx_h Cxx)
    peff       = sum(d ./ (d + h^2));
  end      

endfunction



function [B, S, h, peff] = iridge(Cxx, Cyy, Cxy, dof, options);
%IRIDGE  Individual ridge regressions with generalized cross-validation.
%
%   [B, S, h, peff] = IRIDGE(Cxx, Cyy, Cxy, dof) returns a regularized
%   estimate B of the coefficient matrix for the multivariate multiple
%   regression model Y = X*B + noise(S).  Each column B(:,k) of B is
%   computed by a ridge regression as B(:,k) = Mxx_hk Cxy(:,k), where
%   Mxx_hk is a regularized inverse of Cxx,
%
%             Mxx_h = inv(Cxx + hk^2 * I).
%
%   For each column k of B, an individual regularization parameter
%   ('ridge parameter') hk is selected as the minimizer of the
%   generalized cross-validation function. The matrix Cxx is an
%   estimate of the covariance matrix of the independent variables X,
%   Cyy is an estimate of the covariance matrix of the dependent
%   variables Y, and Cxy is an estimate of the cross-covariance matrix
%   of the independent variables X and the dependent variables Y. The
%   scalar dof is the number of degrees of freedom that were available
%   for the estimation of the covariance matrices.
%
%   The input structure OPTIONS contains optional parameters for the
%   algorithm:
%
%     Field name         Parameter                                   Default
%
%     OPTIONS.relvar_res Minimum relative variance of residuals.       5e-2
%                        From the parameter OPTIONS.relvar_res, a
%                        lower bound for the regularization parameter
%                        is constructed, in order to prevent GCV from
%                        erroneously choosing too small a 
%                        regularization parameter (see GCVRIDGE).
%
%   The OPTIONS structure is also passed to GCVRIDGE.
%     
%   IRIDGE returns an estimate B of the matrix of regression
%   coefficients. Also returned are an estimate S of the residual
%   covariance matrix, a vector h containing the regularization
%   parameters hk for the columns of B, and the scalar peff, an
%   estimate of the effective number of adjustable parameters in each
%   column of B.
% 
%   IRIDGE computes the estimates of the coefficient matrix and of the
%   residual covariance matrix from the covariance matrices Cxx, Cyy,
%   and Cxy by solving the regularized normal equations. The normal
%   equations are solved via an eigendecomposition of the covariance
%   matrix Cxx. However, if the data matrices X and Y are directly
%   available, a method based on a direct factorization of the data
%   matrices will usually be more efficient and more accurate.
%
%   See also: MRIDGE, GCVRIDGE.
  
  error(nargchk(4, 5, nargin))     % check number of input arguments 
  
  px           = size(Cxx, 1);
  py           = size(Cyy, 1);
  if size(Cxx, 2) ~= px | size(Cyy, 2) ~= py | any(size(Cxy) ~= [px, py]) 
    error('Incompatible sizes of covariance matrices.')
  end

  % ==============           process options        ========================
  if nargin < 5 | isempty(options)
    options    = [];
    fopts      = [];
  else
    fopts      = fieldnames(options);
  end
    
  if strmatch('relvar_res', fopts)
    relvar_res = options.relvar_res; 
  else
    relvar_res = 5e-2;
  end
  
  % =================           end options        =========================

  if nargout > 1
    S_out      = 1==1;
  else
    S_out      = 0==1;
  end
  
  if nargout == 4
    peff_out   = 1==1;
  else
    peff_out   = 0==1;
  end
    
  % eigendecomposition of Cxx
  rmax         = min(dof, px);     % maximum possible rank of Cxx
  [V, d, r]    = peigs(Cxx, rmax);
  
  % Fourier coefficients. (The following expression for the Fourier
  % coefficients is only correct if Cxx = X'*X and Cxy = X'*Y for
  % some, possibly scaled and augmented, data matrices X and Y; for
  % general Cxx and Cxy, all eigenvectors V of Cxx must be included,
  % not just those belonging to nonzero eigenvalues.)
  F            = repmat(ones(r, 1)./sqrt(d), 1, px) .* V' * Cxy;

  % Part of residual covariance matrix that does not depend on the
  % regularization parameter h:
  if (dof > r) 
    S0         = Cyy - F'*F;
  else
    S0         = sparse(py, py);
  end
    
  % approximate minimum squared residual
  trSmin       = relvar_res * diag(Cyy);

  % initialize output
  h            = zeros(py, 1);
  B            = zeros(px, py);
  if S_out
    S          = zeros(py, py);
  end
  if peff_out
    peff       = zeros(py, 1);
  end
  
  for k = 1:py                    
    % compute an individual ridge regression for each y-variable
    
    % find regularization parameter that minimizes the GCV object function
    h(k)       = gcvridge(F(:,k), d, S0(k,k), dof, r, trSmin(k), options);
    
    % k-th column of matrix of regression coefficients
    B(:,k)     = V * (sqrt(d)./(d + h(k)^2) .* F(:,k));

    if S_out
      % assemble estimate of covariance matrix of residuals
      for j = 1:k
         diagS  = ( h(j)^2 * h(k)^2 ) ./ ( (d + h(j)^2) .* (d + h(k)^2) );
         S(j,k) = S0(j,k) + F(:,j)' * (diagS .* F(:,k));
         S(k,j) = S(j,k);
      end
    end

    if peff_out 
      % effective number of adjusted parameters in this column
      % of B: peff = trace(Mxx_h Cxx)
      peff(k)  = sum(d ./ (d + h(k)^2));
    end
    
  end

endfunction

