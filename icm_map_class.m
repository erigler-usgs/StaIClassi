%% 
%% Usage:
%%
%% [I_map, P_post, P_cond, P_prior] = icm_map_class (I_chan, I_class, ...
%%                                                   M_class, C_class, W_class, ...
%%                                                   alpha, beta, niter)
%%
%% ICM == "Iterated Conditional Modes"
%%
%% This maximum a posteriori (MAP) multichannel pixel classifier serves as 
%% a wrapper to pixel classification/segmentation functions ml_class.m and 
%% smooth_class.m.  It is a simple iterative algorithm that has been shown
%% to converge to a local maximum in posterior probability space relatively
%% quickly (usually < 10 iterations; Tso & Mather, Classification Methods for
%% Remotely Sensed Data 2nd Ed., CRC Press, 2009).
%%
%%
%% INPUTS:
%%
%% The required input I_chan can be one of two data types, depending on the
%% subsequent input arguments:
%% a) multi-channel 2D image, with channels sorted along the 3rd dimension;
%% b) multi-class probability matrix, with classes sorted along 3rd dimension;
%%
%% If the optional input I_class is passed, use it to seed the ICM algorithm
%% A valid I_class is no more than 2D, and has identical x-y dimensions to
%% I_chan (the only exception is an empty I_class array, in which case the
%% ICM algorithm is seeded with the thematic map that maximizes passed or
%% calculated multi-class conditional probabilities).
%% 
%% The optional M_class, C_class, and W_class, are the multi-channel class 
%% conditioned mean vectors, covariance matrices, and GMM component weight 
%% arrays, respectively. If one is passed, they must all three be passed. 
%% If none are passed, or they are all empty matrices, I_chan will be 
%% treated as a multi-class probability matrix.
%%
%% The optional input alpha is passed to smooth_class.m and holds "clique 1"
%% coefficients used to specify each pixel's class-specific prior.
%%
%% The optional input beta is passed to smooth_class.m and holds "clique 2"
%% coefficient used to specify each pixel's contextual prior.
%%
%% The optional input niter specifies a maximum number of iterations before
%% returning a solution.  If the solution converges first (i.e., the thematic
%% image does not change over two iterations), the algorithm will return.  If
%% niter is not set, it will default to 10 iterations.
%%
%%
%% OUTPUTS:
%%
%% The 2D output matrix I_map is the MAP determination of the most probable 
%% classes for each pixel.
%%
%% The 3D output matrix P_post holds posterior probabilities for each class 
%% and pixel.
%%
%% The 3D output matrix P_cond holds class-conditional probabilities for each
%% class and pixel...this would simply replicate I_chan if I_chan were used
%% to pass the conditional probabilities instead of calculating them with
%% the multichannel observations and GMM parameters.
%%
%% The 3D output matrix P_prior holds contextual prior probabilities derived 
%% from each pixel's nearest neighbors for the final iteration of the ICM.
%% 
%%
%% CONTINUATION WHEN SOLUTION DOESN'T FULLY CONVERGE
%%
%% The ICM algorithm is, mathematically, guaranteed to converge. However, 
%% the numerical reality is that this convergence is dependent on certain
%% assumptions that may not be met, and even if they are, convergence may
%% be slower than anticipated. For these reasons this algorithm uses a 
%% fixed number of iterations instead of a formal convergence criterion.
%% 
%% We described above how the algorithm will exit if full convergence is
%% achieved before niter iterations. However, if the user does not believe
%% the algorithm has achieved an adequate level of convergence after niter
%% iterations, they can pick up exactly where the algorithm left off. For
%% example:
%%
%% octave> rand("state", 1); 
%% octave> [I_map, P_post] = icm_map_class(P_cond, [], [], [], [], 0, 1, 10);
%%
%% ...will give identical end results to:
%%
%% octave> rand("state", 1); 
%% octave> [I_map] = icm_map_class(P_cond, [], [], [], [], 0, 1, 4);
%% octave> [I_map, P_post] = icm_map_class(P_cond, I_map, [], [], [], 0, 1, 6);
%%
%% (rand() simply resets the random number generator to ensure identical results)
%%

%%
%% REFERRENCES:
%%
%% Besag, J. (1986), On the statistical analysis of dirty pictures, Journal
%%   of the Royal Statisical Society, v. 48.3, 259-302.
%% Li, S. Z. (2009), Mathematical MRF models, in Markov Random Field Modeling 
%%   in Image Analysis, 3rd ed., p. 21-48, London.
%% Rigler, E. J., Hill, S. M., Reinard, A. A., and Steenburgh, R. A. (2012),
%%   Solar thematic maps for space weather operations, Space Weather, v. 10, 
%%   online.
%% Tso, B., and Mather, P. (2009), Classification Methods for Remotely Sensed
%%   Data, 2nd Ed., CRC Press.
%% 
%%
%% AUTHOR(S):    E. Joshua Rigler (EJR)
%%
%%
%% CHANGELOG:    17 Jun, 2011 - (EJR) First public release
%%               15 Nov, 2011 - (EJR) Modified to permit initial posterior probabilities
%%                              and explicit conditional probabilities, effectively
%%                              allowing the algorithm to pick up where it left off.
%%               09 Dec, 2011 - (EJR) minor updates to inline documentation
%%               05 Nov, 2012 - (EJR) major updates:
%%                              * no longer allows 3D I_class initial thematic map,
%%                                which makes icm_map_class with other classification
%%                                algorithms like ml_class and gmix_cluster now;
%%                              * no longer calls train_class if training data matrix
%%                                passed in instead of I_class...this was completely
%%                                useless, and led to much code-bloat;
%%                              * no longer interprets I_class as posteriors just
%%                                because it is a floating point array...this was 
%%                                done originally because it seemed necessary when
%%                                the algorithm needed to be restarted, however as
%%                                noted in the help docs, this can be accomplished
%%                                without such unnecessary ambiguity...not sure how
%%                                this ever happened in the first place;
%%                              * now accomodates GMMs, which required adding a new
%%                                input argument W_class;
%%                              * lots of streamlining of pre-processing code and 
%%                                many help documentaion updates.
%%
%%

function [I_map, P_post_all, P_cond_all, P_prior_all] = icm_map_class (I_chan, I_class, ...
                                                                      M_class, C_class, W_class, ...
                                                                      alpha, beta, niter)
   
   %
   % process inputs first
   %
   
   % ensure adequate number of inputs
   if nargin < 1
      error('icm_map_class: at least one input argument is required');
   endif
   
   
   % at least one input
   if nargin > 0
      
      % create default values if parameters not passed
            
      % check existance of I_class, M_class, C_class and W_class inputs
      if ~exist('I_class','var')
         I_class = [];
      endif
      
      if ~exist('M_class','var')
         M_class = [];
      endif
      
      if ~exist('C_class','var')
         C_class = [];
      endif
      
      if ~exist('W_class','var')
         W_class = [];
      endif
      
      
      
      % check existance of alpha, beta, and niter
      if ~exist('alpha','var')
         alpha = []; 
      endif
      
      if ~exist('beta','var')
         beta = [];
      endif
      
      if ~exist('niter', 'var')
         niter = [];
      endif
      
      
   endif
   
   
   % at least 2 inputs
   if nargin > 1
      
      if numel(I_class) ~= 0
         
         % make sure I_class is ONLY 2D
         if ndims(I_class) > 2
            error('icm_map_class: I_class must be nor more than 2D')
         endif
         
         % make sure I_class matches I_chan in first two dimensions
         if size(I_chan,1) ~= size(I_class,1) || ...
            size(I_chan,2) ~= size(I_class,2)
            error('icm_map_class: XY dimension mismatch between I_chan and I_class');
         endif
         
      endif
      
   endif
   
   
   % at least 3-5 inputs
   if nargin > 2
      
      % actually, if M_class is passed, C_class and W_class must be passed too, and be
      % dimensionally consistent with each other
      if ~(nargin > 4)
         error('icm_map_class: M_class, C_class, and W_class must be passed together');
      endif
      
      % determine various GMM distribution parameters from M_class
      nchan = size(M_class,2) * ~isempty(M_class);
      nclass = size(M_class,3) * ~isempty(M_class);
      ncomp = size(M_class,4) * ~isempty(M_class);
      
      % check C_class dimensions
      if size(C_class,1) * ~isempty(C_class) ~= nchan || ...
         size(C_class,2) * ~isempty(C_class) ~= nchan || ...
         size(C_class,3) * ~isempty(C_class) ~= nclass || ...
         size(C_class,4) * ~isempty(C_class) ~= ncomp 
         error('icm_map_class: C_class dimensions are not valid or are inconsistent with M_class');
      endif
      
      % check W_class dimensions
      if size(W_class,3) * ~isempty(W_class) ~= nclass || ...
         size(W_class,4) * ~isempty(W_class) ~= ncomp
         error('icm_map_class: W_class dimensioons are not valid or are inconsistent with M_class');
      endif
      
      % could check W_class XY dimensions, but this will be done in call to ml_class
            
   endif
   
   
   % at this point, it is possible to generate initial posterior probabilities and initial 
   % thematic maps, if they have not already been passed in
   
   
   % two paths, depending on whether GMM parameters were pass or not
   if isempty(M_class)
      
      % if no GMM parameters were passed, I_chan must be conditional probs
      P_cond_all = sum(I_chan, 4);
      
      % if I_class was empty, assume it simply corresonds to the index of the largest
      % probabilty in P_cond_all
      if isempty(I_class)
      
         [P_tmp, I_post_all] = max(P_cond_all, [], 3);
      
         % if the maximum probability is zero, all probabilities must have been
         % zero, which causes max() to return an index of 1...not desirable; 
         % if maximum probabilty is a NaN, then all probabilities must have
         % been a NaN, which indicates that at least one value in the feature
         % vector was a NaN, and therefore unclassifiable unless imputation
         % was used...if imputation *was* used, then there must not hvae been
         % any valid values in the feature vector.
         I_post_all(P_tmp == 0 | isnan(P_tmp)) = 0;
      
      else
                  
         % I_class's dimensions should have already been checked above
         I_post_all = I_class;
      
      endif
      
      % return only integer-valued indices
      I_post_all = int16(I_post_all);
      
      
      % finally, define nclass and ncomp, since they were not defined previously
      % (nchan has no meaning if I_chan is a conditional probability matrix, just set it to zero)
      nchan = 0;
      nclass = size(P_cond_all,3);
      ncomp = size(P_cond_all,4);
      
   else
      

      % if GMM parameters were passed, I_chan must contain multichannel observations,
      % so start by checking dimensions
      if size(I_chan,3) ~= nchan
         error('icm_map_class: multichannel observations and GMM parameters not dimensionally consisent');
      endif
      
      
      % determine initial thematic labels and conditional probabilities from I_chan and GMM
      % parameters; assume all probs are valid, and use imputation to replace missing obs.
      [I_post_all, P_cond_all] = ml_class(I_chan, M_class, C_class, W_class, 1, 1);
      
      
      % if I_class was NOT empty, us it to initialize I_post_all, even if it differs
      % from what was generated by ml_class
      if ~isempty(I_class)
         I_post_all = I_class;
      endif
      
   endif
      
   
   
   % any remaining inputs are non-default MRF/ICM parameters
   
   if isempty(alpha)
      
      alpha = zeros(nclass,1);
      
   else
      
      % check and/or adjust alpha
      if isscalar(alpha)
         alpha = ones(nclass,1) * alpha;
      elseif numel(alpha) ~= nclass
         error('icm_map_class: alpha must be a scalar, or a vector equal in length to the number of classes');
      endif
         
   endif
   
   if isempty(beta)
      
      beta = [0 0 0 0];
   
   else
      
      % check and/or adjust beta
      if isscalar(beta)
         beta = ones(1,4) * beta;
      elseif numel(beta) ~= 4
         error('icm_map_class: beta must be a scalar, or a 4-element vector');
      endif
      
   endif
   
   % at least 8 inputs
   if isempty(niter)
      
      % set default niter
      niter = 10;
      
   endif
   
   
   % ensure that no more inputs were passed than can be used
   if nargin > 8
      error('icm_map_class: too many input arguments passed to funtion');
   endif
   
   
   %
   % begin actual ICM algorithm
   %
   
   % save current copy of I_post_all for later comparisons
   I_map = I_post_all;
   
   % initialize stuff
   P_post_all = P_cond_all / nclass;
   
   % begin ICM loop
   for k=1:niter
      
      % generate whole-image contextual prior probabilities from previous
      % posterior probabilities solution...ICM does not use the labels generated
      % by smooth_class, just the probabilties
      [~, P_prior_all] = smooth_class (I_post_all, alpha, beta, 1);

      
      % combine class-conditional and prior probabilities using Bayes' theorem
      P_post_all = P_cond_all .* P_prior_all;
      
      
      % extract MAP class IDs from P_post_all and compress to 2D for comparison
      [P_max, I_post_all] = max(P_post_all, [], 3);
      
      % convert to integer array
      I_post_all = int16(I_post_all);

      % if the maximum probability is zero, this just means that the pixel was not
      % classified before being passed in, or none of its neighbors were classified
      % before being passed in; either way, its label should remain zero;
      % if the max probability is a NaN, this means that the multichannel pixel
      % must have contained at least one NaN element, so it's label should be zero
      I_post_all(P_max == 0 | isnan(P_max)) = 0;
      
      
      % check if I_post_all differs from I_map
      if (any(I_post_all(:) - I_map(:)))
                  
         % DEBUG DEBUG DEBUG
         disp(sum([I_post_all(:) - I_map(:)]~=0));
         fflush(stdout);
                  
#          cmap_tmp = colormap;
#          figure(k);
#          colormap(cmap_tmp);
#          
#          imagesc(sum(I_post_2d, 3));
#          axis('xy','equal');
#          caxis([-.5 8.5]);
#          %colorbar;
#          %pause;
         
         
         % update I_map 
         I_map = I_post_all;
         
         
         % skip to end of loop
         continue;
         
      else
         
         % no need to update I_map, just break out of loop and return I_map array
         break;
         
      endif
            
            
   end % endfor k=1:niter
      
   
   % finally, copy I_post_all to I_map...this is only because older versions of this
   % function specified I_map as the output, and I didn't want to confuse things
   I_map = I_post_all;
   
   
end % endfunction
