%%
%% Usage:
%%
%% [I_smooth, P_smooth] = smooth_class (I_class, alpha[, beta[, sets]])
%%
%% This function "smooths" the thematic image I_class by maximizing
%% contextual probabilities determined from I_class pixels and their
%% nearest neighbors using the input parameters alpha and beta.
%%
%%
%% INPUTS:
%%
%% The required input I_class is a 2D integer array comprised of class
%% IDs similar to the I_class output of the function ml_class.m.
%%
%% The required input alpha is typicall a vector equal in length to the
%% number of classes, and provides Markov Random Field (MRF) "clique 1" 
%% coeffient(s). This means that a prior probability of Exp(alpha) is 
%% assigned to each class when determining each pixel's posterior  
%% probability.  Alpha may be a scalar integer specifying the number
%% of classes if all classes are to receive equal weight. If the user
%% wishes to assign equal weights that differ from exp(0)==1 (a strange
%% wish, to be sure), they must specify alpha in its vector form.
%%
%% The optional input parameter beta provides MRF "clique 2" coeffient(s).
%% A probability proportional to Exp(beta*nj + beta'*nj' + beta''*nj'' + 
%%  beta'''*nj'''), where nj* is the number of a pixel's nearest neighbors 
%% in the left-right, up-down, down-diagonal, and up-diagonal orientations
%% with a class ID that matches the given pixel, is applied to each class 
%% when determining probabilities of class membership for each pixel.
%%
%% NOTE: the equations coded up in this function assume positive alpha means
%%       weights greater than 1, and positive beta leads to smoother results. 
%%       Negative values are allowed however, implying weights between 0 and 1
%%       for alpha, and something less "smooth" for beta.
%%
%% FIXME_00: we smartly deal with situations where multiple classes get assigned
%%           the maximum probability (this is more common than one might think) 
%%           by randomly assigning one of those classes to the pixel. HOWEVER, 
%%           it would probably be best to leave a pixel label alone if it is 
%%           one of the max probability classes already. This is a relatively 
%%           minor problem in search of a creative solution.
%% FIXME_01: a useful capability would be to specify alpha and beta as XY 
%%           arrays to allow spatially-dependent contextual parameters a la 
%%           Turmon 2010, but this is a fairly major undertaking.
%% FIXME_02: the probabilities generated here are only properly normalized if
%%           beta is a scalar, or all betas are identical (i.e., isotropy is
%%            assumed). This does not adversely impact the smoother, or even 
%%           iterative MAP solvers like ICM, since these are only concerned with
%%           relative probabilities. HOWEVER, this function cannot be used to
%%           iteratively estimate optimal beta parameters unless a more reliable
%%           technique is used to estimate a appropriate normalization constant.
%% FIXME_03: it seems, after reading Besag-1986 ("On the Statistical Analysis of  
%%            Dirty Pictures"), that calculating normalization constants is 
%%           essentially impossible if beta is not isotropic, and maybe even
%%           impossible if it is isotropic, but alpha is not...therefore, we
%%           may as well give up on normalizing these probabilities alltogether,
%%           save the compuation time, and allow ourselves the flexibility to
%%           use not just spatially non-isotropic betas, but also betas that 
%%           describe inter-class relationships (i.e., neighbors of more than 
%%           one pixel type may contribute to probability of class membership,
%%           or two pixel types may have a negative correlation...e.g., filaments
%%           and coronal holes should not exist in adjascent pixels).
%%
%% The optional input parameter sets defines which "coding set" should be 
%% processed. A coding set (see Besag-1986, or even go back to Besag-1974)
%% is a subset of the image where no two pixels are mutual neighbors during a
%% single pass of the algorithm through the set. The most extreme "coding set"
%% is one in which each pixel is passed through the algorithm, adjusted, then
%% the next pixel is processed. This continues sequentially until the last 
%% pixel is processed, and is VERY computationally inefficient. 
%% 
%% A compromise, optionally implemented here, is a simple 4-pixel cycle of 
%% coding sets, where the odd-odd indixed pixels are updated first, the odd-
%% even next, the even-odd after that, and the even-even to finish up. If sets 
%% is a scalar that evaluates true, process all 4 coding sets; if sets is a 
%% 4-vector of scalars, process the coding set if the corresponding element in 
%% sets is true; if sets evaluates false, or is not passed at all, perform a 
%% fully synchronous update of the pixels.
%%
%%
%% OUTPUTS:
%%
%% The output I_smooth is a statistically "smoothed" version of the input 
%% I_class.  In all likelihood, it will NOT be this value that is desired
%% by the user, but rather the probabilities that can be combined with
%% other probabilities in some type of Bayesian pixel classificaiton.
%%
%% The output P_smooth hold the probabilities that each pixel belongs to
%% a class corresponding to the index along the third dimension of P_smooth.
%%

%%
%% REFERENCES:
%%
%% Besag, J. (1986), On the statistical analysis of dirty pictures, Journal
%%   of the Royal Statisical Society, v. 48.3, 259-302.
%% Li, S. Z. (2009), Mathematical MRF models, in Markov Random Field Modeling 
%%   in Image Analysis, 3rd ed., p. 21-48, London.
%% Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., an dTeller, E. (1953),
%%    Equations of state calculations by fast computing machines, Journal of
%%    Chemical Physics, v. 21.6, p. 1087-1092.
%% Rigler, E. J., Hill, S. M., Reinard, A. A., and Steenburgh, R. A. (2012),
%%   Solar thematic maps for space weather operations, Space Weather, v. 10, 
%%   online.
%% Tso, B., and Mather, P. (2009), Classification Methods for Remotely Sensed
%%   Data, 2nd Ed., CRC Press.
%% Turmon, M., Jones, H. P., Malanushenko, O. V., and Pap, J. M. (2010), 
%%   Statistical feature recognition for multidimensional solar imagery, Solar
%%   Physics, v. 262, p. 277-298. (see also Turmon et al. 2002, cited within)%% 
%%
%%
%% AUTHOR(S):     E. Joshua Rigler (EJR)
%%
%%
%% CHANGELOG:     17 Jun, 2011 - (EJR) First public release
%%                04 Dec, 2011 - (EJR) Modified to accept 4-valued beta, even though the
%%                               normalization of probabilities is not yet correct
%%                               if these values are not identical;
%%                               modified to consider non-matching neighbors in
%%                               addition to matching neighbors when calculating
%%                               probability of class membership
%%                09 Dec, 2011 - (EJR) modified to process synchronously (default), or
%%                               semi-synchronously using four "coding sets"; the
%%                               latter is a necessary compromise between a batch
%%                               operation and the pixel-by-pixel processing that
%%                               is typically required to guarantee convergence in
%%                               iterative solvers like ICM or even SA.
%%                10 Apr. 2012 - (EJR) fixed minor indexing bug in coding set #1 that
%%                               choked on arrays with odd-sized dimensions
%%               ~16 May. 2012 - (EJR) now treats invalid labels (i.e., any label that
%%                               is not an integer in the set [0:nclass]) as "missing",
%%                               and replaces it with the most likely label given its
%%                               nearest neighbors; as before, any zero-valued label
%%                               that is passed in will remain a zero-valued label.
%%                05 Nov, 2012 - (EJR) major updates:
%%                               * I_class must now be a 2D array, and all the crazy
%%                                 code-bloat associated with allowing it to be a
%%                                 3D training pixel mask, or a 3D multi-class prob.
%%                                 matrix is gone;
%%                               * scalar alpha is not used to specify the number of
%%                                 classes that will all receive an equal weight of
%%                                 exp(0)==1;
%%                               * substantial cleaning up of documentation.
%%
%%

function [I_smooth, P_smooth] = smooth_class (I_class, alpha, beta, sets)
   
   %
   % Process input arguments
   %
   
   % make sure there are enough inputs
   if (nargin < 2)
      error('smooth_class: At least two input arguments are required.');
   endif
   
   
   % make sure there are not too many inputs
   if (nargin > 4)
      error('smooth_class: No more than four inputs argument are allowed');
   endif
   
      
   % if at least 2 inputs
   if nargin > 1
      
      % check I_class
      if ndims(I_class) > 2
         error('smooth_class: I_class must be a 2D array');
      endif
      
      [ny, nx] = size(I_class);
      
      
      % check and adjust alpha
      if isscalar(alpha)
         
         % it doesn't make much sense to assign identical alphas to each class, since
         % this just weights all classes equally, so use scalar alpha to define the
         % number of classes
         nclass = alpha;
         alpha = zeros(nclass, 1);
      
      elseif isvector(alpha)
         
         nclass = numel(alpha);
      
      else
         
         error('smooth_class: alpha must be a scalar or a vector equal in length to the number of classes');
      
      endif
      
   endif
   
   
   % if at least 3 inputs
   if nargin > 2
      
      if (isscalar(beta))
         
         beta = [beta, beta, beta, beta];
      
      elseif ~(isvector(beta) && length(beta) == 4)
         
         % only accept scalar or vector-4 values for beta
         error('Beta must be a scalar or 4 element vector');
      
      endif
   
   else
      
      % if beta was not passed in, set it equal to zero
      beta = [0,0,0,0];
   
   endif

   % at least 4 inputs...if passed, make sure sets is usable
   if nargin > 3
      
      if isscalar(sets)
         
         if (sets)
            
            sets = [1 1 1 1];
         
         else
            
            sets = [0 0 0 0];
         
         endif
      
      else
         
         if ~(isvector(sets) && numel(sets) == 4)
            
            error('smooth_class.m: input sets must be a 4-vector of logical scalars');
         
         endif
      
      endif
   
   else
      
      % default to fully synchronous update
      sets = [0 0 0 0];
   
   endif
   
     
   
   
#    % check class, structure, and sizes of input arguments
#    if (isinteger(I_class))
#                
#       if (size(I_class,3) == 1)
#          % if 2D integer array is passed as input, we need another way to determine
#          % the number of classes
#          nclass = numel(alpha);
#          
#          if (nclass <= 1 && numel(unique(I_class)) > 1)
#             % alpha either wasn't passed, or it is a scalar, yet more than 1 class 
#             % ID exists in the 2D map, so there is no reliable information to 
#             % determine the number of possible class IDs...throw an error.
#             error('Number of possible class IDs cannot be determined');
#          endif
#          
#          [ny, nx] = size(I_class);
#          I_smooth = I_class;
#          
#       else
#          
#          [ny,nx,nclass] = size(I_class);
#          
#          % convert alpha into vector
#          if (numel(alpha) == 0) 
#             alpha = zeros(nclass, 1);
#          elseif (numel(alpha) == 1)
#             alpha = zeros(nclass, 1) + alpha;
#          endif
#          
#          % compress 3D integer array into a 2D integer array using bin_class.m
#          % to decide between potentially overlapping classes
#          I_smooth = bin_class(I_class, 1, alpha);
#          I_smooth = max(I_smooth, [], 3);
#          
#       endif
#       
#    elseif (isfloat(I_class))
#       
#       % determine the image dimensions and number of class IDs
#       [ny,nx,nclass] = size(I_class);
#                
#       % convert alpha into vector
#       if (numel(alpha) == 0) 
#          alpha = zeros(nclass, 1);
#       elseif (numel(alpha) == 1)
#          alpha = zeros(nclass, 1) + alpha;
#       endif
#       
#       % generate a 2D integer array of class IDs from the index along
#       % the third dimension of I_class corresponding to the highest
#       % probability for that pixel
#       [P_tmp, I_smooth] = max(I_class, [], 3);
#                
#       % convert I_tmp to integer array
#       I_smooth = int16(I_smooth);
#       
#       % if the maximum probability is zero, this just means that the pixel was not
#       % classified before being passed in, so its label should remain zero
#       I_smooth(P_tmp == 0) = 0;
#       
#    else
#       error('I_class is not a valid data type');
#    endif      
# 
#    
#    % if passed, make sure beta is usable
#    if (nargin >= 3)
#       if (isscalar(beta))
#          beta = [beta, beta, beta, beta];
#       elseif ~(isvector(beta) && length(beta) == 4)
#          % only accept scalar or vector-4 values for beta
#          error('Beta must be a scalar or 4 element vector');
#       endif
#    else
#       % if beta was not passed in, set it equal to zero
#       beta = [0,0,0,0];
#    endif
#    
#    
#     % if passed, make sure sets is usable
#     if (nargin == 4)
#       if (isscalar(sets))
#         if (sets)
#           sets = [1 1 1 1];
#         else
#           sets = [0 0 0 0];
#         endif
#       else
#         if ~(isvector(sets) && numel(sets) == 4)
#           error('smooth_class.m: input sets must be a 4-vector of logical scalars');
#         endif
#       endif
#     else
#       % default to fully synchronous update
#       sets = [0 0 0 0];
#     endif
#    
  
      
      
# #
# # The following commented block of code is left intact so that comparisons with
# # the vectorized version below it can be made to confirm the two are equivalent.
# # It should not be used regularly because it is more than 100 times slower.
# # 
#    % loop over each pixel in the thematic map image
#    sy = size(I_smooth,1);
#    sx = size(I_smooth,2);
#    cnt = 0;
#    cnt_end = sy*sx;
#    
#    for ry = 1:sy
#    for rx = 1:sx
#       
#       if I_smooth(ry,rx) == 0
#          
#          % assign zero-probability to all classes if this pixel started out unclassified
#          P_smooth(ry,rx,:) = 0;
#          
#       else
#       
#          % determine if this pixel falls along one of the top, bottom, left, or
#          % right-most edge of the array
#          top = ry == sy;
#          bottom = ry == 1; 
#          left = rx == 1;
#          right = rx == sx;
#          
#          
#          % loop over each class for this pixel
#          for j = 1:nclass
#             
#             nj = 0;
#             
#             
#             % if Beta==0, there is no point in counting up matching nearest neighbors, so
#             % just skip the next if-else block
#             if (beta ~= 0)
#             
#                % determine the number of this pixel's nearest neighbors with class IDs 
#                % that match j (this if-then block is only necessary because pixels that
#                % fall along edges or in corners must be treated specially)
#                if (top && left)
#                   
#                   % this pixel is the upper-left corner and has 3 neighbors
#                   nj = (I_smooth(ry-1,rx) == j) + (I_smooth(ry-1,rx+1) == j) + (I_smooth(ry,rx+1) == j);
#                            
#                elseif (top && right)
#                   
#                   % this pixel is the upper-right corner and has 3 neighbors
#                   nj = (I_smooth(ry-1,rx) == j) + (I_smooth(ry-1,rx-1) == j) + (I_smooth(ry,rx-1) == j);
#             
#                elseif (bottom && left)
#                   
#                   % this pixel is the lower-left corner and has 3 neighbors
#                   nj = (I_smooth(ry,rx+1) == j) + (I_smooth(ry+1,rx+1) == j) + (I_smooth(ry+1,rx) == j);
#                
#                elseif (bottom && right)
#                   
#                   % this pixel is the lower-right corner and has 3 neighbors
#                   nj = (I_smooth(ry,rx-1) == j) + (I_smooth(ry+1,rx-1) == j) + (I_smooth(ry+1,rx) == j);
#             
#                elseif (top)
#                   
#                   % this pixel is along the top edge and has 5 neighbors
#                   nj = (I_smooth(ry,rx-1) == j)            +                  (I_smooth(ry,rx+1) == j) + ...
#                      (I_smooth(ry-1,rx-1) == j) + (I_smooth(ry-1,rx) == j) + (I_smooth(ry-1,rx+1) == j);
#                
#                elseif (bottom)
#                   
#                   % this pixel is along the bottom edge and has 5 neighbors
#                   nj = (I_smooth(ry+1,rx-1) == j) +  (I_smooth(ry+1,rx) == j) + (I_smooth(ry+1,rx+1) == j) + ...
#                      (I_smooth(ry,rx-1) == j)              +                 (I_smooth(ry,rx+1) == j);
#                               
#                elseif (left)
#                
#                   % this pixel is along the left edge and has 5 neighbors
#                   nj = (I_smooth(ry+1,rx) == j) +  (I_smooth(ry+1,rx+1) == j) + ...
#                                                 (I_smooth(ry,rx+1) == j) + ...
#                      (I_smooth(ry-1,rx) == j) +  (I_smooth(ry-1,rx+1) == j);
#                
#                elseif (right)
#                
#                   % this pixel is along the right edge and has 5 neighbors
#                   nj = (I_smooth(ry+1,rx) == j) + (I_smooth(ry+1,rx-1) == j) + ...
#                      (I_smooth(ry,rx-1) == j) +  ...
#                      (I_smooth(ry-1,rx-1) == j) + (I_smooth(ry-1,rx) == j);
#                
#                else
#                   
#                   % this pixel is not special in any way and has 8 neighbors
#                   nj = (I_smooth(ry-1,rx-1) == j) +  (I_smooth(ry-1,rx) == j) + (I_smooth(ry-1,rx+1) == j) + ...
#                      (I_smooth(ry,rx-1) == j)              +                 (I_smooth(ry,rx+1) == j) + ...
#                      (I_smooth(ry+1,rx-1) == j) +  (I_smooth(ry+1,rx) == j) + (I_smooth(ry+1,rx+1) == j);
#                   
#                endif
#             
#             endif
#             
#             
#             % calcuate the un-normalized "smoothness" prior probability using alpha, beta, nj
# #            nj = -nj; % this is necessary to make (beta > 0) imply stronger smoothing
#             energy = -alpha(j) - beta * (nj);            
#             P_smooth(ry,rx,j) = exp(-energy);
#             
#          end % endfor j=1:nclass
#          
#          
#          % generate normalizing constant for the "smoothness" prior
#          P_norm = sum(P_smooth(ry,rx,:),3);
#          
#          % normalize the "smoothness" prior
#          P_smooth(ry,rx,:) = P_smooth(ry,rx,:) ./ P_norm;
#       
#       end % endif I_smooth(ry,rx) == 0
#       
#       
#       % update I_smooth now, or convergence is not guaranteed if being called iteratively (see Besag, 1986, 
#       %  On the Statistical Analysis of Dirty Pictures, J. R. Statist. Soc., Vol. 48, No. 3); note that
#       % though this fully asynchronous approach should converge, it may not give exactly the same answer 
#       % as obtained when partially synchronous techniques using "coding sets" are used used.
#       [P_tmp, I_smooth_tmp] = max(P_smooth(ry, rx, :), [], 3);
#       I_smooth(ry,rx) = int16(I_smooth_tmp);
#       
#    
#       % print 50 dots
#       cnt = cnt + 1;
#       if (~mod(cnt,round(cnt_end/50)))
#          printf('.');
#          fflush(stdout);
#       endif
# 
#    
#    end % endfor rx=1:sx
# 
#    end % endfor ry=1:sy



   % this loop and ugly-looking bit of vectorized code should avoid the big loop over 
   % all pixels commented out above, but still give the same result
   % WOOHOO!!! Doing this sped things up by >~300, and the probabilities were
   %           identical to the old method...not bad
   % BOOHOO!!! Fully synchronous updates of the map violate some conditions that guarantee
   %           convergence in iterative solvers. We need to break this up into "coding sets"
   %           according to Besag-1986, where carefully chosen subsets of the input map are 
   %           updated, then subsequent subsets are updated. No member of a subset may be a
   %           neighbor to another member of the same subset. Therefore, we will lose some 
   %           hard-won performance increase, but hopefully not too much
   
      
   
   % initialize I_smooth and P_smooth
   I_smooth = I_class;
   P_smooth = zeros(ny,nx,nclass);
   
   
   %
   % Coding Set #1:
   %
   if (sets(1))

   % create 8 neighbor pixel arrays for every other pixel (note that these neighbor_* arrays
   % are constructed with 4x the necessary pixels because I could not figure out how to 
   % efficiently construct the neighbor_* array with the coding set's dimensions...the
   % coding set is extracted subsequently to generate the counts for each couple of pair-
   % wise cliques)
   neighbor_lowerleft = 0 * I_smooth;
   neighbor_lowerleft(2:end, 2:end) = I_smooth(1:end-1, 1:end-1); 
   
   neighbor_below = 0 * I_smooth;
   neighbor_below(2:end,:) = I_smooth(1:end-1, :);
   
   neighbor_lowerright = 0 * I_smooth;
   neighbor_lowerright(2:end, 1:end-1) = I_smooth(1:end-1, 2:end); 
   
   neighbor_right = 0 * I_smooth;
   neighbor_right(:, 1:end-1) = I_smooth(:, 2:end);
   
   neighbor_upperright = 0 * I_smooth;
   neighbor_upperright(1:end-1, 1:end-1) = I_smooth(2:end, 2:end);
   
   neighbor_above = 0 * I_smooth;
   neighbor_above(1:end-1,:) = I_smooth(2:end, :);
   
   neighbor_upperleft = 0 * I_smooth;
   neighbor_upperleft(1:end-1, 2:end) = I_smooth(2:end, 1:end-1);
   
   neighbor_left = 0 * I_smooth;
   neighbor_left(:, 2:end) = I_smooth(:, 1:end-1);
      
   
   % cycle through each class and sum up those that match class j
   n1 = [];
   n2 = [];
   n3 = [];
   n4 = [];
   for j=1:nclass
   
      n1(:,:,j) = (neighbor_right(1:2:end,1:2:end) == j) + (neighbor_left(1:2:end,1:2:end) == j);
      n2(:,:,j) = (neighbor_below(1:2:end,1:2:end) == j) + (neighbor_above(1:2:end,1:2:end) == j);
      n3(:,:,j) = (neighbor_lowerleft(1:2:end,1:2:end) == j)  + (neighbor_upperright(1:2:end,1:2:end) == j);
      n4(:,:,j) = (neighbor_lowerright(1:2:end,1:2:end) == j) + (neighbor_upperleft(1:2:end,1:2:end) == j);
      
   end % for j=1:nclass
   

   % multiply n*Beta and use to calculate Coding Set #1 probabilities and normalize
   P_smooth(1:2:end, 1:2:end,:) = exp(-(repmat(-reshape(alpha, [1,1,numel(alpha)]), ...
                                             [size(I_smooth(1:2:end,1:2:end),1), size(I_smooth(1:2:end,1:2:end),2), 1]) - ...
                                       (n1 * beta(1) - (2-n1)*1 * beta(1) + ...
                                        n2 * beta(2) - (2-n2)*1 * beta(2) + ...
                                        n3 * beta(3) - (2-n3)*1 * beta(3) + ...
                                        n4 * beta(4) - (2-n4)*1 * beta(4) ) ) );
   P_norm = sum(P_smooth(1:2:end,1:2:end,:), 3);
   P_smooth(1:2:end,1:2:end,:) = P_smooth(1:2:end,1:2:end,:) ./ repmat(P_norm, [1, 1, nclass]);
   
   
   % update class IDs for Coding Set #1
   [P_tmp, I_tmp] = max(P_smooth(1:2:end,1:2:end,:), [], 3);
   
   
   % often there are multiple "max" values, which leads to a clear bias toward lower indices due to the
   % way Octave's max() function works; therefore we find all matching maxes between P_smooth and P_tmp,
   % then randomly choose one of the labels matching this probability to break the tie
   I_mask = logical(sum([n1+n2+n3+n4], 3)); % ignore pixels that match none of [1:nclass]
   I_tmp = I_tmp .* I_mask;
   max_match = P_smooth(1:2:end, 1:2:end,:) == repmat(P_tmp, [1, 1, nclass]) & ...
               repmat(I_mask, [1, 1, nclass]); % create mask of matching max probabilities
   max_match_count = sum(max_match,3);
   max_match_float = double(max_match);
   max_match_float(max_match) = rand(sum(max_match(:)),1); % reduce # of calls to random number generator
   [rand_max, I_max] = max(max_match_float, [], 3);
   I_tmp(max_match_count>1) = I_max(max_match_count>1);
   
   
   % if the original I_smooth pixels were euql to zero when passed in, they should equal zero when returned;
   % if the max probability is a NaN, the multichannel pixel used to generate conditional probabilities
   % must have contained at least one NaN element, so the label should be zero...not sure this matters now
   I_tmp(I_smooth(1:2:end,1:2:end) == 0 | isnan(P_tmp)) = 0;
   
   % if the original I_smooth pixel labels were invalid, they should have been replaced with valid labels,
   % unless they had no valid neighbors, in which case they are left alone;
   I_tmp(~I_mask) = I_smooth(~I_mask);
   
   % must convert smoothed map to integer array
   I_smooth(1:2:end,1:2:end) = int16(I_tmp);
   
   endif % if (sets(1))
   

   %
   % Coding Set #2
   %
   if (sets(2))

   % create 8 neighbor pixel arrays for every other pixel (note that these neighbor_* arrays
   % are constructed with 4x the necessary pixels because I could not figure out how to 
   % efficiently construct the neighbor_* array with the coding set's dimensions...the
   % coding set is extracted subsequently to generate the counts for each couple of pair-
   % wise cliques)
   neighbor_lowerleft = 0 * I_smooth;
   neighbor_lowerleft(2:end, 2:end) = I_smooth(1:end-1, 1:end-1); 
   
   neighbor_below = 0 * I_smooth;
   neighbor_below(2:end,:) = I_smooth(1:end-1, :);
   
   neighbor_lowerright = 0 * I_smooth;
   neighbor_lowerright(2:end, 1:end-1) = I_smooth(1:end-1, 2:end); 
   
   neighbor_right = 0 * I_smooth;
   neighbor_right(:, 1:end-1) = I_smooth(:, 2:end);
   
   neighbor_upperright = 0 * I_smooth;
   neighbor_upperright(1:end-1, 1:end-1) = I_smooth(2:end, 2:end);
   
   neighbor_above = 0 * I_smooth;
   neighbor_above(1:end-1,:) = I_smooth(2:end, :);
   
   neighbor_upperleft = 0 * I_smooth;
   neighbor_upperleft(1:end-1, 2:end) = I_smooth(2:end, 1:end-1);
   
   neighbor_left = 0 * I_smooth;
   neighbor_left(:, 2:end) = I_smooth(:, 1:end-1);
      
   
   % cycle through each class and sum up those that match class j
   n1 = [];
   n2 = [];
   n3 = [];
   n4 = [];
   for j=1:nclass
   
      n1(:,:,j) = (neighbor_right(1:2:end, 2:2:end) == j) + (neighbor_left(1:2:end, 2:2:end) == j);
      n2(:,:,j) = (neighbor_below(1:2:end, 2:2:end) == j) + (neighbor_above(1:2:end, 2:2:end) == j);
      n3(:,:,j) = (neighbor_lowerleft(1:2:end, 2:2:end) == j)  + (neighbor_upperright(1:2:end, 2:2:end) == j);
      n4(:,:,j) = (neighbor_lowerright(1:2:end, 2:2:end) == j) + (neighbor_upperleft(1:2:end, 2:2:end) == j);
      
   end % for j=1:nclass
   

   % multiply n*Beta and use to calculate Coding Set #2 probabilities and normalize
   P_smooth(1:2:end, 2:2:end,:) = exp(-(repmat(-reshape(alpha, [1,1,numel(alpha)]), ...
                                             [size(I_smooth(1:2:end,2:2:end),1), size(I_smooth(1:2:end,2:2:end),2), 1]) - ...
                                       (n1 * beta(1) - (2-n1)*1 * beta(1) + ...
                                        n2 * beta(2) - (2-n2)*1 * beta(2) + ...
                                        n3 * beta(3) - (2-n3)*1 * beta(3) + ...
                                        n4 * beta(4) - (2-n4)*1 * beta(4) ) ) );
   P_norm = sum(P_smooth(1:2:end,2:2:end,:), 3);
   P_smooth(1:2:end,2:2:end,:) = P_smooth(1:2:end,2:2:end,:) ./ repmat(P_norm, [1, 1, nclass]);
   
   
   % update class IDs for Coding Set #2
   [P_tmp, I_tmp] = max(P_smooth(1:2:end,2:2:end,:), [], 3);


   % often there are multiple "max" values, which leads to a clear bias toward lower indices due to the
   % way Octave's max() function works; therefore we find all matching maxes between P_smooth and P_tmp,
   % then randomly choose one of the labels matching this probability to break the tie
   I_mask = logical(sum([n1+n2+n3+n4], 3)); % ignore pixels that match none of [1:nclass]
   I_tmp = I_tmp .* I_mask;
   max_match = P_smooth(1:2:end,2:2:end,:) == repmat(P_tmp, [1, 1, nclass]) & ...
               repmat(I_mask, [1, 1, nclass]); % create mask of matching max probabilities
   max_match_count = sum(max_match,3);
   max_match_float = double(max_match);
   max_match_float(max_match) = rand(sum(max_match(:)),1); % reduce # of calls to random number generator
   [rand_max, I_max] = max(max_match_float, [], 3);
   I_tmp(max_match_count>1) = I_max(max_match_count>1);
   
   
   % if the original I_smooth pixels were equal to zero when passed in, they should equal zero when returned;
   % if the max probability is a NaN, the multichannel pixel used to generate conditional probabilities
   % must have contained at least one NaN element, so the label should be zero...not sure this matters now
   I_tmp(I_smooth(1:2:end,2:2:end) == 0 | isnan(P_tmp)) = 0;
   
   % if the original I_smooth pixel labels were invalid, they should have been replaced with valid labels,
   % unless they had no valid neighbors, in which case they are left alone;
   I_tmp(~I_mask) = I_smooth(~I_mask);
   
   % must convert smoothed map to integer array
   I_smooth(1:2:end,2:2:end) = int16(I_tmp);
   
   endif % if (sets(2))
   
   
   %
   % Coding Set #3:
   %
   if (sets(3))

   % create 8 neighbor pixel arrays for every other pixel (note that these neighbor_* arrays
   % are constructed with 4x the necessary pixels because I could not figure out how to 
   % efficiently construct the neighbor_* array with the coding set's dimensions...the
   % coding set is extracted subsequently to generate the counts for each couple of pair-
   % wise cliques)
   neighbor_lowerleft = 0 * I_smooth;
   neighbor_lowerleft(2:end, 2:end) = I_smooth(1:end-1, 1:end-1); 
   
   neighbor_below = 0 * I_smooth;
   neighbor_below(2:end,:) = I_smooth(1:end-1, :);
   
   neighbor_lowerright = 0 * I_smooth;
   neighbor_lowerright(2:end, 1:end-1) = I_smooth(1:end-1, 2:end); 
   
   neighbor_right = 0 * I_smooth;
   neighbor_right(:, 1:end-1) = I_smooth(:, 2:end);
   
   neighbor_upperright = 0 * I_smooth;
   neighbor_upperright(1:end-1, 1:end-1) = I_smooth(2:end, 2:end);
   
   neighbor_above = 0 * I_smooth;
   neighbor_above(1:end-1,:) = I_smooth(2:end, :);
   
   neighbor_upperleft = 0 * I_smooth;
   neighbor_upperleft(1:end-1, 2:end) = I_smooth(2:end, 1:end-1);
   
   neighbor_left = 0 * I_smooth;
   neighbor_left(:, 2:end) = I_smooth(:, 1:end-1);
         
   
   % cycle through each class and sum up those that match class j
   n1 = [];
   n2 = [];
   n3 = [];
   n4 = [];
   for j=1:nclass
   
      n1(:,:,j) = (neighbor_right(2:2:end, 1:2:end) == j) + (neighbor_left(2:2:end, 1:2:end) == j);
      n2(:,:,j) = (neighbor_below(2:2:end, 1:2:end) == j) + (neighbor_above(2:2:end, 1:2:end) == j);
      n3(:,:,j) = (neighbor_lowerleft(2:2:end, 1:2:end) == j)  + (neighbor_upperright(2:2:end, 1:2:end) == j);
      n4(:,:,j) = (neighbor_lowerright(2:2:end, 1:2:end) == j) + (neighbor_upperleft(2:2:end, 1:2:end) == j);
      
   end % for j=1:nclass
   

   % multiply n*Beta and use to calculate Coding Set #3 probabilities and normalize
   P_smooth(2:2:end, 1:2:end,:) = exp(-(repmat(-reshape(alpha, [1,1,numel(alpha)]), ...
                                             [size(I_smooth(2:2:end,1:2:end),1), size(I_smooth(2:2:end,1:2:end),2), 1]) - ...
                                       (n1 * beta(1) - (2-n1)*1 * beta(1) + ...
                                        n2 * beta(2) - (2-n2)*1 * beta(2) + ...
                                        n3 * beta(3) - (2-n3)*1 * beta(3) + ...
                                        n4 * beta(4) - (2-n4)*1 * beta(4) ) ) );
   P_norm = sum(P_smooth(2:2:end,1:2:end,:), 3);
   P_smooth(2:2:end,1:2:end,:) = P_smooth(2:2:end,1:2:end,:) ./ repmat(P_norm, [1, 1, nclass]);
   
   
   % update class IDs for Coding Set #3
   [P_tmp, I_tmp] = max(P_smooth(2:2:end,1:2:end,:), [], 3);
   
   
   % often there are multiple "max" values, which leads to a clear bias toward lower indices due to the
   % way Octave's max() function works; therefore we find all matching maxes between P_smooth and P_tmp,
   % then randomly choose one of the labels matching this probability to break the tie
   I_mask = logical(sum([n1+n2+n3+n4], 3)); % ignore pixels that match none of [1:nclass]
   I_tmp = I_tmp .* I_mask;
   max_match = P_smooth(2:2:end,1:2:end,:) == repmat(P_tmp, [1, 1, nclass]) & ...
               repmat(I_mask, [1, 1, nclass]); % create mask of matching max probabilities
   max_match_count = sum(max_match,3);
   max_match_float = double(max_match);
   max_match_float(max_match) = rand(sum(max_match(:)),1); % reduce # of calls to random number generator
   [rand_max, I_max] = max(max_match_float, [], 3);
   I_tmp(max_match_count>1) = I_max(max_match_count>1);


   % if the original I_smooth pixels were euql to zero when passed in, they should equal zero when returned;
   % if the max probability is a NaN, the multichannel pixel used to generate conditional probabilities
   % must have contained at least one NaN element, so the label should be zero...not sure this matters now
   I_tmp(I_smooth(2:2:end,1:2:end) == 0 | isnan(P_tmp)) = 0;
   
   % if the original I_smooth pixel labels were invalid, they should have been replaced with valid labels,
   % unless they had no valid neighbors, in which case they are left alone;
   I_tmp(~I_mask) = I_smooth(~I_mask);

   % must convert smoothed map to integer array
   I_smooth(2:2:end,1:2:end) = int16(I_tmp);

   endif % if (sets(3))

   
   %
   % Coding Set #4:
   %
   if (sets(4))

   % create 8 neighbor pixel arrays for every other pixel (note that these neighbor_* arrays
   % are constructed with 4x the necessary pixels because I could not figure out how to 
   % efficiently construct the neighbor_* array with the coding set's dimensions...the
   % coding set is extracted subsequently to generate the counts for each couple of pair-
   % wise cliques)
   neighbor_lowerleft = 0 * I_smooth;
   neighbor_lowerleft(2:end, 2:end) = I_smooth(1:end-1, 1:end-1); 
   
   neighbor_below = 0 * I_smooth;
   neighbor_below(2:end,:) = I_smooth(1:end-1, :);
   
   neighbor_lowerright = 0 * I_smooth;
   neighbor_lowerright(2:end, 1:end-1) = I_smooth(1:end-1, 2:end); 
   
   neighbor_right = 0 * I_smooth;
   neighbor_right(:, 1:end-1) = I_smooth(:, 2:end);
   
   neighbor_upperright = 0 * I_smooth;
   neighbor_upperright(1:end-1, 1:end-1) = I_smooth(2:end, 2:end);
   
   neighbor_above = 0 * I_smooth;
   neighbor_above(1:end-1,:) = I_smooth(2:end, :);
   
   neighbor_upperleft = 0 * I_smooth;
   neighbor_upperleft(1:end-1, 2:end) = I_smooth(2:end, 1:end-1);
   
   neighbor_left = 0 * I_smooth;
   neighbor_left(:, 2:end) = I_smooth(:, 1:end-1);
         
   
   
   % cycle through each class and sum up those that match class j
   n1 = [];
   n2 = [];
   n3 = [];
   n4 = [];
   for j=1:nclass
   
      n1(:,:,j) = (neighbor_right(2:2:end, 2:2:end) == j) + (neighbor_left(2:2:end, 2:2:end) == j);
      n2(:,:,j) = (neighbor_below(2:2:end, 2:2:end) == j) + (neighbor_above(2:2:end, 2:2:end) == j);
      n3(:,:,j) = (neighbor_lowerleft(2:2:end, 2:2:end) == j)  + (neighbor_upperright(2:2:end, 2:2:end) == j);
      n4(:,:,j) = (neighbor_lowerright(2:2:end, 2:2:end) == j) + (neighbor_upperleft(2:2:end, 2:2:end) == j);
      
   end % for j=1:nclass
   

   % multiply n*Beta and use to calculate Coding Set #4 probabilities and normalize
   P_smooth(2:2:end, 2:2:end,:) = exp(-(repmat(-reshape(alpha, [1,1,numel(alpha)]), ...
                                             [size(I_smooth(2:2:end,2:2:end),1), size(I_smooth(2:2:end,2:2:end),2), 1]) - ...
                                       (n1 * beta(1) - (2-n1)*1 * beta(1) + ...
                                        n2 * beta(2) - (2-n2)*1 * beta(2) + ...
                                        n3 * beta(3) - (2-n3)*1 * beta(3) + ...
                                        n4 * beta(4) - (2-n4)*1 * beta(4) ) ) );
   P_norm = sum(P_smooth(2:2:end,2:2:end,:), 3);
   P_smooth(2:2:end,2:2:end,:) = P_smooth(2:2:end,2:2:end,:) ./ repmat(P_norm, [1, 1, nclass]);
   
   
   % update class IDs for Coding Set #4
   [P_tmp, I_tmp] = max(P_smooth(2:2:end,2:2:end,:), [], 3);
   
   
   % often there are multiple "max" values, which leads to a clear bias toward lower indices due to the
   % way Octave's max() function works; therefore we find all matching maxes between P_smooth and P_tmp, 
   % then randomly choose one of the labels matching this probability to break the tie
   I_mask = logical(sum([n1+n2+n3+n4], 3)); % ignore pixels that match none of [1:nclass]
   I_tmp = I_tmp .* I_mask;
   max_match = P_smooth(2:2:end,2:2:end,:) == repmat(P_tmp, [1, 1, nclass]) & ...
               repmat(I_mask, [1, 1, nclass]); % create mask of matching max probabilities
   max_match_count = sum(max_match,3);
   max_match_float = double(max_match);
   max_match_float(max_match) = rand(sum(max_match(:)),1); % reduce # of calls to random number generator
   [rand_max, I_max] = max(max_match_float, [], 3);
   I_tmp(max_match_count>1) = I_max(max_match_count>1);
   

   % if the original I_smooth pixels were euql to zero when passed in, they should equal zero when returned;
   % if the max probability is a NaN, the multichannel pixel used to generate conditional probabilities
   % must have contained at least one NaN element, so the label should be zero...not sure this matters now
   I_tmp(I_smooth(2:2:end,2:2:end) == 0 | isnan(P_tmp)) = 0;
   
   % if the original I_smooth pixel labels were invalid, they should have been replaced with valid labels,
   % unless they had no valid neighbors, in which case they are left alone;
   I_tmp(~I_mask) = I_smooth(~I_mask);
   
   % must convert smoothed map to integer array
   I_smooth(2:2:end,2:2:end) = int16(I_tmp);

   endif % if (sets(4))
   
   
   %
   % Done with Coding Sets 1-4
   %
      
   %
   % perform a full synchronous update; this is faster, but does not guarantee convergence if used with
   % iterative MAP solvers like ICM or SA
   %
   if (~any(sets))

   % first, create 8 'neighbor' pixel arrays, using zeros to pad the new edges
   neighbor_lowerleft = 0 * I_smooth;
   neighbor_lowerleft(2:end, 2:end) = I_smooth(1:end-1, 1:end-1); 
   
   neighbor_below = 0 * I_smooth;
   neighbor_below(2:end,:) = I_smooth(1:end-1, :);
   
   neighbor_lowerright = 0 * I_smooth;
   neighbor_lowerright(2:end, 1:end-1) = I_smooth(1:end-1, 2:end); 
   
   neighbor_right = 0 * I_smooth;
   neighbor_right(:, 1:end-1) = I_smooth(:, 2:end);
   
   neighbor_upperright = 0 * I_smooth;
   neighbor_upperright(1:end-1, 1:end-1) = I_smooth(2:end, 2:end);
   
   neighbor_above = 0 * I_smooth;
   neighbor_above(1:end-1,:) = I_smooth(2:end, :);
   
   neighbor_upperleft = 0 * I_smooth;
   neighbor_upperleft(1:end-1, 2:end) = I_smooth(2:end, 1:end-1);
   
   neighbor_left = 0 * I_smooth;
   neighbor_left(:, 2:end) = I_smooth(:, 1:end-1);
   

   % next, cycle through classes, and sum up those that match class j
   for j=1:nclass
   
      n1(:,:,j) = (neighbor_right == j) + (neighbor_left == j);
      n2(:,:,j) = (neighbor_below == j) + (neighbor_above == j);
      n3(:,:,j) = (neighbor_lowerleft == j)  + (neighbor_upperright == j);
      n4(:,:,j) = (neighbor_lowerright == j) + (neighbor_upperleft == j);
      
   end % for j=1:nclass


# This probability calculation assumes that only matching neighbors contribute to the probability of class membership,
# and non-matching neighbors don't really contribute at all
# 
#    % finally, multiply n*Beta and use to calculate probabilities and normalize
#    P_smooth = exp(-(repmat(-reshape(alpha,[1,1,numel(alpha)]), [size(I_smooth,1), size(I_smooth,2), 1]) - n * beta));
#    P_norm = sum(P_smooth, 3);
#    P_smooth = P_smooth ./ repmat(P_norm, [1, 1, nclass]);
# 
#    % finally, multiply n*Beta and use to calculate probabilities and normalize
#    P_smooth = exp(-(repmat(-reshape(alpha,[1,1,numel(alpha)]), [size(I_smooth,1), size(I_smooth,2), 1]) - (n1 * beta(1) + ...
#                                                                                                      n2 * beta(2) + ...
#                                                                                                      n3 * beta(3) + ...
#                                                                                                      n4 * beta(4) ) ) );
#    P_norm = sum(P_smooth, 3);
#    P_smooth = P_smooth ./ repmat(P_norm, [1, 1, nclass]);


# This probability calculation assumes that both matching neighbors contribute (positively) to the probability of class
# membership, and non-matching neighbors contribute (negatively) to the probability of class membership...doesn't really change much

#    % finally, multiply n*Beta and use to calculate probabilities and normalize
#    P_smooth = exp(-(repmat(-reshape(alpha,[1,1,numel(alpha)]), [size(I_smooth,1), size(I_smooth,2), 1]) - (n*beta - (8-n)*beta)) );
#    
#    
# # This was an experiment in allowing multi-class betas. It turns out to not be nearly as easy as I first thought to
# # discourage certain classes from being neighbors. If you use a negative beta in the off-diagonal, this just reduces
# # the probability of belonging to the class corresponding to that column. If the off-diagonals are symmetric, then
# # the probability of belonging to both these classes drops. Ultimately, a third class might win out, which is really 
# # NOT what we intended. I need to think about this more, and understand all of Besag-1986's comments about "excluded
# # adjacencies", especially why it is so important to avoid synchronous label updates. -EJR 12/2011
# #    beta = [2  0  0  0  0  0  0  0
# #            0  2  0  0  0  0 -2  0
# #            0  0  2  0  0  0  0  0
# #            0  0  0  2  0  0  0  0
# #            0  0  0  0  2  0  0  0
# #            0  0  0  0  0  2  0  0
# #            0 -2  0  0  0  0  2  0
# #            0  0  0  0  0  0  0  2];
# #    n = n1 + n2 + n3 + n4;
# #    P_smooth = exp(-(repmat(-reshape(alpha,[1,1,numel(alpha)]), [size(I_smooth,1), size(I_smooth,2), 1]) - ...
# #                    sum(repmat(n, [1,1,1,nclass]) .* ...
# #                        repmat(reshape(beta,[1,1,nclass,nclass]), [size(I_smooth,1), size(I_smooth,2), 1] ), 4)   ) );
#                    
#    P_norm = sum(P_smooth, 3);
#    P_smooth = P_smooth ./ repmat(P_norm, [1, 1, nclass]);
   
   
   % finally, multiply n*Beta and use to calculate probabilities and normalize
   P_smooth = exp(-(repmat(-reshape(alpha,[1,1,numel(alpha)]), [size(I_smooth,1), size(I_smooth,2), 1]) -(n1 * beta(1) - (2-n1) * beta(1) + ...
                                                                                                          n2 * beta(2) - (2-n2) * beta(2) + ...
                                                                                                          n3 * beta(3) - (2-n3) * beta(3) + ...
                                                                                                          n4 * beta(4) - (2-n4) * beta(4) ) ) );
   
   P_norm = sum(P_smooth, 3);
   P_smooth = P_smooth ./ repmat(P_norm, [1, 1, nclass]);


   % determine new class IDs
   [P_tmp, I_tmp] = max(P_smooth, [], 3);
   
   
   % often there are multiple "max" values, which leads to a clear bias toward lower indices due to the
   % way Octave's max() function works; therefore we find all matching maxes between P_smooth and P_tmp, 
   % then randomly choose one of the labels matching this probability to break the tie
   I_mask = logical(sum([n1+n2+n3+n4], 3)); % ignore pixels that match none of [1:nclass]
   I_tmp = I_tmp .* I_mask;
   max_match = P_smooth == repmat(P_tmp, [1, 1, nclass]) & ...
               repmat(I_mask, [1, 1, nclass]); % create mask of matching max probabilities
   max_match_count = sum(max_match,3);
   max_match_float = double(max_match);
   max_match_float(max_match) = rand(sum(max_match(:)),1); % reduce # of calls to random number generator
   [rand_max, I_max] = max(max_match_float, [], 3);
   I_tmp(max_match_count>1) = I_max(max_match_count>1);
   
   
   % if the original I_smooth pixel labels were zero when passed in, they should equal zero when returned;
   % if the max probability is a NaN, the multichannel pixel used to generate conditional probabilities
   % must have contained at least one NaN element, so the label should be zero...not sure this matters now
   I_tmp(I_smooth == 0 | isnan(P_tmp)) = 0;
   
   % if the original I_smooth pixel labels were invalid, they should have been replaced with valid labels,
   % unless they had no valid neighbors, in which case they are left alone;
   I_tmp(~I_mask) = I_smooth(~I_mask);
   
   
   % convert smoothed map to integer array
   I_smooth = int16(I_tmp);
   
   
   endif % if (~any(sets))
   
   
   % this shouldn't be necessary, but do it to be safe
   P_smooth(repmat(I_smooth, [1, 1, nclass]) == 0) = 0;      


   % if a 3D I_class was passed in, return a 3D I_smooth
   if (size(I_class,3) > 1)
      
      I_smooth_tmp = I_smooth * 0;
      I_smooth_out = repmat(I_smooth_tmp, [1, 1, nclass]);
      
      for j=1:nclass
         
         % there must be a cleaner way to do this...
         I_smooth_tmp = I_smooth_tmp * 0;
         I_smooth_tmp(I_smooth == j) = j;
         I_smooth_out(:,:,j) = I_smooth_tmp;
         
      end % endfor i=1:nclass
      
      I_smooth = I_smooth_out;
   
   endif
   
      

end % function
