%%
%% Usage:
%%
%% [I_smooth, P_smooth] = metro_map_class(I_class, alpha, beta, niter, T, sets)
%% 
%% This function smooths a thematic map by increasing its joint probability in
%% a manner consistent with MRF parameters alpha and beta (see smooth_class.m). 
%% This is done, as the name suggests, using the famous Metropolis algorithm 
%% (see Metropolis et al., 1953). If a cooling schedule (i.e., monotonically
%% descreasing vector T) is passed in, this function amounts to a simulated
%% annealing (SA) algorithm which, in theory, should escape local minima in 
%% the solution and eventually converge on a global optimum.
%%
%%
%% INPUTS:
%%
%% I_class      - this is the initial map; see smooth_class.m for further 
%%                explanation.
%%
%% alpha        - exp(alpha) define class-specific weights; see smooth_class.m 
%%                for details.
%%
%% beta         - nearest neigbor coupling coefficient; see smooth_class.m 
%%                for details.
%%
%% niter        - number of iterations to perform for each T[emperature]. 
%%
%% T            - "temperature" as used in simulated annealing algorithms. If 
%%                 T is a vector, assume it is a cooling schedule, and run 
%%                 as a simulated annealing algorithm.
%%                 FIXME: To be consistent with SA as described in Tso & Mather 
%%                        (2009), the log of each probability is divided by T, 
%%                        then re-exponentiated to give the scaled probability.
%%                        This is probably needless computation, and should be 
%%                        cleaned up to improve performance.
%%
%% sets         - defines which "coding set" should be processed;
%%                see smooth_class.m for details.
%%
%%
%% OUTPUTS:
%%
%% I_smooth     - MAP smoothed thematic map
%%
%% P_smooth     - Probabilities of class membership.
%%                NOTE: unlike smooth_class.m, the maximum probability for any
%%                      given pixel will not necessarily correspond to the  
%%                      class identified in I_smooth, since there is always a 
%%                      possibility of changing a pixel label even if the label  
%%                      was already the most probable for its local neighborhood 
%%                      of pixels. That is the whole point of this algorithm, 
%%                      since this bumps the solution away from local maxima 
%%                      in probabliity space in the hopes that it will  
%%                      eventually converge on the global maximum.
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
%% CHANGELOG:     09 Dec, 2011 - (EJR) First public release.
%%                28 Dec, 2011 - (EJR) Modified to properly handle coding sets, 
%%                               and not rely on smooth_class.m for this; also 
%%                               removed  rand_state input/output, since resetting 
%%                               the random number generator can be done externally.
%%                07 Nov, 2012 - (EJR) Major updates:
%%                               * now assumes I_class is a 2D thematic map, so 
%%                                 alpha is required to specify number of classes.
%%                               * rearranged input arguments to put niter before T;
%%                                 this is how the SA algorithm is described in most
%%                                 literature, and is consistent with icm_map_class.m;
%%                               * changed name to metro_map_class.m
%%                               * Cleaned up preprocessing and docs accordingly
%%

function [I_smooth, P_smooth] = metro_map_class(I_class, alpha, beta, niter, T, sets)
  
   % just simple input checks here, more complicated checks get performed in smooth_class.m
   if (nargin < 2)
      error('metro_map_class: too few inputs');
   endif
  
   
   % at least two inputs
   if nargin > 1
      
      % define some defaults
      if ~exist('beta', 'var')
         beta = 0;
      endif
      
      if ~exist('niter','var')
         niter = 10;
      endif
      
      if ~exist('T','var')
         T = 1;
      endif
      
      if ~exist('sets','var')
         sets = [1 1 1 1]; % default to all coding sets
      endif
      
      
      % check alpha, extract nclass
      if isscalar(alpha)
         nclass = alpha;
         alpha = zeros(nclass, 1);
      elseif isvector(alpha)
         nclass = numel(alpha);
      else
         error('metro_map_class: alpha must be a scalar, or a vector equal in length to number of classes');
      endif
      
   endif
  
   
   % at least 3 inputs
   if nargin > 2
      
      % check beta
      if (isscalar(beta))
         
         beta = [beta, beta, beta, beta];
      
      elseif ~(isvector(beta) && length(beta) == 4)
         
         % only accept scalar or vector-4 values for beta
         error('Beta must be a scalar or 4 element vector');
      
      endif
      
   endif
   
   
   % at least 4 inputs
   if nargin > 3
      
      % check niter
      if ~isscalar(niter)
         error('metro_map_class: niter must be an integer scalar');
      endif
      
   endif
   
   
   % at least 5 inputs
   if nargin > 4
   
      % if T is a vector, assume it is a cooling schedule, which should always decrease monotonically
      if (isvector(T))
         if (any(diff(T) >= 0))
            error('metro_map_class: T must be monotonically decreasing');
         endif
      else
         error('metro_map_class: T must be a vector');
      endif
   
   endif
   
   
   % at least 6 inputs
   if nargin > 5
      
      % check sets input
      if isscalar(sets)
         if sets
            sets = [1 1 1 1];
         else
            sets = [0 0 0 0];
         endif
      else
         if ~(isvector(sets) && numel(sets) == 4)
            error('metro_map_class.m: input sets must be a 4-vector of logical scalars');
         endif
      endif
   
   endif


# leftover from when 3D I_class was acceptable...delete after testing
#   % use bin_class to flatten I_class and return the number of classes (bin_class.m should check for
#   %  any inconsistent inputs too)
#   [I_class, nclass] = bin_class(I_class, zeros(size(I_class,1), size(I_class,2), 'int16'), exp(alpha));
   
   
   % create a logical array corresponding to zero labels because it will be necessary to change zero
   % labels to ones when indexing the probability arrays inside the loop
   I_class_zeros = I_class == 0;
  
  
   
   % start loops
   for h=1:length(T)
      
      printf('T = %f; iter = 0', T(h));
      fflush(stdout);
      
      for i=1:niter
      
         if (~mod(i,5))
            printf('%d', i);
            fflush(stdout);
         else
            printf('.');
            fflush(stdout);
         endif
         
         
         % a small bias can be introduced due to the order coding sets are processed; 
         % we can minimize this by rotating I_class, permuting beta appropriately, and
         % finally undoing these after a complete pass through the coding sets
         nrot = mod(i,4) - 1;
         I_class = rot90(I_class, nrot);
         
         % swaps up-down/left-right and down-diag/up-diag
         nshift = mod(nrot, 2);
         beta = beta([circshift(beta(1:2), [0 nshift]), circshift(beta(3:4), [0 nshift])]); 
         
      
         % create ij subscripts for calls to sub2ind 
         [jj, ii] = meshgrid([1:columns(I_class)], [1:rows(I_class)]);
         
         
         %
         % CODING SETS: Fully synchronous updates of the map violate some conditions that guarantee
         %              convergence in iterative solvers. We need to break this up into "coding sets"
         %              according to Besag-1986, where carefully chosen subsets of the input map are 
         %              updated, then subsequent subsets are updated. No member of a subset may be a
         %              neighbor to another member of the same subset. Therefore, we will lose some 
         %              hard-won performance increase, but hopefully not too much
         %
         
         
         %
         % Coding Set #1:
         %
         if (sets(1))
            
            % use smooth_class to generate contextual probabilities (ignore smoothed map output;
            %  probs are calculated for all classes, which is more processing than required, but
            %  smooth_class is vectorized and very efficient...for Octave)
            [I_tmp, P_smooth] = smooth_class(I_class, alpha, beta, [1 0 0 0]);
            I_class(I_class_zeros) = 1;
            P_smooth_0 = reshape(P_smooth(1:2:end,1:2:end,:)(sub2ind(size(P_smooth(1:2:end,1:2:end,:)), ...
                                                                     round(ii(1:2:end,1:2:end)(:)/2), ...
                                                                     round(jj(1:2:end,1:2:end)(:)/2), ...
                                                                     I_class(1:2:end,1:2:end)(:)) ), ...
                                 size(I_class(1:2:end,1:2:end)) );
            I_class(I_class_zeros) = 0;


            % generate a corresponding array of randomly assigned labels and determine their probability from
            % P_smooth
            I_class_rand = int16(ceil(rand(size(P_smooth_0)) * nclass));
            P_smooth_1 = reshape(P_smooth(1:2:end,1:2:end,:)(sub2ind(size(P_smooth(1:2:end,1:2:end,:)), ...
                                                                     round(ii(1:2:end,1:2:end)(:)/2), ...
                                                                     round(jj(1:2:end,1:2:end)(:)/2), ...
                                                                     I_class_rand(:))), ...
                                 size(I_class_rand));

            
            % calculate probabilities for replacement 
            P_replace = min(1, exp(log(P_smooth_1)/T(h)) ./ exp(log(P_smooth_0)/T(h)));

            
            % don't replace any pixels with itentical probabilities...this doesn't just prevent
            % "replacement" with the same pixel label, but also replacement with another label
            % with the same probability, which turns out to actually be quite common, and leads
            % to an iterated bias toward lower indexed labels
            P_replace(P_smooth_1 == P_smooth_0) = 0;
            
            
            % do not replace unclassified pixels
            P_replace(I_class_zeros(1:2:end,1:2:end)) = 0;
            
            
            % replace I_class with I_class_rand with probabilty of P_replace
            % (this indexing madness is required for processing coding sets efficiently)
            comp_rand = rand(size(P_replace));
            mask_replace = zeros(size(I_class),'logical');
            mask_replace(1:2:end,1:2:end) = P_replace > comp_rand;
            I_class(mask_replace) = I_class_rand(mask_replace(1:2:end,1:2:end));
                        
         endif % if (sets(1))
         
         
         %
         % Coding Set #2:
         %
         if (sets(2))
            
            % use smooth_class to generate contextual probabilities (ignore smoothed map output;
            %  probs are calculated for all classes, which is more processing than required, but
            %  smooth_class is vectorized and very efficient...for Octave)
            [I_tmp, P_smooth] = smooth_class(I_class, alpha, beta, [0 1 0 0]);
            I_class(I_class_zeros) = 1;
            P_smooth_0 = reshape(P_smooth(1:2:end,2:2:end,:)(sub2ind(size(P_smooth(1:2:end,2:2:end,:)), ...
                                                                     round(ii(1:2:end,2:2:end)(:)/2), ...
                                                                     round(jj(1:2:end,2:2:end)(:)/2), ...
                                                                     I_class(1:2:end,2:2:end)(:))), ...
                                 size(I_class(1:2:end,2:2:end)) );
            I_class(I_class_zeros) = 0;


            % generate a corresponding array of randomly assigned labels and determine their probability from
            % P_smooth
            I_class_rand = int16(ceil(rand(size(P_smooth_0)) * nclass));
            P_smooth_1 = reshape(P_smooth(1:2:end,2:2:end,:)(sub2ind(size(P_smooth(1:2:end,2:2:end,:)), ...
                                                                     round(ii(1:2:end,2:2:end)(:)/2), ...
                                                                     round(jj(1:2:end,2:2:end)(:)/2), ...
                                                                     I_class_rand(:))), ...
                                 size(I_class_rand));

            
            % calculate probabilities for replacement 
            P_replace = min(1, exp(log(P_smooth_1)/T(h)) ./ exp(log(P_smooth_0)/T(h)));

            
            % don't replace any pixels with itentical probabilities...this doesn't just prevent
            % "replacement" with the same pixel label, but also replacement with another label
            % with the same probability, which turns out to actually be quite common, and leads
            % to an iterated bias toward lower indexed labels
            P_replace(P_smooth_1 == P_smooth_0) = 0;
            
            
            % do not replace unclassified pixels
            P_replace(I_class_zeros(1:2:end,2:2:end)) = 0;
            
            
            % replace I_class with I_class_rand with probabilty of P_replace
            % (this indexing madness is required for processing coding sets efficiently)
            comp_rand = rand(size(P_replace));
            mask_replace = zeros(size(I_class),'logical');
            mask_replace(1:2:end,2:2:end) = P_replace > comp_rand;
            I_class(mask_replace) = I_class_rand(mask_replace(1:2:end,2:2:end));      
            
         endif % if (sets(2))
         
         
         %
         % Coding Set #3:
         %
         if (sets(3))
            
            % use smooth_class to generate contextual probabilities (ignore smoothed map output;
            %  probs are calculated for all classes, which is more processing than required, but
            %  smooth_class is vectorized and very efficient...for Octave)
            [I_tmp, P_smooth] = smooth_class(I_class, alpha, beta, [0 0 1 0]);
            I_class(I_class_zeros) = 1;
            P_smooth_0 = reshape(P_smooth(2:2:end,1:2:end,:)(sub2ind(size(P_smooth(2:2:end,1:2:end,:)), ...
                                                                     round(ii(2:2:end,1:2:end)(:)/2), ...
                                                                     round(jj(2:2:end,1:2:end)(:)/2), ...
                                                                     I_class(2:2:end,1:2:end)(:))), ...
                                 size(I_class(2:2:end,1:2:end)) );
            I_class(I_class_zeros) = 0;


            % generate a corresponding array of randomly assigned labels and determine their probability from
            % P_smooth
            I_class_rand = int16(ceil(rand(size(P_smooth_0)) * nclass));
            P_smooth_1 = reshape(P_smooth(2:2:end,1:2:end,:)(sub2ind(size(P_smooth(2:2:end,1:2:end,:)), ...
                                                                     round(ii(2:2:end,1:2:end)(:)/2), ...
                                                                     round(jj(2:2:end,1:2:end)(:)/2), ...
                                                                     I_class_rand(:))), ...
                                 size(I_class_rand));

            
            % calculate probabilities for replacement 
            P_replace = min(1, exp(log(P_smooth_1)/T(h)) ./ exp(log(P_smooth_0)/T(h)));

            
            % don't replace any pixels with itentical probabilities...this doesn't just prevent
            % "replacement" with the same pixel label, but also replacement with another label
            % with the same probability, which turns out to actually be quite common, and leads
            % to an iterated bias toward lower indexed labels
            P_replace(P_smooth_1 == P_smooth_0) = 0;
            
            
            % do not replace unclassified pixels
            P_replace(I_class_zeros(2:2:end,1:2:end)) = 0;
            
            
            % replace I_class with I_class_rand with probabilty of P_replace
            % (this indexing madness is required for processing coding sets efficiently)
            comp_rand = rand(size(P_replace));
            mask_replace = zeros(size(I_class),'logical');
            mask_replace(2:2:end,1:2:end) = P_replace > comp_rand;
            I_class(mask_replace) = I_class_rand(mask_replace(2:2:end,1:2:end));            
            
         endif % if (sets(3))
         
         
         %
         % Coding Set #4:
         %
         if (sets(4))
            
            % use smooth_class to generate contextual probabilities (ignore smoothed map output;
            %  probs are calculated for all classes, which is more processing than required, but
            %  smooth_class is vectorized and very efficient...for Octave)
            [I_tmp, P_smooth] = smooth_class(I_class, alpha, beta, [0 0 0 1]);
            I_class(I_class_zeros) = 1;
            P_smooth_0 = reshape(P_smooth(2:2:end,2:2:end,:)(sub2ind(size(P_smooth(2:2:end,2:2:end,:)), ...
                                                                     round(ii(2:2:end,2:2:end)(:)/2), ...
                                                                     round(jj(2:2:end,2:2:end)(:)/2), ...
                                                                     I_class(2:2:end,2:2:end)(:))), ...
                                 size(I_class(2:2:end,2:2:end)) );
            I_class(I_class_zeros) = 0;


            % generate a corresponding array of randomly assigned labels and determine their probability from
            % P_smooth
            I_class_rand = int16(ceil(rand(size(P_smooth_0)) * nclass));
            P_smooth_1 = reshape(P_smooth(2:2:end,2:2:end,:)(sub2ind(size(P_smooth(2:2:end,2:2:end,:)), ...
                                                                     round(ii(2:2:end,2:2:end)(:)/2), ...
                                                                     round(jj(2:2:end,2:2:end)(:)/2), ...
                                                                     I_class_rand(:))), ...
                                 size(I_class_rand));

            
            % calculate probabilities for replacement 
            P_replace = min(1, exp(log(P_smooth_1)/T(h)) ./ exp(log(P_smooth_0)/T(h)));

            
            % don't replace any pixels with itentical probabilities...this doesn't just prevent
            % "replacement" with the same pixel label, but also replacement with another label
            % with the same probability, which turns out to actually be quite common, and leads
            % to an iterated bias toward lower indexed labels
            P_replace(P_smooth_1 == P_smooth_0) = 0;
            
            
            % do not replace unclassified pixels
            P_replace(I_class_zeros(2:2:end,2:2:end)) = 0;
            
            
            % replace I_class with I_class_rand with probabilty of P_replace
            % (this indexing madness is required for processing coding sets efficiently)
            comp_rand = rand(size(P_replace));
            mask_replace = zeros(size(I_class),'logical');
            mask_replace(2:2:end,2:2:end) = P_replace > comp_rand;
            I_class(mask_replace) = I_class_rand(mask_replace(2:2:end,2:2:end));                  
            
         endif % if (sets(4))
         
         
         %
         % Done with Coding Sets 1-4
         %
            
         %
         % perform a full synchronous update; this is faster, but does not guarantee convergence 
         %
         if (~any(sets))
            
            % use smooth_class to generate contextual probabilities (ignore smoothed map output;
            %  probs are calculated for all classes, which is more processing than required, but
            %  smooth_class is vectorized and very efficient...for Octave)
            [I_tmp, P_smooth] = smooth_class(I_class, alpha, beta);
            I_class(I_class_zeros) = 1;
            P_smooth_0 = reshape(P_smooth(sub2ind(size(P_smooth), ii(:), jj(:), I_class(:))), size(I_class));
            I_class(I_class_zeros) = 0;
            
            
            % generate a corresponding array of randomly assigned labels and determine their probability from
            % P_smooth
            I_class_rand = int16(ceil(rand(size(I_tmp)) * nclass));
            P_smooth_1 = reshape(P_smooth(sub2ind(size(P_smooth), ii(:), jj(:), I_class_rand(:))), size(I_class_rand));

            
            % calculate probabilities for replacement 
            P_replace = min(1, exp(log(P_smooth_1)/T(h)) ./ exp(log(P_smooth_0)/T(h)));

            
            % don't replace any pixels with itentical probabilities...this doesn't just prevent
            % "replacement" with the same pixel label, but also replacement with another label
            % with the same probability, which turns out to actually be quite common, and leads
            % to an iterated bias toward lower indexed labels
            P_replace(P_smooth_1 == P_smooth_0) = 0;
            
            
            % do not replace unclassified pixels
            P_replace(I_class_zeros) = 0;
            
            
            % replace I_class with I_class_rand with probabilty of P_replace
            comp_rand = rand(size(I_class));
            I_class(P_replace > comp_rand) = I_class_rand(P_replace > comp_rand);
            
         endif
         
         % derotate the image that was rotated so that bias arising from the order coding sets 
         % are processed is minimized; un-shift the beta coefficients
         I_class = rot90(I_class, -nrot);
         beta = beta([circshift(beta(1:2), [0 -nshift]), circshift(beta(3:4), [0 -nshift])]); 

         
      endfor % for i=1:niter
      
      printf('\n');
   
   endfor % for h=1:length(T)
   
   
   I_smooth = I_class;
  
  
endfunction