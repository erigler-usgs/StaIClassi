%%
%% Usage:
%%
%% [conf_mtrx, kappa, ratio_mtrx] = confuse(I_class[, T_class[, conf_mtrx[, weights]]])
%% 
%% This function generates a confusion/error matrix describing the accuracy of
%% a given pixel classification when compared to "truth".
%%
%%
%% INPUTS:
%%
%% The first input, I_class, is a matrix of auto-classified pixels.  If it is
%% a 2D matrix of integers, then either a 3D T_class input must be provided, where
%% the size of its 3rd dimension corresponds to the number of possible classes, or
%% a 2D T_class may be provided, along with one or both of the optional arguments 
%% (i.e., conf_mtrx or weights), which can be used to ascertain the total number of 
%% possible classes. This is generally the simplest way to call this routine.
%%
%% If I_class is a 3D matrix of integers, it is assumed that the index for each 
%% slice in the 3rd dimension is the "true" class, while the pixel value is the
%% derived class.  This format is awkward because it allows ambiguous "truth" 
%% (i.e., the same pixel may be two different types).  In this case, the
%% function bin_class.m is called upon to decide what the "truth" really is.
%% In this scenario, T_class should not be provided, or it should be an empty
%% matrix (i.e., T_class=[] if conf_mtrx and/or weights will be provided).
%%
%% The optional input conf_mtrx is a pre-existing confusion matrix to which
%% additional pixel counts will be added.  One important reason for this is so
%% that this routine can be called recursively when cell arrays of I_class
%% and T_class matrices are passed in.  It is initialized as a square matrix
%% of zeros with dimension {nclass+1, nclass+1} (the "+1" is for unclassified
%% pixels). It may be an empty matrix if/when a weights vector is required, but
%% no initial confusion matrix is readily available.
%%
%% The optional input weights is a vector of class weights.  These are required
%% when there might be ambiguous truth or classified pixels (i.e., when I_class
%% or T_class are 3D integer arrays). The weights array may also be used to
%% specify the number of possible classes if/when I_class and T_class are 2D
%% thematic maps.
%%
%%
%% OUTPUTS:
%%
%% The output conf_mtrx is a standard confusion/error matrix as described in
%% Tso & Mather (2009), Classification Methods for Remotely Sensed Data.  In
%% short:
%%
%% 1) diagonal elements of conf_mtrx hold the number of correctly classified
%%    pixels for each class;
%% 2) each element of a column that is not on the diagonal lists the number of
%%    pixels that should have been labeled the same as the diagonal, but were
%%    mislabeled according to the row they are in...this might be considered a
%%    count of false-negatives, or omission error, and is therefore a type II error;
%% 3) each element of a row that is not on the diagonal lists the number of
%%    pixels that were labeled the same as the diagonal, but should have been
%%    labeled according to the column they are in...this might be considered a
%%    count of false positives, or comission error, and is therefore a type I error;
%%
%% The overall accuracy is simply the ratio of the sum of the diagonal elements
%% to the total number of pixels.  A "producer's accuracy" is the ratio of the
%% diagonal element to the sum of its column.  A "user's accuracy" is the ratio
%% of the diagonal element to the sum of its row.
%%
%% The optional output kappa is a widely-accepted metric for specifying the
%% accuracy of the whole confusion matrix, taking into consideration chance
%% assignment of particular labels to pixels.  If all pixels are correctly
%% classified, kappa==1. If one has a confusion matrix already, and simply
%% wants to calculate kappa, pass in empty I_class and T_classs matrices.
%%
%% The optional output ratio_mtrx is the ratio of each element of a column to
%% sum of that column. The sum of each column of ratio_mtrx should be 1. 
%% 
%%
%% AUTHOR(S):     E. Joshua Rigler
%%
%% 
%% CHANGELOG:     20 Jun, 2011 - First public release (EJR)
%%                24 Oct, 2011 - Modified conf_mtrx output to include unclassiffied
%%                               pixels, and added output ratio_mtrx (EJR)
%%
%%

function [conf_mtrx, kappa, ratio_mtrx] = confuse(I_class, T_class, conf_mtrx, weights)
   
   % check inputs for consistency, and assign default initial values where appropriate
   % (this got out of control, there has to be a better way to do this)
   if (nargin < 1)
      
      error('confuse.m: at least one input argument is always required');
   
   elseif (nargin == 1)
      
      T_class = [];
      conf_mtrx = [];
      weights = [];
   
   elseif (nargin == 2)
      
      if (size(I_class,3) > 1)
        if (numel(T_class) > 0)
          error('confuse.m: non-empty T_class not allowed with 3D I_class');
        endif
        #conf_mtrx = zeros(size(I_class,3) + 1);
        #weights = size(I_class,3);
        conf_mtrx = [];
        weights = [];
      else
        if (numel(T_class) == 0 && numel(I_class) ~= 0)
          error('confuse.m: non-empty T_class array required');
        endif
        
        % T_class must have similar first 2 dimensions to I_class 
        if (size(T_class,1) ~= size(I_class,1) ||
            size(T_class,2) ~= size(I_class,2) )
          error('confuse.m: I_class and T_class must match in first 2 dimensions');
        endif
        
        conf_mtrx = [];
        weights = [];
        
      endif            
   
   elseif (nargin == 3)
      
      if (size(I_class,3) > 1)
        
        if (numel(T_class) > 0)
          error('confuse.m: non-empty T_class not allowed with 3D I_class');
        endif        
        
        if (numel(conf_mtrx) == 0)
          % do nothing here, just pass through
        elseif ~(issquare(conf_mtrx) && size(conf_mtrx,1) - 1 == size(I_class,3))
          error('confuse.m: dimensions of I_class and conf_mtrx are not consistent');
        endif
        #weights = size(I_class,3);
        weights = [];
        
      else
        
        if (numel(T_class) == 0 && numel(I_class) ~= 0)
          error('confuse.m: non-empty T_class array required');
        endif
        
        % T_class must have similar first 2 dimensions to I_class 
        if (size(T_class,1) ~= size(I_class,1) ||
            size(T_class,2) ~= size(I_class,2) )
          error('confuse.m: I_class and T_class must match in first 2 dimensions');
        endif

        if (size(T_class,3) > 1)
          if (numel(conf_mtrx) == 0)
            % do nothing here, just pass through
          elseif ~(issquare(conf_mtrx) && size(conf_mtrx,1) - 1 == size(T_class,3))
            error('confuse.m: dimensions of T_class and conf_mtrx are not consistent');
          endif
          #weights = size(T_class,3);
          weights = [];
        else
          if ~(issquare(conf_mtrx))
            error('confuse.m: conf_mtrx must be square');
          endif
          #weights = size(T_class,3);
          weights = [];
        endif
  
      endif
   
   elseif (nargin == 4)
      
      if (size(I_class,3) > 1)
        
        if (numel(T_class) > 0)
          error('confuse.m: non-empty T_class not allowed with 3D I_class');
        endif        
        
        if (numel(conf_mtrx) == 0)
          % do nothing here, just pass through
        elseif ~(issquare(conf_mtrx) && size(conf_mtrx,1) - 1 == size(I_class,3))
          error('confuse.m: dimensions of I_class and conf_mtrx are not consistent');
        endif
        
        if (numel(weights) == 0)
          % do nothing here, just pass through
        elseif (numel(conf_mtrx) == 0)
            % do nothing here, just pass through
        elseif (isvector(weights))
          if ((isscalar(weights) && weights ~= size(I_class,3)) )
            error('confuse.m: I_class and weights do not specify same number of classes');
          elseif (numel(weights) ~= size(I_class,3))
            error('confuse.m: I_class and weights vector do not specify same number of classes');
          endif
        else
          error('confuse.m: weights must be a scalar or a vector');
        endif
        
      else        
        
        if (numel(T_class) == 0 && numel(I_class) ~= 0)
          error('confuse.m: non-empty T_class array required');
        endif
        
        % T_class must have similar first 2 dimensions to I_class 
        if (size(T_class,1) ~= size(I_class,1) ||
            size(T_class,2) ~= size(I_class,2) )
          error('confuse.m: I_class and T_class must match in first 2 dimensions');
        endif
        
        if (size(T_class,3) > 1)
        
          if (numel(conf_mtrx) == 0)
            % do nothing here, just pass through
          elseif ~(issquare(conf_mtrx) && size(conf_mtrx,1) - 1 == size(T_class,3))
            error('confuse.m: dimensions of T_class and conf_mtrx are not consistent');
          endif
          
          if (numel(weights) == 0)
            % do nothing here, just pass through
          elseif (numel(conf_mtrx) == 0)
            % do nothing here, just pass through
          elseif (isvector(weights))
            if ((isscalar(weights) && weights ~= size(T_class,3)) )
              error('confuse.m: T_class and weights do not specify same number of classes');
            elseif (numel(weights) ~= size(T_class,3))
              error('confuse.m: T_class and weights vector do not specify same number of classes');
            endif
          else
            error('confuse.m: weights must be a scalar or a vector');
          endif
          
          
        else
          
          if ~(issquare(conf_mtrx))
            error('confuse.m: conf_mtrx must be square');
          endif
          
          if (numel(weights) == 0)
            % do nothing here, just pass through
          elseif (numel(conf_mtrx) == 0)
            % do nothing here, just pass through
          elseif (isvector(weights))
            if ((isscalar(weights) && weights ~= size(conf_mtrx,1) - 1) )
              error('confuse.m: conf_mtrx and weights do not specify same number of classes');
            elseif (numel(weights) ~= size(conf_mtrx,1) - 1)
              error('confuse.m: conf_mtrx and weights vector do not specify same number of classes');
            endif
          else
            error('confuse.m: weights must be a scalar or a vector');
          endif
          
        endif
   
      endif
   
   else
      
      error('confuse.m: no more than four inputs allowed');
   
   endif

   
#    
#    if (isinteger(I_class))
#       
      if (size(I_class,3) > 1)
         
         % determine nclass from I_class
         [ny, nx, nclass] = size(I_class);
         

         % determine T_class from 3D integer I_class
         T_class = zeros(size(I_class), 'int8');
         T_class(I_class > 0) = 1;
         for i=1:nclass
            T_class(:,:,i) = T_class(:,:,i) * i;
         endfor % endfor i=1:nclass
         
         
         % get rid of ambiguous truth data
         T_class = bin_class(T_class, 1, weights);
         
         % condense T_class to 2D array
         T_class = sum(T_class, 3, 'native');
         
         
         % get rid of any ambiguous classifications
         I_class = bin_class(I_class, 1, weights);
         
         % condence I_class to 2D array
         I_class = sum(I_class, 3, 'native');



      elseif (size(I_class,3) == 1)
         

         if (size(T_class,3) > 1)
            
            [ny, nx, nclass] = size(T_class);
            
            
            if (nclass < numel(unique(I_class(:))) - 1)
              % nclass needs to be at least as large as the number of unique non-zero labels in a flat
              % I_class...this is an imperfect solution to the problem of 2D scenes missing labels
              error('confuse.m: number of classes in T_class is inconsistent with I_class');
            endif
            
            
         else
            
            % if T_class is flat, then conf_mtrx and/or weights must specify the number of classes
            if (numel(conf_mtrx) == 0 && numel(weights) == 0)
              error('confuse.m: the total number of classes cannot be determined');
            endif
            
            if (numel(conf_mtrx) > 0)
                
                nclass = size(conf_mtrx,1) - 1; % don't include unclassified in number of classes
            
                if (nclass < numel(unique(I_class(:))) - 1)
                  % nclass needs to be at least as large as the number of unique non-zero labels in a flat
                  % I_class...this is an imperfect solution to the problem of 2D scenes missing labels
                  error('confuse.m: number of classes specified by conf_mtrx is inconsistent with I_class');
                endif
                
                if (numel(weights) > 0)
                  if (isvector(weights))
                    if (isscalar(weights))
                      if (nclass ~= weights)
                        error('confuse.m: conf_mtrx and weights do not specify same number of classes');
                      endif
                    else
                      if (nclass ~= numel(weights))
                        error('confuse.m: conf_mtrx and weights vector do not specify same number of classes');
                      endif
                    endif
                  else
                    error('confuse.m: weights must be a scalar or a vector');
                  endif
                else
                  weights = nclass;
                endif
            
            else
                
                if (isscalar(weights))
                  nclass = weights;
                elseif (isvector(weights))
                  nclass = numel(weights);
                else
                  error('confuse.m: weights input must be a scalar or a vector');
                endif                  

                if (nclass < numel(unique(I_class(:))) - 1)
                  % nclass needs to be at least as large as the number of unique non-zero labels in a flat
                  % I_class...this is an imperfect solution to the problem of 2D scenes missing labels
                  error('confuse.m: number of classes specified by weights is inconsistent with I_class');
                endif
            
            endif
            
         endif
                  

         % get rid of ambiguous truth data
         T_class = bin_class(T_class, 1, weights);
         
         % condense T_class to 2D array
         T_class = sum(T_class, 3, 'native');
         
      else
         
         error('Bad integer I_class dimensions');
         
      endif
      
## 
## Removed the ability to pass an array of floats representing probability of class membership, it just
## confuses things. -EJR 10/2011
##
#    elseif (isfloat(I_class) && ndims(I_class) <= 3)
#       
#       [ny,nx,nclass] = size(I_class);
#       nclass = numel(unique(I_class(I_class>0))); % not sure what I intended here -EJR 10/2011
#       
#       % I_class has 3rd dimension, and index of slice through 3rd dimension
#       % with maximum value is the classification (not sure why this is allowed);
#       % convert I_class to a 2D array of classes and process normally
#       [P_tmp, I_class] = max(I_class, [], 3);
#       I_class = int8(I_class);
#       I_class(P_tmp == 0) = 0;
#        
#       
#       
#       if (nargin == 1)
#          error('confuse.m: at least two inputs required when I_class is a float array');
#       endif
#       
#       if ~(isscalar(weights) || (isvector(weights) && numel(weights) == nclass))
#          error('confuse.m: number of weights does not match number of classes');
#       endif
#       
#       if ~(ndims(T_class) <=3 &&
#             size(T_class,1) == size(I_class,1) && ...
#             size(T_class,2) == size(I_class,2))
#          error('confuse.m: size mismatch between I_class and T_class');
#       endif
#       
#       
#       % get rid of ambiguous truth data
#       T_class = bin_class(T_class, 1, weights);
#       
#       % condense T_class to 2D array
#       T_class = sum(T_class, 3, 'native');
#       
#       
##
## Removed the ability to pass a cell array of I_classes for now. May renable this later, but hopefully
## in a cleaner way. -EJR 10/2011
##
#    elseif (iscell(I_class))
#       
# #       % come back and fix this
# #       error('confuse.m: cannot handle I_class and T_class cell arrays yet');
#       
#       % do some simple checks, then call confuse.m recursively to fill the
#       % cell array
#       
#       if ~isvector(I_class)
#          error('confuse.m: multiple I_classes must be passed using a cell vector');
#       endif
#       
#       for i=1:numel(I_class)
#          
#          if (nargin == 1)
#             
#             if (isinteger(I_class{1}) && size(I_class{1},3) > 1)
#                I_class_tmp = I_class{i};
#                T_class_tmp = []; % 
#             else
#                error('confuse.m: must have more than 1 argument unless I_class is a 3D integer array');
#             endif
#          
#          elseif (nargin >= 2)
#             % this could probably be just an "else" block
#             
#             if (isinteger(I_class{1}) && size(I_class{1},3) > 1)
#                I_class_tmp = I_class{i};
#                if (numel(T_class{i}) > 0)
#                   error('confuse.m:  non zero-sized T_class not allowed with 3D integer I_class');
#                endif
#                T_class_tmp = []; %
#             else
#                I_class_tmp = I_class{i};
#                if ~(iscell(T_class) && isvector(T_class) && numel(I_class) == numel(T_class))
#                   error('confuse.m: T_class type and/or dimensions do not match I_class');
#                endif
#                T_class_tmp = T_class{i};
#             endif
#          
#          endif
#                   
#          
#          % finally, recursively call confuse.m
#          conf_mtrx = confuse(I_class_tmp, T_class_tmp, conf_mtrx, weights);
#                   
#          
#       endfor % endfor i=1:numel(I_class)
#       
#    elseif (numel(I_class) == 0 && numel(I_class) == 0)
#       
#       % do nothing here, just pass through the empty arrays so that Kappa can be calculated
#       
#    else
#       keyboard
#       error('confuse.m: invalid I_class type');
#    
#    endif
   
   
   if ~(iscell(I_class)) % skip this block if I_class is a cell array (really only necessary if cell array processing is reenabled)
      
      %
      % with I_class and T_class reduced to non-ambiguous 2D arrays of classes, 
      % we can now begin calculating classification statistics
      %
      
      
      % check for and initialize confusion matrix
      if (numel(conf_mtrx) == 0)
         % empty confusion matrix, initialize a new one with room for unexpected zero-labels
         conf_mtrx = zeros(nclass+1,nclass+1);
      elseif (size(conf_mtrx) ~= [nclass+1, nclass+1])
         % incompatible confusion matrix
         error('confuse.m: confusion matrix input is incompatible with I_class/T_class');
      else
         % nothing needs to be done here, the matrix is already prepared
      endif
      
      
      % only compare non-zero pixel classes (is this a good idea?...nope, this causes problems
      %  when a T_class pixel is classified, but the I_class pixel is not; the only solution
      %  would seem to be to rethink this algorithm with 0 as a valid pixel label -EJR 10/2011)
      %valid_class_idx = find((T_class ~= 0) & (I_class ~= 0));
      
      % A zero value in T_class is assigned by some expert, or it is some form of ground truth,
      % but it is always assumed to be a deliberate non-classification and should therefore be
      % ignored. However, it is possible for an I_class pixel to be assigned a zero label, usually
      % because some minimum threshold probability of membership to a class is never achieved. In
      % practice, this means that a non-classified I_class pixel never gets counted as a correct
      % classification, and assumption that may need to be revisited one day. -EJR 10/2011. 
      valid_class_idx = find(T_class ~= 0);
      
      for i=1:numel(valid_class_idx)
         
         if (T_class(valid_class_idx(i)) == I_class(valid_class_idx(i)) )
            
            thisTclass = T_class(valid_class_idx(i)) + 1;
            % increment the diagonal element of conf_mtrx
            conf_mtrx(thisTclass,thisTclass) = conf_mtrx(thisTclass,thisTclass) + 1;
            
         else
            
            thisTclass = T_class(valid_class_idx(i)) + 1;
            thisIclass = I_class(valid_class_idx(i)) + 1;
            
            % increment the off-diagonal elements of conf_mtrx
            conf_mtrx(thisIclass, thisTclass) = conf_mtrx(thisIclass,thisTclass) + 1;
            
            
         endif
         
      end % endfor i=1:numel(valid_class_idx)
      
   endif % endif ~(iscell(I_class))
   
   
   if (nargout >= 2)
            
      % if requested, calculate kappa coefficient following Tso & Mather (2009)
      N = nansum(nansum(conf_mtrx));
      theta1 = nansum(diag(conf_mtrx)) / N;
      theta2 = (nansum(conf_mtrx,1) * nansum(conf_mtrx,2)) / N^2;
      kappa = (theta1-theta2) / (1-theta2);
      
   endif
   
   if (nargout == 3)
      % if requested, calculate a matrix of fractional pixel classifications
      ratio_mtrx = conf_mtrx ./ repmat(sum(conf_mtrx, 1), size(conf_mtrx,1), 1);
      
   endif
   
   
endfunction 