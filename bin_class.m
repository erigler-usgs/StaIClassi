%%
%% Usage:
%%
%% [I_class_out] = bin_class(I_class, binsize|I_class2, weights)
%%
%% This function will bin discrete pixel class arrays such that each new 
%% binned pixel takes on the label of the majority of pixels it includes. 
%% 
%%
%% INPUTS:
%%
%% I_class may be a 2D classification/thematic map like that returned from
%% ml_class, or 3D training pixel array like that returned from train_class.m,
%% If 2D, there is no reliable way to determine the number of possible classes,
%% so the otherwise optional input weights must be passed, with each element 
%% corresponding to a different class. If 3D, each slice along dimension 3
%% corresponds to a different class, and the number of possible classes is
%% equal to the size of the 3rd dimension.  
%%
%% The 2nd input argument normally describes the number of pixels on the edge of
%% a square region that must be classified. For now binsize must be a power of 2.
%% The pixel ID assigned to this region will correspond to the pixel class
%% that occurs most often in that region.
%%
%% The 2nd input argument might also be a 2nd I_class argument, whose first two
%% dimensions must match I_class. The third dimension must either match I_class,
%% or be equal to 1. The dimensions of I_class_out will match I_class2. This is
%% potentially useful for flattening 3D I_class arrays, or expanding 2D I_class
%% arrays to 3D...simply set I_class2 to a zero-array of the proper dimensions.
%%
%% The optional input weights is a vector of weights for each possible pixel
%% class.  This is used to break "ties" when there are equal numbers of 
%% pixel classes in a given region.  It is also used to break ties that
%% occur when I_class is a 3D matrix where non-zero pixels occur on multiple 
%% slices along the 3rd dimension.  If no weights are provided, the pixel 
%% class with the largest index will win all ties.
%%
%%
%% OUTPUTS:
%%
%% I_class_out is a pixel class array whose first two dimension sizes equal  
%% size(I_class)/binsize, and whose 3rd dimension is equal in size to the
%% size of the 3rd dimension of the input I_class, or I_class2.
%%
%% nclass is the number of classes
%%
%% AUTHOR(S):     E. Joshua Rigler (EJR)
%%
%%
%% CHANGELOG:     17 Jun, 2011 - First public release (EJR)
%%                21 Oct, 2011 - Modified to handle 2nd I_class argument 
%%                               instead of binsize as 2nd argument
%%
%%

function [I_class_out, nclass] = bin_class(I_class, binsize, weights)
   
   %
   % Process input arguments first
   %
   
   % make sure there are enough inputs
   if (nargin < 2)
      error('At least two input arguments are required.');
   endif
   
   if (nargin == 2)
      weights = ones(max(size(I_class,3), size(binsize,3)), 1); % not sure this is wise -EJR 11/2011
   endif
   
   if (nargin == 3)
    if (numel(weights) == 1)
      % this is a kludge to handle a scalar weight being used to pass the number of classes
      weights = ones(weights,1) * weights;
    endif
   endif
   
   % make sure there are not too many inputs
   if (nargin > 3)
      error('No more than three input arguments are allowed');
   endif
   
   % check to see if binsize is, instead, a 2nd I_class argument, in which case
   % binsize=1 always.
   if (size(binsize,1) == size(I_class,1) && size(binsize,2) == size(I_class, 2))
      I_class2 = binsize;
      binsize = 1;
   else
      % just set I_class2 equal to 0 array matching dimension of I_class
      I_class2 = I_class * 0;
   endif
   
   % make sure binsize is a power of 2
   if ~(binsize==round(binsize) && log2(binsize)==round(log2(binsize)))
      % probably binsize should be required to divide evenly into both
      % dimensions of the image array, but that would require prior
      % knowledge of these dimensions...a different logic would need
      % to be applied to this entire routine to implement this
      error(['bin_class.m: binsize must be a power of 2']);
   endif
   
   
   % default is to output a 2D class array
   out_3d = false;
   
   if (numel(I_class) == 0)
      I_class_out = []; % just return an empty array...not sure this is a good idea, should probably throw an error
      return;
   endif
   
   % check class, structure, and sizes of input arguments
   if (isinteger(I_class2))
      
      if (size(I_class2,3) > 1)
         
         out_3d = true; % output should be 3D integer matrix
         
         I_map = I_class;  % Note: finish processing I_map below, after 
                           % weights have been figured out
         
         I_map2 = I_class2; % for processing I_class2
         
         [ny,nx,nclass] = size(I_class);
         
         
      elseif (size(I_class2,3) == 1)
         
         % if 2D integer array is passed as input, we need another way to determine
         % the number of classes
         nclass = size(I_class,3);
         if (nclass > 1)
            if (numel(weights) ~= nclass)
                error('Number of classes in I_class does not match number of classes in weights');
            endif
         else   
            nclass = numel(weights);
         endif
         
         
         if (nclass <= 1 && numel(unique(I_class)) > 1)
            % weights either wasn't passed, or it is a scalar, yet more than 1 class 
            % ID exists in the 2D map, so there is no reliable information to 
            % determine the number of possible class IDs...throw an error.
            error('Number of possible class IDs cannot be determined');
            
         endif
         
         
         % copy integer array of class IDs to I_map
         [ny,nx] = size(I_class2);
         I_map = I_class;
         I_map2 = I_class2;

      else
         error('Bad I_class dimensions');
      endif
      
   elseif (isfloat(I_class2))
      
      if (ndims(I_class2) == 3)
         
         % generate a 2D integer array of class IDs from the index along
         % the third dimension of I_class corresponding to the highest
         % probability for that pixel
         [P_tmp, I_map] = max(I_class, [], 3);
         [P_tmp2, I_map2] = max(I_class2, [], 3);
         
         % convert new I_map to integer array
         I_map = int16(I_map);
         I_map2 = int16(I_map2);
   
         % if the maximum probability is zero, this just means that the pixel was not
         % classified before being passed in, so its label should remain zero
         I_map(P_tmp == 0) = 0;
         I_map2(P_tmp2 == 0) = 0;
         
         
         % determine the image dimensions and number of class IDs
         [ny,nx,nclass] = size(I_class2);
         
         
      else         
         error('Bad I_class dimensions');
      endif
      
   else
      error('I_class is not a valid data type');
   endif
   
   
   % make sure weights are usable
   if (~exist('weights','var') || numel(weights) == 0)
      weights = 1; % this will get fixed momentarily
   endif
   if (isscalar(weights))
   
      % weights needs to be a vector for further processing
      weights = zeros(1,nclass) + weights;
   
   elseif (isvector(weights) && numel(weights)~=nclass )
      
      % weights is either NOT a vector, or it is the wrong length
      error('bin_class.m: weights must be a vector of length %d', nclass);
         
   endif
      
   
   % if I_map has a 3rd dimension, and it is equal in length to nclass, it 
   % must be a stacked array of class masks at this point in processing;
   % correct for any overlapping labels, such as might occur if a single
   % pixel were assigned more than one label in training data
   if (size(I_map,3) > 1 && size(I_map,3) == nclass)

      I_map_wgt = I_map * 0;
      for i=1:nclass
         I_map_wgt(I_map == i) = weights(i);
      end % endfor i=1:nclass
            
      [P_tmp, I_map_wgt_idx] = max(I_map_wgt, [], 3);
      I_map_wgt_idx = int16(I_map_wgt_idx);
      I_map_wgt_idx(P_tmp == 0) = 0;
      
      % now, convert 2D I_map_wgt_idx into 3D mask with no overlaps, and 
      % multiply it by I_map to prune overlaps
      for i=1:nclass
         I_map(:,:,i) = I_map(:,:,i) .* (I_map_wgt_idx == i);
      end % endfor i=1:nclass
      
      % finally convert once again to a simple 2D class array
      % (this is simple because there is no longer any overlap by construction)
      I_map = sum(I_map,3);
      
   endif
   
   
   % if I_map2 has a 3rd dimension, and it is equal in length to nclass, it 
   % must be a stacked array of class masks at this point in processing;
   % correct for any overlapping labels, such as might occur if a single
   % pixel were assigned more than one label in training data
   if (size(I_map2,3) > 1 && size(I_map2,3) == nclass)
      
      I_map2_wgt = I_map2 * 0;
      for i=1:nclass
         I_map2_wgt(I_map2 == i) = weights(i);
      end % endfor i=1:nclass
            
      [P_tmp2, I_map2_wgt_idx] = max(I_map2_wgt, [], 3);
      I_map2_wgt_idx = int16(I_map2_wgt_idx);
      I_map2_wgt_idx(P_tmp2 == 0) = 0;
      
      % now, convert 2D I_map_wgt_idx into 3D mask with no overlaps, and 
      % multiply it by I_map to prune overlaps
      for i=1:nclass
         I_map2(:,:,i) = I_map2(:,:,i) .* (I_map2_wgt_idx == i);
      end % endfor i=1:nclass
      
      % finally convert once again to a simple 2D class array
      % (this is simple because there is no longer any overlap by construction)
      I_map2 = sum(I_map2,3);
      
   endif


   % ???   
   I_tmp = zeros(ny,nx,nclass);
   for i=1:nclass
      I_tmp(:,:,i) = (I_map == i);
      I_tmp2(:,:,i) = (I_map2 == i);
   endfor
   
   I_tmp = I_tmp .* repmat(reshape(weights, 1, 1, nclass), ny, nx);
   I_tmp2 = I_tmp2 .* repmat(reshape(weights, 1, 1, nclass), ny, nx);
   
   % bin up (and average) the pixels for each class
   % bin and average pixels to produce lower resolution images
   for i=1:nclass
      
      I_tmp_binned = I_tmp(:,:,i);
      I_tmp2_binned = I_tmp2(:,:,i);
      
      for j=1:log2(binsize)
         
         I_tmp_binned = 1/4 * (I_tmp_binned(1:2:end, 1:2:end) + ...
                               I_tmp_binned(1:2:end, 2:2:end) + ...
                               I_tmp_binned(2:2:end, 1:2:end) + ...
                               I_tmp_binned(2:2:end, 2:2:end));
         
         I_tmp2_binned = 1/4 * (I_tmp2_binned(1:2:end, 1:2:end) + ...
                                I_tmp2_binned(1:2:end, 2:2:end) + ...
                                I_tmp2_binned(2:2:end, 1:2:end) + ...
                                I_tmp2_binned(2:2:end, 2:2:end));
                  
      end % endfor j=1:log2(binsize)
      
      % it would be wise to clean up this code and fix the many ambiguous variable names
      I_binned(:,:,i) = I_tmp_binned; 
      I2_binned(:,:,i) = I_tmp2_binned; 
      
   end % enfor i=1:nclass
   
   
   % generate a 2D integer array of class IDs from the index along
   % the third dimension of I_binned corresponding to the highest
   % "probability" for that pixel
   [P_tmp, I_class] = max(I_binned, [], 3);
   [P_tmp2, I_class2] = max(I2_binned, [], 3);
   
   
   % convert new I_class to integer array
   I_class = int16(I_class);
   I_class(P_tmp == 0) = 0;
   I_class2 = int16(I_class2);
   I_class2(P_tmp2 == 0) = 0;
   
   
   % merge the two 2D integer arrays, and where they differ, assign the label with the
   % higher weight (or higher index if weights are equal...this is embarassingly ugly coding,
   % and extremely slow if there are more than a moderate number mis-matched pixel labels, 
   % excluding non-labels of course...FIXME -EJR 10/2011)
   I_class_merged = 0 * I_class;
   matchedlabel = (I_class == I_class2);
   I_class_merged(matchedlabel) = I_class(matchedlabel);
   onlylabel =  (I_class ~= 0 & I_class2 == 0);
   onlylabel2 = (I_class == 0 & I_class2 ~= 0);
   I_class_merged(onlylabel) = I_class(onlylabel);
   I_class_merged(onlylabel2) = I_class2(onlylabel2);   
   mismatchedlabel = find(~(matchedlabel | onlylabel | onlylabel2));   
   for i=1:numel(mismatchedlabel)
      wgt = wgt2 = 0;
      widx = I_class(mismatchedlabel(i));
      widx2 = I_class2(mismatchedlabel(i));
      if (widx ~= 0)
        wgt = weights(widx);
      endif
      if (widx2 ~= 0)
        wgt2 = weights(widx2);
      endif
      
      % default if weights turn out to be equal
      I_class_merged(mismatchedlabel(i)) = max(I_class(mismatchedlabel(i)), ...
                                               I_class2(mismatchedlabel(i))); 
      % do nothing else if weights are equal
      if (wgt > wgt2)
        I_class_merged(mismatchedlabel(i)) = I_class(mismatchedlabel(i));
      elseif (wgt < wgt2)
        I_class_merged(mismatchedlabel(i)) = I_class2(mismatchedlabel(i));
      endif
      
   endfor
   

   % if the input array was 3D integer (e.g., a Training matrix), convert the
   % 2D I_class into a 3D I_class
   if (out_3d)
      
      I_class_tmp = zeros(size(I_class_merged,1), size(I_class_merged,2), 'int16');
      I_class_out = repmat(I_class_tmp, [1, 1, nclass]);
      
      for i=1:nclass
         
         % there must be a cleaner way to do this...
         I_class_tmp = I_class_tmp * 0;
         I_class_tmp(I_class_merged == i) = i;
         I_class_out(:,:,i) = I_class_tmp;
         
      end % endfor i=1:nclass
            
   else
      I_class_out = I_class_merged;
   endif
   
   
endfunction
