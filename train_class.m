%%
%% Usage:
%%
%% [M_class, C_class, W_class, T_class] = train_class (I_chan, T_class[, select])
%%
%%
%% This function generates multichannel Gaussian Mixture Model (GMM) class
%% statistics that may be used in ml_class.m.  It uses a multichannel image 
%% (I_chan), and specified training pixel masks (T_class).  
%%
%%
%% INPUTS:
%%
%% I_chan is a 3D multichannel image array whose first two dimensions are 
%% y-x coordinates, and whose third dimension corresponds to different 
%% channels. 
%%
%% If T_class is a scalar it specifies the number different classes that
%% must be trained.  If T_class is a 3D array whos y-x dimensions match
%% those of I_chan, each slice along the third dimension should be a mask 
%% for the corresponding pixel class (zeros are non-training pixels,
%% other-than-zeros are training pixels).
%%
%% An optional input argument 'select' can be used to restrict training
%% to a limited subset of pixel classes. It should be a vector of integer
%% values equal in length to the number of classes. The absolute value of
%% each element indicates how many components should be used for the GMM
%% model of each class probability distribution. The sign of each element
%% dictates whether or not to prompt the user to assign more labels to
%% the training data (postive yes, negative no). If select is empty or 
%% zero, skip pixel selection; if select simply is not passed, go ahead
%% and prompt the user to assign more labels to the data.
%%
%%
%% OUTPUTS:
%%
%% M_class holds the class-dependent multichannel means as a row vector. Mean
%% vectors for different classes are sorted by index along the 3rd dimension.
%% GMM components are stored along the 4th dimension.
%%
%% C_class holds the class-dependent multichannel covariance matrices. Matrices
%% for different classes are sorted by index along the 3rd dimension. GMM
%% components are stored along the 4th dimension
%%
%% If requested, the number of training pixels per class W_class, are returned,
%% as well as the updated training pixel array T_class.
%%

%%
%% REFERENCES:
%%
%% de Wit, T. D. (2006), Fast segmentation of solar extreeme ultraviolet 
%%   images, Solar Physics, v. 239.1-2, p. 519-530.
%% McLachlan, G., and Peel, D. (2000), Finite Mixture Models, Wiley and Sons.
%% Rigler, E. J., Hill, S. M., Reinard, A. A., and Steenburgh, R. A. (2012), 
%%   Solar thematic maps for space weather operations, Space Weather, v. 10.
%% Tso, B., and Mather, P. (2009), Classification Methods for Remotely Sensed
%%   Data, 2nd Ed., CRC Press.
%% Turmon, M., Jones, H. P., Malanushenko, O. V., and Pap, J. M. (2010), 
%%   Statistical feature recognition for multidimensional solar imagery, Solar
%%   Physics, v. 262, p. 277-298. (see also Turmon et al. 2002, cited within)
%%
%%
%% AUTHOR(s):     E. Joshua Rigler (EJR)
%%
%%
%% CHANGELOG:     17 Jun, 2011 - First public release (EJR)
%%                11 Mar, 2012 - Modified to call get_coord_ts and get_coord_img,
%%                               depending on dimensions of input image; also to
%%                               allow ithresh input
%%                15 Mar, 2012 - Removed ability to specify colormap and color axis.
%%                               Changed how 'select' works, and other minor mods
%%                               designed to make this function compatible with
%%                               recent changes in get_coord_img.m that allow
%%                               the user to switch channels while training.
%%                18 Oct, 2012 - Modified to train Gaussian mixture models if
%%                               abs(select) is greater than 1; heavy lifting
%%                               now handled by gmix_fit.m; cleaned up code,
%%                               and removed ability to handle cell arrays of 
%%                               input images and training pixel masks, since
%%                               this was code bloat that had not been required
%%                               since merge_class.m was written.
%%
%%

function [M_class, C_class, W_class, T_class] = train_class (I_chan, T_class, select)
	
   %
   % Process input arguments
   %
   
   % ensure adequate inputs
   if (nargin < 2)
      error('At least two input arguments are required.');
   endif
   
   
   % if T_class is a numeric scalar, construct an empty matrix 
   % or cell vector that is compatible with I_chan
   if (isnumeric(T_class) && isscalar(T_class))
            
      nclass = T_class;
      T_class = zeros(size(I_chan,1), size(I_chan,2), nclass, 'int16');
               
   endif
   
   
   % define image dimension, number of channels, number of classes, and 
   % number of GMM components per class
   [ny,nx,nchan] = size(I_chan);
   nclass = size(T_class,3);
   
   
   
   
   % check T_class dimensions
   if size(T_class,1) ~= ny || ...
      size(T_class,2) ~= nx
      
      error('I_chan and T_class are not compatible.');
   
   endif
   
   


   %
   % Define training pixels
   %
   
   
   % initialize cell array for holding pixel values for each channel for each class
   class_chan = cell(1,nclass);
   
   
   if ~exist('select','var')
      % if select was not passed, assume the user wants to select
      % training pixels for single-component GMMs for all classes
      select = true(1,nclass);
   elseif numel(select) == 0
      % if select is empty, assume the user does not want to select
      % any training pixels, and wants single-component GMMs for all classes
      select = false(1,nclass);
   elseif (numel(select) ~= nclass)
      % otherwise select needs to be a logical vector of length equal to nclass
      error('Invalid *select* input argument, *select* must be a vector equal in length to the number of classes');
   endif
   
   
   % cycle over each class, selecting or skipping the class according to 'select'
   for j=1:nclass
      
      % temporary matrix simplifies indexing later in this loop
      T_class_tmp = T_class(:,:,j);
      
      % read in pixels already assigned in T_class
      [idx_tmp] = find(T_class_tmp);         
      
      if select(j) > 0

         if (size(I_chan,1)==1 || size(I_chan,2)==1)
            % a 1 row/column image is probably a time series, and much easier to analyze as a line plot
            [idx_tmp] = get_coord_ts(I_chan, idx_tmp);
         else
            % a 2D array should always be analyzed as an image
            [idx_tmp] = get_coord_img(I_chan, idx_tmp);
         endif
      
      else
                     
         printf('Skipping training pixel selection for class %d\n', j);               
         
      endif
               
      T_class_tmp(:) = 0;
      T_class_tmp(idx_tmp) = j;
      T_class(:,:,j) = T_class_tmp;
      
      % append training pixels to class-channel pixel matrix
      pix_tmp = [];
      if (numel(idx_tmp>0))
         for i=1:nchan
            pix_tmp(:,i) = I_chan(:,:,i)(idx_tmp);
         end % endfor i=1:nchan
      endif
      class_chan{j} = [class_chan{j}; pix_tmp];
      
   end % endfor j=1:nclass


   %
   % Calculate class-based mean and covariance matrices
   %
   
   % M_class should be stacked row vectors
   M_class = zeros(1, nchan, nclass, max(max(select),1));
   
   % C_class should be stacked covariance matrices
   C_class = zeros(nchan, nchan, nclass, max(max(select),1));
   
   for j=1:nclass
      
      
      if ~isempty(class_chan{j}) 
         
         if 0 <= select(j) && select(j) <= 1
         
            % even though covm from the Octave NaN toolbox is supposed to handle NaNs, 
            % experience shows that these covariance matrices are not always valid (i.e., 
            % they are not pos-definite, and have negative determinants); for this reason, 
            % we remove any row in class_chan that contains a NaN
            rows_without_nan = all(~isnan(class_chan{j}),2);
            class_chan{j} = class_chan{j}(rows_without_nan,:);
            
            M_class(1,:,j) = mean(class_chan{j});
            C_class(:,:,j) = covm(class_chan{j}, 'D1'); % 'D1' option removes mean but normalizes by n...
                                                        % technically this is a biased estimate of the
                                                        % covariance (n-1 would provide unbiased est.),
                                                        % but while this might cause problems for very
                                                        % small sample sizes, it is necessary if we
                                                        % wish to be able to use so-called standard
                                                        % mixture reduction (SMR) techniques to later
                                                        % combine means and covariances in a robust
                                                        % manner (see Reece & Roberts, IEEE, 1/2010).
         
            % count up the number of pixels used per class
            W_class(1,1,j) = numel(rows_without_nan);
            
         else
            
            % gmix_fit is slightly slower than simply using mean() and covm() for 
            % single component GMMs, but you can force it with select = -1 to 
            % confirm that it gives identical results.
            [M_class(1,:,j,1:abs(select(j))), ...
             C_class(:,:,j,1:abs(select(j))), ...
             W_class(1,1,j,1:abs(select(j)))] = gmix_fit(class_chan{j}, abs(select(j)));
 
         end % endif 0 <= select(j) && select(j) <=1
         
        
      end % if ~isempty(class_chan{j})
      
      
   end % endfor j=1:nclass


end % function

