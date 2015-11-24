%%
%% Usage:
%%
%% [I_class, P_class, P_comp] = ml_class (I_chan, M_class, C_class[, W_class, pvalue, impute])
%%
%%
%% This function uses the class-dependent multichannel Gaussian Mixture
%% Model (GMM) distribution parameters M_class (the class+component mean 
%% vectors per channel), and C_class (the class+component covariance 
%5 matrices per channel) to assign multichannel image pixels to their
%% most likely class (i.e., the class with the highest class-conditional 
%% probability).
%% 
%%
%% INPUTS:
%%
%% I_chan is a 3D array whose first two dimensions correspond to 
%% multichannel pixel (y-x) coordinates, and whose third dimension
%% corresponds to the channel.
%%
%% M_class is a 4D array with only one 'row'.  Columns hold the
%% mean value for each channel. The 3rd dimension corresponds to
%% different classes of pixel. The 4th dimension corresponds to
%% the different components of each GMM.  M_class may be obtained 
%5 using the train_class.m function.
%%
%% C_class is a 3D array that holds the multi-spectral covariance
%% matrix for each class, class-sorted along the 3rd dimension, 
%% with GMM components stored along the 4th. C_class may be obtained 
%5 using the train_class.m function.
%%
%% W_class is an optional input that should be used to hold GMM
%% component weights. If W_class is a 1x1xNclassxNcomp vector, 
%% each GMM component is weighted by the corresponding prior. If
%% W_class has XY dimensions equal to the XY dimensions of I_chan, 
%% apply these weights to each pixel separately.
%%
%% NOTE: W_class CANNOT be used to specify class weights because GMM
%%       weights are always rescaled to sum to unity (also, the name
%%       of this function is 'ml_class' - max likelihood classify -
%%       which usually implies data-conditioned probabilties, not 
%%       prior-weighted posterior probabilities). If class weights
%%       must be applied, consider them to be Bayesian priors, and
%%       multiply them by the cumulative class probability returned
%%       from this function. A better approach might be to consider
%%       class weights as a 'clique-1' contextual prior (e.g., the
%%       'alpha' parameter to the Markov random field (MRF) model
%%       implemented in smooth_class.m).
%%
%%       An exception to this is when the sum of weights along the 
%%       4th dimension of W_class is zero, in which case this class
%%       will receive a zero weight. This might be useful if the  
%%       user wished to mask certain classes from consideration as
%%       a function of pixel location using the nyXnxXnclassXncomp
%%       W_class array.
%%
%%       If all components for all classes for a given pixel have 
%%       zero weight, that pixel will be considered 'unclassifiable', 
%%       assigned a zero probability for all classes, and assigned a 
%%       zero-label, but otherwise skipped in any further processing. 
%%
%%
%% Pvalue is an optional input. It specifies the area under a Chi- 
%% squared curve for which probabilities are considered valid . It
%% is used to determine a threshold in Mahalanobis Distance (MD)
%% that multichannel pixels in I_chan can be from each component
%% mean. If the MD threshold is exceeded, the probability that
%% the observation belongs to that GMM component is simply set to
%% zero. Pvalue can be a scalar, and  applied to all classes and
%% their components equally. It can also be a vector equal in 
%% length to the number of classes, in which case a separate
%% threshold is used for each class, but applied equally to the
%% class' components. Default pvalue is 1, which implies that all
%% probabilities are considered significant.
%%
%% FIXME: consider modifying to allow a matrix of pvalues, so that 
%%        pvalues can be specified for each component; this should
%%        be simple, since gmix_pdf already accepts such an input.
%%        It might also be nice to make pvalue spatially varying,
%%        but this would require major changes to how gmix_pdf
%%        is designed.
%%
%% Impute is an optional input. If it evaluates TRUE, replace 
%% missing channels with the conditionally expected value  
%% determined from the available channels and class-dependent
%% GMM parameters.
%% 
%% 
%% OUTPUTS:
%%
%% I_class is a 2D array of pixel classes, represented by integer
%% values from 1-Nclass. 
%%
%% The optional 3D output array P_class holds class-conditional
%% probabilities for each pixel (class corresonds to index along
%% 3rd dimension; index of max probability is I_class). P_class
%% is the weighted sum of GMM component probabilities P_comp.
%%
%% The optional 3D output array P_comp holds the un-weighted GMM
%% component probabilities. The element-wise product of P_comp
%% with the nyXnxXnclassXncomp W_class is summed along the 4th
%% dimension to give P_class.
%%
%% NOTE: If probabilities across all classes for a given pixel 
%%       equal zero, this pixel is considered unclassified, and 
%%       its assigned label will be 0.
%%
%%       If probabilities across all classes for a given pixel 
%%       are NaNs, it indicates that one or more channels in a 
%%       multichannel pixel were NaN, and therefore unclassifiable 
%%       without imputation.
%% 
%%       Imputation replaces missing data with the expected value
%%       given the assumed pixel class and available observations.
%%       If no observations are available, the pixel is given a
%%       zero label, and considered unclassifiable.
%%
%%

%%
%% REFERRENCES:
%%
%% de Wit, T. D. (2006), Fast segmentation of solar extreeme ultraviolet 
%%   images, Solar Physics, v. 239.1-2, p. 519-530.
%% Rigler, E. J., Hill, S. M., Reinard, A. A., and Steenburgh, R. A. (2012),
%%   Solar thematic maps for space weather operations, Space Weather, v. 10, 
%%   online.
%% Tso, B., and Mather, P. (2009), Classification Methods for Remotely Sensed
%%   Data, 2nd Ed., CRC Press.
%% Turmon, M., Jones, H. P., Malanushenko, O. V., and Pap, J. M. (2010), 
%%   Statistical feature recognition for multidimensional solar imagery, Solar
%%   Physics, v. 262, p. 277-298. (see also Turmon et al. 2002, cited within)%% 
%%
%% AUTHOR(S):     E. Joshua Rigler (EJR)
%%
%%
%% CHANGELOG:    17 Jun, 2011 - First public release (EJR)
%%               22 May, 2012 - Modified to impute missing values using
%%                              the class-dependent covariance matrices;
%%                              Modified to use pvalue to determine and
%%                              apply significant probability thresholds.
%%               18 Oct, 2012 - Modified to handle multi-class Gaussian
%%                              mixtures; heavy computational lifting 
%%                              now handled by gmix_cluster.m; cleaned
%%                              up code, and removed ability to handle
%%                              cell arrays of input images, since this
%%                              was code bloat that had not been at all
%%                              necessary since merge_class.m was written.
%%               02 Nov, 2012 - very slight modification to call gmix_pdf
%%                              directly instead of gmix_cluster, after
%%                              splitting the clustering and probability
%%                              calculation into two functions.
%%               05 Nov, 2012 - slight modification that makes P_class
%%                              a matrix of class probabilities, while
%%                              adding a 3rd output, P_comp, that holds
%%                              unweigted component probabilities.
%%
%%

function [I_class, P_class, P_comp] = ml_class (I_chan, M_class, C_class, W_class, pvalue, impute)
   
   %
   % Process input arguments
   %
   
   % ensure adequate number of inputs
   if (nargin < 3)
      error('At least three input arguments are required.');
   endif
   
   
   % check for M_class and C_class compatibility
   if size(M_class,2) ~= size(C_class,2) || ...
      size(M_class,3) ~= size(C_class,3) || ...
      size(M_class,4) ~= size(C_class,4)
      
      error('M_class and C_class are not compatible.');
      
   endif
   
   
   % check for I_chan *_class compatibility
   if size(I_chan,3) ~= size(M_class,2)
   
      error('I_chan and M_class/C_class are not compatible.');
      
   endif
   
   
   % define image dimension, number of channels, number of classes, and 
   % number of GMM components per class
   [ny,nx,nchan] = size(I_chan);
   nclass = size(M_class,3);
   ncomp = size(M_class,4);
   
   
   % check W_class dimensions and fix if/when it is sensible 
   
   % W_class empty or not passed in
   if (~(exist('W_class')) || isempty(W_class))
      
      % give equal prior probability to each component of each class
      W_class = ones(size(I_chan,1), size(I_chan,2), size(M_class,3), size(M_class,4));
      
   endif
   
   % W_class is a single set of class GMM weights to be applied to all pixels
   if size(W_class,1) == 1 && ...
      size(W_class,2) == 1 && ...
      size(W_class,3) == nclass && ...
      size(W_class,4) == ncomp
      
      W_class = repmat(W_class, [ny,nx,1,1]);
      
   endif
   
   % W_class is a nyXnxXnclass array of class weights
   if size(W_class,1) == ny && ...
      size(W_class,2) == nx && ...
      size(W_class,3) == nclass && ...
      size(W_class,4) ~= ncomp
      
      % it is not philosphically self-consistent to allow class weights to be
      % applied inside this function; so, just do nothing and allow the catch-
      % all block to exit with error...in the meantime, one may simply combine
      % the probabilities obtained from this function with their own ad-hoc 
      % priors, or apply exponentially-scaled, formal MRF class priors to the 
      % thematic map generated by this function using the smooth_class()  
      % 'alpha' argument.
      
   endif
   
   % if W_class cannot be made compatible with all other inputs, exit with error
   if size(W_class,1) ~= ny || ...
      size(W_class,2) ~= nx || ...
      size(W_class,3) ~= nclass || ...
      size(W_class,4) ~= ncomp
   
      error('I_chan and W_class are not compatible.');
      
   endif
   
   
   
   % normalize W_class so that all(sum(W_class,4)==1)
   W_class = W_class ./ repmat(sum(W_class, 4), [1,1,1,ncomp]);
   % NaNs most likely arose from divide by zero, so convert these all to zero
   W_class(isnan(W_class)) = 0;

      
   
   % set default pvalue if it is not passed in
   % FIXME: what does pvalue mean in context of GMMs?
   % FIXME: consider making pvalue spatially dependent on image XY
   if ~(exist('pvalue')) || numel(pvalue) == 0
      
      pvalue = ones(nclass,1);
      
   elseif isscalar(pvalue)
      
      pvalue = ones(nclass,1) * pvalue;
      
   elseif ~(isvector(pvalue) && numel(pvalue) == nclass)
      
      error('Invalid pvalue specified');
      
   endif

   
   % set default impute if it is not passed in
   if ~(exist('impute'))
      impute = false;
   endif
   
   
   
   %
   % FINALLY, start calculating probabilities and classifying pixels
   %
   
   % copy image array into linear pixel vector
   pix_tmp = reshape(I_chan, ny*nx, nchan);
   
   
   % initialize image class-conditional component probabilities array
   % as columns of linear pixel vectors
   P_comp_tmp = zeros(ny*nx, nclass, ncomp);
   
   
   % generate linear index of multichannel pixels to be classified, designated in W_class
   idx_pix = [];
   for j=1:nclass
   
      idx_pix = unique([idx_pix, find(sum(W_class(:,:,j,:), 4))']);
      
   end % for j=1:nclass
      
   
   for j=1:nclass
      
      % Use gmix_pdf to generate unweighted component probabilities of each class's GMM, 
      % but do NOT use gmix_pdf to apply weights, or determine cumulative probabilities
      % for each pixel, since we apply our own spatially dependent weights below.
      [~, P_comp_tmp(idx_pix,j,:)] = gmix_pdf(pix_tmp(idx_pix,:), ...
                                              reshape(M_class(1,:,j,:), [1,nchan,ncomp]), ...
                                              reshape(C_class(:,:,j,:), [nchan,nchan,ncomp]), ...
                                              [], pvalue(j), impute);
      
   end % endfor j=1:nclass
   
   
   % reshape and copy P_comp_tmp to P_comp
   P_comp = reshape(P_comp_tmp, ny, nx, nclass, ncomp);
   
   
   % apply prior probabilities
   P_class = sum(P_comp .* W_class, 4);
   
               
   % The assigned class will be the index along the third dimension in
   % P_class that corresponds to the maximum probability density mixture
   [P_tmp, I_class] = max(P_class, [], 3);
   
   
   % if the maximum probability is zero, all probabilities must have been
   % zero, which causes max() to return an index of 1...not desirable; 
   % if maximum probabilty is a NaN, then all probabilities must have
   % been a NaN, which indicates that at least one value in the feature
   % vector was a NaN, and therefore unclassifiable unless imputation
   % was used...if imputation *was* used, then there must not hvae been
   % any valid values in the feature vector.
   I_class(P_tmp == 0 | isnan(P_tmp)) = 0;
   
   
   % max() does not return integer-valued indices
   I_class = int16(I_class);
   
   
   % finish timing this loop...
   %printf('Finished classifying image %d in %f seconds\n', i, toc);
   %fflush(stdout);
         
end % function


