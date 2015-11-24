%%
%% Usage:
%%
%% [M_class, C_class, N_class] = merge_class (M_class1, C_class1, N_class1, ...
%%                                            M_class2, C_class2, N_class2);
%%
%%
%% This function combines/merges multichannel class statistical parameters.
%% It might  be used to combine parameters generated from two independent
%% training data sets using train_class. It is mostly a wrapper for the
%% function gmix_merge, which uses a standard mixture reduction (SMR) 
%% technique described in Reece &  Roberts (2010).
%%
%%
%% INPUTS:
%%
%% M_class1 and M_class2 hold class-dependent GMM component mean vectors 
%% to be merged.
%% 
%% C_class1 and C_class2 hold class-dependent GMM component covariance 
%% matrices to be merged.
%% 
%% N_class1 and N_class2 hold, for each class, the number of pixels used 
%% used to generate [M|C]_class1 and [M|C]_class2 respectively. 
%%
%%
%% OUTPUTS:
%%
%% M_class is the merged GMM component means vector given inputs 
%% M_class[1|2] and N_class[1|2].
%%
%% C_class is the merged GMM comonent covariance matrix given inputs 
%% C_class[1|2] and N_class[1|2].
%%
%% N_class is the merged number of pixel locations for each class, or 
%% quite simply, the sum of N_class1 and N_class2.
%%
%% Merging [M|C]_class1 with [M|C}_class2 via SMR will correspond exactly to what 
%% would be obtained if all N_class1+N_class2 pixels had been processed simultaneously 
%% to generate a unified M_class and C_class.  This result is therefore statistically 
%% "consistent" with each set of input class statistics.
%%
%% NOTE1: SMR-merged covariances only match the cumulative covariances exactly if their input
%%        covariances were determined with the sample mean removed, but normalized by N instead 
%%        of (N-1).  This is considered a biased estimate of the sample covariance, or in some
%%        literature, it is referred to as the maximum likelihood covariance estimate.  When
%%        sample sizes are large, none of this really matters, but we explain it here for
%%        completeness.
%%
%% NOTE2: While N_class* *could* be determined from the T_class training data matrix generated
%%        by the function train_class, it really should be obtained from train_class's N_class
%%        output because this corrects for data that may have been pruned prior to calculating
%%        M_class* and C_class* because there were NaNs in it.
%%
%% FIXME: N_class* may be values other than actual pixel counts per class, thereby acting as 
%%        arbitrary weight for each class, however there is no guarantee that the result will 
%%        be "consistent" with the two sets of input class statistics.  To achieve consistency 
%%        with arbitrary class weights, a Covariance Union (CU) algorithm must be implemented
%%        If this is done, the input syntax to this function should not change, except to add
%%        another optional input that designates that CU algorithm should be used, which would
%%        in turn use N_class1 and N_class2 to determine the normalized probabilities necessary
%%        to generate optimal and consistent class statistics with the iterative CU algorithm.
%%        An algorithm for CU can be found in the Reece and Roberts paper mentioned above.
%%

%%
%% REFERENCES:
%%
%% Reece, S., and Roberts, S. (2010), Generalized Covariance Union: a unified
%%   approach to hypothesis merging in tracking, IEEE Trans. Aerospace Elec.
%%   Systems, v. 46, n. 1, p. 207-221.
%%
%%
%% AUTHOR(S):    E. Joshua Rigler (EJR)
%%
%%
%% CHANGELOG:    17 Jun, 2011 - First public release (EJR)
%%               02 Nov, 2012 - Modified to now call gmix_merge, thereby giving it the
%%                              ability to handle multicomponent Gaussian mixtures
%%
%%

function [M_class, C_class, N_class] = merge_class (M_class1, C_class1, N_class1, ...
                                                    M_class2, C_class2, N_class2)
   
   % Fairly simple input argument check for now
   if (nargin ~= 6)
      error ('merge_class: wrong number of input arguments');
   end % if (narging ~= 6)
      
      
   % Check dimensionality of inputs
   % ...this mess needs cleaning -EJR 10/2012
   if (size(N_class1,1) ~= 1 || size(N_class1,2) ~= 1 || ... % not scalars
       size(N_class2,1) ~= 1 || size(N_class2,2) ~= 1 || ... % not scalars
       size(M_class1,1) ~= 1 || size(M_class2,1) ~= 1 || ... % not single column vectors
       size(C_class1,1) ~= size(C_class1,2) || size(C_class2,1) ~= size(C_class2,2) || ... % not square matrices
       size(M_class1,2) ~= size(M_class2,2) || size(M_class1,3) ~= size(M_class2,3) || ...
       size(C_class1,1) ~= size(C_class2,1) || size(C_class1,2) ~= size(C_class2,2) || size(C_class1,3) ~= size(C_class2,3) || ... % not same size inputs
       size(M_class1,2) ~= size(C_class1,2) || size(M_class2,2) ~= size(C_class2,2) || ... % incompatible M_class & C_class
       size(M_class1,3) ~= size(C_class1,3) || size(M_class1,3) ~= size(N_class1,3)) % not same size in 3rd dim.
      
      error ('merge_class: wrong or incompatible dimensions in input arguments');
      
   end % if (size(N_class,1) ~= 1 ...)
   
   
   % Pre-allocate matrices for outputs
   N_class = N_class1 * 0;
   M_class = M_class1 * 0;
   C_class = C_class1 * 0;
   
   % Loop over classes and adjust mean and covariance matrices using gmix_merge
   for k=1:size(N_class1,3)
      
      [M_class_tmp, C_class_tmp, N_class_tmp] = gmix_merge(permute(M_class1(:,:,k,:), [1 2 4 3]), ...
                                                           permute(C_class1(:,:,k,:), [1 2 4 3]), ...
                                                           permute(N_class1(:,:,k,:), [1 2 4 3]), ...
                                                           permute(M_class2(:,:,k,:), [1 2 4 3]), ...
                                                           permute(C_class2(:,:,k,:), [1 2 4 3]), ...
                                                           permute(N_class2(:,:,k,:), [1 2 4 3]));
      
      M_class(:,:,k,:) = M_class_tmp;
      C_class(:,:,k,:) = C_class_tmp;
      N_class(:,:,k,:) = N_class_tmp;      
      
   end % for k=1:numel(N_class1)
   
   
end % function merge_class
