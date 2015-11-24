%%
%% Usage:
%%
%% [I_chan] = gmix_sample (Mu, Sigma, Priors)
%%
%% This function randomly samples a multivariate Gaussian Mixture model (GMM).
%%
%%
%% INPUTS:
%%
%% Mu           - means for GMM, components stored along  3rd dimension.
%%
%% Sigma        - covariances for GMM; components stored along 3rd dimension.
%%
%% Priors       - priors for GMM; components stored along 3rd dimension.
%%
%% N            - number of samples to draw
%% 
%%
%% OUTPUTS
%%
%% I_chan       - multichannel randomly sampled from GMM distribution defined
%%                by Mu, Sigma, and Priors.
%%

%% REFERENCES:
%%
%% ???
%% 
%% 
%% AUTHOR(s):     E. Joshua Rigler (EJR)
%%
%%
%% CHANGELOG:     31 Jan, 2013 - First public release (EJR)
%%
%%

function [I_chan] = gmix_sample (Mu, Sigma, Priors, N)
   
   
   %
   % Check inputs
   %
   
   
   
   % get number of channels
   nchan = size(Mu,2);
   
   % get number of components
   ncomp = size(Mu,3);
   
   % scale Priors to [0 1]
   probs = [0; squeeze(Priors ./ sum(Priors,3))];
   
   
   % initialize output array
   I_chan = zeros(N, nchan);
   
   % create a uniformly-distributed random vector
   U = rand(N,1);
   
   % loop over components
   for k=1:ncomp
      
      lowerprob = sum(probs(1:k));
      upperprob = sum(probs(1:k+1));
      kidx = lowerprob < U & U <= upperprob;
      
      I_chan = I_chan + repmat(kidx, [1, nchan]) .* mvnrnd(Mu(:,:,k), Sigma(:,:,k), N);
      
   endfor
   
%   keyboard
   
endfunction
