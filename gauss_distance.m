%%
%% Usage:
%%
%% [dist] = gauss_distance (Mu1, Sigma1, Mu2, Sigma2, metric);
%%
%% This function calculates various 'distance' measures between two Gaussian 
%% distributions that are fully specified by their mean vectors and positive-
%% definite covariance matrices. It does NOT implement distances that might 
%% require some sort of discretization of a continuous distribution (e.g., the 
%% Earth Mover's distance). Available distances do satisfy the non-negativity,
%% identity, and symmetry metric axioms, but may or may not satisfy the triangle
%% inequality (such cases are noted below).
%%
%% FIXME: The "distances" calculated are all for single multivariate Gaussian
%%         distributions. For many, or even all of them, it is not clear that a
%%         simple wrapper can be applied to calculate distances between Gaussian
%%         mixtures. It may be best to re-do this function as a wrapper to a  
%%         gmix_dist() that deals properly with component weights.
%%
%%
%% INPUTS:
%%
%% Mu1          - 1Xnchan mean vector #1
%%
%% Sigma1       - nchanXnchan covariance matrix #1
%%
%% Mu2          - 1xnchan mean vector #2
%%
%% Sigma2       - nchanXnchan covriance matrix #2
%%
%% metric       - (optional) "metric" to calculate divergence between GMM
%%                components...the default value is is 1, while available
%%                options are:
%%
%%                1 - Bhattacharyya Distance (dB) - this is not a true metric
%%                    because it violates the so-called triangle inequality
%%                    (i.e., d(A,B) + d(B,C) >= d(A,C) is not always true).
%%                    But, it is symmetric, and has a fairly intuitive form
%%                    that separates the distribution divergence into one
%%                    component related to the distance between means, and
%%                    one related to divergence of covariances. Moreover, it
%%                    has been widely used for over half a century to measure
%%                    divergence between probability distributions.
%%                2 - Symmetric Kullback-Leibler Divergence (dKL) - even sym-
%%                    metrized, this is not a true metric, since it too fails
%%                    the triangle inequality. But, it is widely known and used,
%%                    and very well-grounded in information theory (it is often
%%                    referred to as "information divergence", and sometimes
%%                    "relative entropy"). Like dB, dKL has an intuitive form
%%                    that separates divergence into a component related to the
%%                    means, and a component related to the covariance matrices.
%%                3 - Forstner & Moonen metric (dFM, modified by Abou-Moustafa) - 
%%                    this *is* a true metric, and comprised of two terms, one
%%                    related to the distance between means, and one related to 
%%                    divergence of covariances. A downside is that it is not 
%%                    widely known or used (yet).
%%                4 - Cauchy-Schwarz PDF Distance (dCS) - it is not clear if 
%%                    this is a true metric or not. Multiple references suggest
%%                    fails to satisfy the triangle inequality, but its form is
%%                    very similar to so-called Jensen-Shannon divergence(*), 
%%                    whose square root has been proven to fully satisfy all 
%%                    metric-related axioms, including the triangle inequality.
%%                    It is well-grounded in information theory, but does not
%%                    possess the intuitive form noted for dB or dKL.
%%                5 - Hellinger Distance (dH) - this meets all axioms for a true
%%                    metric, including the triangle inequality. It depends on
%%                    the Bhattacharyya Distance by dH = sqrt(1 - exp(-dB)). It
%%                    does not provide a very intuitive closed-form like dB, and
%%                    and while.
%%                6 - Discordance (dDC, i.e., 1-Concordance Correlation) - this
%%                    is neither obviously "information theoretic", nor does it
%%                    obviously satisfy the triangle inequality requirement for
%%                    true metrics. It is symmetric though, and gives a value
%%                    between 0 and 1, not unlike the Hellinger distance.
%%
%%  (*): Ideally we would include the so-called Jensen-Shannon divergence, or
%%       its square root, for a true metric. This amounts to the difference 
%%       between the entropy of a mixture, and a mixture of the entropies. But,
%%       entropy of a continuous PDF mixture has no analytic solution, thereby 
%%       requiring discretized approximations to the PDFs, or computationally
%%       expensive iterative solutions. Several references are offered below  
%%       that describe how this might be done if it is deemed desirable.
%%
%%
%% OUTPUTS:
%%
%% dist         -  "distance" between distributions.
%%

%%
%% REFERENCES:
%%
%% Abou-Moustafa, K. T., De La Toree, F., and Ferrie, F. P. (2010), Designing
%%   a metric for the difference between Gaussian densities, Brain, Body and
%%   Machine in Advances in Intelligent and Soft Computing, Angeles et al., Eds.,
%%   v. 83, p. 57-70, Springer.
%% Endres, D. M., and Schindelin, J. E. (2003), A new metric for probability
%%   distributions, IEEE Transactions on Information Theory, v. 49.7, p 1858-1860.
%% Forstner, W., and moonen, B. (1999), A metric for covariance matrices,
%%   Technical Report, Dept. Geodesy and Geo-Informatics, Stuttgart Univ.
%% Huber, M. F., Bailey, T., Durrant-Whyte, H., and Hanebeck, U. D. (2008),
%%   On entropy approximation for Gaussian mixture random vectors, Proceedings
%%   of the IEEE International Conference on Multisensor Fusion and Integration
%%   for Intelligent Systems, Seoul, Korea.
%% Kampa, K, Hasanbelliu, E., and Principe, J.C. (2011), Closed-form Cauchy-Schwarz
%%   divergence for mixture of Gaussians, 2011 Int. Joint Conference on Neural
%%   Networks (IJCNN), Jul31-Aug5 2011, pp. 2578-2585.
%% Sfikas, G., Constantinopoulos, C., Likas, A., and Galatsanos, N.P. (2005), 
%%   An analytic distance metric for Gaussian mixture models with application 
%%   in image retrieval, Artificial Neural Netorks: Formal Models and Their 
%%   Applications - ICANN2005, Springer-Verlag, pp. 835-840
%%
%% Also, for helpful notes about Shannon-Jensen Divergence, see: 
%%   http://stats.stackexchange.com/questions/8634/jensen-shannon-divergence-for-bivariate-normal-distributions
%%
%%
%% AUTHOR(s):     E. Joshua Rigler (EJR)
%%
%%
%% CHANGELOG:     19 Dec, 2012 - First public release (EJR)
%%

function [dist] = gauss_distance(Mu1, Sigma1, Mu2, Sigma2, metric);
   
   % check minimum number of inputs
   if nargin < 4
      error('gauss_distance: at least 4 input arguments are required');
   endif
   
   
   % assign defaul metric if one is not passed
   if nargin < 5
      metric = 1;
   endif
   
   
   % check for valid Mu1 (Mu2 is checked implicitly if it is same size as Mu1)
   if ndims(Mu1) > 2 || ...
      size(Mu1,1) ~= 1
      error('gauss_distance: Mu1 must be a 1Xnchan row vector');
   else
      nchan = size(Mu1,2);
   endif
   
   
   % check for valid Sigma1 (Sigma2 is checked implicitly if it is same size as Sigma2)
   if ndims(Sigma1) > 2 || ...
      size(Sigma1,1) ~= size(Sigma1,2) || ...
      size(Sigma1,2) ~= nchan
      error('gauss_distance: Sigma1 must be a nchanXnchan square matrix');
   endif
   
   
   % check that first Gaussian is sized identically to second Gaussian
   if ~all(size(Mu1) == size(Mu2)) || ...
      ~all(size(Sigma1) == size(Sigma2))
      error('gauss_distance: multivariate Gaussians must have identical dimensions to merge');
   endif
   
   
  
   
   %
   % begin actual algorithm
   %
   
   
   % calculate distance from distribution #1 to distribution #2
   switch metric
      
      
      case {1}
         
         % Bhattacharyya Distance

% supposedly a more efficient way to scale dMu...test with bigger matrices
%         avgSigma = (Sigma1 + Sigma2) / 2;
%         dMu = (Mu1(1,:) - Mu2(1,:)) / chol(avgSigma);
%         
%         % ...as described on Wikipedia
%         dist = 0.125 * dMu * dMu' + ...
%                0.5 * log(det(avgSigma) / sqrt(det(Sigma1)*det(Sigma2) ) );
         
         
         % testing non-optimized version
         avgSigma = (Sigma1 + Sigma2) / 2;
         dMu = (Mu1(1,:) - Mu2(1,:));
         
         dist = 0.125 * dMu * inv(avgSigma) * dMu' + ...
                0.5 * log(det(avgSigma) / sqrt(det(Sigma1)*det(Sigma2) ) );
      
      
      case 2
      
         % Symmetric Kullback-Leibler Distance
         iSigma1 = inv(Sigma1);
         iSigma2 = inv(Sigma2);
         dMu = (Mu1(1,:) - Mu2(1,:)) * ...
               (iSigma1 + iSigma2) * ...
               (Mu1(1,:) - Mu2(1,:))';
         
         % probably not very speed-optimized
         dist = 0.5 * (dMu + trace(iSigma1*Sigma2 + ...
                                   iSigma2*Sigma1 - ...
                                   2*eye(nchan)) );
         
      
      case 3
         
         % modified Forstner-Moonen Metric

% supposedly a more efficient way to scale dMu...test with bigger matrices
%         avgSigma = (Sigma1 + Sigma2) / 2;
%         dMu = (Mu1(1,:) - Mu2(1,:)) / chol(avgSigma);
%         
%         % probably not very speed-optimized
%         dist = sqrt(dMu * dMu') + ...
%                sqrt(sum(log(eig(Sigma1,Sigma2) ).^2));
         
         avgSigma = (Sigma1 + Sigma2) / 2;
         dMu = (Mu1(1,:) - Mu2(1,:));
 
         dist = sqrt(dMu * inv(avgSigma) * dMu') + ...
                sqrt(sum(log(eig(Sigma1,Sigma2) ).^2));
         
      
      case 4
         
         % Cauchy-Schwarz PDF metric
         
         % This is taken from Kampa, Hasanbelliu, and Principe (2011), which 
         % refers to the "Cauchy-Schwarz" inequality...moreover, Principe seems
         % to be very knowledgable about Renyi Entropy, and dCS is nothing more
         % (or less) than Renyi's relative entropy with alpha=2.
         
         %% NOTE: the z** values are normalization constants that fall out of
         %%       a well-known Gaussian mulitplication identity; it turns out
         %%       that it is not necessary to jump through all the hoops below,
         %%       but that one may simply use the density obtained by treating
         %%       Mu2 as an observation on Mu1, with covariance Sigma1+Sigma2...
         %%       this is precisely what Kampa says...not sure why it didn't
         %%       work the first time I tried it.
         %Sigma12 = inv(inv(Sigma1) + inv(Sigma2));
         %Mu12 = Sigma12 * inv(Sigma1) * Mu1' + Sigma12*inv(Sigma2) * Mu2';
         %z12 =  sqrt(det(Sigma12) / (det(Sigma1) * det(Sigma2) * (2*pi)^nchan) ) * ...
         %       exp(-0.5 * (Mu1*inv(Sigma1)*Mu1' + Mu2*inv(Sigma2)*Mu2' - ...
         %                   Mu12'*inv(Sigma12)*Mu12) );
         %
         %Sigma11 = inv(inv(Sigma1) + inv(Sigma1));
         %Mu11 = Sigma11 * inv(Sigma1) * Mu1' + Sigma11*inv(Sigma1) * Mu1';
         %z11 =  sqrt(det(Sigma11) / (det(Sigma1) * det(Sigma1) * (2*pi)^nchan) ) * ...
         %       exp(-0.5 * (Mu1*inv(Sigma1)*Mu1' + Mu1*inv(Sigma1)*Mu1' - ...
         %                   Mu11'*inv(Sigma11)*Mu11) );
         %
         %Sigma22 = inv(inv(Sigma2) + inv(Sigma2));
         %Mu22 = Sigma22 * inv(Sigma2) * Mu2' + Sigma22*inv(Sigma2) * Mu2';
         %z22 =  sqrt(det(Sigma22) / (det(Sigma2) * det(Sigma2) * (2*pi)^nchan) ) * ...
         %       exp(-0.5 * (Mu2*inv(Sigma2)*Mu2' + Mu2*inv(Sigma2)*Mu2' - ...
         %                   Mu22'*inv(Sigma22)*Mu22) );
         
         z12 = gmix_pdf(Mu1, Mu2, Sigma1+Sigma2);
         z11 = gmix_pdf(Mu1, Mu1, Sigma1+Sigma1);
         z22 = gmix_pdf(Mu2, Mu2, Sigma2+Sigma2);
         
         dist = -log(z12) + 0.5*log(z11) + 0.5*log(z22);
         
         
         keyboard
         
      case 5
         
         % Hellinger Distance
         
         % recursively get Bhattacharyya distance first
         dB = gauss_dist(Mu1, Sigma1, Mu2, Sigma2, 1);
         
         % calculate Hellinger distance
         dist = sqrt(1 - exp(-dB));
         
         
      case 6
         
         % This is taken from Sfikas et al. (2005), which is in turn taken from
         % Ray (2003, dissertation)...it looks a little like dCS at first, but
         % is is not. Sfikas took the -logarithm of the concordance correlation
         % coefficient (CCC), but this doesn't seem as reasonable as using 1-CCC,
         % as suggested by Ray...in fact it seems like a cheap attempt to make
         % the distance appear more information theoretic. 
         % NOTE: this seems to be related (or identical?) to the "concordance 
         %       correlation coefficient" originally devised by Lin (1989), and
         %       most often used for evaluating "inter-rater reliability"...this
         %       is not unlike a continuous analog to Cohen's Kappa, nor is it
         %       functionally much different from "intraclass correlation" that
         %       was originally propsed by none other than Ronald Fisher. In
         %       short, I am not certain it really belongs within this function,
         %       but will leave it until I find a better place -EJR, 12/2012
         % 
         V12 = inv(inv(Sigma1) + inv(Sigma2));
         k12 = Mu1 * inv(Sigma1) * (Mu1 - Mu2)' + ...
               Mu2 * inv(Sigma2) * (Mu2 - Mu1)';
         
         V11 = inv(inv(Sigma1) + inv(Sigma1));
         k11 = Mu1 * inv(Sigma1) * (Mu1 - Mu1)' + ...
               Mu1 * inv(Sigma1) * (Mu1 - Mu1)';
         
         V22 = inv(inv(Sigma2) + inv(Sigma2));
         k22 = Mu2 * inv(Sigma2) * (Mu2 - Mu2)' + ...
               Mu2 * inv(Sigma2) * (Mu2 - Mu2)';
         
         
         dist = 1 - (2*sqrt(det(V12)/(exp(k12)*det(Sigma1)*det(Sigma2))) / ...
                     (sqrt(det(V11)/(exp(k11)*det(Sigma1)*det(Sigma1))) + ...
                      sqrt(det(V22)/(exp(k22)*det(Sigma2)*det(Sigma2))) ) );
         
         
         
      otherwise
         
         error('gauss_distance: unrecognized metric');
   
   endswitch
   
   
endfunction
