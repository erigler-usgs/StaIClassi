%%
%% Usage:
%%
%% [distance] = gmix_distance (Mu1, Sigma1, Prior1, Mu2, Sigma2, Prior2, metric);
%%
%% This function calculates various 'distance' measures between two Gaussian 
%% mixture model (GMM) distributions that possess distinct sets of component
%% means, covariance matrices, and priors. It does NOT implement distances that  
%% require some sort of discretization of a continuous distribution (e.g., the 
%% Earth Mover's distance), but rather ones that are known to have closed-form
%% solutions for GMMs. Available distances should all satisfy non-negativity,
%% identity, and symmetry metric axioms, but do not in general satisfy the
%% triangle inequality.
%%
%%
%% INPUTS:
%%
%% Mu1          - mean vectors for GMM...each mean is a 1Xnchan row vector,
%%                and GMM components are stored along the 3rd dimension.
%%
%% Sigma1       - covariance matrices for GMM; each covariance is
%%                a nchanXnchan matrix, and GMM components are stored along
%%                the 3rd dimension.
%%
%% Prior1       - prior probabilities of GMM component membership;
%%                each prior, or weight, is a 1X1 scalar, and GMM components
%%                are stored along the 3rd dimension.
%%
%% Mu2          - mean vectors for GMM...each mean is a 1Xnchan row vector,
%%                and GMM components are stored along the 3rd dimension.
%%
%% Sigma2       - covariance matrices for GMM; each covariance is
%%                a nchanXnchan matrix, and GMM components are stored along
%%                the 3rd dimension.
%%
%% Prior2       - prior probabilities of GMM component membership;
%%                each prior, or weight, is a 1X1 scalar, and GMM components
%%                are stored along the 3rd dimension.
%%
%% metric       - (optional) "metric" to calculate divergence between GMM
%%                components...the default value is is 1, while available
%%                options are:
%%
%% FIXME: consider limiting options to certain Renyi (i.e., alpha) divergences;
%%         note that alpha = 1/2 is equivalent to Bhattacharyya distance, alpha
%%         --> 1 is equivalent to Kullback-Leibler, and alpha = 2 is the Cauchy-
%%         Schwarz divergence; also note that, at least according to Kampa, no
%%         closed-form distance is possible for alpha<1, but alpha=1/2 would
%%         seem to be an exception...look into this more closely.
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

