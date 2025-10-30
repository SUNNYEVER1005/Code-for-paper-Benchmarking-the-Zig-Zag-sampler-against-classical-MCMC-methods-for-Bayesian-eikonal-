function [samples, logp_samples, epsilon] = runNUTS(f, theta0, varargin)
%  ess_mean, ess_sec_moment
% runNUTS - Sample from a user-defined posterior distribution using NUTS with dual-averaging
%
% Args:
%   f          - function handle that returns [logp, grad], e.g. @(x) deal(logp(x), grad(x))
%   theta0     - initial state (column vector)
%
% Optional Name-Value parameters:
%   'delta'        - target acceptance probability (default 0.8)
%   'n_warmup'     - number of warmup (adaptation) iterations (default 500)
%   'n_mcmc'       - number of post-warmup MCMC iterations (default 2000)
%   'max_tree_depth' - maximum tree depth for NUTS (default 10)
%   'n_updates'    - number of progress prints (default 10)
%
% Returns:
%   samples        - (d x n_mcmc) matrix of posterior samples
%   logp_samples   - (n_mcmc x 1) vector of corresponding log probabilities
%   epsilon        - adapted step size
%   ess_mean       - effective sample size for each parameter's mean
%   ess_sec_moment - effective sample size for each parameter's second moment
%
% Example:
%   f = @(x) deal(-0.5*x'*x, -x);
%   theta0 = randn(5,1);
%   [samples, logp_samples, epsilon, ess_mean] = runNUTS(f, theta0, 'n_mcmc', 1000);

% ------------------------------
% Parse input arguments
p = inputParser;
addParameter(p, 'delta', 0.8);
addParameter(p, 'n_warmup', 500);
addParameter(p, 'n_mcmc', 2000);
addParameter(p, 'max_tree_depth', 10);
addParameter(p, 'n_updates', 10);
parse(p, varargin{:});

delta = p.Results.delta;
n_warmup = p.Results.n_warmup;
n_mcmc = p.Results.n_mcmc;
max_tree_depth = p.Results.max_tree_depth;
n_updates = p.Results.n_updates;

% ------------------------------
% Step-size adaptation (dual-averaging)
[theta, epsilon, epsilon_seq, epsilonbar_seq] = dualAveraging(f, theta0, delta, n_warmup);

figure;
set(0,'defaultAxesFontSize', 14);
plot(epsilon_seq, 'b-'); hold on;
plot(epsilonbar_seq, 'r-');
legend('epsilon sequence','epsilonbar sequence');
title('Dual-Averaging Step-size Adaptation');
xlabel('Iteration'); ylabel('Step-size'); hold off;

% ------------------------------
% NUTS main sampling
samples = zeros(length(theta), n_mcmc);
logp_samples = zeros(n_mcmc, 1);
[samples(:,1), ~, nfevals_total, logp_samples(1)] = NUTS(f, epsilon, theta);

n_itr_per_update = ceil(n_mcmc / n_updates);
for i = 2:n_mcmc
    [samples(:,i), ~, nfevals, logp_samples(i)] = NUTS(f, epsilon, samples(:,i-1));
    nfevals_total = nfevals_total + nfevals;
    if mod(i, n_itr_per_update) == 0
        fprintf('%d iterations completed.\n', i);
    end
end
fprintf('Average gradient evaluations per iteration: %.1f\n', nfevals_total / n_mcmc);

% ------------------------------
% Traceplot
figure;
plot(logp_samples);
xlabel('Iteration'); ylabel('Log probability');
title('Traceplot of log(\pi(\theta))');

end
