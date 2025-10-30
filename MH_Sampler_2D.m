function [samples, potential, acceptance_rate] = MH_Sampler_2D( ...
    thetamin, thetamax, proposal_std, init_sample, num_samples, ...
    timeobs, receiver, v0)

% MH_SAMPLER_2D  Metropolis–Hastings sampler for a 2-D target distribution
%
% Inputs:
%   thetamin, thetamax  - 1×2 vectors, lower and upper bounds for each parameter
%   proposal_std        - scalar, standard deviation of the symmetric Gaussian proposal
%   init_sample         - 2×1 vector, initial sample point
%   num_samples         - scalar, total number of samples to generate
%   timeobs, receiver, v0 - additional arguments passed to the potential-energy routine
%
% Outputs:
%   samples         - 2×num_samples matrix of sampled points
%   potential       - 1×num_samples vector of (unnormalised) log-target values
%   acceptance_rate - scalar, acceptance rate in [0,1]

% -------------------------------------------------------------------------
% Target (unnormalised) density:  π(θ) ∝ exp[ −U(θ) ]
target_dist = @(x) exp(-Potential_energy(x, timeobs, receiver, v0));

dim = size(init_sample);   % 2×1

% -------------------------------------------------------------------------
% Initialise storage
samples  = zeros(dim(1), num_samples);
potential = zeros(1, num_samples);

samples(:,1) = init_sample;
p_curr       = target_dist(init_sample);
accepted     = 0;

% -------------------------------------------------------------------------
% Main sampling loop
for i = 2:num_samples
    i
    % Propose a new candidate (symmetric Gaussian random walk)
    x_prev = samples(:, i-1);
    x_new  = x_prev + proposal_std .* randn(dim(1), 1);

    % Evaluate target density at proposed point
    p_new = target_dist(x_new);

    % Acceptance probability (symmetric proposal ⇒ q-ratio = 1)
    alpha = min(1, p_new / p_curr);

    if rand < alpha
        samples(:, i) = x_new;
        p_curr        = p_new;
        accepted      = accepted + 1;
    else
        samples(:, i) = x_prev;
    end

    % Enforce parameter bounds: reject any component outside [thetamin, thetamax]
    for t = 1:dim(1)
        if samples(t,i) < thetamin(t) || samples(t,i) > thetamax(t)
            samples(t,i) = x_prev(t);
            accepted     = accepted - 1;
        end
    end

    % Store (unnormalised) target value for the current state
    potential(i) = target_dist(samples(:, i));
end

% -------------------------------------------------------------------------
% Final acceptance rate (excluding the initial point)
acceptance_rate = accepted / (num_samples - 1);
end