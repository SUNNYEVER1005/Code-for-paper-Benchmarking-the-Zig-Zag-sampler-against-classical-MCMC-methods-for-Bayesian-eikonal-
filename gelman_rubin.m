function Rhat = gelman_rubin(chain1, chain2, burn_in)
% gelmanRubin  Compute the Gelman–Rubin convergence diagnostic (R_hat)
%              for two MCMC chains.
%
%   Rhat = gelmanRubin(chain1, chain2)
%
%   Inputs:
%     chain1, chain2  —  d×N matrices, where
%                        d = number of parameters,
%                        N = number of samples per chain
%   Output:
%     Rhat            —  d×1 vector whose i-th element is R_hat for the i-th parameter
%
%   Reference formula (J = 2 chains):
%     B_i = N/(J-1) * [ (m1_i - m_i)^2 + (m2_i - m_i)^2 ]
%     W_i = ( s1_i^2 + s2_i^2 ) / J
%     Rhat_i = [ ((N-1)/N)*W_i + B_i/N ] / W_i
%
%   where:
%     m1_i, m2_i     = means of the i-th parameter in chains 1 and 2
%     m_i            = (m1_i + m2_i)/2
%     s1_i^2, s2_i^2 = within-chain variances, unbiased (divided by N-1)
%
%   Example:
%     % Generate two chains, 3 parameters, 1000 samples each
%     c1 = randn(3,1000);
%     c2 = randn(3,1000) + 0.5;   % intentionally shift mean
%     Rhat = gelmanRubin(c1, c2);
%     disp(Rhat);

    chain1 = chain1(:, burn_in:end);
    chain2 = chain2(:, burn_in:end);

    % Check input dimensions
    [d1, N1] = size(chain1);
    [d2, N2] = size(chain2);
    assert(d1 == d2 && N1 == N2, ...
        'chain1 and chain2 must have identical dimensions (d×N)');

    d = d1;
    N = N1;
    J = 2;  % number of chains

    % Compute per-parameter means for each chain (d×1)
    m1 = mean(chain1, 2);
    m2 = mean(chain2, 2);
    % Overall mean across chains
    m = (m1 + m2) / J;

    % Compute between-chain variance B (d×1)
    % For J = 2, denominator J-1 = 1
    B = N / (J-1) * ((m1 - m).^2 + (m2 - m).^2);

    % Compute within-chain variance W (d×1), unbiased (MATLAB var uses N-1)
    s1 = var(chain1, 0, 2);
    s2 = var(chain2, 0, 2);
    W = (s1 + s2) / J;

    % Calculate R-hat (d×1)
    Rhat = (((N-1)/N) * W + B/N) ./ W;
end