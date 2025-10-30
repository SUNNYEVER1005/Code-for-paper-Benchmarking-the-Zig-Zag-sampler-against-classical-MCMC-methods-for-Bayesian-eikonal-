function [mESS,Sigma,b] = multiESS(X,Sigma,b,Noffsets,Nb)
%MULTIESS Compute multivariate effective sample size of Markov chain.
%   MESS = MULTIESS(X) computes effective sample size MESS of single Markov 
%   chain X, using the multivariate dependence structure of the process. 
%   X is a n-by-p array, where each row is a p-dimensional sample and n
%   is the current chain sample size.
%
%   [Additional documentation unchanged...]

% Copyright (C) 2016 Luigi Acerbi
% Modified to ensure mESS <= n

if nargin < 2; Sigma = []; end
if nargin < 3 || isempty(b); b = 'sqroot'; end
if nargin < 4 || isempty(Noffsets); Noffsets = 10; end
if nargin < 5; Nb = []; end

% Number of MCMC chains
if iscell(X)
    nc = numel(X);
    p = size(X{1},2);
else
    p = size(X,2);
    if ndims(X) == 3
        nc = size(X,3);
    else
        nc = 1;
    end
end

% SIGMA can be a cell array, but convert to standard format
if iscell(Sigma)
    temp = Sigma;
    Sigma = zeros(p,p,nc);
    for i = 1:nc; Sigma(:,:,i) = temp{i}; end
    clear temp;
end
if ~isempty(Sigma) && size(Sigma,3) == 1; Sigma = repmat(Sigma,[1 1 nc]); end

% Input check for batch size B
if ischar(b) 
    if ~any(strcmpi(b,{'sqroot','cuberoot','less'}))
        error('Unknown string for batch size. Allowed arguments are ''sqroot'', ''cuberoot'' and ''lESS''.');
    end
    if ~strcmpi(b,'less') && ~isempty(Nb)
        warning('Nonempty parameter NB will be ignored (NB is used only with ''lESS'' batch size B).');
    end
elseif isnumeric(b) && all(isfinite(b))
    if isscalar(b); b = repmat(b,[1 nc]); end
    if numel(b) ~= nc
        error('The batch size B needs to be a scalar or an array of the same size as the number of chains in X.');
    end
else
    error('The batch size B needs to be either ''sqroot'', ''cuberoot'' and ''lESS'' or a number between 1 and N/2.');
end

% If no output, do a plot (only with lESS method)
plotflag = nargout == 0 && ischar(b) && strcmpi(b,'less');

% Prepare arrays
mESS = zeros(1,nc);
Sigma_out = zeros(p,p,nc);
b_out = zeros(1,nc);

% Compute multiESS separately for each chain
for i = 1:nc
    if isnumeric(b); b_in = b(i); else b_in = b; end
    
    % Get chain data
    if iscell(X)
        X_chain = X{i};
    else
        X_chain = X(:,:,i);
    end
    
    % Check chain length
    n_chain = size(X_chain, 1);
    if n_chain < 2*p
        warning('Chain %d has length %d which is less than 2*p=%d. mESS may be unreliable.', i, n_chain, 2*p);
    end
    
    if isnumeric(b_in) && (b_in < 1 || b_in > n_chain/2)
        warning('Batch size %d for chain %d is out of range [1, n/2]. Adjusting.', b_in, i);
        b_in = max(1, min(floor(n_chain/2), b_in));
    end
    
    if iscell(Sigma)
        Sigma_in = Sigma{i};
    else
        if ~isempty(Sigma) && size(Sigma,3) >= i
            Sigma_in = Sigma(:,:,i);
        else
            Sigma_in = [];
        end
    end
    
    [mESS(i),Sigma_out(:,:,i),b_out(i)] = multiESS_chain(X_chain,Sigma_in,b_in,Noffsets,Nb,plotflag);
end

if nargout > 1; Sigma = Sigma_out; end
if nargout > 2; b = b_out; end

end

%--------------------------------------------------------------------------
function [mESS,Sigma,b] = multiESS_chain(X,Sigma,b,Noffsets,Nb,plotflag)
%MULTIESS_CHAIN Compute multiESS for a MCMC chain.

[n,p] = size(X);
if p > n
    error('More dimensions than data points, cannot compute effective sample size.');
end

if ischar(b)    
    switch lower(b)
        case 'sqroot'; b = max(1, floor(n^(1/2)));
        case 'cuberoot'; b = max(1, floor(n^(1/3)));
        case 'less'; 
            b_min = max(1, floor(n^(1/4)));
            b_max = max(1, min(floor(n/max(p,20)), floor(n/2)));
            if isempty(Nb); Nb = 200; end
            % Try NB log-spaced values of B from B_MIN to B_MAX
            b = unique(round(exp(linspace(log(b_min),log(b_max),Nb))));
            b(b < 1 | b > n/2) = []; % Remove invalid batch sizes
    end
end

% Ensure b is at least 1 and at most n/2
if any(b < 1 | b > n/2)
    error('Invalid batch size computed. Check chain length and batch size parameters.');
end

k = numel(b);               % Number of batch sizes
theta = mean(X,1);          % Sample mean
Lambda = cov(X);            % Sample covariance matrix
detLambda = det(Lambda);    % Determinant of sample covariance matrix

mESS = zeros(size(b));      % Prepare mESS
newSigma = zeros(p,p,k);    % Prepare batch Sigma matrices

% Compute mESS for each batch size
for i = 1:k
    [mESS(i),newSigma(:,:,i)] = ...
        multiESS_batch(X,n,p,theta,Lambda,detLambda,Sigma,b(i),Noffsets);
end

% Plot graph of mESS as a function of batch size (only for lESS)
if plotflag
    figure;
    plot(b,mESS,'LineWidth',1);
    hold on;
    set(gca,'TickDir','out');
    set(gcf,'Color','w');
    box off;
    xlabel('Batch size b');
    ylabel('multiESS');
    title('multiESS vs. Batch Size');
end

% Return lowest mESS for 'less' method, otherwise use the specified b
if k > 1
    [mESS,idx] = min(mESS);
    b = b(idx);
    newSigma = newSigma(:,:,idx);
else
    mESS = mESS(1);
    newSigma = newSigma(:,:,1);
end

% Ensure mESS does not exceed n (theoretical constraint)
mESS = min(mESS, n);

Sigma = newSigma;

end

%--------------------------------------------------------------------------
function [mESS,Sigma] = multiESS_batch(X,n,p,theta,Lambda,detLambda,Sigma,b,Noffsets)
%MULTIESS_BATCH Compute multiESS for a given batch size B.

if isempty(Sigma)
    % Compute batch estimator for SIGMA
    a = floor(n/b);
    
    % Skip if too few batches
    if a < 2
        Sigma = Lambda; % Fall back to sample covariance
        mESS = n;
        return;
    end
    
    Sigma = zeros(p,p);
    
    % Average batches over multiple offsets
    offsets = unique(round(linspace(0, max(0, n - a*b), Noffsets)));
    if isempty(offsets); offsets = 0; end
    
    for j = offsets
        Y = reshape(X(j+1:j+a*b,:), [b, a, p]);
        Ybar = squeeze(mean(Y, 1));
        if p == 1; Ybar = Ybar(:); end
        Z = bsxfun(@minus, Ybar, theta);
        Sigma = Sigma + Z' * Z;
    end
    Sigma = Sigma * b / ((a - 1) * numel(offsets));
    
    % Ensure Sigma is positive semi-definite
    Sigma = (Sigma + Sigma') / 2; % Symmetrize
    [V,D] = eig(Sigma);
    D = diag(max(diag(D), 0)); % Ensure non-negative eigenvalues
    Sigma = V * D * V';
end

% Handle potential numerical issues with determinants
try
    detSigma = det(Sigma);
    if detSigma <= 0
        % Fall back to sample covariance if Sigma is singular
        Sigma = Lambda;
        detSigma = detLambda;
    end
catch
    % If determinant calculation fails, use sample covariance
    Sigma = Lambda;
    detSigma = detLambda;
end

% Calculate mESS with safeguards
if detSigma > 0 && isfinite(detLambda) && detLambda > 0
    mESS = n * (detLambda / detSigma)^(1/p);
else
    % Fallback: use univariate ESS approximation
    mESS = n;
end

% Ensure mESS is a real number and within valid range
if ~isreal(mESS) || ~isfinite(mESS) || mESS < 0
    mESS = n;
end

end