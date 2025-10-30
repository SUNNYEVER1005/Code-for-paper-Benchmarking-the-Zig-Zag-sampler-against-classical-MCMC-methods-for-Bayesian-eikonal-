function H = Hellinger(thetamin, thetamax, samples, T_gridip, ...
                       timeobs, Tx_gridip, Tz_gridip)

% Transpose samples to have shape (nSamples, nParams)
samples = samples';

% Grid resolution
nx = 100;
ny = 100;
nz = 100;

% Build evaluation grid
xlin = linspace(thetamin(1), thetamax(1), nx);
ylin = linspace(thetamin(2), thetamax(2), ny);
zlin = linspace(thetamin(3), thetamax(3), nz);
[X, Y, Z] = meshgrid(xlin, ylin, zlin);
pts = [X(:), Y(:), Z(:)];          % (nx*ny*nz, 3)

% Evaluate true (unnormalised) density at grid points
Z_true = arrayfun(@(i) Potential_energy(pts(i,:), T_gridip, ...
                                       timeobs, Tx_gridip, Tz_gridip), ...
                 1:size(pts,1));
Z_true = reshape(Z_true, size(X)); % (nx, ny, nz)

% Kernel-density estimate from MCMC samples
f_emp = mvksdensity(samples, pts); % (nx*ny*nz, 1)
Z_emp = reshape(f_emp, size(X));   % (nx, ny, nz)

% Flatten and convert to probabilities
u = Z_true(:);        % unnormalised negative log-density
p = exp(-u);          % unnormalised density
q = Z_emp(:);         % kernel density values

p = p / sum(p);       % normalise to sum to 1
q = q / sum(q);

% Hellinger distance
H = norm(sqrt(p) - sqrt(q)) / sqrt(2);
end