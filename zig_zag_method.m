% Zig-Zag sampler
function [Time, Xi, Theta, potential, flip_count, m, M, i0, taui0] = ...
    zig_zag_method(initial_point, xmin, xmax, nparameters, nSamples, ...
                   timeobs, sigma_r, gamma, receiver, v0)

flip_count = 0;                           % total number of velocity flips

% ------------------------- storage allocation -------------------------
Time       = zeros(1, nSamples);        % event times
Xi         = zeros(nparameters, nSamples); % positions
Theta      = zeros(nparameters, nSamples); % velocities
potential  = zeros(1, nSamples);        % potential (negative log-target)

% ------------------------- initial state ------------------------------
Xi(:,1)    = initial_point;             % initial position
Theta(:,1) = -1 + 2 * randi([0 1], nparameters, 1); % random ±1 velocities

% ------------------------- sampling loop ------------------------------
for k = 1:nSamples
    k
    % --- current state ---
    pos = Xi(:,k);                       % current position
    vel = Theta(:,k);                    % current velocity

    % --- draw next event time via Cinlar method / bounding rates ---
    [~, tau] = Cinlar_method_and_bound(vel, 0, pos, timeobs, sigma_r, ...
                                       gamma, receiver, v0);

    % --- select earliest event ---
    [tau_i0, i_0] = min(tau);            % minimum time & coordinate index
    Time(k+1)     = Time(k) + tau_i0;
    Xi(:,k+1)     = pos + tau_i0 * vel;  % propagate deterministically

    % --- thinning: compute actual vs upper bound rates ---
    tempm_i0 = rate(tau_i0, pos, vel, timeobs, sigma_r, gamma, receiver, v0);
    m_i0     = tempm_i0(i_0);

    tempM_i0 = Cinlar_method_and_bound(vel, tau_i0, pos, timeobs, sigma_r, ...
                                       gamma, receiver, v0);
    M_i0     = tempM_i0(i_0);

    prob     = m_i0 / M_i0;              % acceptance probability

    % --- store diagnostics ---
    taui0(k) = tau_i0;
    m(k)     = m_i0;
    M(k)     = M_i0;
    i0(k)    = i_0;

    % --- accept / reject flip ---
    if rand < prob
        Theta(:,k+1)        = vel;
        Theta(i_0, k+1)     = -vel(i_0); % flip velocity component
        flip_count          = flip_count + 1;
    else
        Theta(:,k+1)        = vel;       % keep velocity unchanged
    end

    % --- periodic boundary handling (mirror & fold) --------------------
    for i = 1:nparameters
        L = xmax(i) - xmin(i);
        if L <= 0
            error('xmax must exceed xmin for parameter %d', i);
        end

        x   = Xi(i,k+1) - xmin(i);   % offset from lower bound
        n   = floor(x / L);          % number of full spans
        r   = mod(x, L);             % remainder inside one span

        if mod(n,2) == 0             % even → forward direction
            Xi(i,k+1) = xmin(i) + r;
        else                         % odd  → backward direction
            Xi(i,k+1) = xmax(i) - r;
        end

        % flip velocity once for every reflection
        if mod(n,2) == 1
            Theta(i,k+1) = -Theta(i,k+1);
        end
    end

    % --- evaluate potential (negative log-target) -----------------------
    potential(k+1) = Psi([Xi(1,k+1); Xi(2,k+1); Xi(3,k+1)], ...
                         timeobs, sigma_r, receiver, v0);
end

end