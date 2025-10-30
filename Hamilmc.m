%% Leapfrog HMC
function [value1, value2, value3, value4, value5, value6] = ...
         Hamilmc(initial_point, delta, L, nparameters, nSamples, ...
                 timeobs, receiver, M, v0)

% Initialize storage
theta = zeros(nparameters, nSamples);            % parameter samples
Potential_energy_result = zeros(1, nSamples);    % potential energy trace
Trajectory_X = zeros(L, nSamples);               % X-position trace
Trajectory_Y = zeros(L, nSamples);               % Y-position trace
Gradient = zeros(nparameters, nSamples);         % gradient trace

% Set starting point
theta(:, 1) = initial_point;

% Sampling loop
t = 1;
tic;
accept_ratio = 0;

while t < nSamples
    t = t + 1;
    t
    % Current position
    current_theta = theta(:, t-1);
    theta0 = current_theta;

    % Sample momentum from multivariate normal
    p = mvnrnd(zeros(nparameters, 1), M, 1)';

    current_p = p;

    %% Hamiltonian dynamics (leapfrog integrator)

    % Half step for momentum
    p = p - (delta/2) * gradient_P(theta0, timeobs, receiver, v0);

    % Alternate full steps for position and momentum
    for i = 1:L
        % Full step for position
        theta0 = theta0 + delta * (M \ p);
        Trajectory_X(i, t) = theta0(1);
        Trajectory_Y(i, t) = theta0(2);

        % Full step for momentum (except at the end)
        if i ~= L
            p = p - delta * gradient_P(theta0, timeobs, receiver, v0);
        end
    end

    % Half step for momentum at the end
    p = p - (delta/2) * gradient_P(theta0, timeobs, receiver, v0);

    % Negate momentum for symmetric proposal
    p = -p;

    % Store gradient
    Gradient(:, t) = gradient_P(theta0, timeobs, receiver, v0);

    %% Acceptance / rejection step

    % Potential energies
    current_U  = Potential_energy(current_theta, timeobs, receiver, v0);
    proposed_U = Potential_energy(theta0, timeobs, receiver, v0);

    % Kinetic energies
    current_K  = Kinetic_energy(current_p, M);
    proposed_K = Kinetic_energy(p, M);

    % Acceptance probability
    alpha = exp(current_U - proposed_U + current_K - proposed_K);
    u = rand;

    if u < alpha
        theta(:, t) = theta0;
        accept_ratio = accept_ratio + 1;
    else
        theta(:, t) = current_theta;
    end

    % Store potential energy of accepted state
    Potential_energy_result(t) = Potential_energy(theta(:, t), ...
                                                 timeobs, receiver, v0);
end

accept_ratio = accept_ratio / nSamples;

% Return values
value1 = theta;               % samples
value2 = accept_ratio;        % acceptance rate
value3 = Potential_energy_result;
value4 = Trajectory_X;
value5 = Trajectory_Y;
value6 = Gradient;
end