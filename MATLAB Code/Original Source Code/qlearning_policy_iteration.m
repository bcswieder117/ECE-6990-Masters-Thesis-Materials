function result = qlearning_policy_iteration(x,uP,uE,stage_cost,L0,K0,H0,tolerance,max_iterations)
%QLEARNING_POLICY_ITERATION Nominal off-policy quadratic Q-learning PI.
%
% Isolates the sampled Bellman policy-evaluation / saddle-policy-improvement
% core adapted from Rizvi and Lin, Algorithm 3.3, in thesis notation.
%
% Inputs contain the measured transition batch only:
%   x          : nx-by-(N+1) state sequence
%   uP, uE     : realized pursuer and evader actions, each m-by-N
%   stage_cost : N-by-1 signed zero-sum stage costs
%   L0, K0     : initial pursuer and evader gains
%   H0         : initial symmetric action-value matrix (diagnostic baseline)
%
% The analytical plant matrices are intentionally absent. Data generation,
% offset calibration, compensation, simulations, plots, and Monte Carlo
% experiments belong outside this function.

arguments
    x double
    uP double
    uE double
    stage_cost double
    L0 double
    K0 double
    H0 double
    tolerance (1,1) double {mustBePositive} = 1e-10
    max_iterations (1,1) double {mustBeInteger,mustBePositive} = 200
end

[nx,Np1] = size(x);
N = Np1-1;
if size(uP,2) ~= N || size(uE,2) ~= N || numel(stage_cost) ~= N
    error('qlearning_policy_iteration:batchSizeMismatch', ...
        'The state, action, and stage-cost batch sizes are inconsistent.');
end
stage_cost = stage_cost(:);
mP = size(uP,1);
mE = size(uE,1);
nz = nx+mP+mE;
nH = nz*(nz+1)/2;
if ~isequal(size(H0),[nz nz])
    error('qlearning_policy_iteration:invalidH0','H0 must be nz-by-nz.');
end

Phi = zeros(N,nH);
for k = 1:N
    Phi(k,:) = quadratic_features([x(:,k);uP(:,k);uE(:,k)]);
end
if rank(Phi) < nH
    error('qlearning_policy_iteration:insufficientExcitation', ...
        'The measured state-action feature matrix is not full column rank.');
end

H = (H0+H0')/2;
L = L0;
K = K0;

H_change = nan(max_iterations,1);
L_change = nan(max_iterations,1);
K_change = nan(max_iterations,1);
bellman_residual = nan(max_iterations,1);
evaluation_rank = nan(max_iterations,1);
evaluation_condition = nan(max_iterations,1);
game_rcond = nan(max_iterations,1);
pursuer_curvature = nan(max_iterations,1);
evader_curvature = nan(max_iterations,1);
converged = false;

for i = 1:max_iterations
    H_old = H;
    L_old = L;
    K_old = K;

    Psi = zeros(N,nH);
    for k = 1:N
        z_next_policy = [x(:,k+1);L_old*x(:,k+1);K_old*x(:,k+1)];
        Psi(k,:) = Phi(k,:)-quadratic_features(z_next_policy);
    end

    evaluation_rank(i) = rank(Psi);
    if evaluation_rank(i) < nH
        error('qlearning_policy_iteration:rankDeficientEvaluation', ...
            'The fixed-policy Bellman matrix is not full column rank.');
    end

    scale = sqrt(sum(Psi.^2,1));
    scale(scale<eps) = 1;
    Psi_scaled = Psi./scale;
    evaluation_condition(i) = cond(Psi_scaled);

    theta_scaled = Psi_scaled\stage_cost;
    theta = theta_scaled./scale';
    H_new = symmetric_matrix(theta,nz);

    [L_new,K_new,rcond_value,cp,ce] = saddle_gain(H_new,nx,mP,mE);

    H_change(i) = norm(H_new-H_old,'fro');
    L_change(i) = norm(L_new-L_old,'fro');
    K_change(i) = norm(K_new-K_old,'fro');
    bellman_residual(i) = norm(Psi*theta-stage_cost,2)/(1+norm(stage_cost,2));
    game_rcond(i) = rcond_value;
    pursuer_curvature(i) = cp;
    evader_curvature(i) = ce;

    H = H_new;
    L = L_new;
    K = K_new;

    if H_change(i)<tolerance && L_change(i)<tolerance && K_change(i)<tolerance
        converged = true;
        break
    end
end

iterations = i;
result = struct;
result.H = H;
result.L = L;
result.K = K;
result.converged = converged;
result.iterations = iterations;
result.H_change = H_change(1:iterations);
result.L_change = L_change(1:iterations);
result.K_change = K_change(1:iterations);
result.bellman_residual = bellman_residual(1:iterations);
result.evaluation_rank = evaluation_rank(1:iterations);
result.evaluation_condition = evaluation_condition(1:iterations);
result.game_rcond = game_rcond(1:iterations);
result.pursuer_curvature = pursuer_curvature(1:iterations);
result.evader_curvature = evader_curvature(1:iterations);
end

function phi = quadratic_features(z)
n = numel(z);
phi = zeros(1,n*(n+1)/2);
p = 0;
for r = 1:n
    for c = r:n
        p = p+1;
        if r==c
            phi(p) = z(r)*z(c);
        else
            phi(p) = 2*z(r)*z(c);
        end
    end
end
end

function H = symmetric_matrix(theta,n)
H = zeros(n);
p = 0;
for r = 1:n
    for c = r:n
        p = p+1;
        H(r,c) = theta(p);
        H(c,r) = theta(p);
    end
end
H = (H+H')/2;
end

function [L,K,rc,cp,ce] = saddle_gain(H,nx,mP,mE)
idxP = nx+(1:mP);
idxE = nx+mP+(1:mE);
idxU = [idxP idxE];
G = H(idxU,idxU);
F = H(idxU,1:nx);
rc = rcond(G);
if rc <= 1e-13
    error('qlearning_policy_iteration:singularGameMatrix', ...
        'The learned saddle game matrix is numerically singular.');
end
HPP = H(idxP,idxP);
HPE = H(idxP,idxE);
HEP = H(idxE,idxP);
HEE = H(idxE,idxE);
SEE = HEE-HEP*(HPP\HPE);
cp = min(real(eig((HPP+HPP')/2)));
ce = max(real(eig((SEE+SEE')/2)));
if cp <= 0 || ce >= 0
    error('qlearning_policy_iteration:invalidSaddleCurvature', ...
        'The learned Q-function does not define a convex-concave saddle.');
end
gains = -(G\F);
L = gains(1:mP,:);
K = gains(mP+1:mP+mE,:);
end
