function result = gare_policy_iteration(A,BP,BE,Q,RP,RE,L0,K0,tolerance,max_iterations)
%GARE_POLICY_ITERATION Nominal model-based zero-sum policy iteration.
%
% This function isolates the nominal policy-iteration core adapted from
% Rizvi and Lin, Algorithm 3.1, to the thesis pursuer-evader notation:
%
%   x(k+1) = A*x(k) + BP*uP(k) + BE*uE(k)
%   uP(k)  = L*x(k)      minimizing pursuer
%   uE(k)  = K*x(k)      maximizing evader
%
%   stage cost = x'*Q*x + uP'*RP*uP - uE'*RE*uE.
%
% No command offsets, calibration, experiments, plots, Monte Carlo work,
% or file exports are included here.
%
% Inputs
%   A, BP, BE, Q, RP, RE : nominal zero-sum game matrices
%   L0, K0               : admissible initial feedback gains
%   tolerance            : inherited stopping tolerance on ||P_i-P_{i-1}||_F
%   max_iterations       : maximum policy-iteration updates
%
% Output
%   result.P, result.L, result.K and numerical verification diagnostics.

arguments
    A double
    BP double
    BE double
    Q double
    RP double
    RE double
    L0 double
    K0 double
    tolerance (1,1) double {mustBePositive} = 1e-7
    max_iterations (1,1) double {mustBeInteger,mustBePositive} = 2000
end

nx = size(A,1);
validateattributes(A,{'double'},{'2d','square','finite'});
validateattributes(Q,{'double'},{'size',[nx nx],'finite'});
validateattributes(BP,{'double'},{'nrows',nx,'finite'});
validateattributes(BE,{'double'},{'nrows',nx,'finite'});
validateattributes(L0,{'double'},{'ncols',nx,'nrows',size(BP,2),'finite'});
validateattributes(K0,{'double'},{'ncols',nx,'nrows',size(BE,2),'finite'});

Q = (Q+Q')/2;
RP = (RP+RP')/2;
RE = (RE+RE')/2;

L = L0;
K = K0;
P = zeros(nx);

if max(abs(eig(A+BP*L+BE*K))) >= 1
    error('gare_policy_iteration:inadmissibleInitialPolicy', ...
        'The initial pursuer-evader gain pair is not Schur stable.');
end

P_change = nan(max_iterations,1);
L_change = nan(max_iterations,1);
K_change = nan(max_iterations,1);
spectral_radius = nan(max_iterations,1);
game_rcond = nan(max_iterations,1);

converged = false;

for i = 1:max_iterations
    P_previous = P;
    L_previous = L;
    K_previous = K;

    % Source Algorithm 3.1: fixed-policy Lyapunov evaluation.
    Acl_i = A+BP*L_previous+BE*K_previous;
    if max(abs(eig(Acl_i))) >= 1
        error('gare_policy_iteration:unstablePolicy', ...
            'The policy at iteration %d is not Schur stable.',i);
    end

    S_i = Q + L_previous'*RP*L_previous ...
              - K_previous'*RE*K_previous;
    S_i = (S_i+S_i')/2;
    P = dlyap(Acl_i',S_i);
    P = (P+P')/2;

    % Coupled saddle-policy improvement (joint block form, algebraically
    % equivalent to the source Schur-complement equations).
    G = [RP+BP'*P*BP, BP'*P*BE;
         BE'*P*BP,    BE'*P*BE-RE];
    F = [BP'*P*A;
         BE'*P*A];

    if rcond(G) <= 1e-13
        error('gare_policy_iteration:singularGameMatrix', ...
            'The game matrix is numerically singular at iteration %d.',i);
    end

    HPP = RP+BP'*P*BP;
    SEE = BE'*P*BE-RE-(BE'*P*BP)*(HPP\(BP'*P*BE));
    if min(real(eig((HPP+HPP')/2))) <= 0 || ...
       max(real(eig((SEE+SEE')/2))) >= 0
        error('gare_policy_iteration:invalidSaddleCurvature', ...
            'The saddle curvature conditions fail at iteration %d.',i);
    end

    gains = -(G\F);
    mP = size(BP,2);
    L = gains(1:mP,:);
    K = gains(mP+1:end,:);

    P_change(i) = norm(P-P_previous,'fro');
    L_change(i) = norm(L-L_previous,'fro');
    K_change(i) = norm(K-K_previous,'fro');
    spectral_radius(i) = max(abs(eig(A+BP*L+BE*K)));
    game_rcond(i) = rcond(G);

    % Inherited Algorithm 3.1 stopping rule.
    if P_change(i) < tolerance
        converged = true;
        break
    end
end

iterations = i;
Acl = A+BP*L+BE*K;
G = [RP+BP'*P*BP, BP'*P*BE;
     BE'*P*BP,    BE'*P*BE-RE];
F = [BP'*P*A;
     BE'*P*A];
T_P = Q+A'*P*A-F'*(G\F);
lyapunov_rhs = Q+L'*RP*L-K'*RE*K+Acl'*P*Acl;
HPP = RP+BP'*P*BP;
SEE = BE'*P*BE-RE-(BE'*P*BP)*(HPP\(BP'*P*BE));

result = struct;
result.P = P;
result.L = L;
result.K = K;
result.converged = converged;
result.iterations = iterations;
result.P_change = P_change(1:iterations);
result.L_change = L_change(1:iterations);
result.K_change = K_change(1:iterations);
result.spectral_radius_history = spectral_radius(1:iterations);
result.game_rcond_history = game_rcond(1:iterations);
result.closed_loop_spectral_radius = max(abs(eig(Acl)));
result.game_matrix_rcond = rcond(G);
result.pursuer_curvature_minimum = min(real(eig((HPP+HPP')/2)));
result.evader_curvature_maximum = max(real(eig((SEE+SEE')/2)));
result.gare_residual = norm(P-(T_P+T_P')/2,'fro');
result.lyapunov_residual = norm(P-(lyapunov_rhs+lyapunov_rhs')/2,'fro');
end
