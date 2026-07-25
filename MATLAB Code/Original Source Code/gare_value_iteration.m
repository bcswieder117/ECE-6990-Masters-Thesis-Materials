function result = gare_value_iteration(A,BP,BE,Q,RP,RE,P0,tolerance,max_iterations)
%GARE_VALUE_ITERATION Nominal model-based zero-sum value iteration.
%
% Isolates the Riccati value-iteration core adapted from Rizvi and Lin,
% Algorithm 3.2, using the thesis pursuer-evader notation. No disturbance
% experiments, calibration, compensation, plots, or exports are included.

arguments
    A double
    BP double
    BE double
    Q double
    RP double
    RE double
    P0 double
    tolerance (1,1) double {mustBePositive} = 1e-10
    max_iterations (1,1) double {mustBeInteger,mustBePositive} = 5000
end

nx = size(A,1);
validateattributes(A,{'double'},{'2d','square','finite'});
validateattributes(Q,{'double'},{'size',[nx nx],'finite'});
validateattributes(P0,{'double'},{'size',[nx nx],'finite'});
validateattributes(BP,{'double'},{'nrows',nx,'finite'});
validateattributes(BE,{'double'},{'nrows',nx,'finite'});

Q = (Q+Q')/2;
RP = (RP+RP')/2;
RE = (RE+RE')/2;
P = (P0+P0')/2;
mP = size(BP,2);
L = zeros(mP,nx);
K = zeros(size(BE,2),nx);

P_change = nan(max_iterations,1);
L_change = nan(max_iterations,1);
K_change = nan(max_iterations,1);
spectral_radius = nan(max_iterations,1);
game_rcond = nan(max_iterations,1);
converged = false;

for i = 1:max_iterations
    P_old = P;
    L_old = L;
    K_old = K;

    G_old = [RP+BP'*P_old*BP, BP'*P_old*BE;
             BE'*P_old*BP,    BE'*P_old*BE-RE];
    F_old = [BP'*P_old*A;
             BE'*P_old*A];

    check_game_blocks(G_old,BP,BE,P_old,RP,RE,i,'current');

    % Source Algorithm 3.2 Riccati optimality mapping.
    P_new = Q+A'*P_old*A-F_old'*(G_old\F_old);
    P_new = (P_new+P_new')/2;

    G_new = [RP+BP'*P_new*BP, BP'*P_new*BE;
             BE'*P_new*BP,    BE'*P_new*BE-RE];
    F_new = [BP'*P_new*A;
             BE'*P_new*A];

    check_game_blocks(G_new,BP,BE,P_new,RP,RE,i,'updated');

    gains = -(G_new\F_new);
    L_new = gains(1:mP,:);
    K_new = gains(mP+1:end,:);

    P_change(i) = norm(P_new-P_old,'fro');
    L_change(i) = norm(L_new-L_old,'fro');
    K_change(i) = norm(K_new-K_old,'fro');
    spectral_radius(i) = max(abs(eig(A+BP*L_new+BE*K_new)));
    game_rcond(i) = rcond(G_new);

    P = P_new;
    L = L_new;
    K = K_new;

    % Stricter implementation safeguard: matrix and both gains settle.
    if P_change(i) < tolerance && ...
       L_change(i) < tolerance && ...
       K_change(i) < tolerance
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
end

function check_game_blocks(G,BP,BE,P,RP,RE,iteration,label)
if rcond(G) <= 1e-13
    error('gare_value_iteration:singularGameMatrix', ...
        'The %s game matrix is singular at iteration %d.',label,iteration);
end
HPP = RP+BP'*P*BP;
SEE = BE'*P*BE-RE-(BE'*P*BP)*(HPP\(BP'*P*BE));
if min(real(eig((HPP+HPP')/2))) <= 0 || ...
   max(real(eig((SEE+SEE')/2))) >= 0
    error('gare_value_iteration:invalidSaddleCurvature', ...
        'The %s saddle curvature fails at iteration %d.',label,iteration);
end
end
