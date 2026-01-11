% simulate_sampled_errors_fixed.m
clc; clear; close all;

%% Parameters
N = 3; d = 2;
T_final = 100; dt = 0.01;
tvec = 0:dt:T_final; steps = numel(tvec);
sigma = 0.05;                % event-trigger parameter (tune)

% Graph (triangular ring)
A = [0 1 1; 1 0 1; 1 1 0];
L = diag(sum(A,2)) - A;
neighbors = cell(N,1);
for i = 1:N, neighbors{i} = find(A(i,:) == 1); end

% Initial states (same as before)
x = [  20   2;
       3  -15;
      -14   10 ];   % N x 2
z = zeros(N,2);

% Sampled / last-broadcast values (start equal to initial)
x_last = x;
z_last = z;

% Storage
X_hist = zeros(N,d,steps);
Z_hist = zeros(N,d,steps);
trigger_count_x = zeros(N,1);
trigger_count_z = zeros(N,1);

%% --- symbolic gradients for all 3 agents (compute once) ---
% requires Symbolic Math Toolbox
syms sx1 sx2
vars = [sx1 sx2];

% Smooth choices
f1 = (sx1 - 1)^2 + (sx2 + 2)^2;    % quadratic
f2 = exp(sx1 + sx2);               % smooth exp
f3 = sin(sx1) + cos(sx2);          % smooth trig

g1_sym = gradient(f1, vars);
g2_sym = gradient(f2, vars);
g3_sym = gradient(f3, vars);

% create numeric function handles that accept a 1x2 or 2x1 vector
grad1 = matlabFunction(g1_sym, 'Vars', {vars});   % returns 2x1 by default
grad2 = matlabFunction(g2_sym, 'Vars', {vars});
grad3 = matlabFunction(g3_sym, 'Vars', {vars});

%% Simulation (explicit Euler)
for k = 1:steps
    t = tvec(k);
    X_hist(:,:,k) = x;
    Z_hist(:,:,k) = z;
    
    % Precompute some Laplacians
    Lx = L * x;    % Nx2
    Lz = L * z;    % Nx2
    % error Laplacians
    Le_x = L * (x_last - x);   % L(x̂ - x)
    Le_z = L * (z_last - z);   % L(ẑ - z)
    
    % Dynamics
    x_dot = zeros(N,2);
    z_dot = zeros(N,2);
    for i = 1:N
    %     % call the appropriate gradient handle and ensure row vector 1x2
    %     switch i
    %         case 1
    %             gv = grad1(x(i,:));   % 2x1 or 1x2
    %         case 2
    %             gv = grad2(x(i,:));
    %         case 3
    %             gv = grad3(x(i,:));
    %         otherwise
    %             gv = [0;0];
    %     end
        %g = reshape(gv,1,2);% ensure 1x2 row vector
        g = local_grad(i, x(i,:));   
        
        % x_dot = -grad - Lx_i - Lz_i - L(x̂-x)_i - L(ẑ-z)_i
        x_dot(i,:) = -1 * g - Lx(i,:) - Lz(i,:) - 1*Le_x(i,:) - 1*Le_z(i,:);
        
        % z_dot = Lz_i + L(ẑ - z)_i
        z_dot(i,:) = Lx(i,:) + 1*Le_x(i,:);
    end
    
    % Event-trigger checks for x and z (decentralized)
    for i = 1:N
        % x trigger
        e_x = x_last(i,:) - x(i,:);
        local_error_x = dot(e_x,e_x);
        % threshold built from sampled neighbor differences of x_last
        thresh_x = 0;
        for j = neighbors{i}
            diff = x_last(i,:) - x_last(j,:);
            thresh_x = thresh_x + dot(diff,diff);
        end
        if local_error_x > sigma * thresh_x + 1e-12
            x_last(i,:) = x(i,:);
            trigger_count_x(i) = trigger_count_x(i) + 1;
        end
        
        % z trigger
        e_z = z_last(i,:) - z(i,:);
        local_error_z = dot(e_z,e_z);
        thresh_z = 0;
        for j = neighbors{i}
            diffz = z_last(i,:) - z_last(j,:);
            thresh_z = thresh_z + dot(diffz,diffz);
        end
        if local_error_z > sigma * thresh_z + 1e-12
            z_last(i,:) = z(i,:);
            trigger_count_z(i) = trigger_count_z(i) + 1;
        end
    end
    
    % Euler update
    x = x + dt * x_dot;
    z = z + dt * z_dot;
end

%% Pairwise disagreement plot: ||x_i - x_j||
figure; hold on;
pairs = [1 2; 2 3; 1 3];
for p = 1:size(pairs,1)
    i = pairs(p,1);
    j = pairs(p,2);
    diff_ij = squeeze(X_hist(i,:,:) - X_hist(j,:,:));   % 2 x T
    d_ij = vecnorm(diff_ij);                            % 1 x T
    plot(tvec, d_ij, 'LineWidth', 1.4);
end
xlabel('Time'); ylabel('||x_i - x_j||'); title('Pairwise Disagreement');
legend('1-2','2-3','1-3','Location','best'); grid on; hold off;

%% Other plots (kept as before)
figure('Units','normalized','Position',[0.1 0.1 0.6 0.6]);
subplot(3,1,1); hold on;
for i = 1:N, plot(tvec, squeeze(X_hist(i,1,:)), 'LineWidth',1.2); end
xlabel('Time'); ylabel('x_1'); title('Component 1'); legend('A1','A2','A3','Location','best'); grid on; hold off;

subplot(3,1,2); hold on;
for i = 1:N, plot(tvec, squeeze(X_hist(i,2,:)), 'LineWidth',1.2); end
xlabel('Time'); ylabel('x_2'); title('Component 2'); legend('A1','A2','A3','Location','best'); grid on; hold off;

subplot(3,1,3); hold on; colors = lines(N);
for i = 1:N
    plot(squeeze(X_hist(i,1,:)), squeeze(X_hist(i,2,:)), 'Color', colors(i,:), 'LineWidth',1.0);
    plot(X_hist(i,1,1), X_hist(i,2,1), 's', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
    plot(X_hist(i,1,end), X_hist(i,2,end), 'd', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
end
xlabel('x_1'); ylabel('x_2'); title('2D trajectories'); legend('A1','A2','A3','Location','best'); grid on; axis equal; hold off;

%% Print trigger counts
fprintf('\nTriggers for x broadcasts per agent:\n');
for i = 1:N, fprintf('Agent %d: %d\n', i, trigger_count_x(i)); end
fprintf('\nTriggers for z broadcasts per agent:\n');
for i = 1:N, fprintf('Agent %d: %d\n', i, trigger_count_z(i)); end
function g = local_grad(node, v)
    x1 = v(1); x2 = v(2);
    switch node
        case 1
            if x1 > x2, g = [1,0];
            elseif x2 > x1, g = [0,1];
            else g = [0.5,0.5]; end
        case 2
            s = exp(x1 + x2);
            g = [2*x1 + s, 2*x2 + s];
        case 3
            s = sinh(x1 + x2);
            g = [s, s];
        otherwise
            g = [0,0];
    end
end
