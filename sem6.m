clc; clear; close all;

%% Parameters
N = 4; d = 2;
T_final = 100; dt = 0.001;
tvec = 0:dt:T_final; steps = numel(tvec);
gamma1 = 1;
gamma2 = 1;
k_alpha = 1;

% Graph (square)
A = [0 1 0 1;
     1 0 1 0;
     0 1 0 1;
     1 0 1 0];
L = diag(sum(A,2)) - A;

% Initial states
k = 1;
x = k * [ 1   2;
         -2   1;
          3  -4;
         -2   1];
v = zeros(N,2);
alpha = zeros(N,1);

%% Target points for each node (what each agent "wants" to minimize toward)
targets = [ 2   1;   % node 1 wants to go to (2,1)
           -1   3;   % node 2 wants to go to (-1,3)
            1  -1;   % node 3 wants to go to (1,-1)
           -2   2];  % node 4 wants to go to (-2,2)

%% Pack parameters
params = struct('N',N,'dt',dt,'steps',steps,'L',L, ...
    'local_grad',@(node,xi) local_grad(node,xi,targets), ...
    'gamma1',gamma1,'gamma2',gamma2,'k_alpha',k_alpha);

%% Run simulation
out = simulate_euler(@dyn_adaptive, params, x, v, alpha);
X_hist = out.X_hist;

%% Plots
plot_pairwise_disagreement(X_hist, tvec, [1 2; 2 3; 3 4; 4 1; 1 3; 2 4]);
plot_components(X_hist, tvec);
plot_trajectories_2D(X_hist, targets);

%% ================= SIMULATOR =================
function out = simulate_euler(dyn_fun, params, x, v, alpha)
N     = params.N;
dt    = params.dt;
steps = params.steps;
L     = params.L;

X_hist = zeros(N, 2, steps);

for k = 1:steps
    X_hist(:,:,k) = x;
    Lx = L * x;
    [x_dot, v_dot, alpha_dot] = dyn_fun(x, v, alpha, Lx, params);
    x     = x     + dt * x_dot;
    v     = v     + dt * v_dot;
    alpha = alpha + dt * alpha_dot;
end
out.X_hist = X_hist;
end

%% ================= DYNAMICS =================
function [x_dot, v_dot, alpha_dot] = dyn_adaptive(x, v, alpha, Lx, params)
N       = params.N;
gamma1  = params.gamma1;
gamma2  = params.gamma2;
k_alpha = params.k_alpha;
L       = params.L;

Lv = L * v;

x_dot     = zeros(N, 2);
v_dot     = zeros(N, 2);
alpha_dot = zeros(N, 1);

for i = 1:N
    g   = params.local_grad(i, x(i,:));
    e_i = Lx(i,:);
    beta_i = dot(e_i, e_i);
    gain   = alpha(i) + beta_i;

    % Kinematic chain
    x_dot(i,:) = v(i,:);

    % All control in v_dot
    v_dot(i,:) = -gamma2 * g ...
                 - gamma1 * gain * e_i ...
                 - Lv(i,:);

    % Adaptive gain update
    beta_clip    = min(beta_i, 10);
    alpha_dot(i) = -k_alpha * alpha(i) + beta_clip;
end
end

%% ================= GRADIENTS =================
% All cost functions are simple quadratics:
%   f_i(x) = ||x - target_i||^2
%   grad   = 2*(x - target_i)
%
% Node 1: f1 = (x1-2)^2  + (x2-1)^2
% Node 2: f2 = (x1+1)^2  + (x2-3)^2
% Node 3: f3 = (x1-1)^2  + (x2+1)^2
% Node 4: f4 = (x1+2)^2  + (x2-2)^2
%
% Global minimizer = average of targets = (0, 1.25)
function g = local_grad(node, xi, targets)
g = 2 * (xi - targets(node,:));
end

%% ================= PLOT FUNCTIONS =================
function plot_pairwise_disagreement(X_hist, tvec, pairs)
figure; hold on; grid on;
for p = 1:size(pairs,1)
    i = pairs(p,1); j = pairs(p,2);
    dij = squeeze(sum((X_hist(i,:,:) - X_hist(j,:,:)).^2, 2));
    plot(tvec, dij, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%d-%d', i, j));
end
xlabel('Time');
ylabel('||x_i - x_j||^2');
title('Pairwise Disagreement');
legend show;
end

function plot_components(X_hist, tvec)
N = size(X_hist, 1);

% Analytical consensus point = mean of targets
target_mean = [0, 1.25];

figure;
subplot(1,2,1); hold on; grid on;
for i = 1:N
    plot(tvec, squeeze(X_hist(i,1,:)), 'LineWidth', 1.5);
end
yline(target_mean(1), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Consensus x_1');
title('Component x_1'); xlabel('Time'); legend show;

subplot(1,2,2); hold on; grid on;
for i = 1:N
    plot(tvec, squeeze(X_hist(i,2,:)), 'LineWidth', 1.5);
end
yline(target_mean(2), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Consensus x_2');
title('Component x_2'); xlabel('Time'); legend show;
end

function plot_trajectories_2D(X_hist, targets)
N = size(X_hist, 1);
colors = lines(N);
figure; hold on; grid on; axis equal;

% Analytical consensus point
consensus = mean(targets, 1);

for i = 1:N
    traj = squeeze(X_hist(i,:,:))';
    plot(traj(:,1), traj(:,2), 'Color', colors(i,:), 'LineWidth', 1.5);

    % Start marker
    plot(traj(1,1), traj(1,2), 'o', 'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), 'MarkerSize', 7);

    % End marker
    plot(traj(end,1), traj(end,2), 's', 'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), 'MarkerSize', 7);

    % Each node's local target
    plot(targets(i,1), targets(i,2), 'x', 'Color', colors(i,:), ...
        'MarkerSize', 10, 'LineWidth', 2);
end

% Consensus point
plot(consensus(1), consensus(2), 'kp', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 12, 'DisplayName', 'Consensus point');

legend({'','','','','','','','Node targets','','','','Consensus'}, ...
    'Location','best');
title('2D Trajectories  (o=start,  s=end,  x=local target,  p=consensus)');
end