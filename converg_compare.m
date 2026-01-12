% simulate_sampled_errors_fixed.m
clc; clear; close all;

%% Parameters
N = 4; d = 2;
T_final = 100; dt = 0.01;
tvec = 0:dt:T_final; steps = numel(tvec);
sigma = 0.05;

% Graph (triangular ring)
A = [0 1 0 1;
     1 0 1 0;
     0 1 0 1;
     1 0 1 0];   % square graph
L = diag(sum(A,2)) - A;

neighbors = cell(N,1);
for i = 1:N, neighbors{i} = find(A(i,:) == 1); end

% Initial states
x = [20 2; 3 -15; -14 10];
z = zeros(N,2);
x_last = x;
z_last = z;

function dmax = disagreement_metric(X)
    N = size(X,1);
    dmax = 0;
    for i = 1:N
        for j = i+1:N
            dmax = max(dmax, norm(X(i,:) - X(j,:)));
        end
    end
end


% Storage
X_hist = zeros(N,d,steps);
Z_hist = zeros(N,d,steps);
trigger_count_x = zeros(N,1);
trigger_count_z = zeros(N,1);

%% Main simulation
params = struct('N',N,'dt',dt,'steps',steps,'L',L,...
    'x',x,'z',z,'x_last',x_last,'z_last',z_last,...
    'neighbors',{neighbors},'sigma',sigma,...
    'local_grad',@local_grad);

out = simulate_euler(@dyn_example, @trigger_example, params);

% Extract results
X_hist = out.X_hist;
Z_hist = out.Z_hist;
trigger_count_x = out.trigger_count_x;
trigger_count_z = out.trigger_count_z;

%% Plots
% Pairwise disagreement
plot_pairwise_disagreement(X_hist, tvec, [1 2; 2 3; 1 3]);

% Component-wise time plots
plot_components(X_hist, tvec);

% 2D trajectories
plot_trajectories_2D(X_hist);

% Trigger statistics
print_triggers(trigger_count_x, trigger_count_z);

%% ============ FUNCTIONS ============

function out = simulate_euler(dyn_fun, trig_fun, params)
    % Unpack
    N = params.N; dt = params.dt; steps = params.steps;
    x = params.x; z = params.z;
    x_last = params.x_last; z_last = params.z_last;
    L = params.L; neighbors = params.neighbors;
    sigma = params.sigma;
    
    % Initialize storage
    X_hist = zeros(N, 2, steps);
    Z_hist = zeros(N, 2, steps);
    trigger_count_x = zeros(N,1);
    trigger_count_z = zeros(N,1);
    
    for k = 1:steps
        % Store current state
        X_hist(:,:,k) = x;
        Z_hist(:,:,k) = z;
        
        % Laplacians
        Lx = L*x;  Lz = L*z;
        Le_x = L*(x_last - x);
        Le_z = L*(z_last - z);
        
        % Dynamics
        [x_dot, z_dot] = dyn_fun(x,z,x_last,z_last,Lx,Lz,Le_x,Le_z,params);
        
        % Triggering (returns trigger flags too)
        [x_last, z_last, trig_x, trig_z] = trig_fun(x,z,x_last,z_last,neighbors,sigma);
        
        % Count triggers
        trigger_count_x = trigger_count_x + trig_x;
        trigger_count_z = trigger_count_z + trig_z;
        
        % Euler step
        x = x + dt*x_dot;
        z = z + dt*z_dot;
    end
    
    % Pack output
    out.x = x; out.z = z;
    out.x_last = x_last; out.z_last = z_last;
    out.X_hist = X_hist;
    out.Z_hist = Z_hist;
    out.trigger_count_x = trigger_count_x;
    out.trigger_count_z = trigger_count_z;
end

function [x_dot, z_dot] = dyn_example(x,z,~,~,Lx,Lz,Le_x,Le_z,params)
    N = params.N;
    x_dot = zeros(N,2);
    z_dot = zeros(N,2);
    
    for i = 1:N
        g = params.local_grad(i, x(i,:));
        x_dot(i,:) = -g - Lx(i,:) - Lz(i,:) - Le_x(i,:) - Le_z(i,:);
        z_dot(i,:) = Lx(i,:) + Le_x(i,:);
    end
end

function [x_last, z_last, trig_x, trig_z] = trigger_example(x,z,x_last,z_last,neighbors,sigma)
    N = size(x,1);
    trig_x = zeros(N,1);
    trig_z = zeros(N,1);
    
    for i = 1:N
        % x-trigger
        ex = x_last(i,:) - x(i,:);
        thx = 0;
        for j = neighbors{i}
            d = x_last(i,:) - x_last(j,:);
            thx = thx + dot(d,d);
        end
        if dot(ex,ex) > sigma*thx
            x_last(i,:) = x(i,:);
            trig_x(i) = 1;
        end
        
        % z-trigger
        ez = z_last(i,:) - z(i,:);
        thz = 0;
        for j = neighbors{i}
            d = z_last(i,:) - z_last(j,:);
            thz = thz + dot(d,d);
        end
        if dot(ez,ez) > sigma*thz
            z_last(i,:) = z(i,:);
            trig_z(i) = 1;
        end
    end
end

function g = local_grad(node, v)
    x1 = v(1); x2 = v(2);
    switch node
        case 1
            if x1 > x2, g = [1,0];
            elseif x2 > x1, g = [0,1];
            else, g = [0.5,0.5]; end
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

%% ============ PLOTTING FUNCTIONS ============

function plot_pairwise_disagreement(X_hist, tvec, pairs)
    figure('Position', [100 100 800 500]);
    hold on; grid on;
    
    colors = lines(size(pairs,1));
    
    for p = 1:size(pairs,1)
        i = pairs(p,1);
        j = pairs(p,2);
        
        % Compute ||x_i - x_j||^2 over time
        disagreement = squeeze(sum((X_hist(i,:,:) - X_hist(j,:,:)).^2, 2));
        
        plot(tvec, disagreement, 'LineWidth', 2, 'Color', colors(p,:), ...
            'DisplayName', sprintf('Agents %d-%d', i, j));
    end
    
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('||x_i - x_j||^2', 'FontSize', 12);
    title('Pairwise Disagreement over Time', 'FontSize', 14);
    legend('Location', 'best');
    set(gca, 'FontSize', 11);
end

function plot_components(X_hist, tvec)
    N = size(X_hist, 1);
    
    figure('Position', [100 100 1200 400]);
    
    % Component 1
    subplot(1,2,1);
    hold on; grid on;
    for i = 1:N
        plot(tvec, squeeze(X_hist(i,1,:)), 'LineWidth', 2, ...
            'DisplayName', sprintf('Agent %d', i));
    end
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('x_1', 'FontSize', 12);
    title('Component 1 vs Time', 'FontSize', 14);
    legend('Location', 'best');
    
    % Component 2
    subplot(1,2,2);
    hold on; grid on;
    for i = 1:N
        plot(tvec, squeeze(X_hist(i,2,:)), 'LineWidth', 2, ...
            'DisplayName', sprintf('Agent %d', i));
    end
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('x_2', 'FontSize', 12);
    title('Component 2 vs Time', 'FontSize', 14);
    legend('Location', 'best');
end

function plot_trajectories_2D(X_hist)
    N = size(X_hist, 1);
    
    figure('Position', [100 100 700 700]);
    hold on; grid on; axis equal;
    
    colors = lines(N);
    
    for i = 1:N
        traj = squeeze(X_hist(i,:,:))';
        plot(traj(:,1), traj(:,2), 'LineWidth', 2, 'Color', colors(i,:), ...
            'DisplayName', sprintf('Agent %d', i));
        
        % Mark start and end
        plot(traj(1,1), traj(1,2), 'o', 'MarkerSize', 10, ...
            'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', ...
            'HandleVisibility', 'off');
        plot(traj(end,1), traj(end,2), 's', 'MarkerSize', 10, ...
            'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', ...
            'HandleVisibility', 'off');
    end
    
    xlabel('x_1', 'FontSize', 12);
    ylabel('x_2', 'FontSize', 12);
    title('2D Trajectories (o = start, square = end)', 'FontSize', 14);
    legend('Location', 'best');
    set(gca, 'FontSize', 11);
end

function print_triggers(trigger_count_x, trigger_count_z)
    fprintf('\n=== Trigger Statistics ===\n');
    fprintf('Agent | x-triggers | z-triggers\n');
    fprintf('------|------------|------------\n');
    for i = 1:length(trigger_count_x)
        fprintf('  %d   |   %5d    |   %5d\n', i, trigger_count_x(i), trigger_count_z(i));
    end
    fprintf('\nTotal x-triggers: %d\n', sum(trigger_count_x));
    fprintf('Total z-triggers: %d\n', sum(trigger_count_z));
end