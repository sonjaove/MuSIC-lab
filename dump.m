% simulate_sampled_errors.m
% Simulate ẋ = -∇f - Lx - Lz - L(x̂ - x) - L(ẑ - z)
%         ż =  Lz + L(ẑ - z)
% using the same local costs and initial conditions as your original script.
clc; clear; close all;

%% Parameters
N = 3; d = 2;
T_final = 10; dt = 0.01;
tvec = 0:dt:T_final; steps = numel(tvec);
sigma = 0.05;                % event-trigger parameter (tune)

% Graph (triangular ring)
A = [0 1 1; 1 0 1; 1 1 0];
L = diag(sum(A,2)) - A;
neighbors = cell(N,1);
for i = 1:N, neighbors{i} = find(A(i,:) == 1); end

% Initial states (same as before)
x = [  2   2;
       3  -5;
      -4   1 ];   % N x 2
z = zeros(N,2);

% Sampled / last-broadcast values (start equal to initial)
x_last = x;
z_last = z;

% Storage
X_hist = zeros(N,d,steps);
Z_hist = zeros(N,d,steps);
trigger_count_x = zeros(N,1);
trigger_count_z = zeros(N,1);

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
        g = local_grad(i, x(i,:));   % 1x2
        
        % x_dot = -grad - Lx_i - Lz_i - L(x̂-x)_i - L(ẑ-z)_i
        x_dot(i,:) = -g - Lx(i,:) - Lz(i,:) - Le_x(i,:) - Le_z(i,:);
        
        % z_dot = Lz_i + L(ẑ - z)_i
        z_dot(i,:) = Lz(i,:) + Le_z(i,:);
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

%% Simple plots (same style as before)
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

%% Local gradient (kept at file end)
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
