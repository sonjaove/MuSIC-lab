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

% deviding the simulations into different functions.
% euler updates function
function out = simulate_euler(dyn_fun, trig_fun, params)

% unpack
N = params.N; dt = params.dt; steps = params.steps;
x = params.x; z = params.z;
x_last = params.x_last; z_last = params.z_last;
L = params.L; neighbors = params.neighbors;
sigma = params.sigma;

for k = 1:steps
    % Laplacians
    Lx = L*x;  Lz = L*z;
    Le_x = L*(x_last - x);
    Le_z = L*(z_last - z);

    % dynamics (delegated)
    [x_dot, z_dot] = dyn_fun(x,z,x_last,z_last,Lx,Lz,Le_x,Le_z,params);

    % triggering (delegated)
    [x_last,z_last] = trig_fun(x,z,x_last,z_last,neighbors,sigma);

    % Euler
    x = x + dt*x_dot;
    z = z + dt*z_dot;
end

out.x = x; out.z = z;
out.x_last = x_last; out.z_last = z_last;
end
% example dynamics (the one discussed on the board)
function [x_dot,z_dot] = dyn_example(x,z,~,~,Lx,Lz,Le_x,Le_z,params)

N = params.N;
x_dot = zeros(N,2);
z_dot = zeros(N,2);

for i = 1:N
    g = params.local_grad(i,x(i,:));
    x_dot(i,:) = -g - Lx(i,:) - Lz(i,:) - Le_x(i,:) - Le_z(i,:);
    z_dot(i,:) =  Lx(i,:) + Le_x(i,:);
end
end
% trigger example same as the last code.
function [x_last,z_last] = trigger_example(x,z,x_last,z_last,neighbors,sigma)

N = size(x,1);

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
    end
end
end

params = struct('N',N,'dt',dt,'steps',steps,'L',L,...
                'x',x,'z',z,'x_last',x_last,'z_last',z_last,...
                'neighbors',{neighbors},'sigma',sigma,...
                'local_grad',@local_grad);

out = simulate_euler(@dyn_example,@trigger_example,params);



% function to calculate grad.
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
