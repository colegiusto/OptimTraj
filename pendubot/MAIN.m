% MAIN - Pendubot
%

clc; clear;
addpath ../

% Physical parameters of the pendulum
% p.g = 0.5;  % Normalized gravity constant
% p.m = 2;
% p.mtn = @(x)atan(3*cos(3*x));  % Normalized damping constant

load guess.mat

% User-defined dynamics and objective functions
problem.func.dynamics = @(t,x,u)( dynamics(x,u,p) );
problem.func.pathObj = @(t,x,u)( u.^2 );

% Problem bounds
problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = 5;
problem.bounds.finalTime.upp = 30;

problem.bounds.state.low = [-inf; -inf; -inf; -inf];
problem.bounds.state.upp = [inf; inf; inf; inf];
problem.bounds.initialState.low = [-pi/2;0;0;0];
problem.bounds.initialState.upp = [-pi/2;0;0;0];
problem.bounds.finalState.low = [pi/2;0;0;0];
problem.bounds.finalState.upp = [pi/2;0;0;0];

problem.bounds.control.low = -inf; %-inf;
problem.bounds.control.upp = inf; %inf;

% Guess at the initial trajectory



problem.guess = guess;

% problem.guess.time = [0,25];
% problem.guess.state = [-0.5, 0.5; 0, 0];
% problem.guess.control = [0,0];

% t = linspace(0,20,20);
% problem.guess.time = t;
% problem.guess.state = [t.*sin(t*pi/40*9)/20-0.5; 0, diff(t.*sin(t*pi/40*9)/20-0.5)];
% problem.guess.control = [1, diff(problem.guess.state(2,:))];


% Select a solver:
problem.options.method = 'trapezoid';
problem.options.defaultAccuracy = 'high';

% problem.options.nlpOpt = optimoptions('fmincon');
% problem.options.nlpOpt.MaxIterations = 1000;

% Solve the problem
soln = optimTraj(problem);
t = soln.grid.time;
q1 = soln.grid.state(1,:);
q2 = soln.grid.state(2,:);

dq1 = soln.grid.state(3,:);
dq2 = soln.grid.state(4,:);
u = soln.grid.control;

%%

x0 = soln.grid.state(:,1);
xs = zeros(size(soln.grid.state));
xs(:,1) = x0;
for i=2:length(t)
    so = ode45(@(T,x)dynamics(x,u(i)+(u(i)-u(i-1))*(T-t(i-1))/(t(i)-t(i-1)),p), linspace(t(i-1), t(i), 100), xs(:,i-1));
    xs(:,i) = so.y(:, end);
end

% Plot the solution:
figure(1); clf;

subplot(3,1,1)
plot(t,q2, t, xs(2,:))
ylabel('q')
title('Pendulum');

subplot(3,1,2)
plot(t,dq2, t, xs(4,:))
ylabel('dq')

subplot(3,1,3)
plot(t,u)
ylabel('u')


%% animate the result
figure(3)
clf 

hold on 
set(gca, 'YTickLabel', [ ]); 
lblTime = uicontrol('style','text');
lblAction = uicontrol('style','text');
set(lblTime,'Position', [10 20 40 20]);
set(lblAction,'Position', [10 50 40 20]);


for i = 1:length(t)
    clf; hold on;
    animate(t(i), soln.grid.state(:,i), p)


    % set(lblTime,'String', i - 1);
    % set(lblAction,'String', u(i));
    
    drawnow;
    pause(0.05);
end
%% LQR Control

J = @(x,h,F)(F(repmat(x,size(x'))+diag(h))-F(repmat(x,size(x'))))./h';
xsol = soln.grid.state;



K = zeros(1,4,length(t));

Q = diag([0.01, 7, 0.1, 0.9]);
R = 1e-2;

for i = 1:length(t)
    A = J(xsol(:,i), 1e-5*ones(4,1), @(x)dynamics(x, u(i), p));

    B = J(u(i), 1e-5, @(u)dynamics(xsol(:,i), u, p));

    K(:,:, i) = lqr(A, B, Q, R);
end

%% Simulate controlled dynamics

F_closed_loop = @(x, x_star, u_star, K)dynamics(x, clip(u_star-K*(x-x_star), -1, 1), p);

x_cl = zeros(size(soln.grid.state));
x_cl(:,1) = x0;




for i=2:length(t)

    sol = ode45(@(T,x)F_closed_loop(...
        x, ...
        xsol(:,i-1)+(xsol(:,i)-xsol(:,i-1))*(T-t(i-1))/(t(i)-t(i-1)), ...
        u(i-1)+(u(i)-u(i-1))*(T-t(i-1))/(t(i)-t(i-1)), ...
        K(:,:,i-1)) , ...
        [t(i-1), t(i)], ...
        x_cl(:,i-1)...
        );
    x_cl(:,i) = sol.y(:,end);
    
end

plot(t, xsol(2,:), t, x_cl(2,:))

%% Plot error

figure(2)
clf
subplot(2,1,1)
plot(t, xsol);
hold on
plot(t, xs)
plot(t, x_cl)


subplot(2,1,2)

plot(t, xs-xsol)
hold on
plot(t, x_cl-xsol)


%% reanimate

 
figure(3)
clf
for i = 1:length(t)
    clf; hold on;
    animate(t(i), soln.grid.state(:,i), p)
    animate(t(i), x_cl(:,i), p)
    % animate(t(i), xs(:,i), p)



    % set(lblTime,'String', i - 1);
    % set(lblAction,'String', u(i));
    
    drawnow;
    pause(0.05);
end
