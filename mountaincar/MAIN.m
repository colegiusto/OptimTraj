% MAIN - Mountain Car
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

problem.bounds.state.low = [-1.2; -inf];
problem.bounds.state.upp = [0.7; inf];
problem.bounds.initialState.low = [-0.5;0];
problem.bounds.initialState.upp = [-0.5;0];
problem.bounds.finalState.low = [0.5;0];
problem.bounds.finalState.upp = [0.5;0];

problem.bounds.control.low = -1; %-inf;
problem.bounds.control.upp = 1; %inf;

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
q = soln.grid.state(1,:);
dq = soln.grid.state(2,:);
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
plot(t,q, t, xs(1,:))
ylabel('q')
title('Single Pendulum Swing-Up');

subplot(3,1,2)
plot(t,dq, t, xs(2,:))
ylabel('dq')

subplot(3,1,3)
plot(t,u)
ylabel('u')


%% animate the result

clf 



plot(linspace(-1.2, 0.6,100), sin(3*linspace(-1.2, 0.6, 100)))
hold on 
set(gca, 'YTickLabel', [ ]); 
lblTime = uicontrol('style','text');
lblAction = uicontrol('style','text');
set(lblTime,'Position', [10 20 40 20]);
set(lblAction,'Position', [10 50 40 20]);

carReal = plot(0,0, 'ob', 'LineWidth', 4);
car = plot(0,0, 'or', 'LineWidth', 4);

for i = 1:length(t)
    set(carReal, 'XData', xs(1,i));
    set(carReal, 'YData', sin(3 *xs(1,i))); 
    set(car, 'XData', q(i));
    set(car, 'YData', sin(3 *q(i)));

    set(lblTime,'String', i - 1);
    set(lblAction,'String', u(i));
    
    drawnow;
    pause(0.1);
end
%% LQR Contr