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
problem.bounds.finalTime.upp = 20;

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
problem.options.trapezoid.nGrid = 100;
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

Q = diag([0.01, 8, 0.1, 1.3]);
R = 6;

for i = 1:length(t)
    A = J(xsol(:,i), 1e-5*ones(4,1), @(x)dynamics(x, 0.5*(u(i)+u(i)), p));

    B = J(u(i), 1e-5, @(u)dynamics(0.5*(xsol(:,i)+xsol(:,i)), u, p));

    K(:,:, i) = lqr(A, B, Q, R);
end


F_closed_loop = @(x, x_star, u_star, K)dynamics(x, u_star-K*(x-x_star), p);

x_cl = zeros(size(soln.grid.state));
x_cl(:,1) = x0;

for i=2:length(t)
    sol = ode45(@(T,x)F_closed_loop(x, xsol(:,i), u(i)+(u(i)-u(i-1))*(T-t(i-1)), 0.5*(K(:,:,i-1)+K(:, :, i))), [t(i-1), t(i)], x_cl(:,i-1));
    x_cl(:,i) = sol.y(:,end);
    
end
figure(3); clf;
subplot(2,1,1)

plot(t, xsol(1,:), t, x_cl(1,:))
subplot(2,1,2)
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
    animate(t(i), xs(:,i), p)



    % set(lblTime,'String', i - 1);
    % set(lblAction,'String', u(i));
    
    drawnow;
    pause(0.05);
end


%% MPC

start = [pi/2; 0; 0; 0];

T=1e-3;

tmax=20; t=0:T:tmax;
NL=length(t)-1;

A = J(start, 1e-5*ones(4,1), @(x)dynamics(x, 0, p));

B = J(0, 1e-5, @(u)dynamics(start, u, p));

ss1 = ss(A, B, eye(4), 0);
ss1d = c2d(ss1, T);


A=ss1d.A; B=ss1d.B; C=ss1d.C; D=ss1d.D;


Ip = eye(4);
Btilda = [C*B; B];
Itilda = zeros(length(C(:,1))+length(A(:,1)),length(Ip(1,:)));
Itilda(1:length(Ip(:,1)),1:length(Ip(1,:))) = Ip*0.9;

Ftilda = [C*A; A];
Qtilda = zeros(length(Ip(:,1))+length(A(:,1)),...
length(Ip(:,1))+length(A(:,1)));
Qtilda(1:4,1:4) = eye(4)*1;

Atilda = [Itilda Ftilda];
R = 1e-6;

[Ktilda,Ldare,Gdare] = dare(Atilda,Btilda,Qtilda,R);

RKBstuff = (R + Btilda'*Ktilda*Btilda)^-1 * Btilda'*Ktilda;
GI = RKBstuff * Itilda;
GX = RKBstuff * Ftilda;

NLbig=5000;

Gd = zeros(1,NLbig*4); 
Gd(1:4) = -GI;

Actilda = Atilda - (Btilda * RKBstuff * Atilda);
Xtilda(1).mat = -(Actilda' * Ktilda * Itilda);
RKBgd = (R + Btilda'*Ktilda*Btilda)^-1 * Btilda';
Xtilda(NLbig).mat = 0*Xtilda(1).mat; % initializing space in var

for n=2:NLbig
Gd(4*(n-1)+1:4*n) = RKBgd * Xtilda(n-1).mat;
Xtilda(n).mat = Actilda' * Xtilda(n-1).mat;
end

figure(1); clf; subplot(211)
plot(t(2:end),-Gd(1:NL),'LineWidth',2); grid on
xlabel('time (sec)')
ylabel('gain')
title('fig 6')
% axis([0 20 0 1500]);


tstep=0:T:7;
% ystep corresponds to x-directed biped motion of com
tfac = input('Enter a value for tfac. (tfac=1 to match Kajita): ');
if isempty(tfac) || isstr(tfac)
tfac = 1
end
ystep=0*tstep;
ystep=ystep+.3*(tstep>2.6*tfac);
ystep=ystep+.3*(tstep>3.4*tfac);
ystep=ystep+.3*(tstep>4.2*tfac); % desired zmp trajectory
ystep2=0*tstep;
ystep2=ystep2+.1*(tstep>1.8*tfac);
ystep2=ystep2-.2*(tstep>2.6*tfac);
ystep2=ystep2+.2*(tstep>3.4*tfac);
ystep2=ystep2-.2*(tstep>4.2*tfac);
ystep2=ystep2+.1*(tstep>5*tfac); % desired zmp trajectory
% smooth ystep and ystep2 somewhat
for n=1:8
ystep(10:end-10)=.5*(ystep(10:end-10)+ystep(11:end-9));
ystep2(10:end-10)=.5*(ystep2(10:end-10)+ystep2(11:end-9));
end
figure(3); clf
figure(2); clf