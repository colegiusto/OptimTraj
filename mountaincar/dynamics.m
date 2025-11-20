function dx = dynamics(x,u,p)
% dx = dynamics(x,u,p)
%
% Computes the dynamics for the simple pendulum
%

alpha = p.mtn(x(1,:));



q = x(1,:);
dq = x(2,:);

g = p.g; 
m = p.m;

ddq = -sin(alpha).*cos(alpha)*g+ (cos(alpha)/m).*u;

dx = [dq;ddq];

end