function dx = dynamics(x,u,p)
% dx = dynamics(x,u,p)
%
% Computes the dynamics for the pendubot
%


q1 = x(1,:);
q2 = x(2, :);
dq1 = x(3,:);
dq2 = x(4, :);



g = p.g; 
m1 = p.m1;
m2 = p.m2;
I1 = p.I1;
I2 = p.I2;

l1 = p.l1;
l2 = p.l2;

lc1 = p.lc1;
lc2 = p.lc2;

d11 = m1*lc1^2 + m2*(l1^2+lc2^2+2*lc1*lc2.*cos(q2)) + I1+ I2;
d22 = m2*lc2^2 + I2;

d12 = m2*(lc2^2+l1*lc2.*cos(q2)) + I2;
d21 = d12;

h1 = -m2*l1*lc2.*sin(q2).*dq2.^2 - 2*m2*l1*lc2.*sin(q2).*dq2.*dq1;
h2 = m2*l1*lc2.*sin(q2).*dq1.^2;

phi1 = (m1*lc1+m2*l1)*g.*cos(q1)+m2*lc2*g.*cos(q1+q2);
phi2 = m2*lc2*g.*cos(q1+q2);


ddq1 = (d12.*(h2+phi2)-d22.*(h1+phi1-u))./(d11.*d22-d12.*d21);
ddq2 = (d21.*(h1+phi1-u)-d11.*(h2+phi2))./(d11.*d22-d12.*d21);

dx = [dq1; dq2; ddq1; ddq2];

end