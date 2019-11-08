clear; close all; clc;

run('startup.m')
import casadi.*

opti = casadi.Opti();

% T = MX.sym('T', 1); % seconds
T = 10; % seconds

g = 9.81;
l = 1;
m1 = 1;
m2 = 2;

N = 20; % nb colloc nodes

x = opti.variable(4,N);% state (q1,q2,qd1,qd2)
u = opti.variable(1,N);% control (u)

J = 0; % initialization
opti.minimize(J);% minimize objective function 

% opti.set_initial(x(:,1),[0;0;0;0])
opti.set_initial(u,ones(1,N)/100)
opti.subject_to(x(:,1)==[0;0;0;0]);% initial constraint
opti.subject_to(u(1)==0);% initial constraint

opti.subject_to(x(:,end)==[0;pi;0;0]);% terminal constraint on position

opti.subject_to(abs(x(1,:))<5);% constraint on position

% trap colloc

%     q1 = x(1);
q2 = x(2,:);
qd1 = x(3,:);
qd2 = x(4,:);

qdd1 = (l.*m2.*sin(q2).*qd2.^2 + u + m2.*g.*cos(q2).*sin(q2))./(m1 +m2.*(1-cos(q2).^2));
qdd2 = (l.*m2.*cos(q2).*sin(q2).*qd2.^2 + u.*cos(q2) + (m1+m2).*g.*sin(q2))./(l.*m1 +l.*m2.*(1-cos(q2).^2));

fx = [qd1;qd2;qdd1;qdd2];

for n = 1:N-1
    
    J = J + 1/2 * (u(:,n).^2 + u(:,n+1).^2);
    
    opti.subject_to(1/2.*(fx(:,n+1)+fx(:,n)) == x(:,n+1) - x(:,n));
end

opti.solver('ipopt');
sol = opti.solve();
solx = sol.value(x);
solu = sol.value(u);
% plot(solx(1,:))
% figure()
% plot(solx(3,:))
% figure()
% plot(solu)
