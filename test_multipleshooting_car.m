clear; close all; clc;

run('startup.m')
import casadi

opti = casadi.Opti();

% T = MX.sym('T', 1); % seconds
T = 10; % seconds

N = 20; % nb shooting nodes
tn = linspace(0,T,N);% time at nodes
x = opti.variable(3,N);% states (x,y,theta)
u = opti.variable(2,N);% control (a, phi)

Obj = norm(u,'fro');
opti.minimize(Obj);% minimize objective function 

opti.set_initial(u,ones(2,N)/100)
opti.subject_to(x(:,1)==[0;0;0]);% initial constraint
opti.subject_to(u(:,1)==[0;0]);% initial constraint

opti.subject_to(x(1,:)>0);% constraint on position
opti.subject_to(x(2,:)>0);% constraint on position
opti.subject_to(x(:,end)==[1;1;pi]);% terminal constraint on position

opti.subject_to(u(2,:)>0);% constraint on control
% opti.subject_to(abs(u(2,:))<pi/8);% constraint on control

% multiple-shooting (runge-kutta)
for n = 1:N-1
    xe = rg4(x(:,n),u(:,n),N,T);
    opti.subject_to(x(:,n+1)==xe);
end

opti.solver('ipopt');
sol = opti.solve();
solx = sol.value(x);
solu = sol.value(u);
% plot(solx(1,:),solx(2,:))
% axis('equal')
% figure()
% plot(solx(3,:))
% figure()
% plot(solu')

% system's dynamics
function dx = dyn(x, u)
    throttle = u(1);
    steering = u(2); 
%     posx = x(1);% horizontal position
%     posy = x(2);% vertical position
    posth = x(3);% angle
    dx = [throttle*cos(posth),throttle*sin(posth),tan(steering)]';
end

function X_e = rg4(X0,U,N,T)
   M = 4; % RK4 steps per interval
   DT = T/(N-1)/M;
   X_e = X0;
   for j=1:M
       k1 = dyn(X_e, U);
       k2 = dyn(X_e + DT/M/2 * k1, U);
       k3 = dyn(X_e + DT/M/2 * k2, U);
       k4 = dyn(X_e + DT/M * k3, U);
       X_e=X_e+DT/6*(k1 +2*k2 +2*k3 +k4);
   end
end
