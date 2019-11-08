clear; close all; clc;

run('Simulation.m')

run('startup.m')
import casadi.*

opti = casadi.Opti();

% T = MX.sym('T', 1); % seconds
T = 1; % seconds
N = 20*T; % nb colloc nodes

dN = T/N;

NB = 1;

temp_markers = size(model.markers.coordinates);
N_cardinal_coor = temp_markers(1);
N_markers = temp_markers(2);

coor = opti.variable(N, N_cardinal_coor);

J = 0; % initialization

% opti.set_initial(q,zeros(2*model.NB,N))
% opti.subject_to(q(1,:)==[0.25 0.05 -0.012]);% initial constraint
% Estimated_PosMarkers = cell(N, N_cardinal_coor, N_markers);

for n = 1:N
%     Estimated_PosMarkers(n,:,:) = base_referential_coor(model, q(:,n));
    
    for m = 1:1 %N_markers
        distance_between_points = 0;
        for l = 1:N_cardinal_coor
            distance_between_points = distance_between_points + (PosMarkers{ceil(n/N*nb_dt),l,m} - coor(n,l)).^2;
        end
        distance_between_points = distance_between_points;
        J = J + 1/2 * distance_between_points;
    end
end

opti.minimize(J);% minimize objective function

opti.solver('ipopt');
sol = opti.solve();
solx = sol.value(coor);
t_opti = linspace(0.05,T,N);
t_simu = linspace(0,T,1000);

hold on
plot(t_opti,solx,'x');
plot(t_simu,cell2mat(PosMarkers(:,:,1)),'x');
