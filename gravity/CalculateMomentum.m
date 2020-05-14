function [data] = CalculateMomentum(model, data)

T = data.Duration; % secondes
Nint = data.Nint; % nb colloc nodes
dN = T/Nint;

% clear, clc, close all
% close all
% run('../startup.m')

% run('/home/laboratoire/mnt/E/Bureau/Partha/GIT_S2MLib/loadS2MLib_pwd.m');
% addpath('/home/laboratoire/mnt/E/Librairies/biorbd/18_juin_2018/release/wrapper/matlab')
% 
% h = biorbd('new', '../data/DoCi.s2mMod');
% % S2M_rbdl_reader(h)
% 
% Q = zeros(42,1);
% massMatrix = S2M_rbdl('massMatrix', h, Q);

% load('Solutions/Do_822_F3100-3311_U1e-07_N105_IPOPTMA57.mat')

load('../data/Do_822_contact_2_MOD200.00_GenderF_DoCig_Q.mat')
load('../data/Do_822_contact_2_MOD200.00_GenderF_DoCig_V.mat')
load('../data/Do_822_contact_2_MOD200.00_GenderF_DoCig_A.mat')

Q2(15:16, :) = Q2(16:-1:15, :);
Q2(24:25, :) = Q2(25:-1:24, :);
V2(15:16, :) = V2(16:-1:15, :);
V2(24:25, :) = V2(25:-1:24, :);
% new_a(15:16, :) = new_a(16:-1:15, :);
% new_a(24:25, :) = new_a(25:-1:24, :);
% new_tau(15:16, :) = new_tau(16:-1:15, :);
% new_tau(24:25, :) = new_tau(25:-1:24, :);

% Tau_Kalman = S2M_rbdl('inverseDynamics', h, Q2, V2, A2);
% 
% for i = 1:frame_length
%     Tau_spatial_v2(:,i) = ID(model, Q2(:,i), V2(:,i), A2(:,i));
% end

% ret.htot is a vector 6x1, with ret.htot(1:3) = angular momentum, and
% ret.htot(4:6) = linear momentum

htot_full = [];
htot_kalman_full = [];
for i = 1:3750
    ret = EnerMo( model, Q2(:,i), V2(:,i) );
    htot_i = ret.htot;
    % htot_i(1:3) = ret.htot(1:3) - ret.mass * cross(ret.cm, ret.vcm);
    htot_i(1:3) = ret.htot(1:3) - cross(ret.cm, ret.htot(4:6));
    htot_full = [htot_full ret.htot];
    htot_kalman_full = [htot_kalman_full htot_i];
end

htot = [];
htot_kalman = [];
for i = 1:data.Nint+1
    ret = EnerMo( model, data.kalman_q(:,i), data.kalman_v(:,i) );
    htot_i = ret.htot;
    % htot_i(1:3) = ret.htot(1:3) - ret.mass * cross(ret.cm, ret.vcm);
    htot_i(1:3) = ret.htot(1:3) - cross(ret.cm, ret.htot(4:6));
    % htot = [htot ret.htot];
    htot_kalman = [htot_kalman htot_i];
end

htot_estim = [];
for i = 1:data.Nint+1
    ret = EnerMo( model, data.q_opt(:,i), data.v_opt(:,i) );
    htot_i = ret.htot;
    % htot_i(1:3) = ret.htot(1:3) - ret.mass * cross(ret.cm, ret.vcm);
    htot_i(1:3) = ret.htot(1:3) - cross(ret.cm, ret.htot(4:6));
    % htot = [htot ret.htot];
    htot_estim = [htot_estim htot_i];
end

htot_kalman_slope = (htot_kalman(:,2:end) - htot_kalman(:,1:end-1))/ret.mass/dN;
htot_estim_slope = (htot_estim(:,2:end) - htot_estim(:,1:end-1))/ret.mass/dN;

% Linear regression
% lm_x = fitlm(1:data.Nint+1,htot_kalman(4,:),'linear');
% lm_y = fitlm(1:data.Nint+1,htot_kalman(5,:),'linear');
% lm_z = fitlm(1:data.Nint+1,htot_kalman(6,:),'linear');
polylm_x = polyfit(1:data.Nint+1,htot_kalman(4,:),1);
polylm_y = polyfit(1:data.Nint+1,htot_kalman(5,:),1);
polylm_z = polyfit(1:data.Nint+1,htot_kalman(6,:),1);

% data.gravityLinearRegression = [lm_x.Coefficients.Estimate(2); lm_y.Coefficients.Estimate(2); lm_z.Coefficients.Estimate(2)];
data.gravityLinearRegression = [polylm_x(1); polylm_y(1); polylm_z(1)];

colors = [[1, 0, 0]; [0, 0.5, 0]; [0, 0, 1]];
set(groot,'defaultAxesColorOrder', colors);

% figure()
% plot(htot_full(1:3,:)')
% figure()
% plot(htot_kalman_full(1:3,:)')
% figure()
% plot(htot(1:3,:)')
figure()
hold on
plot(htot_kalman(1:3,:)','.-')
plot(htot_estim(1:3,:)','-')
figure()
hold on
plot(htot_estim(4:6,:)')
plot(htot_kalman(4:6,:)','.-')

figure()
hold on
plot(htot_estim_slope(1:3,:)')
plot(htot_kalman_slope(1:3,:)','.-')
figure()
hold on
plot(htot_estim_slope(4:6,:)')
plot(htot_kalman_slope(4:6,:)','.-')

% yline(mean(htot_kalman_slope(4,:)),'-.','Color',colors(1,:));
% yline(mean(htot_kalman_slope(5,:)),'-.','Color',colors(2,:));
% yline(mean(htot_kalman_slope(6,:)),'-.','Color',colors(3,:));
% 
% yline(polylm_x(1),'--','Color',colors(1,:));
% yline(polylm_y(1),'--','Color',colors(2,:));
% yline(polylm_z(1),'--','Color',colors(3,:));

end