run('/home/laboratoire/mnt/E/Bureau/Partha/GIT_S2MLib/loadS2MLib_pwd.m');
addpath('/home/laboratoire/mnt/E/Librairies/biorbd/18_juin_2018/release/wrapper/matlab')

h = biorbd('new', '../data/DoCi.s2mMod');
% S2M_rbdl_reader(h)

Q = zeros(42,1);
massMatrix = S2M_rbdl('massMatrix', h, Q);

load('../data/Do_822_contact_2_MOD200.00_GenderF_DoCig_Q.mat')
load('../data/Do_822_contact_2_MOD200.00_GenderF_DoCig_V.mat')
load('../data/Do_822_contact_2_MOD200.00_GenderF_DoCig_A.mat')

Tau_Kalman = S2M_rbdl('inverseDynamics', h, Q2, V2, A2);

for i = 1:frame_length
    Tau_spatial_v2(:,i) = ID(model, Q2(:,i), V2(:,i), A2(:,i));
end

htot = [];
for i = 1:3750
ret = EnerMo( model, Q2(:,i), V2(:,i) );
htot = [htot ret.htot];
end