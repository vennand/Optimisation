clear, clc, close all
run('startup.m')
import casadi.*

data.nDoF = 42;

% data.Duration = 1; % Time horizon
data.Nint = 100;% number of control nodes
data.odeMethod = 'rk4';
data.NLPMethod = 'MultipleShooting';

data.dataFile = '../data/Do_822_contact_2.c3d';
data.kalmanDataFile_q = '../data/Do_822_contact_2_MOD200.00_GenderF_DoCig_Q_brut.mat';
data.kalmanDataFile_v = '../data/Do_822_contact_2_MOD200.00_GenderF_DoCig_V_brut.mat';
data.kalmanDataFile_a = '../data/Do_822_contact_2_MOD200.00_GenderF_DoCig_A_brut.mat';

% Spécific à Do_822_contact_2.c3d
% Le saut est entre les frames 3050 et 3385
data.frames = 3050:3386;
data.labels = 1:95;

data.realNint = length(data.frames);

[model, data] = GenerateModel(data);

%%
frames = data.frames;

% Assuming there is only one variable, and that it is q, v, a and tau
kalmanData_q = load(data.kalmanDataFile_q);
kalmanData_v = load(data.kalmanDataFile_v);
kalmanData_a = load(data.kalmanDataFile_a);

kalmanData_q = struct2cell(kalmanData_q);
kalmanData_v = struct2cell(kalmanData_v);
kalmanData_a = struct2cell(kalmanData_a);

q = kalmanData_q{1};
v = kalmanData_v{1};
a = kalmanData_a{1};

q = q(:,frames);
v = v(:,frames);
a = a(:,frames);

[dof, frame_length] = size(q);

tau = zeros(dof,frame_length);
for i = 1:frame_length
    tau(:,i) = ID(model, q(:,i), v(:,i), a(:,i));
end

Nint = data.Nint;
realNint = data.realNint;

new_q = zeros(dof, Nint+1);
new_v = zeros(dof, Nint+1);
new_a = zeros(dof, Nint+1);
new_tau = zeros(dof-6, Nint+1);

for old_value = 1:Nint+1
    new_value = range_conversion(old_value, Nint+1, 1, realNint, 1);
    
    new_q(:,old_value) = q(:,round(new_value));
    new_v(:,old_value) = v(:,round(new_value));
    new_a(:,old_value) = a(:,round(new_value));
    new_tau(:,old_value) = tau(7:end,round(new_value));
end

% Storing data
data.kalman_q = new_q;
data.kalman_v = new_v;
data.kalman_a = new_a;
data.kalman_tau = new_tau;

data.kalman_qFull = q;
data.kalman_vFull = v;
data.kalman_aFull = a;
data.kalman_tauFull = tau;


%%
frames = data.frames;
labels = data.labels;

real_data = ezc3dRead(data.dataFile);
frequency = real_data.header.points.frameRate;
data.Duration = length(frames)/frequency;

markers = real_data.data.points(:,labels,frames)/1000; %meter

Nint = data.Nint;
[num_dimension, num_label, num_frames] = size(markers);

markers_reformat = zeros(num_dimension, num_label, Nint+1);

for old_value = 1:Nint+1
    new_value = range_conversion(old_value, Nint+1, 1, num_frames, 1);
    markers_reformat(:,:,old_value) = markers(:,:,round(new_value));
end

markers_reformat = reshape(markers_reformat, num_dimension*num_label, Nint+1)';
%%


kalman_markers = zeros(Nint+1,285);
kalman_q = data.kalman_q;
for i=1:Nint+1
%     kalman_q(1:3,i) = [0;0;0];
    kalman_markers(i,:) = base_referential_coor(model, kalman_q(:,i));
end


%%

x = [markers_reformat(1,1:3) 1]';
y = [markers_reformat(1,4:6) 1]';
z = [markers_reformat(1,7:9) 1]';
zz = [markers_reformat(1,10:12) 1]';
zzz = [markers_reformat(1,13:15) 1]';
X = [x y z zz];

x_prime = [kalman_markers(1,1:3) 1]';
y_prime = [kalman_markers(1,4:6) 1]';
z_prime = [kalman_markers(1,7:9) 1]';
zz_prime = [kalman_markers(1,10:12) 1]';
zzz_prime = [kalman_markers(1,13:15) 1]';
X_prime = [x_prime y_prime z_prime zz_prime];

%%

opti = casadi.Opti();

RT11 = opti.variable();
RT12 = opti.variable();
RT13 = opti.variable();
RT14 = opti.variable();
RT21 = opti.variable();
RT22 = opti.variable();
RT23 = opti.variable();
RT24 = opti.variable();
RT31 = opti.variable();
RT32 = opti.variable();
RT33 = opti.variable();
RT34 = opti.variable();

RT = [RT11 RT12 RT13 RT14; ...
      RT21 RT22 RT23 RT24; ...
      RT31 RT32 RT33 RT34; ...
      0    0    0    1];
  
R = [RT11 RT12 RT13; ...
     RT21 RT22 RT23; ...
     RT31 RT32 RT33];
%%
opti.minimize(norm_fro(X_prime - RT * X)^2);

opti.subject_to(R' == inv(R));
opti.set_initial(RT11, 1)
opti.set_initial(RT22, 1)
opti.set_initial(RT33, 1)

%%

opti.solver('ipopt');
% sol = opti.solve();

%%

scatter3(kalman_markers(1,1:3:end),kalman_markers(1,2:3:end),kalman_markers(1,3:3:end))
hold on
scatter3(markers_reformat(1,1:3:end),markers_reformat(1,2:3:end),markers_reformat(1,3:3:end))
axis equal