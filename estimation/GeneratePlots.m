function [] = GeneratePlots(model, data, q_opt, v_opt, u_opt)

t_estim_x = linspace(0,data.Duration,data.Nint+1);
t_kalman_x = linspace(0,data.Duration,data.Nint+1);

t_estim_u = linspace(0,data.Duration,data.Nint);
t_kalman_tau = linspace(0,data.Duration,data.Nint);

%     42 degrés de liberté, 17 segments
%     Segment			Num	Parent	X	Y	Z
%     Bassin			S1	S0	q4	q5	q6
%     Thorax			S2	S1	−q7	q8	q9
%     Tête			S3	S2	−q10	q11	q12
%     Épaule droite		S4	S2	0	−q13	−q14
%     Bras droit			S5	S4	q15	−q16	−q17
%     Avant-bras droit		S6	S5	q18	0	−q19
%     Main droite		S7	S6	q20	−q21	0
%     Épaule gauche		S8	S2	0	q22	q23
%     Bras gauche		S9	S8	q24	q25	q26
%     Avant-bras gauche		S10	S9	q27	0	q28
%     Main gauche		S11	S10	q29	q30	0
%     Cuisse droite		S12	S1	q31	−q32	−q33
%     Jambe droite		S13	S12	−q34	0	0
%     Pied droit			S14	S13	q35	0	−q36
%     Cuisse gauche		S15	S1	q37	q38	q39
%     Jambe gauche		S16	S15	−q40	0	0
%     Pied gauche		S17	S16	q41	0	q42

base_dof = 6;

pelvis_translation = 1:3;
pelvis_rotation = 4:6;
thorax = 7:9;
head = 10:12;
right_shoulder = 13:14;
right_arm = 15:17;
right_forarm = 18:19;
right_hand = 20:21;
left_shoulder = 22:23;
left_arm = 24:26;
left_forarm = 27:28;
left_hand = 29:30;
right_thigh = 31:33;
right_leg = 34;
right_foot = 35:36;
left_thigh = 37:39;
left_leg = 40;
left_foot = 41:42;

% POSITION PLOT
figure()
subplot(211)
hold on
plot(t_estim_x,q_opt(pelvis_translation,:),'o');
plot(t_kalman_x,data.kalman_q(pelvis_translation,:),'x');
hold off
title('Root translation positions')

subplot(212)
hold on
plot(t_estim_x,q_opt(pelvis_rotation,:),'o');
plot(t_kalman_x,data.kalman_q(pelvis_rotation,:),'x');
hold off
title('Root rotation positions')

figure()
subplot(211)
hold on
plot(t_estim_x,q_opt(thorax,:),'o');
plot(t_kalman_x,data.kalman_q(thorax,:),'x');
hold off
title('Thorax positions')

subplot(212)
hold on
plot(t_estim_x,q_opt(head,:),'o');
plot(t_kalman_x,data.kalman_q(head,:),'x');
hold off
title('Head positions')

figure()
subplot(221)
hold on
plot(t_estim_x,q_opt(right_shoulder,:),'o');
plot(t_kalman_x,data.kalman_q(right_shoulder,:),'x');
hold off
title('Right shoulder positions')

subplot(222)
hold on
plot(t_estim_x,q_opt(right_arm,:),'o');
plot(t_kalman_x,data.kalman_q(right_arm,:),'x');
hold off
title('Right arm positions')

subplot(223)
hold on
plot(t_estim_x,q_opt(right_forarm,:),'o');
plot(t_kalman_x,data.kalman_q(right_forarm,:),'x');
hold off
title('Right forarm positions')

subplot(224)
hold on
plot(t_estim_x,q_opt(right_hand,:),'o');
plot(t_kalman_x,data.kalman_q(right_hand,:),'x');
hold off
title('Right hand positions')

figure()
subplot(221)
hold on
plot(t_estim_x,q_opt(left_shoulder,:),'o');
plot(t_kalman_x,data.kalman_q(left_shoulder,:),'x');
hold off
title('Left shoulder positions')

subplot(222)
hold on
plot(t_estim_x,q_opt(left_arm,:),'o');
plot(t_kalman_x,data.kalman_q(left_arm,:),'x');
hold off
title('Left arm positions')

subplot(223)
hold on
plot(t_estim_x,q_opt(left_forarm,:),'o');
plot(t_kalman_x,data.kalman_q(left_forarm,:),'x');
hold off
title('Left forarm positions')

subplot(224)
hold on
plot(t_estim_x,q_opt(left_hand,:),'o');
plot(t_kalman_x,data.kalman_q(left_hand,:),'x');
hold off
title('Left hand positions')

figure()
subplot(221)
hold on
plot(t_estim_x,q_opt(right_thigh,:),'o');
plot(t_kalman_x,data.kalman_q(right_thigh,:),'x');
hold off
title('Right thigh positions')

subplot(222)
hold on
plot(t_estim_x,q_opt(right_leg,:),'o');
plot(t_kalman_x,data.kalman_q(right_leg,:),'x');
hold off
title('Right leg positions')

subplot(223)
hold on
plot(t_estim_x,q_opt(right_foot,:),'o');
plot(t_kalman_x,data.kalman_q(right_foot,:),'x');
hold off
title('Right foot positions')

figure()
subplot(221)
hold on
plot(t_estim_x,q_opt(left_thigh,:),'o');
plot(t_kalman_x,data.kalman_q(left_thigh,:),'x');
hold off
title('Left thigh positions')

subplot(222)
hold on
plot(t_estim_x,q_opt(left_leg,:),'o');
plot(t_kalman_x,data.kalman_q(left_leg,:),'x');
hold off
title('Left leg positions')

subplot(223)
hold on
plot(t_estim_x,q_opt(left_foot,:),'o');
plot(t_kalman_x,data.kalman_q(left_foot,:),'x');
hold off
title('Left foot positions')

% VELOCITY PLOT

figure()
subplot(211)
hold on
plot(t_estim_x,v_opt(pelvis_translation,:),'o');
plot(t_kalman_x,data.kalman_v(pelvis_translation,:),'x');
hold off
title('Root translation velocity')

subplot(212)
hold on
plot(t_estim_x,v_opt(pelvis_rotation,:),'o');
plot(t_kalman_x,data.kalman_v(pelvis_rotation,:),'x');
hold off
title('Root rotation velocity')

figure()
subplot(211)
hold on
plot(t_estim_x,v_opt(thorax,:),'o');
plot(t_kalman_x,data.kalman_v(thorax,:),'x');
hold off
title('Thorax velocity')

subplot(212)
hold on
plot(t_estim_x,v_opt(head,:),'o');
plot(t_kalman_x,data.kalman_v(head,:),'x');
hold off
title('Head velocity')

figure()
subplot(221)
hold on
plot(t_estim_x,v_opt(right_shoulder,:),'o');
plot(t_kalman_x,data.kalman_v(right_shoulder,:),'x');
hold off
title('Right shoulder velocity')

subplot(222)
hold on
plot(t_estim_x,v_opt(right_arm,:),'o');
plot(t_kalman_x,data.kalman_v(right_arm,:),'x');
hold off
title('Right arm velocity')

subplot(223)
hold on
plot(t_estim_x,v_opt(right_forarm,:),'o');
plot(t_kalman_x,data.kalman_v(right_forarm,:),'x');
hold off
title('Right forarm velocity')

subplot(224)
hold on
plot(t_estim_x,v_opt(right_hand,:),'o');
plot(t_kalman_x,data.kalman_v(right_hand,:),'x');
hold off
title('Right hand velocity')

figure()
subplot(221)
hold on
plot(t_estim_x,v_opt(left_shoulder,:),'o');
plot(t_kalman_x,data.kalman_v(left_shoulder,:),'x');
hold off
title('Left shoulder velocity')

subplot(222)
hold on
plot(t_estim_x,v_opt(left_arm,:),'o');
plot(t_kalman_x,data.kalman_v(left_arm,:),'x');
hold off
title('Left arm velocity')

subplot(223)
hold on
plot(t_estim_x,v_opt(left_forarm,:),'o');
plot(t_kalman_x,data.kalman_v(left_forarm,:),'x');
hold off
title('Left forarm velocity')

subplot(224)
hold on
plot(t_estim_x,v_opt(left_hand,:),'o');
plot(t_kalman_x,data.kalman_v(left_hand,:),'x');
hold off
title('Left hand velocity')

figure()
subplot(221)
hold on
plot(t_estim_x,v_opt(right_thigh,:),'o');
plot(t_kalman_x,data.kalman_v(right_thigh,:),'x');
hold off
title('Right thigh velocity')

subplot(222)
hold on
plot(t_estim_x,v_opt(right_leg,:),'o');
plot(t_kalman_x,data.kalman_v(right_leg,:),'x');
hold off
title('Right leg velocity')

subplot(223)
hold on
plot(t_estim_x,v_opt(right_foot,:),'o');
plot(t_kalman_x,data.kalman_v(right_foot,:),'x');
hold off
title('Right foot velocity')

figure()
subplot(221)
hold on
plot(t_estim_x,v_opt(left_thigh,:),'o');
plot(t_kalman_x,data.kalman_v(left_thigh,:),'x');
hold off
title('Left thigh velocity')

subplot(222)
hold on
plot(t_estim_x,v_opt(left_leg,:),'o');
plot(t_kalman_x,data.kalman_v(left_leg,:),'x');
hold off
title('Left leg velocity')

subplot(223)
hold on
plot(t_estim_x,v_opt(left_foot,:),'o');
plot(t_kalman_x,data.kalman_v(left_foot,:),'x');
hold off
title('Left foot velocity')

% CONTROL PLOT

figure()
subplot(211)
hold on
plot(t_estim_u,u_opt(thorax - base_dof,:),'o');
plot(t_kalman_tau,data.kalman_tau(thorax - base_dof,:),'x');
hold off
title('Thorax control')

subplot(212)
hold on
plot(t_estim_u,u_opt(head - base_dof,:),'o');
plot(t_kalman_tau,data.kalman_tau(head - base_dof,:),'x');
hold off
title('Head control')

figure()
subplot(221)
hold on
plot(t_estim_u,u_opt(right_shoulder - base_dof,:),'o');
plot(t_kalman_tau,data.kalman_tau(right_shoulder - base_dof,:),'x');
hold off
title('Right shoulder control')

subplot(222)
hold on
plot(t_estim_u,u_opt(right_arm - base_dof,:),'o');
plot(t_kalman_tau,data.kalman_tau(right_arm - base_dof,:),'x');
hold off
title('Right arm control')

subplot(223)
hold on
plot(t_estim_u,u_opt(right_forarm - base_dof,:),'o');
plot(t_kalman_tau,data.kalman_tau(right_forarm - base_dof,:),'x');
hold off
title('Right forarm control')

subplot(224)
hold on
plot(t_estim_u,u_opt(right_hand - base_dof,:),'o');
plot(t_kalman_tau,data.kalman_tau(right_hand - base_dof,:),'x');
hold off
title('Right hand control')

figure()
subplot(221)
hold on
plot(t_estim_u,u_opt(left_shoulder - base_dof,:),'o');
plot(t_kalman_tau,data.kalman_tau(left_shoulder - base_dof,:),'x');
hold off
title('Left shoulder control')

subplot(222)
hold on
plot(t_estim_u,u_opt(left_arm - base_dof,:),'o');
plot(t_kalman_tau,data.kalman_tau(left_arm - base_dof,:),'x');
hold off
title('Left arm control')

subplot(223)
hold on
plot(t_estim_u,u_opt(left_forarm - base_dof,:),'o');
plot(t_kalman_tau,data.kalman_tau(left_forarm - base_dof,:),'x');
hold off
title('Left forarm control')

subplot(224)
hold on
plot(t_estim_u,u_opt(left_hand - base_dof,:),'o');
plot(t_kalman_tau,data.kalman_tau(left_hand - base_dof,:),'x');
hold off
title('Left hand control')

figure()
subplot(221)
hold on
plot(t_estim_u,u_opt(right_thigh - base_dof,:),'o');
plot(t_kalman_tau,data.kalman_tau(right_thigh - base_dof,:),'x');
hold off
title('Right thigh control')

subplot(222)
hold on
plot(t_estim_u,u_opt(right_leg - base_dof,:),'o');
plot(t_kalman_tau,data.kalman_tau(right_leg - base_dof,:),'x');
hold off
title('Right leg control')

subplot(223)
hold on
plot(t_estim_u,u_opt(right_foot - base_dof,:),'o');
plot(t_kalman_tau,data.kalman_tau(right_foot - base_dof,:),'x');
hold off
title('Right foot control')

figure()
subplot(221)
hold on
plot(t_estim_u,u_opt(left_thigh - base_dof,:),'o');
plot(t_kalman_tau,data.kalman_tau(left_thigh - base_dof,:),'x');
hold off
title('Left thigh control')

subplot(222)
hold on
plot(t_estim_u,u_opt(left_leg - base_dof,:),'o');
plot(t_kalman_tau,data.kalman_tau(left_leg - base_dof,:),'x');
hold off
title('Left leg control')

subplot(223)
hold on
plot(t_estim_u,u_opt(left_foot - base_dof,:),'o');
plot(t_kalman_tau,data.kalman_tau(left_foot - base_dof,:),'x');
hold off
title('Left foot control')

end
