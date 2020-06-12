% function plotHeatmap(data)
segment_name = {'pelvis', 'thorax', 'right thigh', 'left thigh'};

x1labelpositions = 1:20:201;
x1labelnames = {'20', '40', '60', '80', '100'};
x2labelnames = {'w_0', 'w_{opt}'};

ylabelpositions = [1 4 7 10 13 15 18 20 22 24 27 29 31 34 35 37 40 41];
ylabelnames = {'trans pelvis', 'rot pelvis', 'thorax', 'head',...
             'shoulder r', 'arm r', 'elbow r', 'wrist r', ...
             'shoulder l', 'arm l', 'elbow l', 'wrist l', ...
             'thigh r', 'knee r', 'ankle r', ...
             'thigh l', 'knee l', 'ankle l'};

for l=1:data.nSegment
figure(5*(l-1)+1)
sgtitle(['Mass, segment: ' segment_name{l}])
subplot(121)
imagesc(horzcat(abs(data.stateMassGrad_init_reorganized{l}(1:42,:)), ...
                abs(data.stateMassGrad_opt_reorganized{l}(1:42,:))))
ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position', ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
xticks(ax1, x1labelpositions)
xticklabels(ax1, {'', x1labelnames{:}, x1labelnames{:}})
xlim(ax2, [0 200])
xticks(ax2, [50 150])
xticklabels(ax2, x2labelnames)
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
xline(ax1, 101, 'LineWidth', 2);
title('q')
ax2.YAxis.Visible = 'off';
axes(ax1)
colorbar(ax1)
colorbar(ax2, 'Ticks', [])

subplot(122)
imagesc(horzcat(abs(data.stateMassGrad_init_reorganized{l}(43:end,:)), ...
                abs(data.stateMassGrad_opt_reorganized{l}(43:end,:))))
ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position', ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
xticks(ax1, x1labelpositions)
xticklabels(ax1, {'', x1labelnames{:}, x1labelnames{:}})
xlim(ax2, [0 200])
xticks(ax2, [50 150])
xticklabels(ax2, x2labelnames)
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
xline(ax1, 101, 'LineWidth', 2);
title('qd')
ax2.YAxis.Visible = 'off';
axes(ax1)
colorbar(ax1)
colorbar(ax2, 'Ticks', [])

figure(5*(l-1)+2)
sgtitle(['CoM, segment: ' segment_name{l}])
subplot(121)
imagesc(horzcat(abs(data.stateCoMGrad_init_reorganized{l}(1:42,:)), ...
                abs(data.stateCoMGrad_opt_reorganized{l}(1:42,:))))
ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position', ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
xticks(ax1, 3*x1labelpositions)
xticklabels(ax1, {'', x1labelnames{:}, x1labelnames{:}})
xlim(ax2, [0 600])
% xticks(ax2, [150 450])
% xticklabels(ax2, x2labelnames)
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
xline(ax1, 303.5, 'LineWidth', 2);
for i=1:201
    xline(ax1, i*3+0.5, '--', 'LineWidth', 0.05);
end
title('q')
ax2.YAxis.Visible = 'off';
axes(ax1)
colorbar(ax1)
colorbar(ax2, 'Ticks', [])

subplot(122)
imagesc(horzcat(abs(data.stateCoMGrad_init_reorganized{l}(43:end,:)), ...
                abs(data.stateCoMGrad_opt_reorganized{l}(43:end,:))))
ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position', ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
xticks(ax1, 3*x1labelpositions)
xticklabels(ax1, {'', x1labelnames{:}, x1labelnames{:}})
xlim(ax2, [0 600])
xticks(ax2, [150 450])
xticklabels(ax2, x2labelnames)
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
xline(ax1, 303.5, 'LineWidth', 2);
for i=1:201
    xline(ax1, i*3+0.5, '--', 'LineWidth', 0.05);
end
title('qd')
ax2.YAxis.Visible = 'off';
axes(ax1)
colorbar(ax1)
colorbar(ax2, 'Ticks', [])

figure(5*(l-1)+3)
sgtitle(['Inertia, segment: ' segment_name{l}])
subplot(121)
imagesc(horzcat(abs(data.stateInertiaGrad_init_reorganized{l}(1:42,:)), ...
                abs(data.stateInertiaGrad_opt_reorganized{l}(1:42,:))))
ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position', ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
xticks(ax1, 3*x1labelpositions)
xticklabels(ax1, {'', x1labelnames{:}, x1labelnames{:}})
xlim(ax2, [0 600])
xticks(ax2, [150 450])
xticklabels(ax2, x2labelnames)
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
xline(ax1, 303.5, 'LineWidth', 2);
for i=1:201
    xline(ax1, i*3+0.5, '--', 'LineWidth', 0.05);
end
title('q')
ax2.YAxis.Visible = 'off';
axes(ax1)
colorbar(ax1)
colorbar(ax2, 'Ticks', [])

subplot(122)
imagesc(horzcat(abs(data.stateInertiaGrad_init_reorganized{l}(43:end,:)), ...
                abs(data.stateInertiaGrad_opt_reorganized{l}(43:end,:))))
ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position', ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
xticks(ax1, 3*x1labelpositions)
xticklabels(ax1, {'', x1labelnames{:}, x1labelnames{:}})
xlim(ax2, [0 600])
xticks(ax2, [150 450])
xticklabels(ax2, x2labelnames)
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
xline(ax1, 303.5, 'LineWidth', 2);
for i=1:201
    xline(ax1, i*3+0.5, '--', 'LineWidth', 0.05);
end
title('qd')
ax2.YAxis.Visible = 'off';
axes(ax1)
colorbar(ax1)
colorbar(ax2, 'Ticks', [])


figure(5*(l-1)+4)
sgtitle(['Single shooting simple simulation final state, segment: ' segment_name{l}])

subplot(231)
imagesc(abs(data.simStateMassGrad{l}(1:42,:)))
colorbar
ax1 = gca;
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
ax1.XAxis.Visible = 'off';
title('q / mass')

subplot(234)
imagesc(abs(data.simStateMassGrad{l}(43:end,:)))
colorbar
ax1 = gca;
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
ax1.XAxis.Visible = 'off';
title('qd / mass')

subplot(232)
imagesc(abs(data.simStateCoMGrad{l}(1:42,:)))
colorbar
ax1 = gca;
xticks(ax1, [1 2 3])
xticklabels(ax1, {'x', 'y', 'z'})
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
for i=1:3
    xline(ax1, i+0.5, '--w', 'LineWidth', 0.05);
end
title('q / CoM')

subplot(235)
imagesc(abs(data.simStateCoMGrad{l}(43:end,:)))
colorbar
ax1 = gca;
xticks(ax1, [1 2 3])
xticklabels(ax1, {'x', 'y', 'z'})
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
for i=1:3
    xline(ax1, i+0.5, '--w', 'LineWidth', 0.05);
end
title('qd / CoM')

subplot(233)
imagesc(abs(data.simStateInertiaGrad{l}(1:42,:)))
colorbar
ax1 = gca;
xticks(ax1, [1 2 3])
xticklabels(ax1, {'x', 'y', 'z'})
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
for i=1:3
    xline(ax1, i+0.5, '--w', 'LineWidth', 0.05);
end
title('q / inertia')

subplot(236)
imagesc(abs(data.simStateInertiaGrad{l}(43:end,:)))
colorbar
ax1 = gca;
xticks(ax1, [1 2 3])
xticklabels(ax1, {'x', 'y', 'z'})
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
for i=1:3
    xline(ax1, i+0.5, '--w', 'LineWidth', 0.05);
end
title('qd / inertia')


figure(5*l)
sgtitle(['Single shooting simulation final state, with u_opt, segment: ' segment_name{l}])

subplot(231)
imagesc(abs(data.simStateMassGrad_MX{l}(1:42,:)))
colorbar
ax1 = gca;
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
ax1.XAxis.Visible = 'off';
title('q / mass')

subplot(234)
imagesc(abs(data.simStateMassGrad_MX{l}(43:end,:)))
colorbar
ax1 = gca;
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
ax1.XAxis.Visible = 'off';
title('qd / mass')

subplot(232)
imagesc(abs(data.simStateCoMGrad_MX{l}(1:42,:)))
colorbar
ax1 = gca;
xticks(ax1, [1 2 3])
xticklabels(ax1, {'x', 'y', 'z'})
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
for i=1:3
    xline(ax1, i+0.5, '--w', 'LineWidth', 0.05);
end
title('q / CoM')

subplot(235)
imagesc(abs(data.simStateCoMGrad_MX{l}(43:end,:)))
colorbar
ax1 = gca;
xticks(ax1, [1 2 3])
xticklabels(ax1, {'x', 'y', 'z'})
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
for i=1:3
    xline(ax1, i+0.5, '--w', 'LineWidth', 0.05);
end
title('qd / CoM')

subplot(233)
imagesc(abs(data.simStateInertiaGrad_MX{l}(1:42,:)))
colorbar
ax1 = gca;
xticks(ax1, [1 2 3])
xticklabels(ax1, {'x', 'y', 'z'})
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
for i=1:3
    xline(ax1, i+0.5, '--w', 'LineWidth', 0.05);
end
title('q / inertia')

subplot(236)
imagesc(abs(data.simStateInertiaGrad_MX{l}(43:end,:)))
colorbar
ax1 = gca;
xticks(ax1, [1 2 3])
xticklabels(ax1, {'x', 'y', 'z'})
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
for i=1:3
    xline(ax1, i+0.5, '--w', 'LineWidth', 0.05);
end
title('qd / inertia')
end
