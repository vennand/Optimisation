% function plotHeatmap(data)

x1labelpositions = 1:20:201;
x1labelnames = {'20', '40', '60', '80', '100'};
x2labelnames = {'w_0', 'w_{opt}'};

ylabelpositions = [1 4 7 10 13 15 18 20 22 24 27 29 31 34 35 37 40 41];
ylabelnames = {'trans pelvis', 'rot pelvis', 'thorax', 'head',...
             'shoulder r', 'arm r', 'elbow r', 'wrist r', ...
             'shoulder l', 'arm l', 'elbow l', 'wrist l', ...
             'thigh r', 'knee r', 'ankle r', ...
             'thigh l', 'knee l', 'ankle l'};

figure(1)
sgtitle('Gravity')
subplot(121)
imagesc(horzcat(abs(data.stateGravityGrad_init_reorganized(1:42,:)), ...
                abs(data.stateGravityGrad_opt_reorganized(1:42,:))))
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
imagesc(horzcat(abs(data.stateGravityGrad_init_reorganized(43:end,:)), ...
                abs(data.stateGravityGrad_opt_reorganized(43:end,:))))
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



figure(2)
sgtitle('Single shooting simple simulation final state')

subplot(121)
imagesc(abs(data.simStateGravityGrad(1:42,:)))
colorbar
ax1 = gca;
xticks(ax1, [1 2 3])
xticklabels(ax1, {'x', 'y', 'z'})
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
for i=1:3
    xline(ax1, i+0.5, '--w', 'LineWidth', 0.05);
end
title('q')

subplot(122)
imagesc(abs(data.simStateGravityGrad(43:end,:)))
colorbar
ax1 = gca;
xticks(ax1, [1 2 3])
xticklabels(ax1, {'x', 'y', 'z'})
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
for i=1:3
    xline(ax1, i+0.5, '--w', 'LineWidth', 0.05);
end
title('qd')



figure(3)
sgtitle('Single shooting simulation final state, with u opt')

subplot(121)
imagesc(abs(data.simStateGravityGrad_MX(1:42,:)))
colorbar
ax1 = gca;
xticks(ax1, [1 2 3])
xticklabels(ax1, {'x', 'y', 'z'})
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
for i=1:3
    xline(ax1, i+0.5, '--w', 'LineWidth', 0.05);
end
title('q')

subplot(122)
imagesc(abs(data.simStateGravityGrad_MX(43:end,:)))
colorbar
ax1 = gca;
xticks(ax1, [1 2 3])
xticklabels(ax1, {'x', 'y', 'z'})
yticks(ax1, ylabelpositions)
yticklabels(ax1, ylabelnames)
for i=1:3
    xline(ax1, i+0.5, '--w', 'LineWidth', 0.05);
end
title('qd')
