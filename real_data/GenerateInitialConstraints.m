function [lbw, ubw] = GenerateInitialConstraints(model, data, lbw, ubw)
% To set initial constraint for x, set from 1 to model.nx
% To set initial constraint for u, set from model.x+1 to model.nx+model.nu

x = data.x;

% Initial velocities
% lbw(1:model.nx) = zeros(model.nx,1);
% ubw(1:model.nx) = zeros(model.nx,1);
lbw(1:model.nx) = x(1:model.nx,1)';
ubw(1:model.nx) = x(1:model.nx,1)';
end