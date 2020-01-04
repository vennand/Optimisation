function [lbw, ubw] = GenerateInitialConstraints(model, lbw, ubw)
% To set final constraint for x, set from end-model.nx to end-model.nu
% To set final constraint for u, set from model.nu to end

% Initial velocities
lbw(end-(model.nx+model.nu):end-model.nu) = zeros(model.nx,1);
ubw(1:model.nx) = zeros(model.nx,1);
end