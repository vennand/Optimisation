% function  [q,qd,G,gs] = gamma_q( model, qo, qdo )
% 
%   y = formula for calculating y from qo;
%   q = formula for gamma(y);
% 
%   G = Jacobian of gamma;
% 
%   yd = formula for calculating yd from qo (or y) and qdo;
%   qd = G * yd;  (or a formula equivalent to this expression)
% 
%   g = formula for dG/dt * yd;
% 
%   Tstab = some suitable value, such as 0.1 for humans, and 0.01 for machines faster than humans;
%   gstab = 2/Tstab * (qd - qdo) + 1/Tstab^2 * (q - qo);
% 
%   gs = g + gstab;
% end

function  [q,qd,G,gs] = gamma_q( model, qo, qdo )
  
  y = qo(1:6);
  q = [ y; zeros(length(qo)-6,1) ];
  
  % Matrice 42x6
  G = [ diag(ones(6,1)); ...
        zeros(length(qo)-6,6) ];

  yd = qdo(1:6);
  qd = [ yd; zeros(length(qo)-6,1) ];

  g = [ zeros(6,1); zeros(length(qo)-6,1) ];
  
%% Exemple de gamma_q qui n'a aucun impact
%   y = qo;
%   q = y;
% 
%   G = diag(ones(length(qo),1));
% 
%   yd = qdo;
%   qd = yd;
% 
%   g = zeros(length(qo),1);
%%

  Tstab = 0.1;
  gstab = 2/Tstab * (qd - qdo) + 1/Tstab^2 * (q - qo);

  gs = g + gstab;
end