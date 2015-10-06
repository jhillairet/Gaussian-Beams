% faisceau incident
function [E, H] = gb2D_create(R, param, bool_TE)
%
% Créé un faisceau gaussien 2D se propageant selon z du type :
%
% E =  e^(-jk/2*x^2*Q(z) -jkz) * ey
%
% PARAMÈTRES :
% R : vecteur position [X;Y;Z]   : (3, N)
% param   : structure contenant W0, k, ZW0
% bool_TE : transverse Electrique (E /A) si true. TM (H /A) si false
%
% RETOURNE :
% [E, H] : champ [Ex(X,Y,Z); Ey(X,Y,Z); Ez(X,Y,Z)] (3, N)
%
X = R(1,:);
% y = S(2,:);
Z = R(3,:); % axe de propagation du faisceau
Z0 = 120*pi;

I = sqrt(-1);
W = W(Z, param);

Qz = 1./ (Z - param.ZW0 + I*param.k/2*param.W0^2);
Q0= 1./ (- param.ZW0 + I*param.k/2*param.W0^2);
  
u = sqrt(Qz./Q0) .* exp(-I*param.k/2 * X.^2 .* Qz -I*param.k.*Z);
  
if (bool_TE)
  %%%%%%%%%%%%%%%%% Champ électrique
  % faisceau gaussien généralisé
  E = [0;1;0] * u;

  H = ...
      I / (param.k * Z0) .* (...
	  [0;0;1] * (-I*param.k.*X.*Qz.*u) + ...
	  [1;0;0] * (u.*(...
	      I*param.k + Qz/2 - I*param.k/2*Qz.^2.*X.^2)));
  
else % TM
  H = [0;1;0] * u;
  
  E = Z0 .* (...
      [0;0;1] * (-X .* Qz .* u) +...
      [1;0;0] * (u .* (...
	  1 + Qz ./ (2*I*param.k) - (X.*Qz).^2 ./2)));
end




% fonction W
function W = W(Z, param)
% Fonction Wx, Wy
%
% Paramètres :
% z : axe de propagation du faisceau  (1 x N) 
% param : structure contenant k et W0 
% 
% retourne :
% W (1 x N)
%
% NB : On suppose que le repere local utilise pour calculer ce
% faisceau est positionné au centre du waist, en zW0. Aussi, zW0
% n'apparait pas dans les calculs.

W = param.W0 * sqrt(1 + (Z ./ (param.k/2*param.W0.^2)).^2);