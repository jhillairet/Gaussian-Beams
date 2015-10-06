function F = stackedfield2field(SF, Nx, Ny, Nz)
%
% F = stackedfield2field(SF, Nx, Ny, Nz);
%
% FONCTION :
% transforme une composante d'un champ "stacké" SF calculé à partir d'un stack S
% vers sa représentation 3D type F(X, Y, Z) en vue d'une
% représentation graphique (coupe).
%
% PARAMÈTRES :
% SF : staked field (1, Nx*Ny*Nz)
% Nx, Ny, Nz : taille des dimensions X, Y, Z
%
% RETOURNE :
% F : champ 3D (Nx, Ny, Nz) 
%
% EXEMPLE :
%
% VOIR AUSSI :
% stack : [X, Y, Z] -> S
% unstack : S -> [X, Y, Z]

% verification de la taille des dimensions et de la taille des
% stacks :
if not(isequal(size(SF, 2), Nx*Ny*Nz))
  error('Les tailles des stacks ou des dimensions ne correspondent pas');
end

% recomposition de la matrice 3D du champ. Chacun des points
% (x,y,z) de la matrice correspond au champ calculé en ce point.
% 
F = permute(reshape(SF, Nz, Ny, Nx), [3 2 1]);
