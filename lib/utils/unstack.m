function [X, Y, Z] = unstack(S, Nx, Ny, Nz)
%
% [X, Y, Z] = unstack(S, Nx, Ny, Nz)
%
% dépile le tableau S en 3 tableaux originel X, Y, Z de taille
% respec. Nx, Ny, Nz
%
% paramètres :
% S : vecteur (1, N) où N = Nx * Ny * Nz
% Nx, Ny, Nz : taille des dimensions x,y,z
%
% retourne :
% X, Y, Z : vecteurs de tailles respec Nx, Ny, Nz
%
% VOIR AUSSI :
%  stack : l'opération inverse.
%  stackedfield2field : transformation d'un champ calculé à partir
%  d'un stack (1D) vers un champs 3D

% première vérification : N = Nx * Ny * Nz
if not(isequal(length(S), Nx*Ny*Nz))
    error('la taille du tableau S ne correspond pas avec la taille des dimensions x,y,z !');
end

X = downsample(S(1,:), Nz*Ny);    
Y = downsample(S(2,1:Ny*Nz), Nz);  
Z = S(3, 1:Nz);
