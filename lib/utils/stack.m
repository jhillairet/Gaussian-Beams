function S = stack(X, Y, Z)
% STACK
%
% S = stack(X, Y, Z)
%
% empile 3 vecteurs de tailles diff�rentes vers un seul vecteur
% du type P(n) = [x(n); y(n); z(n)] avec n=[1:Nx*Ny*Nz]
%
% PARAM�TRES :
% X : (1, Nx)
% Y : (1, Ny)
% Z : (1, Nz)
%
% RETOURNE :
% S : tableau (3, Nx*Ny*Nz)
%
% EXEMPLE :
%  >> x = [1:4]; y=[1:3]; z=[1:2]; S=stack(x,y,z)
%  
%  S =
%  
%    Columns 1 through 12 
%  
%       1     1     1     1     1     1     2     2     2     2     2     2
%       1     1     2     2     3     3     1     1     2     2     3     3
%       1     2     1     2     1     2     1     2     1     2     1     2
%  
%    Columns 13 through 24 
%  
%       3     3     3     3     3     3     4     4     4     4     4     4
%       1     1     2     2     3     3     1     1     2     2     3     3
%       1     2     1     2     1     2     1     2     1     2     1     2
% 
% >> isequal(length(x)*length(y)*length(z), length(S))
% 
% ans =
% 
%      1
%
% VOIR AUSSI :
%  destack : l'op�ration inverse.
%  stackedfield2field : transformation d'un champ calcul� � partir
%  d'un stack (1D) vers un champs 3D

% protection contre les tableaux vides : 
if (isempty(X))
  X = 0;
end
if (isempty(Y))
  Y = 0;
end
if (isempty(Z))
  Z = 0;
end

% tailles des vecteurs d'entr�e
Nx = length(X);
Ny = length(Y);
Nz = length(Z);

% creation des lignes du vecteur final. Cf l'exemple ci-dessus pour
% comprendre le principe du vecteur final. Apr�s, c'est de l'astuce
% :)
% tailles : (1, Nx*Ny*Nz)
XX = reshape(ones(Ny*Nz, 1) * X, 1, Ny*Nz*Nx);
YY = repmat(reshape(ones(Nz, 1) * Y, 1, Ny*Nz), 1, Nx);
ZZ = repmat(Z, 1, Nx*Ny); 

% creation du vecteur final 
% taille : (3, Nx*Ny*Nz)
S = [XX; YY; ZZ];
