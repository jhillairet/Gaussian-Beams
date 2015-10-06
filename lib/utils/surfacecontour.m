function surfacecontour(X,Y,Z, varargin)
% SURFACECONTOUR(X,Y,Z, varargin)
%
% FONCTION :
% Produit la fonction SURFACE avec la superposition de la
% fonction CONTOUR avec les memes parametres
%
% PARAMÈTRES :
% X, Y : abscisse et ordonnées
% Z : valeur aux points (X, Y)
% varargin : options supplémentaire pour la fonction "surface"
%
hold on;

surface(X,Y,Z), shading interp;
if (length(varargin) == 0)
    [C,H] = contour(X,Y,Z);
else
    [C,H] = contour(X,Y,Z, varargin);
end