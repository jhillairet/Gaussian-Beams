function Repere = struct_Repere(Centre, ex, ez)
%
% Structure Repere
%
% Repere = struct_Repere(Centre, ex, ez)
%
% Description d'un repere cart�sien quelconque. 
% Ce rep�re est d�fini par les coordonn�es de son centre 
% et de ses trois vecteurs unitaires exprim-A�s dans le rep�re absolu.-b
%
% PARAMETRES :
%  Centre : coordonn-A�es du centre du rep�re            (3 x 1)-b
%  ex     : vecteur unitaire ex du rep-A�re              (3 x 1)-b
%  ez     : vecteur unitaire ez du rep-A�re              (3 x 1)-b
%
% RETOURNE
%  Rep-A�re : structure Rep�re-b

Repere = struct(...
    'Centre', Centre, ...
    'ex', ex, ...
    'ez', ez);
