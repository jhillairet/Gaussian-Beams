function Repere = struct_Repere(Centre, ex, ez)
%
% Structure Repere
%
% Repere = struct_Repere(Centre, ex, ez)
%
% Description d'un repere cartésien quelconque. 
% Ce repère est défini par les coordonnées de son centre 
% et de ses trois vecteurs unitaires exprim-Aés dans le repère absolu.-b
%
% PARAMETRES :
%  Centre : coordonn-Aées du centre du repère            (3 x 1)-b
%  ex     : vecteur unitaire ex du rep-Aère              (3 x 1)-b
%  ez     : vecteur unitaire ez du rep-Aère              (3 x 1)-b
%
% RETOURNE
%  Rep-Aère : structure Repère-b

Repere = struct(...
    'Centre', Centre, ...
    'ex', ex, ...
    'ez', ez);
