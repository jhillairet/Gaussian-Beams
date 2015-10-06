% RELATIVE ERROR ROOT MEAN SQUARE
function RMS = rms(ref, comp)
%
% FONCTION :
% Calcule la moyenne des moindres carrés de l'erreur relative entre
% une référence et son comparé
%
% ARGUMENTS :
%  ref : référence (1 x N)
%  comp: comparé   (1 x N)
%
% RETOURNE :
%  RMS : scalaire 
%
Er_rel = abs(ref - comp) ./ abs(comp);


RMS = sqrt(mean(Er_rel).^2 + std(Er_rel).^2);
