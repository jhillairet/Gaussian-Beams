% RELATIVE ERROR ROOT MEAN SQUARE
function RMS = rms(ref, comp)
%
% FONCTION :
% Calcule la moyenne des moindres carr�s de l'erreur relative entre
% une r�f�rence et son compar�
%
% ARGUMENTS :
%  ref : r�f�rence (1 x N)
%  comp: compar�   (1 x N)
%
% RETOURNE :
%  RMS : scalaire 
%
Er_rel = abs(ref - comp) ./ abs(comp);


RMS = sqrt(mean(Er_rel).^2 + std(Er_rel).^2);
