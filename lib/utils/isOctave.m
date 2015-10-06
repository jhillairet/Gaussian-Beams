function Bool = isOctave()
%
% Bool = isOctave()
%
% Étant donné les différences qui peuvent exister entre MATLAB
% et Octave (incompatibilités, fonctione nommées différemment, 
% etc...), si l'on souhaite réaliser des scripts portables 
% (i.e. qui tournent à la fois sous MATLAB et Octave, 
% ce qui est conseillé !) on peut implémenter du code 
% conditionnel relatif à chacun de ces environnement 
% en réalisant, par exemple, le test sur l'existence 
% d'une variable ou constante builtin n'existant que 
% dans l'un ou l'autre de ces 2 environnements.
% 
% Ex: on test ici l'existence de la constante 
% OCTAVE_VERSION (n'existant que sous Octave)
%
% source : http://enacit1.epfl.ch/cours_matlab/mfiles.html
%
% AUCUN PARAMETRE
%
% RETOURNE :
%  Bool : booléen : True si c'est Octave qui lance le script, False sinon. 
% 
if ~ exist('OCTAVE_VERSION')  % MATLAB
    Bool = false;
else                          % Octave
    Bool = true;
end
