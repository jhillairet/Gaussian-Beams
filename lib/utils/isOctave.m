function Bool = isOctave()
%
% Bool = isOctave()
%
% �tant donn� les diff�rences qui peuvent exister entre MATLAB
% et Octave (incompatibilit�s, fonctione nomm�es diff�remment, 
% etc...), si l'on souhaite r�aliser des scripts portables 
% (i.e. qui tournent � la fois sous MATLAB et Octave, 
% ce qui est conseill� !) on peut impl�menter du code 
% conditionnel relatif � chacun de ces environnement 
% en r�alisant, par exemple, le test sur l'existence 
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
%  Bool : bool�en : True si c'est Octave qui lance le script, False sinon. 
% 
if ~ exist('OCTAVE_VERSION')  % MATLAB
    Bool = false;
else                          % Octave
    Bool = true;
end
