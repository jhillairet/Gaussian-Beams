function Champ = struct_champ(Centre, ex, ez, Coefficients, ...
			      Courbure, Sens)
% Champ = struct_champ(Centre, ex, ez, Coefficients, Courbure)    
%
% Cr-A�ation d'une structure Champ. -b
%
% Le champ est exprim-A� comme une somme de N FG.-b
%
% PARAM-A�TRES-b
%  Centre : Coordonn-A�es des centres des FG        (3 x N)-b
%  ex : vecteur ex des rep-A�res des FG             (3 x N)-b
%  ez : vecteur ez des rep-A�res des FG             (3 x N)-b
%  Coefficients : Amplitudes des FG selon x et y  (2 x N)
%  Courbure : -A�lements des matrices de courbure des FG dans la forme [Q11; Q22; Q12]  -b
%  Sens : (option) Sens de propagation : +1 (defaut) ou -1
%                                                 (3 x N)
%
% RETOURNE
%  Champ : structure Champ
%
% VOIR AUSSI
%  Concat_Champ : pour fusionner deux structures champs
%  getChamp : pour obtenir le Nieme element Champ
%  
%

%  Si aucun argument on initialise une structure vide
if (nargin == 0)
  Champ = struct(...
      'Centre', [], ...
      'ex', [], ...
      'ez', [], ...
      'Coefficients', [], ...
      'Courbure', []);
elseif (nargin == 5)
  Champ = struct(...
      'Centre', Centre, ...
      'ex', ex, ...
      'ez', ez, ...
      'Coefficients', Coefficients, ...
      'Courbure', Courbure, ...
      'Sens', +1);
elseif (nargin == 6)
  Champ = struct(...
      'Centre', Centre, ...
      'ex', ex, ...
      'ez', ez, ...
      'Coefficients', Coefficients, ...
      'Courbure', Courbure, ...
      'Sens', Sens);
else
  error('Nombre d''arguments incorrect');
end