function Plaque = struct_Plaque(Sommets, Centre, N, Epsr)
%
% Structure de description d'un objet quadrilataire convexe de type rectangle.
% 
% Cette structure contient les informations nécessaire sur cet objet, soit:
% ARGUMENTS :
%  Sommets : position des 4 sommets du quadrilatere         (3x4)
%            dans le repere ABSOLU :
%              1____________2
%              |           |
%              |           |
%             4|___________|3
%            Une erreur est retournee si le sommet 4 n'appartient pas
%            au plan definit par les 3 premiers sommets.
%  Centre : centre de la plaque (3x1)
%  
%  Epsr : Epsilon relatif de l'objet 
%  
% RETOURNE :
%  Plaque : structure Plaque
% 



%  structure vide si aucun argument
if (nargin==0)
  Sommets = [];
   Centre = [];
        N = [];
     Epsr = [];

else
  %  ###### tests de validite des arguments
  %  les coordonnees passees sont elles 
  %  conformes a celles d'un rectangle ?
  %  cad les 4 points appartiennent ils au meme plan ?
  %  Pour ce faire on utilise l'equation du plan donnee
  %  a partir des 3 premiers points pour verifier que 
  %  le quatrieme point verifie l'equation
  %  http://en.wikipedia.org/wiki/Plane_%28mathematics%29#Define_a_plane_through_three_points
  x1 = Sommets(1,1); y1 = Sommets(2,1); z1 = Sommets(3,1);
  x2 = Sommets(1,2); y2 = Sommets(2,2); z2 = Sommets(3,2);
  x3 = Sommets(1,3); y3 = Sommets(2,3); z3 = Sommets(3,3);
  x4 = Sommets(1,4); y4 = Sommets(2,4); z4 = Sommets(3,4);
   a = y1*z2-y1*z3-z1*y2+z1*y3-z2*y3+y2*z3;
   b = -x1*z2+x1*z3-x2*z3+x3*z2-x3*z1+x2*z1;
   c = x2*y3-x2*y1-x1*y3+x3*y1+x1*y2-x3*y2;
   d = x1*z2*y3+x3*z1*y2-x3*y1*z2-x2*z1*y3-x1*y2*z3+x2*y1*z3;
  test = a*x4 + b*y4 + c*z4 + d;

  if (test ~= 0)
    error('ERREUR D''ARGUMENT : Le sommet numero 4 n''appartient pas au plan definit par les 3 premiers sommets.');
  end

  % TODO : calculer automatiquement le centre du rectangle
  %  par defaut, le centre de la plaque correspond aux milieux des deux cotes



end
%  ====================================
%  Creation de la structure
%  ====================================
Plaque = struct('Sommets', Sommets, ...
                'Centre', Centre, ...
                'N', N, ....
                'Epsr', Epsr);
