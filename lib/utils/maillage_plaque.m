function [Q,Nx,Ny,Nz] = maillage_plaque(Plaque, dx, dy)
% maillage d'un rectangle a la precision dx,dy
% 
% PARAMETRES :
%  Plaque : structure Plaque contenant les champs
%      dx : echantillonnage requis en x
%      dy : echantillonnage requis en y
%      
% RETOURNE :
%       Q : vecteur contenant les points echantillonne [Qx;Qy;Qz]
%  
%  TODO : 
%  pour l'instant, la fonctione ne genere les points que pour les plaques
%  situee dans les plan z=Plaque.Centre(3)=Cte 
%  et uniquement pour les rectangles, dont le centre est au centre.
%  


X1 = Plaque.Centre(1) + abs(Plaque.Sommets(1,1) - Plaque.Sommets(1,2))/2;
X2 = Plaque.Centre(1) - abs(Plaque.Sommets(1,1) - Plaque.Sommets(1,2))/2;
Y1 = Plaque.Centre(2) + abs(Plaque.Sommets(2,2) - Plaque.Sommets(2,3))/2;
Y2 = Plaque.Centre(2) - abs(Plaque.Sommets(2,2) - Plaque.Sommets(2,3))/2;

X1 = max(Plaque.Sommets(1,:));
X2 = min(Plaque.Sommets(1,:));
Y1 = max(Plaque.Sommets(2,:));
Y2 = min(Plaque.Sommets(2,:));

% on genere tout d'abord une plaque de longueur (X1-X2) 
% et de largeur (Y1-Y2)
x_Q = [X2:dx:X1];
y_Q = [Y2:dy:Y1];
x_Q = [X2:dx:X1];
y_Q = [Y2:dy:Y1];
Q = stack(x_Q, y_Q, Plaque.Centre(3));

% puis on effectue un changement de repere pour mettre le centre de la plaque
% au bon endroit 
Q = Plaque.Centre*ones(1,length(Q)) + Q ;

% taille des composantes de la pile
Nx = length(x_Q);
Ny = length(y_Q);
Nz = 1;