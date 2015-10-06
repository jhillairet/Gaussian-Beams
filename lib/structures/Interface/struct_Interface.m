function Interface = struct_Interface(...
    D, Eps, Epsleft, Epsright, Pleft, Pright, Nleft2right)
% Interface = struct_Interface(D, Eps, Epsleft, Epsright, Pleft, Pright, Nleft2right)
%
% Structure Interface
%
% Cette structure permet de d-Aécrire une interface multicouche de n couches.-b
%
% PARAM-AÈTRES-b
%  D : -Aépaisseur de la nième interface [m]                          (1 x 1)-b
%  Eps : permittivit-Aé relative de la nième couche                   (1 x 1)   -b
%  Epsleft : permittivit-Aé relative de la nième couche (gauche)      (1 x 1)-b
%  Epsright: permittivit-Aé relative de la nième couche (droite)      (1 x 1)-b
%  Pleft : vecteur position des points de la ni-Aème couche gauche    (3 x NP)-b
%  Pright: vecteur position des points de la ni-Aème couche droite    (3 x NP)-b
%  Nleft2right: vecteur unitaire directeur externe -Aà la paroi (gauche->droite)-b
%
%
% RETOURNE : 
%  structure Interface <struct>

Interface = struct(...
    'D', D,...
    'Eps', Eps, ...
    'Epsleft', Epsleft, ...
    'Epsright', Epsright, ...
    'Pleft', Pleft, ...
    'Pright', Pright,...
    'Nleft2right', Nleft2right);
