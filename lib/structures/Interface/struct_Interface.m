function Interface = struct_Interface(...
    D, Eps, Epsleft, Epsright, Pleft, Pright, Nleft2right)
% Interface = struct_Interface(D, Eps, Epsleft, Epsright, Pleft, Pright, Nleft2right)
%
% Structure Interface
%
% Cette structure permet de d-A�crire une interface multicouche de n couches.-b
%
% PARAM-A�TRES-b
%  D : -A�paisseur de la ni�me interface [m]                          (1 x 1)-b
%  Eps : permittivit-A� relative de la ni�me couche                   (1 x 1)   -b
%  Epsleft : permittivit-A� relative de la ni�me couche (gauche)      (1 x 1)-b
%  Epsright: permittivit-A� relative de la ni�me couche (droite)      (1 x 1)-b
%  Pleft : vecteur position des points de la ni-A�me couche gauche    (3 x NP)-b
%  Pright: vecteur position des points de la ni-A�me couche droite    (3 x NP)-b
%  Nleft2right: vecteur unitaire directeur externe -A� la paroi (gauche->droite)-b
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
