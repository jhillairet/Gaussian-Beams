function Objet = struct_Objet(Type, Axes, Z_Interface, Z_Limite, Epsr)
%
% Structure de description d'un objet diélectrique multicouche (lentille,
% radomes,...). 
%
% Objet = struct_Objet(Type, Axes, Z_Interface, Z_Limite, Epsr)
%
% Cette stucture correspond simplement à une succession d'interfaces
% ellipsoidales centrées sur l'axe Oz du repère absolu séparant des milieux
% diélectriques.
% Les milieux et les interfaces sont numérotés de façon croissante selon
% l'axe Oz. Le premier milieu est le milieu incident, le dernier est le
% milieu transmis. 
%
% Pour plus de détails sur les caractéristiques des Objets, cf la thèse
% d'Alexandre Chabory, Annexe G.
%
% PARAMETRES :
%  Type : Type de l'objet 'e' pour ellipsoidal, 'o' pour ogival   (char)
%  Axes : demi grands axes des ellipsoides selon x, y et z        (3 x Ni)
%  Z_Interface : Intersection des ellipsoides avec l'axe Oz       (1 x Ni)
%  Z_Limite : Limite des ellipsoides avec l'axe OZ                (1 x Ni)
%  Epsr : Constantes dielectriques relatives des milieux successifs
%                                                                 (1 x Ni+1)
%
% RETOURNE :
%  Objet : structure Objet

if (nargin == 0)
  Objet = struct(...
      'Type', [], ...
      'Axes', [], ...
      'Z_Interface', [], ...
      'Z_Limite', [], ...
      'Epsr', []);
elseif (nargin == 5)
  Objet = struct(...
      'Type', Type, ...
      'Axes', Axes, ...
      'Z_Interface', Z_Interface, ...
      'Z_Limite', Z_Limite, ...
      'Epsr', Epsr);
end