function ChampFGC = struct_ChampFGC(Centre, ex, ez, CbF, CbS, Beta, Vect, CoeffJ, CoeffM)
%
% Structure de champ exprim� comme une somme de Faisceau Gaussien Conforme
% (FGC)
%
%  ChampFGC = struct_ChampFGC(...
%     Centre, ...
%     ex, ...
%     ez, ...
%     CbF, ...
%     CbS, ...
%     Beta, ...
%     Vect, ...
%     CoeffJ, ...
%     CoeffM);
%
% ARGUMENTS :
%  Centre : Coordonn�es des centres des FGC    (3 x Nf)
%  ex     : vecteur ex des reperes des FGC     (3 x Nf)
%  ez     : vecteur ez des reperes des FGC     (3 x Nf)
%  CbF    : Matrices de courbure des FGC       (3 x Nf)
%  CbS    : Matrices de courbure de la surface (3 x Nf)
%  Beta   : Termes de phase lin�aire           (2 x Nf)
%  Vect   : Orientation des courants           (3 x Nf)
%  CoeffJ : Amplitudes des courants �lectriques [perp;prll](2 x Nf)
%  CoeffM : Amplitudes des courants magn�tiques [perp;prll](2 x Nf)
%
% RETOURNE :
%  ChampFGC : structure ChampFGC

ChampFGC = struct( ...
    'Centre', Centre,...
    'ex', ex, ...
    'ez', ez, ...
    'CbF', CbF, ...
    'CbS', CbS, ...
    'Beta', Beta, ...
    'Vect', Vect, ...
    'CoeffJ', CoeffJ, ...
    'CoeffM', CoeffM);
        