

% Fonction   :    Concat_Champ
% Auteur     :    A .Chabory, J. Hillairet
% Etat       :    ok
%
%
% Concatenation de  2 structures Champ FG
%
% [Champ]=Concat_Champ(Champ1,Champ2)
%
% - Entrees
%    Champ1 : Premier champ   (struct Champ)
%    Champ2 : Deuxieme champ  (struct Champ)
%
% - Sorties
%    Champ  : Champ concatene (struct Champ)

function [Champ]=Concat_Champ(Champ1,Champ2)

if isempty(Champ1)
  Champ.Centre=[Champ2.Centre];
  Champ.ex=[Champ2.ex];
  Champ.ez=[Champ2.ez];
  Champ.Courbure=[Champ2.Courbure];
  Champ.Coefficients=[Champ2.Coefficients];
elseif isempty(Champ2)
  Champ.Centre=[Champ1.Centre];
  Champ.ex=[Champ1.ex];
  Champ.ez=[Champ1.ez];
  Champ.Courbure=[Champ1.Courbure];
  Champ.Coefficients=[Champ1.Coefficients];
else  
  Champ.Centre=[Champ1.Centre Champ2.Centre];
  Champ.ex=[Champ1.ex Champ2.ex];
  Champ.ez=[Champ1.ez Champ2.ez];
  Champ.Courbure=[Champ1.Courbure Champ2.Courbure];
  Champ.Coefficients=[Champ1.Coefficients Champ2.Coefficients];
end

