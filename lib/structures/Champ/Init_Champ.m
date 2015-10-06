

% Fonction   :    Init_Champ
% Date       :    30/06/03
% Auteur     :    A .Chabory
% Etat       :    ok
%
%
% Initialisation ou masquage dun champ 
%
% [Champ]=Init_Champ
% [Champ]=Init_Champ(Champ1,Mq)
%
% - Entrees
%    Champ1 : (Optionnel)     (struct Champ)
%             Champ initial 
%    Mq     : (Optionnel)
%             Indices a conserver 
%
% - Sorties
%    Champ  : Champ concatene (struct Champ)


function [Champ]=Init_Champ(Champ1,Mq)

switch nargin
case 0
   Champ.Centre=[];
   Champ.ex=[];
   Champ.ez=[];
   Champ.Courbure=[];
   Champ.Coefficients=[];
   
case 2
   Champ.Centre=Champ1.Centre(:,Mq);
   Champ.ex=Champ1.ex(:,Mq);
   Champ.ez=Champ1.ez(:,Mq);
   Champ.Courbure=Champ1.Courbure(:,Mq);
   Champ.Coefficients=Champ1.Coefficients(:,Mq);

end



