% Fonction   :    LFG_Inter
% Date       :    22/08/03
% Auteur     :    A .Chabory
% Etat       :    pas ok
%
%
% Lancer de FG pour une interface dielectrique 
%                
% [Champ_Ref,Champ_Tr,Champ_NT]=LFG_Inter(Champ_Inc,Freq,Objet,Interface,Sens)
%
% - Entrees
%    Champ_Inc : Champ incident             (struct Champ)            
%    Freq      : Frequence                  (1)
%    Objet     : Objet                      (stuct Objet)
%    Interface : Interface eclairee         (1)
%    Sens      : Sens de l'interaction      (1) 
%                +1 : Interface->Interface+1
%                -1 : Interface+1->Interface
%
% - Sorties          
%    Champ_Ref : Champ reflechi             (struct Champ)
%    Champ_Tr  : Champ transmis             (struct Champ)
%    Champ_NT  : Champ Inc non traitee      (struct Champ)

function [Champ_Ref,Champ_Tr,Champ_NT]=LFG_Inter(Champ_Inc,Freq,Objet,Interface,Sens)

Champ_Ref=Init_Champ;
Champ_Tr=Init_Champ;
Champ_NT=Init_Champ;

for ind=1:size(Champ_Inc.Coefficients,2)
   Fsc_Inc=Init_Champ(Champ_Inc,ind);
   [Fsc_Ref,Fsc_Tr]=LFG_Elem(Fsc_Inc,Freq,Objet,Interface,Sens);
   if isempty(Fsc_Ref.Coefficients)
      Champ_NT=Concat_Champ(Champ_NT,Fsc_Inc);
   else
      Champ_Ref=Concat_Champ(Champ_Ref,Fsc_Ref);
      Champ_Tr=Concat_Champ(Champ_Tr,Fsc_Tr);
   end
end









   









   


