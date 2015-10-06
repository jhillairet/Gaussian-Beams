% Fonction :   EH_Sigma_FG                           
% Date     :   13/01/03                                    
% Auteur   :   A .Chabory                                  
% Etat     :   OK                          
%                                                          
%                                                          
%  Champs electrique et magnetique rayonne par une somme de FG :
%  Formule analytique approximation paraxiale                                  
%                                                          
%  [E, H] = EH_Sigma_FG( P, Champ, Freq, Epsr)                                                         
%
% - Entrees                                                
%    P         : Points d'observations                (3 x NP)                    
%    Champ     : Description du champ en fsc gassien  (struct Champ)         
%    Freq      : Frequence                            (1)        
%    Epsr      : Permitivite relative du milieu       (1)         
%                                                         
% - Sorties
%    E     :   Champ electtique rayonne en P          (3 x NP)          
%    H     :   Champ magnetique rayonne en P          (3 x NP)          



function [E,H] = EH_Sigma_FG(P, Champ, Freq, Epsr)
   
E = zeros(3,size(P,2));
H = zeros(3,size(P,2));

%Progress Bar
PB = progress('init', 'Sommation des faisceaux gaussiens');

for ind=1:size(Champ.Coefficients,2)
   % Matrice de passage dans le repere associe au fsc
   Mfsc=[Champ.ex(:,ind) cross(Champ.ez(:,ind),Champ.ex(:,ind)) Champ.ez(:,ind)];
   
   % Coordonnees de Q dans le repere fsc
   Pfsc=Mfsc.'*(P-Champ.Centre(:,ind)*ones(1,length(P)));
   
   [E0,H0] = EH_FG(Pfsc,getChamp(Champ,ind),Freq,Epsr); 
   E = E + Mfsc*E0;
   H = H + Mfsc*H0;
   
   % MAJ Progress Bar
   PB = progress(PB, ind/size(Champ.Coefficients,2));
end

