% Fonction :   EH_Sigma_FG_chplointain                           
% Date     :   02/2007                                
% Auteur   :   J. Hillairet                        
% Etat     :   dev                          
%                                                          
%                                                          
%  Champs electrique et magnetique rayonne par une somme de FG :
%  Formule analytique approximation champ lointain                                  
%                                                          
%  [E, H] = EH_Sigma_FG_chplointain( P, Champ, Freq, Epsr)                                                         
%
% PARAMETRES :                                            
%    P         : Points d'observations                (3 x NP)                    
%    Champ     : Description du champ en fsc gassien  (struct Champ)         
%    Freq      : Frequence                            (1)        
%    Epsr      : Permitivite relative du milieu       (1)         
%                                                         
% RETOURNE :
%    E     :   Champ electtique rayonne en P          (3 x NP)          
%    H     :   Champ magnetique rayonne en P          (3 x NP)          
function [E,H] = EH_Sigma_FG_chplointain(P, Champ, Freq, Epsr)
   
E = zeros(3,size(P,2));
H = zeros(3,size(P,2));

  c = 2.997925e8;
 k0 = 2*pi/c*Freq;

%Progress Bar
PB = progress('init', 'Sommation des faisceaux gaussiens : approx. chp. lointain');

for ind=1:size(Champ.Coefficients,2)
%    %  Pour utilser la formulation champ lointain  
%    %  la matrice de courbure complexe dans le plan du waist
%    %  doit etre imaginaire pure (zW0 = 0)
%    %  
%    %  il convient donc de faire attention
%    %  a bien se placer dans ce plan.
%    Q = Champ.Courbure(:,ind);
%    %   si on est pas dans le plan du waist (zW0 ~= 0) 
%    if (real(Q(1)) ~= 0)
%      %  La premiere etape est d'eventuellement diagonaliser la matrice 
%      %  de courbure si elle n'est pas deja diagonale.
%      Q = [Q(1), Q(3); Q(3), Q(2)];
%      if (Q(3) ~= 0)
%        [V,Q] = eig(Q);
%      end
%      % cela permet de remonter aux informations sur zW0
%      zW0_x = -real(1./Q(1,1));
%      zW0_y = -real(1./Q(2,2));
%      % on modifie alors le centre du FG 
%      % en "reculant" le centre du repere propre 
%      % pour qu'il corresponde au centre "physique" du FG
%      Champ.Centre(:,ind) = Champ.Centre(:,ind) + max(zW0_x,zW0_y)*Champ.ez(:,ind);
%  
%  %      % il faut aussi modifier la matrice de courbure
%  %      Q0(1,1) = 1/(j*imag(1/Q(1,1)));
%  %      Q0(2,2) = 1/(j*imag(1/Q(2,2)));
%  %      if exist('V')
%  %        Q0 = V*Q0*inv(V);
%  %      end
%  %      Champ.Courbure(:,ind) = [Q0(1,1);Q0(2,2);Q0(1,2)];
%  
%    end


   % Matrice de passage dans le repere associe au fsc
   Mfsc=[Champ.ex(:,ind) cross(Champ.ez(:,ind),Champ.ex(:,ind)) Champ.ez(:,ind)];
   
   % Coordonnees de Q dans le repere fsc
   Pfsc=Mfsc.'*(P-Champ.Centre(:,ind)*ones(1,length(P)));
   
   [E0,H0] = EH_FG_chplointain(Pfsc,getChamp(Champ,ind),Freq,Epsr); 
   E = E + Mfsc*E0;
   H = H + Mfsc*H0;
   
   % MAJ Progress Bar
   PB = progress(PB, ind/size(Champ.Coefficients,2));
end

