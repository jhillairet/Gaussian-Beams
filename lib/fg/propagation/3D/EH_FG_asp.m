

% Fonction :   EH_FG_asp                           
% Date     :   20/06/03                                    
% Auteur   :   A .Chabory                                  
% Etat     :   OK                          
%                                                          
%                                                          
% champ E et H rayonne par un  FG avec formule analytique     
% asymptotique                                  
%                                                          
%                                                          
% - Entrees                                                
%   P            : Points d'observations                (3 x NP)                    
%   Courbure     : Matrice de courbure du fsc           (3 x NP)  
%   Coefficients : Coefficients en co et cross polar    (2)
%   Freq         : Frequence                            (1)        
%   Epsr         : Permitivite relative du milieu       (1)         
%                                                         
% - Sorties                                                
%   E            : Champ electrique rayonne en P        (3 x NP)          



function EH=EH_FG_asp(P,Courbure,Coefficients,Freq,Epsr)
lam0=2.997925e8./Freq;
km=2.*pi./lam0.*sqrt(Epsr);
Zn=120.*pi./sqrt(Epsr);

NP=size(P,2);

% passage en coord spheriques
r = sqrt(P(1,:).^2+P(2,:).^2+P(3,:).^2);
theta=atan2(sqrt(P(1,:).^2+P(2,:).^2),P(3,:));
phi = atan2(P(2,:),P(1,:));

% matrice de courbure inverse
invCb0=inv([Courbure(1) Courbure(3);...
            Courbure(3) Courbure(2)]);

% champ E
EH(1:3,:) =  ( ( ones(3,1)*(exp(-j.*km.*r+...
   j.*km.*sin(theta).^2./2.*(invCb0(1,1).*cos(phi).^2 + invCb0(2,2).*sin(phi).^2+...
   2.*invCb0(1,2).*cos(phi).*sin(phi) ) )./(sqrt(-1./det(invCb0)).*r))).* ...
   ( [cos(theta).*Coefficients(1);...
      cos(theta).*Coefficients(2);...
      -sin(theta).*(cos(phi).*Coefficients(1)+...
      sin(phi).*Coefficients(2))]  )); 

EH(4:6,:) = cross( P.*(ones(3,1)*(1./r)) ,EH(1:3,:))./Zn;
