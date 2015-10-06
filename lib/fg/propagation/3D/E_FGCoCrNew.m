

% Fonction :   E_FG                           
% Date     :   20/06/03                                    
% Auteur   :   A .Chabory                                  
% Etat     :   ok                          
%                                                          
%                                                          
%  Champ electrique rayonne par un  FG :
%  Formule analytique approximation paraxiale                                  
% 
%  [Eco,Ecr] = E_FGCoCr(P,Courbure,Coefficients,Freq,Epsr)
%                                                          
% - Entrees                                                
%    P            : Points d'observations                (3 x NP)                    
%    Courbure     : Matrice de courbure du fsc           (3 x NP)  
%    Coefficients : Coefficients en co et cross polar    (2)
%    Freq         : Frequence                            (1)        
%    Epsr         : Permitivite relative du milieu       (1)         
%                                                         
% - Sorties                                                
%    Eco, Ecr     : Champ electrique rayonne en P,
%                   en co- et cross-polar                (3 x NP)          



function [Eco,Ecr] = E_FGCoCr(P,InvCourbure,Coefficients,Freq,Epsr)
lam0=2.997925e8./Freq;
km=2.*pi./lam0.*sqrt(Epsr);

NP=size(P,2);

% matrice de courbure sur tous les points
detInvCb = InvCourbure(1,:).*InvCourbure(2,:)-InvCourbure(3,:).^2;
Cb11=InvCourbure(1,:)+P(3,:);
Cb22=InvCourbure(2,:)+P(3,:);
detCb=1./( Cb22.*Cb11-InvCourbure(3,:).^2 );  
Cb11=Cb11.*detCb;
Cb22=Cb22.*detCb;
Cb12=-InvCourbure(3,:).*detCb;

CoefA = exp(-j.*km./2.*( Cb11.*P(1,:).^2+Cb22.*P(2,:).^2+2.*Cb12.*P(1,:).*P(2,:)  )...
   -j.*km.*P(3,:));
CoefB = ones(3,1)* (sqrt(-detCb)./sqrt(-1./detInvCb).*CoefA);

Eco = (([Coefficients(1);0;0]*ones(1,NP))+...
   [zeros(1,NP);zeros(1,NP);-1.*(Coefficients(1).*(Cb11.*P(1,:)+Cb12.*P(2,:)))]).*CoefB;

Ecr = (([0;Coefficients(2);0]*ones(1,NP))+...
   [zeros(1,NP);zeros(1,NP);-1.*(Coefficients(2).*(Cb22.*P(2,:)+Cb12.*P(1,:)))]).*CoefB;

