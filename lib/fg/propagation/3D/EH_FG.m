% Fonction :   EH_FG                           
% Date     :   2003-2007                                   
% Auteur   :   A .Chabory    / J.Hillairet                              
% Etat     :   ok                          
%                                                          
%                                                          
%  Champ electrique et magnetique rayonne par un  FG :
%  Formule analytique approximation paraxiale                                  
% 
%  [E,H]=EH_FG(P,FG,Freq,Epsr)
%                                                          
% - Entrees                                                
%    P            : Points d'observations                (3 x NP)
%    FG           : Structure Champ FG                   <struct>
%    Freq         : Frequence                            (1)        
%    Epsr         : Permitivite relative du milieu       (1)         
%                                                         
% - Sorties                                                
%    E,H            : Champ electrique et magnetique en P  (3 x NP)        
function [E,H] = EH_FG(P,FG,Freq,Epsr)

lam0=2.997925e8./Freq;
km=2.*pi./lam0.*sqrt(Epsr);
Zn=120.*pi./sqrt(Epsr);

NP=size(P,2);

    Courbure = FG.Courbure;
Coefficients = FG.Coefficients;

% matrice de courbure sur tous les points
invCb0=inv([Courbure(1) Courbure(3);...
      Courbure(3) Courbure(2)]);
detCb=1./( (invCb0(2,2)+P(3,:)).*(invCb0(1,1)+P(3,:))-invCb0(1,2).^2 );  
Cb11=(invCb0(2,2)+P(3,:)).*detCb;
Cb22=(invCb0(1,1)+P(3,:)).*detCb;
Cb12=-invCb0(1,2).*detCb;
Coef1 = Cb11.*P(1,:)+Cb12.*P(2,:);
Coef2 = Cb22.*P(2,:)+Cb12.*P(1,:);

u=ones(3,1)* (sqrt(-detCb)./sqrt(-1./det(invCb0)).*...
   exp(-j.*km./2.*(Coef1.*P(1,:)+Coef2.*P(2,:))...
   -j.*km.*P(3,:)));

E = (([Coefficients(1:2);0]*ones(1,NP))+...
   [zeros(1,NP);zeros(1,NP);-(Coefficients(1).*Coef1...
      +Coefficients(2).*Coef2)]).*u;

H = (([-Coefficients(2);Coefficients(1);0]*ones(1,NP))+...
         [zeros(1,NP);zeros(1,NP);...
            -(-Coefficients(2).*Coef1...
            +Coefficients(1).*Coef2)]).*u./Zn;


