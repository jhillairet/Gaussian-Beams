

% Fonction :   EPcc_FG                           
% Date     :   20/06/03                                    
% Auteur   :   A .Chabory                                  
% Etat     :   Amplitude du vct de Poynt a verifier                          
%                                                          
%                                                          
%  Champ electrique rayonne par un  FG unitaire 
%  en co et cross polar
%  Formule analytique approximation paraxiale
%  Calcul du vecteur de poynting
%                                                          
%  [E,Ec,Pt]=EPcc_FG(P,Courbure,Freq,Epsr)
%
% - Entrees                                                
%   P            : Points d'observations                (3 x NP)                    
%   Courbure     : Matrice de courbure du fsc           (3 x NP)  
%   Freq         : Frequence                            (1)        
%   Epsr         : Permitivite relative du milieu       (1)         
%                                                         
% - Sorties                                                
%   E            : Champ electrique en co polar
%   Ec           : Champ electrique en cr polar
%   Pt           : vecteur de poynting en co polar



function [E,Ec,Pt]=EPcc_FG(P,Courbure,Freq,Epsr)

lam0=2.997925e8./Freq;
km=2.*pi./lam0.*sqrt(Epsr);
Zn=120.*pi./sqrt(Epsr);

NP=size(P,2);

% matrice de courbure sur tous les points
invCb0=inv([Courbure(1) Courbure(3);...
      Courbure(3) Courbure(2)]);
detCb=1./( (invCb0(2,2)+P(3,:)).*(invCb0(1,1)+P(3,:))-invCb0(1,2).^2 );  
Cb11=(invCb0(2,2)+P(3,:)).*detCb;
Cb22=(invCb0(1,1)+P(3,:)).*detCb;
Cb12=-invCb0(1,2).*detCb;

u=ones(3,1)* (sqrt(-detCb)./sqrt(-1./det(invCb0)).*...
   exp(-j.*km./2.*(Cb11.*P(1,:).^2+Cb22.*P(2,:).^2+2.*Cb12.*P(1,:).*P(2,:)  )...
   -j.*km.*P(3,:)));
vectz=-(Cb11.*P(1,:)+Cb12.*P(2,:));
vectz_c=-(Cb22.*P(2,:)+Cb12.*P(1,:));
E = [ones(1,NP);zeros(1,NP);vectz].*u;
Ec =[zeros(1,NP);ones(1,NP);vectz_c].*u;

Pt=0.5.*real([-vectz;-vectz_c;ones(1,NP)]);
%Pt=Pt./(ones(3,1)*sum(abs(Pt).^2,1));

