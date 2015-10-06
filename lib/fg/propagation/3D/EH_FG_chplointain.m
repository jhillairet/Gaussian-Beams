% Fonction :   EH_FG_chplointain                           
% Date     :   02/2007
% Auteur   :   Julien Hillairet                          
% Etat     :   ok            
%                                                          
%                                                          
%  Champ electrique et magnetique rayonne par un  FG :
%  Formule analytique approximation champ lointain
% 
%  EH=EH_FG_chplointain(P,FG,Freq,Epsr)
%                                                          
% - Entrees                                                
%    P            : Points d'observations                (3 x NP)                   
%    FG           : Structure Champ FG                   <struct>
%    Freq         : Frequence                            (1)        
%    Epsr         : Permitivite relative du milieu       (1)         
%                                                         
% - Sorties                                                
%    EH            : Champ electrique et magnetiqueen P  (6 x NP)          
function [E,H] = EH_FG_chplointain(P,FG,Freq,Epsr)

lambda0 = 2.997925e8./Freq;
 k0 = 2.*pi./lambda0;
 km = 2.*pi./lambda0.*sqrt(Epsr);
 Zn = 120.*pi./sqrt(Epsr);
NbP = size(P,2);
  R = sqrt(sum(P.^2));

    Q0 = FG.Courbure;
Coeffs = FG.Coefficients;

% matrice de courbure sur tous les points
detQ0 = Q0(1)*Q0(2)-Q0(3).^2;
invQ0 = [Q0(2);Q0(1);-Q0(3)]./(ones(3,1)*detQ0);

u = ones(3,1)*(...
    P(3,:)./R ... % z/R=cos(theta)
    ./nsqrt(detQ0) ...
    .* exp(-j*km*R) ./ R ...
    .* exp(j*km/2./R.^2.*(...
          invQ0(1).*P(1,:).^2 ...
        + invQ0(2).*P(2,:).^2 ...
        + 2*invQ0(3).*P(1,:).*P(2,:) )));



PolE = [Coeffs(1:2)*ones(1,NbP); ...
        -(Coeffs(1)*P(1,:)./P(3,:) + Coeffs(2)*P(2,:)./P(3,:))];

%  NB : Les termes de second ordre dans polH sont negligeables en zone lointaine. 
%  On pourrait les supprimer à priori.
%  PolH = 1./Zn ...
%        .* [- Coeffs(2).*P(3,:)./R ...
%              + (P(2,:).^2*Coeffs(2) - P(1,:).*P(2,:).*Coeffs(1))./(R.*P(3,:)); ...
%            + Coeffs(1).*P(3,:)./R ...
%              - (P(1,:).^2*Coeffs(1) - P(1,:).*P(2,:).*Coeffs(2))./(R.*P(3,:)); ...
%            (P(1,:)*Coeffs(2) - P(2,:)*Coeffs(1))./R];

  a_h_x = -(P(1,:).*P(2,:).*Coeffs(1) + (P(2,:).^2+P(3,:).^2).*Coeffs(2))./(P(3,:).*R);
  a_h_y = +(P(1,:).*P(2,:).*Coeffs(2) + (P(1,:).^2+P(3,:).^2).*Coeffs(1))./(P(3,:).*R);
  a_h_z = +P(1,:)./R.*Coeffs(2) - P(2,:)./R.*Coeffs(1);

PolH = 1./Zn .* [a_h_x;a_h_y;a_h_z];


E = PolE.*u;
H = PolH.*u;


