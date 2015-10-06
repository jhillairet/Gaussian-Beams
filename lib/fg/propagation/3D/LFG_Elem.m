% Fonction   :    LFG_Elem
% Date       :    22/08/03
% Auteur     :    A .Chabory
% Etat       :    Ok
%
%
% Interaction elementaire du LFG sur une interface 
% dielectrique ellipsoidale
%                
% [Fsc_Ref,Fsc_Tr]=TFG_Elem(Fsc,Freq,Objet,Interface,Sens);
%
% - Entrees
%    Fsc       : Faisceau incident          (struct Champ)            
%    Freq      : Frequence                  (1)
%    Objet     : Objet                      (stuct Objet)
%    Interface : Interface eclairee         (1)
%    Sens      : Sens de l'interaction      (1) 
%                +1 : Interface->Interface+1
%                -1 : Interface-1->Interface
%
% - Sorties          
%    Fsc_Ref   : Faisceau reflechi          (struct Champ)            
%    Fsc_Tr    : Faisceau transmis          (struct Champ)            


function [Fsc_Ref,Fsc_Tr]=LFG_Elem(Fsc,Freq,Objet,Interface,Sens);

Fsc_Ref=Init_Champ;
Fsc_Tr=Init_Champ;

Milieu=Interface+(1-Sens)./2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul des parametres initiaux du faisceau
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Sens.*Milieu==-1 | Sens.*Milieu==length(Objet.Epsr)  
   return;
end


Of=Fsc.Centre;
Zf=Fsc.ez;
Xf=Fsc.ex;
Qf=Fsc.Courbure;
Af=Fsc.Coefficients;

Epsr_Inc=Objet.Epsr(Milieu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul du point I																						 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parametres de l interface suivante rencontree

Milieu_Tr=Milieu+Sens;
%Interface=Milieu+(Sens-1)./2;
A=Objet.Axes(1,Interface);
B=Objet.Axes(2,Interface);
C=Objet.Axes(3,Interface);

z_int=Objet.Z_Interface(Interface);
z_lim=Objet.Z_Limite(Interface);

Epsr_Tr=Objet.Epsr(Milieu_Tr);


I=Intersection(Zf,Of,Objet,Interface);

%  error('test')
if isempty(I)
   return
end

% Angles dincidence de transmission et de reflection
nI=Normale(I,Objet,Interface,Sens);

nI=2.*[I(1)./A.^2;I(2)./B.^2;(I(3)-C-z_int)./C.^2];
nI=Sens.*sign(I(3)-C-z_int).*nI./norm(nI);

thetai=acos(dot(nI,Zf));
thetar=thetai;
thetat=asin(sqrt(Epsr_Inc./Epsr_Tr).*sin(thetai));

if imag(thetat)~=0
   disp('reflexion totale');
   return
end

if abs(thetai)>1e-4
   v= cross(Zf,nI);
   v=v./norm(v);
else
   v=cross(nI,Fsc.ex);
   v=v./norm(v);   
end
u=cross(v,nI);

Fsc_Tr.Centre=I;
Fsc_Ref.Centre=I;
Fsc_Tr.ez=-sin(thetat).*u+cos(thetat).*nI;
Fsc_Ref.ez=-sin(thetai).*u-cos(thetai).*nI;
Fsc_Tr.ex=cos(thetat).*u+sin(thetat).*nI;
Fsc_Ref.ex=-cos(thetai).*u+sin(thetai).*nI;

% Matrice de courbure du fsc inc 

MFsc=[cross(v,Zf) v].' * [Fsc.ex  cross(Fsc.ez,Fsc.ex)];
dist=norm(Fsc.Centre-I);
Qinc0=MFsc.'*[Fsc.Courbure(1) Fsc.Courbure(3);Fsc.Courbure(3) Fsc.Courbure(2)]*MFsc;
Qinc=inv(  inv(Qinc0) + dist.*diag([1 1]));


% Matrice de courbure de la surface en I

Qsig11=2.*(u(1).^2./A.^2+u(2).^2./B.^2+u(3).^2./C.^2);
Qsig22=2.*(v(1).^2./A.^2+v(2).^2./B.^2+v(3).^2./C.^2);
Qsig12=u(1).*v(1)./A.^2+u(2).*v(2)./B.^2+u(3).*v(3)./C.^2;
Qsig= -sign(C).*Sens./(2.*sqrt(I(1).^2./A.^4+I(2).^2./B.^4+(I(3)-C-z_int).^2./C.^4)).*...
   [Qsig11 Qsig12;Qsig12 Qsig22];


h=sqrt(Epsr_Tr).*cos(thetat)-sqrt(Epsr_Inc).*cos(thetai);

Fsc_Tr.Courbure(1,1)=1./(sqrt(Epsr_Tr).*cos(thetat).^2).*...
   ( sqrt(Epsr_Inc).*cos(thetai).^2.*Qinc(1,1)+h.*Qsig(1,1));

Fsc_Tr.Courbure(3,1)=1./(sqrt(Epsr_Tr).*cos(thetat)).*...
   (cos(thetai).*Qinc(1,2)+h.*Qsig(1,2));


Fsc_Tr.Courbure(2,1)=1./sqrt(Epsr_Tr).*(sqrt(Epsr_Inc).*Qinc(2,2)+h*Qsig(2,2));

Fsc_Ref.Courbure=[Qinc(1,1)-2*Qsig(1,1)./cos(thetai);...
      Qinc(2,2)-2*Qsig(2,2).*cos(thetai);2.*Qsig(1,2)-Qinc(1,2)];

           
% Coefficients de Fresnel
           
T_TE=2.*cos(thetai)./(cos(thetai)+sqrt(Epsr_Tr./Epsr_Inc).*cos(thetat));
T_TM=2.*cos(thetai)./(cos(thetat)+sqrt(Epsr_Tr./Epsr_Inc).*cos(thetai));
R_TE=T_TE-1;
R_TM=sqrt(Epsr_Tr./Epsr_Inc).*T_TM-1;

lambdam=2.997925e8./(Freq*sqrt(Epsr_Inc));
km=2.*pi./lambdam;
Amp_Inc=sqrt(-det(Qinc))./sqrt(-det(Qinc0)).*exp(-j.*km.*dist).* (MFsc*Fsc.Coefficients);
Fsc_Tr.Coefficients=diag([T_TM T_TE])*Amp_Inc;
Fsc_Ref.Coefficients=diag([R_TM R_TE])*Amp_Inc;










   









   


