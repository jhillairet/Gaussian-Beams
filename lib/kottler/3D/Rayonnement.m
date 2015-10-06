% Fonction  :   Rayonnement                               
% Date      :   13/01/03                                  
% Auteur    :   A .Chabory                                
% Etat      :   ok                                                       
%                                                         
% Rayonnement de J et M a partir des formules de Kottler  
%                                                         
% [E,H]=Rayonnement(P,Q,dS,J,M,Freq,Epsr,Option)
%
% - Entrees                                               
%   P      :   Points d'observations             (3 x NP)                 
%   Q      :   Points courants sur la surface    (3 x NQ)         
%   dS     :   Elements de surface au points Q   (1 x NQ)   
%   J      :   Courant electrique en Q           (3 x NQ)      
%   M      :   Courant magnetique en Q           (3 x NQ)      
%   Freq   :   Frequence                         (1)      
%   Epsr   :   Permitivite relative du milieu    (1)       
%   Option :   1 : E et H                        (1)       
%              2 : E seul                               
%              3 : E seul champ lointain                
%                                                         
% - Sorties                                               
%   E    :   Champ electrique rayonne en P       (3 x NP)            
%   H    :   Champ magnetique rayonne en P       (3 x NP)        
%                                                         

function [E,H]=Rayonnement(P,Q,dS,J,M,Freq,Epsr,Option)

% Constantes
lam0=2.997925e8./Freq;
lamm=lam0./sqrt(Epsr);
k0=2.*pi./lam0;
km=2.*pi./lamm;
Z0=120.*pi;

JdB=10.*log10(sum(abs(J).^2,1));
JdB=JdB-max(JdB);
MdB=10.*log10(sum(abs(M).^2,1));
MdB=MdB-max(MdB);
Mq=find(JdB>-80 | MdB>-80);
Q=Q(:,Mq);
dS=dS(:,Mq);
J=J(:,Mq);
M=M(:,Mq);



switch Option
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Calcul de E et H valable partout %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 1
   % Estimation du temps de calcul
   tic
   r = ones(3,1)*(...
       sqrt((P(1,1) - Q(1,:)).^2 + (P(2,1) - Q(2,:)).^2 + (P(3,1) - Q(3,:)).^2)...
       );
   r1 = [P(1,1)-Q(1,:);P(2,1)-Q(2,:);P(3,1)-Q(3,:)]./r;
   aa = 1 - j./(km.*r) - 1./(km.*r).^2;
   bb = -1 + 3.*j./(km.*r) + 3./(km.*r).^2;
   cc = 1 - j./(km.*r);
   
   E(:,1) = km./(4.*pi.*j) .* ...
       sum(...
       (Z0./sqrt(Epsr) * ...
       (aa.*J + bb.*(ones(3,1)*dot(r1,J)).*r1) - ...
       cc.*cross(r1,M)) .* ...
       exp(-j.*km.*r) ./ r .*...
       (ones(3,1)*dS),2);
   
   H(:,1)=km./(4.*pi.*j).*sum((sqrt(Epsr)./Z0.*(aa.*M+bb.*(ones(3,1)*dot(r1,M)).*r1)+...
      cc.*cross(r1,J)).*exp(-j.*km.*r)./r.*(ones(3,1)*dS),2);
   temps=toc.*(length(P)-1);
   heures=floor(temps./3600);
   temps=mod(temps,3600);
   minu=floor(temps./60);
   sec=mod(temps,60);
   
   h=waitbar(0,['Integrale de rayonnement  ' ...
         'duree estimee : ' num2str(heures) ' heures ' ...
         num2str(minu) ' min ' num2str(sec) ' sec']);
   

   for ind = 2:size(P,2)
      r = ones(3,1)*(...
          sqrt((P(1,ind)-Q(1,:)).^2+(P(2,ind)-Q(2,:)).^2+(P(3,ind)-Q(3,:)).^2)...
          );
      r1 = [ P(1,ind)-Q(1,:) ; P(2,ind)-Q(2,:) ; P(3,ind)-Q(3,:)] ./ r;
      aa = 1 - j./(km.*r) - 1./(km.*r).^2;
      bb = -1 + 3.*j./(km.*r) + 3./(km.*r).^2;
      cc = 1 - j./(km.*r);
     
      E(:,ind) = km./(4.*pi.*j) .* sum(...
          (Z0 ./ sqrt(Epsr) * (aa.*J+bb.*(ones(3,1)*dot(r1,J)).*r1) - ...
          cc .* cross(r1,M) ) .* ...
          exp(-j.*km.*r) ./ r .*...
          (ones(3,1)*dS), 2);
     
      H(:,ind)=km./(4.*pi.*j).*sum((sqrt(Epsr)./Z0.*(aa.*M+bb.*(ones(3,1)*dot(r1,M)).*r1)+...
         cc.*cross(r1,J)).*exp(-j.*km.*r)./r.*(ones(3,1)*dS),2);
     
     waitbar(ind./length(P),h);
   end
   close(h);
   
case 2
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Calcul de E valable partout %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
      % Estimation du temps de calcul
   tic
   r=ones(3,1)*(sqrt((P(1,1)-Q(1,:)).^2+(P(2,1)-Q(2,:)).^2+(P(3,1)-Q(3,:)).^2));
   r1=[P(1,1)-Q(1,:);P(2,1)-Q(2,:);P(3,1)-Q(3,:)]./r;
   aa=1-j./(km.*r)-1./(km.*r).^2;
   bb=-1+3.*j./(km.*r)+3./(km.*r).^2;
   cc=1-j./(km.*r);
   E(:,1)=km./(4.*pi.*j).*sum((Z0./sqrt(Epsr)*(aa.*J+bb.*(ones(3,1)*dot(r1,J)).*r1)-...
      cc.*cross(r1,M)).*exp(-j.*km.*r)./r.*(ones(3,1)*dS),2);
   temps=toc.*(length(P)-1);
   heures=floor(temps./3600);
   temps=mod(temps,3600);
   minu=floor(temps./60);
   sec=mod(temps,60);
   
   h=waitbar(0,['Integrale de rayonnement  ' ...
         'duree estimee : ' num2str(heures) ' heures ' ...
         num2str(minu) ' min ' num2str(sec) ' sec']);
   

   for ind=2:size(P,2)
      r=ones(3,1)*(sqrt((P(1,ind)-Q(1,:)).^2+(P(2,ind)-Q(2,:)).^2+(P(3,ind)-Q(3,:)).^2));
      r1=[P(1,ind)-Q(1,:);P(2,ind)-Q(2,:);P(3,ind)-Q(3,:)]./r;
      aa=1-j./(km.*r)-1./(km.*r).^2;
      bb=-1+3.*j./(km.*r)+3./(km.*r).^2;
      cc=1-j./(km.*r);
      E(:,ind)=km./(4.*pi.*j).*sum((Z0./sqrt(Epsr)*(aa.*J+bb.*(ones(3,1)*dot(r1,J)).*r1)-...
         cc.*cross(r1,M)).*exp(-j.*km.*r)./r.*(ones(3,1)*dS),2);
      waitbar(ind./length(P),h);
   end
   H=[];
   close(h);
   
case 3
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Calcul de E champ lointain  %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   R=ones(3,1)*(sqrt(P(1,:).^2+P(2,:).^2+P(3,:).^2));
   R1=P./R;
   
   tic
   aa(:,1)=sum(J.*(ones(3,1)*...
      (exp(-j.*km.*sqrt((P(1,1)-Q(1,:)).^2+...
      (P(2,1)-Q(2,:)).^2+...
      (P(3,1)-Q(3,:)).^2)).*dS)),2);
   bb(:,1)=sum(M.*(ones(3,1)*...
      (exp(-j.*km.*sqrt((P(1,1)-Q(1,:)).^2+...
      (P(2,1)-Q(2,:)).^2+...
      (P(3,1)-Q(3,:)).^2)).*dS)),2);
   
   temps=toc.*(length(P)-1);
   heures=floor(temps./3600);
   temps=mod(temps,3600);
   minu=floor(temps./60);
   sec=mod(temps,60);
   
   h=waitbar(0,['Integrale de rayonnement  ' ...
         'duree estimee : ' num2str(heures) ' heures ' ...
         num2str(minu) ' min ' num2str(sec) ' sec']);
         
   for ind=2:size(P,2)
      aa(:,ind)=sum(J.*(ones(3,1)*...
         (exp(-j.*km.*sqrt((P(1,ind)-Q(1,:)).^2+...
         (P(2,ind)-Q(2,:)).^2+...
         (P(3,ind)-Q(3,:)).^2)).*dS)),2);
      bb(:,ind)=sum(M.*(ones(3,1)*...
         (exp(-j.*km.*sqrt((P(1,ind)-Q(1,:)).^2+...
         (P(2,ind)-Q(2,:)).^2+...
         (P(3,ind)-Q(3,:)).^2)).*dS)),2);
      waitbar(ind./length(P),h);
   end
   close(h);
   E=-km./(4.*pi.*j.*R).*cross(R1,Z0./sqrt(Epsr).*cross(R1,aa)+bb);
   H=[];
end


