% Fonction  :   RayonnementKottler                               
% Date      :   13/01/03                                  
% Auteur    :   A .Chabory                                
% Révision  :   10/12/04 - S. Bolioli                                  
% Etat      :   OK                                                        
%                                                         
% Rayonnement de J et M a partir des formules de Kottler  
%                                                         
% [E,H]=RayonnementKottler(P,Q,dS,J,M,Freq,Epsr,Option)
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
%   TPs  :   Temps de calcul, en s               (1)
%                                                         

function [E, H, Tps] = RayonnementKottler(P, Q, dS, J, M, Freq, Epsr, Option)

%%%%%%%%% tests des arguments
NP = size(P, 2);
NQ = size(Q, 2);

if (not(isequal(size(P), [3 NP])))
  error('La Taille du vecteur P est erronée');
elseif (not(isequal(size(Q), [3 NQ])))
  error('La Taille du vecteur Q est erronée');
elseif (not(isequal(size(dS), [1 NQ])))
  error('La taille du vecteur dS est erronée');
elseif (not(isequal(size(J), [3 NQ])))
  error('La taille du vecteur J est erronée');
elseif (not(isequal(size(M), [3 NQ])))
  error('La taille du vecteur M est erronée');
end

%%%%%%%%% affichage des informations du calcul
disp('Calcul de champs rayonnés par des courants J et M par les formules de Kottler');
disp(['  ', num2str(size(P)), ' points d''observation P']);
disp(['  ', num2str(size(Q)), ' points courants sur la surface Q']);
disp(['  ', num2str(size(dS)), ' points d''élements de surface dS au point Q']);
disp(['  ', num2str(size(J)), ' points de courant électrique en Q']);
disp(['  ', num2str(size(M)), ' points de courant magnétique en Q']);
disp(['  ', 'Fréquence f=', num2str(Freq), ' Hz']);
disp(['  ', 'Permittivité relative du milieu : Epsr=', num2str(Epsr)]);
switch (Option)
    case 1
        disp(['   ', 'Calcul de E et H']);
    case 2
        disp(['   ', 'Calcul de E seul']);
    case 3
        disp(['   ', 'Calcul de E seul en champ lointain']);
    otherwise
        error('Mauvaise argument : option.');
end


%%%%%%%%%%% Constantes
lam0=2.997925e8./Freq;
lamm=lam0./sqrt(Epsr);
k0=2.*pi./lam0;
km=2.*pi./lamm;
Z0=120.*pi;

%  % Seuillage des valeurs non significatives
%  JdB=10.*log10(sum(abs(J).^2,1));
%  JdB=JdB-max(JdB);
%  MdB=10.*log10(sum(abs(M).^2,1));
%  MdB=MdB-max(MdB);
%  Mq=find(JdB>-80 | MdB>-80);
%  Q=Q(:,Mq);
%  dS=dS(:,Mq);
%  J=J(:,Mq);
%  M=M(:,Mq);

% Progress Bar
PB = progress('init', 'Calcul de l''intégrale de rayonnement');   

switch Option
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calcul de E et H valable partout %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  case 1
    if size(P,2) > 10
      tic
      for ind=1:10
        PB=progress(PB, ind./length(P));
        %,h, ...
        %  ['Integrale de rayonnement                                 ' ...
        %    num2str(ind) '/' num2str(size(P,2))]);
        QP = P(:,ind)*ones(1,size(Q,2))-Q;
        r=ones(3,1)*sqrt(sum(QP.^2,1));
        r1=QP./r;
        kmr = km*r;
        kmr_1 = 1./kmr;
        kmr_2 = kmr_1.^2;
        aa=1-j*kmr_1-kmr_2;
        bb=-1+3*j*kmr_1+3*kmr_2;
        cc=1-j*kmr_1;
        E(:,ind)=km./(4.*pi.*j).*sum((Z0./sqrt(Epsr)*(aa.*J+bb.*(ones(3,1)*dot(r1,J)).*r1)+...
          cc.*cross(r1,-M)).*exp(-j.*km.*r)./r.*(ones(3,1)*dS),2);
        H(:,ind)=-km./(4.*pi.*j).*sum((sqrt(Epsr)./Z0.*(-aa.*M+bb.*(ones(3,1)*dot(r1,-M)).*r1)-...
          cc.*cross(r1,J)).*exp(-j.*km.*r)./r.*(ones(3,1)*dS),2);
      end
      Tps = toc;
      Ideb = 11;
      TempsEstime = Tps/10*size(P,2);
      heures=floor(TempsEstime/3600);
      TempsEstime=mod(TempsEstime,3600);
      minu=floor(TempsEstime/60);
      sec=mod(TempsEstime,60);
      TempsEstime = [num2str(heures) 'h' num2str(minu) 'm' num2str(sec,'%3.0d')];
    else
      Ideb = 1;
      TempsEstime = [];
    end
    
    for ind=Ideb:size(P,2)
      PB = progress(PB, ind/length(P));
      %,h, ...
      %  ['Integrale de rayonnement ... Temps estimé : ' TempsEstime  ...
      %    ' (' num2str(ind) '/' num2str(size(P,2)) ')']);
      QP = P(:,ind)*ones(1,size(Q,2))-Q;
      r=ones(3,1)*sqrt(sum(QP.^2,1));
      r1=QP./r;
      kmr = km*r;
      kmr_1 = 1./kmr;
      kmr_2 = kmr_1.^2;
      aa=1-j*kmr_1-kmr_2;
      bb=-1+3*j*kmr_1+3*kmr_2;
      cc=1-j*kmr_1;
      E(:,ind)=km./(4.*pi.*j).*sum((Z0./sqrt(Epsr)*(aa.*J+bb.*(ones(3,1)*dot(r1,J)).*r1)+...
        cc.*cross(r1,-M)).*exp(-j.*km.*r)./r.*(ones(3,1)*dS),2);
      H(:,ind)=-km./(4.*pi.*j).*sum((sqrt(Epsr)./Z0.*(-aa.*M+bb.*(ones(3,1)*dot(r1,-M)).*r1)-...
        cc.*cross(r1,J)).*exp(-j.*km.*r)./r.*(ones(3,1)*dS),2);
    end
%    close(h);
    Tps = toc;
    
  case 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calcul de E valable partout %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %h=waitbar(0,'Integrale de rayonnement                                           ');
    tic
    if size(P,2) > 10
        tic
        for ind=1:10
            PB = progress(PB, ind./length(P));
            %         waitbar(ind./length(P),h, ...
            %           ['Integrale de rayonnement                                 ' ...
            %             num2str(ind) '/' num2str(size(P,2))]);
            QP = P(:,ind)*ones(1,size(Q,2))-Q;
            r=ones(3,1)*sqrt(sum(QP.^2,1));
            r1=QP./r;
            kmr = km*r;
            kmr_1 = 1./kmr;
            kmr_2 = kmr_1.^2;
            aa=1-j*kmr_1-kmr_2;
            bb=-1+3*j*kmr_1+3*kmr_2;
            cc=1-j*kmr_1;
            E(:,ind)=km./(4.*pi.*j).*sum((Z0./sqrt(Epsr)*(aa.*J+bb.*(ones(3,1)*dot(r1,J)).*r1)+...
                cc.*cross(r1,-M)).*exp(-j.*km.*r)./r.*(ones(3,1)*dS),2);
        end
        Tps = toc;
        Ideb = 11;
        TempsEstime = Tps/10*size(P,2);
        heures=floor(TempsEstime/3600);
        TempsEstime=mod(TempsEstime,3600);
        minu=floor(TempsEstime/60);
        sec=mod(TempsEstime,60);
        TempsEstime = [num2str(heures) 'h' num2str(minu) 'm' num2str(sec,'%3.0d')];
    else
        Ideb = 1;
        TempsEstime = [];
    end

    for ind=Ideb:size(P,2)
        PB = progress(PB, ind./length(P));
        %        waitbar(ind/length(P),h, ...
        %         ['Integrale de rayonnement ... Temps estimé : ' TempsEstime  ...
        %           ' (' num2str(ind) '/' num2str(size(P,2)) ')']);
        QP = P(:,ind)*ones(1,size(Q,2))-Q;
        r=ones(3,1)*sqrt(sum(QP.^2,1));
        r1=QP./r;
        kmr = km*r;
        kmr_1 = 1./kmr;
        kmr_2 = kmr_1.^2;
        aa=1-j*kmr_1-kmr_2;
        bb=-1+3*j*kmr_1+3*kmr_2;
        cc=1-j*kmr_1;
        E(:,ind)=km./(4.*pi.*j).*sum((Z0./sqrt(Epsr)*(aa.*J+bb.*(ones(3,1)*dot(r1,J)).*r1)+...
            cc.*cross(r1,-M)).*exp(-j.*km.*r)./r.*(ones(3,1)*dS),2);
    end
    H = [];
    %     close(h);
    Tps = toc;
    
    case 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calcul de E champ lointain  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    R=ones(3,1)*(sqrt(P(1,:).^2+P(2,:).^2+P(3,:).^2));
    R1=P./R;
    
%    h=waitbar(0,'Integrale de rayonnement                            ');
    if size(P,2) > 10
      tic
      for ind=1:10
          PB = progress(PB, ind./length(P));
%         waitbar(ind./length(P),h, ...
%           ['Integrale de rayonnement                                 ' ...
%             num2str(ind) '/' num2str(size(P,2))]);
        QP = P(:,ind)*ones(1,size(Q,2))-Q;
        r=sqrt(sum(QP.^2,1));
        U = ones(3,1)*(exp(-j.*km.*r).*dS);
        aa(:,ind)=sum(J.*U,2);
        bb(:,ind)=sum(M.*U,2);
      end
      Tps = toc;
      Ideb = 11;
      TempsEstime = Tps/10*size(P,2);
      heures=floor(TempsEstime/3600);
      TempsEstime=mod(TempsEstime,3600);
      minu=floor(TempsEstime/60);
      sec=mod(TempsEstime,60);
      TempsEstime = [num2str(heures) 'h' num2str(minu) 'm' num2str(sec,'%3.0d')];
    else
      Ideb = 1;
      TempsEstime = [];
    end
    for ind=Ideb:size(P,2)
        PB = progress(PB, ind./length(P));
%       waitbar(ind/length(P),h, ...
%         ['Integrale de rayonnement ... Temps estimé : ' TempsEstime  ...
%           ' (' num2str(ind) '/' num2str(size(P,2)) ')']);
      QP = P(:,ind)*ones(1,size(Q,2))-Q;
      r=sqrt(sum(QP.^2,1));
      U = ones(3,1)*(exp(-j.*km.*r).*dS);
      aa(:,ind)=sum(J.*U,2);
      bb(:,ind)=sum(M.*U,2);
    end
%    close(h);
    E=-km./(4.*pi.*j.*R).*cross(R1,Z0./sqrt(Epsr).*cross(R1,aa)+bb);
    H=[];
end


