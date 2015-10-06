function [E,H, Epc1, Epc2, Epc3] = EH_FG_OP_Plaque_conductrice_uniforme(Champ, Plaque, P, f, epsr)
% Calcul du champ EM rayonne par une Plaque conductrice 
% dans l'hypothese de l'optique physique
% par un faisceau gaussien
% 
% Les effets de la finitude sont pris en compte
% Ainsi que les effets "de coin" de maniere cette fois uniforme
% 
% PARAMETRES
%  Champ : structure Champ faisceau gaussien, contenant :
%           {Centre, ex, ez, Coefficients, Courbure}
%  Plaque: structure Plaque, contenant :
%           {centre, X1, X2, Y1, Y2, n, cart, dx, dy}
%  P : points d'observations dans le repere absolu (3 x NbP)
%  f : frequence
%  epsr : epsilon relatif du milieu de propagation
% 
% RETOURNE
%  E : Champ Electrique aux points P (3 x NbP)
%  H : Champ Magnetique aux points P (3 x NbP)
%  Epc1 : Champ issu des points critiques du premier ordre
%  Epc2 : Champ issu des points critiques du second ordre
%  Epc3 : Champ issu des points critiques du troisieme ordre
global Axx_l Axy_l Ayy_l Bx_l By_l C_l;


% ===========================================
%      CONSTANTES et PARAMETRES
% ===========================================
     c = 2.997925e8;
    Z0 = 120*pi;
lambda = c/f;
     k = 2*pi/lambda;
   NbP = length(P);


X1 = max(Plaque.Sommets(1,:));
X2 = min(Plaque.Sommets(1,:));
Y1 = max(Plaque.Sommets(2,:));
Y2 = min(Plaque.Sommets(2,:));

%  Calcul du point d'intersection du faisceau gaussien
%  avec le plan. Ce point correspond à l'origine du repere local
%  dans lequel on va faire les calculs.
%  
%  Distance à parcourir entre le centre du repere faisceau
%  et le plan :
d =  (Plaque.Centre(3) - Champ.Centre(3))./Champ.ez(3);
%  coordonnees du point d'intersection :
O_l = Champ.Centre + d*Champ.ez;
%  Changements de variable
 P_l = P  - O_l*ones(1,NbP);
X1_l = X1 - O_l(1);
X2_l = X2 - O_l(1);
Y1_l = Y1 - O_l(2);
Y2_l = Y2 - O_l(2);
   R = sqrt(sum(P.^2));
 R_l = sqrt(sum(P_l.^2)); 

%  % matrice de courbure en z_fc = 0
%  Q0    = Champ.Courbure;
%  detQ0 = Q0(1)*Q0(2)-Q0(3)^2;  
%  invQ0 = [Q0(2); Q0(1); -Q0(3)] ./ (ones(3,1)*detQ0);
%  % matrice de courbure propagée en zi=d, soit sur le plan metallique
%  % Hypothese : elle est consideree (a tord) comme constante sur toute
%  % la plaque
%      inv_Qp = [invQ0(1)+d; invQ0(2)+d; invQ0(3)];
%  det_inv_Qp = inv_Qp(1,:).*inv_Qp(2,:) - inv_Qp(3,:).^2;
%          Qp = [inv_Qp(2,:);inv_Qp(1,:);-inv_Qp(3,:)]...
%                ./(ones(3,1)*det_inv_Qp);

[Q0, detQ0] = subfun_Qz(0, Champ);  
[Qp, det_Qp] = subfun_Qz(d, Champ);  

% Matrice de passage pour passer du repere faisceau
% au repere absolu
R_Rfsc2Rabs = [Champ.ex cross(Champ.ez,Champ.ex) Champ.ez];
% raccourcis de notation de la matrice de rotation
% Matrice de passage pour passer du repere 
R_Rabs2Rfsc = transpose(R_Rfsc2Rabs);
r11 = R_Rabs2Rfsc(1,1);
r12 = R_Rabs2Rfsc(1,2);
r21 = R_Rabs2Rfsc(2,1);
r22 = R_Rabs2Rfsc(2,2);
r31 = R_Rabs2Rfsc(3,1);
r32 = R_Rabs2Rfsc(3,2);

% Coordonnees de P dans le repere fsc
P_l_fsc = R_Rabs2Rfsc*(P_l-Champ.Centre*ones(1,NbP)-O_l*ones(1,NbP));

% position du centre de la Plaque dans le repere faisceau
%  PlaqueCentre_l_fsc = R_Rabs2Rfsc*(Plaque.Centre - O_l - Champ.Centre);





% centre du faisceau exprimé dans le repere faisceau
r_l_i_Oi = R_Rabs2Rfsc*(Champ.Centre - O_l);
x_l_i_Oi = r_l_i_Oi(1);
y_l_i_Oi = r_l_i_Oi(2);
z_l_i_Oi = r_l_i_Oi(3);

ri_Oi = R_Rabs2Rfsc*Champ.Centre;
xi_Oi = ri_Oi(1);
yi_Oi = ri_Oi(2);
zi_Oi = ri_Oi(3);

M_fsc = R_Rfsc2Rabs;

% elements de la notation matricielle de la phase
Axx_l = j*(1./R_l - P_l(1,:).^2./R_l.^3 ...
  + Qp(1).*r11.^2 ...
  + Qp(2).*r21.^2 ...
  + 2*Qp(3).*r11.*r21);

Ayy_l =  j*(1./R_l - P_l(2,:).^2./R_l.^3 ...
  + Qp(2).*r22.^2 ...
  + Qp(1).*r12.^2 ...
  + 2*Qp(3).*r12.*r22);

Axy_l = j*( - P_l(1,:).*P_l(2,:)./R_l.^3 ...
  + Qp(3).*(r11.*r22 + r12.*r21) ...
  + Qp(2).*r21.*r22 ...
  + Qp(1).*r11.*r12)

Bx_l = j*(- r31 ...
  + P_l(1,:)./R_l ...
  + Qp(3).*(x_l_i_Oi.*r21 + r11.*y_l_i_Oi) ...
  + Qp(2).*r21.*y_l_i_Oi ...
  + Qp(1).*r11.*x_l_i_Oi);

By_l = j*(- r32 ...
  + P_l(2,:)./R_l ...
  + Qp(3).*(x_l_i_Oi.*r22 + r12.*y_l_i_Oi) ...
  + Qp(2).*r22.*y_l_i_Oi ...
  + Qp(1).*r12.*x_l_i_Oi);

C_l = j*(R_l ...
  - z_l_i_Oi ...
  +1/2*Qp(1).*x_l_i_Oi.^2 ...
  +1/2*Qp(2).*y_l_i_Oi.^2 ...
  +Qp(3).*x_l_i_Oi.*y_l_i_Oi);

detA_l = Axx_l.*Ayy_l - Axy_l.^2;
iA = [Ayy_l;Axx_l;-Axy_l]./(ones(3,1)*detA_l);



% ======================================================
% Calcul du point col (point critique du premier ordre)
% ======================================================
%  point col dans le repere absolu modifie
xQ_l_s = iA(1,:).*Bx_l + iA(3,:).*By_l;
yQ_l_s = iA(3,:).*Bx_l + iA(2,:).*By_l;
rQ_l_s = [xQ_l_s; yQ_l_s; zeros(1,NbP)];
%  point col dans le repere absolu 
xQ_s = xQ_l_s + O_l(1);
yQ_s = yQ_l_s + O_l(2);
rQ_s = [xQ_s; yQ_s; Plaque.Centre(3)*ones(1,NbP)];

% point col dans le repere fsc
rQ_fsc_s = R_Rabs2Rfsc*(rQ_s  - Champ.Centre*ones(1,NbP) );

%  rQ_l_fsc_s = R_Rabs2Rfsc*(rQ_l_s  - Champ.Centre*ones(1,NbP) );
% matrice de courbure evaluee au point col dans le repere faisceau
%  iQz_fsc_s = [invQ0(1)+rQ_fsc_s(3,:); ...
%             invQ0(2)+rQ_fsc_s(3,:); ...
%             invQ0(3)*ones(1,NbP)];
%  det_iQz_fsc_s = iQz_fsc_s(1,:).*iQz_fsc_s(2,:) - iQz_fsc_s(3,:).^2;
%  Qz_fsc_s = [iQz_fsc_s(2,:); ...
%              iQz_fsc_s(1,:); ...
%             -iQz_fsc_s(3,:)]./(ones(3,1)*det_iQz_fsc_s);
%  det_Qz_fsc_s = Qz_fsc_s(1,:).*Qz_fsc_s(2,:) - Qz_fsc_s(3,:).^2;

% vecteur distance au point col
    r_s = P_l - rQ_l_s;
sum_r_s = sqrt(sum(abs(r_s).^2));
   rr_s = r_s./(ones(3,1)*sum_r_s);

% fonction de phase au point col
%  g_l = 1/2*(xQ_l_s.^2.*Axx_l + 2*xQ_l_s.*yQ_l_s.*Axy_l + yQ_l_s.^2.*Ayy_l) ...
%      - Bx_l.*xQ_l_s - By_l.*yQ_l_s ...
%      + C_l;

%  matrice de courbure evaluee au point critique du premier ordre
[Qz_s, det_Qz_fsc_s] = subfun_Qz(rQ_fsc_s(3,:), Champ);  
%  fonction de phase evaluee au point critique du premier ordre
    g_s_l = subfun_g_l(xQ_l_s, yQ_l_s);
%  fonction de polarisation dans le repere faisceau incident
%  evaluee au point critique du premier ordre
h_fsc_s = subfun_h(rQ_fsc_s, Champ);
%  fonction de polarisation dans le repere absolu
%  evaluee au point critique du premier ordre
h_s = R_Rfsc2Rabs*h_fsc_s;

% fonction d'amplitude au point col
ampl_s =  j*sqrt(epsr) ...
        .* sqrt(det_Qz_fsc_s./detQ0./detA_l) ...
        ./ sum_r_s;

% application de la technique du point col
saddle_point =  ampl_s .* exp(-k*g_s_l);


%  h_l_s_x = -Champ.Coefficients(2).*ones(1,NbP);
%  h_l_s_y =  Champ.Coefficients(1).*ones(1,NbP);
%       Qx = Qz_fsc_s(1,:).*rQ_fsc_s(1,:) + Qz_fsc_s(3,:).*rQ_fsc_s(2,:);
%       Qy = Qz_fsc_s(3,:).*rQ_fsc_s(1,:) + Qz_fsc_s(2,:).*rQ_fsc_s(2,:);
%  h_l_s_z = -Champ.Coefficients(1).*Qy + Champ.Coefficients(2).*Qx;
%    h_l_s = [h_l_s_x;h_l_s_y;h_l_s_z];
%    h_s   = M_fsc*h_l_s;



% =====================================================
% point critique de premier ordre (point col)
% =====================================================
% on applique le point col que lorsque le point appartient au domaine d'integration
%  U1= (real(rQ_s(1,:)) > X2) ...
%    & (real(rQ_s(1,:)) < X1) ...
%    & (real(rQ_s(2,:)) > Y2) ...
%    & (real(rQ_s(2,:)) < Y1);
U1 = ones(1,NbP);

%  Si le centre d'impact du faisceau correspond touche l'une
%  des bordures de la plaque, il faut prendre la demi contribution
%  du point critique du premier ordre (car dans ce cas le point col
%  correspond à une des limites d'integration)
%  U2=1-1/2*(O_l(1) == X1 ...
%         || O_l(1) == X2 ...
%         || O_l(2) == Y2 ...
%         || O_l(2) == Y2);
U2 = ones(1,NbP);


Epc1 = cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), h_s))) .* ...
      (ones(3,1)*(U1.*U2.*saddle_point));


  
% =====================================================
% Contribution des points critiques du second ordre
% =====================================================
% points critiques du second ordre
x_Y1_l = (Bx_l-Axy_l.*Y1_l)./Axx_l;
x_Y2_l = (Bx_l-Axy_l.*Y2_l)./Axx_l;
y_X1_l = (By_l-Axy_l.*X1_l)./Ayy_l;
y_X2_l = (By_l-Axy_l.*X2_l)./Ayy_l;

x_Y1 = x_Y1_l + O_l(1);
x_Y2 = x_Y2_l + O_l(1);
y_X1 = y_X1_l + O_l(2);
y_X2 = y_X2_l + O_l(2);

%  %  Trace des points critiques
%  figure(1)
%    hold on
%    plot3(real(rQ_s(3,:)), ...
%          real(-rQ_s(2,:)), ...
%          real(rQ_s(1,:)), 'r.')  
%  
%    plot3(X1*ones(size(y_X1))/lambda, real(y_X1)/lambda, zeros(size(y_X1)), 'g.');
%    plot3(X2*ones(size(y_X1))/lambda, real(y_X2)/lambda, zeros(size(y_X1)), 'g.');
%    plot3(real(x_Y1)/lambda, Y1*ones(size(x_Y1))/lambda, zeros(size(y_X1)), 'b.');
%    plot3(real(x_Y2)/lambda, Y2*ones(size(x_Y2))/lambda, zeros(size(y_X1)), 'b.')
%    hold off


% vecteur distance aux points critiques du second ordre

   rQ_X1 = [X1*ones(1,NbP); y_X1;Plaque.Centre(3)*ones(1,NbP)]; 
 rQ_X1_l = [X1_l*ones(1,NbP); y_X1_l;Plaque.Centre(3)*ones(1,NbP)]; 
    r_X1 = P_l - rQ_X1_l;
sum_r_X1 = sqrt(sum(abs(r_X1.^2)));
   rr_X1 = r_X1./(ones(3,1)*sum_r_X1);

   rQ_X2 = [X2*ones(1,NbP); y_X2;Plaque.Centre(3)*ones(1,NbP)]; 
 rQ_X2_l = [X2_l*ones(1,NbP); y_X2_l;Plaque.Centre(3)*ones(1,NbP)]; 
    r_X2 = P_l - rQ_X2_l;
sum_r_X2 = sqrt(sum(abs(r_X2.^2)));
   rr_X2 = r_X2./(ones(3,1)*sum_r_X2);

   rQ_Y1 = [x_Y1; Y1*ones(1,NbP); Plaque.Centre(3)*ones(1,NbP)]; 
 rQ_Y1_l = [x_Y1_l; Y1_l*ones(1,NbP); Plaque.Centre(3)*ones(1,NbP)]; 
    r_Y1 = P_l - rQ_Y1_l;
sum_r_Y1 = sqrt(sum(abs(r_Y1.^2)));
   rr_Y1 = r_Y1./(ones(3,1)*sum_r_Y1);

   rQ_Y2 = [x_Y2; Y2*ones(1,NbP); Plaque.Centre(3)*ones(1,NbP)]; 
 rQ_Y2_l = [x_Y2_l; Y2_l*ones(1,NbP); Plaque.Centre(3)*ones(1,NbP)]; 
    r_Y2 = P_l - rQ_Y2_l;
sum_r_Y2 = sqrt(sum(abs(r_Y2).^2));
   rr_Y2 = r_Y2./(ones(3,1)*sum_r_Y2);

%  figure(1)
%    hold on
%    plot3(real(rQ_X1(3,:)), -real(rQ_X1(2,:)), real(rQ_X1(1,:)), 'g.');
%    plot3(real(rQ_X2(3,:)), -real(rQ_X2(2,:)), real(rQ_X2(1,:)), 'g.');
%    plot3(real(rQ_Y1(3,:)), -real(rQ_Y1(2,:)), real(rQ_Y1(1,:)), 'b.');
%    plot3(real(rQ_Y2(3,:)), -real(rQ_Y2(2,:)), real(rQ_Y2(1,:)), 'b.');
%    hold off

% point critiques du second ordre dans le repere du faisceau incident
rQ_fsc_X1 = R_Rabs2Rfsc*(rQ_X1  - Champ.Centre*ones(1,NbP));
rQ_fsc_X2 = R_Rabs2Rfsc*(rQ_X2  - Champ.Centre*ones(1,NbP));
rQ_fsc_Y1 = R_Rabs2Rfsc*(rQ_Y1  - Champ.Centre*ones(1,NbP));
rQ_fsc_Y2 = R_Rabs2Rfsc*(rQ_Y2  - Champ.Centre*ones(1,NbP));

%  xQ_fsc_X1 = r11.*rQ_X1(1,:) + r12.*rQ_X1(2,:) - xi_Oi;
%  yQ_fsc_X1 = r21.*rQ_X1(1,:) + r22.*rQ_X1(2,:) - yi_Oi;
%  zQ_fsc_X1 = r31.*rQ_X1(1,:) + r32.*rQ_X1(2,:) - zi_Oi;
%  
%  xQ_fsc_X2 = r11.*rQ_X2(1,:) + r12.*rQ_X2(2,:) - xi_Oi;
%  yQ_fsc_X2 = r21.*rQ_X2(1,:) + r22.*rQ_X2(2,:) - yi_Oi;
%  zQ_fsc_X2 = r31.*rQ_X2(1,:) + r32.*rQ_X2(2,:) - zi_Oi;
%  
%  xQ_fsc_Y1 = r11.*rQ_Y1(1,:) + r12.*rQ_Y1(2,:) - xi_Oi;
%  yQ_fsc_Y1 = r21.*rQ_Y1(1,:) + r22.*rQ_Y1(2,:) - yi_Oi;
%  zQ_fsc_Y1 = r31.*rQ_Y1(1,:) + r32.*rQ_Y1(2,:) - zi_Oi;
%  
%  xQ_fsc_Y2 = r11.*rQ_Y2(1,:) + r12.*rQ_Y2(2,:) - xi_Oi;
%  yQ_fsc_Y2 = r21.*rQ_Y2(1,:) + r22.*rQ_Y2(2,:) - yi_Oi;
%  zQ_fsc_Y2 = r31.*rQ_Y2(1,:) + r32.*rQ_Y2(2,:) - zi_Oi;



%  % matrice de courbure dans le repere incident
%  iQz_fsc_X1 = [invQ0(1)+rQ_fsc_X1(3,:); ...
%                invQ0(2)+rQ_fsc_X1(3,:); ...
%                invQ0(3)*ones(1,NbP)];
%  det_iQz_fsc_X1 = iQz_fsc_X1(1,:).*iQz_fsc_X1(2,:) - iQz_fsc_X1(3,:).^2;
%  Qz_fsc_X1 = [iQz_fsc_X1(2,:); ...
%               iQz_fsc_X1(1,:); ...
%              -iQz_fsc_X1(3,:)]./(ones(3,1)*det_iQz_fsc_X1);
%  det_Qz_fsc_X1 = Qz_fsc_X1(1,:).*Qz_fsc_X1(2,:) - Qz_fsc_X1(3,:).^2;
%  
%  iQz_fsc_X2 = [invQ0(1)+rQ_fsc_X2(3,:); ...
%                invQ0(2)+rQ_fsc_X2(3,:); ...
%                invQ0(3)*ones(1,NbP)];
%  det_iQz_fsc_X2 = iQz_fsc_X2(1,:).*iQz_fsc_X2(2,:) - iQz_fsc_X2(3,:).^2;
%  Qz_fsc_X2 = [iQz_fsc_X2(2,:); ...
%               iQz_fsc_X2(1,:); ...
%              -iQz_fsc_X2(3,:)]./(ones(3,1)*det_iQz_fsc_X2);
%  det_Qz_fsc_X2 = Qz_fsc_X2(1,:).*Qz_fsc_X2(2,:) - Qz_fsc_X2(3,:).^2;
%  
%  iQz_fsc_Y1 = [invQ0(1)+rQ_fsc_Y1(3,:); ...
%                invQ0(2)+rQ_fsc_Y1(3,:); ...
%                invQ0(3)*ones(1,NbP)];
%  det_iQz_fsc_Y1 = iQz_fsc_Y1(1,:).*iQz_fsc_Y1(2,:) - iQz_fsc_Y1(3,:).^2;
%  Qz_fsc_Y1 = [iQz_fsc_Y1(2,:); ...
%            iQz_fsc_Y1(1,:); ...
%           -iQz_fsc_Y1(3,:)]./(ones(3,1)*det_iQz_fsc_Y1);
%  det_Qz_fsc_Y1 = Qz_fsc_Y1(1,:).*Qz_fsc_Y1(2,:) - Qz_fsc_Y1(3,:).^2;
%  
%  iQz_fsc_Y2 = [invQ0(1)+rQ_fsc_Y2(3,:); ...
%                invQ0(2)+rQ_fsc_Y2(3,:); ...
%                invQ0(3)*ones(1,NbP)];
%  det_iQz_fsc_Y2 = iQz_fsc_Y2(1,:).*iQz_fsc_Y2(2,:) - iQz_fsc_Y2(3,:).^2;
%  Qz_fsc_Y2 = [iQz_fsc_Y2(2,:); ...
%            iQz_fsc_Y2(1,:); ...
%           -iQz_fsc_Y2(3,:)]./(ones(3,1)*det_iQz_fsc_Y2);
%  det_Qz_fsc_Y2 = Qz_fsc_Y2(1,:).*Qz_fsc_Y2(2,:) - Qz_fsc_Y2(3,:).^2;

%  Matrices de courbures aux points critiques du second ordre
[Qz_fsc_X1, det_Qz_fsc_X1] = subfun_Qz(rQ_fsc_X1(3,:), Champ);
[Qz_fsc_X2, det_Qz_fsc_X2] = subfun_Qz(rQ_fsc_X2(3,:), Champ);
[Qz_fsc_Y1, det_Qz_fsc_Y1] = subfun_Qz(rQ_fsc_Y1(3,:), Champ);
[Qz_fsc_Y2, det_Qz_fsc_Y2] = subfun_Qz(rQ_fsc_Y2(3,:), Champ);

%  % vecteur polarisation aux points critiques du second ordre
%  h_l_x = -Champ.Coefficients(2).*ones(1,NbP);
%  h_l_y =  Champ.Coefficients(1).*ones(1,NbP);
%  
%  h_l_X1_z = -Champ.Coefficients(1).*(rQ_fsc_X1(1,:).*Qz_fsc_X1(3,:) + rQ_fsc_X1(2,:).*Qz_fsc_X1(2,:)) ...
%             +Champ.Coefficients(2).*(rQ_fsc_X1(1,:).*Qz_fsc_X1(1,:) + rQ_fsc_X1(2,:).*Qz_fsc_X1(3,:));
%  h_l_X1 = [h_l_x;h_l_y;h_l_X1_z];
%  h_X1   = M_fsc*h_l_X1;
%  
%  h_l_X2_z = -Champ.Coefficients(1).*(rQ_fsc_X2(1,:).*Qz_fsc_X2(3,:) + rQ_fsc_X2(2,:).*Qz_fsc_X2(2,:)) ...
%             +Champ.Coefficients(2).*(rQ_fsc_X2(1,:).*Qz_fsc_X2(1,:) + rQ_fsc_X2(2,:).*Qz_fsc_X2(3,:));
%  h_l_X2 = [h_l_x;h_l_y;h_l_X2_z];
%  h_X2   = M_fsc*h_l_X2;
%  
%  h_l_Y1_z = -Champ.Coefficients(1).*(rQ_fsc_Y1(1,:).*Qz_fsc_Y1(3,:) + rQ_fsc_Y1(2,:).*Qz_fsc_Y1(2,:)) ...
%             +Champ.Coefficients(2).*(rQ_fsc_Y1(1,:).*Qz_fsc_Y1(1,:) + rQ_fsc_Y1(2,:).*Qz_fsc_Y1(3,:));
%  h_l_Y1 = [h_l_x;h_l_y;h_l_Y1_z];
%  h_Y1   = M_fsc*h_l_Y1;
%  
%  h_l_Y2_z = -Champ.Coefficients(1).*(rQ_fsc_Y2(1,:).*Qz_fsc_Y2(3,:) + rQ_fsc_Y2(2,:).*Qz_fsc_Y2(2,:)) ...
%             +Champ.Coefficients(2).*(rQ_fsc_Y2(1,:).*Qz_fsc_Y2(1,:) + rQ_fsc_Y2(2,:).*Qz_fsc_Y2(3,:));
%  h_l_Y2 = [h_l_x;h_l_y;h_l_Y2_z];
%  h_Y2   = M_fsc*h_l_Y2;

%  Fonctions de polarisations aux points critiques du second ordre
h_X1 = M_fsc*subfun_h(rQ_fsc_X1, Champ);
h_X2 = M_fsc*subfun_h(rQ_fsc_X2, Champ);
h_Y1 = M_fsc*subfun_h(rQ_fsc_Y1, Champ);
h_Y2 = M_fsc*subfun_h(rQ_fsc_Y2, Champ);

%  % fonction de phase aux points critiques
%  g_X1 = 1/2*(rQ_X1_l(1,:).^2.*Axx_l ...
%              + 2*rQ_X1_l(1,:).*rQ_X1_l(2,:).*Axy_l ...
%              + rQ_X1_l(2,:).^2.*Ayy_l) ...
%      - Bx_l.*rQ_X1_l(1,:) - By_l.*rQ_X1_l(2,:) ...
%      + C_l;
%  
%  g_X2 = 1/2*(rQ_X2_l(1,:).^2.*Axx_l ...
%              + 2*rQ_X2_l(1,:).*rQ_X2_l(2,:).*Axy_l ...
%              + rQ_X2_l(2,:).^2.*Ayy_l) ...
%      - Bx_l.*rQ_X2_l(1,:) - By_l.*rQ_X2_l(2,:) ...
%      + C_l ;
%  
%  g_Y1 = 1/2*(rQ_Y1_l(1,:).^2.*Axx_l ...
%              + 2*rQ_Y1_l(1,:).*rQ_Y1_l(2,:).*Axy_l ...
%              + rQ_Y1_l(2,:).^2.*Ayy_l) ...
%      - Bx_l.*rQ_Y1_l(1,:) - By_l.*rQ_Y1_l(2,:) ...
%      + C_l;
%  
%  g_Y2 = 1/2*(rQ_Y2_l(1,:).^2.*Axx_l ...
%              + 2*rQ_Y2_l(1,:).*rQ_Y2_l(2,:).*Axy_l ...
%              + rQ_Y2_l(2,:).^2.*Ayy_l) ...
%      - Bx_l.*rQ_Y2_l(1,:) - By_l.*rQ_Y2_l(2,:) ...
%      + C_l;

%  Fonctions de phase aux points critiques du second ordre
g_X1 = subfun_g_l(rQ_X1_l(1,:), rQ_X1_l(2,:));
g_X2 = subfun_g_l(rQ_X2_l(1,:), rQ_X2_l(2,:));
g_Y1 = subfun_g_l(rQ_Y1_l(1,:), rQ_Y1_l(2,:));
g_Y2 = subfun_g_l(rQ_Y2_l(1,:), rQ_Y2_l(2,:));

% derivees de la fonction de phase aux points critiques
%  dgdx_X1_l = rQ_X1_l(1,:).*Axx_l + rQ_X1_l(2,:).*Axy_l - Bx_l;
%  dgdx_X2_l = rQ_X2_l(1,:).*Axx_l + rQ_X2_l(2,:).*Axy_l - Bx_l;
%  
%  dgdy_Y1_l = rQ_Y1_l(1,:).*Axy_l + rQ_Y1_l(2,:).*Ayy_l - By_l;
%  dgdy_Y2_l = rQ_Y2_l(1,:).*Axy_l + rQ_Y2_l(2,:).*Ayy_l - By_l;
dgdx_X1_l = (Axx_l-Axy_l.^2./Ayy_l).*X1_l-Bx_l+By_l.*Axy_l./Ayy_l;
dgdx_X2_l = (Axx_l-Axy_l.^2./Ayy_l).*X2_l-Bx_l+By_l.*Axy_l./Ayy_l;

dgdy_Y1_l = (-Axy_l.^2./Axx_l+Ayy_l).*Y1_l+Bx_l.*Axy_l./Axx_l-By_l;
dgdy_Y2_l = (-Axy_l.^2./Axx_l+Ayy_l).*Y2_l+Bx_l.*Axy_l./Axx_l-By_l;

%  % fonction d'amplitude au point col
%  ampl_X1 = 1./sum_r_X1;
%  ampl_X2 = 1./sum_r_X2;
%  ampl_Y1 = 1./sum_r_Y1;
%  ampl_Y2 = 1./sum_r_Y2;
%  
%  % application de la formule 
%  Cte = j*sqrt(epsr)*k/2/pi ;
%  
%  U_X1 = (real(y_X1)>Y2) & (real(y_X1)<Y1);
%  U_X2 = (real(y_X2)>Y2) & (real(y_X2)<Y1);
%  U_Y1 = (real(x_Y1)>X2) & (real(x_Y1)<X1);
%  U_Y2 = (real(x_Y2)>X2) & (real(x_Y2)<X1);
%    
%  %  U_X1 = ones(1,NbP);
%  %  U_X2 = ones(1,NbP);
%  %  U_Y1 = ones(1,NbP);
%  %  U_Y2 = ones(1,NbP);
%  
%  %  ---------- Formulation non uniforme
%  E_X1 = cross(rr_X1, cross(rr_X1, cross(Plaque.N*ones(1,NbP), h_X1))) .* ...
%        (ones(3,1)*(U_X1.*Cte .* ...
%        sqrt(2*pi).*(-1).*... 
%        ampl_X1 .* exp(-k*g_X1) ...
%        ./k^(3/2) ...
%        .*sqrt(det_Qz_fsc_X1./detQ0./Ayy_l)./dgdx_X1_l ));
%  
%  E_X2 = cross(rr_X2, cross(rr_X2, cross(Plaque.N*ones(1,NbP), h_X2))) .* ...
%        (ones(3,1)*(U_X2.*Cte .* ...
%        sqrt(2*pi)*(+1).*...
%        ampl_X2 .* exp(-k*g_X2) ...
%        ./k^(3/2) ...
%        .*sqrt(det_Qz_fsc_X1./detQ0./Ayy_l)./dgdx_X2_l ));
%  
%  E_Y1 = cross(rr_Y1, cross(rr_Y1, cross(Plaque.N*ones(1,NbP), h_Y1))) .* ...
%        (ones(3,1)*(U_Y1.*Cte .* ...
%        sqrt(2*pi)*(-1).*...
%        ampl_Y1 .* exp(-k*g_Y1) ...
%        ./k^(3/2) ...
%        .*sqrt(det_Qz_fsc_X1./detQ0./Axx_l)./dgdy_Y1_l ));
%  
%  E_Y2 = cross(rr_Y2, cross(rr_Y2, cross(Plaque.N*ones(1,NbP), h_Y2))) .* ...
%        (ones(3,1)*(U_Y2.*Cte .* ...
%        sqrt(2*pi)*(+1).*...
%        ampl_Y2 .* exp(-k*g_Y2) ...
%        ./k^(3/2) ...
%        .*sqrt(det_Qz_fsc_X1./detQ0./Axx_l)./dgdy_Y2_l ));


%  Formulations uniformes
Cte = j*sqrt(epsr)*k/2/pi ;

hh_s_X  = sqrt(2./(Axx_l-Axy_l.^2./Ayy_l));
hh_s_Y  = sqrt(2./(-Axy_l.^2./Axx_l+Ayy_l));


  %  Le signe de la racine carrée est défini tel que l'on doit avoir
  %  h_x -> h_s lorsque s -> 0, cad lq 
  %  premiere version ()
%    s_X1 = sign2(O_l(1)).*sign2(real(hh_s_X./(X1_l-xQ_l_s))).*sqrt(g_X1 - g_s_l);
%    s_X2 = - sign2(O_l(1)).*sign2(real(hh_s_X./(X2_l-xQ_l_s))).*sqrt(g_X2 - g_s_l);
%    s_Y1 = sign2(O_l(2)).*sign2(real(hh_s_Y./(Y1_l-yQ_l_s))).*sqrt(g_Y1 - g_s_l);
%    s_Y2 = - sign2(O_l(2)).*sign2(real(hh_s_Y./(Y2_l-yQ_l_s))).*sqrt(g_Y2 - g_s_l);

  s_X1 = sign(real(O_l(1).*hh_s_X./(X1_l-xQ_l_s))).*sqrt(g_X1 - g_s_l);
  s_X2 = - sign(real(O_l(1).*hh_s_X./(X2_l-xQ_l_s))).*sqrt(g_X2 - g_s_l);
  s_Y1 = sign(real(O_l(2).*hh_s_Y./(Y1_l-yQ_l_s))).*sqrt(g_Y1 - g_s_l);
  s_Y2 = - sign(real(O_l(2).*hh_s_Y./(Y2_l-yQ_l_s))).*sqrt(g_Y2 - g_s_l);

%  Lorsque le FG est a gde distance de la plaque :
%  * changer les signes de s_X/Y en +sign2...
%  * changer les signes - dans les expressions de E_X/Y en +
%  * supprimer le Epc1 a la fin

%  %  Moins bon que le premier
%  s_X1 = sign2(O_l(1)).*sqrt((g_X1 - g_s_l)).*sign(real((X1-xQ_s))) ;
%  s_X2 = -sign2(O_l(1)).*sqrt((g_X2 - g_s_l)).*sign(real((X2-xQ_s))) ;
%  s_Y1 = sign2(O_l(2)).*sqrt((g_Y1 - g_s_l)).*sign(real((Y1-yQ_s))) ;
%  s_Y2 = -sign2(O_l(2)).*sqrt((g_Y2 - g_s_l)).*sign(real((Y2-yQ_s))) ;

%  s_Y1 = sqrt((g_Y1 - g_s_l));
%  s_Y1(Y1-yQ_s>0) = -s_Y1(Y1-yQ_s>0);
%  s_Y2 = sqrt((g_Y2 - g_s_l));
%  s_Y2(Y2-yQ_s<0) = -s_Y2(Y2-yQ_s<0);

%  hh_X1 = 2*s_X1./dgdx_X1_l;
%  hh_s  = sqrt(2./(Axx_l-Axy_l.^2./Ayy_l));
%  
%  f_X1 = cross(rr_X1, cross(rr_X1, cross(Plaque.N*ones(1,NbP), h_X1))) .* ...
%        (ones(3,1)*(Cte./ sum_r_X1.* sqrt(det_Qz_fsc_X1./detQ0)));
%  f_s = cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), h_s))) .* ...
%        (ones(3,1)*(Cte./ sum_r_s.* sqrt(det_Qz_fsc_s./detQ0)));
%  
%  E_X1 =  ...
%        f_s.*(ones(3,1)*(exp(-k*g_s_l).*sqrt(2*pi/k./Ayy_l)/sqrt(k).*hh_s.*subfun_Transfun(sqrt(k).*s_X1))) ...
%         +(f_X1.*(ones(3,1)*(hh_X1.*sqrt(2*pi/k./Ayy_l))) ...
%          - f_s.*(ones(3,1)*(hh_s.*sqrt(2*pi/k./Ayy_l))))./(ones(3,1)*(2*k*s_X1));
%  
%  E_X1 = -E_X1;

E_X1 = ...
      Epc1.*sqrt(1/pi).*(ones(3,1)*subfun_Transfun(sqrt(k).*s_X1)) ...
      + ...
     cross(rr_X1, cross(rr_X1, cross(Plaque.N*ones(1,NbP), h_X1))) .* ...
      (ones(3,1)*(Cte*sqrt(2*pi)/k^(3/2)  ...
      ./ sum_r_X1 ./dgdx_X1_l .* exp(-k*g_X1)...
      .* sqrt(det_Qz_fsc_X1./detQ0./Ayy_l))) ...
      - ...
      cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), h_s))) .* ...
      (ones(3,1)*(Cte.*sqrt(pi)/k^(3/2) ...
      ./ sum_r_s  ./s_X1.* exp(-k*g_X1) ...
      .* sqrt(det_Qz_fsc_s./detQ0./detA_l)));
    
E_X1 = -E_X1;

E_X2 = ...
      Epc1.*sqrt(1/pi).*(ones(3,1)*(subfun_Transfun(sqrt(k).*s_X2))) ...
      + ...
     cross(rr_X2, cross(rr_X2, cross(Plaque.N*ones(1,NbP), h_X2))) .* ...
      (ones(3,1)*(Cte*sqrt(2*pi)/k^(3/2)  ...
      ./ sum_r_X2 ./dgdx_X2_l .* exp(-k*g_X2)...
      .* sqrt(det_Qz_fsc_X2./detQ0./Ayy_l))) ...
      - ...
      cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), h_s))) .* ...
      (ones(3,1)*(Cte.*sqrt(pi)/k^(3/2) ...
      ./ sum_r_s  ./s_X2.* exp(-k*g_X2) ...
      .* sqrt(det_Qz_fsc_s./detQ0./detA_l)));
    
E_X2 = +E_X2;

E_Y1 = ...
      Epc1.*sqrt(1/pi).*(ones(3,1)*subfun_Transfun(sqrt(k).*s_Y1)) ...
      + ...
     cross(rr_Y1, cross(rr_Y1, cross(Plaque.N*ones(1,NbP), h_Y1))) .* ...
      (ones(3,1)*(Cte*sqrt(2*pi)/k^(3/2)  ...
      ./ sum_r_Y1 ./dgdy_Y1_l .* exp(-k*g_Y1)...
      .* sqrt(det_Qz_fsc_Y1./detQ0./Axx_l))) ...
      - ...
      cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), h_s))) .* ...
      (ones(3,1)*(Cte.*sqrt(pi)/k^(3/2) ...
      ./ sum_r_s ./s_Y1 .* exp(-k*g_Y1) ...
      .* sqrt(det_Qz_fsc_s./detQ0./detA_l)));
    
E_Y1 = -E_Y1;

E_Y2 = ...
      Epc1.*sqrt(1/pi).*(ones(3,1)*(subfun_Transfun(sqrt(k).*s_Y2))) ...
      + ...
     cross(rr_Y2, cross(rr_Y2, cross(Plaque.N*ones(1,NbP), h_Y2))) .* ...
      (ones(3,1)*(Cte*sqrt(2*pi)/k^(3/2)  ...
      ./ sum_r_Y2 ./dgdy_Y2_l .* exp(-k*g_Y2)...
      .* sqrt(det_Qz_fsc_Y2./detQ0./Axx_l))) ...
      - ...
      cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), h_s))) .* ...
      (ones(3,1)*(Cte.*sqrt(pi)/k^(3/2) ...
      ./ sum_r_s ./s_Y2 .* exp(-k*g_Y2) ...
      .* sqrt(det_Qz_fsc_s./detQ0./detA_l)));
    
E_Y2 = +E_Y2;


%  error('test')

%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

%  %  %  ----------  Formulation uniforme 
%  nu_X1 = -dgdx_X1_l.*sqrt(k/2./Axx_l)
%  %  nu_X2 = dgdx_X2_l.*sqrt(k/2./Axx_l);
%  %  nu_Y1 = dgdy_Y1_l.*sqrt(k/2./Ayy_l);
%  %  nu_Y2 = dgdy_Y2_l.*sqrt(k/2./Ayy_l);
%  %  nu_X1 = abs(dgdx_X1).*nsqrt(k./2./abs(-j*Axx));
%  %  nu_X2 = abs(dgdx_X2).*nsqrt(k./2./abs(-j*Axx));
%  %  nu_Y1 = abs(dgdy_Y1).*nsqrt(k./2./abs(-j*Ayy));
%  %  nu_Y2 = abs(dgdy_Y2).*nsqrt(k./2./abs(-j*Ayy));
%  %  
%  %  %  % Lorsque le point critique du premier ordre n'appartient plus 
%  %  %  % au domaine d'integration, on supprime le coefficient 1/2
%  %  %  % car le point critique du premier ordre devient alors un 
%  %  %  % point du second ordre
%  %  Cte = j*k/2/pi ;%.* (1 + U1);
%  %  
%  E_X1 = cross(rr_X1, cross(rr_X1, cross(Plaque.N*ones(1,NbP), h_X1))) .* ...
%        (ones(3,1)*(U_X1 .* Cte .* ...
%        (+1)*2/k ...
%        .* ampl_X1 ...
%        .* exp(-k*g_X1) ...
%        .* sqrt(pi.*det_Qz_fsc_X1./detQ0./(Axx_l)./(Ayy_l)) ...
%        .* Fmin(nu_X1) ...
%        .* exp(sign(-real(Axx_l)).*nu_X1.^2) ...
%        .*sign(-real(dgdx_X1_l))));


%  
%  E_X2 = cross(rr_X2, cross(rr_X2, cross(Plaque.N*ones(1,NbP), h_X2))) .* ...
%        (ones(3,1)*(U_X2 .* Cte .* ...
%        (-1)*2/k ...
%        .* ampl_X2 ...
%        .* exp(-k*g_X2) ...
%        .* sqrt(pi.*det_Qz_fsc_X1./detQ0./(Axx_l)./(Ayy_l)) ...
%        .* Fmin(nu_X2) ...
%        .* exp(sign(real(Axx_l)).*nu_X2.^2).*sign(-imag(dgdx_X2_l))));
%  
%  E_Y1 = cross(rr_Y1, cross(rr_Y1, cross(Plaque.N*ones(1,NbP), h_Y1))) .* ...
%        (ones(3,1)*(U_Y1 .* Cte .* ...
%        (+1)*2/k ...
%        .* ampl_Y1 ...
%        .* exp(-k*g_Y1) ...
%        .* sqrt(pi.*det_Qz_fsc_X1./detQ0./(Axx_l)./(Ayy_l)) ...
%        .* Fmin(nu_Y1) ...
%        .* exp(sign(real(Ayy_l)).*nu_Y1.^2).*sign(real(dgdy_Y1_l))    ));
%  
%  E_Y2 = cross(rr_Y2, cross(rr_Y2, cross(Plaque.N*ones(1,NbP), h_Y2))) .* ...
%        (ones(3,1)*(U_Y2 .* Cte .* ...
%        (-1)*2/k ...
%        .* ampl_Y2 ...
%        .* exp(-k*g_Y2) ...
%        .* sqrt(pi.*det_Qz_fsc_X1./detQ0./(Axx_l)./(Ayy_l)) ...
%        .* Fmin(nu_Y2) ...
%        .* exp(sign(real(Ayy_l)).*nu_Y2.^2).*sign(real(dgdy_Y2_l))  ));


%  
%  % =======================================================
%  % Points critiques du troisieme ordre (coins)
%  % =======================================================
%  
%  % vecteur distance aux points critiques du troisieme ordre
%     rQ_X1Y1 = [X1*ones(1,NbP); Y1*ones(1,NbP);zeros(1,NbP)]; 
%      r_X1Y1 = P - rQ_X1Y1;
%  sum_r_X1Y1 = nsqrt(sum(r_X1Y1.^2));
%     rr_X1Y1 = r_X1Y1./(ones(3,1)*sum_r_X1Y1);
%  
%     rQ_X1Y2 = [X1*ones(1,NbP); Y2*ones(1,NbP);zeros(1,NbP)]; 
%      r_X1Y2 = P - rQ_X1Y2;
%  sum_r_X1Y2 = nsqrt(sum(r_X1Y2.^2));
%     rr_X1Y2 = r_X1Y2./(ones(3,1)*sum_r_X1Y2);
%  
%     rQ_X2Y1 = [X2*ones(1,NbP); Y1*ones(1,NbP);zeros(1,NbP)]; 
%      r_X2Y1 = P - rQ_X2Y1;
%  sum_r_X2Y1 = nsqrt(sum(r_X2Y1.^2));
%     rr_X2Y1 = r_X2Y1./(ones(3,1)*sum_r_X2Y1);
%  
%     rQ_X2Y2 = [X2*ones(1,NbP); Y2*ones(1,NbP);zeros(1,NbP)]; 
%      r_X2Y2 = P - rQ_X2Y2;
%  sum_r_X2Y2 = nsqrt(sum(r_X2Y2.^2));
%     rr_X2Y2 = r_X2Y2./(ones(3,1)*sum_r_X2Y2);
%  
%  % point col dans le repere incident
%  xQ_fsc_X1Y1 = r11.*rQ_X1Y1(1,:) + r12.*rQ_X1Y1(2,:) - xi_Oi;
%  yQ_fsc_X1Y1 = r21.*rQ_X1Y1(1,:) + r22.*rQ_X1Y1(2,:) - yi_Oi;
%  zQ_fsc_X1Y1 = r31.*rQ_X1Y1(1,:) + r32.*rQ_X1Y1(2,:) - zi_Oi;
%  
%  xQ_fsc_X1Y2 = r11.*rQ_X1Y2(1,:) + r12.*rQ_X1Y2(2,:) - xi_Oi;
%  yQ_fsc_X1Y2 = r21.*rQ_X1Y2(1,:) + r22.*rQ_X1Y2(2,:) - yi_Oi;
%  zQ_fsc_X1Y2 = r31.*rQ_X1Y2(1,:) + r32.*rQ_X1Y2(2,:) - zi_Oi;
%  
%  xQ_fsc_X2Y1 = r11.*rQ_X2Y1(1,:) + r12.*rQ_X2Y1(2,:) - xi_Oi;
%  yQ_fsc_X2Y1 = r21.*rQ_X2Y1(1,:) + r22.*rQ_X2Y1(2,:) - yi_Oi;
%  zQ_fsc_X2Y1 = r31.*rQ_X2Y1(1,:) + r32.*rQ_X2Y1(2,:) - zi_Oi;
%  
%  xQ_fsc_X2Y2 = r11.*rQ_X2Y2(1,:) + r12.*rQ_X2Y2(2,:) - xi_Oi;
%  yQ_fsc_X2Y2 = r21.*rQ_X2Y2(1,:) + r22.*rQ_X2Y2(2,:) - yi_Oi;
%  zQ_fsc_X2Y2 = r31.*rQ_X2Y2(1,:) + r32.*rQ_X2Y2(2,:) - zi_Oi;
%  
%  
%  
%  % matrice de courbure dans le repere incident
%  iQz_fsc_X1Y1 = [invQ0(1)+zQ_fsc_X1Y1; ...
%                  invQ0(2)+zQ_fsc_X1Y1; ...
%                  invQ0(3)*ones(1,NbP)];
%  det_iQz_fsc_X1Y1 = iQz_fsc_X1Y1(1,:).*iQz_fsc_X1Y1(2,:) - iQz_fsc_X1Y1(3,:).^2;
%  Qz_fsc_X1Y1 = [iQz_fsc_X1Y1(2,:); ...
%                 iQz_fsc_X1Y1(1,:); ...
%                -iQz_fsc_X1Y1(3,:)]./(ones(3,1)*det_iQz_fsc_X1Y1);
%  det_Qz_fsc_X1Y1 = Qz_fsc_X1Y1(1,:).*Qz_fsc_X1Y1(2,:) - Qz_fsc_X1Y1(3,:).^2;
%  
%  iQz_fsc_X1Y2 = [invQ0(1)+zQ_fsc_X1Y2; ...
%                  invQ0(2)+zQ_fsc_X1Y2; ...
%                  invQ0(3)*ones(1,NbP)];
%  det_iQz_fsc_X1Y2 = iQz_fsc_X1Y2(1,:).*iQz_fsc_X1Y2(2,:) - iQz_fsc_X1Y2(3,:).^2;
%  Qz_fsc_X1Y2 = [iQz_fsc_X1Y2(2,:); ...
%                 iQz_fsc_X1Y2(1,:); ...
%                -iQz_fsc_X1Y2(3,:)]./(ones(3,1)*det_iQz_fsc_X1Y2);
%  det_Qz_fsc_X1Y2 = Qz_fsc_X1Y2(1,:).*Qz_fsc_X1Y2(2,:) - Qz_fsc_X1Y2(3,:).^2;
%  
%  iQz_fsc_X2Y1 = [invQ0(1)+zQ_fsc_X2Y1; ...
%                  invQ0(2)+zQ_fsc_X2Y1; ...
%                  invQ0(3)*ones(1,NbP)];
%  det_iQz_fsc_X2Y1 = iQz_fsc_X2Y1(1,:).*iQz_fsc_X2Y1(2,:) - iQz_fsc_X2Y1(3,:).^2;
%  Qz_fsc_X2Y1 = [iQz_fsc_X2Y1(2,:); ...
%                 iQz_fsc_X2Y1(1,:); ...
%                -iQz_fsc_X2Y1(3,:)]./(ones(3,1)*det_iQz_fsc_X2Y1);
%  det_Qz_fsc_X2Y1 = Qz_fsc_X2Y1(1,:).*Qz_fsc_X2Y1(2,:) - Qz_fsc_X2Y1(3,:).^2;
%  
%  iQz_fsc_X2Y2 = [invQ0(1)+zQ_fsc_X2Y2; ...
%                  invQ0(2)+zQ_fsc_X2Y2; ...
%                  invQ0(3)*ones(1,NbP)];
%  det_iQz_fsc_X2Y2 = iQz_fsc_X2Y2(1,:).*iQz_fsc_X2Y2(2,:) - iQz_fsc_X2Y2(3,:).^2;
%  Qz_fsc_X2Y2 = [iQz_fsc_X2Y2(2,:); ...
%                 iQz_fsc_X2Y2(1,:); ...
%                -iQz_fsc_X2Y2(3,:)]./(ones(3,1)*det_iQz_fsc_X2Y2);
%  det_Qz_fsc_X2Y2 = Qz_fsc_X2Y2(1,:).*Qz_fsc_X2Y2(2,:) - Qz_fsc_X2Y2(3,:).^2;
%  
%  % vecteur polarisation aux points critiques du second ordre
%  h_l_x = -Champ.Coefficients(2).*ones(1,NbP);
%  h_l_y =  Champ.Coefficients(1).*ones(1,NbP);
%  
%  h_l_X1Y1_z = -Champ.Coefficients(1).*(xQ_fsc_X1Y1.*Qz_fsc_X1Y1(3,:) + yQ_fsc_X1Y1.*Qz_fsc_X1Y1(2,:)) ...
%               +Champ.Coefficients(2).*(xQ_fsc_X1Y1.*Qz_fsc_X1Y1(1,:) + yQ_fsc_X1Y1.*Qz_fsc_X1Y1(3,:));
%  h_l_X1Y1 = [h_l_x;h_l_y;h_l_X1Y1_z];
%  h_X1Y1   = M_fsc*h_l_X1Y1;
%  
%  h_l_X1Y2_z = -Champ.Coefficients(1).*(xQ_fsc_X1Y2.*Qz_fsc_X1Y2(3,:) + yQ_fsc_X1Y2.*Qz_fsc_X1Y2(2,:)) ...
%               +Champ.Coefficients(2).*(xQ_fsc_X1Y2.*Qz_fsc_X1Y2(1,:) + yQ_fsc_X1Y2.*Qz_fsc_X1Y2(3,:));
%  h_l_X1Y2 = [h_l_x;h_l_y;h_l_X1Y2_z];
%  h_X1Y2   = M_fsc*h_l_X1Y2;
%  
%  h_l_X2Y1_z = -Champ.Coefficients(1).*(xQ_fsc_X2Y1.*Qz_fsc_X2Y1(3,:) + yQ_fsc_X2Y1.*Qz_fsc_X2Y1(2,:)) ...
%               +Champ.Coefficients(2).*(xQ_fsc_X2Y1.*Qz_fsc_X2Y1(1,:) + yQ_fsc_X2Y1.*Qz_fsc_X2Y1(3,:));
%  h_l_X2Y1 = [h_l_x;h_l_y;h_l_X2Y1_z];
%  h_X2Y1   = M_fsc*h_l_X2Y1;
%  
%  h_l_X2Y2_z = -Champ.Coefficients(1).*(xQ_fsc_X2Y2.*Qz_fsc_X2Y2(3,:) + yQ_fsc_X2Y2.*Qz_fsc_X2Y2(2,:)) ...
%               +Champ.Coefficients(2).*(xQ_fsc_X2Y2.*Qz_fsc_X2Y2(1,:) + yQ_fsc_X2Y2.*Qz_fsc_X2Y2(3,:));
%  h_l_X2Y2 = [h_l_x;h_l_y;h_l_X2Y2_z];
%  h_X2Y2   = M_fsc*h_l_X2Y2;
%  
%  % fonction de phase aux points critiques
%  g_X1Y1 = -j/2*(rQ_X1Y1(1,:).^2.*Axx ...
%              + 2*rQ_X1Y1(1,:).*rQ_X1Y1(2,:).*Axy ...
%              + rQ_X1Y1(2,:).^2.*Ayy) ...
%      + j*Bx.*rQ_X1Y1(1,:) + j*By.*rQ_X1Y1(2,:) ...
%      + j*C;
%  
%  g_X1Y2 = -j/2*(rQ_X1Y2(1,:).^2.*Axx ...
%              + 2*rQ_X1Y2(1,:).*rQ_X1Y2(2,:).*Axy ...
%              + rQ_X1Y2(2,:).^2.*Ayy) ...
%      + j*Bx.*rQ_X1Y2(1,:) + j*By.*rQ_X1Y2(2,:) ...
%      + j*C;
%  
%  g_X2Y1 = -j/2*(rQ_X2Y1(1,:).^2.*Axx ...
%              + 2*rQ_X2Y1(1,:).*rQ_X2Y1(2,:).*Axy ...
%              + rQ_X2Y1(2,:).^2.*Ayy) ...
%      + j*Bx.*rQ_X2Y1(1,:) + j*By.*rQ_X2Y1(2,:) ...
%      + j*C;
%  
%  g_X2Y2 = -j/2*(rQ_X2Y2(1,:).^2.*Axx ...
%              + 2*rQ_X2Y2(1,:).*rQ_X2Y2(2,:).*Axy ...
%              + rQ_X2Y2(2,:).^2.*Ayy) ...
%      + j*Bx.*rQ_X2Y2(1,:) + j*By.*rQ_X2Y2(2,:) ...
%      + j*C;
%  
%  % derivees de la fonction de phase aux points critiques
%  dgdx_X1Y1 = -j*rQ_X1Y1(1,:).*Axx - j*rQ_X1Y1(2,:).*Axy + j*Bx;
%  dgdy_X1Y1 = -j*rQ_X1Y1(1,:).*Axy - j*rQ_X1Y1(2,:).*Ayy + j*By;
%  
%  dgdx_X1Y2 = -j*rQ_X1Y2(1,:).*Axx - j*rQ_X1Y2(2,:).*Axy + j*Bx;
%  dgdy_X1Y2 = -j*rQ_X1Y2(1,:).*Axy - j*rQ_X1Y2(2,:).*Ayy + j*By;
%  
%  dgdx_X2Y1 = -j*rQ_X2Y1(1,:).*Axx - j*rQ_X2Y1(2,:).*Axy + j*Bx;
%  dgdy_X2Y1 = -j*rQ_X2Y1(1,:).*Axy - j*rQ_X2Y1(2,:).*Ayy + j*By;
%  
%  dgdx_X2Y2 = -j*rQ_X2Y2(1,:).*Axx - j*rQ_X2Y2(2,:).*Axy + j*Bx;
%  dgdy_X2Y2 = -j*rQ_X2Y2(1,:).*Axy - j*rQ_X2Y2(2,:).*Ayy + j*By;
%  
%  % fonction d'amplitude au point col
%  ampl_X1Y1 = sqrt(det_Qz_fsc_X1Y1./detQ0) ...
%              ./sum_r_X1Y1;
%  
%  ampl_X1Y2 = sqrt(det_Qz_fsc_X1Y2./detQ0) ...
%              ./sum_r_X1Y2;
%  
%  ampl_X2Y1 = sqrt(det_Qz_fsc_X2Y1./detQ0) ...
%              ./sum_r_X2Y1;
%  
%  ampl_X2Y2 = sqrt(det_Qz_fsc_X2Y2./detQ0) ...
%              ./sum_r_X2Y2;
%  
%  %  % application de la formule non uniforme
%  Cte = j*k/2/pi;
%  E_X1Y1 = - cross(rr_X1Y1, cross(rr_X1Y1, cross(Plaque.N*ones(1,NbP), h_X1Y1))) .* ...
%        (ones(3,1)*(Cte ...
%          .* ampl_X1Y1 ...
%          .* exp(-k*g_X1Y1) ...
%          ./ k^2 ./ dgdx_X1Y1 ./ dgdy_X1Y1));
%  E_X1Y2 = + cross(rr_X1Y2, cross(rr_X1Y2, cross(Plaque.N*ones(1,NbP), h_X1Y2))) .* ...
%        (ones(3,1)*(Cte ...
%          .* ampl_X1Y2 ...
%          .* exp(-k*g_X1Y2) ...
%          ./ k^2 ./ dgdx_X1Y2 ./ dgdy_X1Y2));
%  E_X2Y1 = + cross(rr_X2Y1, cross(rr_X2Y1, cross(Plaque.N*ones(1,NbP), h_X2Y1))) .* ...
%        (ones(3,1)*(Cte ...
%          .* ampl_X2Y1 ...
%          .* exp(-k*g_X2Y1) ...
%          ./ k^2 ./ dgdx_X2Y1 ./ dgdy_X2Y1));
%  E_X2Y2 = - cross(rr_X2Y2, cross(rr_X2Y2, cross(Plaque.N*ones(1,NbP), h_X2Y2))) .* ...
%        (ones(3,1)*(Cte ...
%          .* ampl_X2Y2 ...
%          .* exp(-k*g_X2Y2) ...
%          ./ k^2 ./ dgdx_X2Y2 ./ dgdy_X2Y2));
%  
%  %  % ===================================
%  %  %  version uniforme
%  %  % ===================================
%  %  nu_X1 = abs(dgdx_X1).*nsqrt(k/2./abs(-Axx));
%  %  nu_X2 = abs(dgdx_X2).*nsqrt(k/2./abs(-Axx));
%  %  nu_Y1 = abs(dgdy_Y1).*nsqrt(k/2./abs(-Ayy));
%  %  nu_Y2 = abs(dgdy_Y2).*nsqrt(k/2./abs(-Ayy));
%  %  T_X1 = sign(real(-Axx)) .* Fmin(nu_X1,7) .* exp(-j*sign(real(-Axx)).*nu_X1.^2);
%  %  T_X2 = sign(real(-Axx)) .* Fmin(nu_X2,7) .* exp(-j*sign(real(-Axx)).*nu_X2.^2);
%  %  T_Y1 = sign(real(-Ayy)) .* Fmin(nu_Y1,7) .* exp(-j*sign(real(-Ayy)).*nu_Y1.^2);
%  %  T_Y2 = sign(real(-Ayy)) .* Fmin(nu_Y2,7) .* exp(-j*sign(real(-Ayy)).*nu_Y2.^2);
%  %  
%  %  E_X1Y1 = cross(rr_X1Y1, cross(rr_X1Y1, cross(Plaque.N*ones(1,NbP), h_X1Y1))) .* ...
%  %        (ones(3,1)*(Cte ...
%  %          .* (-1) .* 2/k ...
%  %          .* ampl_X1Y1 ...
%  %          .* exp(j*k*g_X1Y1) ...
%  %          .* unSurSqrtAxx .* unSurSqrtAyy ...
%  %          ...%.* sign(real(dgdx_X1)).* sign(real(dgdy_Y1)) ...
%  %          .* T_X1 .* T_Y1        ));
%  %  
%  %  E_X1Y2 = cross(rr_X1Y2, cross(rr_X1Y2, cross(Plaque.N*ones(1,NbP), h_X1Y2))) .* ...
%  %        (ones(3,1)*(Cte ...
%  %          .* (+1) .* 2/k ...
%  %          .* ampl_X1Y1 ...
%  %          .* exp(j*k*g_X1Y2) ...
%  %          .* unSurSqrtAxx .* unSurSqrtAyy ...
%  %          .* sign(real(dgdx_X1)).* sign(real(dgdy_Y2)) ...
%  %          .* T_X1 .* T_Y2        ));
%  %  
%  %  E_X2Y1 = cross(rr_X2Y1, cross(rr_X2Y1, cross(Plaque.N*ones(1,NbP), h_X2Y1))) .* ...
%  %        (ones(3,1)*(Cte ...
%  %          .* (+1) .* 2/k ...
%  %          .* ampl_X2Y1 ...
%  %          .* exp(j*k*g_X2Y1) ...
%  %          .* unSurSqrtAxx .* unSurSqrtAyy ...
%  %          .* sign(real(dgdx_X2)).* sign(real(dgdy_Y1)) ...
%  %          .* T_X2 .* T_Y1        ));
%  %  
%  %  E_X2Y2 = cross(rr_X2Y2, cross(rr_X2Y2, cross(Plaque.N*ones(1,NbP), h_X2Y2))) .* ...
%  %        (ones(3,1)*(Cte ...
%  %          .* (-1) .* 2/k ...
%  %          .* ampl_X2Y1 ...
%  %          .* exp(j*k*g_X2Y2) ...
%  %          .* unSurSqrtAxx .* unSurSqrtAyy ...
%  %          .* sign(real(dgdx_X2)).* sign(real(dgdy_Y2)) ...
%  %          .* T_X2 .* T_Y2        ));

% ==========================================================
% champ total
% ==========================================================
%  %  filtrage des NaN
%  E_X1(isnan(E_X1)) = eps;
%  E_X2(isnan(E_X2)) = eps;
%  E_Y1(isnan(E_Y1)) = eps;
%  E_Y2(isnan(E_Y2)) = eps;
%  Epc1(isnan(Epc1)) = eps;

Epc2 =  E_X1 + E_X2  + E_Y1 + E_Y2;
Epc3 = [];%E_X1Y1 +  E_X1Y2 + E_X2Y1 + E_X2Y2;

  E = Epc1+Epc2 ;%

H = zeros(size(E));








%  Fonction de phase
function g = subfun_g_l(xQ_l,yQ_l)
global Axx_l Axy_l Ayy_l Bx_l By_l C_l;  
  g = 1/2*(xQ_l.^2.*Axx_l + 2*xQ_l.*yQ_l.*Axy_l + yQ_l.^2.*Ayy_l) ...
      - Bx_l.*xQ_l - By_l.*yQ_l ...
      + C_l;

function [Qz, det_Qz] = subfun_Qz(z_fsc, Champ)
%  Determine les valeurs de la matrice de courbure d'un faisceau gaussien
%  propage au point P
%  
%  ENTREE
%  z_fsc : distance d'evaluation (exprimes dans le repere faisceau)
%  Champ : structure Champ
%  
%  RETOURNE :
%  Qz : matrice de courbure propagee sur la distance z_fsc
%  det_Qz : determinant de cette matrice de courbure a la distance z_fsc

  % matrice de courbure en z_fc = 0
  Q0    = Champ.Courbure;
  detQ0 = Q0(1)*Q0(2)-Q0(3)^2;  
  invQ0 = [Q0(2); Q0(1); -Q0(3)] ./ (ones(3,1)*detQ0);
  % matrice de courbure en z_fc
  iQz = [invQ0(1)+z_fsc; ...
         invQ0(2)+z_fsc; ...
         invQ0(3)*ones(1,size(z_fsc,2))];
  det_iQz = iQz(1,:).*iQz(2,:) - iQz(3,:).^2;
  Qz = [iQz(2,:); iQz(1,:); -iQz(3,:)]./(ones(3,1)*det_iQz);
  % determinant de la matrice de courbure en z_fc
  det_Qz = Qz(1,:).*Qz(2,:) - Qz(3,:).^2;

function h = subfun_h(P, Champ)
%  fonction de polarisation dans le repere faisceau
%  
%  ENTREES :
%  P : points d'evaluation (exprimes dans le repere faisceau)
%  Champ : structure Champ
%  
%  RETOURNE :
%  h : vecteur polarisation (3xNbP)
%  det_Qz : determinant de la matrice de courbure evaluee au point P
%  
%  ATTENTION : ne plas oulier de changer le repere du vecteur h
%  dans le repere ansolu par la suite.
  NbP = size(P,2);
  
  [Qz, det_Qz] = subfun_Qz(P(3,:), Champ);  
  
  h_x = -Champ.Coefficients(2).*ones(1,NbP);
  h_y =  Champ.Coefficients(1).*ones(1,NbP);
      Qx = Qz(1,:).*P(1,:) + Qz(3,:).*P(2,:);
      Qy = Qz(3,:).*P(1,:) + Qz(2,:).*P(2,:);
  h_z = -Champ.Coefficients(1).*Qy + Champ.Coefficients(2).*Qx;
    h = [h_x;h_y;h_z];


function F = subfun_Transfun(z)
%  Fonction de transition pour dev. asyptotique uniforme
F = sqrt(pi)/2.*(1-erfz(z));

%  fonction signe sans le 0 si x = 0
function s = sign2(x)
s = +1*(x>=0) -1*(x<0);  

  