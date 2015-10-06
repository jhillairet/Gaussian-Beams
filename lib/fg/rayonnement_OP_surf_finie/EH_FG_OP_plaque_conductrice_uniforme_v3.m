function [E,H, Epc1, Epc2, Epc3] = EH_FG_OP_Plaque_conductrice_uniforme(Champ, Plaque, P, f, epsr)
% Calcul du champ EM rayonne par une Plaque conductrice 
% dans l'hypothese de l'optique physique
% par un faisceau gaussien
% 
% Les effets de la finitude sont pris en compte
% Ainsi que les effets "de coin" de maniere cette fois uniforme
% 
% La normale de la plaque est consideree comme etant l'axe z.
% 
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
global R_Rfsc2Rabs epsr detQ0 detA_l;


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
%  avec le plan. Ce point correspond � l'origine du repere local
%  dans lequel on va faire les calculs.
%  
%  Distance � parcourir entre le centre du repere faisceau
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

%  matrice de courbure
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

% centre du faisceau exprim� dans le repere faisceau
r_l_i_Oi = R_Rabs2Rfsc*(Champ.Centre - O_l);
x_l_i_Oi = r_l_i_Oi(1);
y_l_i_Oi = r_l_i_Oi(2);
z_l_i_Oi = r_l_i_Oi(3);

ri_Oi = R_Rabs2Rfsc*Champ.Centre;
xi_Oi = ri_Oi(1);
yi_Oi = ri_Oi(2);
zi_Oi = ri_Oi(3);

M_fsc = R_Rfsc2Rabs;

% elements de la forme quadratique (QF)
% dans sa notation matricielle
QF.Axx_l = j*(1./R_l - P_l(1,:).^2./R_l.^3 ...
  + Qp(1).*r11.^2 ...
  + Qp(2).*r21.^2 ...
  + 2*Qp(3).*r11.*r21);

QF.Ayy_l =  j*(1./R_l - P_l(2,:).^2./R_l.^3 ...
  + Qp(2).*r22.^2 ...
  + Qp(1).*r12.^2 ...
  + 2*Qp(3).*r12.*r22);

QF.Axy_l = j*( - P_l(1,:).*P_l(2,:)./R_l.^3 ...
  + Qp(3).*(r11.*r22 + r12.*r21) ...
  + Qp(2).*r21.*r22 ...
  + Qp(1).*r11.*r12);

QF.Bx_l = j*(- r31 ...
  + P_l(1,:)./R_l ...
  + Qp(3).*(x_l_i_Oi.*r21 + r11.*y_l_i_Oi) ...
  + Qp(2).*r21.*y_l_i_Oi ...
  + Qp(1).*r11.*x_l_i_Oi);

QF.By_l = j*(- r32 ...
  + P_l(2,:)./R_l ...
  + Qp(3).*(x_l_i_Oi.*r22 + r12.*y_l_i_Oi) ...
  + Qp(2).*r22.*y_l_i_Oi ...
  + Qp(1).*r12.*x_l_i_Oi);

QF.C_l = j*(R_l ...
  - z_l_i_Oi ...
  +1/2*Qp(1).*x_l_i_Oi.^2 ...
  +1/2*Qp(2).*y_l_i_Oi.^2 ...
  +Qp(3).*x_l_i_Oi.*y_l_i_Oi);

detA_l = QF.Axx_l.*QF.Ayy_l - QF.Axy_l.^2;
iA = [QF.Ayy_l;QF.Axx_l;-QF.Axy_l]./(ones(3,1)*detA_l);



% ======================================================
% Calcul du point col (point critique du premier ordre)
% ======================================================
%  point col dans le repere absolu modifie
xQ_l_s = iA(1,:).*QF.Bx_l + iA(3,:).*QF.By_l;
yQ_l_s = iA(3,:).*QF.Bx_l + iA(2,:).*QF.By_l;
rQ_l_s = [xQ_l_s; yQ_l_s; zeros(1,NbP)];
%  point col dans le repere absolu 
xQ_s = xQ_l_s + O_l(1);
yQ_s = yQ_l_s + O_l(2);
rQ_s = [xQ_s; yQ_s; Plaque.Centre(3)*ones(1,NbP)];

% point col dans le repere fsc
rQ_fsc_s = R_Rabs2Rfsc*(rQ_s  - Champ.Centre*ones(1,NbP) );

% vecteur distance au point col
    r_s = P_l - rQ_l_s;
sum_r_s = sqrt(sum((r_s).^2));
   rr_s = r_s./(ones(3,1)*sum_r_s);

%  matrice de courbure evaluee au point critique du premier ordre
[Qz_s, det_Qz_fsc_s] = subfun_Qz(rQ_fsc_s(3,:), Champ);  
%  fonction de phase evaluee au point critique du premier ordre
    g_s_l = subfun_g_l(xQ_l_s, yQ_l_s, QF);
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

% =====================================================
% point critique de premier ordre (point col)
% =====================================================
Epc1 = cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), h_s))) .* ...
      (ones(3,1)*(ampl_s .* exp(-k*g_s_l)));


  
% =====================================================
% Contribution des points critiques du second ordre
% =====================================================
  % points critiques du second ordre
  % dans le repere modifie correspondant a l intersection
  % du FG avec la plaque
  x_Y1_l = (QF.Bx_l-QF.Axy_l.*Y1_l)./QF.Axx_l;
  x_Y2_l = (QF.Bx_l-QF.Axy_l.*Y2_l)./QF.Axx_l;
  y_X1_l = (QF.By_l-QF.Axy_l.*X1_l)./QF.Ayy_l;
  y_X2_l = (QF.By_l-QF.Axy_l.*X2_l)./QF.Ayy_l;

  % points critiques du second ordre
  % dans le repere absolu lie a la plaque
  x_Y1 = x_Y1_l + O_l(1);
  x_Y2 = x_Y2_l + O_l(1);
  y_X1 = y_X1_l + O_l(2);
  y_X2 = y_X2_l + O_l(2);

  % vecteur distance aux points critiques du second ordre
    rQ_X1 = [X1*ones(1,NbP); y_X1;Plaque.Centre(3)*ones(1,NbP)]; 
  rQ_X1_l = [X1_l*ones(1,NbP); y_X1_l;Plaque.Centre(3)*ones(1,NbP)]; 
      r_X1 = P_l - rQ_X1_l;
  sum_r_X1 = sqrt(sum((r_X1.^2)));
    rr_X1 = r_X1./(ones(3,1)*sum_r_X1);
  
    rQ_X2 = [X2*ones(1,NbP); y_X2;Plaque.Centre(3)*ones(1,NbP)]; 
  rQ_X2_l = [X2_l*ones(1,NbP); y_X2_l;Plaque.Centre(3)*ones(1,NbP)]; 
      r_X2 = P_l - rQ_X2_l;
  sum_r_X2 = sqrt(sum((r_X2.^2)));
    rr_X2 = r_X2./(ones(3,1)*sum_r_X2);
  
    rQ_Y1 = [x_Y1; Y1*ones(1,NbP); Plaque.Centre(3)*ones(1,NbP)]; 
  rQ_Y1_l = [x_Y1_l; Y1_l*ones(1,NbP); Plaque.Centre(3)*ones(1,NbP)]; 
      r_Y1 = P_l - rQ_Y1_l;
  sum_r_Y1 = sqrt(sum((r_Y1.^2)));
    rr_Y1 = r_Y1./(ones(3,1)*sum_r_Y1);
  
    rQ_Y2 = [x_Y2; Y2*ones(1,NbP); Plaque.Centre(3)*ones(1,NbP)]; 
  rQ_Y2_l = [x_Y2_l; Y2_l*ones(1,NbP); Plaque.Centre(3)*ones(1,NbP)]; 
      r_Y2 = P_l - rQ_Y2_l;
  sum_r_Y2 = sqrt(sum((r_Y2).^2));
    rr_Y2 = r_Y2./(ones(3,1)*sum_r_Y2);


  % point critiques du second ordre dans le repere du faisceau incident
  rQ_fsc_X1 = R_Rabs2Rfsc*(rQ_X1  - Champ.Centre*ones(1,NbP));
  rQ_fsc_X2 = R_Rabs2Rfsc*(rQ_X2  - Champ.Centre*ones(1,NbP));
  rQ_fsc_Y1 = R_Rabs2Rfsc*(rQ_Y1  - Champ.Centre*ones(1,NbP));
  rQ_fsc_Y2 = R_Rabs2Rfsc*(rQ_Y2  - Champ.Centre*ones(1,NbP));

  %  Matrices de courbures aux points critiques du second ordre
  [Qz_fsc_X1, det_Qz_fsc_X1] = subfun_Qz(rQ_fsc_X1(3,:), Champ);
  [Qz_fsc_X2, det_Qz_fsc_X2] = subfun_Qz(rQ_fsc_X2(3,:), Champ);
  [Qz_fsc_Y1, det_Qz_fsc_Y1] = subfun_Qz(rQ_fsc_Y1(3,:), Champ);
  [Qz_fsc_Y2, det_Qz_fsc_Y2] = subfun_Qz(rQ_fsc_Y2(3,:), Champ);

  %  Fonctions de polarisations aux points critiques du second ordre
  h_X1 = M_fsc*subfun_h(rQ_fsc_X1, Champ);
  h_X2 = M_fsc*subfun_h(rQ_fsc_X2, Champ);
  h_Y1 = M_fsc*subfun_h(rQ_fsc_Y1, Champ);
  h_Y2 = M_fsc*subfun_h(rQ_fsc_Y2, Champ);

  %  Fonctions de phase aux points critiques du second ordre
  g_X1 = subfun_g_l(rQ_X1_l(1,:), rQ_X1_l(2,:), QF);
  g_X2 = subfun_g_l(rQ_X2_l(1,:), rQ_X2_l(2,:), QF);
  g_Y1 = subfun_g_l(rQ_Y1_l(1,:), rQ_Y1_l(2,:), QF);
  g_Y2 = subfun_g_l(rQ_Y2_l(1,:), rQ_Y2_l(2,:), QF);

  % derivees de la fonction de phase aux points critiques
  dgdx_X1_l = (QF.Axx_l-QF.Axy_l.^2./QF.Ayy_l).*X1_l-QF.Bx_l+QF.By_l.*QF.Axy_l./QF.Ayy_l;
  dgdx_X2_l = (QF.Axx_l-QF.Axy_l.^2./QF.Ayy_l).*X2_l-QF.Bx_l+QF.By_l.*QF.Axy_l./QF.Ayy_l;
  
  dgdy_Y1_l = (-QF.Axy_l.^2./QF.Axx_l+QF.Ayy_l).*Y1_l+QF.Bx_l.*QF.Axy_l./QF.Axx_l-QF.By_l;
  dgdy_Y2_l = (-QF.Axy_l.^2./QF.Axx_l+QF.Ayy_l).*Y2_l+QF.Bx_l.*QF.Axy_l./QF.Axx_l-QF.By_l;


  %  Constante multiplicative
  Cte = j*k/2/pi*sqrt(epsr) ;

  %  terme de changement de variable dans le passage 
  %  aux coordonnees "col" : dx/ds
  hh_s_X  = sqrt(2./(QF.Axx_l-QF.Axy_l.^2./QF.Ayy_l));
  hh_s_Y  = sqrt(2./(-QF.Axy_l.^2./QF.Axx_l+QF.Ayy_l));


  %  Le signe de la racine carr�e est d�fini tel que l'on doit avoir
  %  h_x -> h_s lorsque s -> 0, cad lq 
  s_X1 = sign(real(hh_s_X./(X1_l-xQ_l_s))).*sqrt(g_X1 - g_s_l);
  s_X2 = sign(real(hh_s_X./(X2_l-xQ_l_s))).*sqrt(g_X2 - g_s_l);
  s_Y1 = sign(real(hh_s_Y./(Y1_l-yQ_l_s))).*sqrt(g_Y1 - g_s_l);
  s_Y2 = sign(real(hh_s_Y./(Y2_l-yQ_l_s))).*sqrt(g_Y2 - g_s_l);

%  s_X1 = sign(real(X1_l-xQ_l_s)).*sqrt((g_X1 - g_s_l));

%  Expressions uniformes des 4 contributions
E_X1 = ...
      Epc1.*sqrt(1/pi).*(ones(3,1)*subfun_Transfun(sqrt(k).*s_X1)) ...
      + ...
     cross(rr_X1, cross(rr_X1, cross(Plaque.N*ones(1,NbP), h_X1))) .* ...
      (ones(3,1)*(Cte*sqrt(2*pi)/k^(3/2)  ...
      ./ sum_r_X1 ./dgdx_X1_l .* exp(-k*g_X1)...
      .* sqrt(det_Qz_fsc_X1./detQ0./QF.Ayy_l))) ...
      - ...
      cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), h_s))) .* ...
      (ones(3,1)*(Cte.*sqrt(pi)/k^(3/2) ...
      ./ sum_r_s  ./s_X1.* exp(-k*g_X1) ...
      .* sqrt(det_Qz_fsc_s./detQ0./detA_l)));
    

E_X2 = ...
      Epc1.*sqrt(1/pi).*(ones(3,1)*(sqrt(pi)-subfun_Transfun(sqrt(k).*s_X2))) ...
      - ...
     cross(rr_X2, cross(rr_X2, cross(Plaque.N*ones(1,NbP), h_X2))) .* ...
      (ones(3,1)*(Cte*sqrt(2*pi)/k^(3/2)  ...
      ./ sum_r_X2 ./dgdx_X2_l .* exp(-k*g_X2)...
      .* sqrt(det_Qz_fsc_X2./detQ0./QF.Ayy_l))) ...
      + ...
      cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), h_s))) .* ...
      (ones(3,1)*(Cte.*sqrt(pi)/k^(3/2) ...
      ./ sum_r_s  ./s_X2.* exp(-k*g_X2) ...
      .* sqrt(det_Qz_fsc_s./detQ0./detA_l)));
    

E_Y1 = ...
      Epc1.*sqrt(1/pi).*(ones(3,1)*subfun_Transfun(sqrt(k).*s_Y1)) ...
      + ...
     cross(rr_Y1, cross(rr_Y1, cross(Plaque.N*ones(1,NbP), h_Y1))) .* ...
      (ones(3,1)*(Cte*sqrt(2*pi)/k^(3/2)  ...
      ./ sum_r_Y1 ./dgdy_Y1_l .* exp(-k*g_Y1)...
      .* sqrt(det_Qz_fsc_Y1./detQ0./QF.Axx_l))) ...
      - ...
      cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), h_s))) .* ...
      (ones(3,1)*(Cte.*sqrt(pi)/k^(3/2) ...
      ./ sum_r_s ./s_Y1 .* exp(-k*g_Y1) ...
      .* sqrt(det_Qz_fsc_s./detQ0./detA_l)));
    

E_Y2 = ...
      Epc1.*sqrt(1/pi).*(ones(3,1)*(sqrt(pi)-subfun_Transfun(sqrt(k).*s_Y2))) ...
      - ...
     cross(rr_Y2, cross(rr_Y2, cross(Plaque.N*ones(1,NbP), h_Y2))) .* ...
      (ones(3,1)*(Cte*sqrt(2*pi)/k^(3/2)  ...
      ./ sum_r_Y2 ./dgdy_Y2_l .* exp(-k*g_Y2)...
      .* sqrt(det_Qz_fsc_Y2./detQ0./QF.Axx_l))) ...
      + ...
      cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), h_s))) .* ...
      (ones(3,1)*(Cte.*sqrt(pi)/k^(3/2) ...
      ./ sum_r_s ./s_Y2 .* exp(-k*g_Y2) ...
      .* sqrt(det_Qz_fsc_s./detQ0./detA_l)));
    



Epc2 = - E_X1 - E_X2 - E_Y1 - E_Y2;

%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  



%  %  =======================================================
%  %  Points critiques du troisieme ordre (coins)
%  %  =======================================================
%  
%  %  points critiques du troisieme ordre
%     rQ_X1Y1 = [X1*ones(1,NbP); Y1*ones(1,NbP); Plaque.Centre(3)*ones(1,NbP)]; 
%   rQ_l_X1Y1 = [X1_l*ones(1,NbP); Y1_l*ones(1,NbP); Plaque.Centre(3)*ones(1,NbP)]; 
%      r_X1Y1 = P_l - rQ_l_X1Y1;
%  sum_r_X1Y1 = sqrt(sum(r_X1Y1.^2));
%     rr_X1Y1 = r_X1Y1./(ones(3,1)*sum_r_X1Y1);
%  
%  %  matrice de courbure evaluee au point critique du troisieme ordre
%  [Qz_X1Y1, det_Qz_fsc_X1Y1] = subfun_Qz(rQ_l_X1Y1(3,:), Champ);  
%  
%  %  fonction de phase evaluee au point critique du troisieme ordre
%  g_X1Y1_l = subfun_g_l(rQ_l_X1Y1(1), rQ_l_X1Y1(2));
%  
%  %  fonction de polarisation dans le repere faisceau incident
%  %  evaluee au point critique du premier ordre
%  h_fsc_X1Y1 = subfun_h(rQ_l_X1Y1, Champ);
%  
%  %  fonction de polarisation dans le repere absolu
%  %  evaluee au point critique du premier ordre
%  h_X1Y1 = R_Rfsc2Rabs*h_fsc_X1Y1;
%  
%  % fonction d'amplitude au point critique du troisieme ordre :
%  ampl_X1Y1 = j*sqrt(epsr) ...
%          .* sqrt(det_Qz_fsc_X1Y1./detQ0./detA_l) ...
%          ./ sum_r_X1Y1;
%  
%  dgdx_X1Y1 = rQ_X1Y1(1,:).*Axx_l + rQ_X1Y1(2,:).*Axy_l - Bx_l ;
%  dgdy_X1Y1 = rQ_X1Y1(2,:).*Ayy_l + rQ_X1Y1(1,:).*Axy_l - By_l ;
%  
%  uu =  sign(real(hh_s_X./(X1_l-xQ_l_s)));
%  vv = sign(real(hh_s_Y./(Y1_l-yQ_l_s)));
%  ww = sign(real(hh_s_X.*hh_s_Y./((X1_l-xQ_l_s).*(Y1_l-yQ_l_s))));
%  u_X1 = uu.*(sqrt(subfun_g_l(X1_l, rQ_l_s(2,:)) - g_s_l));
%  u_Y1 = vv.*(sqrt(subfun_g_l(rQ_l_s(1,:), Y1_l) - g_s_l));
%  
%    %  Constante multiplicative
%    Cte = j*sqrt(epsr)*k/2/pi ;
%  %  %  =====================================================
%  %  %  point critique de premier ordre (point col)
%  %  %  =====================================================
%  %  J1 = cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), h_s))) .* ...
%  %        (ones(3,1)*(ampl_s .* exp(-k*g_s_l) .*Cte*2/k.* ...
%  %          subfun_Transfun(sqrt(k).*u_Y1).*subfun_Transfun(sqrt(k).*u_X1)));
%  %  
%  %  J2 =  - cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), h_s))) .* ...
%  %        (ones(3,1)*(ampl_s .* exp(-k*g_s_l) .* Cte/k^(3/2) .* ...
%  %        subfun_Transfun(sqrt(k).*u_Y1)./u_X1 .* exp(-k.*u_X1.^2) ));
%  %  %  
%  %  %     rQ_X1Ys = [X1*ones(1,NbP); rQ_s(2,:); Plaque.Centre(3)*ones(1,NbP)]; 
%  %  %   rQ_l_X1Ys = [X1_l*ones(1,NbP); rQ_l_s(2,:); Plaque.Centre(3)*ones(1,NbP)]; 
%  %  %      r_X1Ys = P_l - rQ_l_X1Ys;
%  %  %  sum_r_X1Ys = sqrt(sum(r_X1Ys.^2));
%  %  %     rr_X1Ys = r_X1Ys./(ones(3,1)*sum_r_X1Ys);
%  %  %  h_fsc_X1Ys = subfun_h(rQ_l_X1Ys, Champ);
%  %  %      h_X1Ys = R_Rfsc2Rabs*h_fsc_X1Ys;
%  %  %  
%  %  %  u_X1Y1 = uu.*(sqrt(subfun_g_l(X1_l, Y1_l) - g_s_l));
%  %  %  
%  %  %  J3 = cross(rr_X1Ys, cross(rr_X1Ys, cross(Plaque.N*ones(1,NbP), h_X1Ys))) .* ...
%  %  %        (ones(3,1)*(Cte*sqrt(2*pi)/k^(3/2) .* exp(-k.*g_s_l).* exp(-k.*u_X1.^2) .* ...
%  %  %        sqrt(det_Qz_fsc_Y2./detQ0./Axx_l) ./ dgdx_X1_l .* ...
%  %  %        subfun_Transfun(sqrt(k).*u_X1Y1)));
%  %  %  
%  %  %  
%  %  %  
%  %    E_X1Y1 = J1+J2;
%  
%  l1 = Axx_l./2;
%  l2 = detA_l./Axx_l./2;
%  
%  Bu = Bx_l;
%  Bv = By_l - Axy_l./Axx_l.*Bx_l;
%  
%  AA = (-1+erfz(1/2*(Axx_l.*X1_l+Axy_l.*Y1_l-Bx_l).*2^(1/2)./Axx_l.^(1/2))) ...
%    .* ((1/2*exp(-By_l.*Axy_l.*Bx_l./detA_l) ...
%    .* pi ...
%    .* (exp(By_l.*Axy_l.*Bx_l./detA_l) ...
%        .* exp(1/2.*(Bx_l.^2.*Ayy_l + Axx_l.*By_l.^2 -2*By_l.*Axy_l.*Bx_l)./detA_l) ...
%        -  erfz(1/2.*2^(1/2).*(By_l.*Axx_l-Bx_l.*Axy_l-Y1_l.*Ayy_l.*Axx_l+Y1_l.*Axy_l.^2).*sqrt(1./(Axx_l.*detA_l))) ...
%      .*exp(1/2*(Bx_l.^2.*Ayy_l+Axx_l.*By_l.^2)./detA_l)).*sqrt(1./detA_l)));
%  
%  
%  
%  E_X1Y1 = cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), h_s))) .* ...
%          (ones(3,1)*(ampl_s .* exp(-k*C_l) * pi/2/k.*sqrt(1./(l1.*l2)) .* ...
%          exp(k/2.*(Bu.^2./l1+Bv.^2./l2))  .* ...
%          (1-erfz((k*l1.*(X1_l+Axy_l./Axx_l.*Bv./l2)-k*Bu).*sqrt(1./(2*k*l1)))) .* ...
%          (1-erfz((k*l2.*Y1_l - k.*Bv).*sqrt(1./(2*k*l2))))));
%  %  E_X1Y1 = cross(rr_X1Y1, cross(rr_X1Y1, cross(Plaque.N*ones(1,NbP), h_X1Y1))) .* ...
%  %          (ones(3,1)*(ampl_X1Y1 *1/k^(3/2).* exp(-k*C_l) .* AA));

%  E_X1Y1 = integrate_I_XiYj(P, P_l, Champ, Plaque, QF, X1, X2, Y1, Y2, O_l, k);

E_X1Y1 = zeros(size(Epc1));

phi_X1 = 1/2*(-QF.Axy_l.^2./QF.Ayy_l+QF.Axx_l).*X1_l^2 ...
        -1*(QF.Bx_l-QF.By_l.*QF.Axy_l./QF.Ayy_l).*X1_l ...
        -1/2*QF.By_l.^2./QF.Ayy_l+QF.C_l;
phi_X2 = 1/2*(-QF.Axy_l.^2./QF.Ayy_l+QF.Axx_l).*X2_l^2 ...
        -1*(QF.Bx_l-QF.By_l.*QF.Axy_l./QF.Ayy_l).*X2_l ...
        -1/2*QF.By_l.^2./QF.Ayy_l+QF.C_l;

phi_Y1 = 1/2*(-QF.Axy_l.^2./QF.Axx_l+QF.Ayy_l).*Y1_l^2 ...
        -1*(+QF.By_l-QF.Bx_l.*QF.Axy_l./QF.Axx_l).*Y1_l ...
        -1/2*QF.Bx_l.^2./QF.Axx_l+QF.C_l;
phi_Y2 = 1/2*(-QF.Axy_l.^2./QF.Axx_l+QF.Ayy_l).*Y2_l^2 ...
        -1*(+QF.By_l-QF.Bx_l.*QF.Axy_l./QF.Axx_l).*Y2_l ...
        -1/2*QF.Bx_l.^2./QF.Axx_l+QF.C_l;

s_phi_X1 = sign(real(sqrt(2./(detA_l./QF.Ayy_l))./(X1_l - xQ_l_s))) ...
          .* sqrt(phi_X1 - g_s_l);
s_phi_X2 = sign(real(sqrt(2./(detA_l./QF.Ayy_l))./(X2_l - xQ_l_s))) ...
          .* sqrt(phi_X2 - g_s_l);
s_phi_Y1 = sign(real(sqrt(2./(detA_l./QF.Axx_l))./(Y1_l - yQ_l_s))) ...
          .* sqrt(phi_Y1 - g_s_l);
s_phi_Y2 = sign(real(sqrt(2./(detA_l./QF.Axx_l))./(Y2_l - yQ_l_s))) ...
          .* sqrt(phi_Y2 - g_s_l);

E_X1Y1 = (Epc1)./pi.*(ones(3,1) * (...
    subfun_Transfun(sqrt(k).*s_phi_X1) ... 
    .*subfun_Transfun(sqrt(k).*s_phi_Y1) )) ;
E_X2Y1 = (Epc1)./pi.*(ones(3,1) * (...
    (sqrt(pi)-subfun_Transfun(sqrt(k).*s_phi_X2)) ... 
    .*subfun_Transfun(sqrt(k).*s_phi_Y1) )) ;

E_X1Y2 = (Epc1)./pi.*(ones(3,1) * (...
    subfun_Transfun(sqrt(k).*s_phi_X1) ... 
    .*(sqrt(pi) -subfun_Transfun(sqrt(k).*s_phi_Y2)) )) ;
E_X2Y2 = (Epc1)./pi.*(ones(3,1) * (...
    (sqrt(pi)-subfun_Transfun(sqrt(k).*s_phi_X2)) ... 
    .*(sqrt(pi)-subfun_Transfun(sqrt(k).*s_phi_Y2)) )) ;

%  % vecteur distance aux points critiques du troisieme ordre
%     rQ_X1Y1 = [X1*ones(1,Ndg_dy_x_Y_idbP); Y1*ones(1,NbP);zeros(1,NbP)]; 
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
Epc3 = E_X1Y1 + E_X2Y1 + E_X1Y2 + E_X2Y2;

%  astuce pour supprimer la mauvaise évaluation
%  des integrales "de coin" : 
%  on ne conserve que la zone proche du spéculaire
Epc3 = Epc3.*(ones(3,1)*exp(-k*g_s_l));

% ==========================================================
% champ total
% ==========================================================
%  %  filtrage des NaN
%  E_X1(isnan(E_X1)) = eps;
%  E_X2(isnan(E_X2)) = eps;
%  E_Y1(isnan(E_Y1)) = eps;
%  E_Y2(isnan(E_Y2)) = eps;
%  Epc1(isnan(Epc1)) = eps;



E = Epc1+Epc2+Epc3;%
H = zeros(size(E));



%  Fonction de phase
function g = subfun_g_l(xQ_l,yQ_l,QF)
  g = 1/2*(xQ_l.^2.*QF.Axx_l + 2*xQ_l.*yQ_l.*QF.Axy_l + yQ_l.^2.*QF.Ayy_l) ...
      - QF.Bx_l.*xQ_l - QF.By_l.*yQ_l ...
      + QF.C_l;

%  %  Fonction de phase
%  function g = subfun_g_l2(xQ_l,yQ_l,id,QF)
%    g = 1/2*(xQ_l.^2.*QF.Axx_l(id) + 2*xQ_l.*yQ_l.*QF.Axy_l(id) + yQ_l.^2.*QF.Ayy_l(id)) ...
%        - QF.Bx_l(id).*xQ_l - QF.By_l(id).*yQ_l ...
%        + QF.C_l(id);

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




%  ==============================================
%  ==============================================
%  ==============================================
%   integration numerique des integrales de coins
%  ==============================================
%  ==============================================
%  ==============================================
function Inte = integrate_I_XiYj(P, P_l, Champ, Plaque, QF, X1, X2, Y1, Y2, O_l, k)
global R_Rfsc2Rabs epsr;
NbP=size(P,2);
Y_stop = 5*Y1;
X_stop = 5*X1;

N = 20;
W=sqrt(imag(2./k./Champ.Courbure(1)))  ;
if (O_l(1) - N*W) < X1
  BorneInfX1 = X1;
  BorneSupX1 = N*W;
else
  BorneInfX1 = O_l(1) - N*W;
  BorneSupX1 = O_l(1) + N*W;
end


if (O_l(1) + N*W) > X2
  BorneInfX2 = N*X2;
  BorneSupX2 = X2;
else
  BorneInfX2 = O_l(1) - N*W;
  BorneSupX2 = O_l(1) + N*W;
end

Inte = zeros(3,NbP);
%  
%  a = X1;
%  b = X_stop;
%  c = Y1;
%  d = Y_stop;
%  
%  N = 1e3; % au carr�
%  u=a+(b-a)*rand(1,N);
%  v=c+(d-c)*rand(1,N);
%  
%  [uu,vv] = ndgrid(u,v);
%  uu = reshape(uu, 1,N*N);
%  vv = reshape(vv, 1,N*N);

PB = progress('init', ['Sommation sur ', num2str(NbP), ' points entre X=', num2str(X1), ' et ', num2str(X_stop), ' et Y=', num2str(Y1), ' et ', num2str(Y_stop)]);
%  Pour tous les points de visualisation
for id = 1:NbP
%    %  =============================  
%    %  integration purement numerique
%    %  ============================= 
%    Inte(:,id) = dblquad(@(x,y)integrande2D(x,y,P,id,O_l,Champ,Plaque,k,epsr), X, X_stop, Y, Y_stop);

%    %  =============================  
%    %  integration par Monte-Carlo
%    %  ============================= 
%    %  nombre de samples 
%   
%    F = integrande2D(uu,vv,P,id,O_l,Champ,Plaque,k,epsr);
%    Inte(:,id)=(b-a)*(d-c)*mean(F,2);


%    %  =============================  
%    %  utilisation d'algo "rapide"
%    %  =============================
%    Nby = 50;
%    y = linspace(Y1, 3*Y1, Nby);
%      
%    for idy=1:size(y,2)
%      f = @(x)integrande2D(x, y(idy), P(:,id), O_l, Champ, Plaque, getQF(QF,id), k, epsr);
%  %  Integration par la methode des trapezes 
%  %      Inte(:,id) = Inte(:,id) + trapz(x,f(x),2);
%  
%  %  Integration par la methode de la quadrature de Gauss
%  %      Inte(:,id) = Inte(:,id) + gaussquadrature(f, X, X_stop,20);   
%  
%  %  Integration numerique rapide utilisant la quadrature de Clenshaw-Curtis   
%      Inte(:,id) = Inte(:,id) + clenshawcurtisquadrature(f, X1, 3*X1, 100);
%    end
%    Inte(:,id) = Inte(:,id)./Nby;

  %  ===================================
  %  Integration de l'integrale "simple"
  %  ===================================
%    % integration numerique
%    Inte(:,id) = quad(...
%                @(x)integrande1D(x, Y, P(:,id), O_l, Champ, Plaque, getQF(QF,id), k, epsr), ...
%                X, X_stop, 1e-9);

%    %  trapeze
%    x = [X1:2*pi/k/50:3*X1];
%    Inte(:,id) = trapz(x,integrande1D(x, Y1, P(:,id), O_l, Champ, Plaque, getQF(QF,id), k, epsr),2);


  Inte(:,id) = clenshawcurtisquadrature(...
                @(x)integrande_x(x, Y1, P(:,id), O_l, Champ, Plaque, getQF(QF,id), k, epsr, +1), ...
                BorneInfX1, BorneSupX1, N*10)  ...
              + ...
              clenshawcurtisquadrature(...
                @(x)integrande_x(x, Y1, P(:,id), O_l, Champ, Plaque, getQF(QF,id), k, epsr, +1), ...
                BorneInfX2, BorneSupX2, N*10)  ...
              + ...
              clenshawcurtisquadrature(...
                @(x)integrande_x(x, Y2, P(:,id), O_l, Champ, Plaque, getQF(QF,id), k, epsr, -1), ...
                BorneInfX1, BorneSupX1, N*10)  ...
              + ...
              clenshawcurtisquadrature(...
                @(x)integrande_x(x, Y2, P(:,id), O_l, Champ, Plaque, getQF(QF,id), k, epsr, -1), ...
                BorneInfX2, BorneSupX2, N*10);

  % MAJ Progress Bar
  PB = progress(PB, id/NbP);
end




%  %  integrande pour integrale bidim.
%  %  x vecteur, y scalaire
%  %  size(P) = (3,1)
%  function intgd = integrande2D(x,y,P,O_l,Champ,Plaque,QF,k,epsr)
%  global R_Rfsc2Rabs;
%  
%  x_l = x - O_l(1);
%  y_l = y - O_l(2);
%  
%    %  constante amplitude
%    Cte = j*k*sqrt(epsr)./(2*pi);
%  
%    %  vecteur position
%      rQ = [x; y*ones(size(x)); Plaque.Centre(3)*ones(size(x))]; 
%    rQ_l = [x_l; y_l*ones(size(x)); Plaque.Centre(3)*ones(size(x))]; 
%       r = (P-O_l)*ones(size(x)) - rQ_l;
%   sum_r = sqrt(sum(r.^2));
%      rr = r./(ones(3,1)*sum_r);
%  
%    % point dans le repere fsc
%    zQ_fsc = dot([0;0;1]*ones(size(x)), R_Rfsc2Rabs.'*(rQ  - Champ.Centre*ones(size(x))));
%  
%    [Qz, det_Qz] = subfun_Qz(zQ_fsc, Champ);
%    [Q0, det_Q0] = subfun_Qz(0, Champ);
%  
%    %  fonction de polarisation dans le repere faisceau incident
%    %  evaluee au point critique du premier ordre
%    h = R_Rfsc2Rabs*subfun_h(rQ, Champ);
%  
%    
%    g = subfun_g_l(x_l, y_l, QF);
%    f_s = Cte./sum_r.*sqrt(det_Qz./det_Q0);
%    f_v = cross(rr, cross(rr, cross(Plaque.N*ones(size(x)), h))); 
%  
%    intgd = f_v .* (ones(3,1)*(f_s.*exp(-k*g)));



%  integrande pour l'integrale simple de X a X_stop
function intgd = integrande_x(x, Y, P, O_l, Champ, Plaque, QF, k, epsr, signe)
global R_Rfsc2Rabs;
detA_l = QF.Axx_l.*QF.Ayy_l - QF.Axy_l.^2;

P_l = P - O_l;
x_l = x - O_l(1);
Y_l = Y - O_l(2);

NbX = size(x,2);
NbP = size(P,2);
  %  ========================================================
  %  Definition des fonctions
%    %  ========================================================
      ys_l = (QF.By_l - x_l.*QF.Axy_l)./QF.Ayy_l;
      ys = ys_l + O_l(2);
  
  g_x_ys_l = subfun_g_l(x_l, ys_l, QF);
   g_x_Y_l = subfun_g_l(x_l, Y_l, QF);
  
    rQ_x_Y = [x; Y*ones(1,NbX); Plaque.Centre(3)*ones(1,NbX)]; 
  rQ_l_x_Y = [x_l; Y_l*ones(1,NbX); Plaque.Centre(3)*ones(1,NbX)]; 
     r_x_Y = P_l*ones(1,NbX) - rQ_l_x_Y;
 sum_r_x_Y = sqrt(sum(r_x_Y.^2));
    rr_x_Y = r_x_Y./(ones(3,1)*sum_r_x_Y);
  
    rQ_x_ys = [x; ys; Plaque.Centre(3)*ones(1,NbX)]; 
  rQ_l_x_ys = [x_l; ys_l; Plaque.Centre(3)*ones(1,NbX)]; 
     r_x_ys = P_l*ones(1,NbX) - rQ_l_x_ys;
 sum_r_x_ys = sqrt(sum(r_x_ys.^2));
    rr_x_ys = r_x_ys./(ones(3,1)*sum_r_x_ys);

  % point dans le repere fsc
  rQ_fsc_x_ys = R_Rfsc2Rabs.'*(rQ_x_ys  - Champ.Centre*ones(1,NbX));
   rQ_fsc_x_Y = R_Rfsc2Rabs.'*(rQ_x_Y  - Champ.Centre*ones(1,NbX));

    %  matrice de courbure evaluee au point critique du troisieme ordre
    [Q_x_Y, det_Qz_fsc_x_Y] = subfun_Qz(rQ_fsc_x_Y(3,:), Champ);  
  [Q_x_ys, det_Qz_fsc_x_ys] = subfun_Qz(rQ_fsc_x_ys(3,:), Champ);  
      [Q0, det_Q0] = subfun_Qz(0, Champ);

%    aux coordonnees "col" : dx/ds
%    hh_s_X  = sqrt(2./(QF.Axx_l-QF.Axy_l.^2./QF.Ayy_l));
  hh_s_Y  = sqrt(2./(-QF.Axy_l.^2./QF.Axx_l + QF.Ayy_l));
      s_Y = sign(real(hh_s_Y./(Y_l-rQ_l_x_ys(2,:)))).*sqrt(g_x_Y_l - g_x_ys_l);

  %  fonction de polarisation dans le repere faisceau incident
  %  evaluee au point critique du premier ordre
  h_x_Y = R_Rfsc2Rabs*subfun_h(rQ_fsc_x_Y, Champ);
 h_x_ys = R_Rfsc2Rabs*subfun_h(rQ_fsc_x_ys, Champ);

  %  derivee de la fonction de phase au point (x,Y)
  dgdy_x_Y = x_l*QF.Axy_l + Y_l*QF.Ayy_l - QF.By_l;

  if signe == +1
    Q_s = subfun_Transfun(sqrt(k).*s_Y);
  elseif signe == -1
    Q_s = sqrt(pi) - subfun_Transfun(sqrt(k).*s_Y);
  end

  intgd1 = cross(rr_x_ys, cross(rr_x_ys, cross(Plaque.N*ones(1,NbX), h_x_ys))) ...
      .* (ones(3,1)*(j*k*sqrt(epsr)/2/pi ...
      ./ sum_r_x_ys ...
      .* sqrt(det_Qz_fsc_x_Y./det_Q0.*2./QF.Ayy_l/k) ...
      .*  exp(-k*g_x_ys_l).*Q_s));

  intgd2 = cross(rr_x_Y, cross(rr_x_Y, cross(Plaque.N*ones(1,NbX), h_x_Y))) ...
      .* (ones(3,1)*(j*k*sqrt(epsr)/2/pi ....
      ./ sum_r_x_Y ...
      .* sqrt(det_Qz_fsc_x_Y./det_Q0) ...
      ./ k ./dgdy_x_Y ...
      .* exp(-k*g_x_Y_l)));

  intgd3 = - cross(rr_x_ys, cross(rr_x_ys, cross(Plaque.N*ones(1,NbX), h_x_ys))) ...
      .* (ones(3,1)*(j*k*sqrt(epsr)/2/pi ...
      ./ sum_r_x_ys ...
      .* sqrt(det_Qz_fsc_x_Y./det_Q0.*2./QF.Ayy_l) ...
      .* 1./(2*k*s_Y) ...
      .* exp(-k*g_x_Y_l))) ;

  intgd = intgd1 + intgd2 + intgd3;








function QF = getQF(QF, id)
  QF.Axx_l = QF.Axx_l(id);
  QF.Axy_l = QF.Axy_l(id);
  QF.Ayy_l = QF.Ayy_l(id);
   QF.Bx_l = QF.Bx_l(id);
   QF.By_l = QF.By_l(id);
    QF.C_l = QF.C_l(id);
