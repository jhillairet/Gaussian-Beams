function [E,H, Epc1, Epc2, Epc3] = EH_FG_OP_Plaque_conductrice_uniforme(Champ, Plaque, P, f, epsr, varargin)
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
%  
%  ATTENTION : la methode utilisee pour calculer les fontions
%  d'erreur peut faire varier tres sensiblement le resultat
%  (cf fonction F = subfun_Transfun(z) )
%  
%  
%  
debug = false;

%  Accelerer le calcul ? 
%  Pour accelerer le calcul, on regarde ou se situe 
%  le point d'impact du faisceau sur la plaque. 
%  Si celui-ci est loint des bords (à l'interieur)
%  alors on ne prendra en compte que le speculaire.
if debug
  bool_accelerer = true;
else
  bool_accelerer = false;
end

% ===========================================
%      CONSTANTES et PARAMETRES
% ===========================================
     c = 2.997925e8;
    Z0 = 120*pi;
lambda = c/f;
     k = 2*pi/lambda;
   NbP = length(P);






%  Coordonnees des sommets
X1 = max(Plaque.Sommets(1,:));
X2 = min(Plaque.Sommets(1,:));
Y1 = max(Plaque.Sommets(2,:));
Y2 = min(Plaque.Sommets(2,:));

%  Calcul du point d'intersection du faisceau gaussien
%  avec le plan. Ce point correspond ï¿½ l'origine du repere local
%  dans lequel on va faire les calculs.
%  
%  Distance ï¿½ parcourir entre le centre du repere faisceau
%  et le plan :
d =  (Plaque.Centre(3) - Champ.Centre(3))./Champ.ez(3);
%  coordonnees du point d'intersection :
O_l = Champ.Centre + d*Champ.ez;
%  Si l'argument optionnel n'est pas vide, il d'agit du 
%  point d'intersection corrige.
if ~isempty(varargin)
  if ~isempty(varargin{1})
    O_l = (varargin{1})
  end
end

if debug
  %  DEBUG
  %  Si le point d'intersection est en dehors de la plaque
  %  alors on prend comme origine le point de la plaque le 
  %  plus proche
  if ( O_l(1) > X1 ...
    | O_l(1) < X2 ...
    | O_l(2) > Y1 ...
    | O_l(2) < Y2)
    disp('point inter en dehors')
    if O_l(1) > X1
      O_l(1) = X1;
    elseif O_l(1) < X2
      O_l(1) = X2;
    end
    if O_l(2) > Y1
      O_l(2) = Y1;
    elseif O_l(2) < Y2
      O_l(2) = Y2;
    end
  end
end




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
m11 = R_Rabs2Rfsc(1,1);
m12 = R_Rabs2Rfsc(1,2);
m13 = R_Rabs2Rfsc(1,3);
m21 = R_Rabs2Rfsc(2,1);
m22 = R_Rabs2Rfsc(2,2);
m23 = R_Rabs2Rfsc(2,3);
m31 = R_Rabs2Rfsc(3,1);
m32 = R_Rabs2Rfsc(3,2);
m33 = R_Rabs2Rfsc(3,3);

% Coordonnees de P dans le repere fsc
P_l_fsc = R_Rabs2Rfsc*(P_l-Champ.Centre*ones(1,NbP)-O_l*ones(1,NbP));

% centre du faisceau exprimï¿½ dans le repere faisceau
r_l_i_Oi = R_Rabs2Rfsc*(Champ.Centre - O_l);
x_l_i_Oi = r_l_i_Oi(1);
y_l_i_Oi = r_l_i_Oi(2);
z_l_i_Oi = r_l_i_Oi(3);

ri_Oi = R_Rabs2Rfsc*Champ.Centre;
xi_Oi = ri_Oi(1);
yi_Oi = ri_Oi(2);
zi_Oi = ri_Oi(3);

M_fsc = R_Rfsc2Rabs;

%  % elements de la forme quadratique (QF)
%  % dans sa notation matricielle
%  QF.Axx_l = j*(1./R_l - P_l(1,:).^2./R_l.^3 ...
%    + Qp(1).*m11.^2 ...
%    + Qp(2).*m21.^2 ...
%    + 2*Qp(3).*m11.*m21);
%  
%  QF.Ayy_l =  j*(1./R_l - P_l(2,:).^2./R_l.^3 ...
%    + Qp(2).*m22.^2 ...
%    + Qp(1).*m12.^2 ...
%    + 2*Qp(3).*m12.*m22);
%  
%  QF.Axy_l = j*( - P_l(1,:).*P_l(2,:)./R_l.^3 ...
%    + Qp(3).*(m11.*m22 + m12.*m21)/2 ...
%    + Qp(2).*m21.*m22 ...
%    + Qp(1).*m11.*m12);
%  
%  QF.Bx_l = j*(- m31 ...
%    + P_l(1,:)./R_l ...
%    + Qp(3).*(x_l_i_Oi.*m21 + m11.*y_l_i_Oi) ...
%    + Qp(2).*m21.*y_l_i_Oi ...
%    + Qp(1).*m11.*x_l_i_Oi);
%  
%  QF.By_l = j*(- m32 ...
%    + P_l(2,:)./R_l ...
%    + Qp(3).*(x_l_i_Oi.*m22 + m12.*y_l_i_Oi) ...
%    + Qp(2).*m22.*y_l_i_Oi ...
%    + Qp(1).*m12.*x_l_i_Oi);
%  
%  QF.C_l = j*(R_l ...
%    - z_l_i_Oi ...
%    +1/2*Qp(1).*x_l_i_Oi.^2 ...
%    +1/2*Qp(2).*y_l_i_Oi.^2 ...
%    +Qp(3).*x_l_i_Oi.*y_l_i_Oi);


%  ==================================================
%  Forme quadratique avec les notations de la these
%  ==================================================
alpha1 = m11*(O_l(1) - Champ.Centre(1)) + m12*(O_l(2) - Champ.Centre(2)) - m13*Champ.Centre(3);
alpha2 = m21*(O_l(1) - Champ.Centre(1)) + m22*(O_l(2) - Champ.Centre(2)) - m23*Champ.Centre(3);
alpha3 = m31*(O_l(1) - Champ.Centre(1)) + m32*(O_l(2) - Champ.Centre(2)) - m33*Champ.Centre(3);

Mij = [m11, m12; m21, m22];
QO_l= [Qp(1), Qp(3); Qp(3), Qp(2)];
tMijQ0_lMij = transpose(Mij)*QO_l*Mij;
tMijQO_lalpha = transpose(Mij)*QO_l*[alpha1;alpha2];

% on rajoute une partie imaginaire (fictive)
% pour eviter la division par zero
R_l = R_l+j*eps; 

P11 = 1./R_l - P_l(1,:).^2./R_l.^3;
P22 = 1./R_l - P_l(2,:).^2./R_l.^3;
P12 = - P_l(1,:).*P_l(2,:)./R_l.^3;

QF.Axx_l = j*(P11 + tMijQ0_lMij(1,1));
QF.Ayy_l = j*(P22 + tMijQ0_lMij(2,2));
QF.Axy_l = j*(P12 + tMijQ0_lMij(1,2));

QF.Bx_l = j*(P_l(1,:)./R_l - tMijQO_lalpha(1) - m31);
QF.By_l = j*(P_l(2,:)./R_l - tMijQO_lalpha(2) - m32);

QF.C_l = j*(R_l + alpha3 + 1/2*transpose([alpha1;alpha2])*QO_l*[alpha1;alpha2]);


%  %  ====================================================
%  %  Forme developpee avec MAPLE (resultat egal)
%  %  ====================================================
%  QF.Axx_l = j./R_l-j.*P_l(1,:).^2./R_l.^3+j.*m21.^2.*Qp(2)+j.*m11.^2.*Qp(1)+2*j.*m11.*m21.*Qp(3);
%  
%  QF.Ayy_l = 2*j.*m12.*m22.*Qp(3)+j.*m12.^2.*Qp(1)+j.*m22.^2.*Qp(2)+j./R_l-j.*P_l(2,:).^2./R_l.^3;
%  
%  QF.Axy_l = j.*Qp(2).*m21.*m22+j.*Qp(3).*m12.*m21+j.*Qp(1).*m11.*m12-j.*P_l(2,:).*P_l(1,:)./R_l.^3+j.*Qp(3).*m11.*m22;
%  
%  QF.Bx_l = -( j.*Qp(2).*m21.*m22.*O_l(2)-j.*Qp(2).*m21.*m22.*Champ.Centre(2)-j.*Qp(2).*m21.^2.*Champ.Centre(1)-j.*Qp(1).*m11.*m12.*Champ.Centre(2)+2.*j.*Qp(3).*m11.*m21.*O_l(1)-j.*Qp(3).*m12.*Champ.Centre(2).*m21+j.*Qp(1).*m11.*m12.*O_l(2)-j.*Qp(3).*m13.*Champ.Centre(3).*m21+j.*m31+j.*Qp(1).*m11.^2.*O_l(1)+j.*Qp(2).*m21.^2.*O_l(1)-j.*Qp(1).*m11.*m13.*Champ.Centre(3)+j.*Qp(3).*m12.*O_l(2).*m21-j.*Qp(2).*m21.*m23.*Champ.Centre(3)+j.*Qp(3).*m11.*m22.*O_l(2)-j.*P_l(1,:)./R_l-2.*j.*Qp(3).*m11.*m21.*Champ.Centre(1)-j.*Qp(1).*m11.^2.*Champ.Centre(1)-j.*Qp(3).*m11.*m22.*Champ.Centre(2)-j.*Qp(3).*m11.*m23.*Champ.Centre(3));
%  
%  QF.By_l = -( -j.*Qp(2).*m21.*Champ.Centre(1).*m22-j.*Qp(3).*m13.*Champ.Centre(3).*m22+j.*Qp(2).*m22.^2.*O_l(2)+j.*Qp(1).*m12.^2.*O_l(2)+j.*Qp(3).*m12.*m21.*O_l(1)-j.*Qp(1).*m12.^2.*Champ.Centre(2)+j.*Qp(1).*m11.*O_l(1).*m12-j.*Qp(3).*m11.*Champ.Centre(1).*m22-j.*Qp(1).*m11.*Champ.Centre(1).*m12-j.*Qp(2).*m22.^2.*Champ.Centre(2)+j.*Qp(3).*m11.*O_l(1).*m22+j.*Qp(2).*m21.*O_l(1).*m22-j.*Qp(3).*m12.*m21.*Champ.Centre(1)-j.*P_l(2,:)./R_l-j.*Qp(1).*m12.*m13.*Champ.Centre(3)-2.*j.*Qp(3).*m12.*m22.*Champ.Centre(2)-j.*Qp(2).*m22.*m23.*Champ.Centre(3)+j.*m32+2.*j.*Qp(3).*m12.*m22.*O_l(2)-j.*Qp(3).*m12.*m23.*Champ.Centre(3));
%  
%  QF.C_l = j.*R_l+j.*(m31.*(O_l(1)-Champ.Centre(1))+m32.*(O_l(2)-Champ.Centre(2))-m33.*Champ.Centre(3))+1./2.*j.*((m11.*(O_l(1)-Champ.Centre(1))+m12.*(O_l(2)-Champ.Centre(2))-m13.*Champ.Centre(3)).^2.*Qp(1)+(m21.*(O_l(1)-Champ.Centre(1))+m22.*(O_l(2)-Champ.Centre(2))-m23.*Champ.Centre(3)).^2.*Qp(2)+2.*(m11.*(O_l(1)-Champ.Centre(1))+m12.*(O_l(2)-Champ.Centre(2))-m13.*Champ.Centre(3)).*(m21.*(O_l(1)-Champ.Centre(1))+m22.*(O_l(2)-Champ.Centre(2))-m23.*Champ.Centre(3)).*Qp(3));


detA_l = QF.Axx_l.*QF.Ayy_l - QF.Axy_l.^2;
iA = [QF.Ayy_l;QF.Axx_l;-QF.Axy_l]./(ones(3,1)*detA_l);

% ======================================================
% Calcul du point col (point critique du premier ordre)
% ======================================================
%  point col dans le repere absolu modifie
xQ_l_s = iA(1,:).*QF.Bx_l + iA(3,:).*QF.By_l;
yQ_l_s = iA(3,:).*QF.Bx_l + iA(2,:).*QF.By_l;
rQ_l_s = [xQ_l_s; yQ_l_s; Plaque.Centre(3)*ones(1,NbP)];
%  point col dans le repere absolu 
xQ_s = xQ_l_s + O_l(1);
yQ_s = yQ_l_s + O_l(2);
rQ_s = [xQ_s; yQ_s; Plaque.Centre(3)*ones(1,NbP)];

% point col dans le repere fsc
rQ_fsc_s = R_Rabs2Rfsc*(rQ_s  - Champ.Centre*ones(1,NbP) );

% vecteur distance au point col
    r_s = P_l - rQ_l_s;
sum_r_s = sqrt(sum((r_s).^2)) + j*eps;
   rr_s = r_s./(ones(3,1)*sum_r_s);

%  matrice de courbure evaluee au point critique du premier ordre
[Qz_s, det_Qz_fsc_s] = subfun_Qz(rQ_fsc_s(3,:), Champ);  
%  meme dans les fonctions d'amplitudes
%  on fait l'approximation Qz = Qp
%  cela permet d'eviter une singularite contenue dans det_Qz_fsc_s
det_Qz_fsc_s = det_Qp;

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


if bool_accelerer
  DD = 10*lambda;

  test =  O_l(1) < X1-DD ...
        & O_l(1) > X2+DD ...
        & O_l(2) < Y1-DD ...
        & O_l(2) > Y2+DD;
  
  if test
    disp('Impact dans la zone centrale --> Diffraction negligee')
    E = Epc1;
    H = zeros(3, size(E,2));
    Epc2 = zeros(3, size(E,2));
    Epc3 = zeros(3, size(E,2));
    return;
  end
end
  
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
  sum_r_X1 = sqrt(sum((r_X1.^2))) + j*eps;
    rr_X1 = r_X1./(ones(3,1)*sum_r_X1);
  
    rQ_X2 = [X2*ones(1,NbP); y_X2;Plaque.Centre(3)*ones(1,NbP)]; 
  rQ_X2_l = [X2_l*ones(1,NbP); y_X2_l;Plaque.Centre(3)*ones(1,NbP)]; 
      r_X2 = P_l - rQ_X2_l;
  sum_r_X2 = sqrt(sum((r_X2.^2))) + j*eps;
    rr_X2 = r_X2./(ones(3,1)*sum_r_X2);
  
    rQ_Y1 = [x_Y1; Y1*ones(1,NbP); Plaque.Centre(3)*ones(1,NbP)]; 
  rQ_Y1_l = [x_Y1_l; Y1_l*ones(1,NbP); Plaque.Centre(3)*ones(1,NbP)]; 
      r_Y1 = P_l - rQ_Y1_l;
  sum_r_Y1 = sqrt(sum((r_Y1.^2))) + j*eps;
    rr_Y1 = r_Y1./(ones(3,1)*sum_r_Y1);
  
    rQ_Y2 = [x_Y2; Y2*ones(1,NbP); Plaque.Centre(3)*ones(1,NbP)]; 
  rQ_Y2_l = [x_Y2_l; Y2_l*ones(1,NbP); Plaque.Centre(3)*ones(1,NbP)]; 
      r_Y2 = P_l - rQ_Y2_l;
  sum_r_Y2 = sqrt(sum((r_Y2).^2)) + j*eps;
    rr_Y2 = r_Y2./(ones(3,1)*sum_r_Y2);

%  %  affichage des points critiques
%  figure(1)
%  hold on
%    line(real(rQ_X1(1,:))/lambda, real(rQ_X1(2,:))/lambda, real(rQ_X1(3,:))/lambda, 'Color', 'b');
%    line(real(rQ_X2(1,:))/lambda, real(rQ_X2(2,:))/lambda, real(rQ_X2(3,:))/lambda, 'Color', 'b');
%    line(real(rQ_Y1(1,:))/lambda, real(rQ_Y1(2,:))/lambda, real(rQ_Y1(3,:))/lambda, 'Color', 'r');
%    line(real(rQ_Y2(1,:))/lambda, real(rQ_Y2(2,:))/lambda, real(rQ_Y2(3,:))/lambda, 'Color', 'r');
%  hold off


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
  %  meme remaque que pour le point col
  det_Qz_fsc_X1 = det_Qp;
  det_Qz_fsc_X2 = det_Qp;
  det_Qz_fsc_Y1 = det_Qp;
  det_Qz_fsc_Y2 = det_Qp;

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
  d2gdx2 = QF.Axx_l-QF.Axy_l.^2./QF.Ayy_l;
  d2gdy2 = -QF.Axy_l.^2./QF.Axx_l+QF.Ayy_l;

  hh_s_X  = sqrt(2./d2gdx2);
  hh_s_Y  = sqrt(2./d2gdy2);


  %  Le signe de la racine carrï¿½e est dï¿½fini tel que l'on doit avoir
  %  h_x -> h_s lorsque s -> 0, cad lq 
%      s_X1 = sign(real(hh_s_X./(X1_l-xQ_l_s))).*sqrt(g_X1 - g_s_l);
%      s_X2 = sign(real(hh_s_X./(X2_l-xQ_l_s))).*sqrt(g_X2 - g_s_l);
%      s_Y1 = sign(real(hh_s_Y./(Y1_l-yQ_l_s))).*sqrt(g_Y1 - g_s_l);
%      s_Y2 = sign(real(hh_s_Y./(Y2_l-yQ_l_s))).*sqrt(g_Y2 - g_s_l);
  %  Methode Felsen
  s_X1 = exp(j*angle((X1_l-xQ_l_s).*sqrt(d2gdx2/2))).*sqrt(abs(g_X1 - g_s_l));
  s_X2 = exp(j*angle((X2_l-xQ_l_s).*sqrt(d2gdx2/2))).*sqrt(abs(g_X2 - g_s_l));
  s_Y1 = exp(j*angle((Y1_l-yQ_l_s).*sqrt(d2gdy2/2))).*sqrt(abs(g_Y1 - g_s_l));
  s_Y2 = exp(j*angle((Y2_l-yQ_l_s).*sqrt(d2gdy2/2))).*sqrt(abs(g_Y2 - g_s_l));

%    %  Methode James
%    s_X1 = exp(j*angle(dgdx_X1_l.*sqrt(2./d2gdx2))).*sqrt(abs(g_X1 - g_s_l));
%    s_X2 = exp(j*angle(dgdx_X2_l.*sqrt(2./d2gdx2))).*sqrt(abs(g_X2 - g_s_l));
%    s_Y1 = exp(j*angle(dgdy_Y1_l.*sqrt(2./d2gdy2))).*sqrt(abs(g_Y1 - g_s_l));
%    s_Y2 = exp(j*angle(dgdy_Y2_l.*sqrt(2./d2gdy2))).*sqrt(abs(g_Y2 - g_s_l));

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
    

%  on filtre les NaN
E_X1(find(isnan(E_X1)))=0;
E_X2(find(isnan(E_X2)))=0;
E_Y1(find(isnan(E_Y1)))=0;
E_Y2(find(isnan(E_Y2)))=0;

Epc2 = - E_X1 - E_X2 - E_Y1 - E_Y2;
 

%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  



% =======================================================
% Integrale du troisieme type
% =======================================================
%  phi_X1 = 1/2*(-QF.Axy_l.^2./QF.Ayy_l+QF.Axx_l).*X1_l^2 ...
%          -1*(QF.Bx_l-QF.By_l.*QF.Axy_l./QF.Ayy_l).*X1_l ...
%          -1/2*QF.By_l.^2./QF.Ayy_l+QF.C_l;
%  phi_X2 = 1/2*(-QF.Axy_l.^2./QF.Ayy_l+QF.Axx_l).*X2_l^2 ...
%          -1*(QF.Bx_l-QF.By_l.*QF.Axy_l./QF.Ayy_l).*X2_l ...
%          -1/2*QF.By_l.^2./QF.Ayy_l+QF.C_l;
%  
%  phi_Y1 = 1/2*(-QF.Axy_l.^2./QF.Axx_l+QF.Ayy_l).*Y1_l^2 ...
%          -1*(+QF.By_l-QF.Bx_l.*QF.Axy_l./QF.Axx_l).*Y1_l ...
%          -1/2*QF.Bx_l.^2./QF.Axx_l+QF.C_l;
%  phi_Y2 = 1/2*(-QF.Axy_l.^2./QF.Axx_l+QF.Ayy_l).*Y2_l^2 ...
%          -1*(+QF.By_l-QF.Bx_l.*QF.Axy_l./QF.Axx_l).*Y2_l ...
%          -1/2*QF.Bx_l.^2./QF.Axx_l+QF.C_l;
%  
%  s_phi_X1 = sign(real(sqrt(2./(detA_l./QF.Ayy_l))./(X1_l - xQ_l_s))) ...
%            .* sqrt((phi_X1 - g_s_l));
%  s_phi_X2 = sign(real(sqrt(2./(detA_l./QF.Ayy_l))./(X2_l - xQ_l_s))) ...
%            .* sqrt((phi_X2 - g_s_l));
%  s_phi_Y1 = sign(real(sqrt(2./(detA_l./QF.Axx_l))./(Y1_l - yQ_l_s))) ...
%            .* sqrt((phi_Y1 - g_s_l));
%  s_phi_Y2 = sign(real(sqrt(2./(detA_l./QF.Axx_l))./(Y2_l - yQ_l_s))) ...
%            .* sqrt((phi_Y2 - g_s_l));

%  E_X1Y1 = (Epc1)./pi.*(ones(3,1) * (...
%      subfun_Transfun(sqrt(k).*s_phi_X1) ... 
%      .*subfun_Transfun(sqrt(k).*s_phi_Y1)  )) ;
%  
%  E_X2Y1 = (Epc1)./pi.*(ones(3,1) * (...
%      (sqrt(pi)-subfun_Transfun(sqrt(k).*s_phi_X2)) ... 
%      .*subfun_Transfun(sqrt(k).*s_phi_Y1) )) ;
%  
%  E_X1Y2 = (Epc1)./pi.*(ones(3,1) * (...
%      subfun_Transfun(sqrt(k).*s_phi_X1) ... 
%      .*(sqrt(pi) -subfun_Transfun(sqrt(k).*s_phi_Y2)) )) ;
%  E_X2Y2 = (Epc1)./pi.*(ones(3,1) * (...
%      (sqrt(pi)-subfun_Transfun(sqrt(k).*s_phi_X2)) ... 
%      .*(sqrt(pi)-subfun_Transfun(sqrt(k).*s_phi_Y2)) )) ;


% vecteur distance aux coins
    rQ_X1Y1 = [X1; Y1;0]; 
  rQ_X1Y1_l = [X1_l; Y1_l;0]; 
     r_X1Y1 = P_l - rQ_X1Y1_l*ones(1,NbP);
 sum_r_X1Y1 = sqrt(sum((r_X1Y1.^2))) + j*eps;
    rr_X1Y1 = r_X1Y1./(ones(3,1)*sum_r_X1Y1);
%  dans le repere du faisceau incident
    rQ_fsc_X1Y1 = R_Rabs2Rfsc*(rQ_X1Y1  - Champ.Centre)*ones(1,NbP);

    rQ_X2Y1 = [X2; Y1;0]; 
  rQ_X2Y1_l = [X2_l; Y1_l;0]; 
     r_X2Y1 = P_l - rQ_X2Y1_l*ones(1,NbP);
 sum_r_X2Y1 = sqrt(sum((r_X2Y1.^2))) + j*eps;
    rr_X2Y1 = r_X2Y1./(ones(3,1)*sum_r_X2Y1);
%  dans le repere du faisceau incident
    rQ_fsc_X2Y1 = R_Rabs2Rfsc*(rQ_X2Y1  - Champ.Centre)*ones(1,NbP);

    rQ_X1Y2 = [X1; Y2;0]; 
  rQ_X1Y2_l = [X1_l; Y2_l;0]; 
     r_X1Y2 = P_l - rQ_X1Y2_l*ones(1,NbP);
 sum_r_X1Y2 = sqrt(sum((r_X1Y2.^2))) + j*eps;
    rr_X1Y2 = r_X1Y2./(ones(3,1)*sum_r_X1Y2);
%  dans le repere du faisceau incident
    rQ_fsc_X1Y2 = R_Rabs2Rfsc*(rQ_X1Y2  - Champ.Centre)*ones(1,NbP);

    rQ_X2Y2 = [X2; Y2;0]; 
  rQ_X2Y2_l = [X2_l; Y2_l;0]; 
     r_X2Y2 = P_l - rQ_X2Y2_l*ones(1,NbP);
 sum_r_X2Y2 = sqrt(sum((r_X2Y2.^2))) + j*eps;
    rr_X2Y2 = r_X2Y2./(ones(3,1)*sum_r_X2Y2);
%  dans le repere du faisceau incident
    rQ_fsc_X2Y2 = R_Rabs2Rfsc*(rQ_X2Y2  - Champ.Centre)*ones(1,NbP);

%  matrices de courbure
%  [Qz_fsc_X1Y1, det_Qz_fsc_X1Y1] = subfun_Qz(rQ_fsc_X1Y1(3,:), Champ);
%  [Qz_fsc_X2Y1, det_Qz_fsc_X2Y1] = subfun_Qz(rQ_fsc_X2Y1(3,:), Champ);
%  [Qz_fsc_X1Y2, det_Qz_fsc_X1Y2] = subfun_Qz(rQ_fsc_X1Y2(3,:), Champ);
%  [Qz_fsc_X2Y2, det_Qz_fsc_X2Y2] = subfun_Qz(rQ_fsc_X2Y2(3,:), Champ);
%  meme remarque que pour le point col
[Qz_fsc_X1Y1, det_Qz_fsc_X1Y1] = subfun_Qz(d, Champ);
[Qz_fsc_X2Y1, det_Qz_fsc_X2Y1] = subfun_Qz(d, Champ);
[Qz_fsc_X1Y2, det_Qz_fsc_X1Y2] = subfun_Qz(d, Champ);
[Qz_fsc_X2Y2, det_Qz_fsc_X2Y2] = subfun_Qz(d, Champ);
det_Qz_fsc_X1Y1 = det_Qp;
det_Qz_fsc_X2Y1 = det_Qp;
det_Qz_fsc_X1Y2 = det_Qp;
det_Qz_fsc_X2Y2 = det_Qp;

%  fonction d'amplitude aux coins
ampl_X1Y1 =  j*sqrt(epsr) ...
        .* sqrt(det_Qz_fsc_X1Y1./detQ0) ...
        ./ sum_r_X1Y1;
ampl_X2Y1 =  j*sqrt(epsr) ...
        .* sqrt(det_Qz_fsc_X2Y1./detQ0) ...
        ./ sum_r_X2Y1;
ampl_X1Y2 =  j*sqrt(epsr) ...
        .* sqrt(det_Qz_fsc_X1Y2./detQ0) ...
        ./ sum_r_X1Y2;
ampl_X2Y2 =  j*sqrt(epsr) ...
        .* sqrt(det_Qz_fsc_X2Y2./detQ0) ...
        ./ sum_r_X2Y2;
%  Fonctions de phase aux coins
g_X1Y1 = subfun_g_l(rQ_X1Y1_l(1,:), rQ_X1Y1_l(2,:), QF);
g_X2Y1 = subfun_g_l(rQ_X2Y1_l(1,:), rQ_X2Y1_l(2,:), QF);
g_X1Y2 = subfun_g_l(rQ_X1Y2_l(1,:), rQ_X1Y2_l(2,:), QF);
g_X2Y2 = subfun_g_l(rQ_X2Y2_l(1,:), rQ_X2Y2_l(2,:), QF);
%  Fonctions de polarisations aux coins
h_X1Y1 = M_fsc*subfun_h(rQ_fsc_X1Y1, Champ);
h_X2Y1 = M_fsc*subfun_h(rQ_fsc_X2Y1, Champ);
h_X1Y2 = M_fsc*subfun_h(rQ_fsc_X1Y2, Champ);
h_X2Y2 = M_fsc*subfun_h(rQ_fsc_X2Y2, Champ);
%  derivee simples aux coins
dgdx_X1Y1_l = QF.Axx_l.*X1_l + QF.Axy_l.*Y1_l - QF.Bx_l;
dgdy_X1Y1_l = QF.Axy_l.*X1_l + QF.Ayy_l.*Y1_l - QF.By_l;

dgdx_X2Y1_l = QF.Axx_l.*X2_l + QF.Axy_l.*Y1_l - QF.Bx_l;
dgdy_X2Y1_l = QF.Axy_l.*X2_l + QF.Ayy_l.*Y1_l - QF.By_l;

dgdx_X1Y2_l = QF.Axx_l.*X1_l + QF.Axy_l.*Y2_l - QF.Bx_l;
dgdy_X1Y2_l = QF.Axy_l.*X1_l + QF.Ayy_l.*Y2_l - QF.By_l;

dgdx_X2Y2_l = QF.Axx_l.*X2_l + QF.Axy_l.*Y2_l - QF.Bx_l;
dgdy_X2Y2_l = QF.Axy_l.*X2_l + QF.Ayy_l.*Y2_l - QF.By_l;

  %  terme de changement de variable dans le passage 
  %  aux coordonnees "col" : dx/ds
  hh_s_X  = sqrt(2./(QF.Axx_l-QF.Axy_l.^2./QF.Ayy_l));
  hh_s_Y  = sqrt(2./(-QF.Axy_l.^2./QF.Axx_l+QF.Ayy_l));


%  nu_X1Y1 = dgdx_X1Y1_l.*sqrt(k./(2.*QF.Axx_l));
%  mu_X1Y1 = dgdy_X1Y1_l.*sqrt(k./(2.*QF.Ayy_l));
%  
%  nu_X2Y1 = dgdx_X2Y1_l.*sqrt(k./(2.*QF.Axx_l));
%  mu_X2Y1 = dgdy_X2Y1_l.*sqrt(k./(2.*QF.Ayy_l));
%  
%  nu_X1Y2 = dgdx_X1Y2_l.*sqrt(k./(2.*QF.Axx_l));
%  mu_X1Y2 = dgdy_X1Y2_l.*sqrt(k./(2.*QF.Ayy_l));
%  
%  nu_X2Y2 = dgdx_X2Y2_l.*sqrt(k./(2.*QF.Axx_l));
%  mu_X2Y2 = dgdy_X2Y2_l.*sqrt(k./(2.*QF.Ayy_l));
%  
%  
%  nu_X1Y1 = dgdx_X1Y1_l.*sqrt(k./(2.*d2gdx2));
%  mu_X1Y1 = dgdy_X1Y1_l.*sqrt(k./(2.*d2gdy2));
%  
%  nu_X2Y1 = dgdx_X2Y1_l.*sqrt(k./(2.*d2gdx2));
%  mu_X2Y1 = dgdy_X2Y1_l.*sqrt(k./(2.*d2gdy2));
%  
%  nu_X1Y2 = dgdx_X1Y2_l.*sqrt(k./(2.*d2gdx2));
%  mu_X1Y2 = dgdy_X1Y2_l.*sqrt(k./(2.*d2gdy2));
%  
%  nu_X2Y2 = dgdx_X2Y2_l.*sqrt(k./(2.*d2gdx2));
%  mu_X2Y2 = dgdy_X2Y2_l.*sqrt(k./(2.*d2gdy2));

%  Methode Felsen
g_X1 = 1/2*QF.Axx_l.*X1_l.^2 - QF.Bx_l.*X1_l;
g_X2 = 1/2*QF.Axx_l.*X2_l.^2 - QF.Bx_l.*X2_l;
g_Y1 = 1/2*QF.Ayy_l.*Y1_l.^2 - QF.By_l.*Y1_l;
g_Y2 = 1/2*QF.Ayy_l.*Y2_l.^2 - QF.By_l.*Y2_l;

g_sx = 1/2*QF.Axx_l.*xQ_l_s.^2 - QF.Bx_l.*xQ_l_s;
g_sy = 1/2*QF.Ayy_l.*yQ_l_s.^2 - QF.By_l.*yQ_l_s;

s_X1 = exp(j*angle((X1_l-xQ_l_s).*sqrt(QF.Axx_l/2))).*sqrt(abs(g_X1 - g_sx));
s_Y1 = exp(j*angle((Y1_l-yQ_l_s).*sqrt(QF.Ayy_l/2))).*sqrt(abs(g_Y1 - g_sy));
s_X2 = exp(j*angle((X2_l-xQ_l_s).*sqrt(QF.Axx_l/2))).*sqrt(abs(g_X2 - g_sx));
s_Y2 = exp(j*angle((Y2_l-yQ_l_s).*sqrt(QF.Ayy_l/2))).*sqrt(abs(g_Y2 - g_sy));

nu_X1Y1 = sqrt(k).*s_X1;
mu_X1Y1 = sqrt(k).*s_Y1;

nu_X2Y1 = sqrt(k).*s_X2;
mu_X2Y1 = sqrt(k).*s_Y1;

nu_X1Y2 = sqrt(k).*s_X1;
mu_X1Y2 = sqrt(k).*s_Y2;

nu_X2Y2 = sqrt(k).*s_X2;
mu_X2Y2 = sqrt(k).*s_Y2;



A11A22 = QF.Axx_l.*QF.Ayy_l;
%  A11A22 = detA_l;

E_X1Y1 = cross(rr_X1Y1, cross(rr_X1Y1, cross(Plaque.N*ones(1,NbP), h_X1Y1))) .* ...
      (ones(3,1)*( 1/pi*ampl_X1Y1 .* exp(-k*g_X1Y1) ...
                  .* sqrt(1./(A11A22)) ...
                  .* subfun_Transfun(nu_X1Y1) ...
                  .* subfun_Transfun(mu_X1Y1) ...
                  .* exp(nu_X1Y1.^2) ...
                  .* exp(mu_X1Y1.^2) ));

E_X2Y1 = cross(rr_X2Y1, cross(rr_X2Y1, cross(Plaque.N*ones(1,NbP), h_X2Y1))) .* ...
      (ones(3,1)*( 1/pi*ampl_X2Y1 .* exp(-k*g_X2Y1) ...
                  .* sqrt(1./(A11A22)) ...
                  .* (sqrt(pi) - subfun_Transfun(nu_X2Y1)) ...
                  .* subfun_Transfun(mu_X2Y1) ...
                  .* exp(nu_X2Y1.^2) ...
                  .* exp(mu_X2Y1.^2) ));

E_X1Y2 = cross(rr_X1Y2, cross(rr_X1Y2, cross(Plaque.N*ones(1,NbP), h_X1Y2))) .* ...
      (ones(3,1)*( 1/pi*ampl_X1Y2 .* exp(-k*g_X1Y2) ...
                  .* sqrt(1./(A11A22)) ...
                  .* subfun_Transfun(nu_X1Y2) ...
                  .* (sqrt(pi) - subfun_Transfun(mu_X1Y2)) ...
                  .* exp(nu_X1Y2.^2) ...
                  .* exp(mu_X1Y2.^2) ));

E_X2Y2 = cross(rr_X2Y2, cross(rr_X2Y2, cross(Plaque.N*ones(1,NbP), h_X2Y2))) .* ...
      (ones(3,1)*( 1/pi*ampl_X2Y1 .* exp(-k*g_X2Y2) ...
                  .* sqrt(1./(A11A22)) ...
                  .* (sqrt(pi) - subfun_Transfun(nu_X2Y2)) ...
                  .* (sqrt(pi) - subfun_Transfun(mu_X2Y2)) ...
                  .* exp(nu_X2Y2.^2) ...
                  .* exp(mu_X2Y2.^2) ));


%  on filtre les NaN
E_X1Y1(find(isnan(E_X1Y1)))=0;
E_X1Y2(find(isnan(E_X1Y2)))=0;
E_X2Y1(find(isnan(E_X2Y1)))=0;
E_X2Y2(find(isnan(E_X2Y2)))=0;


Epc3 = + E_X1Y1 + E_X2Y2 + E_X2Y1 + E_X1Y2 ;



%  figure
%    plot(P(1,:), angle(([Epc1(1,:)+Epc2(1,:);Epc3(1,:)])))

% ==========================================================
% champ total
% ==========================================================

E = Epc1+Epc2+Epc3;
H = zeros(size(E));

%  if isnan(sum(E(:)))
%    dbstop;
%  end


%  Fonction de phase
function g = subfun_g_l(xQ_l,yQ_l,QF)
  g = 1/2*(xQ_l.^2.*QF.Axx_l + 2*xQ_l.*yQ_l.*QF.Axy_l + yQ_l.^2.*QF.Ayy_l) ...
      - QF.Bx_l.*xQ_l - QF.By_l.*yQ_l ...
      + QF.C_l;

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
  h_z = Champ.Coefficients(1).*Qy - Champ.Coefficients(2).*Qx;
    h = [h_x;h_y;h_z];


function F = subfun_Transfun(z)
%  Fonction de transition pour dev. asChamp.Centre(1)totique uniforme

% sans approx 
cerr_fun = 1-erfz(z);

%  %  Approximation legere
%  cerr_fun = (abs(z)<20).*(1-erfz(z)) + ...
%             (abs(z)>20).*exp(-z.^2)/sqrt(pi).*(1./z-1/2./z.^3);
%  %  Approximation grossiere (uniquement grand arguments)
%  cerr_fun = exp(-z.^2)/sqrt(pi).*(1./z-1/2./z.^3);


%  cerr_fun = exp(-z.^2).*W(i*z);
%  cerr_fun = exp(-z.^2).*faddeeva(i*z);

F = sqrt(pi)/2.*cerr_fun ;

