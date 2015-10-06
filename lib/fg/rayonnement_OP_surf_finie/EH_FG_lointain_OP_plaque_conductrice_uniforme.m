function [E,H, Epc1, Epc2, Epc3] = EH_FG_lointain_OP_Plaque_conductrice_uniforme(Champ, Plaque, P, f, epsr, varargin)
% Calcul du champ EM rayonne par une Plaque conductrice 
% dans l'hypothese de l'optique physique
% par un faisceau gaussien approximation champ lointain
% 
% Les effets de la finitude sont pris en compte
% Ainsi que les effets "de coin" de maniere uniforme
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

debug = false;

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
%  coordonnees du point d'intersection
R_inter_Ofg = Champ.Centre + d*Champ.ez
%  Si l'argument optionnel n'est pas vide, il d'agit du 
%  point d'intersection corrige.
if ~isempty(varargin)
  R_inter_Ofg = (varargin{1})
end

if debug
  %  Si le point d'intersection en en dehors de la plaque,
  %  alors on calcule le point le sur le bord de la plaque 
  %  qui est le plus proche
  %  DEBUG
  if ( R_inter_Ofg(1) > X1 ...
    | R_inter_Ofg(1) < X2 ...
    | R_inter_Ofg(2) > Y1 ...
    | R_inter_Ofg(2) < Y2)
    disp('point inter en dehors')
    if R_inter_Ofg(1) > X1
      R_inter_Ofg(1) = X1;
    elseif R_inter_Ofg(1) < X2
      R_inter_Ofg(1) = X2;
    end
  
    if R_inter_Ofg(2) > Y1
      R_inter_Ofg(2) = Y1;
    elseif R_inter_Ofg(2) < Y2
      R_inter_Ofg(2) = Y2;
    end
  end
end

%  Points dans l'origine translatee
 P_l = P  - R_inter_Ofg*ones(1,NbP);
X1_l = X1 - R_inter_Ofg(1);
X2_l = X2 - R_inter_Ofg(1);
Y1_l = Y1 - R_inter_Ofg(2);
Y2_l = Y2 - R_inter_Ofg(2);
%  Distance de l'origine au point obs
   R = sqrt(sum(P.^2));
  rr = P./(ones(3,1)*R);

%  Distance de l'origine translatee au point d'obs   
 R_l = sqrt(sum(P_l.^2)); 
rr_l = P_l./(ones(3,1)*R_l);

% Matrice de passage pour passer du repere faisceau
% au repere absolu
R_Bfsc2B = [Champ.ex cross(Champ.ez,Champ.ex) Champ.ez];

% raccourcis de notation de la matrice de rotation
% Matrice de passage pour passer du repere 
R_B2Bfsc = transpose(R_Bfsc2B);
     r11 = R_B2Bfsc(1,1);
     r12 = R_B2Bfsc(1,2);
     r13 = R_B2Bfsc(1,3);
     r21 = R_B2Bfsc(2,1);
     r22 = R_B2Bfsc(2,2);
     r23 = R_B2Bfsc(2,3);

%  %  Origine translatée du FG
%     P_l_O_fg = Champ.Centre - R_inter_Ofg;
%  sum_R_l_Ofg = sqrt(sum(P_l_O_fg.^2));
%  sum_P_fsc_O_fg = sqrt(sum(P_l_O_fg.^2));

%  distance parcourue par le fsc 
%  dans le repere faisceau
P_fg_Op = R_B2Bfsc*(R_inter_Ofg-Champ.Centre);
r_fg_Op = d;


%  distance parcourue par le faisceau pour 
%  couper le plan z=0
r_fg_Op = - Champ.Centre(3)./(dot(Champ.ez, [0;0;1]));


%  elements de la matrice K, correspondant à l'inverse de la
%  matrice Q0
detQ0 = Champ.Courbure(1).*Champ.Courbure(2) - Champ.Courbure(3).^2;
invQ0 = [Champ.Courbure(2); Champ.Courbure(1); -Champ.Courbure(3)]./(ones(3,1)*detQ0);
K11 = invQ0(1);
K22 = invQ0(2);
K12 = invQ0(3);

%  Raccourcis de notation utile pour la forme quadratique
d1 = r11*(R_inter_Ofg(1) - Champ.Centre(1)) ...
    +r12*(R_inter_Ofg(2) - Champ.Centre(2)) ...
    -r13*Champ.Centre(3);

d2 = r21*(R_inter_Ofg(1) - Champ.Centre(1)) ...
    +r22*(R_inter_Ofg(2) - Champ.Centre(2)) ...
    -r23*Champ.Centre(3);

%  % ============================================
%  % elements de la forme quadratique (QF)
%  % dans sa notation matricielle
%  % ============================================
%  rot = [r11,r12;r21,r22];
%    d = [d1;d2];
%   kk = [K11,K12;K12,K22];
%  
%  %  A
%  PP = [1 - P_l(1,:).^2./R_l.^2; ...
%        1 - P_l(2,:).^2./R_l.^2; ...
%        - P_l(1,:).*P_l(2,:)./R_l.^2]./(ones(3,1)*R_l);
%  
%  PP_l_O_fg = [1 - P_l_O_fg(1).^2./sum_R_l_Ofg.^2; ...
%         1 - P_l_O_fg(2).^2./sum_R_l_Ofg.^2; ...
%        - P_l_O_fg(1).*P_l_O_fg(2)./sum_R_l_Ofg.^2]./(ones(3,1)*sum_R_l_Ofg);
%  
%  TMP = transpose(rot)*kk*rot./sum_R_l_Ofg.^2;
%    
%  A = PP ...
%      + PP_l_O_fg*ones(1,NbP) ...
%      - [TMP(1,1);TMP(2,2);TMP(1,2)]*ones(1,NbP);
%  
%  %  B
%  TMP = transpose(rot)*kk*d./sum_R_l_Ofg.^2;
%  
%  B = [P_l(1,:); P_l(2,:)]./(ones(2,1)*R_l) ...
%      + [P_l_O_fg(1); P_l_O_fg(2)]./sum_R_l_Ofg*ones(1,NbP) ...
%      + [TMP(1); TMP(2)]*ones(1,NbP);
%  
%  QF.A11 = A(1,:);
%  QF.A22 = A(2,:);
%  QF.A12 = A(3,:);
%   QF.B1 = B(1,:);
%   QF.B2 = B(2,:);
%    QF.C = R_l + sum_R_l_Ofg - transpose(d)*kk*d./(2*sum_R_l_Ofg.^2);

% ============================================
% elements de la forme quadratique (QF)
% dans sa notation matricielle
% VERSION THESE
% ============================================

%  qq. raccourcis de notations
 M = [r11,r12;r12,r22];
 m11 = r11;
 m12 = r12;
 m21 = r21;
 m22 = r22;
 m13 = r13; 
 m23 = r23;

 d1 = m11*(R_inter_Ofg(1) - Champ.Centre(1)) ...
     +m12*(R_inter_Ofg(2) - Champ.Centre(2)) ...
     -m13*Champ.Centre(3);
 d2 = m21*(R_inter_Ofg(1) - Champ.Centre(1)) ...
     +m22*(R_inter_Ofg(2) - Champ.Centre(2)) ...
     -m23*Champ.Centre(3);

  
 tMiQ0 = transpose(M)*[invQ0(1),invQ0(3);invQ0(3),invQ0(2)];
 tMiQ0M = tMiQ0*M;

%  coefficients de la forme quadratique canonique finale
QF.A11 = 1./R_l - (P_l(1,:).^2)./(R_l.^3) ...
      + 1./r_fg_Op - (R_inter_Ofg(1) - Champ.Centre(1)).^2./(r_fg_Op.^3) ...
      - tMiQ0M(1,1)./(r_fg_Op.^2) ;
QF.A22 = 1./R_l - (P_l(2,:).^2)./(R_l.^3) ...
      + 1./r_fg_Op - (R_inter_Ofg(2) - Champ.Centre(2)).^2./(r_fg_Op.^3) ...
      - tMiQ0M(2,2)./(r_fg_Op.^2);
QF.A12 = - (P_l(1,:).*P_l(2,:))./(R_l.^3) ...
      - (R_inter_Ofg(1) - Champ.Centre(1)).*(R_inter_Ofg(2) - Champ.Centre(2))./(r_fg_Op.^3) ...
      - tMiQ0M(1,2)./(r_fg_Op.^2);

QF.B1 = P_l(1,:)./R_l ...
      - (R_inter_Ofg(1) - Champ.Centre(1))./r_fg_Op ...
      - tMiQ0(1,1)./(r_fg_Op.^2).*d1 ...
      - tMiQ0(1,2)./(r_fg_Op.^2).*d2;
QF.B2 = P_l(2,:)./R_l ...
      - (R_inter_Ofg(2) - Champ.Centre(2))./r_fg_Op ...
      - tMiQ0(2,1)./(r_fg_Op.^2).*d1 ...
      - tMiQ0(2,2)./(r_fg_Op.^2).*d2;

QF.C = r_fg_Op ...
      + R_l ...
      - 1./(2*r_fg_Op.^2).*(invQ0(1).*d1.^2 + invQ0(2).*d2.^2 + 2*invQ0(3).*d1.*d2);




% ======================================================
% Calcul du point col (point critique du premier ordre)
% et des fonctions evaluees en ce point
% ======================================================
detA = (QF.A11.*QF.A22 - QF.A12.^2);
  iA = [QF.A22;QF.A11;-QF.A12]./(ones(3,1)*detA);
%  point col dans le repere absolu translate
xQ_l_s = iA(1,:).*QF.B1 + iA(3,:).*QF.B2;
yQ_l_s = iA(3,:).*QF.B1 + iA(2,:).*QF.B2;
rQ_l_s = [xQ_l_s; yQ_l_s; zeros(1,NbP)];
%  point col dans le repere absolu 
xQ_s = xQ_l_s + R_inter_Ofg(1);
yQ_s = yQ_l_s + R_inter_Ofg(2);
rQ_s = [xQ_s; yQ_s; Plaque.Centre(3)*ones(1,NbP)];

% point col dans le repere fsc
    rQ_fsc_s = R_B2Bfsc*(rQ_s  - Champ.Centre*ones(1,NbP));
sum_rQ_fsc_s = sqrt(sum(rQ_fsc_s.^2));


% vecteur distance au point col
    r_s = P - rQ_s;
sum_r_s = sqrt(sum(r_s.^2));
   rr_s = r_s./(ones(3,1)*sum_r_s);

%  fonction de phase evaluee au point critique du premier ordre
g_s_l = subfun_g_l(xQ_l_s, yQ_l_s, QF);

%  fonction de polarisation dans le repere faisceau incident
%  evaluee au point critique du premier ordre
ah_fsc_s = subfun_ah(rQ_fsc_s, Champ);

%  fonction de polarisation dans le repere absolu
%  evaluee au point critique du premier ordre
ah_s = R_Bfsc2B*ah_fsc_s;

ah_Op = R_Bfsc2B*subfun_ah(P_fg_Op, Champ);

% fonction d'amplitude au point col
% Le signe (-1) provient du signe de la racine carree complexe : 
% (z)^(1/2) = -j * sqrt(-z)
ampl_s =  -sqrt(epsr./(detQ0.*detA)) ...
        ./ sum_r_s ...
        .* rQ_fsc_s(3,:)./sum_rQ_fsc_s.^2;
%          .* P_fg_Op(3)./r_fg_Op.^2;

% =====================================================
% point critique de premier ordre (point col)
% =====================================================
Epc1 = cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), ah_s))) .* ...
      (ones(3,1)*(ampl_s .* exp(-k*g_s_l)));




% ======================================================
% Calcul des points critiques du second ordre
% et des fonctions evaluees en ces points
% ======================================================
  % points critiques du second ordre
  % dans le repere modifie correspondant a l intersection
  % du FG avec la plaque
  x_Y1_l = (QF.B1-QF.A12.*Y1_l)./QF.A11;
  x_Y2_l = (QF.B1-QF.A12.*Y2_l)./QF.A11;
  y_X1_l = (QF.B2-QF.A12.*X1_l)./QF.A22;
  y_X2_l = (QF.B2-QF.A12.*X2_l)./QF.A22;
%  
%    % points critiques du second ordre
%    % dans le repere absolu lie a la plaque
%    x_Y1 = x_Y1_l + R_inter_Ofg(1);
%    x_Y2 = x_Y2_l + R_inter_Ofg(1);
%    y_X1 = y_X1_l + R_inter_Ofg(2);
%    y_X2 = y_X2_l + R_inter_Ofg(2);
%  
%    % vecteur distance aux points critiques du second ordre
%      rQ_X1 = [X1*ones(1,NbP); y_X1;Plaque.Centre(3)*ones(1,NbP)]; 
  rQ_X1_l = [X1_l*ones(1,NbP); y_X1_l;Plaque.Centre(3)*ones(1,NbP)]; 
%       r_X1 = P_l - rQ_X1_l;
%   sum_r_X1 = sqrt(sum((r_X1.^2)));
%      rr_X1 = r_X1./(ones(3,1)*sum_r_X1);
%      rQ_fsc_X1 = R_B2Bfsc*(rQ_X1 - Champ.Centre*ones(1,NbP));
%    sum_rQ_fsc_X1 = sqrt(sum(rQ_fsc_X1.^2));
%      
%    
%      rQ_X2 = [X2*ones(1,NbP); y_X2;Plaque.Centre(3)*ones(1,NbP)]; 
  rQ_X2_l = [X2_l*ones(1,NbP); y_X2_l;Plaque.Centre(3)*ones(1,NbP)]; 
%       r_X2 = P_l - rQ_X2_l;
%   sum_r_X2 = sqrt(sum((r_X2.^2)));
%      rr_X2 = r_X2./(ones(3,1)*sum_r_X2);
%        rQ_fsc_X2 = R_B2Bfsc*(rQ_X2 - Champ.Centre*ones(1,NbP));
%    sum_rQ_fsc_X2 = sqrt(sum(rQ_fsc_X2.^2));
%  
%  
%      rQ_Y1 = [x_Y1; Y1*ones(1,NbP); Plaque.Centre(3)*ones(1,NbP)]; 
  rQ_Y1_l = [x_Y1_l; Y1_l*ones(1,NbP); Plaque.Centre(3)*ones(1,NbP)]; 
%       r_Y1 = P_l - rQ_Y1_l;
%   sum_r_Y1 = sqrt(sum((r_Y1.^2)));
%      rr_Y1 = r_Y1./(ones(3,1)*sum_r_Y1);
%          rQ_fsc_Y1 = R_B2Bfsc*(rQ_Y1 - Champ.Centre*ones(1,NbP));
%    sum_rQ_fsc_Y1 = sqrt(sum(rQ_fsc_Y1.^2));
%  
%      rQ_Y2 = [x_Y2; Y2*ones(1,NbP); Plaque.Centre(3)*ones(1,NbP)]; 
  rQ_Y2_l = [x_Y2_l; Y2_l*ones(1,NbP); Plaque.Centre(3)*ones(1,NbP)]; 
%       r_Y2 = P_l - rQ_Y2_l;
%   sum_r_Y2 = sqrt(sum((r_Y2).^2));
%      rr_Y2 = r_Y2./(ones(3,1)*sum_r_Y2);
%        rQ_fsc_Y2 = R_B2Bfsc*(rQ_Y2 - Champ.Centre*ones(1,NbP));
%    sum_rQ_fsc_Y2 = sqrt(sum(rQ_fsc_Y2.^2));
%  
%  
%    % point critiques du second ordre dans le repere du faisceau incident
%    rQ_fsc_X1 = R_B2Bfsc*(rQ_X1  - Champ.Centre*ones(1,NbP));
%    rQ_fsc_X2 = R_B2Bfsc*(rQ_X2  - Champ.Centre*ones(1,NbP));
%    rQ_fsc_Y1 = R_B2Bfsc*(rQ_Y1  - Champ.Centre*ones(1,NbP));
%    rQ_fsc_Y2 = R_B2Bfsc*(rQ_Y2  - Champ.Centre*ones(1,NbP));
%  
%    %  Fonctions de polarisations aux points critiques du second ordre
%    ah_X1 = R_Bfsc2B*subfun_ah(rQ_fsc_X1, Champ);
%    ah_X2 = R_Bfsc2B*subfun_ah(rQ_fsc_X2, Champ);
%    ah_Y1 = R_Bfsc2B*subfun_ah(rQ_fsc_Y1, Champ);
%    ah_Y2 = R_Bfsc2B*subfun_ah(rQ_fsc_Y2, Champ);
%  
%  
  %  Fonctions de phase aux points critiques du second ordre
  g_X1 = subfun_g_l(rQ_X1_l(1,:), rQ_X1_l(2,:), QF);
  g_X2 = subfun_g_l(rQ_X2_l(1,:), rQ_X2_l(2,:), QF);
  g_Y1 = subfun_g_l(rQ_Y1_l(1,:), rQ_Y1_l(2,:), QF);
  g_Y2 = subfun_g_l(rQ_Y2_l(1,:), rQ_Y2_l(2,:), QF);
%  
%  
%  
%    % derivees de la fonction de phase aux points critiques
%    dgdx_X1_l = j*(QF.A11-QF.A12.^2./QF.A22).*X1_l-j*QF.B1+j*QF.B2.*QF.A12./QF.A22;
%    dgdx_X2_l = j*(QF.A11-QF.A12.^2./QF.A22).*X2_l-j*QF.B1+j*QF.B2.*QF.A12./QF.A22;
%    
%    dgdy_Y1_l = j*(-QF.A12.^2./QF.A11+QF.A22).*Y1_l+j*QF.B1.*QF.A12./QF.A11-j*QF.B2;
%    dgdy_Y2_l = j*(-QF.A12.^2./QF.A11+QF.A22).*Y2_l+j*QF.B1.*QF.A12./QF.A11-j*QF.B2;

  %  terme de changement de variable dans le passage 
  %  aux coordonnees "col" : dx/ds
  d2gdx2 =  j*QF.A11-j*QF.A12.^2./QF.A22;
  d2gdy2 = -j*QF.A12.^2./QF.A11+j*QF.A22;
  hh_s_X  = sqrt(2./d2gdx2);
  hh_s_Y  = sqrt(2./d2gdy2);


%  %    Le signe de la racine carrée est défini tel que l'on doit avoir
%  %    h_x -> h_s lorsque s -> 0, cad lq 
%    s_X1 = sign(real(hh_s_X./(X1_l-xQ_l_s))).*sqrt(g_X1 - g_s_l);
%    s_X2 = sign(real(hh_s_X./(X2_l-xQ_l_s))).*sqrt(g_X2 - g_s_l);
%    s_Y1 = sign(real(hh_s_Y./(Y1_l-yQ_l_s))).*sqrt(g_Y1 - g_s_l);
%    s_Y2 = sign(real(hh_s_Y./(Y2_l-yQ_l_s))).*sqrt(g_Y2 - g_s_l);

   %  Methode Felsen
  s_X1 = exp(j*angle((X1_l-xQ_l_s).*sqrt(d2gdx2/2))).*sqrt(abs(g_X1 - g_s_l));
  s_X2 = exp(j*angle((X2_l-xQ_l_s).*sqrt(d2gdx2/2))).*sqrt(abs(g_X2 - g_s_l));
  s_Y1 = exp(j*angle((Y1_l-yQ_l_s).*sqrt(d2gdy2/2))).*sqrt(abs(g_Y1 - g_s_l));
  s_Y2 = exp(j*angle((Y2_l-yQ_l_s).*sqrt(d2gdy2/2))).*sqrt(abs(g_Y2 - g_s_l));

% =====================================================
% point critique du second ordre
% =====================================================
%  E_X1 = ...
%        Epc1.*sqrt(1/pi).*(ones(3,1)*subfun_Transfun(sqrt(k).*s_X1)) ...
%        - ...
%       cross(rr_X1, cross(rr_X1, cross(Plaque.N*ones(1,NbP), ah_X1))) .* ...
%        (ones(3,1)*( ...
%          1 .* sqrt(j*epsr./(2*pi*k*detQ0*QF.A22))  ...
%          .* P_fg_Op(3)./r_fg_Op.^2 ...
%          ./ sum_r_X1 ...
%          ./ dgdx_X1_l ...
%          .* exp(-k*g_X1))) ...
%        + ...
%        cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), ah_X1))) .* ...
%        (ones(3,1)*( ...
%          1/2./s_X1.*sqrt(epsr./(pi*k*detQ0*detA)) ...
%          .* P_fg_Op(3)./r_fg_Op.^2 ...
%          ./ sum_r_s ...
%          .* exp(-k*g_X1) ));
%  
%  
%  
%  E_X2 = ...
%        Epc1.*sqrt(1/pi).*(ones(3,1)*(sqrt(pi)-subfun_Transfun(sqrt(k).*s_X2))) ...
%        + ...
%       cross(rr_X2, cross(rr_X2, cross(Plaque.N*ones(1,NbP), ah_X2))) .* ...
%        (ones(3,1)*( ...
%          1 .* sqrt(j*epsr./(2*pi*k*detQ0*QF.A22))  ...
%          .* P_fg_Op(3)./r_fg_Op.^2 ...
%          ./ sum_r_X2 ...
%          ./ dgdx_X2_l ...
%          .* exp(-k*g_X2))) ...
%        - ...
%        cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), ah_X2))) .* ...
%        (ones(3,1)*( ...
%          1/2./s_X2.*sqrt(epsr./(pi*k*detQ0*detA)) ...
%          .* P_fg_Op(3)./r_fg_Op.^2 ...
%          ./ sum_r_s ...
%          .* exp(-k*g_X2) ));
%  
%  E_Y1 = ...
%        Epc1.*sqrt(1/pi).*(ones(3,1)*subfun_Transfun(sqrt(k).*s_Y1)) ...
%        - ...
%       cross(rr_Y1, cross(rr_Y1, cross(Plaque.N*ones(1,NbP), ah_Y1))) .* ...
%        (ones(3,1)*( ...
%          1 .* sqrt(j*epsr./(2*pi*k*detQ0*QF.A11))  ...
%          .* P_fg_Op(3)./r_fg_Op.^2 ...
%          ./ sum_r_Y1 ...
%          ./ dgdy_Y1_l ...
%          .* exp(-k*g_Y1))) ...
%        + ...
%        cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), ah_Y1))) .* ...
%        (ones(3,1)*( ...
%          1/2./s_Y1.*sqrt(epsr./(pi*k*detQ0*detA)) ...
%          .* P_fg_Op(3)./r_fg_Op.^2 ...
%          ./ sum_r_s ...
%          .* exp(-k*g_Y1) ));
%  
%  E_Y2 = ...
%        Epc1.*sqrt(1/pi).*(ones(3,1)*(sqrt(pi)-subfun_Transfun(sqrt(k).*s_Y2))) ...
%        + ...
%       cross(rr_Y2, cross(rr_Y2, cross(Plaque.N*ones(1,NbP), ah_Y2))) .* ...
%        (ones(3,1)*( ...
%          1 .* sqrt(j*epsr./(2*pi*k*detQ0*QF.A11))  ...
%          .* P_fg_Op(3)./r_fg_Op.^2 ...
%          ./ sum_r_Y2 ...
%          ./ dgdy_Y2_l ...
%          .* exp(-k*g_Y2))) ...
%        - ...
%        cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), ah_Y2))) .* ...
%        (ones(3,1)*( ...
%          1/2./s_Y2.*sqrt(epsr./(pi*k*detQ0*detA)) ...
%          .* P_fg_Op(3)./r_fg_Op.^2 ...
%          ./ sum_r_s ...
%          .* exp(-k*g_Y2) ));

%  VERSION "COURTE"
E_X1 = Epc1.*sqrt(1/pi).*(ones(3,1)*subfun_Transfun(sqrt(k).*s_X1));
E_X2 = Epc1.*sqrt(1/pi).*(ones(3,1)*(sqrt(pi)-subfun_Transfun(sqrt(k).*s_X2)));
E_Y1 = Epc1.*sqrt(1/pi).*(ones(3,1)*subfun_Transfun(sqrt(k).*s_Y1));
E_Y2 = Epc1.*sqrt(1/pi).*(ones(3,1)*(sqrt(pi)-subfun_Transfun(sqrt(k).*s_Y2)));

%  on filtre les NaN
E_X1(find(isnan(E_X1)))=0;
E_X2(find(isnan(E_X2)))=0;
E_Y1(find(isnan(E_Y1)))=0;
E_Y2(find(isnan(E_Y2)))=0;

Epc2 =  - E_X1 - E_X2 - E_Y1 - E_Y2;

%  figure(10)
%    plot(P(1,:), 20*log10(abs([E_X1(1,:);E_X2(1,:);E_Y1(1,:);E_Y2(1,:)])));
%    legend('X_1','X_2','Y_1','Y_2')
%    axis([-29 29 -150 -40])

%  ==============================================================
%  Coins
%  ==============================================================
%  
%  % vecteur distance aux coins
%      rQ_X1Y1 = [X1; Y1;0]; 
%    rQ_X1Y1_l = [X1_l; Y1_l;0]; 
%       r_X1Y1 = P_l - rQ_X1Y1_l*ones(1,NbP);
%   sum_r_X1Y1 = sqrt(sum((r_X1Y1.^2)));
%      rr_X1Y1 = r_X1Y1./(ones(3,1)*sum_r_X1Y1);
%  %  dans le repere du faisceau incident
%      rQ_fsc_X1Y1 = R_B2Bfsc*(rQ_X1Y1  - Champ.Centre)*ones(1,NbP);
%  sum_rQ_fsc_X1Y1 = sqrt(sum(rQ_fsc_X1Y1.^2));
%  
%      rQ_X2Y1 = [X2; Y1;0]; 
%    rQ_X2Y1_l = [X2_l; Y1_l;0]; 
%       r_X2Y1 = P_l - rQ_X2Y1_l*ones(1,NbP);
%   sum_r_X2Y1 = sqrt(sum((r_X2Y1.^2)));
%      rr_X2Y1 = r_X2Y1./(ones(3,1)*sum_r_X2Y1);
%  %  dans le repere du faisceau incident
%      rQ_fsc_X2Y1 = R_B2Bfsc*(rQ_X2Y1  - Champ.Centre)*ones(1,NbP);
%  sum_rQ_fsc_X2Y1 = sqrt(sum(rQ_fsc_X2Y1.^2));
%  
%      rQ_X1Y2 = [X1; Y2;0]; 
%    rQ_X1Y2_l = [X1_l; Y2_l;0]; 
%       r_X1Y2 = P_l - rQ_X1Y2_l*ones(1,NbP);
%   sum_r_X1Y2 = sqrt(sum((r_X1Y2.^2)));
%      rr_X1Y2 = r_X1Y2./(ones(3,1)*sum_r_X1Y2);
%  %  dans le repere du faisceau incident
%      rQ_fsc_X1Y2 = R_B2Bfsc*(rQ_X1Y2  - Champ.Centre)*ones(1,NbP);
%  sum_rQ_fsc_X1Y2 = sqrt(sum(rQ_fsc_X1Y2.^2));
%  
%      rQ_X2Y2 = [X2; Y2;0]; 
%    rQ_X2Y2_l = [X2_l; Y2_l;0]; 
%       r_X2Y2 = P_l - rQ_X2Y2_l*ones(1,NbP);
%   sum_r_X2Y2 = sqrt(sum((r_X2Y2.^2)));
%      rr_X2Y2 = r_X2Y2./(ones(3,1)*sum_r_X2Y2);
%  %  dans le repere du faisceau incident
%      rQ_fsc_X2Y2 = R_B2Bfsc*(rQ_X2Y2  - Champ.Centre)*ones(1,NbP);
%  sum_rQ_fsc_X2Y2 = sqrt(sum(rQ_fsc_X2Y2.^2));
%  
%  
%  %  fonction d'amplitude aux coins
%  ampl_X1Y1 =  j*k/2/pi*sqrt(1./detQ0)  ...
%          .* P_fg_Op(3)./r_fg_Op.^2 ...
%          ./ sum_r_X1Y1;
%  ampl_X2Y1 =  j*k/2/pi*sqrt(1./detQ0)  ...
%          .* P_fg_Op(3)./r_fg_Op.^2 ...
%          ./ sum_r_X2Y1;
%  ampl_X1Y2 =  j*k/2/pi*sqrt(1./detQ0)  ...
%          .* P_fg_Op(3)./r_fg_Op.^2 ...
%          ./ sum_r_X1Y2;
%  ampl_X2Y2 =  j*k/2/pi*sqrt(1./detQ0)  ...
%          .* P_fg_Op(3)./r_fg_Op.^2 ...
%          ./ sum_r_X2Y2;
%  
%  %  ampl_X1Y1 =  j*k/2/pi*sqrt(1./detQ0)  ...
%  %          .* rQ_fsc_X1Y1(3)./sum_rQ_fsc_X1Y1.^2 ...
%  %          ./ sum_r_X1Y1;
%  %  ampl_X2Y1 =  j*k/2/pi*sqrt(1./detQ0)  ...
%  %          .* rQ_fsc_X2Y1(3)./sum_rQ_fsc_X2Y1.^2 ...
%  %          ./ sum_r_X2Y1;
%  %  ampl_X1Y2 =  j*k/2/pi*sqrt(1./detQ0)  ...
%  %          .* rQ_fsc_X1Y2(3)./sum_rQ_fsc_X1Y2.^2 ...
%  %          ./ sum_r_X1Y2;
%  %  ampl_X2Y2 =  j*k/2/pi*sqrt(1./detQ0)  ...
%  %          .* rQ_fsc_X2Y2(3)./sum_rQ_fsc_X2Y2.^2 ...
%  %          ./ sum_r_X2Y2;
%  
%  
%  %  Fonctions de phase aux coins
%  g_X1Y1 = subfun_g_l(rQ_X1Y1_l(1,:), rQ_X1Y1_l(2,:), QF);
%  g_X2Y1 = subfun_g_l(rQ_X2Y1_l(1,:), rQ_X2Y1_l(2,:), QF);
%  g_X1Y2 = subfun_g_l(rQ_X1Y2_l(1,:), rQ_X1Y2_l(2,:), QF);
%  g_X2Y2 = subfun_g_l(rQ_X2Y2_l(1,:), rQ_X2Y2_l(2,:), QF);
%  
%  %  Fonctions de polarisations aux coins
%  h_X1Y1 = R_Bfsc2B*subfun_ah(rQ_fsc_X1Y1, Champ);
%  h_X2Y1 = R_Bfsc2B*subfun_ah(rQ_fsc_X2Y1, Champ);
%  h_X1Y2 = R_Bfsc2B*subfun_ah(rQ_fsc_X1Y2, Champ);
%  h_X2Y2 = R_Bfsc2B*subfun_ah(rQ_fsc_X2Y2, Champ);
%  
%  %  %  derivee simples aux coins
%  %  dgdx_X1Y1_l = j*QF.A11.*X1_l + j*QF.A12.*Y1_l - j*QF.B1;
%  %  dgdy_X1Y1_l = j*QF.A12.*X1_l + j*QF.A22.*Y1_l - j*QF.B2;
%  %  
%  %  dgdx_X2Y1_l = j*QF.A11.*X2_l + j*QF.A12.*Y1_l - j*QF.B1;
%  %  dgdy_X2Y1_l = j*QF.A12.*X2_l + j*QF.A22.*Y1_l - j*QF.B2;
%  %  
%  %  dgdx_X1Y2_l = j*QF.A11.*X1_l + j*QF.A12.*Y2_l - j*QF.B1;
%  %  dgdy_X1Y2_l = j*QF.A12.*X1_l + j*QF.A22.*Y2_l - j*QF.B2;
%  %  
%  %  dgdx_X2Y2_l = j*QF.A11.*X2_l + j*QF.A12.*Y2_l - j*QF.B1;
%  %  dgdy_X2Y2_l = j*QF.A12.*X2_l + j*QF.A22.*Y2_l - j*QF.B2;
%  
%  Methode Felsen
%  hypothese : terme croise nul
g_X1 = j/2*QF.A11.*X1_l.^2 - j*QF.B1.*X1_l;
g_X2 = j/2*QF.A11.*X2_l.^2 - j*QF.B1.*X2_l;
g_Y1 = j/2*QF.A22.*Y1_l.^2 - j*QF.B2.*Y1_l;
g_Y2 = j/2*QF.A22.*Y2_l.^2 - j*QF.B2.*Y2_l;

g_sx = j/2*QF.A11.*xQ_l_s.^2 - j*QF.B1.*xQ_l_s;
g_sy = j/2*QF.A22.*yQ_l_s.^2 - j*QF.B2.*yQ_l_s;

s_X1 = exp(j*angle((X1_l-xQ_l_s).*sqrt(j*QF.A11/2))).*sqrt(abs(g_X1 - g_sx));
s_Y1 = exp(j*angle((Y1_l-yQ_l_s).*sqrt(j*QF.A22/2))).*sqrt(abs(g_Y1 - g_sy));
s_X2 = exp(j*angle((X2_l-xQ_l_s).*sqrt(j*QF.A11/2))).*sqrt(abs(g_X2 - g_sx));
s_Y2 = exp(j*angle((Y2_l-yQ_l_s).*sqrt(j*QF.A22/2))).*sqrt(abs(g_Y2 - g_sy));

nu_X1Y1 = sqrt(k).*s_X1;
mu_X1Y1 = sqrt(k).*s_Y1;

nu_X2Y1 = sqrt(k).*s_X2;
mu_X2Y1 = sqrt(k).*s_Y1;

nu_X1Y2 = sqrt(k).*s_X1;
mu_X1Y2 = sqrt(k).*s_Y2;

nu_X2Y2 = sqrt(k).*s_X2;
mu_X2Y2 = sqrt(k).*s_Y2;
%  
%  
%  
%  
%  A11A22 = -QF.A11.*QF.A22;
%  
%  
%  E_X1Y1 = cross(rr_X1Y1, cross(rr_X1Y1, cross(Plaque.N*ones(1,NbP), h_X1Y1))) .* ...
%        (ones(3,1)*( 2/k*ampl_X1Y1 .* exp(-k*g_X1Y1) ...
%                    .* sqrt(1./(A11A22)) ...
%                    .* subfun_Transfun(nu_X1Y1) ...
%                    .* subfun_Transfun(mu_X1Y1) ...
%                    .* exp(nu_X1Y1.^2) ...
%                    .* exp(mu_X1Y1.^2) ));
%  
%  
%  E_X2Y1 = cross(rr_X2Y1, cross(rr_X2Y1, cross(Plaque.N*ones(1,NbP), h_X2Y1))) .* ...
%        (ones(3,1)*( 2/k*ampl_X2Y1 .* exp(-k*g_X2Y1) ...
%                    .* sqrt(1./(A11A22)) ...
%                    .* (sqrt(pi) - subfun_Transfun(nu_X2Y1)) ...
%                    .* subfun_Transfun(mu_X2Y1) ...
%                    .* exp(nu_X2Y1.^2) ...
%                    .* exp(mu_X2Y1.^2) ));
%  
%  E_X1Y2 = cross(rr_X1Y2, cross(rr_X1Y2, cross(Plaque.N*ones(1,NbP), h_X1Y2))) .* ...
%        (ones(3,1)*( 2/k*ampl_X1Y2 .* exp(-k*g_X1Y2) ...
%                    .* sqrt(1./(A11A22)) ...
%                    .* subfun_Transfun(nu_X1Y2) ...
%                    .* (sqrt(pi) - subfun_Transfun(mu_X1Y2)) ...
%                    .* exp(nu_X1Y2.^2) ...
%                    .* exp(mu_X1Y2.^2) ));
%  
%  E_X2Y2 = cross(rr_X2Y2, cross(rr_X2Y2, cross(Plaque.N*ones(1,NbP), h_X2Y2))) .* ...
%        (ones(3,1)*( 2/k*ampl_X2Y1 .* exp(-k*g_X2Y2) ...
%                    .* sqrt(1./(A11A22)) ...
%                    .* (sqrt(pi) - subfun_Transfun(nu_X2Y2)) ...
%                    .* (sqrt(pi) - subfun_Transfun(mu_X2Y2)) ...
%                    .* exp(nu_X2Y2.^2) ...
%                    .* exp(mu_X2Y2.^2) ));

%  VERSION "COURTE"
E_X1Y1 = Epc1/pi.*(ones(3,1)*(subfun_Transfun(nu_X1Y1).*subfun_Transfun(mu_X1Y1)));
E_X2Y1 = Epc1/pi.*(ones(3,1)*((sqrt(pi)-subfun_Transfun(nu_X2Y1)).*subfun_Transfun(mu_X2Y1)));
E_X1Y2 = Epc1/pi.*(ones(3,1)*((sqrt(pi)-subfun_Transfun(mu_X1Y2)).*subfun_Transfun(nu_X1Y2)));
E_X2Y2 = Epc1/pi.*(ones(3,1)*((sqrt(pi)-subfun_Transfun(nu_X2Y2)).*(sqrt(pi)-subfun_Transfun(mu_X2Y2))));


%  on filtre les NaN
E_X1Y1(find(isnan(E_X1Y1)))=0;
E_X1Y2(find(isnan(E_X1Y2)))=0;
E_X2Y1(find(isnan(E_X2Y1)))=0;
E_X2Y2(find(isnan(E_X2Y2)))=0;


Epc3 = + E_X1Y1 + E_X2Y2 + E_X2Y1 + E_X1Y2 ;


%  figure(6)
%    plot(P(1,:), 20*log10(abs([E_X1Y1(1,:);E_X2Y2(1,:);E_X2Y1(1,:);E_X1Y2(1,:)])))
%    legend('X_1 Y_1','X_2 Y_2','X_2 Y_1','X_1 Y_2')
%    axis([-22.4 22.4 -150 -20])
% ==========================================================
% champ total
% ==========================================================
E = Epc1+Epc2+Epc3;

H = zeros(size(E));





%  Fonction de phase
function g = subfun_g_l(xQ_l,yQ_l,QF)
  g = j/2*(xQ_l.^2.*QF.A11 + 2*xQ_l.*yQ_l.*QF.A12 + yQ_l.^2.*QF.A22) ...
      - j*QF.B1.*xQ_l - j*QF.B2.*yQ_l ...
      + j*QF.C;


function a_h = subfun_ah(P, Champ)
%  fonction de polarisation du FG lointain magnetique 
%  dans son repere faisceau
%  
%  ENTREES :
%  P : points d'evaluation (exprimes dans le repere faisceau)
%  Champ : structure Champ
%  
%  RETOURNE :
%  a_h : vecteur polarisation (3xNbP)
%  
%  ATTENTION : ne plas oulier de changer le repere du vecteur h
%  dans le repere absolu par la suite.
  NbP = size(P,2);

  R = sqrt(sum(P.^2));

  a_h_x = -(P(1,:).*P(2,:).*Champ.Coefficients(1) + (P(2,:).^2+P(3,:).^2).*Champ.Coefficients(2))./(P(3,:).*R);
  a_h_y = +(P(1,:).*P(2,:).*Champ.Coefficients(2) + (P(1,:).^2+P(3,:).^2).*Champ.Coefficients(1))./(P(3,:).*R);
  a_h_z =  +P(1,:)./R.*Champ.Coefficients(2) - P(2,:)./R.*Champ.Coefficients(1);

  a_h = [a_h_x;a_h_y;a_h_z];


function F = subfun_Transfun(z)
%  Fonction de transition pour dev. asyptotique uniforme
cerr_fun = 1-erfz(z);
%  cerr_fun = exp(-z.^2).*W(i*z);
%  cerr_fun = exp(-z.^2).*faddeeva(i*z);

F = sqrt(pi)/2.*cerr_fun;

