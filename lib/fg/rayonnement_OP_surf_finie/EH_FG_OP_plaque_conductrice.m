function [E,H, Epc1, Epc2, Epc3] = EH_FG_OP_Plaque_conductrice(Champ, Plaque, P, f, epsr)
% Calcul du champ EM rayonne par une Plaque conductrice 
% dans l'hypothese de l'optique physique
% par un faisceau gaussien
% 
% Les effets de la finitude sont pris en compte
% Ainsi que les effets "de coin"
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

% ===========================================
%      CONSTANTES et PARAMETRES
% ===========================================
     c = 2.997925e8;
    Z0 = 120*pi;
lambda = c/f;
     k = 2*pi/lambda;
   NbP = length(P);

%  code en cours de debug ? 
%  (affiche la position des points critiques)
 debug = false; 

Plaque.X1 = max(Plaque.Sommets(1,:));
Plaque.X2 = min(Plaque.Sommets(1,:));
Plaque.Y1 = max(Plaque.Sommets(2,:));
Plaque.Y2 = min(Plaque.Sommets(2,:));

% matrice de courbure en z_fc = 0
Q0    = Champ.Courbure;
detQ0 = Q0(1)*Q0(2)-Q0(3)^2;  
invQ0 = [Q0(2); Q0(1); -Q0(3)] ./ (ones(3,1)*detQ0);

% Matrice de passage dans le repere associe au fsc
M_fsc = [Champ.ex cross(Champ.ez,Champ.ex) Champ.ez];
% raccourcis de notation de la matrice de rotation
R_3 = transpose(M_fsc);
r11 = R_3(1,1);
r12 = R_3(1,2);
r21 = R_3(2,1);
r22 = R_3(2,2);
r31 = R_3(3,1);
r32 = R_3(3,2);

% Coordonnees de P dans le repere fsc
P_fsc = R_3*(P-Champ.Centre*ones(1,NbP));

% position du centre de la Plaque dans le repere faisceau
PlaqueCentre_fsc = R_3*(Plaque.Centre - Champ.Centre);

% TODO : position du "centre de phase" sur la Plaque dans le rep faisceau
% Il devrait s'agir de l'intersection de la droite 
% definit par le vecteur directeur ez du repere fsc 
% avec le plan definit par les 3 premiers points de la Plaque
% ou encore mieux, du point correspondant 
% a la plus courte distance entre la Plaque et le FG 
%
% en attendant, ce point correspond au centre de la Plaque  
Centre_Phase = Plaque.Centre;
CentrePhase_fsc = PlaqueCentre_fsc;

%  distance entre le centre du repere faisceau et l'intersection de la droite de 
%  cecteur directeur ez_fsc avec le plan
d = - Champ.Centre(3)./Champ.ez(3);
%  coordonnees du point d'intersection :
CentrePhase_fsc = Champ.Centre + d*Champ.ez;

%  %  TODO : centre de phase 
%  % on trace le prolongement du vecteur directeur e_s^i jusqu'au plan z=0
%  Centre_Phase_x = Champ.ez(1)*(-Champ.Centre(3)/Champ.ez(3))+Champ.Centre(1);
%  Centre_Phase_y = Champ.ez(2)*(-Champ.centre(3)/Champ.ez(3))+Champ.Centre(2);

%  CentrePhase_fsc = transpose(M_fsc)*([Centre_Phase_x;Centre_Phase_y;0] - Champ.Centre);

        z_Q_fsc = d;

%  %  R = sqrt(sum((P-CentrePhase_fsc*ones(1,size(P,2))).^2));
R = sqrt(sum(P.^2));

% inverse de la matrice de courbure et matrice de courbure
% hypothese : zi = cte et vaut la distance pour rejoindre 
% "le centre de phase".
    inv_Qp = [invQ0(1)+z_Q_fsc; invQ0(2)+z_Q_fsc; invQ0(3)];
det_inv_Qp = inv_Qp(1,:).*inv_Qp(2,:) - inv_Qp(3,:).^2;
        Qp = [inv_Qp(2,:);inv_Qp(1,:);-inv_Qp(3,:)]...
              ./(ones(3,1)*det_inv_Qp);

% centre du faisceau exprim� dans le repere faisceau
ri_Oi = R_3*Champ.Centre;
xi_Oi = ri_Oi(1);
yi_Oi = ri_Oi(2);
zi_Oi = ri_Oi(3);



% elements de la notation matricielle de la phase
Axx = - 1./R ...
  + P(1,:).^2./R.^3 ...
  - Qp(1).*r11.^2 ...
  - Qp(2).*r21.^2 ...
  - 2*Qp(3).*r11.*r21;

Ayy = - 1./R ...
  + P(2,:).^2./R.^3 ...
  - Qp(2).*r22.^2 ...
  - Qp(1).*r12.^2 ...
  - 2*Qp(3).*r12.*r22;

Axy = + P(1,:).*P(2,:)./R.^3 ...
  - Qp(3).*(r11.*r22 + r12.*r21) ...
  - Qp(2).*r21.*r22 ...
  - Qp(1).*r11.*r12;

Bx = + r31 ...
  - P(1,:)./R ...
  - Qp(3).*(xi_Oi.*r21 + r11.*yi_Oi) ...
  - Qp(2).*r21.*yi_Oi ...
  - Qp(1).*r11.*xi_Oi;

By = + r32 ...
  - P(2,:)./R ...
  - Qp(3).*(xi_Oi.*r22 + r12.*yi_Oi) ...
  - Qp(2).*r22.*yi_Oi ...
  - Qp(1).*r12.*xi_Oi;

C = + R ...
  - zi_Oi ...
  +1/2*Qp(1).*xi_Oi.^2 ...
  +1/2*Qp(2).*yi_Oi.^2 ...
  +Qp(3).*xi_Oi.*yi_Oi;

detA = Axx.*Ayy - Axy.^2;
iA = [Ayy;Axx;-Axy]./(ones(3,1)*detA);

% ======================================================
% Calcul du point col (point critique du premier ordre)
% ======================================================
%  point col dans le repere absolu
xQ_s = iA(1,:).*Bx + iA(3,:).*By;
yQ_s = iA(3,:).*Bx + iA(2,:).*By;
rQ_s = [xQ_s; yQ_s; zeros(1,NbP)];

if debug
  %  affichage du point critique
  %  du premier ordre par rapport 
  %  au domaine d'integration
  figure(100)
    set(gca, 'FontSize', 15);
    fill([Plaque.X2 Plaque.X1 Plaque.X1 Plaque.X2], ...
         [Plaque.Y1 Plaque.Y1 Plaque.Y2 Plaque.Y2], ...
         [.7 .7 .7]);
    grid on;
    axis([Plaque.X2 - 10, Plaque.X1+10, Plaque.Y2-10, Plaque.Y1+10]);
    line(real(xQ_s), real(yQ_s), 'Color', 'r', 'LineWidth', 2);
    line(imag(xQ_s), imag(yQ_s), 'Color', 'b', 'LineWidth', 2);
    legend('Surface S', 'Re (x^{(1)}_Q, y^{(1)}_Q)', 'Im (x^{(1)}_Q, y^{(1)}_Q)', ...
           'Location', 'SouthOutside');
    xlabel('x (\lambda)');
    ylabel('y (\lambda)');
  
    tt = atan2(P(1,:), P(3,:));

    id_X2 = find(real(xQ_s) >= Plaque.X2, 1);
    id_X1 = find(real(xQ_s) >= Plaque.X1, 1);
    line(real(xQ_s(id_X2)), real(yQ_s(id_X2)), 'Color', 'k', 'Marker', '.', 'MarkerSize', 30);
    line(real(xQ_s(id_X1)), real(yQ_s(id_X1)), 'Color', 'k', 'Marker', '.', 'MarkerSize', 30);
    text(real(xQ_s(id_X2))+1, real(yQ_s(id_X2))+3, ['\theta=', num2str(180/pi*tt(id_X2)), '^\circ'], 'FontSize', 15);
    text(real(xQ_s(id_X1))+1, real(yQ_s(id_X1))+3, ['\theta=', num2str(180/pi*tt(id_X1)), '^\circ'], 'FontSize', 15);
    text(Plaque.X2+1, Plaque.Y2+2, 'S', 'FontSize', 20);
    axis equal
%      saveas(figure(100), 'diff_3D_evolution_point_col_par_rapport_S.png');
end


% point col dans le repere fsc
rQ_fsc_s = R_3*(rQ_s  - Champ.Centre*ones(1,NbP) );

% matrice de courbure evaluee au point col dans le repere faisceau
iQz_fsc_s = [invQ0(1)+rQ_fsc_s(3,:); ...
           invQ0(2)+rQ_fsc_s(3,:); ...
           invQ0(3)*ones(1,NbP)];
det_iQz_fsc_s = iQz_fsc_s(1,:).*iQz_fsc_s(2,:) - iQz_fsc_s(3,:).^2;
Qz_fsc_s = [iQz_fsc_s(2,:); ...
            iQz_fsc_s(1,:); ...
           -iQz_fsc_s(3,:)]./(ones(3,1)*det_iQz_fsc_s);
det_Qz_fsc_s = Qz_fsc_s(1,:).*Qz_fsc_s(2,:) - Qz_fsc_s(3,:).^2;

% vecteur distance au point col
    r_s = P - rQ_s;
sum_r_s = sqrt(sum(abs(r_s).^2));
   rr_s = r_s./(ones(3,1)*sum_r_s);

% fonction de phase au point col
g = -j/2*(xQ_s.^2.*Axx + 2*xQ_s.*yQ_s.*Axy + yQ_s.^2.*Ayy) ...
    + j*Bx.*xQ_s + j*By.*yQ_s ...
    + j*C;

% fonction d'amplitude au point col
% det. de la matrice Hessienne au point col : detHs = -detA
ampl =  j*sqrt(epsr) ...
        .* sqrt(det_Qz_fsc_s./detQ0./(-detA)) ...
        ./ sum_r_s;

% application de la technique du point col
saddle_point =  ampl .* exp(-k*g);


h_l_s_x = -Champ.Coefficients(2).*ones(1,NbP);
h_l_s_y =  Champ.Coefficients(1).*ones(1,NbP);
     Qx = Qz_fsc_s(1,:).*rQ_fsc_s(1,:) + Qz_fsc_s(3,:).*rQ_fsc_s(2,:);
     Qy = Qz_fsc_s(3,:).*rQ_fsc_s(1,:) + Qz_fsc_s(2,:).*rQ_fsc_s(2,:);
h_l_s_z = -Champ.Coefficients(1).*Qy + Champ.Coefficients(2).*Qx;
  h_l_s = [h_l_s_x;h_l_s_y;h_l_s_z];
  h_s   = M_fsc*h_l_s;



% =====================================================
% point critique de premier ordre (point col)
% =====================================================
% on applique le point col que lorsque le point appartient au domaine d'integration
U1= (abs(rQ_s(1,:)) >= Plaque.X2) ...
  & (abs(rQ_s(1,:)) <= Plaque.X1) ...
  & (abs(rQ_s(2,:)) >= Plaque.Y2) ...
  & (abs(rQ_s(2,:)) <= Plaque.Y1);
%  U1 = ones(1,NbP);

Epc1 = cross(rr_s, cross(rr_s, cross(Plaque.N*ones(1,NbP), h_s))) .* ...
      (ones(3,1)*(U1.*saddle_point));

%  %  error('test')
%  figure(1)
%    hold on; 
%    plot3(abs(rQ_s(1,:))/lambda, ...
%          abs(rQ_s(2,:))/lambda, ...
%          abs(rQ_s(3,:))/lambda, 'r.')  
%    hold off;
  
% =====================================================
% Contribution des points critiques du second ordre
% =====================================================
% points critiques du second ordre
x_Y1 = (Bx-Axy.*Plaque.Y1)./Axx;
x_Y2 = (Bx-Axy.*Plaque.Y2)./Axx;
y_X1 = (By-Axy.*Plaque.X1)./Ayy;
y_X2 = (By-Axy.*Plaque.X2)./Ayy;

% vecteur distance aux points critiques du second ordre
   rQ_X1 = [Plaque.X1*ones(1,NbP); y_X1;zeros(1,NbP)]; 
    r_X1 = P - rQ_X1;
sum_r_X1 = sqrt(sum((r_X1.^2)));
   rr_X1 = r_X1./(ones(3,1)*sum_r_X1);

   rQ_X2 = [Plaque.X2*ones(1,NbP); y_X2;zeros(1,NbP)]; 
    r_X2 = P - rQ_X2;
sum_r_X2 = sqrt(sum((r_X2.^2)));
   rr_X2 = r_X2./(ones(3,1)*sum_r_X2);

   rQ_Y1 = [x_Y1; Plaque.Y1*ones(1,NbP); zeros(1,NbP)]; 
    r_Y1 = P - rQ_Y1;
sum_r_Y1 = sqrt(sum((r_Y1.^2)));
   rr_Y1 = r_Y1./(ones(3,1)*sum_r_Y1);

   rQ_Y2 = [x_Y2; Plaque.Y2*ones(1,NbP); zeros(1,NbP)]; 
    r_Y2 = P - rQ_Y2;
sum_r_Y2 = sqrt(sum((r_Y2).^2));
   rr_Y2 = r_Y2./(ones(3,1)*sum_r_Y2);

% point col dans le repere incident
xQ_fsc_X1 = r11.*rQ_X1(1,:) + r12.*rQ_X1(2,:) - xi_Oi;
yQ_fsc_X1 = r21.*rQ_X1(1,:) + r22.*rQ_X1(2,:) - yi_Oi;
zQ_fsc_X1 = r31.*rQ_X1(1,:) + r32.*rQ_X1(2,:) - zi_Oi;

xQ_fsc_X2 = r11.*rQ_X2(1,:) + r12.*rQ_X2(2,:) - xi_Oi;
yQ_fsc_X2 = r21.*rQ_X2(1,:) + r22.*rQ_X2(2,:) - yi_Oi;
zQ_fsc_X2 = r31.*rQ_X2(1,:) + r32.*rQ_X2(2,:) - zi_Oi;

xQ_fsc_Y1 = r11.*rQ_Y1(1,:) + r12.*rQ_Y1(2,:) - xi_Oi;
yQ_fsc_Y1 = r21.*rQ_Y1(1,:) + r22.*rQ_Y1(2,:) - yi_Oi;
zQ_fsc_Y1 = r31.*rQ_Y1(1,:) + r32.*rQ_Y1(2,:) - zi_Oi;

xQ_fsc_Y2 = r11.*rQ_Y2(1,:) + r12.*rQ_Y2(2,:) - xi_Oi;
yQ_fsc_Y2 = r21.*rQ_Y2(1,:) + r22.*rQ_Y2(2,:) - yi_Oi;
zQ_fsc_Y2 = r31.*rQ_Y2(1,:) + r32.*rQ_Y2(2,:) - zi_Oi;

% matrice de courbure dans le repere incident
iQz_fsc_X1 = [invQ0(1)+zQ_fsc_X1; ...
              invQ0(2)+zQ_fsc_X1; ...
              invQ0(3)*ones(1,NbP)];
det_iQz_fsc_X1 = iQz_fsc_X1(1,:).*iQz_fsc_X1(2,:) - iQz_fsc_X1(3,:).^2;
Qz_fsc_X1 = [iQz_fsc_X1(2,:); ...
             iQz_fsc_X1(1,:); ...
            -iQz_fsc_X1(3,:)]./(ones(3,1)*det_iQz_fsc_X1);
det_Qz_fsc_X1 = Qz_fsc_X1(1,:).*Qz_fsc_X1(2,:) - Qz_fsc_X1(3,:).^2;

iQz_fsc_X2 = [invQ0(1)+zQ_fsc_X2; ...
              invQ0(2)+zQ_fsc_X2; ...
              invQ0(3)*ones(1,NbP)];
det_iQz_fsc_X2 = iQz_fsc_X2(1,:).*iQz_fsc_X2(2,:) - iQz_fsc_X2(3,:).^2;
Qz_fsc_X2 = [iQz_fsc_X2(2,:); ...
             iQz_fsc_X2(1,:); ...
            -iQz_fsc_X2(3,:)]./(ones(3,1)*det_iQz_fsc_X2);
det_Qz_fsc_X2 = Qz_fsc_X2(1,:).*Qz_fsc_X2(2,:) - Qz_fsc_X2(3,:).^2;

iQz_fsc_Y1 = [invQ0(1)+zQ_fsc_Y1; ...
              invQ0(2)+zQ_fsc_Y1; ...
              invQ0(3)*ones(1,NbP)];
det_iQz_fsc_Y1 = iQz_fsc_Y1(1,:).*iQz_fsc_Y1(2,:) - iQz_fsc_Y1(3,:).^2;
Qz_fsc_Y1 = [iQz_fsc_Y1(2,:); ...
          iQz_fsc_Y1(1,:); ...
         -iQz_fsc_Y1(3,:)]./(ones(3,1)*det_iQz_fsc_Y1);
det_Qz_fsc_Y1 = Qz_fsc_Y1(1,:).*Qz_fsc_Y1(2,:) - Qz_fsc_Y1(3,:).^2;

iQz_fsc_Y2 = [invQ0(1)+zQ_fsc_Y2; ...
              invQ0(2)+zQ_fsc_Y2; ...
              invQ0(3)*ones(1,NbP)];
det_iQz_fsc_Y2 = iQz_fsc_Y2(1,:).*iQz_fsc_Y2(2,:) - iQz_fsc_Y2(3,:).^2;
Qz_fsc_Y2 = [iQz_fsc_Y2(2,:); ...
          iQz_fsc_Y2(1,:); ...
         -iQz_fsc_Y2(3,:)]./(ones(3,1)*det_iQz_fsc_Y2);
det_Qz_fsc_Y2 = Qz_fsc_Y2(1,:).*Qz_fsc_Y2(2,:) - Qz_fsc_Y2(3,:).^2;

% vecteur polarisation aux points critiques du second ordre
h_l_x = -Champ.Coefficients(2).*ones(1,NbP);
h_l_y =  Champ.Coefficients(1).*ones(1,NbP);

h_l_X1_z = -Champ.Coefficients(1).*(xQ_fsc_X1.*Qz_fsc_X1(3,:) + yQ_fsc_X1.*Qz_fsc_X1(2,:)) ...
          +Champ.Coefficients(2).*(xQ_fsc_X1.*Qz_fsc_X1(1,:) + yQ_fsc_X1.*Qz_fsc_X1(3,:));
h_l_X1 = [h_l_x;h_l_y;h_l_X1_z];
h_X1   = M_fsc*h_l_X1;

h_l_X2_z = -Champ.Coefficients(1).*(xQ_fsc_X2.*Qz_fsc_X2(3,:) + yQ_fsc_X2.*Qz_fsc_X2(2,:)) ...
          +Champ.Coefficients(2).*(xQ_fsc_X2.*Qz_fsc_X2(1,:) + yQ_fsc_X2.*Qz_fsc_X2(3,:));
h_l_X2 = [h_l_x;h_l_y;h_l_X2_z];
h_X2   = M_fsc*h_l_X2;

h_l_Y1_z = -Champ.Coefficients(1).*(xQ_fsc_Y1.*Qz_fsc_Y1(3,:) + yQ_fsc_Y1.*Qz_fsc_Y1(2,:)) ...
          +Champ.Coefficients(2).*(xQ_fsc_Y1.*Qz_fsc_Y1(1,:) + yQ_fsc_Y1.*Qz_fsc_Y1(3,:));
h_l_Y1 = [h_l_x;h_l_y;h_l_Y1_z];
h_Y1   = M_fsc*h_l_Y1;

h_l_Y2_z = -Champ.Coefficients(1).*(xQ_fsc_Y2.*Qz_fsc_Y2(3,:) + yQ_fsc_Y2.*Qz_fsc_Y2(2,:)) ...
          +Champ.Coefficients(2).*(xQ_fsc_Y2.*Qz_fsc_Y2(1,:) + yQ_fsc_Y2.*Qz_fsc_Y2(3,:));
h_l_Y2 = [h_l_x;h_l_y;h_l_Y2_z];
h_Y2   = M_fsc*h_l_Y2;

% fonction de phase aux points critiques
g_X1 = -j/2*(rQ_X1(1,:).^2.*Axx + 2*rQ_X1(1,:).*rQ_X1(2,:).*Axy + rQ_X1(2,:).^2.*Ayy) ...
    + j*Bx.*rQ_X1(1,:) + j*By.*rQ_X1(2,:) ...
    + j*C;

g_X2 = -j/2*(rQ_X2(1,:).^2.*Axx + 2*rQ_X2(1,:).*rQ_X2(2,:).*Axy + rQ_X2(2,:).^2.*Ayy) ...
    + j*Bx.*rQ_X2(1,:) + j*By.*rQ_X2(2,:) ...
    + j*C ;

g_Y1 = -j/2*(rQ_Y1(1,:).^2.*Axx + 2*rQ_Y1(1,:).*rQ_Y1(2,:).*Axy + rQ_Y1(2,:).^2.*Ayy) ...
    + j*Bx.*rQ_Y1(1,:) + j*By.*rQ_Y1(2,:) ...
    + j*C;

g_Y2 = -j/2*(rQ_Y2(1,:).^2.*Axx + 2*rQ_Y2(1,:).*rQ_Y2(2,:).*Axy + rQ_Y2(2,:).^2.*Ayy) ...
    + j*Bx.*rQ_Y2(1,:) + j*By.*rQ_Y2(2,:) ...
    + j*C;


% derivees de la fonction de phase aux points critiques
dgdx_X1 = - j*rQ_X1(1,:).*Axx - j*rQ_X1(2,:).*Axy + j*Bx;
dgdx_X2 = - j*rQ_X2(1,:).*Axx - j*rQ_X2(2,:).*Axy + j*Bx;

dgdy_Y1 = - j*rQ_Y1(1,:).*Axy - j*rQ_Y1(2,:).*Ayy + j*By;
dgdy_Y2 = - j*rQ_Y2(1,:).*Axy - j*rQ_Y2(2,:).*Ayy + j*By;

% fonction d'amplitude au point col
ampl_X1 = 1./sum_r_X1;
ampl_X2 = 1./sum_r_X2;
ampl_Y1 = 1./sum_r_Y1;
ampl_Y2 = 1./sum_r_Y2;

% application de la formule 

% Lorsque le point critique du premier ordre n'appartient plus 
% au domaine d'integration, on supprime le coefficient 1/2
% car le point critique du premier ordre devient alors un 
% point du second ordre
Uid = (Plaque.X2<Centre_Phase(1)) & ...
      (Plaque.X1>Centre_Phase(1)) & ...
      (Plaque.Y2<Centre_Phase(2)) & ...
      (Plaque.Y1>Centre_Phase(2));

Cte = j*k/2/pi ;%.* (1 + Uid);

%  U_X1 = (real(y_X1)>Plaque.Y2) & (real(y_X1)<Plaque.Y1);
%  U_X2 = (real(y_X2)>Plaque.Y2) & (real(y_X2)<Plaque.Y1);
%  U_Y1 = (real(x_Y1)>Plaque.X2) & (real(x_Y1)<Plaque.X1);
%  U_Y2 = (real(x_Y2)>Plaque.X2) & (real(x_Y2)<Plaque.X1);
%  
U_X1 = (abs(y_X1)>Plaque.Y2) & (abs(y_X1)<Plaque.Y1);
U_X2 = (abs(y_X2)>Plaque.Y2) & (abs(y_X2)<Plaque.Y1);
U_Y1 = (abs(x_Y1)>Plaque.X2) & (abs(x_Y1)<Plaque.X1);
U_Y2 = (abs(x_Y2)>Plaque.X2) & (abs(x_Y2)<Plaque.X1);
%  
%  %  U_X1 = ones(1,NbP);
%  %  U_X2 = ones(1,NbP);
%  %  U_Y1 = ones(1,NbP);
%  %  U_Y2 = ones(1,NbP);
%  
%  ---------- Formulation non uniforme
E_X1 = cross(rr_X1, cross(rr_X1, cross(Plaque.N*ones(1,NbP), h_X1))) .* ...
      (ones(3,1)*(U_X1.*Cte .* ...
      sqrt(2*pi).*(-1).*... 
      ampl_X1 .* exp(-k*g_X1) ...
      ./k^(3/2) ...
      .*sqrt(det_Qz_fsc_X1./detQ0./(-j*Ayy))./dgdx_X1 ));

E_X2 = cross(rr_X2, cross(rr_X2, cross(Plaque.N*ones(1,NbP), h_X2))) .* ...
      (ones(3,1)*(U_X2.*Cte .* ...
      sqrt(2*pi)*(+1).*...
      ampl_X2 .* exp(-k*g_X2) ...
      ./k^(3/2) ...
      .*sqrt(det_Qz_fsc_X1./detQ0./(-j*Ayy))./dgdx_X2 ));

E_Y1 = cross(rr_Y1, cross(rr_Y1, cross(Plaque.N*ones(1,NbP), h_Y1))) .* ...
      (ones(3,1)*(U_Y1.*Cte .* ...
      sqrt(2*pi)*(-1).*...
      ampl_Y1 .* exp(-k*g_Y1) ...
      ./k^(3/2) ...
      .*sqrt(det_Qz_fsc_X1./detQ0./(-j*Axx))./dgdy_Y1 ));

E_Y2 = cross(rr_Y2, cross(rr_Y2, cross(Plaque.N*ones(1,NbP), h_Y2))) .* ...
      (ones(3,1)*(U_Y2.*Cte .* ...
      sqrt(2*pi)*(+1).*...
      ampl_Y2 .* exp(-k*g_Y2) ...
      ./k^(3/2) ...
      .*sqrt(det_Qz_fsc_X1./detQ0./(-j*Axx))./dgdy_Y2 ));

%  %  ----------  Formulation uniforme (valide pour la phase stat?)
%  nu_X1 = dgdx_X1.*sqrt(k/2./j*Axx);
%  nu_X2 = dgdx_X2.*sqrt(k/2./j*Axx);
%  nu_Y1 = dgdy_Y1.*sqrt(k/2./j*Ayy);
%  nu_Y2 = dgdy_Y2.*sqrt(k/2./j*Ayy);
%  %  nu_X1 = abs(dgdx_X1).*nsqrt(k./2./abs(-j*Axx));
%  %  nu_X2 = abs(dgdx_X2).*nsqrt(k./2./abs(-j*Axx));
%  %  nu_Y1 = abs(dgdy_Y1).*nsqrt(k./2./abs(-j*Ayy));
%  %  nu_Y2 = abs(dgdy_Y2).*nsqrt(k./2./abs(-j*Ayy));
%  
%  % Lorsque le point critique du premier ordre n'appartient plus 
%  % au domaine d'integration, on supprime le coefficient 1/2
%  % car le point critique du premier ordre devient alors un 
%  % point du second ordre
%  Cte = j*k/2/pi .* (1 + U1);
%  
%  E_X1 = cross(rr_X1, cross(rr_X1, cross(Plaque.N*ones(1,NbP), h_X1))) .* ...
%        (ones(3,1)*(U_X1 .* Cte .* ...
%        (+1)*2/k ...
%        .* ampl_X1 ...
%        .* exp(-k*g_X1) ...
%        .* sqrt(pi.*det_Qz_fsc_X1./detQ0./(-j*Axx)./(-j*Ayy)) ...
%        .* Fmin(nu_X1) ...
%        .* exp(-j*sign(real(Axx)).*nu_X1.^2) .*sign(-imag(dgdx_X1))));
%  
%  E_X2 = cross(rr_X2, cross(rr_X2, cross(Plaque.N*ones(1,NbP), h_X2))) .* ...
%        (ones(3,1)*(U_X2 .* Cte .* ...
%        (-1)*2/k ...
%        .* ampl_X2 ...
%        .* exp(-k*g_X2) ...
%        .* sqrt(pi.*det_Qz_fsc_X1./detQ0./(-j*Axx)./(-j*Ayy)) ...
%        .* Fmin(nu_X2) ...
%        .* exp(-j*sign(real(Axx)).*nu_X2.^2).*sign(-imag(dgdx_X2))));
%  
%  E_Y1 = cross(rr_Y1, cross(rr_Y1, cross(Plaque.N*ones(1,NbP), h_Y1))) .* ...
%        (ones(3,1)*(U_Y1 .* Cte .* ...
%        (+1)*2/k ...
%        .* ampl_Y1 ...
%        .* exp(-k*g_Y1) ...
%        .* sqrt(pi.*det_Qz_fsc_X1./detQ0./(-j*Axx)./(-j*Ayy)) ...
%        .* Fmin(nu_Y1) ...
%        .* exp(-j*sign(real(Ayy)).*nu_Y1.^2).*sign(real(dgdy_Y1))    ));
%  
%  E_Y2 = cross(rr_Y2, cross(rr_Y2, cross(Plaque.N*ones(1,NbP), h_Y2))) .* ...
%        (ones(3,1)*(U_Y2 .* Cte .* ...
%        (-1)*2/k ...
%        .* ampl_Y2 ...
%        .* exp(-k*g_Y2) ...
%        .* sqrt(pi.*det_Qz_fsc_X1./detQ0./(-j*Axx)./(-j*Ayy)) ...
%        .* Fmin(nu_Y2) ...
%        .* exp(-j*sign(real(Ayy)).*nu_Y2.^2).*sign(real(dgdy_Y2))  ));



% =======================================================
% Points critiques du troisieme ordre (coins)
% =======================================================

% vecteur distance aux points critiques du troisieme ordre
   rQ_X1Y1 = [Plaque.X1*ones(1,NbP); Plaque.Y1*ones(1,NbP);zeros(1,NbP)]; 
    r_X1Y1 = P - rQ_X1Y1;
sum_r_X1Y1 = nsqrt(sum(r_X1Y1.^2));
   rr_X1Y1 = r_X1Y1./(ones(3,1)*sum_r_X1Y1);

   rQ_X1Y2 = [Plaque.X1*ones(1,NbP); Plaque.Y2*ones(1,NbP);zeros(1,NbP)]; 
    r_X1Y2 = P - rQ_X1Y2;
sum_r_X1Y2 = nsqrt(sum(r_X1Y2.^2));
   rr_X1Y2 = r_X1Y2./(ones(3,1)*sum_r_X1Y2);

   rQ_X2Y1 = [Plaque.X2*ones(1,NbP); Plaque.Y1*ones(1,NbP);zeros(1,NbP)]; 
    r_X2Y1 = P - rQ_X2Y1;
sum_r_X2Y1 = nsqrt(sum(r_X2Y1.^2));
   rr_X2Y1 = r_X2Y1./(ones(3,1)*sum_r_X2Y1);

   rQ_X2Y2 = [Plaque.X2*ones(1,NbP); Plaque.Y2*ones(1,NbP);zeros(1,NbP)]; 
    r_X2Y2 = P - rQ_X2Y2;
sum_r_X2Y2 = nsqrt(sum(r_X2Y2.^2));
   rr_X2Y2 = r_X2Y2./(ones(3,1)*sum_r_X2Y2);

% point col dans le repere incident
xQ_fsc_X1Y1 = r11.*rQ_X1Y1(1,:) + r12.*rQ_X1Y1(2,:) - xi_Oi;
yQ_fsc_X1Y1 = r21.*rQ_X1Y1(1,:) + r22.*rQ_X1Y1(2,:) - yi_Oi;
zQ_fsc_X1Y1 = r31.*rQ_X1Y1(1,:) + r32.*rQ_X1Y1(2,:) - zi_Oi;

xQ_fsc_X1Y2 = r11.*rQ_X1Y2(1,:) + r12.*rQ_X1Y2(2,:) - xi_Oi;
yQ_fsc_X1Y2 = r21.*rQ_X1Y2(1,:) + r22.*rQ_X1Y2(2,:) - yi_Oi;
zQ_fsc_X1Y2 = r31.*rQ_X1Y2(1,:) + r32.*rQ_X1Y2(2,:) - zi_Oi;

xQ_fsc_X2Y1 = r11.*rQ_X2Y1(1,:) + r12.*rQ_X2Y1(2,:) - xi_Oi;
yQ_fsc_X2Y1 = r21.*rQ_X2Y1(1,:) + r22.*rQ_X2Y1(2,:) - yi_Oi;
zQ_fsc_X2Y1 = r31.*rQ_X2Y1(1,:) + r32.*rQ_X2Y1(2,:) - zi_Oi;

xQ_fsc_X2Y2 = r11.*rQ_X2Y2(1,:) + r12.*rQ_X2Y2(2,:) - xi_Oi;
yQ_fsc_X2Y2 = r21.*rQ_X2Y2(1,:) + r22.*rQ_X2Y2(2,:) - yi_Oi;
zQ_fsc_X2Y2 = r31.*rQ_X2Y2(1,:) + r32.*rQ_X2Y2(2,:) - zi_Oi;



% matrice de courbure dans le repere incident
iQz_fsc_X1Y1 = [invQ0(1)+zQ_fsc_X1Y1; ...
                invQ0(2)+zQ_fsc_X1Y1; ...
                invQ0(3)*ones(1,NbP)];
det_iQz_fsc_X1Y1 = iQz_fsc_X1Y1(1,:).*iQz_fsc_X1Y1(2,:) - iQz_fsc_X1Y1(3,:).^2;
Qz_fsc_X1Y1 = [iQz_fsc_X1Y1(2,:); ...
               iQz_fsc_X1Y1(1,:); ...
              -iQz_fsc_X1Y1(3,:)]./(ones(3,1)*det_iQz_fsc_X1Y1);
det_Qz_fsc_X1Y1 = Qz_fsc_X1Y1(1,:).*Qz_fsc_X1Y1(2,:) - Qz_fsc_X1Y1(3,:).^2;

iQz_fsc_X1Y2 = [invQ0(1)+zQ_fsc_X1Y2; ...
                invQ0(2)+zQ_fsc_X1Y2; ...
                invQ0(3)*ones(1,NbP)];
det_iQz_fsc_X1Y2 = iQz_fsc_X1Y2(1,:).*iQz_fsc_X1Y2(2,:) - iQz_fsc_X1Y2(3,:).^2;
Qz_fsc_X1Y2 = [iQz_fsc_X1Y2(2,:); ...
               iQz_fsc_X1Y2(1,:); ...
              -iQz_fsc_X1Y2(3,:)]./(ones(3,1)*det_iQz_fsc_X1Y2);
det_Qz_fsc_X1Y2 = Qz_fsc_X1Y2(1,:).*Qz_fsc_X1Y2(2,:) - Qz_fsc_X1Y2(3,:).^2;

iQz_fsc_X2Y1 = [invQ0(1)+zQ_fsc_X2Y1; ...
                invQ0(2)+zQ_fsc_X2Y1; ...
                invQ0(3)*ones(1,NbP)];
det_iQz_fsc_X2Y1 = iQz_fsc_X2Y1(1,:).*iQz_fsc_X2Y1(2,:) - iQz_fsc_X2Y1(3,:).^2;
Qz_fsc_X2Y1 = [iQz_fsc_X2Y1(2,:); ...
               iQz_fsc_X2Y1(1,:); ...
              -iQz_fsc_X2Y1(3,:)]./(ones(3,1)*det_iQz_fsc_X2Y1);
det_Qz_fsc_X2Y1 = Qz_fsc_X2Y1(1,:).*Qz_fsc_X2Y1(2,:) - Qz_fsc_X2Y1(3,:).^2;

iQz_fsc_X2Y2 = [invQ0(1)+zQ_fsc_X2Y2; ...
                invQ0(2)+zQ_fsc_X2Y2; ...
                invQ0(3)*ones(1,NbP)];
det_iQz_fsc_X2Y2 = iQz_fsc_X2Y2(1,:).*iQz_fsc_X2Y2(2,:) - iQz_fsc_X2Y2(3,:).^2;
Qz_fsc_X2Y2 = [iQz_fsc_X2Y2(2,:); ...
               iQz_fsc_X2Y2(1,:); ...
              -iQz_fsc_X2Y2(3,:)]./(ones(3,1)*det_iQz_fsc_X2Y2);
det_Qz_fsc_X2Y2 = Qz_fsc_X2Y2(1,:).*Qz_fsc_X2Y2(2,:) - Qz_fsc_X2Y2(3,:).^2;

% vecteur polarisation aux points critiques du second ordre
h_l_x = -Champ.Coefficients(2).*ones(1,NbP);
h_l_y =  Champ.Coefficients(1).*ones(1,NbP);

h_l_X1Y1_z = -Champ.Coefficients(1).*(xQ_fsc_X1Y1.*Qz_fsc_X1Y1(3,:) + yQ_fsc_X1Y1.*Qz_fsc_X1Y1(2,:)) ...
             +Champ.Coefficients(2).*(xQ_fsc_X1Y1.*Qz_fsc_X1Y1(1,:) + yQ_fsc_X1Y1.*Qz_fsc_X1Y1(3,:));
h_l_X1Y1 = [h_l_x;h_l_y;h_l_X1Y1_z];
h_X1Y1   = M_fsc*h_l_X1Y1;

h_l_X1Y2_z = -Champ.Coefficients(1).*(xQ_fsc_X1Y2.*Qz_fsc_X1Y2(3,:) + yQ_fsc_X1Y2.*Qz_fsc_X1Y2(2,:)) ...
             +Champ.Coefficients(2).*(xQ_fsc_X1Y2.*Qz_fsc_X1Y2(1,:) + yQ_fsc_X1Y2.*Qz_fsc_X1Y2(3,:));
h_l_X1Y2 = [h_l_x;h_l_y;h_l_X1Y2_z];
h_X1Y2   = M_fsc*h_l_X1Y2;

h_l_X2Y1_z = -Champ.Coefficients(1).*(xQ_fsc_X2Y1.*Qz_fsc_X2Y1(3,:) + yQ_fsc_X2Y1.*Qz_fsc_X2Y1(2,:)) ...
             +Champ.Coefficients(2).*(xQ_fsc_X2Y1.*Qz_fsc_X2Y1(1,:) + yQ_fsc_X2Y1.*Qz_fsc_X2Y1(3,:));
h_l_X2Y1 = [h_l_x;h_l_y;h_l_X2Y1_z];
h_X2Y1   = M_fsc*h_l_X2Y1;

h_l_X2Y2_z = -Champ.Coefficients(1).*(xQ_fsc_X2Y2.*Qz_fsc_X2Y2(3,:) + yQ_fsc_X2Y2.*Qz_fsc_X2Y2(2,:)) ...
             +Champ.Coefficients(2).*(xQ_fsc_X2Y2.*Qz_fsc_X2Y2(1,:) + yQ_fsc_X2Y2.*Qz_fsc_X2Y2(3,:));
h_l_X2Y2 = [h_l_x;h_l_y;h_l_X2Y2_z];
h_X2Y2   = M_fsc*h_l_X2Y2;

% fonction de phase aux points critiques
g_X1Y1 = -j/2*(rQ_X1Y1(1,:).^2.*Axx ...
            + 2*rQ_X1Y1(1,:).*rQ_X1Y1(2,:).*Axy ...
            + rQ_X1Y1(2,:).^2.*Ayy) ...
    + j*Bx.*rQ_X1Y1(1,:) + j*By.*rQ_X1Y1(2,:) ...
    + j*C;

g_X1Y2 = -j/2*(rQ_X1Y2(1,:).^2.*Axx ...
            + 2*rQ_X1Y2(1,:).*rQ_X1Y2(2,:).*Axy ...
            + rQ_X1Y2(2,:).^2.*Ayy) ...
    + j*Bx.*rQ_X1Y2(1,:) + j*By.*rQ_X1Y2(2,:) ...
    + j*C;

g_X2Y1 = -j/2*(rQ_X2Y1(1,:).^2.*Axx ...
            + 2*rQ_X2Y1(1,:).*rQ_X2Y1(2,:).*Axy ...
            + rQ_X2Y1(2,:).^2.*Ayy) ...
    + j*Bx.*rQ_X2Y1(1,:) + j*By.*rQ_X2Y1(2,:) ...
    + j*C;

g_X2Y2 = -j/2*(rQ_X2Y2(1,:).^2.*Axx ...
            + 2*rQ_X2Y2(1,:).*rQ_X2Y2(2,:).*Axy ...
            + rQ_X2Y2(2,:).^2.*Ayy) ...
    + j*Bx.*rQ_X2Y2(1,:) + j*By.*rQ_X2Y2(2,:) ...
    + j*C;

% derivees de la fonction de phase aux points critiques
dgdx_X1Y1 = -j*rQ_X1Y1(1,:).*Axx - j*rQ_X1Y1(2,:).*Axy + j*Bx;
dgdy_X1Y1 = -j*rQ_X1Y1(1,:).*Axy - j*rQ_X1Y1(2,:).*Ayy + j*By;

dgdx_X1Y2 = -j*rQ_X1Y2(1,:).*Axx - j*rQ_X1Y2(2,:).*Axy + j*Bx;
dgdy_X1Y2 = -j*rQ_X1Y2(1,:).*Axy - j*rQ_X1Y2(2,:).*Ayy + j*By;

dgdx_X2Y1 = -j*rQ_X2Y1(1,:).*Axx - j*rQ_X2Y1(2,:).*Axy + j*Bx;
dgdy_X2Y1 = -j*rQ_X2Y1(1,:).*Axy - j*rQ_X2Y1(2,:).*Ayy + j*By;

dgdx_X2Y2 = -j*rQ_X2Y2(1,:).*Axx - j*rQ_X2Y2(2,:).*Axy + j*Bx;
dgdy_X2Y2 = -j*rQ_X2Y2(1,:).*Axy - j*rQ_X2Y2(2,:).*Ayy + j*By;

% fonction d'amplitude au point col
ampl_X1Y1 = sqrt(det_Qz_fsc_X1Y1./detQ0) ...
            ./sum_r_X1Y1;

ampl_X1Y2 = sqrt(det_Qz_fsc_X1Y2./detQ0) ...
            ./sum_r_X1Y2;

ampl_X2Y1 = sqrt(det_Qz_fsc_X2Y1./detQ0) ...
            ./sum_r_X2Y1;

ampl_X2Y2 = sqrt(det_Qz_fsc_X2Y2./detQ0) ...
            ./sum_r_X2Y2;

%  % application de la formule non uniforme
Cte = j*k/2/pi;
E_X1Y1 = + cross(rr_X1Y1, cross(rr_X1Y1, cross(Plaque.N*ones(1,NbP), h_X1Y1))) .* ...
      (ones(3,1)*(Cte ...
        .* ampl_X1Y1 ...
        .* exp(-k*g_X1Y1) ...
        ./ k^2 ./ dgdx_X1Y1 ./ dgdy_X1Y1));
E_X1Y2 = - cross(rr_X1Y2, cross(rr_X1Y2, cross(Plaque.N*ones(1,NbP), h_X1Y2))) .* ...
      (ones(3,1)*(Cte ...
        .* ampl_X1Y2 ...
        .* exp(-k*g_X1Y2) ...
        ./ k^2 ./ dgdx_X1Y2 ./ dgdy_X1Y2));
E_X2Y1 = - cross(rr_X2Y1, cross(rr_X2Y1, cross(Plaque.N*ones(1,NbP), h_X2Y1))) .* ...
      (ones(3,1)*(Cte ...
        .* ampl_X2Y1 ...
        .* exp(-k*g_X2Y1) ...
        ./ k^2 ./ dgdx_X2Y1 ./ dgdy_X2Y1));
E_X2Y2 = + cross(rr_X2Y2, cross(rr_X2Y2, cross(Plaque.N*ones(1,NbP), h_X2Y2))) .* ...
      (ones(3,1)*(Cte ...
        .* ampl_X2Y2 ...
        .* exp(-k*g_X2Y2) ...
        ./ k^2 ./ dgdx_X2Y2 ./ dgdy_X2Y2));

%  % ===================================
%  %  version uniforme
%  % ===================================
%  nu_X1 = abs(dgdx_X1).*nsqrt(k/2./abs(-Axx));
%  nu_X2 = abs(dgdx_X2).*nsqrt(k/2./abs(-Axx));
%  nu_Y1 = abs(dgdy_Y1).*nsqrt(k/2./abs(-Ayy));
%  nu_Y2 = abs(dgdy_Y2).*nsqrt(k/2./abs(-Ayy));
%  T_X1 = sign(real(-Axx)) .* Fmin(nu_X1,7) .* exp(-j*sign(real(-Axx)).*nu_X1.^2);
%  T_X2 = sign(real(-Axx)) .* Fmin(nu_X2,7) .* exp(-j*sign(real(-Axx)).*nu_X2.^2);
%  T_Y1 = sign(real(-Ayy)) .* Fmin(nu_Y1,7) .* exp(-j*sign(real(-Ayy)).*nu_Y1.^2);
%  T_Y2 = sign(real(-Ayy)) .* Fmin(nu_Y2,7) .* exp(-j*sign(real(-Ayy)).*nu_Y2.^2);
%  
%  E_X1Y1 = cross(rr_X1Y1, cross(rr_X1Y1, cross(Plaque.N*ones(1,NbP), h_X1Y1))) .* ...
%        (ones(3,1)*(Cte ...
%          .* (-1) .* 2/k ...
%          .* ampl_X1Y1 ...
%          .* exp(j*k*g_X1Y1) ...
%          .* unSurSqrtAxx .* unSurSqrtAyy ...
%          ...%.* sign(real(dgdx_X1)).* sign(real(dgdy_Y1)) ...
%          .* T_X1 .* T_Y1        ));
%  
%  E_X1Y2 = cross(rr_X1Y2, cross(rr_X1Y2, cross(Plaque.N*ones(1,NbP), h_X1Y2))) .* ...
%        (ones(3,1)*(Cte ...
%          .* (+1) .* 2/k ...
%          .* ampl_X1Y1 ...
%          .* exp(j*k*g_X1Y2) ...
%          .* unSurSqrtAxx .* unSurSqrtAyy ...
%          .* sign(real(dgdx_X1)).* sign(real(dgdy_Y2)) ...
%          .* T_X1 .* T_Y2        ));
%  
%  E_X2Y1 = cross(rr_X2Y1, cross(rr_X2Y1, cross(Plaque.N*ones(1,NbP), h_X2Y1))) .* ...
%        (ones(3,1)*(Cte ...
%          .* (+1) .* 2/k ...
%          .* ampl_X2Y1 ...
%          .* exp(j*k*g_X2Y1) ...
%          .* unSurSqrtAxx .* unSurSqrtAyy ...
%          .* sign(real(dgdx_X2)).* sign(real(dgdy_Y1)) ...
%          .* T_X2 .* T_Y1        ));
%  
%  E_X2Y2 = cross(rr_X2Y2, cross(rr_X2Y2, cross(Plaque.N*ones(1,NbP), h_X2Y2))) .* ...
%        (ones(3,1)*(Cte ...
%          .* (-1) .* 2/k ...
%          .* ampl_X2Y1 ...
%          .* exp(j*k*g_X2Y2) ...
%          .* unSurSqrtAxx .* unSurSqrtAyy ...
%          .* sign(real(dgdx_X2)).* sign(real(dgdy_Y2)) ...
%          .* T_X2 .* T_Y2        ));

% ==========================================================
% champ total
% ==========================================================
%  filtrage des NaN
E_X1(isnan(E_X1)) = eps;
E_X2(isnan(E_X2)) = eps;
E_Y1(isnan(E_Y1)) = eps;
E_Y2(isnan(E_Y2)) = eps;
E_X1Y1(isnan(E_X1Y1)) = eps;
E_X1Y2(isnan(E_X1Y2)) = eps;
E_X2Y1(isnan(E_X2Y1)) = eps;
E_X2Y2(isnan(E_X2Y2)) = eps;
Epc1(isnan(Epc1)) = eps;

Epc2 = E_X1 + E_X2  + E_Y1 + E_Y2;
Epc3 = E_X1Y1 +  E_X1Y2 + E_X2Y1 + E_X2Y2;
E = Epc1 + Epc2 + Epc3;

H = zeros(size(E));