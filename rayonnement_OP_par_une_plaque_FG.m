% Rayonnement par une plaque conductrice d'un faisceau gaussien incident
% dans l'hypothese de l'optique physique
% Formulation paraxiale du FG incident
clear all;

% =========================================================
% Definition generale
% =========================================================
     c = 2.997925e8;
     f = 10e9;
    Z0 = 120*pi;
lambda = c/f;
     k = 2*pi /lambda;
  epsr = 1;

%  Afficher des champs normalise au max ?
bool_normalisation = false;

%  Quelle polarisation afficher ?
%  1 : Etheta
%  2 : Ephi
%  3 : Er
id_pol = 1;

%  Quelle type de formulation employer ?
%  1 : non uniforme
%  2 : uniforme
type_formulation = 1;

%  faisceau gaussien paraxial ? 
%  true : paraxiale
%  false: champ lointain
bool_paraxial = true;

% repere absolu cartesien : repere de la plaque 
rep_plaque.centre  = [0;0;0]*lambda;
rep_plaque.cart.ex = [1;0;0];
rep_plaque.cart.ey = [0;1;0];
rep_plaque.cart.ez = [0;0;1];

% =========================================================
% Definition des points d'observation
% =========================================================
rep_obs.R     = 1000*lambda;
rep_obs.theta = linspace(-pi/2, pi/2,1000);
rep_obs.phi   = deg2rad(0);
P             = rep_obs.R*[cos(rep_obs.phi)*sin(rep_obs.theta); ...
                         sin(rep_obs.phi)*sin(rep_obs.theta); ...
                         cos(rep_obs.theta)];
NbP = length(P);
R = sqrt(sum(P.^2));



% =========================================================
% Definition des coordonnees et tailles de la plaque
% =========================================================
DemiCote =10*lambda;

Sommets = [-DemiCote, +DemiCote, +DemiCote, -DemiCote;...
            DemiCote, DemiCote, -DemiCote, -DemiCote;...
            0, 0, 0, 0];
Epsr_plaque = 10000000000;

Plaque = struct_Plaque(Sommets, [0;0;Sommets(3,1)]*lambda, [0;0;1], Epsr_plaque);

X1 = max(Sommets(1,:));
X2 = min(Sommets(1,:));
Y1 = max(Sommets(2,:));
Y2 = min(Sommets(2,:));

% maillage de la plaque pour le calcul
% du champ reference (kottler)
dx = lambda/13;
dy = lambda/13;

[Q,Nx,Ny,Nz] = maillage_plaque(Plaque, dx, dy);



% =========================================================
% Definition des propriete du faisceau gaussien incident
% =======================================================
    R_inc = 10*lambda;
theta_inc = deg2rad(0);
  phi_inc = deg2rad(0);
rep_inc.centre  = R_inc*[cos(phi_inc)*sin(theta_inc); ...
                         sin(phi_inc)*sin(theta_inc); ...
                         cos(theta_inc)];

rep_inc.centre  = [10;0;10]*lambda; % centre du repere incident

% angle de rotation (spheriques) pour passer du repere
% plaque au repere incident
theta = pi+deg2rad(0);
phi   = deg2rad(0);

% definition de la matrice de rotation pour passer
% du repere de la plaque au repere incident
R_phi = [cos(phi), sin(phi), 0; ...
        -sin(phi), cos(phi), 0; ...
                0,        0, 1];
R_theta = [cos(theta), 0, -sin(theta); ...
                    0, 1, 0 ;...
           sin(theta), 0, cos(theta)];
R_theta_phi = R_theta*R_phi;

% creation du repere fscal au faisceau
rep_inc.cart.ex = transpose(R_theta_phi)*rep_plaque.cart.ex; % e_prll^i
rep_inc.cart.ey = transpose(R_theta_phi)*rep_plaque.cart.ey; % e_perp^i
rep_inc.cart.ez = transpose(R_theta_phi)*rep_plaque.cart.ez; % e_s^i

% proprietes du faisceau
         kW0 = 4*pi;
%            W0 = 2*lambda;%2/(k*deg2rad(15))
          W0 = kW0/k;
Courbure     = [1;1;0]/(j*k*W0^2/2);
Coefficients = [-1;0]; % co et cross polar

% creation de la structure Champ
Champ = struct_Champ(rep_inc.centre, ...
                     rep_inc.cart.ex, ...
                     rep_inc.cart.ez, ...
                     Coefficients, ...
                     Courbure);

% =========================================================
% Visualisation de la geometrie de la scene
% =========================================================
clf(figure(1))
  h=plot3(Q(1,:)/lambda, Q(2,:)/lambda, Q(3,:)/lambda-.7, '.', 'Color', [.9 .9 .9])
  xlabel('x (\lambda)', 'FontSize', 20);
  ylabel('y (\lambda)', 'FontSize', 20);
  zlabel('z (\lambda)', 'FontSize', 20);
%    title('Geometrie', 'FontSize', 20)
  view(-30,25);
  grid on;
  
  hold on

%  representation de la plaque dans l'espace
  line([Plaque.Sommets(1,:), Plaque.Sommets(1,1)]/lambda, ...
      [Plaque.Sommets(2,:), Plaque.Sommets(2,1)]/lambda, ...
      [Plaque.Sommets(3,:), Plaque.Sommets(3,1)]/lambda);

% centre de la plaque
  line(Plaque.Centre(1)/lambda, Plaque.Centre(2)/lambda, Plaque.Centre(3)/lambda, ...
      'color', 'r', 'markersize', 20, 'marker', '.'); 

% centre du repere incident
  line(rep_inc.centre(1)/lambda, rep_inc.centre(2)/lambda, rep_inc.centre(3)/lambda, ...
      'color', 'b', 'markersize', 20, 'marker', '.'); 

%  affiche un cercle de rayon W0 dans le plan fsc 
  alpha=linspace(0,2*pi,100);
  cercle = transpose(R_theta_phi)*[W0.*cos(alpha); W0.*sin(alpha) ; zeros(size(alpha))] ...
          + rep_inc.centre*ones(size(alpha))  ;
  line(cercle(1,:)/lambda, cercle(2,:)/lambda, cercle(3,:)/lambda, 'Color', 'r')
  text(rep_inc.centre(1)/lambda, rep_inc.centre(2)/lambda, rep_inc.centre(3)/lambda+5, ['W_0=', num2str(W0/lambda), ' \lambda'], 'Color', 'r');

% vecteur unitaires du repere fscal (dans le repere absolu)
hold on;
  Cte = 10;
  if (Coefficients(1)~=0) % e parallel
    C = Cte * Coefficients(1);
    quiver3(...
          rep_inc.centre(1)/lambda, rep_inc.centre(2)/lambda, rep_inc.centre(3)/lambda,...
          C*rep_inc.cart.ex(1), C*rep_inc.cart.ex(2), C*rep_inc.cart.ex(3), 'LineWidth', 2, 'Color', 'r');
    text(rep_inc.centre(1)/lambda+C*rep_inc.cart.ex(1), ...
        rep_inc.centre(2)/lambda+C*rep_inc.cart.ex(2), ...
        rep_inc.centre(3)/lambda+C*rep_inc.cart.ex(3)+2, 'ê_{fg;||}', 'FontSize', 18, 'Color', 'r');
  end
  if (Coefficients(2)~=0) % e perpendiculaire
    C = Cte * Coefficients(2);
    quiver3(...
          rep_inc.centre(1)/lambda, rep_inc.centre(2)/lambda, rep_inc.centre(3)/lambda,...
          C*rep_inc.cart.ey(1), C*rep_inc.cart.ey(2), C*rep_inc.cart.ey(3), 'LineWidth', 2, 'Color', 'r');
    text(rep_inc.centre(1)/lambda+C*rep_inc.cart.ey(1), ...
        rep_inc.centre(2)/lambda+C*rep_inc.cart.ey(2), ...
        rep_inc.centre(3)/lambda+C*rep_inc.cart.ey(3)+2, 'ê_{fg;\perp}', 'FontSize', 18, 'Color', 'r');
  end
  % e_s (direction de propa FG)
  quiver3(...
        rep_inc.centre(1)/lambda, rep_inc.centre(2)/lambda, rep_inc.centre(3)/lambda,...
        Cte*rep_inc.cart.ez(1), Cte*rep_inc.cart.ez(2), Cte*rep_inc.cart.ez(3), 'LineWidth', 2);
  text(rep_inc.centre(1)/lambda+Cte*rep_inc.cart.ez(1)+1, ...
      rep_inc.centre(2)/lambda+Cte*rep_inc.cart.ez(2)+1, ...
      rep_inc.centre(3)/lambda+Cte*rep_inc.cart.ez(3)+2, 'ê_{fg;z}', 'FontSize', 18, 'Color', 'r');
% on trace le prolongement du vecteur directeur e_s^i jusqu'au plan z=0
  Centre_Phase_x = rep_inc.cart.ez(1)*(-rep_inc.centre(3)/rep_inc.cart.ez(3))+rep_inc.centre(1);
  Centre_Phase_y = rep_inc.cart.ez(2)*(-rep_inc.centre(3)/rep_inc.cart.ez(3))+rep_inc.centre(2);
  line([rep_inc.centre(1)/lambda Centre_Phase_x/lambda], ...
      [rep_inc.centre(2)/lambda Centre_Phase_y/lambda], ...
      [rep_inc.centre(3)/lambda 0], 'Color', 'r', 'LineStyle', ':');


% vecteur unitaires du repere absolu
  quiver3(...
        rep_plaque.centre(1)/lambda, rep_plaque.centre(2)/lambda, rep_plaque.centre(3)/lambda,...
        Cte*rep_plaque.cart.ex(1), Cte*rep_plaque.cart.ex(2), Cte*rep_plaque.cart.ex(3), 'LineWidth', 2, 'Color', 'k');
  text(Cte*rep_plaque.cart.ex(1), Cte*rep_plaque.cart.ex(2), Cte*rep_plaque.cart.ex(3)+2, 'ê_x', 'FontSize', 18);
  quiver3(...
        rep_plaque.centre(1)/lambda, rep_plaque.centre(2)/lambda, rep_plaque.centre(3)/lambda,...
        Cte*rep_plaque.cart.ey(1), Cte*rep_plaque.cart.ey(2), Cte*rep_plaque.cart.ey(3), 'LineWidth', 2, 'Color', 'k');
  text(Cte*rep_plaque.cart.ey(1), Cte*rep_plaque.cart.ey(2), Cte*rep_plaque.cart.ey(3)+2, 'ê_y', 'FontSize', 18);
  quiver3(...
        rep_plaque.centre(1)/lambda, rep_plaque.centre(2)/lambda, rep_plaque.centre(3)/lambda,...
        Cte*rep_plaque.cart.ez(1), Cte*rep_plaque.cart.ez(2), Cte*rep_plaque.cart.ez(3), 'LineWidth', 2, 'MaxHeadSize', 1, 'Color', 'k');
  text(Cte*rep_plaque.cart.ez(1)+1, Cte*rep_plaque.cart.ez(2)+1, Cte*rep_plaque.cart.ez(3)-4, 'ê_z', 'FontSize', 18);
  
  %  Visualisation des points d'observation en zone lointaine
  Cte2 = 30;
  hP=line(Cte2*cos(rep_obs.phi)*sin(rep_obs.theta), Cte2*sin(rep_obs.phi)*sin(rep_obs.theta), Cte2*cos(rep_obs.theta), 'Color', 'k');
  
  hold off;
  axis equal

%  %  on sauvegarde la geometrie
%  saveas(figure(1), 'geometrie_rayonnement_plaque.png');



%  #############################################################
%  #############################################################
%  #############################################################
%  #########################  CALCULS  #########################
%  #############################################################
%  #############################################################
%  #############################################################
disp('##### Geometrie correcte ? (appuyer pour continuer) #####')
pause



%  ==================================================================
%  Calcul du champ asymptotique rayonne par les courants electriques
%  dans l'approximation de l'optique physique
%  ==================================================================
if type_formulation == 1
  [Er_fg,Hr_fg, Epc1, Epc2, Epc3] = EH_FG_OP_plaque_conductrice(Champ, Plaque, P, f, epsr);
elseif type_formulation == 2
  if bool_paraxial
  [Er_fg,Hr_fg, Epc1, Epc2, Epc3] = EH_FG_OP_plaque_conductrice_uniforme(Champ, Plaque, P, f, epsr);
%    [Er_fg,Hr_fg, Epc1, Epc2, Epc3] = EH_FG_OP_plaque_conductrice_uniforme_v2(Champ, Plaque, P, f, epsr);
  else
  [Er_fg,Hr_fg, Epc1, Epc2, Epc3] = EH_FG_lointain_OP_plaque_conductrice_uniforme(Champ, Plaque, P, f, epsr);
  end
end

%  =========================================================
%  Calcul du champ incident sur la plaque 
%  (pour kottler et pour la visualisation)
%  =========================================================
if bool_paraxial 
  [Ei_Q, Hi_Q] = EH_Sigma_FG(Q, Champ, f, epsr);
else
  [Ei_Q, Hi_Q] = EH_Sigma_FG_chplointain(Q, Champ, f, epsr);
end

% =========================================================
% Visualisation des champs sur la plaque
% =========================================================
Qx = stackedfield2field(Q(1,:), Nx,Ny,Nz);
Qy = stackedfield2field(Q(2,:), Nx,Ny,Nz);

Ei_Qx = stackedfield2field(Ei_Q(1,:), Nx,Ny,Nz);
Ei_Qy = stackedfield2field(Ei_Q(2,:), Nx,Ny,Nz);
Ei_Qz = stackedfield2field(Ei_Q(3,:), Nx,Ny,Nz);
Ei_Qi = {20*log10(abs(Ei_Qx)), 20*log10(abs(Ei_Qy)), 20*log10(abs(Ei_Qz)),...
        180/pi*angle(Ei_Qx) , 180/pi*angle(Ei_Qy) , 180/pi*angle(Ei_Qz)};

figure(1)
  set(gca, 'FontSize', 15);
  delete(h)
  hold on
  pcolor(Qx/lambda, Qy/lambda, Ei_Qi{1})
  shading interp, caxis([-80 0]);
  h=colorbar
  set(get(h, 'Ylabel'), 'String', 'Amplitude du champ en dB', 'FontSize', 15);
  
  %  visualisation du point central d'impact du faisceau sur la plaque
  line(Centre_Phase_x/lambda,Centre_Phase_y/lambda,0+.1,'Color', 'k', 'Marker', '+', 'MarkerSize', 20);

if ~exist('save.EHkottler.mat')
  % =========================================================
  % calcul du courant electrique sur la plaque
  % =========================================================
  J = 2*cross(Plaque.N*ones(1,length(Hi_Q)), Hi_Q);

  % =========================================================
  % calcul du champ rayonne
  % =========================================================
        dS = dx*dy*ones(1,length(Q));
  [Er_k] = RayonnementKottler(P,Q,dS,J,zeros(size(J)),f,epsr,3);
  Hr_k = [];

  save('save.EHkottler.mat', 'Er_k', 'Hr_k', 'J');
else
  load('save.EHkottler.mat')
end

%  =========================================================
%    Calcul du champ rayonne derriere la plaque (total)
%  =========================================================
% Calcul du champ incident aux points d'observation
% pour la recomposition avec l'OP
[Ei,Hi] = EH_Sigma_FG(P, Champ, f, epsr);

%  En champ lointain, on ajoute le champ incident
%  uniquement derriï¿½re la plaque pour ne pas avoir
%  le champ incident superposï¿½ au champ reflechi.
u = ones(3,1)*(P(3,:)<0);

 Et_k = Er_k  + u.*Ei;
Et_fg = Er_fg + u.*Ei;


%  =========================================================
%  Calcul des composantes /theta ou phi selon le choix
%  =========================================================
switch id_pol 
  case 1 % Etheta
    e_theta = [cos(rep_obs.theta).*cos(rep_obs.phi); ...
              cos(rep_obs.theta).*sin(rep_obs.phi); ...
              -sin(rep_obs.theta)];
    
    plot_1 = dot(Et_k, e_theta);
    plot_2 = dot(Et_fg,e_theta);
    plot_i = dot(Ei, e_theta);

    plot_pc1 = dot(Epc1, e_theta);
    plot_pc2 = dot(Epc2, e_theta);
    plot_pc3 = dot(Epc3, e_theta);
    legend1 = 'OP Numerique E_\theta';
    legend2 = 'OP Asymptotique E_\theta (FG inc paraxial))';

  case 2 % Phi
    e_phi = [-sin(rep_obs.phi);cos(rep_obs.phi);0]*ones(1,length(rep_obs.theta));
    
    plot_1 = dot(Et_k, e_phi);
    plot_2 = dot(Et_fg,e_phi);
    plot_i = dot(Ei, e_phi);

    plot_pc1 = dot(Epc1, e_phi);
    plot_pc2 = dot(Epc2, e_phi);
    plot_pc3 = dot(Epc3, e_phi);
    legend1 = 'OP Numerique E_\phi';
    legend2 = 'OP Asymptotique E_\phi (FG inc paraxial))';

  case 3 % R
    e_r = [sin(rep_obs.theta).*cos(rep_obs.phi); ...
           sin(rep_obs.theta).*sin(rep_obs.phi); ...
           cos(rep_obs.theta)];
    
    plot_1 = dot(Et_k, e_r);
    plot_2 = dot(Et_fg,e_r);
    plot_i = dot(Ei, e_r);

    plot_pc1 = dot(Epc1, e_r);
    plot_pc2 = dot(Epc2, e_r);
    plot_pc3 = dot(Epc3, e_r);
    legend1 = 'OP Numerique E_\r';
    legend2 = 'OP Asymptotique E_\r (FG inc paraxial))';    
end

%  coefficient de normalisation
if bool_normalisation
  no = 1./max(plot1);
else
  no = 1;
end


%  =========================================================
%  REPRESENTATION
%  =========================================================
conf_axis = [min(rad2deg(rep_obs.theta)) max(rad2deg(rep_obs.theta)) -140 -30];


%  AMPLITUDE
%  
clf(figure(2))
  set(gca, 'FontSize', 15);
  plot(rad2deg(rep_obs.theta), 20*log10(abs(no.*plot_1)), ...
       rad2deg(rep_obs.theta), 20*log10(abs(no.*plot_2)), 'LineWidth', 2) ;
  hold on;
  plot(rad2deg(rep_obs.theta), 20*log10(abs(no.*plot_1 - no.*plot_2)), 'r--');
  hold off;
  legend(legend1, legend2, 'Difference', 'Location', 'SouthOutside');
  xlabel('\theta (deg)');
  ylabel('amplitude (dB)');
  grid on;
  axis(conf_axis);

%  PHASE
%  
clf(figure(3))
  set(gca, 'FontSize', 15);
  plot(rad2deg(rep_obs.theta), 180/pi*(angle(plot_1)), ...
       rad2deg(rep_obs.theta), 180/pi*(angle(plot_2)), 'LineWidth', 2);
  legend(legend1, legend2, 'Location', 'SouthOutside');
  xlabel('\theta (deg)');
  ylabel('phase (deg)');
  grid on;
  axis([min(rad2deg(rep_obs.theta)) max(rad2deg(rep_obs.theta)) -180 180]);


%  Illustration des apports de chacun des point critiques
%  ou de chacunes des types d'integrale 
%  selon la formulation utilisee

clf(figure(4))
  set(gca, 'FontSize', 15);
  plot(rad2deg(rep_obs.theta), 20*log10(abs(plot_1)), '-', ...
  rad2deg(rep_obs.theta),  20*log10(abs([plot_pc1;plot_pc2;plot_pc3])), 'LineWidth', 2);
  grid on;
if type_formulation == 1
  legend('Reference', 'Contr. pts. prem. ordre', 'Contr. pts. snd. ordre', 'Contr. pts. trois. ordre', 'Location', 'SouthOutside');
elseif type_formulation == 2
  legend('Reference', 'Contr. I_0', 'Contr. I_{1,2}', 'Contr. I_{ij}', 'Location', 'SouthOutside');
end
  xlabel('\theta (deg)');
  ylabel('Amplitude (dB)');
  grid on;
  axis([min(rad2deg(rep_obs.theta)) max(rad2deg(rep_obs.theta)) -150 -30]);

%  %  représentation de l'erreur
%  clf(figure(5))
%    set(gca, 'FontSize', 15);
%    plot(rad2deg(rep_obs.theta), 20*log10(abs(plot_1 - plot_2)./max(abs(plot_1))));
%    legend('Difference', 'Location', 'SouthOutside');
%    xlabel('\theta (deg)');
%    ylabel('amplitude (dB)');
%    grid on;
%    axis([min(rad2deg(rep_obs.theta)) max(rad2deg(rep_obs.theta)) -120 -40]);




