% Rayonnement par une plaque conductrice d'un faisceau gaussien incident
% dans l'hypothese de l'optique physique
% Formulation champ lointain du FG incident
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
%  index de polarisation pour l'affichage
  id_pol = 1;

%  normalisation des resultats au max
  bool_normalisation = false;

% repere absolu cartesien : repere de la plaque 
rep_plaque.centre  = [0;0;0]*lambda;
rep_plaque.cart.ex = [1;0;0];
rep_plaque.cart.ey = [0;1;0];
rep_plaque.cart.ez = [0;0;1];


% =========================================================
% Definition des points d'observation
% =========================================================
rep_obs.R     = 1000*lambda;
rep_obs.theta = linspace(-pi/2,pi/2,1000);
rep_obs.phi   = deg2rad(0);
P             = rep_obs.R*[cos(rep_obs.phi)*sin(rep_obs.theta); ...
                         sin(rep_obs.phi)*sin(rep_obs.theta); ...
                         cos(rep_obs.theta)];
NbP = length(P);
R = sqrt(sum(P.^2));



% =========================================================
% Definition des coordonnees et tailles de la plaque
% =========================================================
Sommets = [-10, +10, +10, -10;...
            10, 10, -10, -10;...
            0, 0, 0, 0]*lambda;


Epsr_plaque = 10000000000;

Plaque = struct_Plaque(Sommets, [0;0;Sommets(3,1)]*lambda, [0;0;1], Epsr_plaque);

X1 = max(Sommets(1,:));
X2 = min(Sommets(1,:));
Y1 = max(Sommets(2,:));
Y2 = min(Sommets(2,:));

% maillage de la plaque pour le calcul
% du champ reference (kottler)
dx = lambda/12;
dy = lambda/12;

[Q,Nx,Ny,Nz] = maillage_plaque(Plaque, dx, dy);



% =========================================================
% Definition des propriete du faisceau gaussien incident
% =======================================================
rep_inc.centre  = [0;0;10]*lambda; % centre du repere incident

% angle de rotation (spheriques) pour passer du repere
% plaque au repere incident
theta = pi ;
phi   = 0;

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
W0           = 4*lambda;%2/(k*deg2rad(15))
Courbure     = [1;1;0]/(j*k*W0^2/2);
Coefficients = [1;0]; % co et cross polar

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
  title('Géometrie de la scène', 'FontSize', 20)
  view(-30,25);
  grid on;

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
%  %  ligne reliant le centre du repere fsc avec le centre du repere plaque
%  line([rep_inc.centre(1) rep_plaque.centre(1)]/lambda, ...
%       [rep_inc.centre(2) rep_plaque.centre(2)]/lambda, ...
%       [rep_inc.centre(3) rep_plaque.centre(3)]/lambda, ...
%      'Color', 'k');

%  affiche un cercle de rayon W0 dans le plan fsc 
alpha=linspace(0,2*pi,100);
cercle = transpose(R_theta_phi)*[W0.*cos(alpha); W0.*sin(alpha) ; zeros(size(alpha))] ...
        + rep_inc.centre*ones(size(alpha))  ;
line(cercle(1,:)/lambda, cercle(2,:)/lambda, cercle(3,:)/lambda, 'Color', 'r')
text(rep_inc.centre(1)/lambda, rep_inc.centre(2)/lambda, rep_inc.centre(3)/lambda+2, ['W_0=', num2str(W0/lambda), ' \lambda'], 'Color', 'r');

% vecteur unitaires du repere fscal (dans le repere absolu)
hold on;
Cte = 10;
if (Coefficients(1)~=0) % e parallel
  C = Cte * Coefficients(1);
  quiver3(...
        rep_inc.centre(1)/lambda, rep_inc.centre(2)/lambda, rep_inc.centre(3)/lambda,...
        C*rep_inc.cart.ex(1), C*rep_inc.cart.ex(2), C*rep_inc.cart.ex(3), 'LineWidth', 2);
  text(rep_inc.centre(1)/lambda+C*rep_inc.cart.ex(1), ...
      rep_inc.centre(2)/lambda+C*rep_inc.cart.ex(2), ...
      rep_inc.centre(3)/lambda+C*rep_inc.cart.ex(3), 'e_{||}^i', 'FontSize', 18);
end
if (Coefficients(2)~=0) % e perpendiculaire
  C = Cte * Coefficients(2);
  quiver3(...
        rep_inc.centre(1)/lambda, rep_inc.centre(2)/lambda, rep_inc.centre(3)/lambda,...
        C*rep_inc.cart.ey(1), C*rep_inc.cart.ey(2), C*rep_inc.cart.ey(3), 'LineWidth', 2);
  text(rep_inc.centre(1)/lambda+C*rep_inc.cart.ey(1), ...
      rep_inc.centre(2)/lambda+C*rep_inc.cart.ey(2), ...
      rep_inc.centre(3)/lambda+C*rep_inc.cart.ey(3), 'e_{\perp}^i', 'FontSize', 18);
end
% e_s (direction de propa FG)
quiver3(...
      rep_inc.centre(1)/lambda, rep_inc.centre(2)/lambda, rep_inc.centre(3)/lambda,...
      Cte*rep_inc.cart.ez(1), Cte*rep_inc.cart.ez(2), Cte*rep_inc.cart.ez(3), 'LineWidth', 2);
text(rep_inc.centre(1)/lambda+Cte*rep_inc.cart.ez(1), ...
     rep_inc.centre(2)/lambda+Cte*rep_inc.cart.ez(2), ...
     rep_inc.centre(3)/lambda+Cte*rep_inc.cart.ez(3), 'e_s^i', 'FontSize', 18);
% on trace le prolongement du vecteur directeur e_s^i jusqu'au plan z=0
Centre_Phase_x = rep_inc.cart.ez(1)*(-rep_inc.centre(3)/rep_inc.cart.ez(3))+rep_inc.centre(1);
Centre_Phase_y = rep_inc.cart.ez(2)*(-rep_inc.centre(3)/rep_inc.cart.ez(3))+rep_inc.centre(2);
line([rep_inc.centre(1)/lambda Centre_Phase_x/lambda], ...
     [rep_inc.centre(2)/lambda Centre_Phase_y/lambda], ...
     [rep_inc.centre(3)/lambda 0], 'Color', 'r');


% vecteur unitaires du repere absolu
quiver3(...
      rep_plaque.centre(1)/lambda, rep_plaque.centre(2)/lambda, rep_plaque.centre(3)/lambda,...
      Cte*rep_plaque.cart.ex(1), Cte*rep_plaque.cart.ex(2), Cte*rep_plaque.cart.ex(3), 'LineWidth', 2);
text(Cte*rep_plaque.cart.ex(1), Cte*rep_plaque.cart.ex(2), Cte*rep_plaque.cart.ex(3)+2, 'e_x', 'FontSize', 18);
quiver3(...
      rep_plaque.centre(1)/lambda, rep_plaque.centre(2)/lambda, rep_plaque.centre(3)/lambda,...
      Cte*rep_plaque.cart.ey(1), Cte*rep_plaque.cart.ey(2), Cte*rep_plaque.cart.ey(3), 'LineWidth', 2);
text(Cte*rep_plaque.cart.ey(1), Cte*rep_plaque.cart.ey(2), Cte*rep_plaque.cart.ey(3)+2, 'e_y', 'FontSize', 18);
quiver3(...
      rep_plaque.centre(1)/lambda, rep_plaque.centre(2)/lambda, rep_plaque.centre(3)/lambda,...
      Cte*rep_plaque.cart.ez(1), Cte*rep_plaque.cart.ez(2), Cte*rep_plaque.cart.ez(3), 'LineWidth', 2, 'MaxHeadSize', 1);
text(Cte*rep_plaque.cart.ez(1), Cte*rep_plaque.cart.ez(2), Cte*rep_plaque.cart.ez(3), 'e_z', 'FontSize', 18);

%  Visualisation des points d'observation en zone lointaine
Cte2 = 30;
hP=line(Cte2*cos(rep_obs.phi)*sin(rep_obs.theta), Cte2*sin(rep_obs.phi)*sin(rep_obs.theta), Cte2*cos(rep_obs.theta), 'Color', 'k');

hold off;
axis equal
 



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
tic;
[Er_fg,Hr_fg, Epc1, Epc2, Epc3] = EH_FG_lointain_OP_plaque_conductrice_uniforme(Champ, Plaque, P, f, epsr);
t_ana = toc

%  [Er_fg2,Hr_fg2, Epc12, Epc22, Epc32] = EH_FG_OP_plaque_conductrice_uniforme(Champ, Plaque, P, f, epsr);
%  =========================================================
%  Calcul du champ incident sur la plaque 
%  (pour kottler et pour la visualisation)
%  =========================================================
[Ei_Q, Hi_Q] = EH_Sigma_FG_chplointain(Q, Champ, f, epsr);


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
  delete(h)
  hold on
  pcolor(Qx/lambda, Qy/lambda, Ei_Qi{1})
  shading interp, caxis([-80 0]);
  colorbar
  %  visualisation du point central d'impact du faisceau sur la plaque
  line(Centre_Phase_x/lambda,Centre_Phase_y/lambda,0+.1,'Color', 'k', 'Marker', '+', 'MarkerSize', 20);

if ~exist('save.EHkottler_loinntain.mat')
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

  save('save.EHkottler_loinntain.mat', 'Er_k', 'Hr_k', 'J');
else
  load('save.EHkottler_loinntain.mat')
end




conf_axis = [min(rad2deg(rep_obs.theta)) max(rad2deg(rep_obs.theta)) -150 0];

%  =========================================================
%    Calcul du champ incident
%  =========================================================
[Ei,Hi] = EH_Sigma_FG_chplointain(P, Champ, f, epsr);

%  =========================================================
%    Calcul du champ rayonne derriere la plaque (total)
%  =========================================================
%  En champ lointain, on ajoute le champ incident
%  uniquement derrière la plaque pour ne pas avoir
%  le champ incident superposé au champ reflechi.
u = ones(3,1)*(P(3,:)<0);

 Et_k = Er_k  + u.*Ei;
Et_fg = Er_fg + u.*Ei;
%  Et_fg2 = Er_fg2 + u.*Ei;

%  AMPLITUDE
clf(figure(6))
plot(rad2deg(rep_obs.theta), 20*log10(abs(Et_k)), '--' ,...
     rad2deg(rep_obs.theta), 20*log10(abs(Et_fg)), '-' ,...
     rad2deg(rep_obs.theta), 20*log10(abs(Et_k - Et_fg)), '-.');
xlabel('\theta (deg)');
title('Champ rayonne total')
grid on;
axis(conf_axis);

legend('numerical E_x', 'numerical E_y', 'numerical E_z', ...
       'asymptotic E_x', 'asymptotic E_y', 'asymptotic E_z', ...
       '|\Delta_x|', '|\Delta_y|', '|\Delta_z|');

% figure pour presentation
clf(figure(20))
set(gca, 'FontSize', 20)
plot(rad2deg(rep_obs.theta), 20*log10(abs(Er_k(1,:))), ...
     rad2deg(rep_obs.theta), 20*log10(abs(Epc1(1,:))), ...
     rad2deg(rep_obs.theta), 20*log10(abs(Epc2(1,:))), ...
     rad2deg(rep_obs.theta), 20*log10(abs(Epc3(1,:))));
legend('Integration numerique des courants : E_x', ...
        'Contribution des points critiques 1er ordre : E_x', ...
        'Contribution des points critiques 2nd ordre : E_x', ...
        'Contribution des points critiques 3eme ordre : E_x');
xlabel('\theta (deg)'); ylabel('Amplitude (dB)');
title(sprintf('Comparaison entre integration numerique \n et expression asymptotique'));
grid on;
axis(conf_axis);







if id_pol == 1
  %  Etheta
  e_theta = [cos(rep_obs.theta).*cos(rep_obs.phi); ...
            cos(rep_obs.theta).*sin(rep_obs.phi); ...
            -sin(rep_obs.theta)];
  
  plot_1 = dot(Et_k, e_theta);
  plot_2 = dot(Et_fg,e_theta);
%    plot_3 = dot(Et_fg2,e_theta);
  legend1 = 'PO Numerical integration E_\theta';
  legend2 = 'PO Asymptotic E_\theta (farfield GB)';
%    legend3 = 'PO Asymptotic E_\theta (paraxial GB)';
elseif id_pol == 2
  %  Phi
  e_phi = [-sin(rep_obs.phi);cos(rep_obs.phi);0]*ones(1,length(rep_obs.theta));
  
  plot_1 = dot(Et_k, e_phi);
  plot_2 = dot(Et_fg,e_phi);
%    plot_3 = dot(Et_fg2,e_phi);
  legend1 = 'PO Numerical integration E_\phi';
  legend2 = 'PO Asymptotic E_\phi (farfield GB)';
%    legend3 = 'PO Asymptotic E_\phi (paraxial GB)';
end

%  coefficient de normalisation
if bool_normalisation
  no = 1./max(plot1);
else
  no = 1;
end


% figure pour publication
clf(figure(30))
plot(rad2deg(rep_obs.theta), 20*log10(abs(no.*plot_1)), '--', ...
     rad2deg(rep_obs.theta), 20*log10(abs(no.*plot_2)), ...
...%       rad2deg(rep_obs.theta), 20*log10(abs(no.*plot_3)), ...
     rad2deg(rep_obs.theta), 20*log10(abs(no.*plot_1 - no.*plot_2)), 'g:') ...
%       rad2deg(rep_obs.theta), 20*log10(abs(no.*plot_1 - no.*plot_3)), 'r:'  );
set(gca, 'FontSize', 20);
legend(legend1, legend2, ...%legend3, ...
        'Difference (farfield)', ...%'Difference (paraxial)', ...
       'Location', 'SouthOutside');
xlabel('\theta (deg)');
ylabel('amplitude (dB)');
title(sprintf('Comparison between numerical integration\n and asymptotic expression'));
grid on;
axis(conf_axis);

