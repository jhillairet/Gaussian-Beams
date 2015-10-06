% Calcul des courants equivalents electriques et magnetique sur une ouverture rayonnante 
% d'un guide d'onde rectangulaire de largeur a et de hauteur b.
%
% ENTRÉES :
% 	a : rayon du guide (m)
%	drho,dphi : pas d'echantillonage de la surface rayonnante
% 	m,n : indice du mode 
% 	bool_TE : 1 pour mode TE, 0 pour mode TM (pas de mode TEM)
% 	f : frequence
% 
% SORTIES : 
% 	Q : point d'echantillonnage des courants (3,Nb)
%      Qx = [0:dx:a]; (1,Nx)
%      Qy = [0:dy:b]; (1,Ny)
%  		 Qz = 0;        (1,1)
%  		 Q = stack(Qx,Qy,Qz); (3,Nx*Ny)=(3,Nb)
% 	J : vecteur densite de courant électrique (3,Nb)
% 	M : vecteur densite de courant magnetique (3,Nb)
%		E : champ électrique sur la surface
% 	H : champ magnetique sur la surface
function [Q,J,M,E,H] = JM_ouverture_rayonnante_circulaire(a,drho,dphi,m,n,bool_TE,f);

% constantes et parametres utiles
    c  = 2.997925e8;
 omega = 2*pi*f;
lambda = c/f;
    k0 = 2*pi/lambda;
   mu0 = 4*pi*1e-7;
 epsr0 = 8.85*1e-12 ;

% points de la surface rayonnante
    Qrho = [0:drho:a];
    Qphi = [0:dphi:2*pi];
      Qz = 0;

% on met tous les points de la surface dans un vecteur unique de taille (3x(X*Y*Z))
		 Q_cyl = stack(Qrho,Qphi,Qz);
     Q = [Q_cyl(1,:).*cos(Q_cyl(2,:)); Q_cyl(1,:).*sin(Q_cyl(2,:)); Q_cyl(3,:)];



%  m-ieme zeros de la derivee de la fonction de Bessel J_n
pmn = [3.832 7.016 10.174; ...% n = 0
       1.841 5.331 8.536; ... % n = 1
       3.054 6.706 9.970]; % n= 2

%  nombre d'onde de coupure (cutoff)
kc = pmn(m+1,n+1)/a;

%  beta
beta = sqrt(k0^2 - kc^2);

% pre allocation
% E et H sont des vecteurs dont les composantes
% sont exprimees dans le repere cylindrique : [Erho;Ephi;Ez]
		E = zeros(3,size(Q,2));
		H = zeros(3,size(Q,2));


dbesselJ = @(x,n) -besselj(n + 1, x) + n./x.*besselj(n, x);



if bool_TE
	H(1,:) = -j*beta/kc.*sin(n.*Q_cyl(2,:)).*dbesselJ(kc.*Q_cyl(1,:),n);
	H(2,:) = -j*beta.*n/kc^2./Q_cyl(1,:).*cos(n.*Q_cyl(2,:)).*besselj(n, kc.*Q_cyl(1,:));
	H(3,:) = 	sin(n.*Q_cyl(2,:)).*besselj(n, kc.*Q_cyl(1,:));
	E(1,:) = - j*omega*mu0*n/kc^2./Q_cyl(1,:).*cos(n.*Q_cyl(2,:)).*besselj(n, kc.*Q_cyl(1,:));
	E(2,:) = j*omega*mu0/kc.*sin(n.*Q_cyl(2,:)).*dbesselJ(kc.*Q_cyl(1,:),n);
	% Ez nul en TE
else

	% Hz nul en TM	
end
%  passage en coordonnees cartesiennes
Ex = E(1,:).*cos(Q(2,:)) - E(2,:).*sin(Q(2,:));
Ey = E(1,:).*sin(Q(2,:)) + E(2,:).*cos(Q(2,:));
E = [Ex;Ey;E(3,:)];

Hx = H(1,:).*cos(Q(2,:)) - H(2,:).*sin(Q(2,:));
Hy = H(1,:).*sin(Q(2,:)) + H(2,:).*cos(Q(2,:));
H = [Hx;Hy;H(3,:)];
%  % normalisation des champs à un champ electrique unite :
%  E = E./(ones(3,1)*sqrt(sum(E.^2)));
%  H = H./(ones(3,1)*sqrt(sum(E.^2)));
% calcul des courants equivalents sur la surface :
	n = [0;0;1]*ones(1,size(Q,2));	
	J = cross(n,H);
	M = cross(E,n);