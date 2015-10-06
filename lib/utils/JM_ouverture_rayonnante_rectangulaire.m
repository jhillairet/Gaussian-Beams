% Calcul des courants equivalents electriques et magnetique sur une ouverture rayonnante 
% d'un guide d'onde rectangulaire de largeur a et de hauteur b.
%
% ENTRÉES :
% 	a,b : largeur et hauteur du guide (m)
%	dx,dy : pas d'echantillonage de la surface rayonnante
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
function [Q,J,M,E,H] = JM_ouverture_rayonnante_rectangulaire(a,b,dx,dy,m,n,bool_TE,f);

% constantes et parametres utiles
    c  = 2.997925e8;
 omega = 2*pi*f;
lambda = c/f;
    k0 = 2*pi/lambda;
   mu0 = 4*pi*1e-7;
 epsr0 = 8.85*1e-12 ;

		kx = m*pi/a;
		ky = n*pi/b;
%  nombre d'onde de coupure (cutoff)
    kc = sqrt(kx^2 + ky^2);
%  nombre d'onde guide
  beta = sqrt(k0^2-kc^2);

% points de la surface rayonnante
     Qx = [0:dx:a];
     Qy = [0:dy:b];
		 Qz = 0;
% on met tous les points de la surface dans un vecteur unique de taille (3x(X*Y*Z))
		 Q = stack(Qx,Qy,Qz);

% pre allocation
		E = zeros(3,size(Q,2));
		H = zeros(3,size(Q,2));

if bool_TE
  E(1,:) = j*omega*mu0*n*pi/kc^2/b.*cos(kx.*Q(1,:)).*sin(ky.*Q(2,:));
  E(2,:) = -j*omega*mu0*m*pi/kc^2/a.*sin(kx.*Q(1,:)).*cos(ky.*Q(2,:));
  H(1,:) = j*beta*m*pi/kc^2/a.*sin(kx.*Q(1,:)).*cos(ky.*Q(2,:));
  H(2,:) = j*beta*n*pi/kc^2/b.*cos(kx.*Q(1,:)).*sin(ky.*Q(2,:));
  H(3,:) = cos(kx.*Q(1,:)).*sin(ky.*Q(2,:));
	% Ez nul en TE
else
	E(1,:) = j*k0*kx./(kx.^2+ky.^2).*cos(kx.*Q(1,:)).*sin(ky.*Q(2,:));
	E(2,:) = j*k0*ky./(kx.^2+ky.^2).*sin(kx.*Q(1,:)).*cos(ky.*Q(2,:));
	E(3,:) = sin(kx.*Q(1,:)).*sin(ky.*Q(2,:));	
	H(1,:) = + omega/c/k0.*E(2,:);
	H(2,:) = - omega/c/k0.*E(1,:);
	% Hz nul en TM	
end

%  % normalisation des champs à un champ electrique unite :
%  E = E./(ones(3,1)*sqrt(sum(E.^2)));
%  H = H./(ones(3,1)*sqrt(sum(E.^2)));
% calcul des courants equivalents sur la surface :
	n = [0;0;1]*ones(1,size(Q,2));	
	J = cross(n,H);
	M = cross(E,n);