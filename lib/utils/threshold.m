function A = threshold(A, n)
%  A = threshold(A, n)
%
% Seuillage du vecteur A � la valeur min n. 
% Tout ce qui est en dessous de n est d�clar� nul (0)
% cad : si A < n, A = 0.
%
% PARAMETRES :
%  A : vecteur  (1 x NbA)
%  n : seuil (threshold)  (1 x 1) ou (1 x NbA)
%
% RETOURNE :
%  A : vecteur A seuill�
%

A = A .*  (abs(A) > n);
