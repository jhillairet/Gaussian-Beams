function d=nsqrt(c)
%  % calcul de la racine carree avec determination 
%  % de la partie imaginaire negative
%  % Source : wikipedia
%  %  http://fr.wikipedia.org/wiki/Racine_carr%C3%A9e
%  u = sign(-imag(c)) + (imag(c)==0);
%  d = sqrt((abs(c)+real(c))/2)-j*sqrt((abs(c)-real(c))/2); 

%  Modification de la définition usuelle de la racine carree
%  complexe. Il s'agit en utilisant la definition suivante, de
%  toujours obtenir la partie imaginaire de d negative. 
%  
%  Contrairement à la définition usuelle, la partie réelle
%  peut etre positive ET negative.
d = -i*sqrt(-c);


%  %  Autre definition courant de la racine complexe :
%  theta = atan2(imag(c), real(c));
%  d = exp(log(abs(c))/2).*exp(i*theta/2+i*pi);

%  % definition de la racine pour que la partie imaginaire soit negative
%  d = exp((log(-1)-log(-1./c))/2);



%  http://mathworld.wolfram.com/SquareRoot.html
%  d = sqrt(abs(c)).*(cos((atan2(imag(c),real(c)))/2)+j*sin((atan2(imag(c),real(c)))/2));

%  %  Mathematics: Complex functions
%  %  Marcel Leutenegger
%  %  Ecole Polytechnique Federale de Lausanne
%  %  Laboratoire d Optique Biomedicale
%  %  January 2, 2005
%  d1 = sqrt((abs(c)+real(c))/2) + j.*imag(c)/2.*sqrt(2./(abs(c)+real(c))); ...
%  d2 = abs(imag(c))/2.*sqrt(2./(abs(c)-real(c))) - j*sqrt((abs(c)-real(c))/2);
%   u = real(c) >= 0;
