function textwaitbar(x, textbar)
%
% waitbar en mode texte à la manière d'Octave
%
% ARGUMENTS :
%  x : avancement (entre 0 et 1) (1 x 1)
%  textbar (optionnel) : description du calcul (string)
%

% ==========
% CONSTANTES
% ==========
NB_CHAR = 40;
TITLE = ['##[  ', textbar, '  ]##'];
NB_TITLE = length(TITLE);

% On efface tout
clc;

% On affiche le titre du calcul
if nargin>=2
    if isstr(textbar)
        disp(TITLE);
    end
end

% On affiche la barre de progression
pourcentage = x*100;
S = '[';
for q = 1:NB_CHAR
    if (q < ceil(pourcentage*NB_CHAR/100))
        S = [S, '='];
    else
        S = [S, ' '];
    end
end
S=[S, '] ', num2str(pourcentage), ' %'];
disp(S);
