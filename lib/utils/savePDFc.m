function [status, pdfname] = savePDFc(h, filename)
% SAVEEPSC
%
% enregistre une figure au format PDF couleur à partir d'une transformation
% d'un fichier EPS couleur. Pour cela, on appelle un programme exterieur,
% epstopdf, qui convertit correctement les EPS couleur de matlab em PDF.
%
% Parametres :
%   h       :   figure handler
%   filename:   nom du fichier 
%
% Il n'est pas utile de prêciser l'extension du fichier, le programme la
% rajoute automatiquement. Si elle est precisée tout de meme, le programme
% ne la rajoute pas.
%
% Retourne ; 
%   status  : 0 si tout c'est  bien passé, autre sinon.
%   pdfname : nom du fichier PDF créé

%%%%%%%%%%%%%%%%%%%%%%%%%%%% test des parametres
if ~ishandle(h)
    error('Le parametre h n est pas un handler de figure !');
end

if ~isstr(filename)
    error('Le parametre filename n est pas une chaine de caractere !');
end

status = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% supression de l'extension si elle est fournise
i = strfind(filename, '.');
% si l'index i est un entier, cela signifie que l'extension a été donnée.
if i~=0
    % on supprime les caracteres a partir de l'index i
    filename = filename(1:i-1);
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% creation d'un fichier EPS couleur
saveas(h, strcat(filename, '.eps'), 'epsc');
%%%%%%%%%%%%%%%%%%%%%%%%%%%% conversion EPS -> PDF
if exist(strcat(filename, '.eps'))
    [status, w] = system(sprintf('epstopdf %s', strcat(filename, '.eps')));

    if status ~= 0
        error('Erreur dans l appel de la fonction epstopdf. Est elle presente sur le systeme ?');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% suppression du fichier temporaire EPS si tout
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% c'est bien passé
    if exist(strcat(filename, '.pdf'), 'file')
        delete(strcat(filename, '.eps'));
        pdfname = strcat(filename, '.pdf');
    else
        error('Aucun fichier pdf n a été créé !');
    end
else
    error('Aucun fichier .eps créé !');
end