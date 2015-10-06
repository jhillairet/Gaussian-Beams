function addpath_recursive(dir_name)
% ADDPATH_RECURSIVE(dir_name)
%
% Ajoute un répertoire ainsi que tous ses sous-répertoire 
% au PATH de matlab (répertoire absolu)
%
% utilisation :
%  addpath_recursive('nom_du_repertoire')
%
% arguments
%   dir_name : nom du répertoire (relatif) à ajouter au path
%
if (isdir(dir_name))
    abs_dir = [pwd, '/', dir_name];
    if (isdir(abs_dir))
        addpath(genpath(abs_dir));
    else
        error(['"',abs_dir,'"', ' n''est pas un répertoire valide']);
    end
else
    error(['"', dir_name, '"', ' n''est pas un répertoire valide']);
end 
