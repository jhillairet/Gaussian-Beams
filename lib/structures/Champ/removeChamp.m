function ChampRest = removeChamp(Champ, ind)
%  Supprime les elements de la structure Champ pour 
%  les indexes ind
%  
%  ARGUMENTS :
%    Champ : structure champ
%    ind : index(es)
%  
%  RETOURNE :
%    ChampInd : structure champ restante
%    

if ind > size(Champ.Centre,2)
  error(['Erreur : on essaye d''acceder à l''index ', num2str(ind), ...
        ' d''une structure Champ de taille ', num2str(size(Champ.Centre,2))]);
end



Champ.Centre(:,ind) = [];
Champ.ex(:,ind) = [];
Champ.ez(:,ind) = [];
Champ.Coefficients(:,ind) = [];
Champ.Courbure(:,ind) = [];

ChampRest = Champ;

