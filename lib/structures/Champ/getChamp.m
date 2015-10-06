function ChampInd = getChamp(Champ, ind)
%  Recupere les elements de la structure Champ pour 
%  les indexes ind
%  
%  ARGUMENTS :
%    Champ : structure champ
%    ind : index(es)
%  
%  RETOURNE :
%    ChampInd : structure champ aux index 
%    
ChampInd = struct_Champ(...
              Champ.Centre(:,ind), ...
              Champ.ex(:,ind), ...
              Champ.ez(:,ind), ...
              Champ.Coefficients(:,ind), ...
              Champ.Courbure(:,ind));
