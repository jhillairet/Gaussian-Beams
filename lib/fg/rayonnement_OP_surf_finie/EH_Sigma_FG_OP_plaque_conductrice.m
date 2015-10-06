% Fonction :   EH_Sigma_FG_OP_plaque_conductrice                          
% Date     :   Juin 2007            
% Auteur   :   J.Hillairet                 
%                                                          
%                                                                                        
%                                                          
function [E,H] = EH_Sigma_FG_OP_plaque_conductrice(P, ChampFG, normale_plaque, Plaque, Freq, Epsr, W0, varargin)

%  Doit on effectuer un lancer de FG pour eliminer
%  les faisceaux ne rencontrent pas la plaque ?   
%  EXPERIMENTAL
bool_lancerFG = true;

%  doit on inverser le signe du champ rayonne 
%  par les faisceaux dont le centre est 
%  situe sous le sol (y<0) ?
%  EXPERIMENTAL
bool_signeimage = false;

if bool_lancerFG
  %  ====================================================================
  %  filtrage des FG
  %  ====================================================================
  %  On recupere les index des faisceaux qui vont contribuer
  %  au champ rayonne. (ceux qui touchent le plan ou son voisinage) 
  %  Les autres faisceaux ne sont pas pris
  %  en compte
  index = [];
  for ind=1:size(ChampFG.Centre,2)
    % faisceau actuel
    FG_ind=getChamp(ChampFG,ind);
    % distance parcourue pour atteindre le plan de la plaque
    d =  (Plaque.Centre(3) - FG_ind.Centre(3))./FG_ind.ez(3);
    %  coordonnees du point d'intersection :
    O_l = FG_ind.Centre + d*FG_ind.ez;
    %  est ce que le point d'intersection 
    %  appartient au domaine de la plaque (+ 2*largeurFG) ?
    delta = 2*sqrt(W0^2.*(1+d^2));
  
    test =  O_l(1) > min(Plaque.Sommets(1,:))-delta ...
          & O_l(1) < max(Plaque.Sommets(1,:))+delta ...
          & O_l(2) > min(Plaque.Sommets(2,:))-delta ...
          & O_l(2) < max(Plaque.Sommets(2,:))+delta;
    if test 
      index=[index, ind];
      figure(1)
      hold on
        plot3(-FG_ind.Centre(1), FG_ind.Centre(3), FG_ind.Centre(2), 'g.');
      hold off
    else
%        disp('Pas d''intersection : faisceau elimine');
      figure(1)
      hold on
        plot3(-FG_ind.Centre(1), FG_ind.Centre(3), FG_ind.Centre(2), 'r.');
      hold off
    end
  end
  
  disp([num2str(size(ChampFG.Centre,2) - length(index)), 'Faisceaux elimines']);
else
  disp(['Prise en compte de TOUS les faisceaux']);
  index=1:size(ChampFG.Centre,2);
end


%  ====================================================================
%  Changement de repere
%  la normale au plan rectangulaire devient le nouvel axe z
%  et le barycentre du rectangle la nouvelle origine absolue
%  ====================================================================
  %  vecteur normal a la plaque
  %  Il s'agit de la normale aux diagonales du rectangle
  N = cross(Plaque.Sommets(:,4)-Plaque.Sommets(:,2), Plaque.Sommets(:,3)-Plaque.Sommets(:,1));
  N = N/norm(N)  


  %  On teste si tous les Sommets appartiennent bien au meme plan
  %  Si ce n'est pas le cas, c'est une erreur !
  % TODO
  
  %  Centre (barycentre) du rectangle :
  %  il s'agit de la nouvelle origine du repere absolu 
  CentreGrav =  sum(Plaque.Sommets,2)/4;
  
  
  %  deux vecteurs appartenant au plan du rectangle
  e0 = [1;0;0];
  e1 = Plaque.Sommets(:,2)-Plaque.Sommets(:,1);
  e1 = e1/norm(e1);
  e2 = Plaque.Sommets(:,1)-Plaque.Sommets(:,4);
  e2 = e2/norm(e2);
  
  MRot = [e1 e2 N];
  P_plaque = MRot.'*(P - CentreGrav*ones(1,size(P,2)));

  Plaque.Centre = [0;0;0];
  Plaque.N = N;


  Plaque.Sommets = MRot.'*(Plaque.Sommets - CentreGrav*ones(1,size(Plaque.Sommets,2))) ;
  if ~isempty(varargin)
    NewCentre = MRot.'*(varargin{1} - CentreGrav) ;
  else
    NewCentre = [];
  end
%    MRot = eye(3);
%    P_plaque = P;







%Progress Bar

PB = progress('init', ['Sommation des ', num2str(length(index)),  ' faisceaux gaussiens']);

warning off;
E = zeros(3,size(P,2));
H = zeros(3,size(P,2));
ind2=0;
for ind=index  
   % faisceau actuel
    FG_ind = getChamp(ChampFG,ind);

    FG_ind.Centre = MRot.'*(FG_ind.Centre - CentreGrav);
    FG_ind.ex = MRot.'*FG_ind.ex;
    FG_ind.ez = MRot.'*FG_ind.ez;
    FG_ind.Coefficients = [FG_ind.Coefficients(1)*dot([1;0;0], MRot.'*[1;0;0]); ...
                           FG_ind.Coefficients(2)*dot([0;1;0], MRot.'*[0;1;0])]   ;

    [Et,Ht, Epc1, Epc2, Epc3] = EH_FG_OP_plaque_conductrice_uniforme_v2(FG_ind, Plaque, P_plaque, Freq, Epsr, NewCentre);
%      [Et,Ht, Epc1, Epc2, Epc3] = EH_FG_lointain_OP_plaque_conductrice_uniforme(FG_ind, Plaque, P_plaque, Freq, Epsr, NewCentre);
  if bool_signeimage
    TMP = getChamp(ChampFG,ind);
    if TMP.Centre(2) < 0
      Et = diag([-1;1;-1])*Et;
      Ht = diag([-1;1;-1])*Ht;
    end
  end


    E = E + MRot*Et;
    H = H + MRot*Ht;

  % MAJ Progress Bar
  ind2=ind2+1;
  PB = progress(PB, ind2/length(index));
end

warning on;