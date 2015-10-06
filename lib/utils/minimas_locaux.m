% P = minimas_locaux(V)
%
% renvoie la position P des minimas locaux d'un vecteur V
%
function P = minimas_locaux(V)
D = diff(V); 
P=[];
for k=1:length(D)-1
    if ((D(k) < 0) & (D(k+1)>0))
        P(k) = k;
    end
end
% on enleve les points nuls (qui ne satisfaisait pas aux "if")
P = nonzeros(P);