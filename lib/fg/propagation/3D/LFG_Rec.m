% Fonction   :    LFG_Rec
% Date       :    22/08/03
% Auteur     :    A .Chabory
% Etat       :    Ok
%
%
% Lancer de FG recursif
%                
% [Champ_Ref,Champ_Tr,Champ]=LFG_Rec(Champ_Inc,Freq,Objet)
%
% - Entrees
%    Champ Inc :  Champ incident (struct Champ)
%    Freq      :  Frequence      (1)
%    Objet     :  Objet          (struct Objet)
%
% - Sorties          
%    Champ_Ref :  Champ reflechi (struct Champ) 
%    Champ_Tr  :  Champ transmis (struct Champ)
%    Champ     :  Champ total    (liste de struct Champ)

function [Champ_Ref,Champ_Tr,Champ]=LFG_Rec(Champ_Inc,Freq,Objet,Interface,Sens,Champ_Ref,Champ_Tr,Champ,Puiss_Min)

lam0=2.997925e8./Freq;
k0=2*pi./lam0;
Z0=120*pi;

switch nargin
case 3
   Interface=1;
   Sens=1;
   Champ_Ref=Init_Champ;
   Champ_Tr=Init_Champ;
   for ind=1:size(Objet.Epsr,2)
      Champ{ind}{1}=Init_Champ;
      Champ{ind}{2}=Init_Champ;
   end
   Puiss_Inc=  sum(pi./(2*k0*Z0).*sum(abs(Champ_Inc.Coefficients).^2)./...
      sqrt(imag(Champ_Inc.Courbure(1,:)).*imag(Champ_Inc.Courbure(2,:))-imag(Champ_Inc.Courbure(3,:)).^2));
   Puiss_Min=1e-10.*Puiss_Inc;
case 5
   Champ_Ref=Init_Champ;
   Champ_Tr=Init_Champ;
   for ind=1:size(Objet.Epsr,2)
      Champ{ind}{1}=Init_Champ;
      Champ{ind}{2}=Init_Champ;
   end
end



if isempty(Champ_Inc.Coefficients)
   return
end


[Tmp_Ref,Tmp_Tr,Tmp_NT]=LFG_Inter(Champ_Inc,Freq,Objet,Interface,Sens);

if ~isempty(Tmp_Ref.Coefficients)
    Puiss_Tmp_Ref=  pi./(2*k0*Z0).*sum(abs(Tmp_Ref.Coefficients).^2)./...
        sqrt(imag(Tmp_Ref.Courbure(1,:)).*imag(Tmp_Ref.Courbure(2,:))-imag(Tmp_Ref.Courbure(3,:)).^2);
    Mq_Ref=find(Puiss_Tmp_Ref>Puiss_Min);
    Tmp_Ref=Init_Champ(Tmp_Ref,Mq_Ref);
end
if ~isempty(Tmp_Tr.Coefficients)
    Puiss_Tmp_Tr=  pi./(2*k0*Z0).*sum(abs(Tmp_Tr.Coefficients).^2)./...
        sqrt(imag(Tmp_Tr.Courbure(1,:)).*imag(Tmp_Tr.Courbure(2,:))-imag(Tmp_Tr.Courbure(3,:)).^2);
    Mq_Tr=find(Puiss_Tmp_Tr>Puiss_Min);
    Tmp_Tr=Init_Champ(Tmp_Tr,Mq_Tr);
end

%  Rajout JH car on partait trop loin 
%  dans la recursion
%  
%  Utilite ???
if Interface>=length(Objet.Epsr)
 return
end

Milieu=Interface+(1-Sens)./2;

Champ{Milieu}{(3+Sens)./2}=Concat_Champ(Champ{Milieu}{(3+Sens)./2},Tmp_Ref);
Champ{Milieu+Sens}{(3-Sens)./2}=Concat_Champ(Champ{Milieu+Sens}{(3-Sens)./2},Tmp_Tr);

if Milieu==1
   Champ_Ref=Concat_Champ(Champ_Ref,Tmp_Ref);
   [Champ_Ref,Champ_Tr,Champ]=LFG_Rec(Tmp_Tr,Freq,Objet,Interface+Sens,Sens,Champ_Ref,Champ_Tr,Champ,Puiss_Min); 
elseif (Milieu+Sens)==1
   Champ_Ref=Concat_Champ(Champ_Ref,Tmp_Tr);
   [Champ_Ref,Champ_Tr,Champ]=LFG_Rec(Tmp_Ref,Freq,Objet,Interface-Sens,-Sens,Champ_Ref,Champ_Tr,Champ,Puiss_Min); 
elseif (Milieu+Sens)==length(Objet.Epsr)
   Champ_Tr=Concat_Champ(Champ_Tr,Tmp_Tr);
   [Champ_Ref,Champ_Tr,Champ]=LFG_Rec(Tmp_Ref,Freq,Objet,Interface-Sens,-Sens,Champ_Ref,Champ_Tr,Champ,Puiss_Min);  
else
   [Champ_Ref,Champ_Tr,Champ]=LFG_Rec(Tmp_Tr,Freq,Objet,Interface+Sens,Sens,Champ_Ref,Champ_Tr,Champ,Puiss_Min); 
   [Champ_Ref,Champ_Tr,Champ]=LFG_Rec(Tmp_Ref,Freq,Objet,Interface-Sens,-Sens,Champ_Ref,Champ_Tr,Champ,Puiss_Min);  
end
