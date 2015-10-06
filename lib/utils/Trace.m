% Fonction : Trace	
% Rep      : Focntions/Visualisation
% Date     : 27/01/03
% Auteur   : A .Chabory
% Etat     : OK
%
%  Trace un champ
%
%  Trace(U,Axes)
%
%  - Entrees
%    U         : Valeurs du champ
%    Axes      : Axes du trace
%    varargin  : (Optionnels)
%                 Definition des preferences
%                 'caxis'     : caxis ou dynamique
%                 'type'      : valeurs a tracer dB,abs,nat,phase
%                 'titre'     : {'titre','xlabel','ylabel','zlabel'}
%                 'legende '  : legende




function Trace(U,Axes,varargin)


type='nat';
titre={'','','',''};
cax=[];
legende={};
for ind=1:2:length(varargin)-1
   argument=varargin{ind};
   switch argument
   case 'titre'
      titre=varargin{ind+1};
   case 'caxis'
      cax=varargin{ind+1};
   case 'type'
      type=varargin{ind+1};
   case 'legende'
      legende=varargin{ind+1};
   end
end

%set(0,'DefaultAxesFontSize',18);
%set(0,'DefaultAxesFontWeight','normal');
%set(0,'DefaultAxesFontName','helvetica');


% mise en forme pour que E soit (Naxes1,Naxes2)
if strcmp(class(U),'double')
   U={U};
end
if strcmp(class(Axes),'double')
   Axes={Axes};
end

switch type
case 'dB'
   U{1}=10.*log10(sum( abs(U{1}).^2 ,1) );
case 'abs'
   U{1}=sqrt(sum(abs(U{1}).^2 ,1));
case 'nat'
   U{1}=U{1}(1,:);
case 'phase'
   U{1}=angle(U{1}(1,:)).*180./pi;
end

maxi=max(max(max(U{1})));
mini=min(min(min(U{1})));

for ind=2:length(U)
   % calcul des valeurs a tracer
   switch type
   case 'dB'
      U{ind}=10.*log10(sum( abs(U{ind}).^2 ,1 ) );
   case 'abs'
      U{ind}=sqrt(sum(abs(U{ind}).^2 ,1));
   case 'nat'
      U{ind}=U{ind}(1,:);
   case 'phase'
      U{ind}=angle(U{ind}(1,:)).*180./pi;
   end
   maxi=max([max(max(U{ind})) maxi]);
   mini=min([min(min(U{ind})) mini]);
end

switch length(cax)
case 0
   cax=[mini maxi];
case 1
   cax=[maxi-cax maxi];
end

if length(Axes)==2
   for ind=1:length(U)
      U{ind}=reshape(U{ind},length(Axes{2}),length(Axes{1}));
      figure
      pcolor(Axes{1},Axes{2},U{ind});
      caxis(cax);
      shading interp;
      colorbar;
      axis image;
      title(titre{1});
      xlabel(titre{2}); 
      ylabel(titre{3}); 
   end
elseif length(Axes)==3
   for ind=1:length(U)
      U{ind}=reshape(U{ind},length(Axes{1}),length(Axes{2}));
      figure
      [tetg,phig]=meshgrid(Axes{1},Axes{2});
       pcolor(sin(tetg).*cos(phig),...
         sin(tetg).*sin(phig),U{ind});
      caxis(cax);
      shading interp;
      colorbar;
      axis image;
      title(titre{1});
      xlabel(titre{2}); 
      ylabel(titre{3}); 
   end
else
   if size(Axes{1},1)==1
      style={'b','r','k','b--','r--','k--'};
      figure;
      hold on;
      for ind=1:length(U)
         plot(Axes{1},U{ind},style{ind});
      end
      axis([Axes{1}(1),Axes{1}(end),cax]);
      grid on;
      title(titre{1});
      xlabel(titre{2}); 
      ylabel(titre{3});
      if ~isempty(legende)
         legend(legende);
      end
   end
   if size(Axes{1},1)==3
      Tri=delaunay(Axes{1}(1,:),Axes{1}(2,:));
      for ind=1:length(U)
         figure
         trisurf(Tri,Axes{1}(1,:),Axes{1}(2,:),Axes{1}(3,:),U{ind});
         caxis(cax);
         shading interp;
         colorbar;
         axis image;
         title(titre{1});
         xlabel(titre{2}); 
         ylabel(titre{3});
         zlabel(titre{4});
      end
   end
end



         
   



   

