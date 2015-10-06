%  Modifie les attribue des colorbar 
%  d'un ensemble de fichiers
%  
function status = modif_colorbar()
%  on liste tous les fichiers de figure contenant "geometrie"
 list_fig = dir('*geometrie*.fig');

for idx=1:size(list_fig,1)
  open(list_fig(idx).name);
  [path,name,ext,ver]=fileparts(list_fig(idx).name); 
  % TEMPORAIRE   
%    h=colorbar, set(get(h, 'Title'), 'String','dB','FontSize', 14), set(h, 'FontSize', 14);   
  h=colorbar; set(h,'FontSize', 12), set(get(h, 'YLabel'), 'String', 'Amplitude du champ en dB', 'FontSize', 16)   
  saveas(gcf, [name,'.png']);
%    print(gcf, '-dpng', '-r300', [name,'.png']) 
  close(gcf);
end

cd(actual_dir);

