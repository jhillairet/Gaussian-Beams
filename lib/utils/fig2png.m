%  transforme des images .fig 
%  dans le repertoire rep
%  en image bitmap au format png 
%  
function status = fig2png(rep)
if ~exist(rep, 'dir')
  status = 'repertoire non existant';
  return;
end

%  on se place dans le repertoire rep
actual_dir = pwd;
cd(rep);

%  on liste tous les fichiers .fig
 list_fig = dir('*.fig');

for idx=1:size(list_fig,1)
  open(list_fig(idx).name);
  [path,name,ext,ver]=fileparts(list_fig(idx).name); 
  % TEMPORAIRE   
  %    h=colorbar, set(get(h, 'Title'), 'String','dB','FontSize', 14), set(h, 'FontSize', 14);   
    
  saveas(gcf, [name,'.png']);
%    print(gcf, '-dpng', '-r300', [name,'.png']) 
  close(gcf);
end

cd(actual_dir);

