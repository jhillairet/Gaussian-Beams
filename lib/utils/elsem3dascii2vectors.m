%  Adapte le format ASCII provenant d'ELSEM3D 
%  en donnees vectorielles
%  
%  PARAMETRES : 
%   - ascii_filename : nom du fichier au format ASCII / PTDO
%   - Q : (3xNbQ) points relatifs aux champs
%  RETOURNE :
%   - E,H (3xNbQ)
function [E,H] = elsem3dascii2vectors(ascii_filename)
  % 1ere etape : on determine le numero de la ligne ou l'on a :
  %   "Param1, Param2, Param3, Param4, Tableau"
  fid = fopen(ascii_filename, 'r');
  if fid==-1
    error('File not found or permission denied');
  end
  
  line_str = [];
  line_Nb = 0;
  while (isempty(findstr(line_str, 'Param1, Param2')) | feof(fid) == 1)
    line_str = fgetl(fid);
    line_Nb = line_Nb + 1;
  end
    
  % puis on lit le fichier ASCII a partir de la ligne suivante  
%    C = textscan(fid, '%f %d8 %d8 %s %n %n', 'headerlines', line_Nb+1);
  C = textscan(fid, '%f %d8 %d8 %s %n %n');

  fclose(fid);
  % tous les champs dans l'ordre :
  % A=[Ex1;Ey1;Ez1;Hx1;Hy1;Hz1;Ex2;Ey2...]  
  A = C{5} + j*C{6};
  E = [A(1:6:end).'; A(2:6:end).'; A(3:6:end).'];
  H = [A(4:6:end).'; A(5:6:end).'; A(6:6:end).'];
