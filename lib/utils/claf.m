function claf()
%
% Clear All Figures
%
% Author : JH
% 
%

if ~isOctave % Instructions pour MATLAB
% get all figures handles
  handles = get(0, 'Children');

% clear all figures
  for ind=1:length(handles)
    clf(handles(ind));
  end
else % Instructions pour OCTAVE
  close all;
end
