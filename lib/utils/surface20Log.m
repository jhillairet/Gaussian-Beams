% 
function surface20Log(varargin)
%
% raccourci vers la commande :
% surface(20*log10(Z))
if (length(varargin) > 1)
    surface(varargin{1}, varargin{2}, 20*log10(varargin{3}));
else
    surface(20*log10(varargin{1}));
end