function h = imshow2(varargin)
varargin = [varargin, 'InitialMagnification', 'fit'];

hh = imshow(varargin{:}); 
titlestr = inputname(1);
titlestr = regexprep(titlestr,'_','\\_');
title(titlestr);

if (nargout > 0)
% Only return handle if caller requested it.
h = hh;
end