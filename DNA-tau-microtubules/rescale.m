function varargout = rescale(varargin)
%RESCALE Rescales the range of data
%   R = RESCALE(A)
%   R = RESCALE(A,B,C)
%   R = RESCALE(...,'InputMin',IMIN)
%   R = RESCALE(...,'InputMax',IMAX)
%   
%   Example:
%       % Clip all entries to [2,4] and then rescale all entries to [-1,1]
%       a = distributed([1;2;3;4;5]);
%       r = rescale(a,-1,1,'InputMin',2,'InputMax',4)
%   
%   See also RESCALE, DISTRIBUTED.


%   Copyright 2017-2019 The MathWorks, Inc.

[varargout{1:max(nargout, 1)}] = distributedutil.distributedSpmdWrapper( @rescale, varargin{:} );
