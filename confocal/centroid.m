function [X, Y] = centroid(pshape, I)
%MATLAB Code Generation Library Function
% CENTROID Get centroid of the polyshape or of single boundary
% inside the polyshape

%   Copyright 2022 The MathWorks, Inc.

%#codegen

coder.internal.polyshape.checkConsistency(pshape, nargin);

if (nargin == 1)
    n = coder.internal.polyshape.checkArray(pshape);
    np = numel(pshape);
    pm = zeros(np, 2);
    for i=1:np
        if pshape(i).isEmptyShape()
            pm(i, :) = [NaN NaN];
        else
            pm(i, :) = pshape(i).polyImpl.getCentroid();
        end
    end
    X = reshape(pm(:, 1), n);
    Y = reshape(pm(:, 2), n);
else
    coder.internal.polyshape.checkScalar(pshape);
    coder.internal.polyshape.checkEmpty(pshape);
    II = coder.internal.polyshape.checkIndex(pshape, I);
    ni = length(II);
    pm = zeros(ni, 2);
    for i=1:ni
        pm(i, :) = pshape.polyImpl.getBoundaryCentroid(II(i));
    end
    X = pm(:, 1);
    Y = pm(:, 2);
end