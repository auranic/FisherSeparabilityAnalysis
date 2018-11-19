function [separ_fraction,py] = checkSeparability(xy,alpha)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% xy is gram matrix : xy = X*X'

dxy = diag(xy);
sm = bsxfun(@rdivide,xy,dxy);
sm = sm - diag(diag(sm));
sm = sm>alpha;
py = sum(sm');
py = py/length(py);
separ_fraction = sum(py==0)/length(py);

end

