function [separ_fraction,py,unsep_map] = checkSeparability(xy,alpha)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% xy is gram matrix : xy = X*X'

npoints = size(xy,1);

dxy = diag(xy);
sm = bsxfun(@rdivide,xy,dxy);
sm = sm - diag(diag(sm));
sm = sm>alpha;
py = sum(sm');
py = py/length(py);
separ_fraction = sum(py==0)/length(py);

unsep_map = containers.Map('KeyType','int32','ValueType','any');
for i=1:npoints
    unsep_map(i) = find(sm(i,:));
end

end

