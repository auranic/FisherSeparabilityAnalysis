function [Xp] = preprocessing(X,center,dimred,whiten,projectonsphere,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

eigval_retaining_factor = 10;
for i=1:length(varargin)
    if(strcmp(varargin{i},'EigValueRetainingFactor'))
        eigval_retaining_factor = varargin{i+1};
    end
end

Xp = X;

% centering

nobjects = size(X,1);
sampleMean = mean(X,1);
if center
    Xp = bsxfun(@minus, X, sampleMean);
    %Xp = (X - repmat(sampleMean,nobjects,1));
end

% dimension reduction

PCAcomputed = 0;

if dimred|whiten
    [v,u,s] = pca(Xp);
    PCAcomputed = 1;
    sc = s/s(1);
    ind = find(sc>1/eigval_retaining_factor);
    Xp = Xp*v(:,ind);
    disp(sprintf('%i components is retained using factor %2.2f',length(ind),eigval_retaining_factor));
end

% whitening

if whiten
    Xp = u(:,ind);
    st = std(Xp);
    Xp = Xp./repmat(st,nobjects,1);
end

% project on sphere (scale each vector to unit length)
if projectonsphere
    st = sqrt(sum(Xp.^2,2));
    Xp = Xp./repmat(st,1,size(Xp,2));
end

end

