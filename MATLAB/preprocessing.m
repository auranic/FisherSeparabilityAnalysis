function X = preprocessing(X,center,dimred,whiten,projectonsphere,varargin)
%preprocessing form preprocessed dataset
%
%Inputs
%   X is n-by-d data matrix with n d-dimensional datapoints.
%   center is boolean. True means subtraction of mean vector.
%   dimred is boolean. True means applying of dimensionality reduction with
%       PCA. Number of used PCs is defined by ConditionalNumber argument.
%   whiten is boolean. True means applying of whitenning. True whiten
%       automatically caused true dimred.
%   projectonsphere is boolean. True means projecting data onto unit sphere
%   varargin contains Name Value pairs. One possible value can be:
%       'ConditionalNumber' - a positive real value used to select the top
%           princinpal components. We consider only PCs with eigen values
%           which are not less than the maximal eigenvalue divided by
%           ConditionalNumber Default value is 10. 
%
%Outputs:
%   X is preprocessed data matrix.

    % Check the optional parameter
    ConditionalNumber = 10;
    for i=1:length(varargin)
        if(strcmp(varargin{i},'ConditionalNumber'))
            ConditionalNumber = varargin{i+1};
        end
    end

    % Centering
    if center
        X = bsxfun(@minus, X, mean(X));
    end

    % dimensionality reduction if requested dimensionality reduction or
    % whitenning
    if dimred || whiten
        [v,u,s] = pca(X);
        sc = s/s(1);
        ind = find(sc>1/ConditionalNumber);
        X = X * v(:,ind);
        %disp(sprintf('%i components is retained using factor %2.2f',length(ind),ConditionalNumber));
    end

    % whitening if requested
    if whiten
        X = u(:,ind);
        X = bsxfun(@rdivide, X, std(X));
    end

    % project on sphere (scale each vector to unit length)
    if projectonsphere
        X = bsxfun(@rdivide, X, sqrt(sum(X.^2,2)));
    end
end

