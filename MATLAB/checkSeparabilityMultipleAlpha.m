function [separ_fraction, py] = checkSeparabilityMultipleAlpha(data, alpha)
%checkSeparabilityMultipleAlpha calculate fraction of points inseparable
%for each alpha and fraction of points which are inseparable from each
%point for different alpha.
%
%Inputs:
%   data is data matrix to calculate separability. Each row contains one
%       data point.
%   alpha is array of alphas to test separability.
%
%Outputs:
%   separ_fraction fraction of points inseparable from at least one point.
%       Fraction is calculated for each alpha.
%   py is n-by-m matrix. py(i,j) is fraction of points which are
%       inseparable from point data(i, :) for alphas(j).

    %Number of points per 1 loop. 20k assumes approx 3.2GB
    nP = 20000;

    %Normalise alphas
    if(size(alpha, 1)) > 1
        alpha = alpha(:)';
    end
    
    addedone = 0;
    if max(alpha)<1
        alpha = [alpha 1];
        addedone = 1;
    end
    
    alpha = [-Inf, alpha, Inf];

    n = size(data, 1);
    counts = zeros(n, length(alpha));
    leng = zeros(n, 1);
    for k = 1:nP:n
        e = k + nP - 1;
        if e > n
            e = n;
        end
        
        % Calculate diagonal part, divide each row by diagonal element
        xy = data(k:e, :) * data(k:e, :)';
        leng(k:e) = diag(xy);
        xy = xy - diag(leng(k:e));
        xy = bsxfun(@rdivide, xy, leng(k:e));
        counts(k:e, :) = counts(k:e, :) + histc(xy', alpha)';
        % Calculate nondiagonal part
        for kk = 1:nP:n
            %Ignore diagonal part
            if k == kk
                continue;
            end
            ee = kk + nP - 1;
            if ee > n
                ee = n;
            end
            xy = data(k:e, :) * data(kk:ee, :)';
            xy = bsxfun(@rdivide, xy, leng(k:e));
            counts(k:e, :) = counts(k:e, :) + histc(xy', alpha)';
        end
    end
    
    % Calculate comulative sum
    counts = cumsum(counts, 2, 'reverse');
    % 
    py = counts / (n - 1);
    py = py';
     if addedone
         py = py(2:end - 2, :);
     else
         py = py(2:end - 1, :);
     end
    separ_fraction = sum(py == 0) / length(py);
end
