function [separ_fraction,py] = checkSeparabilityMultipleAlphaBig(data, alpha)
% data is data matrix to compare separability
% alpha is array of alphas to test

    %Number of points per 1 loop
    nP = 10000;

    %Normalise alphas
    if(size(alpha,1))>1
        alpha = alpha';
    end
    addedone = 0;
    if max(alpha)<1
        alpha = [alpha 1];
        addedone = 1;
    end
    
    alpha = [-Inf alpha Inf];

    n = size(data, 1);
    counts = zeros(n, length(alpha));
    leng = zeros(n, 1);
    for k = 1:nP:n
        e = k + nP - 1;
        if e > n
            e = n;
        end
        
        %Calculate diagonal part
        xy = data(k:e, :) * data(k:e, :)';
        leng(k:e) = diag(xy);
        xy = xy - diag(leng(k:e));
        xy = bsxfun(@rdivide, xy, leng(k:e));
        counts(k:e, :) = counts(k:e, :) + histc(xy', alpha)';
                
        for kk = 1:nP:n
            ee = kk + nP - 1;
            if ee > n
                ee = n;
            end
            %Ignore diagonal part
            if k == kk
                continue;
            end
            xy = data(k:e, :) * data(kk:ee, :)';
            xy = bsxfun(@rdivide, xy, leng(k:e));
            counts(k:e, :) = counts(k:e, :) + histc(xy', alpha)';
        end
    end
    
    na=length(alpha); 
    for i=1:na-1 
        counts(:,na-i)=counts(:,na-i)+counts(:,na-i+1); 
    end; 

    py = counts / (n - 1);
    py = py';
     if addedone
         py = py(2:end-2,:);
     else
         py = py(2:end-1,:);
     end
    separ_fraction = sum(py==0) / length(py);
end
