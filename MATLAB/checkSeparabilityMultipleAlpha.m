function [separ_fraction,py] = checkSeparabilityMultipleAlpha(xy,alpha)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% xy is gram matrix : xy = X*X'

dxy = diag(xy);
sm = bsxfun(@rdivide,xy,dxy);


  if(size(alpha,1))>1 
      alpha = alpha';
  end
  addedone = 0;
  if max(alpha)<1
      alpha = [alpha 1];
      addedone = 1;
  end
  
  alpha = [-Inf alpha Inf];
  
%  
%  if min(alpha)>-1
%      alpha = [-1 alpha];
%  end
%  if max(alpha)<1
%      alpha = [alpha 1];
%  end
%  alpha = [-Inf alpha Inf];

sm = sm - diag(diag(sm));
counts = histc(sm',alpha)';
na=length(alpha); 
for i=1:na-1 
    counts(:,na-i)=counts(:,na-i)+counts(:,na-i+1); 
end; 

py = counts/size(sm,2);
py = py';
 if addedone
     py = py(2:end-2,:);
 else
     py = py(2:end-1,:);
 end
separ_fraction = sum(py==0)/length(py);

end

