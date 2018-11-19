function [n,n_single_estimate,alfa_single_estimate] = dimension_uniform_sphere(py,alpha,number_of_data_points)
% Gives an estimation of the dimension of uniformly sampled n-sphere
% corresponding to the average probability of being unseparable 
% and a margin value
%   Arguments:
%   py - average probabilities of a data point to be UNseparable from the
%   rest of data points
%   alpha - set of values (margins), must be in the range (0;1)
%   number_of_data_points (optional) - number of data points used to
%   estimate py
% It is assumed that the length of py and alpha vectors must be of the
% same length
% 

if length(py)~=length(alpha)
    error(sprintf('ERROR: length of py (%i) and alpha (%i) does not match',length(py),length(alpha)));
    return;
end

n = zeros(1,length(alpha));
for i=1:length(alpha)
    if alpha(i)>=1
        error('ERROR: alpha is >=1');
    end 
    if alpha(i)<=0    
        error('ERROR: alpha is <=0');
    end
    if py(i)==0
        %py(i) = 1/number_of_data_points;
        n(i) = NaN;
    else
        p = py(i);
        a = alpha(i);
        a2 = alpha(i)^2;
        n(i) = lambertw(-(log(1-a2)/(2*pi*p*p*a*a*(1-a2))))/(-log(1-a2));
        %napprox(i) = log(-log(1-a2))/(-log(1-a2))-log(2*pi*p*p*a*a*(1-a2))/(-log(1-a2));
        %disp(sprintf('Argument for W(x) %f n(i)=%2.2f approx n(i)=%2.2f',-(log(1-a2)/(2*pi*p*p*a*a*(1-a2))),n(i),napprox(i)));
    end
end

inds = find(~isnan(n));
%k = floor(length(inds)/2+0.5);
alpha_max = max(alpha(inds));
alpha_ref = alpha_max*0.9;
k = find(abs(alpha(inds)-alpha_ref)==min(abs(alpha-alpha_ref)));
alfa_single_estimate = alpha(inds(k));
        p = py(inds(k));
        a = alfa_single_estimate;
        a2 = alfa_single_estimate^2;        
n_single_estimate = lambertw(-(log(1-a2)/(2*pi*p*p*a*a*(1-a2))))/(-log(1-a2));
