function [n, n_single_estimate, alfa_single_estimate]...
    = dimension_uniform_sphere(py, alphas)
%Gives an estimation of the dimension of uniformly sampled n-sphere
%corresponding to the average probability of being unseparable and a margin
%value 
%
%Inputs:
%   py - average fraction of data points which are INseparable.
%   alphas - set of values (margins), must be in the range (0;1)
% It is assumed that the length of py and alpha vectors must be of the
% same.
%
%Outputs:
%   n - effective dimension profile as a function of alpha
%   n_single_estimate - a single estimate for the effective dimension 
%   alfa_single_estimate is alpha for n_single_estimate.
%

    % Sanity check of arguments
    if length(py)~=length(alphas)
        error('ERROR: length of py (%i) and alpha (%i) does not match',...
            length(py), length(alphas));
    end
    if ~isnumeric(alphas) || sum(alphas <= 0) > 0 || sum(alphas >= 1) > 0
        error(['"Alphas" must be a real vector, with alpha range,'...
            ' the values must be within (0,1) interval']);
    end
    

    % Calculate dimension for each alpha
    n = zeros(1,length(alphas));
    for i=1:length(alphas)
        if py(i) == 0
            %All points are separable. Nothing to do and is not interesting
            n(i) = NaN;
        else
            p = py(i);
            a2 = alphas(i)^2;
            w = log(1 - a2);
            n(i) = lambertw(-(w / (2 * pi * p * p * a2 *(1 - a2)))) / (-w);
        end
    end
    % Find indices of alphas which ar not completely separable 
    inds = find(~isnan(n));
    % Find the maximal value of such alpha
    alpha_max = max(alphas(inds));
    % The reference alpha is the closest to 90 of maximal partially
    % separable alpha
    alpha_ref = alpha_max * 0.9;
    [~, k] = min(abs(alphas-alpha_ref));
    % Get corresponding values
    alfa_single_estimate = alphas(inds(k));
    n_single_estimate = n(inds(k));
end