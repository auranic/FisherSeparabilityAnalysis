function p = probability_unseparable_sphere(alpha,n)
%probability_unseparable_sphere calculate theoretical probability for point
%to be inseparable for dimension n
%
%Inputs:
%   alpha is 1-by-d vector of possible alphas. Must be row vector or scalar
%   n is c-by-1 vector of dimnesions. Must be column vector or scalar.
%
%Outputs:
%   p is c-by-d matrix of probabilities.

    p = bsxfun(@power, 1-alpha .^ 2, (n - 1) / 2)...
        ./ bsxfun(@times, alpha, sqrt(2 * pi * n));
end

