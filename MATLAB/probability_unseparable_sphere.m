function [p] = probability_unseparable_sphere(alpha,n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%p = power((1-alpha*alpha),(n-1)/2)/(alpha*sqrt(pi*(n-2)));
p = power((1-alpha*alpha),(n-1)/2)/(alpha*sqrt(2*pi*n));
%p = power((1-alpha*alpha),(n-1)/2)/(alpha*power(40*pi*(n-2),1/3));
%p = power((1-alpha*alpha),(n-1)/2)/(alpha*power(200*pi*(n-2),1/4));
%p = power((1-alpha*alpha),(n-1)/2)/(alpha*power(20*pi*(n-1),1/3));

end

