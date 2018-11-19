npoints = 1000;
dimension = 5;

% case 1: uniformly sampled sphere
X = randsphere(npoints,dimension,1);

% case 2: iris dataset
% load fisheriris.mat;
% X = zscore(meas);

% case 3: macrophages
%X = importdata('macrophages.txt');
%X = X.data;
%X = X';

%case 4: single cell data
%X = importdata('neuro_ul.txt');
%X = X.data;
%mn = mean(X,2);
%X = X - repmat(mn,1,size(X,2));

[n_alpha,n_single,p_alpha,alphas,separable_fraction] = SeparabilityAnalysis(X);
