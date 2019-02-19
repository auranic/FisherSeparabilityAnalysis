function [n_alpha,n_single,p_alpha,alphas,separable_fraction] = SeparabilityAnalysis(X,varargin)
% Performs standard analysis of separability and
%   produces standard plots. 
% output - 
% 
% n_alpha - effective dimension profile as a function of alpha
% n_single - a single estimate for the effective dimension 
% p_alpha - p_alpha distributions as a function of alpha,
%           matrix with columns corresponding to the alpha values,
%           and with rows corresponding to objects
% separable_fraction - separable fraction of data points as a function of
%           alpha
% alphas  - alpha values
% 
% Input 
% X  - is a data matrix
% Optional arguments in varargin
% 'ColinearityControl' - a real value used to select the top princinpal components (default: 10)
% 'ProjectOnSphere' - a boolean value indicating if projecting on a sphere should be performed (default: 1)
% 'Alphas' - a real vector, with alpha range, the values must be within (0,1) interval (default: [0.6,0.62,...,0.98] )
% 'ProducePlots' - a boolean value indicating if the standard plots needs
%                   to be drawn (default: 1)

ColinearityControl = 10;
ProjectOnSphere = 1;
alphas = 0.5:0.02:0.98;
ProducePlots = 1;

for i=1:length(varargin)
   if strcmpi(varargin{i},'ColinearityControl')
        ColinearityControl = varargin{i + 1};
   end
   if strcmpi(varargin{i},'ProjectOnSphere')
        ProjectOnSphere = varargin{i + 1};
   end
   if strcmpi(varargin{i},'Alphas')
        alphas = varargin{i + 1};
   end
   if strcmpi(varargin{i},'ProducePlots')
        alphas = varargin{i + 1};
   end
end

npoints = size(X,1);

Xp = preprocessing(X,1,1,1,ProjectOnSphere,'EigValueRetainingFactor',ColinearityControl);
if size(Xp, 1) > 10000
    [separable_fraction,p_alpha] = checkSeparabilityMultipleAlphaBig(Xp,alphas);
else
    xy = Xp*Xp';
    [separable_fraction,p_alpha] = checkSeparabilityMultipleAlpha(xy,alphas);
end
py_mean = mean(p_alpha,2);

[n_alpha,n_single] = dimension_uniform_sphere(py_mean,alphas,npoints);

alpha_ind_selected = find(n_single==n_alpha);

if ProducePlots

n_min = floor(min(n_alpha));
n_max = floor(max(n_alpha)+0.8);
if n_min==0
    n_min = 1;
end
if n_min>1
    n_min = n_min-1;
end
if n_min>1
    n_min = n_min-1;
end
n_max = n_max+2;
ns = n_min:n_max;
    
subplot(1,3,1);
plot(alphas,n_alpha,'ko-'); hold on; plot(alphas(alpha_ind_selected),n_single,'rx','MarkerSize',20);
xlabel('Alpha','FontSize',16); ylabel('Effective dimension','FontSize',16);

nbins = floor(npoints/200);
if nbins<20
    nbins = 20;
end
subplot(1,3,2);
hist(p_alpha(alpha_ind_selected,:),nbins); 
xlabel(sprintf('unseparability prob.p for alpha=%2.2f',alphas(alpha_ind_selected)),'FontSize',16); ylabel('Number of values','FontSize',16);

subplot(1,3,3);
semilogy(alphas,py_mean,'bo-','LineWidth',3); hold on;
xlabel('Alpha','FontSize',16); ylabel('Mean unseparability prob.','FontSize',16);
title(sprintf('Theor.curves for n=%i:%i',n_min,n_max),'FontSize',16);

pteor = zeros(length(ns),length(alphas));
for k=1:length(ns)
    for j=1:length(alphas)
        pteor(k,j) = probability_unseparable_sphere(alphas(j),ns(k));
    end
end
semilogy(alphas,pteor,'-');

set(gcf,'Position',[105         161        1115         344]);

end

end
