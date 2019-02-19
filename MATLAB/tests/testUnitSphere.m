%clear all;
npoints = 10000;

dim_min = 30;
dim_max = 30;

for dimension=dim_min:dim_max

disp(dimension);
    
X = randsphere(npoints,dimension,1);
%X = spheredist(npoints,dimension);

Xp = preprocessing(X,1,1,1,1,'EigValueRetainingFactor',10);
xy = Xp*Xp';

nbins = 10;
lowalfa = 0.6;

clear sf py pteor alfa;

for i=1:nbins
alfa(i)=lowalfa+(1-lowalfa)*(i-1)/nbins;
pteor(i) = probability_unseparable_sphere(alfa(i),dimension);
end

tic;
[sf,py] = checkSeparabilityMultipleAlpha(xy,alfa);
toc;
%for i=1:length(alfa)
%    [sf(i),py(i,:)] = checkSeparability(xy,alfa(i));
%end

py_mean = mean(py,2);
[dims(dimension-dim_min+1,:),dims_single(dimension-dim_min+1),alfa_single(dimension-dim_min+1)] = dimension_uniform_sphere(py_mean,alfa,npoints);

%plot(alfa,1-sf,'b-'); hold on;
%plot(alfa,fracteor,'r-');



%plot(alfa,fract_emp,'b-'); hold on;
%plot(alfa,fracteor,'r-');

semilogy(alfa,py_mean,'bo-','MarkerSize',3); hold on;
semilogy(alfa,pteor,'r-');

%loglog(alfa,py_mean,'bo-','MarkerSize',3); hold on;
%loglog(alfa,pteor,'r-');


%  lpy_mean = log10(py_mean);
%  x = log10(1-alfa.*alfa); y = lpy_mean+log10(alfa);
%  x = x'; y = y';
%  inds = find(~isinf(y));
%  x = x(inds); y = y(inds);
%  Xc = [ones(length(x),1) x];
%  b(dimension-dim_min+1,:) = Xc\y;

end

% xx = dim_min:dim_max; yy = b(:,2); xx = xx'; yy = yy';
% xxc = [ones(length(xx),1) xx];
% b2 = xxc\yy';
% disp(sprintf('A(n)=%f+%fn',b2(1),b2(2)));
% 
% xx = dim_min:dim_max; xx=log10(xx); yy = b(:,1); xx = xx'; yy = yy';
% xxc = [ones(length(xx),1) xx];
% b1 = xxc\yy';
% disp(sprintf('B(n)=%f+%flog10(n)',b1(1),b1(2)));


title('Non-separability p_y vs alfa'); 
xlabel('alfa','FontSize',14); 
%ylabel('Unseparable fraction of points','FontSize',14);
ylabel('Unseparability probability p_y','FontSize',14);

figure;
plot(alfa,dims,'ko-'); hold on;
%plot(alfa_single,dims_single,'rx','MarkerSize',10);
for i=1:length(alfa_single)
    text(alfa_single(i),dims_single(i),sprintf(' %i',floor(dims_single(i)+0.5)),'Color','r','FontSize',10);
end
ylabel('Estimated effective dimension','FontSize',14);
xlabel('alpha','FontSize',14);