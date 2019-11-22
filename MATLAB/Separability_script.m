%tab = importdata('_input.txt');
%x = tab.data;
x = load('_input.txt');

[n_alpha, n_single, p_alpha, alphas, separable_fraction,Xp] = SeparabilityAnalysis(x,'ConditionalNumber',1000000,'ncomp',1,'ProducePlots',0);

fid = fopen('_palpha.txt','w');
for i=1:size(p_alpha,2)
    for j=1:size(p_alpha,1)
        fprintf(fid,'%f\t',p_alpha(j,i));
    end
    fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('_nalpha.txt','w');
fprintf(fid,'%f\n',n_alpha);
fclose(fid);

fid = fopen('_nalpha_single.txt','w');
fprintf(fid,'%f\n',n_single);
fclose(fid);

exit();