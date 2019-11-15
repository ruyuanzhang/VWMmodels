clear all;close all;clc;



%% data load and prepration
tmp = load('SUBJ02behavior.mat');
tmp = tmp.a;

tmp.allWMload = flatten(cell2mat(tmp.WMload));
tmp.alltarg = flatten(cell2mat(tmp.targ));
tmp.allresp = flatten(cell2mat(tmp.resp));

idx = ~isnan(tmp.allWMload);
tmp.allWMload = tmp.allWMload(idx); 


tmp.alltarg(tmp.alltarg<=0) = tmp.alltarg(tmp.alltarg<=0) + 180;
tmp.alltarg = (tmp.alltarg - 90)/180*pi; %[-pi/2,pi/2] 
tmp.allresp(tmp.allresp<=0) = tmp.alltarg(tmp.allresp<=0) + 180;
tmp.allresp = (tmp.allresp - 90)/180*pi; %[-pi/2,pi/2] 
tmp.allerror = circulardiff(tmp.allresp, tmp.alltarg, pi);

tmp.allerror = tmp.allerror(idx); 
a = tmp;
save('SUBJ02behavior.mat','a')

data.N = tmp.allWMload;
data.error_vec = tmp.allerror;
%%
N_vec = [1 2];
parnames = {'J1','power','tau','kappa_r'};

[fitpars, max_lh, AIC, BIC] = fit_VPA_model_rz(data);

fprintf('\nDone. Computing predicted response distribution...\n');
data_fit = gen_fake_VPA_data(fitpars,1e5,N_vec);
    
fprintf('\n-------------RESULTS---------------\n');
% show input parameters
fprintf('Results from fitting these data with VPA model:\n');
for ii=1:length(parnames)
    fprintf(' %s = %2.2f\n',parnames{ii},fitpars(ii));
end
fprintf('\nmax log likelihood=%2.2f, AIC=%2.2f, BIC=%2.2f\n',max_lh,AIC,BIC);
fprintf('-----------------------------------\n');

% plot fit
figure
X = linspace(-pi,pi,52);
X = X(1:end-1)+diff(X(1:2))/2;
Y_emp = hist(data.error_vec,X);
Y_emp = Y_emp/sum(Y_emp)/diff(X(1:2));
Y_fit = hist(data_fit.error_vec,X);
Y_fit = Y_fit/sum(Y_fit)/diff(X(1:2));
bar(X,Y_emp,'k');
hold on
plot(X,Y_fit,'r-','Linewidth',3)
legend('Data','Fit')
xlabel('Response error');
ylabel('Probability');
xlim([-pi pi]);
title(['AIC=' num2str(AIC,4)]);