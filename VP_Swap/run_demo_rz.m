% function run_demo()
%
% This demo generates synthetic VPA data and then fits the VPA model to it.
% 
% After fitting the model, the maximum-likelihood parameter estimates, AIC, 
% and BIC estimates will be shown. In addition, a plot will be generated 
% that shows the estimation error distribution of the synthetic data and 
% model fit (in this plot, errors and fit are pooled across set sizes).
%
% This code accompanies the paper "Conceptualizing and testing working 
% memory models in a three-dimensional model space" by Van den Berg, Awh,
% and Ma, published in Psychological Review, 2014.
%
% For questions/bug reports/etc, please email nronaldvdberg@gmail.com

function run_demo_rz

%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Generate synthetic dataset  %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

% set number of trials and set size(s)
n_trials = 25000; % number of trials to generate
N_vec = [1 2 3 4 6 8]; % set sizes included in synthetic data

% draw random parameter values
J1bar = rand*20+40;
power = -2-rand;
tau = rand*10+20;
kappa_r = rand*20+40; 

% J1bar = 40;
% power = -2;
% tau = 20;
% kappa_r = 40; 

if numel(N_vec)==1
    power = NaN;
end

pars_in = [J1bar, power, tau, kappa_r];
parnames = {'J1','power','tau','kappa_r'};

% show input parameters
fprintf('\nGenerating VPA data with the following parameter values:\n');
for ii=1:length(parnames)
    fprintf('%s = %2.2f\n',parnames{ii},pars_in(ii));
end

data = gen_fake_VPA_data(pars_in,n_trials, N_vec);
fprintf('\nDone.\n\n');

%-%-%-%-%-%-%
% Fit model %
%-%-%-%-%-%-%
fprintf('Fitting VPA model to data');
opt.x0 = [40, -2, 20, 40]; % initial guess
opt.PLB = [0,-10, 0, 0];
opt.PUB = [200, 2, 200, 200];
opt.LB = [0, -20, 0, 0];
opt.UB = [500, 5, 500, 500];
opt.options = bads('defaults');
opt.options.MaxIter = '500*nvars';
opt.data = data;

[fitpars, max_lh, AIC, BIC] = fit_VPA_model_rz(opt);
fprintf('\nDone. Computing predicted response distribution...\n');
data_fit = gen_fake_VPA_data(fitpars,1e5,N_vec);
    
fprintf('\n-------------RESULTS---------------\n');
% show input parameters
fprintf('Data were generated from VPA model with the following parameter values:\n');
for ii=1:length(parnames)
    fprintf(' %s = %2.2f\n',parnames{ii},pars_in(ii));
end
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