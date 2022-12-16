function c = VWM_VP_config(varparams)

% varparams is a struct that accept input from outside functions

nFit = varparams.nFit;
optimizer = varparams.optimizer;

% some setting
maxJ1bar = 700;

%% define optimization variable
opt.PLB = [0, 0, 0, 0];
opt.PUB = [maxJ1bar, 5, maxJ1bar, maxJ1bar];
opt.LB = [0, 0, 0, 0];
opt.UB = [maxJ1bar, 5, maxJ1bar, maxJ1bar];
opt.paramLabels = {'J1bar', 'power', 'tau', 'kappa_r'};

% 
nvars = length(opt.PLB);
opt.options = bads('defaults');
opt.options.MaxIter = maxJ1bar*nvars;

x0 = [opt.PLB(1)+(opt.PUB(1)-opt.PLB(1))*rand(1,nFit);
    opt.PLB(2)+(opt.PUB(2)-opt.PLB(2))*rand(1,nFit);
    opt.PLB(3)+(opt.PUB(3)-opt.PLB(3))*rand(1,nFit);
    opt.PLB(4)+(opt.PUB(4)-opt.PLB(4))*rand(1,nFit)];

opt.x0=x0';  %x0,  nFit x n params

opt.nFit = nFit;
opt.nvars = nvars;
opt.optimizer = optimizer;

c.opt=opt;
%% define fitting function
c.fitFun = @VWM_VP_fit;

%% define negative loglikelihood function
c.negLogLikeliFun=@VWM_VP_nll;

end