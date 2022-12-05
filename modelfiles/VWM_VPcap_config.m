function c = VWM_VPcap_config(varparams)

% varparams is a struct that accept input from outside functions
nFit = varparams.nFit;
optimizer = varparams.optimizer;

% some setting
maxJ1bar = 700;

%% define optimization variable
opt.PLB = [0, 0, 0, 0, 0];
opt.PUB = [maxJ1bar, 10, maxJ1bar, maxJ1bar, 20];
opt.LB = [0, 0, 0, 0, 0];
opt.UB = [maxJ1bar, 20, maxJ1bar, maxJ1bar, 20];
opt.paramLabels = {'J1bar', 'power', 'var', 'kappa_r', 'K'};

% 
nvars = length(opt.PLB);
opt.options = bads('defaults');
opt.options.MaxIter = maxJ1bar*nvars;

x0 = [opt.PLB(1)+eps:(opt.PUB(1)-opt.PLB(1))/(nFit-1):opt.PUB(1)-eps;
    opt.PLB(2)+eps:(opt.PUB(2)-opt.PLB(2))/(nFit-1):opt.PUB(2)-eps;
    opt.PLB(3)+eps:(opt.PUB(3)-opt.PLB(3))/(nFit-1):opt.PUB(3)-eps;
    opt.PLB(4)+eps:(opt.PUB(4)-opt.PLB(4))/(nFit-1):opt.PUB(4)-eps;
    opt.PLB(5)+eps:(opt.PUB(5)-opt.PLB(5))/(nFit-1):opt.PUB(5)-eps];
opt.x0=x0';  %x0,  nFit x n params

opt.nFit = nFit;
opt.nvars = nvars;
opt.optimizer = optimizer;

c.opt=opt;
%% define fitting function
c.fitFun = @VWM_VPcap_fit;

%% define negative loglikelihood function
c.negLogLikeliFun=@VWM_VPcap_nll;

end