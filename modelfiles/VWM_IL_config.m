function c = VWM_IL_config(varparams)

% varparams is a struct that accept input from outside functions

nFit = varparams.nFit;
optimizer = varparams.optimizer;

% some setting
maxJ1bar = 700;

%% define optimization variable
opt.PLB = [0, 0];
opt.PUB = [maxJ1bar, 10];
opt.LB = [0, 0];
opt.UB = [maxJ1bar, 20];
opt.paramLabels = {'kappa_r', 'K'};

nvars = length(opt.PLB);
opt.options = bads('defaults');
opt.options.MaxIter = maxJ1bar*nvars;
% initial guess
opt.x0 = [opt.PLB(1):(opt.PUB(1)-opt.PLB(1))/(nFit-1):opt.PUB(1);
    opt.PLB(2):(opt.PUB(2)-opt.PLB(2))/(nFit-1):opt.PUB(2);
    ];
opt.x0=opt.x0';

opt.nFit = nFit;
opt.nvars = nvars;
opt.optimizer = optimizer;
c.opt=opt;
%% define fitting function
c.fitFun = @VWM_IL_fit;

%% define negative loglikelihood function
c.negLogLikeliFun=@VWM_IL_nll;


end