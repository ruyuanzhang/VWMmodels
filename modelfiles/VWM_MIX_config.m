function c = VWM_MIX_config(varparams)

% varparams is a struct that accept input from outside functions

nFit = varparams.nFit;
optimizer = varparams.optimizer;
nSetSize = varparams.nSetSize;

% some setting
maxJ1bar = 700;

%% define optimization variable
opt.PLB = [zeros(1, nSetSize), 0];
opt.PUB = [maxJ1bar*ones(1, nSetSize), 10];
opt.LB = opt.PLB;
opt.UB = [maxJ1bar*ones(1, nSetSize), 20];
opt.paramLabels = {'kappa_r', 'K'};

nvars = length(opt.PLB);
opt.options = bads('defaults');
opt.options.MaxIter = maxJ1bar*nvars;
% initial guess
x0 = nan(nSetSize+1, nFit);
for i=1:nSetSize+1
    x0(i,:) = opt.PLB(i):(opt.PUB(i)-opt.PLB(i))/(nFit-1):opt.PUB(i);
end
opt.x0=x0';

opt.nFit = nFit;
opt.nvars = nvars;
opt.optimizer = optimizer;

c.opt=opt;
%% define fitting function
c.fitFun = @VWM_MIX_fit;

%% define negative loglikelihood function
c.negLogLikeliFun=@VWM_MIX_nll;


end