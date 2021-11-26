function c = VWM_MIX_config(nSetSize, nFit)

% 
if notDefined('nSetSize')
    error('For MIX model, you need to specifiy level of set size levels')
end
if notDefined('nFit')
    nFit = 20;
end

% some setting
maxJ1bar = 700;

%% define optimization variable
opt.PLB = [zeros(1, nSetSize), 0, 0];
opt.PUB = [maxJ1bar*ones(1, nSetSize), 20, 1];
opt.LB = opt.PLB;
opt.UB = [maxJ1bar*ones(1, nSetSize), 20, 1];
opt.paramLabels = {'kappa_r', 'K', 'pNT'};

nvars = length(opt.PLB);
opt.options = bads('defaults');
opt.options.MaxIter = maxJ1bar*nvars;
% initial guess
x0 = nan(nSetSize+2, nFit);
for i=1:nSetSize+2
    x0(i,:) = opt.PLB(i):(opt.PUB(i)-opt.PLB(i))/(nFit-1):opt.PUB(i);
end
opt.x0=x0';

opt.nFit = nFit;
opt.nvars = nvars;

c.opt=opt;
%% define fitting function
c.fitFun = @VWM_MIXnt_fit;

%% define negative loglikelihood function
c.negLogLikeliFun=@VWM_MIXnt_nll;


end