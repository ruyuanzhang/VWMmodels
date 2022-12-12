function data = prepVWMdatant(data)
% prepare data to fit swap error
%
data.unique_N = unique(data.N);

% if the input error range is [-90, 89]
% convert to radians to [-pi, pi)
data.error = data.error * 2 * pi/180;
data.distrError = cellfun(@(x) x*2*pi/180, data.distrError, 'UniformOutput', false); % [-pi,pi]

% discretesize of the error space
% this step is necessary if errors are continous...
error_range = linspace(-pi,pi, 181); %
error_range = error_range(1:end-1) + diff(error_range(1:2))/2; %
interpFun = @(x,y,z) interp1(x,y,z,'nearest', 'extrap');
for ii=1:length(data.unique_N) %loop set size
    trial_idx = find(data.N==data.unique_N(ii));
    data.error_idx{ii} = interp1(error_range, 1:length(error_range), data.error(trial_idx),'nearest','extrap');
    
    % index non-target swap error
    if data.unique_N(ii)<=1
        data.distrError_idx{ii} = {};
    else
        data.distrError_idx{ii} = cellfun(interpFun, repmat({error_range}, 1, length(trial_idx)), repmat({1:length(error_range)},1,length(trial_idx)),...
            data.distrError(trial_idx), 'UniformOutput', false);
    end
end

data.gvar.error_range = error_range;
return