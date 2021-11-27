%%
cl;
d = load('LYF_sham_Day3_set1__set3__set5__set8_2110201359.mat');
d = d.results;

%% 
data.probe = d.probeInd;
data.resp = d.respInd;
data.error = d.error;
data.N = d.stimNum;

% extractor
for i = 1:length(data.probe)
    data.distr{i} = deleteel(d.stimuliInd(i,:), find((data.probe(i)==d.stimuliInd(i,:)|isnan(d.stimuliInd(i,:)))));
    if ~isempty(data.distr{i})
            data.distrError{i} = circulardiff(data.resp(i),data.distr{i}, 180); % error range, -90~89;
    end
end

%%
cd ..
save('sampledata2','data');